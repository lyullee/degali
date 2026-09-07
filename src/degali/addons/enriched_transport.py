"""Gauge-free, weak modal transport with newly reconstructed mass and stress.

A scalar mixing profile is REQUIRED. This provides conditional local rates,
not a downstream eddy-mixing closure or a full curved RANS model.
"""

import math
import copy
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import legvander, legval, legint
from .phase_radial_quadrature import gauss_rule
from .edge_conservative_refit import FaceSplitSquareMoments
from .reservoir_thermal import _PhaseForceView
from .energy_crosswind import IndependentEnergyCrosswind
from .buoyancy_profile import BuoyancyConstrainedEnthalpySection
from .transverse_mixing import ConservativeTransverseMixing


def advective_moments(projection, parameters, *, order=8, angular_order=8):
    """Independent physical flux/modal values, without any rate reconstruction."""
    result=np.zeros(5+projection.count)
    theta=projection.mixing.state[3]
    for prepared in FaceSplitSquareMoments(projection).rule(parameters,order=order,angular_order=angular_order):
        d=projection.fields(parameters,prepared)
        rho,c,h,u,w=d["rho"],d["c"],d["h"],prepared["u"],prepared["weights"]
        result[:5] += w@np.column_stack([rho*u,c*u,rho*u*u*math.cos(theta),rho*u*u*math.sin(theta),h*u+.5*rho*u**3])
        size=projection.basis.size
        result[5:5+size] += (w*c*u)@prepared["psi"]
        result[5+size:] += (w*h*u)@prepared["psi"]
    return result


def shifted_advective_moments(projection, parameters, rates, ds, *, order=8, angular_order=8):
    """Verification-only perturbed fields; no valid downstream state is implied."""
    old=projection.mixing.state
    encoded=np.array([math.log(old[0]),math.log(old[0]*old[1]),math.log(old[2]),old[3],math.log(old[4]),old[5],old[6]])
    new=encoded+ds*rates[:7]
    state=np.array([math.exp(new[0]),math.exp(new[1]-new[0]),math.exp(new[2]),new[3],math.exp(new[4]),new[5],new[6]])
    original=projection.section
    section=BuoyancyConstrainedEnthalpySection(original.jetplume,original.thermodynamics,
        thermal_width_ratio=original.thermal_width_ratio,quadrature_points=original.quadrature_points)
    section._quadrature_cache=original._quadrature_cache
    m=ConservativeTransverseMixing(section,state)
    view=copy.copy(projection)
    view.section,view.mixing=section,m
    view.hc=section.phase_inverse.enthalpy_and_slope(state[0],state[0]*state[1])[0]
    # Only prepare/fields/rule are used on this local view. Its old baseline
    # is NOT evaluated or represented as a newly valid weak transport state.
    return advective_moments(view,parameters+ds*rates[7:],order=order,angular_order=angular_order)


@lru_cache(maxsize=8)
def panel_antiderivative(order):
    x, w = gauss_rule(order)
    v = legvander(2*x-1., order-1)
    primitive = []
    for j in range(order):
        coefficients = np.zeros(order)
        coefficients[j] = 1.
        anti = legint(coefficients)
        primitive.append(.5*(legval(2*x-1., anti)-legval(-1., anti)))
    return np.column_stack(primitive)@np.linalg.inv(v)


def panel_cumulative(values, knots, order):
    """Within-panel indefinite integration, plus exact Gauss panel prefixes."""
    values = np.asarray(values)
    _, w = gauss_rule(order)
    shaped = values.reshape((len(knots)-1, order)+values.shape[1:])
    width = np.diff(knots)
    full = np.tensordot(shaped, w, axes=(1, 0))*width.reshape((-1,)+(1,)*(values.ndim-1))
    prefix = np.concatenate([np.zeros_like(full[:1]), np.cumsum(full, axis=0)], axis=0)
    partial = np.einsum("ij,pj...->pi...", panel_antiderivative(order), shaped)
    partial *= width.reshape((-1, 1)+(1,)*(values.ndim-1))
    return (partial+prefix[:-1, None]).reshape(values.shape), prefix[-1]


class PrescribedRadialMixing:
    """Explicit piecewise mixing input; no inference from the enriched flow."""
    def __init__(self, knots, values):
        self.knots = np.array(knots, float)
        values = np.asarray(values, float)
        if (values.ndim != 3 or values.shape[0] != len(knots)-1 or values.shape[2] != 2
                or np.any(np.diff(self.knots) <= 0.) or not np.all(np.isfinite(values))):
            raise ValueError("supply ordered phase panels and chi/old-mass Gauss values")
        self.order = values.shape[1]
        x, _ = gauss_rule(self.order)
        v = legvander(2*x-1., self.order-1)
        self.coefficients = np.einsum("ij,pjk->pik", np.linalg.inv(v), values)
        if np.any(values[:, :, 0] <= 0.):
            raise ValueError("prescribed mixing must be strictly positive")

    @classmethod
    def from_weak_baseline(cls, projection, *, order):
        m, baseline = projection.mixing, projection.baseline
        q, _ = m.partition.nodes(order)
        data = m.mixing_family(q.ravel(), baseline["family"], order=16)
        at = np.array([1., baseline["gamma"]])
        values = np.column_stack([data["chi"]@at, data["mass"]@at]).reshape(*q.shape, 2)
        return cls(m.partition.knots, values)

    def evaluate(self, q):
        q = np.asarray(q, float)
        if q.ndim != 1 or not np.all(np.isfinite(q)) or np.any(q < self.knots[0]) or np.any(q > self.knots[-1]*(1+1e-14)):
            raise ValueError("mixing samples are outside the supplied radial domain")
        index = np.clip(np.searchsorted(self.knots, q, side="right")-1, 0, len(self.knots)-2)
        x = 2*(q-self.knots[index])/(self.knots[index+1]-self.knots[index])-1.
        out = np.column_stack([legval(x, self.coefficients[index, :, k].T, tensor=False) for k in (0, 1)])
        if np.any(out[:, 0] <= 0.):
            raise ValueError("interpolated mixing is not positive; refine or reject input")
        return out


class EnrichedModalTransport:
    def __init__(self, projection, parameters, *, scalar_mixing, thermal_species_ratio, mechanical_work):
        if not isinstance(scalar_mixing, PrescribedRadialMixing):
            raise TypeError("supply an explicit phase-resolved scalar mixing input")
        if not math.isfinite(thermal_species_ratio) or thermal_species_ratio <= 0.:
            raise ValueError("supply a positive thermal/species diffusivity ratio")
        if mechanical_work != "reduced_buoyancy_work_immediate_shear_heat":
            raise ValueError("explicit reduced work and immediate shear-heating assumptions required")
        self.projection, self.parameters = projection, np.array(parameters, float)
        self.mixing, self.ratio = scalar_mixing, thermal_species_ratio
        self.quadrature = FaceSplitSquareMoments(projection)
        self.size, self.count = projection.basis.size, 7+projection.count
        self.area = projection.mixing.area
        self.active = np.r_[np.arange(5), np.arange(7, self.count)]
        # q lies exactly in the centre-zero polynomial span. This is a gauge
        # identity, not a fitted thermal transport coefficient.
        x, w = gauss_rule(16)
        a, b = np.meshgrid(x, x, indexing="ij")
        prep = projection.prepare(a.ravel(), b.ravel())
        self.q_mode = prep["psi"].T@(np.outer(w, w).ravel()*prep["q"])

    def old_rate_embedding(self, old_rates):
        rates = np.zeros(self.count)
        rates[:7] = old_rates[:7]
        rates[7+self.size:] = 2*self.q_mode/self.projection.section.thermal_width_ratio**2*old_rates[7]
        return rates

    def local(self, a, b):
        p, m = self.projection, self.projection.mixing
        prep = p.prepare(a, b)
        d = p.fields(self.parameters, prep)
        rho, c, h, u = d["rho"], d["c"], d["h"], prep["u"]
        n, size = len(a), self.size
        dc, dh, du = [np.zeros((n, self.count)) for _ in range(3)]
        dc[:, 1] = c
        dc[:, 7:7+size] = c[:, None]*prep["psi"]
        rc, yc = m.state[:2]
        hr0, hc0 = p.section._phase_partials(np.asarray(rc), np.asarray(rc*yc))
        dh[:, 0], dh[:, 1] = hr0*rc*h/p.hc, hc0*rc*yc*h/p.hc
        dh[:, 7+size:] = h[:, None]*prep["psi"]
        drho = (dh-d["hc"][:, None]*dc)/d["hr"][:, None]
        du[:, :7] = m.wind_partials[:7]*math.cos(m.state[3])
        du[:, 3] -= m.wind*math.sin(m.state[3])
        du[:, 4] += m.state[4]*np.exp(-p.section.velocity_shape_exponent*prep["q"])
        coordinate = np.column_stack([a, b])
        grad_psi = np.stack([prep["da"], prep["db"]], axis=1)
        grad_c = c[:, None]*(-2*p.q0*coordinate+np.einsum("nik,k->ni", grad_psi, self.parameters[:size]))
        grad_h = h[:, None]*(-2*p.q0/p.section.thermal_width_ratio**2*coordinate
                            +np.einsum("nik,k->ni", grad_psi, self.parameters[size:]))
        grad_rho = (grad_h-d["hc"][:, None]*grad_c)/d["hr"][:, None]
        grad_y = (grad_c-d["y"][:, None]*grad_rho)/rho[:, None]
        specific_h = h/rho
        grad_specific_h = (grad_h-specific_h[:, None]*grad_rho)/rho[:, None]
        uq = -p.section.velocity_shape_exponent*m.state[4]*np.exp(-p.section.velocity_shape_exponent*prep["q"])
        bm = self.area*(drho*u[:, None]+rho[:, None]*du)
        bc = self.area*(dc*u[:, None]+c[:, None]*du)
        bp = self.area*(drho*(u*u)[:, None]+(2*rho*u)[:, None]*du)
        bh = self.area*(dh*u[:, None]+h[:, None]*du)
        bk = .5*self.area*(drho*(u**3)[:, None]+(3*rho*u*u)[:, None]*du)
        for result, value in ((bm,rho*u),(bc,c*u),(bp,rho*u*u),(bh,h*u),(bk,.5*rho*u**3)):
            result[:, 2] += self.area*value
        force = 9.81*self.area*(p.section.rhoa-rho)*math.sin(m.state[3])
        return dict(**prep, **d, dc=dc, dh=dh, drho=drho, du=du, coordinate=coordinate,
            grad_psi=grad_psi, grad_y=grad_y, grad_h=grad_specific_h, specific_h=specific_h,
            uq=uq, bm=bm, bc=bc, bp=bp, bh=bh, bk=bk, force=force)

    def _ray(self, angle, knots, order):
        x, w = gauss_rule(order)
        q = (knots[:-1,None]+np.diff(knots)[:,None]*x).ravel()
        radius = np.sqrt(q/self.projection.q0)
        d = self.local(np.minimum(radius*math.cos(angle),1.), np.minimum(radius*math.sin(angle),1.))
        integrand = np.zeros((len(q), 2, self.count+1))
        integrand[:, 0, 1:] = -d["bm"]
        integrand[:, 1, 0], integrand[:, 1, 1:] = d["force"], -d["bp"]
        cumulative, endpoint = panel_cumulative(integrand, knots, order)
        f = cumulative/(2*q[:,None,None])
        edge = endpoint/(2*knots[-1])
        return d, f[:,0], f[:,1], edge, (np.diff(knots)[:,None]*w).ravel()

    def physical_velocity(self,a,b,rates,mass_flux):
        d=self.local(a,b)
        m=self.projection.mixing
        ly,ln=math.sqrt(2*self.projection.q0)*np.array([m.sy,m.sn])
        metric=1.-rates[3]*ln*b
        if np.any(metric<=0.):
            raise ValueError("curved physical section metric is not positive")
        common=np.asarray(mass_flux)/(self.area*d["rho"])
        ay=m.log_sy_partials[:7]@rates[:7]
        an=m.log_sn_partials[:7]@rates[:7]
        return ly*a/metric*(d["u"]*ay+common),ln*b/metric*(d["u"]*an+common)

    def assemble(self, *, order=4, angular_order=8, probes=1025):
        p, m, size, count = self.projection, self.projection.mixing, self.size, self.count
        ak = self.quadrature.angular_knots(self.parameters, probes=probes)
        crossed = self.mixing.knots[(self.mixing.knots>p.q0)&(self.mixing.knots<2*p.q0)]
        ak = np.sort(np.r_[ak, np.arccos(np.sqrt(p.q0/crossed))])
        ak = ak[np.r_[True, np.diff(ak)>2e-13]]
        x, w = gauss_rule(angular_order)
        angles = (ak[:-1,None]+np.diff(ak)[:,None]*x).ravel()
        aw = (np.diff(ak)[:,None]*w).ravel()
        jf = np.zeros((5,count))
        jm = np.zeros((2*size,count))
        rhs = np.zeros((2*size,count+1))
        flux = np.zeros(5+2*size)
        force, work = 0.,0.
        heat0 = np.zeros(count+1)
        heat_axial,heat_production,heat_boundary=[np.zeros(count+1) for _ in range(3)]
        boundary_mass = np.zeros(count+1)
        samples, edges = [], []
        for first in range(0,len(angles),8):
            subset = angles[first:first+8]
            partitions = self.quadrature.partitions(self.parameters,subset)
            for j, (angle,k) in enumerate(zip(subset,partitions)):
                extra = self.mixing.knots[(self.mixing.knots>0.)&(self.mixing.knots<k[-1])]
                k = np.sort(np.r_[k,extra])
                k = k[np.r_[True,np.diff(k)>k[-1]*2e-13]]
                d,fm,fp,fe,rw = self._ray(angle,k,order)
                weight = 8*aw[first+j]*rw
                rho,c,h,u,q = [d[s] for s in ("rho","c","h","u","q")]
                psi = d["psi"]
                chi = self.mixing.evaluate(q)[:,0]
                rp = fp-u[:,None]*fm
                production = -(2*q*d["uq"])[:,None]*rp
                dot = np.einsum("nik,ni->nk",d["grad_psi"],d["coordinate"])
                diffusion_c = np.einsum("nik,ni->nk",d["grad_psi"],d["grad_y"])
                diffusion_h = np.einsum("nik,ni->nk",d["grad_psi"],d["grad_h"])
                rhs[:size] += np.einsum("n,nk,nr->kr",weight*d["y"],dot,fm)
                rhs[size:] += np.einsum("n,nk,nr->kr",weight*d["specific_h"],dot,fm)
                rhs[size:] += np.einsum("n,nk,nr->kr",weight,psi,production)
                rhs[:size,0] -= (weight*self.area*rho*chi/(2*p.q0))@diffusion_c
                rhs[size:,0] -= (weight*self.area*rho*chi*self.ratio/(2*p.q0))@diffusion_h
                jm[:size] += np.einsum("n,nk,nr->kr",weight,psi,d["bc"])
                jm[size:] += np.einsum("n,nk,nr->kr",weight,psi,d["bh"])
                theta = m.state[3]
                px,pz = d["bp"]*math.cos(theta),d["bp"]*math.sin(theta)
                px[:,3] -= self.area*rho*u*u*math.sin(theta)
                pz[:,3] += self.area*rho*u*u*math.cos(theta)
                jf += np.einsum("n,nkr->kr",weight,np.stack([d["bm"],d["bc"],px,pz,d["bh"]+d["bk"]],axis=1))
                flux[:5] += self.area*weight@np.column_stack([rho*u,c*u,rho*u*u*math.cos(theta),rho*u*u*math.sin(theta),h*u+.5*rho*u**3])
                flux[5:5+size] += (self.area*weight*c*u)@psi
                flux[5+size:] += (self.area*weight*h*u)@psi
                force += float(weight@(9.81*self.area*(p.section.rhoa-rho)))
                work += float(weight@(d["force"]*u))
                heat0[1:] += weight@d["bh"]
                heat0 -= weight@production
                heat_axial[1:] += weight@d["bh"]
                heat_production += weight@production
                t = math.tan(angle)
                ed = self.local(np.array([1.]),np.array([t]))
                ue = ed["u"][0]
                fk = ue*fe[1]-.5*ue*ue*fe[0]
                fh = .5*m.wind*m.wind*fe[0]-fk
                ew = 16*p.q0*aw[first+j]/math.cos(angle)**2
                rhs[size:] -= ew*ed["psi"][0,:,None]*fh
                heat0 += ew*fh
                heat_boundary += ew*fh
                boundary_mass += ew*fe[0]
                idx = np.unique(np.linspace(0,len(q)-1,17,dtype=int))
                samples.append(dict(q=q[idx],rho=rho[idx],uq=d["uq"][idx],rp=rp[idx],fm=fm[idx]))
                edges.append(dict(t=t, fm=fe[0],fp=fe[1],h_target=fh,
                    rho=ed["rho"][0],y=ed["y"][0],specific_h=ed["specific_h"][0],u=ue,
                    grad_y=ed["grad_y"][0,0],grad_h=ed["grad_h"][0,0]))
        source = IndependentEnergyCrosswind.source_terms(_PhaseForceView(m,force),m.state)
        source[4] += work
        matrix = np.vstack([jf,jm-rhs[:,1:]])
        right = np.r_[source,rhs[:,0]]
        rates = np.zeros(count)
        rates[5:7] = math.cos(m.state[3]),math.sin(m.state[3])
        scales = np.maximum(np.max(abs(matrix[:,self.active]),axis=1),1e-12)
        scaled = matrix[:,self.active]/scales[:,None]
        singular = np.linalg.svd(scaled,compute_uv=False)
        if singular[-1] <= singular[0]*1e-12:
            raise ValueError("gauge-free modal transport matrix is numerically rank deficient")
        rates[self.active] = np.linalg.solve(scaled,(right-matrix@rates)/scales)
        residual = float(max(abs(matrix@rates-right)/np.maximum(abs(right),1.)))
        at = np.r_[1.,rates]
        sampled_chi = np.concatenate([-s["rp"]@at/(self.area*s["rho"]*s["uq"]) for s in samples])
        fm = np.array([e["fm"]@at for e in edges])
        fp = np.array([e["fp"]@at for e in edges])
        targets = np.array([e["h_target"]@at for e in edges])
        t = np.array([e["t"] for e in edges])
        q = p.q0*(1+t*t)
        chi,old_fm = self.mixing.evaluate(q).T
        old = m.local(q)
        ec = np.array([e["y"] for e in edges])*fm-self.area*np.array([e["rho"]*e["grad_y"] for e in edges])*chi/(2*p.q0)
        eh = np.array([e["specific_h"] for e in edges])*fm-self.area*np.array([e["rho"]*e["grad_h"] for e in edges])*chi*self.ratio/(2*p.q0)-targets
        ep = fp-m.wind*math.cos(m.state[3])*fm
        defects = dict(hydrogen=float(max(abs(ec)/np.maximum(abs(old["y"]*old_fm),1e-12))),
            heat=float(max(abs(eh)/np.maximum(abs(old["h"]/old["rho"]*old_fm),1.))),
            momentum=float(max(abs(ep)/np.maximum(abs(old["u"]*old_fm),1.))))
        curvature = abs(rates[3])*math.sqrt(2*p.q0)*m.sn
        heat_terms=dict(axial=float(heat_axial@at),production=float(heat_production@at),boundary=float(heat_boundary@at))
        heat_scaled=abs(float(heat0@at))/max(1.,*(abs(v) for v in heat_terms.values()))
        old_rates=self.old_rate_embedding(p.baseline["rates"])
        return dict(rates=rates, matrix=matrix, right=right, moment_jacobian=np.vstack([jf,jm]),
            fluxes=flux, source=source, force=force, work=work, matrix_condition=float(singular[0]/singular[-1]),
            linear_scaled_error=residual, weak_heat_zeroth=float(heat0@at),
            weak_heat_scaled_error=heat_scaled,heat_terms=heat_terms,
            frozen_rate_scaled_defect=float(max(abs(matrix@old_rates-right)/np.maximum(abs(right),1.))),
            mass_boundary_error=float(boundary_mass@at+source[0]),
            minimum_chi_momentum=float(min(sampled_chi)), maximum_outward_mass=float(max(fm)),
            curvature_half_width=float(curvature), edge_samples=len(edges), edge_defects=defects,
            constitutive_passed=bool(min(sampled_chi)>=0. and max(fm)<0. and curvature<.1
                and max(defects["hydrogen"],defects["heat"])<=.05))
