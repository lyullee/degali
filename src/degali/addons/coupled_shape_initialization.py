"""Research-only simultaneous shape/rate initialization on a fixed trial mesh.

The mesh accelerates nonlinear iterations, NOT scientific acceptance. Every
evaluation solves the full gauge-free rate system for the current shape.
Independent moving-phase quadrature is required before adopting any result.
"""

import math
import numpy as np
from .phase_radial_quadrature import gauss_rule
from .edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit
from .enriched_transport import EnrichedModalTransport, panel_cumulative
from .reservoir_thermal import _PhaseForceView
from .energy_crosswind import IndependentEnergyCrosswind


class FixedMeshCoupledTransport:
    """Vectorized counterpart of EnrichedModalTransport on a frozen fit mesh.

    Phase cuts are found on the seed, but phase inversion is evaluated on the
    CURRENT parameters. Moving phase events need independent verification.
    Extra face rays have zero volume weight and cannot alter the rate solve.
    """
    def __init__(self, projection, seed, *, scalar_mixing, thermal_species_ratio,
                 mechanical_work, order=4, angular_order=24, face_samples=33,
                 split_angles=False):
        self.prototype = EnrichedModalTransport(projection, seed,
            scalar_mixing=scalar_mixing, thermal_species_ratio=thermal_species_ratio,
            mechanical_work=mechanical_work)
        self.projection, self.mixing = projection, scalar_mixing
        self.ratio, self.mechanical_work = thermal_species_ratio, mechanical_work
        self.order, self.angular_order = order, angular_order
        quad = FaceSplitSquareMoments(projection)
        if split_angles:
            ak = quad.angular_knots(seed)
            q0 = projection.q0
            crossed = scalar_mixing.knots[(scalar_mixing.knots > q0) & (scalar_mixing.knots < 2*q0)]
            ak = np.sort(np.r_[ak, np.arccos(np.sqrt(q0/crossed))])
            ak = ak[np.r_[True, np.diff(ak) > 2e-13]]
        else:
            ak = np.array([0., math.pi/4])
        ax, aw = gauss_rule(angular_order)
        angles = (ak[:-1, None]+np.diff(ak)[:, None]*ax).ravel()
        weights = (np.diff(ak)[:, None]*aw).ravel()
        self.volume_rays = len(angles)
        if face_samples < 3:
            raise ValueError("at least three independent face rays required")
        t = .5*(1-np.cos(np.linspace(0., math.pi, face_samples)))
        self.face_slice = slice(len(angles), len(angles)+face_samples)
        self.angles = np.r_[angles, np.arctan(t)]
        self.angular_weights = np.r_[weights, np.zeros_like(t)]
        self.face_t = np.tan(self.angles)
        self.face_t[-1] = 1.
        self.knots, self.slices = [], []
        qs, phis, ws = [], [], []
        x, w = gauss_rule(order)
        offset = 0
        for start in range(0, len(self.angles), 8):
            subset = self.angles[start:start+8]
            partitions = quad.partitions(seed, subset)
            for j, k in enumerate(partitions):
                extra = scalar_mixing.knots[(scalar_mixing.knots > 0.) & (scalar_mixing.knots < k[-1])]
                k = np.sort(np.r_[k, extra])
                k = k[np.r_[True, np.diff(k) > k[-1]*2e-13]]
                q = (k[:-1, None]+np.diff(k)[:, None]*x).ravel()
                qs.append(q)
                phis.append(np.full(len(q), subset[j]))
                ws.append(8*self.angular_weights[start+j]*(np.diff(k)[:, None]*w).ravel())
                self.knots.append(k)
                self.slices.append(slice(offset, offset+len(q)))
                offset += len(q)
        self.q, phi, self.weight = np.concatenate(qs), np.concatenate(phis), np.concatenate(ws)
        radius = np.sqrt(self.q/projection.q0)
        self.a = np.minimum(radius*np.cos(phi), 1.)
        self.b = np.minimum(radius*np.sin(phi), 1.)
        self.chi = scalar_mixing.evaluate(self.q)[:, 0]
        self.edge_chi, old_fm = scalar_mixing.evaluate(projection.q0*(1+self.face_t**2)).T
        old = projection.mixing.local(projection.q0*(1+self.face_t**2))
        self.edge_scales = np.column_stack([
            np.maximum(abs(old["y"]*old_fm), 1e-12),
            np.maximum(abs(old["h"]/old["rho"]*old_fm), 1.),
            np.maximum(abs(old["u"]*old_fm), 1.)])
        self.sample_indices = np.concatenate([
            s.start+np.unique(np.linspace(0, s.stop-s.start-1, 17, dtype=int)) for s in self.slices])
        # Same finite polynomial basis and geometry on every nonlinear call.
        self.prepared = projection.prepare(self.a, self.b)
        self.face_prepared = projection.prepare(np.ones_like(self.face_t), self.face_t)

    def evaluate(self, parameters):
        p, m, w = self.projection, self.projection.mixing, self.weight
        model = EnrichedModalTransport(p, parameters, scalar_mixing=self.mixing,
            thermal_species_ratio=self.ratio, mechanical_work=self.mechanical_work)
        size, count, area = model.size, model.count, model.area
        d = model.local(self.a, self.b)
        ed = model.local(np.ones_like(self.face_t), self.face_t)
        rho, c, h, u, q = [d[k] for k in ("rho", "c", "h", "u", "q")]
        integrand = np.zeros((len(q), 2, count+1))
        integrand[:, 0, 1:] = -d["bm"]
        integrand[:, 1, 0], integrand[:, 1, 1:] = d["force"], -d["bp"]
        f, fe = np.empty_like(integrand), np.empty((len(self.angles), 2, count+1))
        for i, (s, k) in enumerate(zip(self.slices, self.knots)):
            cumulative, endpoint = panel_cumulative(integrand[s], k, self.order)
            f[s] = cumulative/(2*q[s, None, None])
            fe[i] = endpoint/(2*k[-1])
        fm, fp = f[:, 0], f[:, 1]
        rp = fp-u[:, None]*fm
        production = -(2*q*d["uq"])[:, None]*rp
        dot = np.einsum("nik,ni->nk", d["grad_psi"], d["coordinate"])
        diff_c = np.einsum("nik,ni->nk", d["grad_psi"], d["grad_y"])
        diff_h = np.einsum("nik,ni->nk", d["grad_psi"], d["grad_h"])
        rhs = np.zeros((2*size, count+1))
        rhs[:size] = np.einsum("n,nk,nr->kr", w*d["y"], dot, fm, optimize=True)
        rhs[size:] = np.einsum("n,nk,nr->kr", w*d["specific_h"], dot, fm, optimize=True)
        rhs[size:] += np.einsum("n,nk,nr->kr", w, d["psi"], production, optimize=True)
        rhs[:size, 0] -= (w*area*rho*self.chi/(2*p.q0))@diff_c
        rhs[size:, 0] -= (w*area*rho*self.chi*self.ratio/(2*p.q0))@diff_h
        jm = np.vstack([np.einsum("n,nk,nr->kr", w, d["psi"], d[k], optimize=True) for k in ("bc", "bh")])
        theta = m.state[3]
        px, pz = d["bp"]*math.cos(theta), d["bp"]*math.sin(theta)
        px[:, 3] -= area*rho*u*u*math.sin(theta)
        pz[:, 3] += area*rho*u*u*math.cos(theta)
        jf = np.stack([w@v for v in (d["bm"], d["bc"], px, pz, d["bh"]+d["bk"])])
        ue = ed["u"]
        fk = ue[:, None]*fe[:, 1]-.5*ue[:, None]**2*fe[:, 0]
        fh = .5*m.wind**2*fe[:, 0]-fk
        ew = 16*p.q0*self.angular_weights/np.cos(self.angles)**2
        rhs[size:] -= np.einsum("n,nk,nr->kr", ew, ed["psi"], fh, optimize=True)
        force = float(w@(9.81*area*(p.section.rhoa-rho)))
        work = float(w@(d["force"]*u))
        source = IndependentEnergyCrosswind.source_terms(_PhaseForceView(m, force), m.state)
        source[4] += work
        matrix, right = np.vstack([jf, jm-rhs[:, 1:]]), np.r_[source, rhs[:, 0]]
        rates = np.zeros(count)
        rates[5:7] = math.cos(theta), math.sin(theta)
        scales = np.maximum(np.max(abs(matrix[:, model.active]), axis=1), 1e-12)
        scaled = matrix[:, model.active]/scales[:, None]
        singular = np.linalg.svd(scaled, compute_uv=False)
        if singular[-1] <= 1e-12*singular[0]:
            raise ValueError("fixed-mesh coupled rates lost rank")
        rates[model.active] = np.linalg.solve(scaled, (right-matrix@rates)/scales)
        at = np.r_[1., rates]
        em, ep, eh = fe[:, 0]@at, fe[:, 1]@at, fh@at
        ec = ed["y"]*em-area*ed["rho"]*ed["grad_y"][:, 0]*self.edge_chi/(2*p.q0)
        heat = ed["specific_h"]*em-area*ed["rho"]*ed["grad_h"][:, 0]*self.edge_chi*self.ratio/(2*p.q0)-eh
        momentum = ep-m.wind*math.cos(theta)*em
        edge = np.column_stack([ec, heat, momentum])/self.edge_scales
        chi_p = -(rp@at)/(area*rho*d["uq"])
        edge_chi_p = -(ep-ue*em)/(area*ed["rho"]*ed["uq"])
        r2 = (m.sy*m.sy+m.sn*m.sn)*q
        moments = area*w@np.column_stack([rho*u, c*u, rho*u*u, h*u+.5*rho*u**3,
            9.81*(p.section.rhoa-rho), r2*h*u])
        dr, dc, dh = [d[k][:, 7:] for k in ("drho", "dc", "dh")]
        moment_jacobian = area*np.stack([(w*u)@dr, (w*u)@dc, (w*u*u)@dr,
            (w*u)@dh+(.5*w*u**3)@dr, -9.81*w@dr, (w*r2*u)@dh])
        si = self.sample_indices
        stress_penalty = np.minimum(np.r_[chi_p[si]/self.chi[si], edge_chi_p/self.edge_chi], 0.)
        # Equal dimensionless block weighting; these are numerical objectives,
        # not fitted physical diffusivities. No negative stress is clipped in
        # the equations, diagnostic values or acceptance test.
        fit_edge = edge[self.face_slice]
        objective = np.r_[fit_edge.T.ravel()/math.sqrt(len(fit_edge)),
            stress_penalty/math.sqrt(len(stress_penalty))]
        heat_axial = w@d["bh"]@rates
        heat_prod, heat_boundary = w@production@at, ew@fh@at
        return dict(rates=rates, matrix=matrix, right=right, source=source,
            moment_jacobian=moment_jacobian, moments=moments, edge_residual=edge,
            edge_defects=dict(zip(("hydrogen", "heat", "momentum"), np.max(abs(edge), axis=0))),
            objective=objective, minimum_chi_momentum=float(min(np.min(chi_p), np.min(edge_chi_p))),
            maximum_outward_mass=float(max(em)), matrix_condition=float(singular[0]/singular[-1]),
            linear_scaled_error=float(max(abs(matrix@rates-right)/np.maximum(abs(right), 1.))),
            weak_heat_scaled_error=float(abs(heat_axial-heat_prod+heat_boundary)/max(1., abs(heat_axial), abs(heat_prod), abs(heat_boundary))),
            mass_boundary_error=float(ew@em+source[0]),
            curvature_half_width=float(abs(rates[3])*math.sqrt(2*p.q0)*m.sn))


class SimultaneousShapeInitializer:
    """Constrained SQP, recomputing all rates on every trial shape/FD sample."""
    def __init__(self, mesh):
        self.mesh, self.projection = mesh, mesh.projection
        self.constraint = ConservativeEdgeRefit(self.projection)

    def retract(self, parameters, maximum_iterations=8):
        par = np.array(parameters, float)
        p = self.projection
        for _ in range(maximum_iterations):
            out = self.mesh.evaluate(par)
            r = (out["moments"]-p.target)/p.scales
            if max(abs(r)) <= 1e-9:
                return par, out
            step, _ = self.constraint.step(par, r, out["moment_jacobian"]/p.scales[:, None],
                np.zeros(p.count), np.eye(p.count))
            for factor in (1., .5, .25, .125, .0625):
                try:
                    new = self.mesh.evaluate(par+factor*step)
                    if np.linalg.norm((new["moments"]-p.target)/p.scales) < np.linalg.norm(r):
                        par += factor*step
                        break
                except (ValueError, RuntimeError):
                    continue
            else:
                raise ValueError("coupled fit-mesh conservation retraction failed")
        raise ValueError("coupled fit-mesh conservation retraction iteration limit")

    def fit(self, initial, *, maximum_iterations=6, callback=None):
        p, history = self.projection, []
        par, out = self.retract(initial)
        for iteration in range(maximum_iterations):
            residual = (out["moments"]-p.target)/p.scales
            row = dict(iteration=iteration, parameters=par.copy(), edge_defects=out["edge_defects"],
                moment_error=float(max(abs(residual))), minimum_chi_momentum=out["minimum_chi_momentum"],
                objective_norm=float(np.linalg.norm(out["objective"])))
            history.append(row)
            if callback:
                callback(row)
            if (max(out["edge_defects"].values()) <= .04 and out["minimum_chi_momentum"] >= 0.
                    and out["maximum_outward_mass"] < 0. and out["curvature_half_width"] < .1):
                row["stopped"] = "fit_targets_met_pending_independent_verification"
                break
            columns = []
            constraint_values = self.constraint.bounds_matrix@par
            margin = .1-abs(constraint_values)
            for j in range(p.count):
                speed = abs(self.constraint.bounds_matrix[:, j])
                mask = speed > 1e-12
                step = min(1e-6, float(np.min(.2*margin[mask]/speed[mask])))
                if step < 1e-10:
                    raise ValueError("no resolved two-sided shape difference inside trust region")
                e = np.zeros(p.count)
                e[j] = step
                plus, minus = self.mesh.evaluate(par+e), self.mesh.evaluate(par-e)
                columns.append((plus["objective"]-minus["objective"])/(2*step))
            jac = np.column_stack(columns)
            dx, row["qp"] = self.constraint.step(par, residual,
                out["moment_jacobian"]/p.scales[:, None], out["objective"], jac)
            for factor in (1., .5, .25, .125, .0625):
                try:
                    candidate, new = self.retract(par+factor*dx)
                    if np.linalg.norm(new["objective"]) < np.linalg.norm(out["objective"])*(1-1e-5*factor):
                        par, out = candidate, new
                        row["accepted_factor"] = factor
                        break
                except (ValueError, RuntimeError):
                    continue
            else:
                row["stopped"] = "no_acceptable_self_consistent_step"
                break
        return dict(parameters=par, history=history, output=out)
