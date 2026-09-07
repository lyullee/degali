"""Opt-in reduced finite-TKE transport, NOT a validated LH2 closure.

Q=rho*k is separate from thermodynamic enthalpy. All old scalar weak rows
are retained; dissipation heats H, while production first enters Q. The
caller supplies Q, ambient k, positive mixing/time fields and four fixed
circulation amplitudes. These inputs are not inferred or calibrated here.
Normal Reynolds stresses, pressure transport and buoyancy correlations are
still omitted in this reduced mean-momentum approximation.
"""

import math
import numpy as np
from .enriched_transport import EnrichedModalTransport, panel_cumulative
from .aligned_flux_witness import modes
from .phase_radial_quadrature import gauss_rule
from .reservoir_thermal import _PhaseForceView
from .energy_crosswind import IndependentEnergyCrosswind
from .stress_realizability import shear_tke_lower_bound


def _positive_field(supplied, data, name):
    """Fields depend on the current state, never on unknown state rates."""
    raw = supplied(data) if callable(supplied) else supplied
    value = np.asarray(raw, float)
    if value.ndim == 0:
        value = np.full(len(data['q']), float(value))
    if (value.shape != data['q'].shape or not np.all(np.isfinite(value))
            or np.any(value <= 0.)):
        raise ValueError(f'{name} must be a positive finite scalar or matching field')
    return value


def _ratio_field(supplied, data, name):
    raw = supplied(data) if callable(supplied) else supplied
    value = np.asarray(raw, float)
    if value.ndim == 0:
        value = np.full(len(data['q']), float(value))
    if (
        value.shape != data['q'].shape or not np.all(np.isfinite(value))
        or np.any(value < 0.) or np.any(value > 2.)
    ):
        raise ValueError(f'{name} must be a finite scalar or matching field in [0, 2]')
    return value


class FiniteTkeModalTransport:
    def __init__(self, base_model, tke_parameters, *, ambient_tke,
                 tke_diffusivity, dissipation_time, circulation_amplitudes,
                 normal_stress_ratio=None):
        if not isinstance(base_model, EnrichedModalTransport):
            raise TypeError('an explicit enriched mean-state transport model is required')
        self.base = base_model
        self.size, self.mean_count = base_model.size, base_model.count
        self.tke_count = self.size + 1
        self.count = self.mean_count + self.tke_count
        self.parameters = np.array(tke_parameters, float)
        if (self.parameters.shape != (self.tke_count,)
                or not np.all(np.isfinite(self.parameters))):
            raise ValueError('finite log-Q center and one coefficient per square mode required')
        if not np.isscalar(ambient_tke) or not np.isfinite(ambient_tke) or ambient_tke < 0.:
            raise ValueError('explicit nonnegative finite ambient TKE required')
        self.ambient_tke = float(ambient_tke)
        self.diffusivity, self.dissipation_time = tke_diffusivity, dissipation_time
        self.amplitudes = np.array(circulation_amplitudes, float)
        if self.amplitudes.shape != (4,) or not np.all(np.isfinite(self.amplitudes)):
            raise ValueError('four finite fixed circulation amplitudes required; no closure is inferred')
        self.normal_stress_ratio = normal_stress_ratio
        self.active = np.r_[base_model.active, np.arange(self.mean_count, self.count)]
        # Also reject invalid scalar/callable inputs at construction, not only
        # after an expensive quadrature. Full positivity remains pointwise.
        self.local(np.array([0.]), np.array([0.]))

    def _augment(self, data):
        d = dict(data)
        b, p = self.base, self.base.projection
        correction = d['psi'] @ self.parameters[1:]
        if np.any(abs(correction) > .1):
            raise ValueError('Q shape exceeds the declared log-shape trust region')
        exponent = p.section.velocity_shape_exponent
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            q_density = np.exp(self.parameters[0] - exponent*d['q'] + correction)
        if not np.all(np.isfinite(q_density)) or np.any(q_density <= 0.):
            raise ValueError('Q must stay strictly positive and finite; no clipping is used')
        theta = np.column_stack([np.ones(len(q_density)), d['psi']])
        grad_theta = np.concatenate([np.zeros((len(q_density), 2, 1)), d['grad_psi']], axis=2)
        grad_q = q_density[:, None] * (-2*p.q0*exponent*d['coordinate']
            + np.einsum('nik,k->ni', d['grad_psi'], self.parameters[1:]))
        grad_c = d['c'][:, None] * (-2*p.q0*d['coordinate']
            + np.einsum('nik,k->ni', d['grad_psi'], b.parameters[:self.size]))
        grad_h = d['h'][:, None] * (-2*p.q0/p.section.thermal_width_ratio**2*d['coordinate']
            + np.einsum('nik,k->ni', d['grad_psi'], b.parameters[self.size:]))
        grad_rho = (grad_h - d['hc'][:, None]*grad_c)/d['hr'][:, None]
        tke = q_density/d['rho']
        grad_k = (grad_q - tke[:, None]*grad_rho)/d['rho'][:, None]
        bq = np.zeros((len(q_density), self.count))
        bq[:, :self.mean_count] = b.area*q_density[:, None]*d['du']
        bq[:, 2] += b.area*q_density*d['u']
        bq[:, self.mean_count:] = (b.area*q_density*d['u'])[:, None]*theta
        d.update(q_density=q_density, tke=tke, grad_q=grad_q, grad_k=grad_k,
                 theta_q=theta, grad_theta_q=grad_theta, bq=bq)
        d['chi_k'] = _positive_field(self.diffusivity, d, 'TKE diffusivity')
        d['tau'] = _positive_field(self.dissipation_time, d, 'dissipation time')
        if self.normal_stress_ratio is None:
            d['normal_stress_ratio'] = np.zeros(len(q_density))
        else:
            d['normal_stress_ratio'] = _ratio_field(
                self.normal_stress_ratio, d, 'normal stress ratio'
            )
        d['dissipation'] = b.area*q_density/d['tau']
        if not np.all(np.isfinite(d['dissipation'])):
            raise ValueError('dissipation overflow; supplied state/time is inadmissible')
        return d

    def local(self, a, b):
        return self._augment(self.base.local(a, b))

    def _circulation(self, d):
        if not np.any(self.amplitudes):
            return np.zeros((len(d['q']), 2)), np.zeros((len(d['q']), 2))
        v = modes(self.base, d)
        return (np.einsum('nik,k->ni', v['mass_modes'], self.amplitudes),
                np.einsum('nik,k->ni', v['momentum_modes'], self.amplitudes))

    def assemble(self, *, order=8, angular_order=16, split_angles=True,
                 batch_size=8, solve=True, callback=None,
                 include_work=False, include_thermal_relaxation=False,
                 include_phase_transition=False, include_turbulent_heat_exchange=False):
        """Assemble all conservation rows with bounded-memory phase quadrature.

        k/tau/chi inputs are frozen at this state. Passing solve=False returns
        the same affine operator without accepting any candidate trajectory.
        Natural Q boundary data are used here, NOT a claimed gradient match.
        """
        if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
            raise ValueError('positive integer batch size required')
        b, n, old, size = self.base, self.count, self.mean_count, self.size
        p, m = b.projection, b.projection.mixing
        if split_angles:
            ak = b.quadrature.angular_knots(b.parameters)
            crossed = b.mixing.knots[(b.mixing.knots > p.q0) & (b.mixing.knots < 2*p.q0)]
            ak = np.sort(np.r_[ak, np.arccos(np.sqrt(p.q0/crossed))])
            ak = ak[np.r_[True, np.diff(ak) > 2e-13]]
        else:
            ak = np.array([0., math.pi/4])
        ax, aw0 = gauss_rule(angular_order)
        angles = (ak[:-1, None] + np.diff(ak)[:, None]*ax).ravel()
        aws = (np.diff(ak)[:, None]*aw0).ravel()
        x, w = gauss_rule(order)
        jf, jch, jq = np.zeros((5, n)), np.zeros((2*size, n)), np.zeros((size+1, n))
        rch, rq = np.zeros((2*size, n+1)), np.zeros((size+1, n+1))
        names = ('heat_axial', 'mean_kinetic_axial', 'tke_axial', 'production',
                 'heat_boundary', 'mean_kinetic_boundary', 'tke_boundary',
                 'mass_boundary', 'normal_stress_boundary')
        ledger = {key: np.zeros(n+1) for key in names}
        force, work, dissipation, q_flux, peak = 0., 0., 0., 0., 0
        for start in range(0, len(angles), batch_size):
            phi, aw = angles[start:start+batch_size], aws[start:start+batch_size]
            knots, qs, weights, phis, slices = [], [], [], [], []
            offset = 0
            for angle, weight, c in zip(phi, aw, b.quadrature.partitions(b.parameters, phi)):
                crossed = b.mixing.knots[(b.mixing.knots > 0.) & (b.mixing.knots < c[-1])]
                c = np.sort(np.r_[c, crossed])
                c = c[np.r_[True, np.diff(c) > c[-1]*2e-13]]
                q = (c[:-1, None] + np.diff(c)[:, None]*x).ravel()
                knots.append(c); qs.append(q); phis.append(np.full(len(q), angle))
                weights.append((8*weight*np.diff(c)[:, None]*w).ravel())
                slices.append(slice(offset, offset+len(q))); offset += len(q)
            q, weight, ph = np.concatenate(qs), np.concatenate(weights), np.concatenate(phis)
            peak = max(peak, len(q))
            radius = np.sqrt(q/p.q0)
            d = self.local(np.minimum(radius*np.cos(ph), 1.), np.minimum(radius*np.sin(ph), 1.))
            ed = self.local(np.ones(len(phi)), np.tan(phi))
            integrand = np.zeros((len(q), 2, n+1))
            integrand[:, 0, 1:old+1] = -d['bm']
            integrand[:, 1, 0], integrand[:, 1, 1:old+1] = d['force'], -d['bp']
            f, fe = np.empty_like(integrand), np.empty((len(phi), 2, n+1))
            for j, (s, c) in enumerate(zip(slices, knots)):
                cum, end = panel_cumulative(integrand[s], c, order)
                f[s], fe[j] = cum/(2*q[s, None, None]), end/(2*c[-1])
            fm, fp = f[:, 0], f[:, 1]
            cm, cp = self._circulation(d)
            ecm, ecp = self._circulation(ed)
            production = -(2*q*d['uq'])[:, None]*(fp-d['u'][:, None]*fm)
            grad_u = (2*p.q0*d['uq'])[:, None]*d['coordinate']
            production[:, 0] -= np.sum(grad_u*(cp-d['u'][:, None]*cm), axis=1)
            dot = np.einsum('nik,ni->nk', d['grad_psi'], d['coordinate'])
            cdot = np.einsum('nik,ni->nk', d['grad_psi'], cm)
            for rows, scalar, gradient, chi in (
                (slice(0, size), d['y'], d['grad_y'], b.mixing.evaluate(q)[:, 0]),
                (slice(size, 2*size), d['specific_h'], d['grad_h'], b.mixing.evaluate(q)[:, 0]*b.ratio),
            ):
                rch[rows] += np.einsum('n,nk,nr->kr', weight*scalar, dot, fm, optimize=True)
                rch[rows, 0] += (weight*scalar) @ cdot
                diff = np.einsum('nik,ni->nk', d['grad_psi'], gradient)
                rch[rows, 0] -= (weight*b.area*d['rho']*chi/(2*p.q0)) @ diff
            rch[size:, 0] += (weight*d['dissipation']) @ d['psi']
            # Q zeroth and all shape tests: P-D, not immediate P heating.
            dotq = np.column_stack([np.zeros(len(q)), dot])
            cdotq = np.column_stack([np.zeros(len(q)), cdot])
            rq += np.einsum('n,nk,nr->kr', weight*d['tke'], dotq, fm, optimize=True)
            rq[:, 0] += (weight*d['tke']) @ cdotq
            dq = np.einsum('nik,ni->nk', d['grad_theta_q'], d['grad_k'])
            rq[:, 0] -= (weight*b.area*d['rho']*d['chi_k']/(2*p.q0)) @ dq
            rq += np.einsum('n,nk,nr->kr', weight, d['theta_q'], production, optimize=True)
            rq[:, 0] -= (weight*d['dissipation']) @ d['theta_q']
            jch[:size, :old] += np.einsum('n,nk,nr->kr', weight, d['psi'], d['bc'], optimize=True)
            jch[size:, :old] += np.einsum('n,nk,nr->kr', weight, d['psi'], d['bh'], optimize=True)
            jq += np.einsum('n,nk,nr->kr', weight, d['theta_q'], d['bq'], optimize=True)
            theta = m.state[3]
            normal = b.area*d['normal_stress_ratio']*d['q_density']
            normal_modal = np.zeros((len(q), n))
            normal_modal[:, old:] = normal[:, None]*d['theta_q']
            normal_work = normal_modal*d['u'][:, None]
            px = np.zeros((len(q), n))
            pz = np.zeros((len(q), n))
            px[:, :old] = d['bp']*math.cos(theta)
            pz[:, :old] = d['bp']*math.sin(theta)
            px[:, 3] -= b.area*d['rho']*d['u']**2*math.sin(theta)
            pz[:, 3] += b.area*d['rho']*d['u']**2*math.cos(theta)
            px[:, old:] += normal_modal[:, old:]*math.cos(theta)
            pz[:, old:] += normal_modal[:, old:]*math.sin(theta)
            bm = np.zeros((len(q), n))
            bc = np.zeros((len(q), n))
            thermal = np.zeros((len(q), n))
            bm[:, :old] = d['bm']
            bc[:, :old] = d['bc']
            thermal[:, :old] = d['bh'] + d['bk']
            thermal[:, old:] = normal_work[:, old:]
            jf += np.stack([weight @ v for v in (
                bm,
                bc,
                px,
                pz,
                thermal,
            )])
            jf[4] += weight @ d['bq']
            force += float(weight @ (9.81*b.area*(p.section.rhoa-d['rho'])))
            work += float(weight @ (d['force']*d['u']))
            dissipation += float(weight @ d['dissipation'])
            q_flux += float(weight @ (b.area*d['q_density']*d['u']))
            em, ep = fe[:, 0].copy(), fe[:, 1].copy()
            em[:, 0] += ecm[:, 0]; ep[:, 0] += ecp[:, 0]
            en = b.area*ed['normal_stress_ratio']*ed['q_density']
            normal_edge = np.zeros((len(phi), n+1))
            normal_edge[:, old+1:] = en[:, None]*ed['theta_q']
            normal_work_edge = normal_edge*ed['u'][:, None]
            normal_work_lateral = normal_work_edge*math.sin(theta)
            fq = self.ambient_tke*em
            fk = ed['u'][:, None]*ep - .5*ed['u'][:, None]**2*em
            if np.any(d['normal_stress_ratio']):
                fk += normal_work_lateral
            fh = .5*m.wind**2*em - fk
            ew = 16*p.q0*aw/np.cos(phi)**2
            rch[size:] -= np.einsum('n,nk,nr->kr', ew, ed['psi'], fh, optimize=True)
            rq -= np.einsum('n,nk,nr->kr', ew, ed['theta_q'], fq, optimize=True)
            ledger['heat_axial'][1:old+1] += weight @ d['bh']
            ledger['mean_kinetic_axial'][1:old+1] += weight @ d['bk']
            ledger['tke_axial'][1:] += weight @ d['bq']
            ledger['production'] += weight @ production
            for key, value in (('heat_boundary', fh), ('mean_kinetic_boundary', fk),
                               ('tke_boundary', fq), ('mass_boundary', em),
                               ('normal_stress_boundary', normal_work_edge)):
                ledger[key] += ew @ value
            if callback:
                callback(dict(completed_angles=min(start+batch_size, len(angles)),
                              total_angles=len(angles), peak_batch_nodes=peak))
        source = IndependentEnergyCrosswind.source_terms(
            _PhaseForceView(m, force),
            m.state,
            include_work=include_work,
            include_thermal_relaxation=include_thermal_relaxation,
            include_phase_transition=include_phase_transition,
            include_turbulent_heat_exchange=include_turbulent_heat_exchange,
        )
        source[4] += work + self.ambient_tke*source[0]
        thermal_source = ['dissipation_only']
        if include_work:
            thermal_source.append('buoyancy_work')
        if include_thermal_relaxation:
            thermal_source.append('thermal_relaxation')
        if include_phase_transition:
            thermal_source.append('phase_transition')
        if include_turbulent_heat_exchange:
            thermal_source.append('turbulent_heat_exchange')
        matrix = np.vstack([jf, jch-rch[:, 1:], jq-rq[:, 1:]])
        right = np.r_[source, rch[:, 0], rq[:, 0]]
        out = dict(matrix=matrix, right=right, source=source,
                   advective_jacobian=np.vstack([jf, jch, jq]), ledger_coefficients=ledger,
                   dissipation=dissipation, force=force, work=work, tke_flux=q_flux,
                   order=order, angular_order=angular_order, angles=len(angles),
                   peak_batch_nodes=peak, batch_size=batch_size, phase_split_angles=split_angles,
                   adopted=False, field_scored=False, closure_inputs_validated=False,
                   thermal_source='+'.join(thermal_source), tke_source='production_minus_dissipation')
        if solve:
            rates = np.zeros(n)
            rates[5:7] = math.cos(m.state[3]), math.sin(m.state[3])
            scales = np.maximum(np.max(abs(matrix[:, self.active]), axis=1), 1e-12)
            scaled = matrix[:, self.active]/scales[:, None]
            singular = np.linalg.svd(scaled, compute_uv=False)
            if singular[-1] <= 1e-12*singular[0]:
                raise ValueError('finite-TKE coupled matrix is numerically rank deficient')
            rates[self.active] = np.linalg.solve(scaled, (right-matrix@rates)/scales)
            out.update(rates=rates, matrix_condition=float(singular[0]/singular[-1]),
                       linear_scaled_error=float(max(abs(matrix@rates-right)/np.maximum(abs(right), 1.))))
            out['ledgers'] = self.energy_ledgers(out, rates)
        return out

    def energy_ledgers(self, operator, rates):
        rates = np.asarray(rates, float)
        if rates.shape != (self.count,) or not np.all(np.isfinite(rates)):
            raise ValueError('matching finite state rates required')
        values = {key: float(value @ np.r_[1., rates])
                  for key, value in operator['ledger_coefficients'].items()}
        d, w = operator['dissipation'], operator['work']
        values.update(dissipation=d, work=w)
        groups = dict(
            heat=(values['heat_axial'], values['heat_boundary'], -d),
            tke=(values['tke_axial'], values['tke_boundary'], -values['production'], d),
            mean_kinetic=(values['mean_kinetic_axial'], values['mean_kinetic_boundary'], values['production'], -w),
        )
        residual = {key: math.fsum(terms) for key, terms in groups.items()}
        scaled = {key: abs(residual[key])/max(1., *(abs(v) for v in terms)) for key, terms in groups.items()}
        ambient_energy = .5*self.base.projection.mixing.wind**2 + self.ambient_tke
        boundary = math.fsum(values[key] for key in ('heat_boundary', 'mean_kinetic_boundary', 'tke_boundary'))
        return dict(terms=values, residuals=residual, scaled_errors=scaled,
                    total_budget=math.fsum(residual.values()),
                    natural_boundary_identity=boundary-ambient_energy*values['mass_boundary'],
                    mass_boundary_error=values['mass_boundary']+operator['source'][0])

    def independent_diagnostics(self, rates, *, angles=65, order=16):
        """Fresh rays, including endpoints: necessary conditions only, no clipping."""
        rates = np.asarray(rates, float)
        if rates.shape != (self.count,) or not np.all(np.isfinite(rates)):
            raise ValueError('matching finite state rates required')
        if not isinstance(angles, int) or isinstance(angles, bool) or angles < 2:
            raise ValueError('at least two independent rays required')
        b, p, m = self.base, self.base.projection, self.base.projection.mixing
        phi = np.linspace(0., math.pi/4, angles)
        at = np.r_[1., rates[:self.mean_count]]
        edge_max = np.zeros(3)
        q_edge, q_edge_raw = 0., 0.
        minimum_k, margin, minimum_prod, minimum_chi = np.inf, np.inf, np.inf, np.inf
        max_mass, nonradial = -np.inf, 0.
        lengths = math.sqrt(2*p.q0)*np.array([m.sy, m.sn])
        for angle, cuts in zip(phi, b.quadrature.partitions(b.parameters, phi)):
            crossed = b.mixing.knots[(b.mixing.knots > 0.) & (b.mixing.knots < cuts[-1])]
            cuts = np.sort(np.r_[cuts, crossed])
            cuts = cuts[np.r_[True, np.diff(cuts) > cuts[-1]*2e-13]]
            raw, fm, fp, fe, _ = b._ray(angle, cuts, order)
            d = self._augment(raw)
            ed = self.local(np.array([1.]), np.array([math.tan(angle)]))
            for local, mass_coeff, momentum_coeff in ((d, fm, fp), (ed, fe[None, 0], fe[None, 1])):
                cm, cp = self._circulation(local)
                mass = local['coordinate']*(mass_coeff @ at)[:, None] + cm
                momentum = local['coordinate']*(momentum_coeff @ at)[:, None] + cp
                stress = momentum-local['u'][:, None]*mass
                coord = local['coordinate']
                radial = np.sum(coord*stress, axis=1)/np.sum(coord**2, axis=1)
                fraction = np.linalg.norm(stress-coord*radial[:, None], axis=1)/np.maximum(np.linalg.norm(stress, axis=1), 1e-30)
                production = -2*p.q0*local['uq']*np.sum(coord*stress, axis=1)
                chi = -radial/(b.area*local['rho']*local['uq'])
                k_min = shear_tke_lower_bound(stress, b.area, local['rho'], lengths)['minimum_tke']
                minimum_k = min(minimum_k, float(min(local['tke'])))
                margin = min(margin, float(min(local['tke']-k_min)))
                minimum_prod = min(minimum_prod, float(min(production)))
                minimum_chi = min(minimum_chi, float(min(chi)))
                nonradial = max(nonradial, float(max(fraction)))
                if local is ed:
                    em, ep = mass[:, 0], momentum[:, 0]
                    ce, oldfm = b.mixing.evaluate(ed['q']).T
                    old_fields = m.local(ed['q'])
                    scale = np.column_stack([np.maximum(abs(old_fields['y']*oldfm), 1e-12),
                        np.maximum(abs(old_fields['h']/old_fields['rho']*oldfm), 1.),
                        np.maximum(abs(old_fields['u']*oldfm), 1.)])
                    common = b.area*ed['rho']/(2*p.q0)
                    fh = .5*(m.wind**2+ed['u']**2)*em-ed['u']*ep
                    defects = np.column_stack([ed['y']*em-common*ce*ed['grad_y'][:, 0],
                        ed['specific_h']*em-common*ce*b.ratio*ed['grad_h'][:, 0]-fh,
                        ep-m.wind*math.cos(m.state[3])*em])
                    edge_max = np.maximum(edge_max, np.max(abs(defects)/scale, axis=0))
                    q_constitutive = ed['tke']*em-common*ed['chi_k']*ed['grad_k'][:, 0]
                    q_natural = self.ambient_tke*em
                    q_error = abs(q_constitutive-q_natural)
                    q_scale = np.maximum.reduce([abs(ed['tke']*em), abs(q_natural), np.full(len(em), 1e-12)])
                    q_edge = max(q_edge, float(max(q_error/q_scale)))
                    q_edge_raw = max(q_edge_raw, float(max(q_error)))
                    max_mass = max(max_mass, float(max(em)))
        return dict(edge_defects=dict(zip(('hydrogen', 'heat', 'momentum'), edge_max)),
                    tke_edge_scaled_error=q_edge, tke_edge_maximum_raw_error=q_edge_raw,
                    minimum_tke=minimum_k, minimum_tke_above_shear_bound=margin,
                    minimum_shear_production=minimum_prod, minimum_momentum_diffusivity=minimum_chi,
                    maximum_outward_mass=max_mass, maximum_nonradial_stress_fraction=nonradial,
                    curvature_half_width=float(abs(rates[3])*lengths[1]), angles=angles, order=order,
                    physical_closure_passed=False, field_scored=False)
