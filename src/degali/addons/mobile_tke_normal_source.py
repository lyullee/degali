"""Opt-in source matching with mobile centers/area and prescribed Q/Rss.

Only six source matching degrees of freedom are solved. No observations,
turbulence inputs, transport rates or pressure compensation are fitted.
The inherited mixing object is used solely to construct a local field helper;
it is never evaluated or certified for the new source state.
"""

import math
import numpy as np
from .enriched_segments import encode_enriched, decode_enriched
from .enriched_transport import EnrichedModalTransport
from .exact_transverse_geometry import exact_geometry_field_view
from .edge_conservative_refit import FaceSplitSquareMoments


class MobileTkeNormalSource:
    variable_names = ('log_rho_center', 'log_C_center', 'log_area',
                      'log_excess_velocity', 'C_log_q_shape', 'H_log_q_shape')

    def __init__(self, base_model, tke_parameters, *, axial_variance_fraction,
                 pressure_assumption):
        if not isinstance(base_model, EnrichedModalTransport):
            raise TypeError('an explicit enriched mean field helper is required')
        ratio = float(axial_variance_fraction)
        if not math.isfinite(ratio) or not 0. <= ratio <= 2.:
            raise ValueError('explicit Rss/k must lie in [0,2]')
        if pressure_assumption != 'ambient_pressure_no_compensation':
            raise ValueError('explicit ambient pressure with no compensation is required')
        p = base_model.projection
        qpar = np.array(tke_parameters, float)
        if qpar.shape != (base_model.size+1,) or not np.all(np.isfinite(qpar)):
            raise ValueError('explicit finite log-Q center and square modes required')
        self.base, self.q_parameters = base_model, qpar
        self.ratio, self.pressure_assumption = ratio, pressure_assumption
        self.origin = encode_enriched(p, base_model.parameters)
        self.mapping = np.zeros((base_model.count, 6))
        self.mapping[[0, 1, 2, 4], np.arange(4)] = 1.
        self.mapping[7:7+base_model.size, 4] = base_model.q_mode
        self.mapping[7+base_model.size:, 5] = base_model.q_mode
        self.lower, self.upper = np.full(6, -.1), np.full(6, .1)
        t = np.linspace(0., 1., 33)
        a, b = np.meshgrid(t, t, indexing='ij')
        prep = p.prepare(a.ravel(), b.ravel())
        if max(abs(prep['psi'] @ base_model.q_mode-prep['q'])) > 1e-11:
            raise ValueError('q direction is not exactly represented in the scalar basis')
        self._q_density(p, prep)
        positive = prep['q'] > 0.
        for column, start in ((4, 0), (5, base_model.size)):
            initial = prep['psi'] @ base_model.parameters[start:start+base_model.size]
            if max(abs(initial)) > .0999:
                raise ValueError('initial scalar shape is outside the sampled trust region')
            self.lower[column] = np.max((-.0999-initial[positive])/prep['q'][positive])
            self.upper[column] = np.min((.0999-initial[positive])/prep['q'][positive])

    def _q_density(self, projection, prepared):
        correction = prepared['psi'] @ self.q_parameters[1:]
        if np.any(abs(correction) > .1):
            raise ValueError('prescribed Q shape leaves the log-shape trust region')
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            value = np.exp(self.q_parameters[0]
                           -projection.section.velocity_shape_exponent*prepared['q']+correction)
        if np.any(value <= 0.) or not np.all(np.isfinite(value)) or np.any(prepared['u'] <= 0.):
            raise ValueError('positive finite Q and forward source velocity required')
        return value

    def field_view(self, change):
        change = np.asarray(change, float)
        if change.shape != (6,) or not np.all(np.isfinite(change)):
            raise ValueError('six finite mobile-source changes required')
        if np.any(change < self.lower) or np.any(change > self.upper):
            raise ValueError('mobile source leaves the registered trust region')
        encoded = self.origin+self.mapping @ change
        p, par, geometry = exact_geometry_field_view(self.base.projection, encoded)
        local_model = EnrichedModalTransport(p, par, scalar_mixing=self.base.mixing,
            thermal_species_ratio=self.base.ratio,
            mechanical_work='reduced_buoyancy_work_immediate_shear_heat')
        return local_model, encoded, geometry

    def evaluate(self, change, *, order=4, angular_order=8, jacobian=False):
        """Actual six moments, optionally their analytic moving-state Jacobian.

        B and thermal Hsecond are matching conditions, not conserved fluxes.
        The physical Hsecond weight changes with both exact section widths.
        Q is held as a prescribed volumetric field in normalized coordinates.
        """
        model, encoded, geometry = self.field_view(change)
        p, mixing = model.projection, model.projection.mixing
        values, derivative, energy = np.zeros(6), np.zeros((6, 6)), np.zeros(4)
        normal_momentum = 0.
        for prep in FaceSplitSquareMoments(p).rule(model.parameters, order=order, angular_order=angular_order):
            d = model.local(prep['a'], prep['b']) if jacobian else p.fields(model.parameters, prep)
            rho, c, h, u = d['rho'], d['c'], d['h'], prep['u']
            q, w, r2 = self._q_density(p, prep), prep['weights'], prep['r2']
            thermal, kinetic, turbulent = w @ (h*u), .5*w @ (rho*u**3), w @ (q*u)
            energy += [thermal, kinetic, turbulent, self.ratio*turbulent]
            normal_momentum += float(self.ratio*w @ q)
            values += [w @ (rho*u), w @ (c*u), w @ (rho*u*u+self.ratio*q),
                       thermal+kinetic+(1+self.ratio)*turbulent,
                       9.81*w @ (p.section.rhoa-rho), w @ (r2*h*u)]
            if jacobian:
                area = mixing.area
                bq = area*q[:, None]*d['du']
                bq[:, 2] += area*q*u
                bn = np.zeros_like(bq)
                bn[:, 2] = area*self.ratio*q
                buoyancy = -9.81*area*d['drho']
                buoyancy[:, 2] += 9.81*area*(p.section.rhoa-rho)
                dr2 = np.zeros_like(bq)
                dr2[:, :7] = 2*prep['q'][:, None]*(mixing.sy**2*mixing.log_sy_partials[:7]
                                                       +mixing.sn**2*mixing.log_sn_partials[:7])
                second = r2[:, None]*d['bh']+(area*h*u)[:, None]*dr2
                full = np.stack([(w/area) @ v for v in
                    (d['bm'], d['bc'], d['bp']+bn, d['bh']+d['bk']+(1+self.ratio)*bq,
                     buoyancy, second)])
                derivative += full @ self.mapping
        if not np.all(np.isfinite(values)) or (jacobian and not np.all(np.isfinite(derivative))):
            raise ValueError('nonfinite source moments or Jacobian')
        return dict(moments=values, jacobian=derivative if jacobian else None,
            energy_terms=dict(zip(('thermal', 'mean_kinetic', 'tke', 'normal_stress_work'), energy)),
            axial_normal_momentum=normal_momentum, encoded_state=encoded,
            parameters=model.parameters.copy(), geometry=geometry)

    def solve(self, *, total_moment_target, callback=None, maximum_iterations=20):
        target = np.array(total_moment_target, float)
        if target.shape != (6,) or not np.all(np.isfinite(target)) or np.any(target[:3] <= 0.):
            raise ValueError('six finite targets with positive mass/H2/axial momentum required')
        if isinstance(maximum_iterations, bool) or int(maximum_iterations) != maximum_iterations or not 1 <= maximum_iterations <= 20:
            raise ValueError('registered iteration limit is an integer from 1 to 20')
        scales = np.maximum(abs(target), [1e-12, 1e-12, 1., 1., 1e-3, 1.])
        change, history, updates = np.zeros(6), [], 0
        for order, angular in ((4, 8), (8, 16)):
            current = self.evaluate(change, order=order, angular_order=angular, jacobian=True)
            while True:
                residual = (current['moments']-target)/scales
                row = dict(updates=updates, order=order, angular_order=angular,
                           maximum_error=float(max(abs(residual))), changes=change.copy())
                history.append(row)
                if callback:
                    callback(row)
                if max(abs(residual)) <= 1e-10:
                    break
                if updates >= maximum_iterations:
                    raise ValueError('mobile source exceeded registered iteration limit')
                jac = current['jacobian']/scales[:, None]
                column_scale = np.max(abs(jac), axis=0)
                if np.any(column_scale <= 0.):
                    raise ValueError('mobile source has a zero matching direction')
                scaled = jac/column_scale
                singular = np.linalg.svd(scaled, compute_uv=False)
                if singular[-1] <= 1e-12*singular[0]:
                    raise ValueError('mobile source matching Jacobian is rank deficient')
                row['scaled_condition'] = float(singular[0]/singular[-1])
                step = np.linalg.solve(scaled, -residual)/column_scale
                rejected = []
                for factor in (1., .5, .25, .125, .0625, .03125, .015625):
                    candidate = change+factor*step
                    try:
                        trial = self.evaluate(candidate, order=order, angular_order=angular, jacobian=True)
                        norm = np.linalg.norm((trial['moments']-target)/scales)
                        if norm < np.linalg.norm(residual):
                            change, current = candidate, trial
                            row['accepted_factor'] = factor
                            break
                        rejected.append(dict(factor=factor, reason='no residual decrease'))
                    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                        rejected.append(dict(factor=factor, reason=str(exc)))
                else:
                    if callback:
                        callback(dict(failure='no acceptable mobile-source Newton step', rejected=rejected))
                    raise ValueError('no acceptable mobile-source Newton step within registered bounds')
                row['rejected'] = rejected
                updates += 1
        fine = current
        coarse = self.evaluate(change, order=4, angular_order=8)
        errors = abs(fine['moments']-target)/scales
        refinement = abs(fine['moments']-coarse['moments'])/scales
        state, _ = decode_enriched(fine['encoded_state'], self.base.projection.count)
        return dict(**fine, changes=change, physical_state=state, variable_names=self.variable_names,
            tke_parameters=self.q_parameters.copy(), axial_variance_fraction=self.ratio,
            pressure_assumption=self.pressure_assumption, total_moment_target=target,
            moment_errors=errors, moment_refinement=refinement, history=history,
            numerical_passed=bool(max(errors) <= 1e-8 and max(refinement) <= 1e-5),
            physical_initialization_passed=False, normal_transport_closed=False,
            scalar_mixing_recomputed=False, q_was_fitted=False, variance_ratio_was_fitted=False,
            centers_geometry_velocity_fixed=False, old_rates_reused=False, adopted=False, field_scored=False)
