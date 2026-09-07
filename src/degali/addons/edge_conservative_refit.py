"""Face/ray phase-cell quadrature and constrained fixed-transport refitting.

These are research cross-sections, not a transport law for the shape modes.
The earlier failed experiment and its numerical implementation are preserved.
"""

import math
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from .phase_radial_quadrature import gauss_rule
from .edge_phase_quadrature import PhaseSplitSquareMoments


class FaceSplitSquareMoments(PhaseSplitSquareMoments):
    def face_coordinates(self, parameters, angles):
        q = self.projection.q0/np.cos(angles)**2
        return self.coordinates(parameters, q, angles)

    def angular_knots(self, parameters, *, probes=1025):
        if isinstance(probes, bool) or int(probes) != probes or probes < 33:
            raise ValueError("face screening needs an integer >=33")
        angles = np.linspace(0., math.pi/4, int(probes))
        ti, y = self.face_coordinates(parameters, angles)
        inv = self.projection.section.phase_inverse
        lows, highs, targets, types = [], [], [], []
        for values, grid, is_t in ((ti, inv.t_grid, True), (y, inv.y_grid, False)):
            for target in grid[(grid >= min(values)) & (grid <= max(values))]:
                r = values-target
                crossing = np.flatnonzero((r[:-1] < 0.) != (r[1:] < 0.))
                lows.extend(angles[crossing])
                highs.extend(angles[crossing+1])
                targets.extend([target]*len(crossing))
                types.extend([is_t]*len(crossing))
        lo, hi, targets, types = np.asarray(lows), np.asarray(highs), np.asarray(targets), np.asarray(types, bool)
        if len(lo):
            tl, yl = self.face_coordinates(parameters, lo)
            negative_lo = np.where(types, tl, yl) < targets
            for _ in range(36):
                mid = .5*(lo+hi)
                tm, ym = self.face_coordinates(parameters, mid)
                same_side = (np.where(types, tm, ym) < targets) == negative_lo
                lo, hi = np.where(same_side, mid, lo), np.where(same_side, hi, mid)
        knots = np.sort(np.r_[0., .5*(lo+hi), math.pi/4])
        knots = knots[np.r_[True, np.diff(knots) > 2e-13]]
        knots[-1] = math.pi/4
        return knots

    def rule(self, parameters, *, order=8, angular_order=8, probes=1025):
        x, w = gauss_rule(order)
        ax, aw = gauss_rule(angular_order)
        ak = self.angular_knots(parameters, probes=probes)
        angles = (ak[:-1, None]+np.diff(ak)[:, None]*ax).ravel()
        angular_weights = (np.diff(ak)[:, None]*aw).ravel()
        for first in range(0, len(angles), 8):
            subset = angles[first:first+8]
            knots = self.partitions(parameters, subset)
            qs, phis, weights = [], [], []
            for j, k in enumerate(knots):
                nodes = k[:-1, None]+np.diff(k)[:, None]*x
                qs.append(nodes.ravel())
                phis.append(np.full(nodes.size, subset[j]))
                weights.append((np.diff(k)[:, None]*w).ravel()*angular_weights[first+j])
            q, phi = np.concatenate(qs), np.concatenate(phis)
            prepared = self.prepare_ray(q, phi)
            prepared["weights"] = 8*self.projection.mixing.area*np.concatenate(weights)
            prepared["r2"] = (self.projection.mixing.sy**2+self.projection.mixing.sn**2)*q
            yield prepared

    def moments(self, parameters, *, order=8, angular_order=8, probes=1025, jacobian=False):
        p = self.projection
        result = np.zeros(6)
        derivative = np.zeros((6, p.count))
        for prepared in self.rule(parameters, order=order, angular_order=angular_order, probes=probes):
            d = p.fields(parameters, prepared)
            rho, c, h, u = d["rho"], d["c"], d["h"], prepared["u"]
            w, r2 = prepared["weights"], prepared["r2"]
            result += np.array([w@(rho*u), w@(c*u), w@(rho*u*u), w@(h*u+.5*rho*u**3),
                9.81*w@(p.section.rhoa-rho), w@(r2*h*u)])
            if jacobian:
                psi = prepared["psi"]
                zero = np.zeros_like(psi)
                dc = np.column_stack([c[:, None]*psi, zero])
                dh = np.column_stack([zero, h[:, None]*psi])
                dr = (dh-d["hc"][:, None]*dc)/d["hr"][:, None]
                derivative += np.stack([(w*u)@dr, (w*u)@dc, (w*u*u)@dr,
                    (w*u)@dh+(.5*w*u**3)@dr, -9.81*w@dr, (w*r2*u)@dh])
        return (result, derivative) if jacobian else result


class ConservativeEdgeRefit:
    def __init__(self, projection):
        self.projection = projection
        self.quadrature = FaceSplitSquareMoments(projection)
        t = np.linspace(0., 1., 33)
        a, b = np.meshgrid(t, t, indexing="ij")
        psi = projection.prepare(a.ravel(), b.ravel())["psi"]
        zero = np.zeros_like(psi)
        self.bounds_matrix = np.block([[psi, zero], [zero, psi]])

    def step(self, parameters, residual, jacobian, edge, edge_jacobian, *, radius=.02):
        """QP with hard linearized six moments and sampled scalar bounds."""
        p = self.projection
        je = edge_jacobian
        hessian = je.T@je+1e-4*np.eye(p.count)
        gradient = je.T@edge
        constraints = [LinearConstraint(jacobian, -residual, -residual),
            LinearConstraint(self.bounds_matrix, -.0999-self.bounds_matrix@parameters,
                             .0999-self.bounds_matrix@parameters)]
        result = minimize(lambda x: .5*x@hessian@x+gradient@x, np.zeros(p.count),
            jac=lambda x: hessian@x+gradient, method="SLSQP", constraints=constraints,
            bounds=[(-radius, radius)]*p.count, options=dict(ftol=1e-12, maxiter=150))
        dx = result.x
        violation = max(float(max(abs(jacobian@dx+residual))),
                        float(max(abs(self.bounds_matrix@(parameters+dx)))-.0999))
        if violation > 1e-8 or not np.all(np.isfinite(dx)):
            raise ValueError(f"constrained step infeasible: {result.message}, violation={violation:.3e}")
        return dx, dict(success=bool(result.success), message=str(result.message), violation=violation)

    def retract(self, parameters, *, callback=None, maximum_iterations=8):
        p = self.projection
        par = np.array(parameters, float)
        history = []
        for iteration in range(maximum_iterations):
            values, j = self.quadrature.moments(par, jacobian=True)
            r, j = (values-p.target)/p.scales, j/p.scales[:, None]
            history.append(dict(iteration=iteration, maximum_error=float(max(abs(r)))))
            if callback:
                callback(history[-1])
            if max(abs(r)) <= 1e-9:
                return par, history
            step, info = self.step(par, r, j, np.zeros(p.count), np.eye(p.count))
            history[-1]["step"] = info
            for factor in (1., .5, .25, .125, .0625):
                candidate = par+factor*step
                try:
                    new = self.quadrature.moments(candidate)
                    p.edge_residual(candidate)
                    if np.linalg.norm((new-p.target)/p.scales) < np.linalg.norm(r):
                        par = candidate
                        break
                except (ValueError, RuntimeError):
                    continue
            else:
                raise ValueError("phase-accurate conservation retraction did not improve")
        raise ValueError("phase-accurate conservation retraction exceeded iteration limit")

    def fit(self, initial, *, callback=None, maximum_iterations=20):
        p = self.projection
        parameters, initial_retraction = self.retract(initial, callback=callback)
        history = []
        for iteration in range(maximum_iterations):
            values, j = self.quadrature.moments(parameters, jacobian=True)
            r, j = (values-p.target)/p.scales, j/p.scales[:, None]
            edge = p.edge_residual(parameters)
            row = dict(iteration=iteration, moment_error=float(max(abs(r))), edge_max=float(max(abs(edge))),
                       edge_rms=float(np.sqrt(np.mean(edge*edge))))
            history.append(row)
            if callback:
                callback(row)
            if max(abs(edge)) <= .04:
                row["stopped"] = "edge_target_reached_pending_independent_verification"
                break
            step, row["qp"] = self.step(parameters, r, j, edge, p.edge_jacobian(parameters))
            for factor in (1., .5, .25, .125, .0625):
                try:
                    candidate, _ = self.retract(parameters+factor*step)
                    candidate_edge = p.edge_residual(candidate)
                    if np.linalg.norm(candidate_edge) < np.linalg.norm(edge)*(1-1e-5*factor):
                        parameters = candidate
                        break
                except (ValueError, RuntimeError):
                    continue
            else:
                row["stopped"] = "no_acceptable_step"
                break
        return dict(parameters=parameters, initial_retraction=initial_retraction, history=history)
