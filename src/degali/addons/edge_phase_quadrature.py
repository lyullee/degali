"""Phase-cell-split polar integration of square-symmetric enriched scalars.

This is a fixed-section integration diagnostic, not new physics or transport.
Radial monotonicity is screened separately for each angle and new profile.
"""

import math
import numpy as np
from .phase_radial_quadrature import gauss_rule


class PhaseSplitSquareMoments:
    def __init__(self, projection):
        self.projection = projection

    def prepare_ray(self, q, phi):
        radius = np.sqrt(q/self.projection.q0)
        # Roundoff at the known geometric face can exceed one by one ulp.
        a = np.minimum(radius*np.cos(phi), 1.)
        b = np.minimum(radius*np.sin(phi), 1.)
        return self.projection.prepare(a, b)

    def coordinates(self, parameters, q, phi):
        d = self.projection.fields(parameters, self.prepare_ray(q, phi))
        inv = self.projection.section.phase_inverse
        return inv.a/(d["rho"]+inv.k*d["c"]), d["y"]

    def partitions(self, parameters, angles):
        qmax = self.projection.q0/np.cos(angles)**2
        fraction = np.linspace(0., 1., 513)
        q = qmax[:, None]*fraction
        phi = np.broadcast_to(angles[:, None], q.shape)
        ti, y = [v.reshape(q.shape) for v in self.coordinates(parameters, q.ravel(), phi.ravel())]
        if np.any(np.diff(ti, axis=1) < -1e-9) or np.any(np.diff(y, axis=1) > 1e-11):
            raise ValueError("enriched phase ray is not on the supported monotonic branch")
        inv = self.projection.section.phase_inverse
        targets, tags, indices = [], [], []
        for i in range(len(angles)):
            ts = inv.t_grid[(inv.t_grid > ti[i, 0]) & (inv.t_grid < ti[i, -1])]
            ys = inv.y_grid[(inv.y_grid < y[i, 0]) & (inv.y_grid > y[i, -1])]
            targets.extend(np.r_[ts, ys])
            tags.extend(np.r_[np.ones(len(ts), bool), np.zeros(len(ys), bool)])
            indices.extend([i]*(len(ts)+len(ys)))
        targets, tags, indices = np.asarray(targets), np.asarray(tags, bool), np.asarray(indices, int)
        lo, hi = np.zeros(len(targets)), qmax[indices].copy()
        for _ in range(43):
            if not len(targets):
                break
            mid = .5*(lo+hi)
            tm, ym = self.coordinates(parameters, mid, angles[indices])
            below = np.where(tags, tm < targets, ym > targets)
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        roots = .5*(lo+hi)
        partitions = []
        for i in range(len(angles)):
            knots = np.sort(np.r_[0., roots[indices == i], qmax[i]])
            knots = knots[np.r_[True, np.diff(knots) > qmax[i]*2e-13]]
            knots[-1] = qmax[i]
            partitions.append(knots)
        return partitions

    def moments(self, parameters, *, order=16, angular_order=64):
        x, w = gauss_rule(order)
        ax, aw = gauss_rule(angular_order)
        angles, angular_weights = math.pi/4*ax, math.pi/4*aw
        p, m = self.projection, self.projection.mixing
        result = np.zeros(6)
        for first in range(0, len(angles), 8):
            subset = angles[first:first+8]
            knots = self.partitions(parameters, subset)
            qs, phis, weights = [], [], []
            for j, k in enumerate(knots):
                nodes = k[:-1, None]+np.diff(k)[:, None]*x
                qw = np.diff(k)[:, None]*w
                qs.append(nodes.ravel())
                phis.append(np.full(nodes.size, subset[j]))
                weights.append(qw.ravel()*angular_weights[first+j])
            q, phi, weight = np.concatenate(qs), np.concatenate(phis), np.concatenate(weights)
            prepared = self.prepare_ray(q, phi)
            d = p.fields(parameters, prepared)
            rho, c, h, u = d["rho"], d["c"], d["h"], prepared["u"]
            values = np.column_stack([rho*u, c*u, rho*u*u, h*u+.5*rho*u**3,
                9.81*(p.section.rhoa-rho), (m.sy**2+m.sn**2)*q*h*u])
            result += 8*m.area*weight@values
        return result
