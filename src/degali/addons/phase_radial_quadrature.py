"""Vectorized, phase-cell-split radial quadrature for a fixed section.

Only monotonically increasing ideal temperature and decreasing H2 fraction
are supported. Splitting table knots improves numerical integration, not the
EOS or physical closure. Order refinement is still required by callers.
"""

from functools import lru_cache
import math

import numpy as np
from scipy.special import roots_legendre


@lru_cache(maxsize=12)
def gauss_rule(order):
    if isinstance(order, bool) or int(order) != order or order < 4:
        raise ValueError("radial quadrature order must be an integer >=4")
    x, w = roots_legendre(int(order))
    return .5*(x+1.), .5*w


class PhaseRadialQuadrature:
    def __init__(self, section, state, *, probes=513):
        if section.ground_interaction != "free":
            raise NotImplementedError("phase radial quadrature requires a free section")
        section._physical(state)
        self.section, self.state = section, np.array(state, float)
        self.beta = section.thermal_width_ratio
        self.q0 = math.pi*section.k.delta**2/8.
        self.qmax = 2.*self.q0
        if isinstance(probes, bool) or int(probes) != probes or probes < 33:
            raise ValueError("monotonicity screen requires an integer >=33")
        self.probes = int(probes)
        q = np.linspace(0., self.qmax, self.probes)
        ti, y = self.coordinates(q)
        if np.any(np.diff(ti) < -1e-9) or np.any(np.diff(y) > 1e-11):
            raise ValueError("phase radial path is not on the supported monotonic branch")
        inv = section.phase_inverse
        t_targets = inv.t_grid[(inv.t_grid > ti[0]) & (inv.t_grid < ti[-1])]
        y_targets = inv.y_grid[(inv.y_grid < y[0]) & (inv.y_grid > y[-1])]
        targets = np.r_[t_targets, y_targets]
        is_t = np.arange(len(targets)) < len(t_targets)
        lo, hi = np.zeros(len(targets)), np.full(len(targets), self.qmax)
        # Simultaneous bisection across every crossed table edge.
        for _ in range(43):
            if not len(targets):
                break
            mid = .5*(lo+hi)
            tm, ym = self.coordinates(mid)
            below = np.where(is_t, tm < targets, ym > targets)
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        knots = np.sort(np.r_[0., self.q0, .5*(lo+hi), self.qmax])
        # Only remove arithmetic duplicate boundaries; no grid coarsening.
        self.knots = knots[np.r_[True, np.diff(knots) > self.qmax*2e-13]]
        if self.knots[-1] != self.qmax:
            self.knots[-1] = self.qmax
        self._rules = {}

    def check(self):
        if self.section.thermal_width_ratio != self.beta or self.section.ground_interaction != "free":
            raise RuntimeError("fixed phase partition cannot be reused after section changes")

    def coordinates(self, q):
        self.check()
        rho, y, _, _ = self.section.thermodynamic_profile(self.state, np.exp(-np.asarray(q)))
        c = rho*y
        return self.section.phase_inverse.a/(rho+self.section.phase_inverse.k*c), y

    def nodes(self, order=16, *, square=False):
        self.check()
        key = (order, square)
        if key in self._rules:
            return self._rules[key]
        x, w = gauss_rule(order)
        left, right = self.knots[:-1], self.knots[1:]
        nodes = left[:, None]+(right-left)[:, None]*x
        weights = (right-left)[:, None]*np.broadcast_to(w, nodes.shape)
        if square:
            inner = right <= self.q0
            weights[inner] *= 2.*math.pi
            outer = ~inner
            tl = np.sqrt(np.maximum(left[outer]/self.q0-1., 0.))
            tr = np.sqrt(right[outer]/self.q0-1.)
            t = tl[:, None]+(tr-tl)[:, None]*x
            nodes[outer] = self.q0*(1.+t*t)
            weights[outer] = ((tr-tl)[:, None]*w*2.*self.q0*t*(2.*math.pi-8.*np.arctan(t)))
        self._rules[key] = (nodes, weights)
        return nodes, weights

    def integrate(self, function, order=16, *, square=False):
        q, w = self.nodes(order, square=square)
        values = np.asarray(function(q.ravel()), float)
        return np.tensordot(w.ravel(), values, axes=(0, 0))

    def edge_integrate(self, function, order=16):
        """Integral_0^1 function(q0*(1+t²)) dt, split at phase cells."""
        self.check()
        bounds = np.r_[self.q0, self.knots[self.knots > self.q0]]
        t = np.sqrt(np.maximum(bounds/self.q0-1., 0.))
        x, w = gauss_rule(order)
        nodes = t[:-1, None]+np.diff(t)[:, None]*x
        weights = np.diff(t)[:, None]*w
        values = function(self.q0*(1.+nodes.ravel()**2))
        return np.tensordot(weights.ravel(), values, axes=(0, 0))

    def cumulative(self, function, q, order=16):
        """Integral_0^q of vector-valued function with exact panel boundaries."""
        self.check()
        q = np.asarray(q, float)
        if q.ndim != 1 or not np.all(np.isfinite(q)) or np.any(q < 0.) or np.any(q > self.qmax):
            raise ValueError("cumulative radii must be a finite 1-D vector in the section")
        nodes, weights = self.nodes(order)
        values = np.asarray(function(nodes.ravel()), float)
        trailing = values.shape[1:]
        panels = np.sum(values.reshape(nodes.shape+trailing)*weights.reshape(weights.shape+(1,)*len(trailing)), axis=1)
        prefix = np.concatenate([np.zeros((1,)+trailing), np.cumsum(panels, axis=0)])
        index = np.clip(np.searchsorted(self.knots, q, side="right")-1, 0, len(self.knots)-2)
        left = self.knots[index]
        x, w = gauss_rule(order)
        partial_nodes = left[:, None]+(q-left)[:, None]*x
        partial = np.asarray(function(partial_nodes.ravel())).reshape(partial_nodes.shape+trailing)
        partial_weight = ((q-left)[:, None]*w).reshape(partial_nodes.shape+(1,)*len(trailing))
        return prefix[index]+np.sum(partial*partial_weight, axis=1)
