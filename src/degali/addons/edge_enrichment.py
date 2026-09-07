"""Conservative fixed-transport, two-dimensional edge-profile projection.

No transport law for the new coefficients is assumed. This representation
must not be passed to the existing radial-only plume/phase quadratures.
"""

import math
import numpy as np
from scipy.linalg import null_space
from scipy.special import roots_legendre
from numpy.polynomial.legendre import legval, legder


class SquareEvenBasis:
    """Square-symmetric, centre-zero even polynomials; L2 mean norm is one."""
    def __init__(self, degree):
        if isinstance(degree, bool) or int(degree) != degree or degree < 2:
            raise ValueError("supply an integer square-basis degree >=2")
        self.pairs = [(i, j) for i in range(degree+1) for j in range(i, degree+1)
                      if 0 < i+j <= degree]
        self.degree, self.size = int(degree), len(self.pairs)
        v = self._raw(np.array([0.]), np.array([0.]))[0][0]
        self.center = v
        norm2 = v@v
        self.whitener = np.eye(self.size)+((1/math.sqrt(1+norm2)-1)/norm2)*np.outer(v, v)

    def _raw(self, a, b):
        ap, bp, ad, bd = [], [], [], []
        for i in range(self.degree+1):
            coef = np.zeros(2*i+1)
            coef[-1] = math.sqrt(4*i+1)
            ap.append(legval(a, coef))
            bp.append(legval(b, coef))
            ad.append(legval(a, legder(coef)))
            bd.append(legval(b, legder(coef)))
        values, da, db = [], [], []
        for i, j in self.pairs:
            if i == j:
                values.append(ap[i]*bp[j])
                da.append(ad[i]*bp[j])
                db.append(ap[i]*bd[j])
            else:
                values.append((ap[i]*bp[j]+ap[j]*bp[i])/math.sqrt(2))
                da.append((ad[i]*bp[j]+ad[j]*bp[i])/math.sqrt(2))
                db.append((ap[i]*bd[j]+ap[j]*bd[i])/math.sqrt(2))
        return tuple(np.column_stack(x) for x in (values, da, db))

    def values(self, a, b):
        a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
        if a.ndim != 1 or not np.all(np.isfinite(a+b)) or np.any(abs(a) > 1.) or np.any(abs(b) > 1.):
            raise ValueError("basis coordinates must be finite vectors inside the square")
        v, da, db = self._raw(a, b)
        return (v-self.center)@self.whitener, da@self.whitener, db@self.whitener


class FixedTransportEdgeProjection:
    """Change scalar shapes at fixed geometry/velocity/transport, preserving moments."""
    def __init__(self, reservoir, baseline, *, degree):
        if not baseline["valid"]:
            raise ValueError("edge projection requires a passed weak baseline")
        self.reservoir, self.baseline = reservoir, baseline
        self.mixing = reservoir.mixing
        self.section = self.mixing.section
        self.basis = SquareEvenBasis(degree)
        self.count = 2*self.basis.size
        self.q0 = self.mixing.partition.q0
        self.hc = self.section.phase_inverse.enthalpy_and_slope(
            self.mixing.state[0], self.mixing.state[0]*self.mixing.state[1])[0]
        if self.hc == 0.:
            raise ValueError("zero excess enthalpy needs a different shape coordinate")
        self._grids, self._edges = {}, {}
        m = self.mixing
        def target(q):
            d = m.local(q)
            rho, c, h, u = [d[k] for k in ("rho", "c", "h", "u")]
            return m.area*np.column_stack([rho*u, c*u, rho*u*u, h*u+.5*rho*u**3,
                9.81*(self.section.rhoa-rho), (m.sy*m.sy+m.sn*m.sn)*q*h*u])
        self.target = m.partition.integrate(target, 16, square=True)
        self.scales = np.maximum(abs(self.target), [1e-12, 1e-12, 1., 1., 1e-3, 1.])

    def prepare(self, a, b):
        psi, da, db = self.basis.values(a, b)
        q = self.q0*(a*a+b*b)
        state, m = self.mixing.state, self.mixing
        return dict(a=a, b=b, q=q, psi=psi, da=da, db=db,
            c0=state[0]*state[1]*np.exp(-q),
            h0=self.hc*np.exp(-q/self.section.thermal_width_ratio**2),
            u=m.wind*math.cos(state[3])+state[4]*np.exp(-self.section.velocity_shape_exponent*q))

    def fields(self, parameters, prepared):
        par = np.asarray(parameters, float)
        if par.shape != (self.count,) or not np.all(np.isfinite(par)):
            raise ValueError("invalid scalar shape parameters")
        size, p = self.basis.size, prepared
        lc, lh = p["psi"]@par[:size], p["psi"]@par[size:]
        if max(np.max(abs(lc), initial=0.), np.max(abs(lh), initial=0.)) > .1000000001:
            raise ValueError("shape correction exceeds the sampled 0.1 log trust region")
        c, h = p["c0"]*np.exp(lc), p["h0"]*np.exp(lh)
        rho, y, temperature = self.section.phase_inverse.state(c, h)
        hr, hc = self.section._phase_partials(rho, c)
        if np.any(hr >= 0.):
            raise ValueError("phase inverse must have a negative density slope")
        return dict(c=c, h=h, rho=rho, y=y, temperature=temperature, hr=hr, hc=hc, log_c=lc, log_h=lh)

    def grid(self, order):
        if order not in self._grids:
            x, w = roots_legendre(order)
            a, b = np.meshgrid(.5*(x+1), .5*(x+1), indexing="ij")
            p = self.prepare(a.ravel(), b.ravel())
            p["weights"] = (2*self.q0*self.mixing.area*np.outer(w, w)).ravel()
            p["r2"] = 2*self.q0*(self.mixing.sy**2*p["a"]**2+self.mixing.sn**2*p["b"]**2)
            self._grids[order] = p
        return self._grids[order]

    def moments(self, parameters, order=96, *, jacobian=False):
        p = self.grid(order)
        d = self.fields(parameters, p)
        rho, c, h, u, w = d["rho"], d["c"], d["h"], p["u"], p["weights"]
        values = np.array([w@(rho*u), w@(c*u), w@(rho*u*u), w@(h*u+.5*rho*u**3),
                           9.81*w@(self.section.rhoa-rho), w@(p["r2"]*h*u)])
        if not jacobian:
            return values
        psi = p["psi"]
        zero = np.zeros_like(psi)
        dc = np.column_stack([c[:, None]*psi, zero])
        dh = np.column_stack([zero, h[:, None]*psi])
        dr = (dh-d["hc"][:, None]*dc)/d["hr"][:, None]
        j = np.stack([(w*u)@dr, (w*u)@dc, (w*u*u)@dr,
                      (w*u)@dh+(.5*w*u**3)@dr, -9.81*w@dr, (w*p["r2"]*u)@dh])
        return values, j

    def edge(self, samples):
        if samples not in self._edges:
            # Fit at cosine-spaced locations, verify at a separate dense uniform grid.
            t = .5*(1-np.cos(np.linspace(0., math.pi, samples))) if samples == 65 else np.linspace(0., 1., samples)
            p = self.prepare(np.ones_like(t), t)
            base = self.reservoir.edge_fields(p["q"], self.baseline["family"])
            at = np.array([1., self.baseline["gamma"]])
            p["fm"], p["chi"], p["target_h"] = [base[k]@at for k in ("mass_flux", "chi_species", "reservoir_enthalpy")]
            if np.any(p["fm"] >= 0.) or np.any(p["chi"] <= 0.):
                raise ValueError("projection requires inflow and positive frozen mixing")
            p["scale_c"] = np.maximum(abs(base["local"]["y"]*p["fm"]), 1e-12)
            p["scale_h"] = np.maximum(abs(base["local"]["h"]/base["local"]["rho"]*p["fm"]), 1.)
            self._edges[samples] = p
        return self._edges[samples]

    def edge_residual(self, parameters, samples=65):
        p = self.edge(samples)
        d = self.fields(parameters, p)
        size = self.basis.size
        ca = d["c"]*(-2*self.q0+p["da"]@parameters[:size])
        ha = d["h"]*(-2*self.q0/self.section.thermal_width_ratio**2+p["da"]@parameters[size:])
        ra = (ha-d["hc"]*ca)/d["hr"]
        specific_h = d["h"]/d["rho"]
        ya, hsa = (ca-d["y"]*ra)/d["rho"], (ha-specific_h*ra)/d["rho"]
        conductance = self.mixing.area*d["rho"]*p["chi"]/(2*self.q0)
        fc = d["y"]*p["fm"]-conductance*ya
        fh = specific_h*p["fm"]-conductance*hsa
        return np.r_[fc/p["scale_c"], (fh-p["target_h"])/p["scale_h"]]

    def edge_jacobian(self, parameters):
        step = 1e-6
        eye = np.eye(self.count)*step
        return np.column_stack([(self.edge_residual(parameters+e)-self.edge_residual(parameters-e))/(2*step) for e in eye])

    def retract(self, parameters, order):
        par = np.array(parameters, float)
        for _ in range(10):
            values, j = self.moments(par, order, jacobian=True)
            residual = (values-self.target)/self.scales
            if max(abs(residual)) <= 1e-9:
                self.edge_residual(par)  # Include faces in the trust-region check.
                return par
            j = j/self.scales[:, None]
            step, _, rank, _ = np.linalg.lstsq(j, -residual, rcond=1e-12)
            if rank != 6:
                raise ValueError("six independent moments lost numerical rank")
            for factor in (.999, .5, .25, .125, .0625, .03125, .015625):
                try:
                    candidate = par+factor*step
                    error = (self.moments(candidate, order)-self.target)/self.scales
                    self.edge_residual(candidate)
                    if np.linalg.norm(error) < np.linalg.norm(residual):
                        par = candidate
                        break
                except (ValueError, RuntimeError):
                    pass
            else:
                raise ValueError("nonlinear moment retraction did not decrease its residual")
        raise ValueError("nonlinear moment retraction exceeded ten iterations")

    def fit(self, *, callback=None, maximum_iterations=25):
        parameters, history = np.zeros(self.count), []
        for order in (48, 96):
            parameters = self.retract(parameters, order)
            for iteration in range(maximum_iterations):
                r = self.edge_residual(parameters)
                values, jm = self.moments(parameters, order, jacobian=True)
                j = jm/self.scales[:, None]
                singular = np.linalg.svd(j, compute_uv=False)
                if singular[-1] <= 1e-12*singular[0]:
                    raise ValueError("moment tangent constraints are rank deficient")
                tangent = null_space(j, rcond=1e-12)
                je = self.edge_jacobian(parameters)
                step = tangent@np.linalg.lstsq(je@tangent, -r, rcond=1e-10)[0]
                norm = float(np.linalg.norm(r))
                row = dict(order=order, iteration=iteration, edge_rms=float(np.sqrt(np.mean(r*r))),
                    edge_max=float(max(abs(r))), moment_error=float(max(abs((values-self.target)/self.scales))),
                    moment_condition=float(singular[0]/singular[-1]))
                history.append(row)
                if callback is not None and iteration % 5 == 0:
                    callback(row)
                accepted = False
                for factor in (1., .5, .25, .125, .0625, .03125, .015625, .0078125):
                    try:
                        candidate = self.retract(parameters+factor*step, order)
                        candidate_r = self.edge_residual(candidate)
                        if np.linalg.norm(candidate_r) < norm*(1-1e-5*factor):
                            parameters, accepted = candidate, True
                            break
                    except (ValueError, RuntimeError):
                        pass
                if not accepted or max(abs(r)) < .005 or np.max(abs(step)) < 1e-10:
                    row["stopped"] = "no_acceptable_step" if not accepted else "converged_or_small_step"
                    break
        return dict(parameters=parameters, history=history)
