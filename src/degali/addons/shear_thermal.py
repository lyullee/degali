"""Reduced axial shear-work identity and an EXPLICIT thermal hypothesis test.

Not a closed curved Favre-RANS model or a downstream solver. Force density
and thermal/species diffusivity ratio are required. Shear production is not
automatically heat: only ``equilibrium_budgets`` tests that extra hypothesis.
"""

import math
import numpy as np


def affine_root(coefficients):
    """Solve a+b*gamma=0 without replacing a missing equation by gamma=0."""
    pair = np.asarray(coefficients, float)
    if pair.shape != (2,) or not np.all(np.isfinite(pair)):
        raise ValueError("an affine equation requires two finite coefficients")
    if pair[1] == 0.:
        raise ValueError("affine equation has no unique thermal-width rate")
    root = float(-pair[0]/pair[1])
    if not math.isfinite(root):
        raise ValueError("affine root is not finite")
    return root


class ReducedShearThermal:
    """Radial ALE stress reconstruction in the existing axial-velocity model.

    ``axial_force_density(q, local_data)`` supplies force per physical volume
    [N/m³] in the reduced even axial equation. No default buoyancy, pressure,
    gravity-energy treatment, or instantaneous conversion to heat is chosen.
    """

    def __init__(self, mixing, *, axial_force_density):
        if not callable(axial_force_density):
            raise TypeError("supply the reduced axial force density explicitly")
        self.mixing = mixing
        self.force_density = axial_force_density

    def axial_fields(self, q, family):
        m = self.mixing
        d = m.local(q)
        u, rho, h = d["u"], d["rho"], d["h"]
        directions = np.column_stack([family.origin, family.response])
        area_rate = directions[2]
        drho, dc, dh, du = [d[k]@directions for k in ("drho", "dc", "dh", "du")]
        mass = m.area*(drho*u[:, None]+rho[:, None]*du+(rho*u)[:, None]*area_rate)
        species = m.area*(dc*u[:, None]+d["c"][:, None]*du+(d["c"]*u)[:, None]*area_rate)
        momentum = m.area*(drho*(u*u)[:, None]+(2*rho*u)[:, None]*du+(rho*u*u)[:, None]*area_rate)
        enthalpy = m.area*(dh*u[:, None]+h[:, None]*du+(h*u)[:, None]*area_rate)
        kinetic = .5*m.area*(drho*(u**3)[:, None]+(3*rho*u*u)[:, None]*du+(rho*u**3)[:, None]*area_rate)
        force = np.broadcast_to(np.asarray(self.force_density(q, d), float), u.shape)
        if not np.all(np.isfinite(force)):
            raise ValueError("axial force density must be finite")
        uq = -m.section.velocity_shape_exponent*m.state[4]*np.exp(-m.section.velocity_shape_exponent*q)
        return dict(local=d, mass=mass, species=species, momentum=momentum,
                    enthalpy=enthalpy, kinetic=kinetic, force=m.area*force, uq=uq)

    def radial_fields(self, q, family, *, thermal_species_ratio, order=16):
        if not math.isfinite(thermal_species_ratio) or thermal_species_ratio <= 0.:
            raise ValueError("supply a positive finite thermal/species diffusivity ratio")
        q = np.asarray(q, float)
        m = self.mixing
        def integrand(qq):
            a = self.axial_fields(qq, family)
            bp = a["momentum"].copy()
            bp[:, 0] -= a["force"]
            return np.stack([a["mass"], a["species"], bp], axis=1)
        cumulative = m.partition.cumulative(integrand, q, order)
        flux = np.empty_like(cumulative)
        nonzero = q > 0.
        flux[nonzero] = -cumulative[nonzero]/(2*q[nonzero, None, None])
        if np.any(~nonzero):
            flux[~nonzero] = -.5*integrand(np.array([0.]))[0]
        a = self.axial_fields(q, family)
        d, uq = a["local"], a["uq"]
        fm, fc, fp = flux[:, 0], flux[:, 1], flux[:, 2]
        rp = fp-d["u"][:, None]*fm
        denominator = m.area*d["rho"]*d["yq"]
        if np.any(abs(denominator) < 1e-14) or np.any(uq == 0.):
            raise ValueError("gradients do not identify finite scalar/shear diffusivities")
        chi_c = -(fc-d["y"][:, None]*fm)/denominator[:, None]
        chi_p = -rp/(m.area*d["rho"]*uq)[:, None]
        fh = (d["h"]/d["rho"])[:, None]*fm-(m.area*d["rho"]*thermal_species_ratio*d["hq"])[:, None]*chi_c
        fk = d["u"][:, None]*fp-.5*(d["u"]**2)[:, None]*fm
        production = -(2*q*uq)[:, None]*rp
        work = np.column_stack([a["force"]*d["u"], np.zeros(len(q))])
        return dict(**a, mass_flux=fm, species_flux=fc, momentum_flux=fp,
                    stress=rp, chi_species=chi_c, chi_momentum=chi_p,
                    enthalpy_flux=fh, kinetic_flux=fk, production=production, work=work)

    def equilibrium_budgets(self, family, *, thermal_species_ratio, order=16):
        """Affine R0/R2 for the ADDITIONAL Q_H=P local-equilibrium hypothesis.

        The kinetic identity is independently meaningful without Q_H=P.
        Negative inferred diffusion is retained for rejection by the caller.
        Boundary fluxes are signed outward and already exclude face motion.
        """
        m, p = self.mixing, self.mixing.partition
        spread = m.sy*m.sy+m.sn*m.sn
        def fields(q):
            return self.radial_fields(q, family, thermal_species_ratio=thermal_species_ratio, order=order)
        def volume(q):
            d = fields(q)
            return np.stack([d["enthalpy"], d["kinetic"], d["production"], d["work"],
                spread*q[:, None]*d["enthalpy"], spread*q[:, None]*d["production"],
                2*spread*q[:, None]*d["enthalpy_flux"]], axis=1)
        def edge(q):
            d = fields(q)
            return np.stack([d["enthalpy_flux"], d["kinetic_flux"],
                             spread*q[:, None]*d["enthalpy_flux"]], axis=1)
        v = p.integrate(volume, order, square=True)
        e = 16*p.q0*p.edge_integrate(edge, order)
        result = dict(enthalpy_axial=v[0], kinetic_axial=v[1], production=v[2], body_work=v[3],
            weighted_enthalpy_axial=v[4], weighted_production=v[5], weighted_transverse=v[6],
            enthalpy_outward=e[0], kinetic_outward=e[1], weighted_enthalpy_outward=e[2])
        result["zeroth_residual"] = v[0]+e[0]-v[2]
        result["second_residual"] = v[4]+e[2]-v[6]-v[5]
        result["kinetic_identity"] = v[1]+e[1]-v[3]+v[2]
        result["total_energy_discrepancy"] = np.array([family.sources[4], 0.])-v[3]+e[1]+e[0]
        return result
