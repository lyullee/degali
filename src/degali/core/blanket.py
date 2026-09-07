"""The secondary source blanket.

Ports ``DEG1.for`` (the ``SRC1`` derivative routine and its ``SRC1O`` output
routine) and ``NOBL.FOR``.

What the model says
-------------------
A ground-level release competes with the atmosphere.  The atmosphere can only
take up so much contaminant from a patch of given streamwise extent; that
ceiling is

.. math::
    q^*_{max} = C_c \\, \\kappa\\, u_* (1+\\alpha)
                \\frac{\\delta_{lay}}{\\delta_{lay}-1} \\frac{1}{\\hat\\Phi(\\rho, L)}

If the source flux stays below it, the gas is swept away as fast as it is
produced and there is no blanket: :func:`no_blanket` handles that case.  If the
flux exceeds it, gas accumulates into a dense pool -- the *secondary source*
blanket -- which spreads sideways under its own weight until the extra area
brings take-up back into balance.

The blanket is a cylinder of radius :math:`R` and depth :math:`H`, centred on
the source, with spatially uniform properties.  Six states are integrated in
time:

===========  ================================================================
``R``        blanket radius, m
``m``        total mass, kg
``m_c``      contaminant mass, kg
``m_a``      dry-air mass, kg
``E``        total enthalpy, J
``P``        radial momentum, kg m/s (only when the van Ulden balance is on)
===========  ================================================================

Spreading
---------
Two regimes, following van Ulden (1983):

*Gravity slumping*, the quasi-steady limit,

.. math:: u_f = C_E \\sqrt{g \\frac{\\rho - \\rho_a}{\\rho_a} H}

with :math:`C_E = 1.15` from laboratory measurements.

*Momentum balance*, used while an instantaneous release is still accelerating
from rest.  Equation (V.9) of the user's guide balances the static pressure
force against form drag and the acceleration reaction of the displaced air,
with the cloud split into a head of depth :math:`H_h` and a tail of depth
:math:`H_t`.  The frontal velocity appears on both sides, so ``SRC1`` solves it
by a bracketed golden-section iteration; that iteration is reproduced here
because the convergence criterion (``rcrit = 0.002`` relative) leaves a
visible signature.  The balance is switched off, permanently, the first time
the momentum solution overtakes the slumping velocity.

Air entrainment
---------------
Air is drawn in at the advancing front at a rate set by the frontal Richardson
number :math:`Ri = g' H / u_f^2`:

.. math:: \\dot m_a = 2\\pi R H u_f \\rho_a \\frac{\\epsilon}{Ri}

with :math:`\\epsilon = 0.59`.  Entrainment through the *top* is not a separate
term: it is the take-up rate ``totrteout``, which leaves the blanket and
becomes the input to the downwind model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .atmosphere import friction_velocity
from .constants import GG, PI, VKC
from .entrainment import phi_hat, surface_exchange
from .thermo import MixtureState, Thermo

#: Index of each state in the integration vector, matching ``SRC1``'s
#: ``DATA iR/1/, mass/2/, massc/3/, massa/4/, iebal/5/, mbal/6/``.
I_R, I_M, I_MC, I_MA, I_E, I_P = range(6)


@dataclass
class BlanketParameters:
    """Numerical and closure constants, from the ``ER1`` parameter file."""

    ce: float = 1.15  #: gravity slumping coefficient
    epsilon: float = 0.59  #: frontal air entrainment coefficient
    delrmn: float = 0.0  #: stop spreading below this density excess
    dellay: float = 2.15  #: ratio of layer thickness to effective depth
    iphifl: int = 3  #: entrainment prescription
    srccut: float = 1.0e-5  #: minimum blanket depth, m
    srcss: float = 5.2  #: minimum time before declaring steady state, s
    srcoer: float = 0.007  #: output recording criterion
    #: van Ulden momentum-balance constants ``a_v``, ``b_v``, ``e_v``, ``d_v``.
    vua: float = 1.3
    vub: float = 1.2
    vue: float = 20.0
    vud: float = 0.64
    #: Relative convergence for the frontal-velocity iteration in ``SRC1``.
    rcrit: float = 0.002
    #: Half-width of the finite difference used for ``dR_p/dt``, s.
    delt: float = 0.1
    #: Offset added when clamping the blanket radius to the pool radius.
    #: ``SRC1`` writes ``Y(iR) = RADP + zero`` with ``zero = 1e-10``; the
    #: mutation is visible to the integrator, so it is reproduced.
    zero: float = 1.0e-10


@dataclass
class BlanketState:
    """Diagnostics evaluated alongside the derivatives, ``PRMT(6..23)``."""

    time: float
    radius: float
    height: float
    qstar: float  #: take-up flux, kg/(m**2 s)
    sz: float  #: equivalent vertical dispersion parameter, m
    mixture: MixtureState
    richardson: float
    velocity: float  #: frontal spreading velocity, m/s
    ht: float = 0.0  #: tail depth (van Ulden), m
    hh: float = 0.0  #: head depth (van Ulden), m
    rate: float = 0.0  #: primary source rate at this time, kg/s


@dataclass
class BlanketResult:
    """Output of a blanket integration."""

    history: list[BlanketState] = field(default_factory=list)
    #: Conditions handed to the downwind model (``/comss/`` in the Fortran).
    ess: float = 0.0  #: source strength, kg/s
    outcc: float = 0.0  #: secondary source concentration, kg/m**3
    outsz: float = 0.0  #: secondary source sigma_z, m
    outb: float = 0.0  #: secondary source half-width, m
    outl: float = 0.0  #: secondary source streamwise length, m
    swcl: float = 0.0  #: contaminant mass fraction
    swal: float = 0.0  #: air mass fraction
    senl: float = 0.0  #: enthalpy, J/kg
    srhl: float = 0.0  #: density, kg/m**3
    rm: float = 0.0  #: radius at peak take-up, m
    szm: float = 0.0  #: sigma_z at peak take-up, m
    emax: float = 0.0  #: peak take-up rate, kg/s
    rmax: float = 0.0  #: largest radius reached, m
    aleph: float = 0.0  #: wind-speed scaling constant used downwind
    tend: float = 0.0
    formed: bool = True  #: False when the atmosphere kept up and no blanket grew

    def as_array(self) -> np.ndarray:
        """History as the columns ``b9.scl`` prints."""
        return np.array(
            [
                [
                    s.time, s.radius, s.height, s.qstar, s.sz,
                    s.mixture.yc, s.mixture.rho, s.mixture.temp, s.richardson,
                ]
                for s in self.history
            ]
        )


class Blanket:
    """The ``DEG1`` secondary source model."""

    def __init__(self, case, thermo: Thermo, params: BlanketParameters | None = None):
        self.case = case
        self.th = thermo
        self.p = params or BlanketParameters()

        self.ustar = friction_velocity(case.u0, case.z0, case.zr, case.rml)
        self.alpha = None  # set by the caller, usually from fit_alpha
        self.rhoa = thermo.table.rhoa if thermo.table else None

    # -- closures ----------------------------------------------------------

    def max_takeup_flux(self, cc: float, rho: float, fetch: float) -> float:
        r"""``QSTRMX``: the atmospheric take-up flux over a fetch, kg/(m**2 s).

        Zero in a windless simulation -- with no wind there is nothing to
        carry the gas away, and the blanket only loses mass by dilution.
        """
        if self.case.u0 == 0.0:
            return 0.0
        alpha1 = self.alpha + 1.0
        ph = phi_hat(
            rho, fetch,
            rhoa=self.rhoa, alpha=self.alpha, ustar=self.ustar,
            u0=self.case.u0, z0=self.case.z0,
            gammaf=self.gammaf, dellay=self.p.dellay,
        )
        return (
            cc * VKC * self.ustar * alpha1
            * self.p.dellay / (self.p.dellay - 1.0) / ph
        )

    def sigma_z(self, qstar: float, fetch: float, cc: float) -> float:
        """``SZ``: the vertical dispersion parameter equivalent to the take-up.

        Inverts ``u_eff = q* L / C_c`` against the power-law profile, giving
        ``sz = z0 * (u_eff (1+alpha) / (u0 z0))**(1/(1+alpha))``.
        """
        if self.case.u0 == 0.0 or cc <= 0.0:
            return 0.0
        alpha1 = self.alpha + 1.0
        uheff = qstar * fetch / cc
        return (uheff * alpha1 / self.case.u0 / self.case.z0) ** (1.0 / alpha1) * self.case.z0

    # -- the frontal velocity ---------------------------------------------

    def _slumping_velocity(self, gprime: float) -> float:
        """Quasi-steady gravity current, ``CE * sqrt(g' H)``."""
        return self.p.ce * math.sqrt(gprime)

    def _momentum_velocity(
        self,
        y: np.ndarray,
        rho: float,
        hei: float,
        delrho: float,
        gprime: float,
        slump: float,
        vel_prev: float,
        head_developed: bool,
    ) -> tuple[float, float, float, float]:
        r"""Frontal velocity from the van Ulden momentum balance.

        Returns ``(velocity, dP/dt, H_t, H_h)``.

        The head depth follows from the velocity itself,
        :math:`H_h = u_f^2 / (C_E^2 g \Delta\rho/\rho_a)`, and the tail depth
        from volume conservation, so the momentum
        :math:`P = \int \rho u \, dV` closes only for the correct
        :math:`u_f`.  ``SRC1`` iterates on that with a golden-section bracket;
        the ratio 0.382 is reproduced because the iterate it converges to,
        and therefore the recorded output points, depend on it.
        """
        p = self.p
        R = y[I_R]
        vel = vel_prev
        velmin = 0.0
        velmax = max(slump, 0.1, vel)

        for _ in range(41):
            hh = vel * vel / p.ce / p.ce / GG / (delrho / self.rhoa)
            rh = R - p.vua * p.vub * hh
            value = R**2 / rh**2

            if head_developed:
                ht = 2.0 * (value * hei - p.vua * hh * (value - 1.0)) - hh
                denom = (
                    0.4 * PI * rho * (2.0 / 3.0 * ht + hh) * rh**3 / R
                    + 2.0 / 3.0 * PI * p.vua * rho * hh * (R**2 - rh**3 / R)
                    + p.vue * PI * R * hei**2 * self.rhoa
                )
                velc = y[I_P] / denom
                dpdt = (
                    PI * GG * delrho * R * ht**2
                    - p.vua * p.vud * PI * self.rhoa * R * hh * vel**2
                )
            else:
                ht = value * hei - p.vua * hh * (value - 1.0)
                denom = (
                    2.0 / 3.0 * PI * rho * ht * rh**3 / R
                    + 2.0 / 3.0 * PI * p.vua * rho * hh * (R**2 - rh**3 / R)
                    + p.vue * PI * R * hei**2 * self.rhoa
                )
                velc = y[I_P] / denom
                dpdt = (
                    PI * GG * delrho * (rh * ht**2 + p.vua * p.vub * hh**3)
                    - p.vua * p.vud * PI * self.rhoa * R * hh * vel**2
                )

            dif = abs(vel - velc)
            total = abs(vel) + abs(velc) + 1.0e-10
            if dif / total <= p.rcrit:
                return (vel + velc) / 2.0, dpdt, ht, hh

            if velc < velmin:
                velmin = max(velc, 0.0)
            if velc > velmax:
                velmax = velc
            if vel - velc > 0.0:
                velmax = vel
                vel = 0.382 * (velmax - velmin) + velmin
            else:
                velmin = velc
                vel = (1.0 - 0.382) * (velmax - velmin) + velmin

        # 40 iterations without convergence.  The Fortran accepts the bracket
        # if the cloud has accelerated past the slumping velocity -- that
        # happens when the integrator steps over the crossover -- and aborts
        # otherwise, with "STOP SRC1 velocity loop".
        #
        # That abort is reachable on real data.  A Thorney Island release is a
        # 14 m cylinder of Freon-air with a height-to-diameter ratio near 1,
        # and the momentum balance does not converge for it: the original
        # stops in exactly the same place, at the same step, so this is a
        # limit of the model rather than of the port.  It is raised rather
        # than papered over, because a frontal velocity that cannot be found
        # is not a small numerical inconvenience -- it is the whole gravity
        # current.
        vel = min(velmin, velmax)
        if vel > slump:
            return slump, 0.0, hei, hei
        raise RuntimeError(
            "SRC1 frontal-velocity iteration did not converge; the van Ulden "
            "momentum balance has no solution for this cloud shape "
            f"(H/D = {hei / (2.0 * R):.2f}). DEGADIS 2.1 stops here too."
        )

    # -- derivatives -------------------------------------------------------

    def derivatives(
        self, time: float, y: np.ndarray, ctx: dict
    ) -> tuple[np.ndarray, BlanketState]:
        """Port of ``SRC1``.

        ``ctx`` carries the mutable state the Fortran kept in ``PRMT``: the
        previous frontal velocity (the iteration's starting guess), the head
        and tail depths, and the latch that disables the momentum balance.
        """
        p = self.p
        src = self.case.source
        d = np.zeros(6)

        # --- composition and thermodynamics -------------------------------
        if y[I_M] <= 0.0:
            wc = max(ctx.get("wc_last", 1.0e-10), 1.0e-10)
            if wc > 1.0:
                wc = 1.0e-10
            wa = 1.0 - wc
            enth = wc * self.th.hmrte  # entrained air contributes nothing
        else:
            wc = y[I_MC] / y[I_M]
            wa = y[I_MA] / y[I_M]
            enth = y[I_E] / y[I_M]

        self.th.ambient.humsrc = (
            1.0 - wc - wa * (1.0 + self.th.ambient.humid)
        ) / wc
        mix = self.th.properties(wc, wa, enth, ifl=1)

        radp = src.radius_at(time)
        hei = max(y[I_M] / PI / y[I_R] / y[I_R] / mix.rho, 0.0)
        delrho = max(mix.rho - self.rhoa, 0.0)
        gprime = GG * delrho / self.rhoa * hei

        # --- spreading -----------------------------------------------------
        vel = 0.0
        airrte = 0.0
        ri = 0.0
        ht = hh = 0.0

        if gprime > 0.0:
            slump = self._slumping_velocity(gprime)
            if ctx.get("vuflag", False):
                vel, d[I_P], ht, hh = self._momentum_velocity(
                    y, mix.rho, hei, delrho, gprime, slump,
                    ctx.get("vel", 0.0), ctx.get("hh_last", 0.0) >= ctx.get("ht_last", 0.0),
                )
                ctx["vel"] = vel
                if vel >= slump:
                    # the momentum solution has caught up; latch it off
                    ctx["slump_reached"] = True
            else:
                vel = slump
                ht = hh = hei
            if vel > 0.0:
                ri = gprime / vel**2
                airrte = 2.0 * PI * p.epsilon / ri * self.rhoa * y[I_R] * hei * vel
                d[I_R] = vel

        # --- constraints on the radius ------------------------------------
        if delrho < p.delrmn and self.case.u0 != 0.0:
            d[I_R] = 0.0
            airrte = 0.0

        area = PI * (y[I_R] ** 2 - radp**2)
        if y[I_R] <= radp:
            # the blanket may never be smaller than the pool feeding it; it
            # then grows at the pool's own rate, by central difference
            area = 0.0
            y[I_R] = radp + p.zero  # SRC1 mutates the state vector here
            h = p.delt / 2.0
            if time > h:
                d[I_R] = max(
                    0.0,
                    (src.radius_at(time + h) - src.radius_at(time - h)) / p.delt,
                )
            else:
                d[I_R] = max(0.0, (src.radius_at(h) - src.radius_at(0.0)) / h)

        # --- fluxes ---------------------------------------------------------
        erate = src.rate_at(time)
        pwcp = src.wc_at(time)
        totprt = erate / pwcp
        hprim = src.enthalpy_at(time)
        fetch = 2.0 * y[I_R]

        qstrmx = self.max_takeup_flux(mix.cc, mix.rho, fetch)
        qstrll = qstrmx * PI * fetch * fetch / 4.0
        totrteout = qstrll / wc

        watrte, qrte = surface_exchange(
            mix.temp, hei, mix.rho, mix.wm, mix.cp, mix.yw,
            tsurf=self.th.ambient.tsurf, pamb=self.th.ambient.pamb,
            zr=self.case.zr, u0=self.case.u0, z0=self.case.z0,
            alpha=self.alpha, ustar=self.ustar,
            isofl=self.th.ambient.isofl, ihtfl=self.th.ambient.ihtfl,
            iwtfl=self.th.ambient.iwtfl, htco=self.th.ambient.htco,
            wtco=self.th.ambient.wtco,
            vapour_pressure=self.th.backend.water_vapour_pressure,
        )
        qrte = max(area * qrte, 0.0)  # never let the ground cool the cloud
        watrte = area * watrte

        totrtein = airrte + totprt + watrte

        # --- shrinking blanket ---------------------------------------------
        if totrtein < totrteout and not self.case.instantaneous:
            d[I_R] = 0.0
            if hei > p.srccut and y[I_R] > p.srccut:
                dhdt = (totrtein - totrteout) / 3.0 / PI / y[I_R] ** 2 / mix.rho
                d[I_R] = y[I_R] / hei * dhdt
            if y[I_R] <= 0.01 * ctx.get("rmax", y[I_R]) and erate == 0.0:
                d[I_R] = 0.0
            if y[I_R] <= hei and erate == 0.0:
                d[I_R] = 0.0

        # --- mass and energy balances ---------------------------------------
        humid = self.th.ambient.humid
        d[I_M] = totrtein - totrteout
        d[I_MC] = erate - qstrll
        d[I_MA] = (airrte + erate * (1.0 / pwcp - 1.0)) / (1.0 + humid) - (
            wa / wc * qstrll
        )
        d[I_E] = 0.0
        if self.th.ambient.ihtfl != 0:
            # for ihtfl == 0 the enthalpy follows the adiabatic mixing line
            # and TPROP supplies it, so no balance is integrated
            d[I_E] = (
                hprim * totprt
                + self.th.harte * airrte
                + self.th.hwrte * watrte
                - mix.enthalpy * totrteout
                + qrte
            )

        sz = self.sigma_z(qstrmx, fetch, mix.cc)

        ctx["wc_last"] = wc
        ctx["ht_last"], ctx["hh_last"] = ht, hh

        state = BlanketState(
            time=time, radius=float(y[I_R]), height=hei, qstar=qstrmx, sz=sz,
            mixture=mix, richardson=ri, velocity=vel, ht=ht, hh=hh, rate=erate,
        )
        return d, state
