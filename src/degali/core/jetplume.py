"""Elevated and pressurised releases: a port of ``JETPLU``.

Ports ``JETPLU.FOR`` (the ``MODEL``, ``MODOUT`` and ``SETJET`` routines),
``ELLIPS.FOR`` and the matrix solve in ``SIMUL.FOR``.

Everything else in this package assumes the contaminant arrives at ground
level with no momentum of its own.  A pressurised release does not: it leaves
as a jet, and where it goes is decided by its own momentum, by buoyancy, and
by the wind bending it over.  ``JETPLU`` is an integral model of that
trajectory in the style of Ooms, and it hands over to ``DEG1`` at the point
the plume touches down -- or reports that it never does.

What is integrated
------------------
Six quantities against arc length ``s`` along the plume axis:

=========  =============================================================
``cc``     centreline concentration, kg/m**3
``sysz``   product of the cross-section's Gaussian half-widths
``theta``  trajectory angle above horizontal
``uc``     excess centreline velocity over the local wind
``x``      downwind distance
``zj``     centreline elevation
=========  =============================================================

The first four come from four coupled balances -- contaminant mass, total
mass, and momentum in ``x`` and ``z`` -- which are *not* separable: each
balance involves the derivatives of all four states.  So the right-hand side
is assembled as a 4x4 linear system and solved at every step. That is what
``SIMUL`` did with Gauss elimination and what :func:`numpy.linalg.solve` does
here.

``x`` and ``z`` follow directly from the angle, ``dx/ds = cos(theta)`` and
``dz/ds = sin(theta)``.

The cross-section
-----------------
The plume cross-section is an ellipse whose *area* is fixed by ``sysz`` but
whose aspect ratio is not.  DEGADIS closes it by requiring that the plume's
excess spread over the ambient be isotropic,

.. math:: S_y^2 - S_z^2 = \\sigma_{ya}^2 - \\sigma_{za}^2

so as the ambient dispersion becomes anisotropic the plume follows.  With
``sysz`` known this is one equation in one unknown, solved by root-finding at
every derivative call.

Entrainment has three parts: shear from the excess velocity (``ALFA1``),
cross-flow from the wind blowing across an inclined plume (``ALFA2``), and
the growth of the ambient sigmas themselves. Form drag on the inclined
cross-section appears in both momentum balances.

Termination
-----------
Two ways out. The plume descends to ``z = 0``, in which case the state is
interpolated back to touchdown, the concentration is doubled for the ground
reflection, and ``DEGBRIDG`` takes over. Or the maximum ground-level mole
fraction falls below the level of concern while still airborne, in which case
the release never produces a ground-level hazard and the run stops there.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import ellipe, ellipeinc, erf

from .atmosphere import psi
from .constants import PI, RT2, VKC
from .entrainment import phi
from .numerics import limit, zbrent
from .rkgst import Control, rkgst

#: State indices in the jet integration.
J_CC, J_SYSZ, J_THETA, J_UC, J_X, J_Z = range(6)


#: Classical drag coefficients, for the shape interpolation in
#: :attr:`JetCoefficients.shape_drag`: a flat plate normal to the flow, and a
#: circular cylinder in crossflow.
CD_PLATE = 2.0
CD_CIRCULAR = 1.2


@dataclass
class JetCoefficients:
    """Closure constants, from the ``DATA`` statements in ``JETPLU``."""

    alfa1: float = 0.057  #: shear entrainment coefficient (read from the deck)
    #: Use Sandia's source-momentum entrainment in place of JETPLU's local
    #: perimeter/shear term.  HyRAM+ defines
    #: ``E_mom = beta_A sqrt(A0 rho0 v0^2 / rho_a)`` with ``beta_A=0.28``.
    #: It is a constant volumetric entrainment flux fixed by the atmospheric
    #: source plane.  Zero preserves JETPLU exactly.
    momentum_entrainment_beta: float = 0.0
    #: Move the shear coefficient from its jet value towards its plume value
    #: as buoyancy takes over.
    #:
    #: Papanicolaou and List (1988) measured 0.0545 for a pure jet and 0.0875
    #: for a pure plume; ``JETPLU`` uses 0.057 throughout, at the jet end. A
    #: flashing liquid hydrogen release is buoyant within a metre of the
    #: nozzle, so it spends almost all of its measured range in the plume
    #: regime while being entrained at the jet rate.
    #:
    #: The transition uses the same experiment's pure-plume Richardson number,
    #: 0.716. The often-quoted 0.0533, 0.0833 and 0.557 are Fischer et al.'s
    #: earlier proposed values, quoted for comparison in that paper rather
    #: than measured by Papanicolaou and List (see ``docs/references.md``).
    #: Zero leaves ``alfa1`` alone, which is what every parity test needs.
    plume_transition: bool = False
    #: Scale the shear entrainment by the square root of the density ratio.
    #:
    #: Ricou and Spalding (1961) measured that the entrained mass flux of a
    #: jet scales as ``sqrt(rho_ambient / rho_jet)`` -- *inversely* with the
    #: square root of the jet density, as Panda and Hecht restate for
    #: cryogenic hydrogen: "entrained mass flow increases with distance from
    #: the nozzle exit, scaling inversely by the square root of the density of
    #: the jet at the nozzle".
    #:
    #: The direction matters and is easy to invert. A jet *lighter* than
    #: ambient entrains more, not less, and a flashing liquid hydrogen jet
    #: becomes lighter than air within a metre of the nozzle -- so this raises
    #: the entrainment over most of the measured range, while lowering it in
    #: the dense first metre. Applying the reciprocal by mistake moved the
    #: vertical spread the wrong way, from 0.70 of the measured value to 0.65.
    #:
    #: ``JETPLU`` has no such scaling: it was written for releases close to
    #: ambient density. Zero leaves the coefficient alone.
    density_scaled_entrainment: bool = False
    #: Extra drag on the rising plume, opposing the vertical velocity.
    #:
    #: Mack and co-workers report that "for strongly buoyant plumes, it was
    #: observed that vertical plume speeds are overpredicted by some integral
    #: models", and that the net buoyancy force is hard to model because it is
    #: formed from averaged plume variables: the volume may not contribute in
    #: full because the density varies along the axis, and the entrainment
    #: into large vortical structures is not represented.
    #:
    #: Some models correct this with an added-mass term. EFFECTS does not --
    #: it adds a drag that depends on the plume's *vertical velocity and
    #: shape* rather than on the buoyant volume, which is the distinction that
    #: matters: added mass has been pre-registered and falsified elsewhere on
    #: wind-tunnel data across six source widths, and Tickle reports the same
    #: from DRIFT.
    #:
    #: ``JETPLU``'s own drag acts along the trajectory and is proportional to
    #: ``(u sin theta)**2``; this term is separate and acts only on the rise.
    #: Zero recovers the shipped behaviour.
    rise_drag: float = 0.0
    #: Include the vertical velocity in the shear entrainment velocity.
    #:
    #: ``JETPLU`` entrains on ``max(u_c, 0)``, the *horizontal* excess over
    #: the wind. Mack and co-workers record making exactly this change to
    #: EFFECTS: "the shear term contribution was extended with the vertical
    #: velocity component which was neglected in dU in the original model.
    #: This was due to the focus on heavy gas applications where horizontal
    #: entrainment due to plume vertical velocity plays a minor role."
    #:
    #: For a buoyant plume it is not minor: the faster it rises, the more air
    #: it drags in, which dilutes the buoyancy that is lifting it. That is a
    #: negative feedback the shipped model has no access to, and it acts on
    #: the rise through the entrainment rather than through a drag force --
    #: which matters here because the drag term is quadratic in a rise
    #: velocity of a few centimetres per second and does almost nothing.
    vertical_shear: bool = False
    #: Coefficient on entrainment by the plume's own buoyant velocity scale;
    #: zero disables it.
    #:
    #: **This is the one mechanism AEA Technology recommended and nobody
    #: tried.** Having found that suppressing rise through the vertical
    #: velocity spoils the dilution, they wrote:
    #:
    #:     Modifying the entrainment formula to depend upon a buoyant velocity
    #:     scale rather than vertical component of velocity may circumvent
    #:     such problems
    #:
    #: The scale is `sqrt(g' H)` with `g'` the reduced gravity and `H` the
    #: plume depth. It is added to the shear velocity in quadrature, so the
    #: larger of the two governs.
    #:
    #: What makes it different from every other mechanism here is **what it
    #: does to the wind dependence**. The shear scale falls with the wind and
    #: the buoyant scale does not, so in a light wind the buoyant term
    #: dominates and dilutes the buoyancy that is doing the lifting, while in
    #: a strong wind it is negligible. Every mechanism tried before this one
    #: changes the rise by roughly the same factor at 2 m/s and at 6 m/s, and
    #: the error does not: the model over-predicts the plume height by 0.66 at
    #: 6.3 m/s and by more than 2.2 at 1.55 m/s.
    buoyant_entrainment: float = 0.0
    #: Add the published Houf--Schefer buoyancy-driven entrainment closure.
    #:
    #: HyRAM adds an entrainment flux based on the local densimetric Froude
    #: number, with a source-Froude coefficient, then limits the combined
    #: momentum-plus-buoyancy entrainment coefficient to the pure-plume value
    #: 0.082. The equations and constants are reproduced here; see
    #: ``docs/prereg-houf-entrainment.md``. False preserves JETPLU exactly.
    houf_buoyant_entrainment: bool = False
    #: Add DEGADIS's atmospheric-surface-layer entrainment through the top
    #: of the part of a plume that is lying on the ground.
    #:
    #: ``ground_effect`` removes the buried curved perimeter, but the original
    #: JETPLU equations have no replacement for the top-surface mixing used by
    #: DEGADIS's ground-layer path.  This switch adds that existing
    #: friction-velocity/Richardson closure over the ground-contact chord.  It
    #: has no effect unless ``JetPlume.ground_effect`` is also true, and is off
    #: by default to preserve the Fortran oracle.
    ground_layer_entrainment: bool = False
    #: Coefficient on the shape-dependent pressure drag; zero disables it.
    #:
    #: `rise_drag` is a fixed coefficient on `w|w|` acting on the perimeter.
    #: Mack records EFFECTS using a drag that depends on the plume's **shape**
    #: instead -- a cylinder in crossflow with aspect ratio `AR = H / B` --
    #: and that dependence is the whole of AEA Technology's wind-tunnel
    #: finding, which is that rise suppression is worst for **wide** sources.
    #: A fixed coefficient cannot express it.
    #:
    #: Moving vertically, a plume that is wide and flat is a bluff body and a
    #: compact one is not. The endpoints are classical: a flat plate normal to
    #: the flow has `Cd` near 2.0 and a circular cylinder near 1.2. The
    #: interpolation between them is a choice made here and is not from a
    #: source, so it is a shape *structure* with a fitted blend rather than a
    #: measured law, and must be described that way.
    #:
    #: The force acts on the plume's width, `2 delta sigma_y`, because that is
    #: what it presents to vertical motion -- not on the perimeter, which is
    #: what `rise_drag` uses.
    shape_drag: float = 0.0
    #: Critical bulk Richardson number for lift-off; zero disables the test.
    #:
    #: **Falsified on the Spadeadam tests and kept for the record.** The
    #: shipped model has no lift-off criterion, and adding Briggs' changes
    #: nothing, because `Ri*` runs 187 to 363 over the whole trajectory
    #: against a threshold of 20 to 30. The criterion agrees the plume should
    #: lift. The fault is the rise *rate*, not its onset.
    #:
    #: **The shipped model has no lift-off criterion at all.** It lets a
    #: buoyant plume rise from the source, and every model in AEA
    #: Technology's survey of the problem has a threshold instead: Briggs
    #: proposed a bulk Richardson number `Lp` with a critical value near 30,
    #: and HGSYSTEM-MMES holds the plume on the ground until `Lp` reaches 20.
    #:
    #: The parameter here is that threshold, on
    #:
    #:     Ri* = g H (rho_a - rho_m) / (rho_a u*^2)
    #:
    #: with `H` the plume depth, `rho_m` the plume density and `u*` the
    #: friction velocity. Below it the plume is held down; above it the
    #: vertical momentum equation runs as shipped.
    #:
    #: This acts on *when* the rise starts, not on how fast, so it touches no
    #: entrainment term -- which is what distinguishes it from the two
    #: mechanisms AEA warns will cost dilution.
    liftoff_richardson: float = 0.0
    alfa2: float = 0.5  #: cross-flow entrainment coefficient
    sc: float = 1.42  #: turbulent Schmidt number
    cd: float = 0.2  #: drag coefficient
    delta: float = 2.15  #: ratio of the ellipse semi-axis to sigma

    #: Richardson number at which the transition reaches the plume value.
    RI_PLUME = 0.716
    ALPHA_JET = 0.0545
    ALPHA_PLUME = 0.0875
    HOUF_ALPHA_LIMIT = 0.082

    @staticmethod
    def houf_buoyancy_coefficient(source_froude: float) -> float:
        """Houf--Schefer buoyancy coefficient from source Froude number."""
        if source_froude < 268.0:
            return (
                17.313
                - 0.11665 * source_froude
                + 2.0771e-4 * source_froude**2
            )
        return 0.97

    def shear_coefficient(
        self, richardson: float, density_ratio: float = 1.0
    ) -> float:
        """``alfa1`` at this local Richardson number and density ratio."""
        alpha = self.alfa1
        if self.plume_transition:
            ratio = min(abs(richardson) / self.RI_PLUME, 1.0)
            alpha = self.ALPHA_JET + (self.ALPHA_PLUME - self.ALPHA_JET) * ratio**2
        if self.density_scaled_entrainment:
            alpha /= math.sqrt(max(density_ratio, 1e-12))
        return alpha

    def profile_constants(self) -> tuple[float, ...]:
        r"""``RK1``..``RK6``: the shape factors of the Gaussian profiles.

        Each is the integral of a product of Gaussians over the elliptical
        cross-section, truncated at ``delta`` sigma. ``RK3`` is the plain
        area factor; the rest carry an :func:`~scipy.special.erf` because the
        truncation is at finite width.
        """
        return (
            self.profile_integral(1.0),
            self.profile_integral(1.0 + self.sc),
            self.profile_integral(0.0),
            self.profile_integral(self.sc),
            self.profile_integral(2.0 * self.sc),
            self.profile_integral(1.0 + 2.0 * self.sc),
        )

    def profile_integral(self, exponent: float) -> float:
        """Area factor for one Gaussian exponent in JETPLU's section."""
        if exponent < 0.0:
            raise ValueError("Gaussian exponent must be non-negative")
        if exponent == 0.0:
            return PI * self.delta**2
        ppp = math.sqrt(PI) / 2.0
        q = erf(self.delta / RT2 * math.sqrt(exponent) * ppp)
        return 2.0 * PI / exponent * q**2


def ellipse(yaxis: float, zaxis: float, zseg: float) -> tuple[float, float]:
    """Port of ``ELLIPS``: curved perimeter and area of an elliptical segment.

    ``zseg`` is how much of the ellipse is present above the cut: equal to
    ``zaxis`` for a whole ellipse, zero for a half. The original evaluates the
    elliptic integrals with the arithmetic-geometric mean iteration from
    Spanier and Oldham; :func:`scipy.special.ellipe` and
    :func:`scipy.special.ellipeinc` are the same functions to machine
    precision.

    Returns ``(perimeter, area)``.
    """
    if zaxis <= yaxis:
        aaa, bbb = yaxis, zaxis / yaxis
    else:
        aaa, bbb = zaxis, yaxis / zaxis
    m = max(1.0 - bbb**2, 0.0)  # scipy takes the parameter m = k^2

    if zseg == zaxis:
        return 4.0 * aaa * float(ellipe(m)), PI * yaxis * zaxis
    if zseg == 0.0:
        return 2.0 * aaa * float(ellipe(m)), PI * yaxis * zaxis / 2.0

    aaa = zaxis
    bbb = yaxis / zaxis
    ratio = zseg / zaxis
    phi = math.asin(ratio)
    area = yaxis * zaxis * (PI / 2.0 + phi + ratio * math.sqrt(1.0 - ratio**2))

    if bbb <= 1.0:
        m = 1.0 - bbb**2
        circum = 2.0 * aaa * (float(ellipeinc(phi, m)) + float(ellipe(m)))
    else:
        m = 1.0 - (1.0 / bbb) ** 2
        phi = PI / 2.0 - phi
        circum = 2.0 * aaa * bbb * (
            2.0 * float(ellipe(m)) - float(ellipeinc(phi, m))
        )
    return circum, area


@dataclass
class JetIntegralFluxes:
    """Integral transport represented by one JETPLU similarity state."""

    total_mass: float  #: kg/s
    contaminant_mass: float  #: kg/s
    momentum_x: float  #: N
    momentum_z: float  #: N
    energy: float = math.nan  #: W, relative to ambient enthalpy


@dataclass
class JetState:
    """Diagnostics evaluated alongside the derivatives."""

    s: float  #: arc length, m
    x: float
    z: float
    cc: float
    sy: float
    sz: float
    theta: float
    uc: float
    ua: float  #: wind speed averaged over the cross-section
    rho: float
    temp: float
    yc: float  #: centreline mole fraction
    rate: float  #: contaminant flux through the cross-section, kg/s


@dataclass
class JetResult:
    """Outcome of a jet/plume integration."""

    rows: np.ndarray = field(default_factory=lambda: np.zeros((0, 12)))
    touchdown: bool = False
    #: Downwind distance at touchdown, m; ``0`` when the plume never lands,
    #: which is the signal ``JETPLU`` writes to stop ``DEGADIS`` continuing.
    distance: float = 0.0
    concentration: float = 0.0  #: ground-level concentration at touchdown
    halfwidth: float = 0.0  #: ``delta * sy`` at touchdown, m


class JetPlume:
    """The ``JETPLU`` model."""

    def __init__(
        self, thermo, *, coefficients: JetCoefficients, u0, z0, zr, rml,
        ustar, rhoa, rhoe, deltay, betay, deltaz, betaz, gammaz, yclow,
        ground_effect: bool = False,
        liquid_fraction: float = 0.0,
        evaporation_ratio: float = 0.0,
        non_boussinesq: bool = False,
        spread_floor: bool = False,
    ):
        self.th = thermo
        self.k = coefficients
        self.u0, self.z0, self.zr, self.rml = u0, z0, zr, rml
        self.ustar = ustar
        self.rhoa, self.rhoe = rhoa, rhoe
        self.deltay, self.betay = deltay, betay
        self.deltaz, self.betaz, self.gammaz = deltaz, betaz, gammaz
        self.yclow = yclow
        self.rk = coefficients.profile_constants()
        #: Whether a plume in contact with the ground is treated as such.
        #:
        #: ``JETPLU`` models a *free* plume until it touches down: it entrains
        #: through its whole perimeter and rises without restraint. A jet
        #: released half a metre up, whose vertical spread has grown past that,
        #: is already lying on the ground -- which blocks entrainment on the
        #: underside and holds the plume down. Leaving that out is what makes
        #: the model lose the PRESLHY plume out of the 0.5 m sensor plane by
        #: 6 m, predicting 7 mole per cent where 32 was measured.
        #:
        #: Off by default, because every parity test against the Fortran
        #: depends on the free-plume behaviour.
        self.ground_effect = ground_effect
        #: Mass fraction of the released contaminant still liquid at the
        #: orifice, and the masses of air needed to evaporate one of it.
        #:
        #: A flashing cryogen leaves as a droplet-laden two-phase fluid. For
        #: liquid hydrogen at 5 barg only a quarter has flashed, so the bulk
        #: density at the orifice is about 5.3 kg/m**3 -- four times ambient --
        #: while the *gas* mixing line, which is all an integral model carries,
        #: tops out at 1.09 times ambient.
        #:
        #: Reading buoyancy off the gas line alone launches a jet four times
        #: denser than air as though it were barely denser, so a horizontal
        #: release climbs where it should fall. On the PRESLHY trials released
        #: at 1.5 m the sensors a metre below the axis measured 43 to 98 vol %
        #: while the model gave 0.0 to 0.2: the real jet reached the ground and
        #: the modelled one climbed away.
        #:
        #: Droplets are carried as extra mass with negligible volume,
        #: evaporating as entrained air supplies their latent heat. Zero
        #: recovers the single-phase behaviour exactly, so every parity test
        #: against the Fortran is unaffected.
        self.liquid_fraction = liquid_fraction
        self.evaporation_ratio = evaporation_ratio
        #: Divide the buoyancy by the parcel density rather than the ambient.
        #:
        #: ``JETPLU`` writes the vertical momentum source linear in the
        #: density difference, which is the Boussinesq form and accelerates a
        #: parcel at ``g (rho_a - rho)/rho_a``. The non-Boussinesq form uses
        #: the parcel's own density, ``g (rho_a - rho)/rho``. Xiao and
        #: co-workers report the Boussinesq approximation failing for hydrogen
        #: when buoyancy is comparable with momentum, which is this regime.
        #:
        #: For a plume *lighter* than ambient that makes the acceleration
        #: larger, not smaller, so it cannot be the fix for a model that
        #: already rises too fast -- see ``docs/prereg-boussinesq.md``, where
        #: that was predicted before it was run. It is here because the term
        #: is worth having correct for releases denser than air, where it
        #: increases the sink rate.
        self.non_boussinesq = non_boussinesq
        #: Require the lateral spread to be at least the passive one.
        #:
        #: ``JETPLU`` splits the product ``sigma_y sigma_z`` by requiring the
        #: plume's excess over ambient to be isotropic,
        #: ``sy**2 - sz**2 = sya**2 - sza**2``. Measured against FLADIS, where
        #: both spreads are resolved on the same arc, that gives
        #: ``sigma_y/sigma_z = 3.03`` against a measured 1.66 to 1.91 -- the
        #: modelled section is too flat.
        #:
        #: Hart and Harper report the same failure in the UDM from a different
        #: route: its heavy-gas spread rate goes as the square root of the
        #: density excess, which is zero for a buoyant cloud, so "any
        #: entrainment of air into the cloud therefore is constrained to make
        #: the cloud taller rather than wider". Their fix is to floor the
        #: spread rate at the passive value, and this is the same idea applied
        #: to the split: ``sigma_y`` may not fall below the ambient lateral
        #: spread.
        self.spread_floor = spread_floor
        #: Filled from a consistent source plane by the initial-condition
        #: routine. The coefficient depends on source, not local, Froude.
        self.houf_source_froude: float | None = None
        self._houf_alpha_buoy: float | None = None
        self._source_momentum_entrainment: float | None = None

    # -- the cross-section --------------------------------------------------

    def _split(self, sysz: float, sya: float, sza: float) -> tuple[float, float]:
        """Split ``sysz`` into ``sy`` and ``sz``.

        Port of the ``SYFUN`` root find: the plume's excess spread over the
        ambient is required to be isotropic, ``sy^2 - sz^2 = sya^2 - sza^2``.
        Written as a difference of products in the original to avoid
        cancellation, and reproduced that way.
        """

        def f(sy: float) -> float:
            sz = sysz / sy
            return (sya - sza) * (sya + sza) + (sz - sy) * (sz + sy)

        start = math.sqrt(sysz)
        hi, lo = limit(f, start, start / 20.0, max(sysz, 10.0), 1.0e-10)
        sy = zbrent(f, lo, hi, 1.0e-4)
        if self.spread_floor and sy < sya:
            # the plume cannot be narrower across the wind than the ambient
            # turbulence alone would make it
            sy = sya
        return sy, sysz / sy

    def _wind_profile(
        self, zj: float, sz: float, ct: float
    ) -> tuple[float, float]:
        """Wind speed averaged over the plume's vertical extent.

        A local power-law exponent is fitted from two points on the log
        profile -- at the plume top and at its centre -- and the resulting
        profile is integrated between top and bottom. The 0.01 m offset on the
        top height is in the original to keep the exponent from dividing by
        zero when the plume is thin.
        """
        # A quiescent ambient is a valid free-jet boundary condition.  The
        # logarithmic profile below is undefined at ``ustar == 0`` because it
        # tries to form ``log(0/0)`` even though both the local and integrated
        # wind speeds are exactly zero.  JETPLU was written for atmospheric
        # cross-wind releases and never met this laboratory limit.
        if self.ustar <= 0.0:
            return 0.0, 0.0

        zj = max(zj, self.zr)
        ztop = zj + self.k.delta * sz * ct + 0.01
        utop = self.ustar / VKC * (
            math.log((ztop + self.zr) / self.zr) - psi(ztop, self.rml)
        )
        umid = self.ustar / VKC * (
            math.log((zj + self.zr) / self.zr) - psi(zj, self.rml)
        )
        alpha = math.log(utop / umid) / math.log(ztop / zj)
        alpha1 = 1.0 + alpha
        zbot = max(zj - self.k.delta * sz * ct, self.zr)
        averaged = (
            umid / (alpha1 * (ztop - zbot) * zj**alpha)
            * (ztop**alpha1 - zbot**alpha1)
        )
        return averaged, alpha

    def _wind(self, zj: float, sz: float, ct: float) -> float:
        """Compatibility wrapper returning only the averaged wind speed."""
        return self._wind_profile(zj, sz, ct)[0]

    def integral_fluxes(
        self, y: np.ndarray, *, include_energy: bool = True,
        energy_quadrature: int = 32,
    ) -> JetIntegralFluxes:
        """Return the conserved fluxes carried by a similarity state.

        The four balance quantities use the exact ``RK1``--``RK6`` shape
        factors used by :meth:`derivatives`. The optional energy diagnostic is
        analytic as well: the adiabatic-line enthalpy flux is proportional to
        contaminant flux, and the cubic kinetic term uses two additional
        Gaussian factors. JETPLU does not transport energy as an ODE state;
        exposing it here makes that omission measurable at a model handoff
        instead of silently accepting it.
        """
        if np.shape(y) != (6,):
            raise ValueError("a JETPLU state must contain six variables")
        cc, sysz, theta, uc, dist, zj = (
            float(y[J_CC]), abs(float(y[J_SYSZ])), float(y[J_THETA]),
            float(y[J_UC]), float(y[J_X]), float(y[J_Z]),
        )
        if cc <= 0.0 or sysz <= 0.0 or not np.all(np.isfinite(y)):
            raise ValueError("JETPLU concentration, area and state must be physical")

        rk1, rk2, rk3, rk4, rk5, rk6 = self.rk
        st, ct = math.sin(theta), math.cos(theta)
        sya = self.deltay * max(dist, 0.0) ** self.betay
        if dist > 0.0:
            sza = (
                self.deltaz * dist**self.betaz
                * math.exp(self.gammaz * math.log(dist) ** 2)
            )
        else:
            sza = 0.0
        sy, sz = self._split(sysz, sya, sza)
        ua = self._wind(zj, sz, ct)
        mix = self.th.table.from_concentration(cc)
        gamma = (mix.rho - self.rhoa) / cc

        axial_ambient = ua * ct
        contaminant = cc * sysz * (
            rk1 * axial_ambient + rk2 * uc
        )
        total_mass = self.rhoa * sysz * (
            rk3 * axial_ambient + rk4 * uc
        ) + gamma * contaminant
        momentum = self.rhoa * sysz * (
            rk3 * axial_ambient**2
            + 2.0 * rk4 * axial_ambient * uc
            + rk5 * uc**2
        ) + gamma * cc * sysz * (
            rk1 * axial_ambient**2
            + 2.0 * rk2 * axial_ambient * uc
            + rk6 * uc**2
        )

        energy = math.nan
        if include_energy:
            table_fraction = self.th.table.cc / self.th.table.rho
            source_fraction = float(table_fraction[-1])
            if source_fraction <= 0.0:
                raise ValueError("mixing table has no contaminant source state")
            positive = table_fraction > max(source_fraction * 1.0e-10, 1.0e-14)
            enthalpy_per_contaminant = (
                (self.th.table.h[positive] - self.th.table.h[0])
                / table_fraction[positive]
            )
            linear_enthalpy = bool(np.allclose(
                enthalpy_per_contaminant,
                enthalpy_per_contaminant[-1],
                rtol=1.0e-8,
                atol=max(abs(float(enthalpy_per_contaminant[-1])) * 1.0e-10, 1.0e-8),
            ))
            if linear_enthalpy:
                enthalpy_flux = float(
                    enthalpy_per_contaminant[-1] * contaminant
                )
            else:
                if energy_quadrature < 8:
                    raise ValueError(
                        "nonlinear energy quadrature requires at least eight points"
                    )
                nodes, weights = np.polynomial.legendre.leggauss(
                    energy_quadrature
                )
                coords = self.k.delta * nodes
                grid_y, grid_z = np.meshgrid(coords, coords, indexing="ij")
                weight_y, weight_z = np.meshgrid(
                    weights, weights, indexing="ij"
                )
                exponent = math.pi / 8.0 * (grid_y**2 + grid_z**2)
                concentration = cc * np.exp(-exponent)
                axial_velocity = (
                    axial_ambient + uc * np.exp(-self.k.sc * exponent)
                )
                density_profile = self.rhoa + gamma * concentration
                enthalpy = np.array([
                    self.th.table.from_concentration(float(value)).enthalpy
                    for value in concentration.ravel()
                ]).reshape(concentration.shape)
                area_weights = (
                    math.pi / 4.0 * sysz * self.k.delta**2
                    * weight_y * weight_z
                )
                enthalpy_flux = float(np.sum(
                    density_profile
                    * (enthalpy - float(self.th.table.h[0]))
                    * axial_velocity * area_weights
                ))
            k3u = self.k.profile_integral(3.0 * self.k.sc)
            k_c3u = self.k.profile_integral(1.0 + 3.0 * self.k.sc)
            kinetic_integral = self.rhoa * sysz * (
                rk3 * axial_ambient**3
                + 3.0 * rk4 * axial_ambient**2 * uc
                + 3.0 * rk5 * axial_ambient * uc**2
                + k3u * uc**3
            ) + gamma * cc * sysz * (
                rk1 * axial_ambient**3
                + 3.0 * rk2 * axial_ambient**2 * uc
                + 3.0 * rk6 * axial_ambient * uc**2
                + k_c3u * uc**3
            )
            energy = float(enthalpy_flux + 0.5 * kinetic_integral)

        return JetIntegralFluxes(
            total_mass=float(total_mass),
            contaminant_mass=float(contaminant),
            momentum_x=float(momentum * ct),
            momentum_z=float(momentum * st),
            energy=energy,
        )

    # -- derivatives --------------------------------------------------------

    def derivatives(self, s: float, y: np.ndarray, dery: np.ndarray) -> JetState:
        """Port of ``MODEL``: the four coupled balances plus the trajectory."""
        k = self.k
        rk1, rk2, rk3, rk4, rk5, rk6 = self.rk
        cc, sysz, theta, uc, dist, zj = (
            y[J_CC], abs(y[J_SYSZ]), y[J_THETA], y[J_UC], y[J_X], y[J_Z]
        )
        uentr = max(uc, 0.0)
        st, ct = math.sin(theta), math.cos(theta)

        mix = self.th.table.from_concentration(cc)
        rho = mix.rho
        # Droplets still in the jet add mass without volume. Their fraction
        # falls as entrained air supplies the latent heat: one mass of liquid
        # needs `evaporation_ratio` masses of air, and the jet has (1/cc - 1)
        # masses of air per mass of contaminant.
        rho_bulk = rho
        if self.liquid_fraction > 0.0 and cc > 0.0:
            air_per_contaminant = max(1.0 / max(cc / rho, 1e-30) - 1.0, 0.0)
            evaporated = (
                air_per_contaminant / self.evaporation_ratio
                if self.evaporation_ratio > 0.0 else 1.0
            )
            still_liquid = self.liquid_fraction * max(1.0 - evaporated, 0.0)
            rho_bulk = rho + still_liquid * cc
        gamma = (rho_bulk - self.rhoa) / cc
        if self.non_boussinesq:
            gamma *= self.rhoa / max(rho_bulk, 1e-9)

        # ambient dispersion at this distance, and its growth rate
        sya = self.deltay * dist**self.betay
        sza = (
            self.deltaz * dist**self.betaz
            * math.exp(self.gammaz * math.log(dist) ** 2)
        )
        dsya = dsza = 0.0
        if dist > 1.0:
            dsya = sya / dist * self.betay * ct
            dsza = (
                sza / dist * (self.betaz + 2.0 * self.gammaz * math.log(dist)) * ct
            )

        sy, sz = self._split(sysz, sya, sza)
        # The cross-section reaches the ground when its lower edge does, at
        # delta sigma_z below the centreline. ELLIPS already knows how to cut
        # an ellipse at a plane, so the free perimeter follows from the same
        # routine that gives the full one; entrainment then acts only on the
        # part still exposed to the air.
        half_depth = k.delta * sz * ct
        ground_width = 0.0
        if self.ground_effect and zj < half_depth:
            segment = max(min(zj / max(ct, 1e-6), k.delta * sz), -k.delta * sz)
            pe, _area = ellipse(sy * k.delta, sz * k.delta, segment)
            ratio = segment / max(k.delta * sz, 1e-30)
            ground_width = (
                2.0 * k.delta * sy
                * math.sqrt(max(1.0 - ratio * ratio, 0.0))
            )
        else:
            pe, _area = ellipse(sy * k.delta, sz * k.delta, sz * k.delta)
        ua, alpha_wind = self._wind_profile(zj, sz, ct)

        ua2, uact, uast = ua * ua, ua * ct, ua * st
        if k.buoyant_entrainment > 0.0:
            # sqrt(g' H): the speed a parcel of this density deficit reaches
            # falling through its own depth. Independent of the wind, which
            # is the whole point.
            depth = 2.0 * k.delta * sz
            reduced_g = 9.81 * abs(rho_bulk - self.rhoa) / self.rhoa
            u_buoyant = math.sqrt(max(reduced_g * depth, 0.0))
            uentr = math.hypot(uentr, k.buoyant_entrainment * u_buoyant)
        if k.vertical_shear:
            # The relative velocity between plume and air is a vector. The
            # shipped model takes only its horizontal part, `uc`; the plume's
            # own vertical velocity is `(uc + ua cos theta) sin theta`, and
            # leaving it out is exactly what EFFECTS records correcting --
            # "the shear term contribution was extended with the vertical
            # velocity component which was neglected in dU in the original
            # model", because the original was aimed at heavy gases where a
            # plume does not rise.
            w_plume = (uc + uact) * st
            uentr = math.hypot(uentr, w_plume)
        uauc, uc2, ucct = ua * uc, uc * uc, uc * ct
        ccsysz = cc * sy * sz
        a = np.zeros((4, 5))

        # contaminant mass
        qqq = rk1 * uact + rk2 * uc
        rate = ccsysz * qqq
        a[0] = [qqq * sysz, qqq * cc, -rk1 * ccsysz * uast, rk2 * ccsysz, 0.0]

        # total mass, and the three entrainment terms
        qqq = rk3 * uact + rk4 * uc
        e3 = rk3 * ua * (sza * dsya + sya * dsza)
        e2 = k.alfa2 * uact * abs(st) * pe
        if k.momentum_entrainment_beta > 0.0:
            if self._source_momentum_entrainment is None:
                raise RuntimeError(
                    "source-momentum entrainment requires source conditions"
                )
            e1 = self._source_momentum_entrainment
        else:
            e1 = k.shear_coefficient(
                9.81 * abs(rho - self.rhoa) * math.sqrt(max(sysz, 0.0))
                / max(self.rhoa * max(uentr, 1e-6) ** 2, 1e-12),
                rho / self.rhoa,
            ) * uentr * pe
        e_ground = 0.0
        if k.ground_layer_entrainment and ground_width > 0.0:
            # DEGADIS's ground-layer top entrainment.  The centreline density
            # excess and the above-ground depth are both converted by delta;
            # equivalently this is the layer-average density over the full
            # depth.  Only a positive density excess is stably stratified.
            layer_depth = max(zj + half_depth, 1e-9)
            ri_layer = (
                9.81 * max(rho_bulk - self.rhoa, 0.0) * layer_depth
                / (
                    self.rhoa
                    * max(self.ustar, 1e-6) ** 2
                    * k.delta
                )
            )
            layer_velocity = (
                VKC * self.ustar * (1.0 + alpha_wind)
                / phi(ri_layer, 0.0, 3)
            )
            e_ground = layer_velocity * ground_width
        if k.houf_buoyant_entrainment:
            if self._houf_alpha_buoy is None:
                raise RuntimeError(
                    "Houf entrainment requires source conditions before integration"
                )
            # HyRAM uses v=v_cl exp(-r^2/B^2), whereas the DEGADIS widths
            # are Gaussian standard deviations. The geometric mean preserves
            # area for the elliptical cross-section.
            width = math.sqrt(max(2.0 * sysz, 1e-30))
            vcl = max(uc + uact, 1e-9)
            density_difference = abs(self.rhoa - rho_bulk)
            e_buoy = 0.0
            if density_difference > 1e-12 and st > 0.0:
                fr_local = (
                    vcl**2 * max(rho_bulk, 1e-12)
                    / (9.81 * width * density_difference)
                )
                e_buoy = (
                    self._houf_alpha_buoy / max(fr_local, 1e-30)
                    * 2.0 * PI * vcl * width * st
                )
            cap = k.HOUF_ALPHA_LIMIT * 2.0 * PI * width * vcl
            e1 = min(e1 + e_buoy, cap)
        a[1] = [
            0.0, self.rhoa * qqq, -self.rhoa * rk3 * sysz * uast,
            self.rhoa * rk4 * sysz,
            self.rhoa * (e1 + e2 + e3 + e_ground),
        ]

        drag = k.cd * pe * self.rhoa / 2.0 * uast**2

        # vertical momentum
        rrr = 2.0 * rk2 * uauc * ct + rk1 * ua2 * ct * ct + rk6 * uc2
        qqq = 2.0 * rk4 * uauc * ct + rk3 * ua2 * ct * ct + rk5 * uc2
        a[2] = [
            gamma * sysz * rrr * st,
            self.rhoa * qqq * st + gamma * cc * rrr * st,
            self.rhoa * sysz * (
                ct * qqq + st * (-2.0 * rk4 * uauc * st - 2.0 * rk3 * ua2 * ct * st)
            ) + gamma * ccsysz * (
                ct * rrr + st * (-2.0 * rk2 * uauc * st - 2.0 * rk1 * ua2 * ct * st)
            ),
            self.rhoa * sysz * st * (2.0 * rk4 * uact + 2.0 * rk5 * uc)
            + gamma * ccsysz * st * (2.0 * rk2 * uact + 2.0 * rk6 * uc),
            -rk1 * 9.81 * gamma * ccsysz - math.copysign(1.0, theta) * drag * ct,
        ]
        if k.shape_drag > 0.0:
            # Aspect ratio of the cross-section: vertical over horizontal.
            # Below one the plume is wider than it is tall and is bluff to
            # its own rise; at one it is circular.
            ar = min(max(sz / max(sy, 1e-9), 0.0), 1.0)
            cd = CD_CIRCULAR + (CD_PLATE - CD_CIRCULAR) * (1.0 - ar)
            w = uact * st
            width = 2.0 * k.delta * sy
            a[2][4] -= (
                k.shape_drag * cd * self.rhoa / 2.0 * w * abs(w) * width
            )
        if k.rise_drag > 0.0:
            # opposes the rise, quadratic in the vertical velocity, acting on
            # the plume's horizontal projection
            w = uact * st
            a[2][4] -= (
                k.rise_drag * self.rhoa / 2.0 * w * abs(w) * pe * abs(ct)
            )
        if k.liftoff_richardson > 0.0 and a[2][4] > 0.0:
            # Hold a grounded plume down until it is buoyant enough to lift.
            #
            # The transition is smoothed over the decade below the threshold
            # rather than switched, because Hall and Walker's wind-tunnel
            # measurements show lift-off "varying continuously with this
            # parameter, rather than occurring at a precise critical value".
            # A step would also stiffen the integration at the crossing.
            depth = max(k.delta * sz, 1e-9)
            if zj < depth:
                deficit = max(self.rhoa - rho_bulk, 0.0)
                ri = 9.81 * depth * deficit / (self.rhoa * max(self.ustar, 1e-6) ** 2)
                held = min(max(ri / k.liftoff_richardson, 0.0), 1.0)
                a[2][4] *= held

        if self.ground_effect and zj < half_depth and a[2][4] > 0.0:
            # A plume lying on the ground cannot lift its whole cross-section
            # at once: the part still in contact is held there. Scale the
            # buoyancy by the fraction that has cleared, which goes smoothly
            # to one as the plume rises free.
            a[2][4] *= min(max(zj / max(half_depth, 1e-9), 0.0), 1.0)

        # streamwise momentum
        rrr = rk1 * ua2 * ct**3 + 2.0 * rk2 * uauc * ct * ct + rk6 * uc2 * ct
        qqq = rk3 * ua2 * ct**3 + 2.0 * rk4 * uauc * ct * ct + rk5 * uc2 * ct
        a[3] = [
            gamma * sysz * rrr,
            self.rhoa * qqq + gamma * cc * rrr,
            self.rhoa * sysz * (
                -3.0 * rk3 * ua2 * ct * ct * st - 4.0 * rk4 * uauc * ct * st
                - rk5 * uc2 * st
            ) + gamma * ccsysz * (
                -3.0 * rk1 * ua2 * ct * ct * st - 4.0 * rk2 * uauc * ct * st
                - rk6 * uc2 * st
            ),
            self.rhoa * sysz * (2.0 * rk4 * ua * ct * ct + 2.0 * rk5 * ucct)
            + gamma * ccsysz * (2.0 * rk2 * ua * ct * ct + 2.0 * rk6 * ucct),
            ua * a[1][4] + drag * abs(st),
        ]

        # SIMUL solved this by Gauss elimination; the system is 4x4 and dense
        sol = np.linalg.solve(a[:, :4], a[:, 4])
        dery[0:4] = sol
        # the original clamps two derivatives that can go the wrong way
        # through round-off: concentration must fall, area must grow
        if dery[J_CC] > 0.0:
            dery[J_CC] = 0.0
        if dery[J_SYSZ] < 0.0:
            dery[J_SYSZ] = 0.0
        dery[J_X] = ct
        dery[J_Z] = st

        return JetState(
            s=s, x=dist, z=zj, cc=cc, sy=sy, sz=sz, theta=theta, uc=uc,
            ua=ua, rho=rho, temp=mix.temp, yc=mix.yc, rate=rate,
        )

    # -- initial conditions -------------------------------------------------

    @staticmethod
    def expanded_source_start(
        *, mass_flow: float, orifice_diameter: float, orifice_density: float,
        expanded_density: float, expanded_fraction: float,
        entrainment_momentum_conserving: bool = False,
        source_specific_momentum: float | None = None,
    ) -> tuple[float, float, float]:
        """Density, contaminant fraction and diameter at the expanded plane.

        The jet leaves the orifice as a droplet-laden two-phase fluid and
        expands, entraining the air its remaining liquid needs to evaporate.
        The plume model should start where that has happened, not at the
        orifice: the state there is what the flash calculation computes, and
        the area follows from continuity.

        The historical candidate scales velocity by
        ``sqrt(rho_orifice / rho_expanded)``.  When
        ``entrainment_momentum_conserving`` is true, the expanded plane is
        treated as a short entrainment/heating zone whose surrounding air
        enters at rest: ``m_H2 u_in = m_total u_out``, hence
        ``u_out = expanded_fraction * u_in``.  ``source_specific_momentum``
        may supply ``u_in`` from a preceding measured-pressure notional
        nozzle; otherwise it follows from orifice continuity as before.
        Pressure thrust belongs to that preceding calculation and is not
        created again by air entrainment.

        Returns ``(density, concentration, diameter)``.
        """
        area = PI * orifice_diameter**2 / 4.0
        if source_specific_momentum is not None and source_specific_momentum <= 0.0:
            raise ValueError("source specific momentum must be positive")
        u_orifice = (
            mass_flow / max(orifice_density * area, 1e-30)
            if source_specific_momentum is None
            else float(source_specific_momentum)
        )
        if entrainment_momentum_conserving:
            u_expanded = expanded_fraction * u_orifice
        else:
            u_expanded = u_orifice * math.sqrt(
                orifice_density / max(expanded_density, 1e-30)
            )
        expanded_area = mass_flow / max(
            expanded_fraction * expanded_density * u_expanded, 1e-30
        )
        return (
            expanded_density,
            expanded_fraction,
            math.sqrt(4.0 * expanded_area / PI),
        )

    def initial_conditions_directed(
        self, *, erate: float, diajet: float, elejet: float, ua: float,
        theta0: float, rho_exit: float | None = None,
        concentration: float = 1.0,
    ) -> np.ndarray:
        """Start a jet released in a given direction.

        ``SETJET`` places the plume using the Kamotani and Greber trajectory
        for a jet issuing *vertically* into a crossflow, which is what
        DEGADIS's test cases are.  A horizontal release does not follow that
        trajectory, so the correlation cannot be used to find the starting
        point.

        What survives is the other half of ``SETJET``: the zone of flow
        development still has to be skipped, and its length is a property of
        the jet rather than of its direction.  Pratte and Baines' development
        length is kept and the plume is started that far along its own axis,
        at the release angle.

        Parameters
        ----------
        theta0
            Release angle above horizontal, radians.  ``0`` for a horizontal
            jet, ``+pi/2`` vertically up, ``-pi/2`` vertically down.
        rho_exit
            Density of the material entering the plume, kg/m**3.  For a
            flashing release this is the density after the flash, not the
            stored liquid.
        concentration
            Contaminant mass fraction entering the plume.  Below one for a
            release that has already entrained air while flashing.

        Notes
        -----
        The three must describe the *same* plane. Taking the density from
        after the flash while keeping the concentration and the area at the
        orifice is the error EPA's 1991 evaluation records the SLAB developer
        objecting to -- "for jet releases the source area should be the
        cross-section of the fully expanded jet rather than the orifice" --
        and it was made here before it was noticed. Correcting it moved the
        modelled plume rise on the PRESLHY trials from 1.07 m to 0.60 m
        against a measured -0.12, without touching a coefficient.

        :func:`expanded_source_start` builds a consistent triple.
        """
        rho = self.rhoe if rho_exit is None else rho_exit
        qinit = erate / concentration / rho
        uj = qinit / PI / diajet**2 * 4.0

        density_difference = abs(self.rhoa - rho)
        if density_difference <= 1e-15:
            source_froude = math.inf
        else:
            source_froude = uj / math.sqrt(
                9.81 * diajet * density_difference / max(rho, 1e-30)
            )
        self.houf_source_froude = source_froude
        self._houf_alpha_buoy = self.k.houf_buoyancy_coefficient(source_froude)
        if self.k.momentum_entrainment_beta > 0.0:
            area = PI * diajet**2 / 4.0
            self._source_momentum_entrainment = (
                self.k.momentum_entrainment_beta
                * math.sqrt(area * rho * uj**2 / self.rhoa)
            )

        # Pratte and Baines: the development length, in orifice diameters
        if ua <= 0.0:
            # Limit of Pratte and Baines' development-length correlation as
            # the co/cross-flow speed tends to zero.
            sod = 7.7
        else:
            sod = 7.7 * (
                1.0 - math.exp(-0.48 * math.sqrt(rho * uj / self.rhoa / ua))
            )
        distance = sod * diajet
        x = distance * math.cos(theta0)
        z = elejet + distance * math.sin(theta0)
        if z <= 0.0:
            # a downward jet that reaches the ground inside its own
            # development length has no established plume to model
            raise ValueError(
                f"the jet reaches the ground within its development length "
                f"({distance:.2f} m from a {elejet:.2f} m release)"
            )

        uc = uj - ua * math.cos(theta0)
        rk1, rk2 = self.rk[0], self.rk[1]
        cc = concentration * rho
        sysz = erate / cc / (rk1 * ua * math.cos(theta0) + rk2 * max(uc, 1e-6))

        y = np.zeros(6)
        y[J_CC] = cc
        y[J_SYSZ] = sysz
        y[J_THETA] = theta0
        y[J_UC] = uc
        y[J_X] = max(x, 1.0e-30)
        y[J_Z] = z
        return y

    def initial_conditions(
        self, *, erate: float, diajet: float, elejet: float, ua: float
    ) -> np.ndarray:
        """Port of ``SETJET``.

        The jet leaves the orifice with a top-hat velocity profile, which the
        similarity model cannot start from. ``SETJET`` skips the zone of flow
        development entirely, using the Kamotani and Greber (NASA CR-2392)
        wind-tunnel correlations to place the plume where the profile has
        become Gaussian: a trajectory ``z/D = a (x/D)^b`` whose coefficients
        depend on the momentum ratio ``J = (rho_e/rho_a)(u_j/u_a)^2`` and, in
        the mid range, on the Froude number. The development length follows
        Pratte and Baines.

        The starting point is where that trajectory reaches an arc length
        ``S/D``, found by Newton-Raphson.
        """
        qinit = erate / self.rhoe
        uj = qinit / PI / diajet**2 * 4.0
        cc = self.rhoe

        if self.k.momentum_entrainment_beta > 0.0:
            area = PI * diajet**2 / 4.0
            self._source_momentum_entrainment = (
                self.k.momentum_entrainment_beta
                * math.sqrt(area * self.rhoe * uj**2 / self.rhoa)
            )

        sod = 7.7 * (1.0 - math.exp(-0.48 * math.sqrt(self.rhoe * uj / self.rhoa / ua)))
        rj = (self.rhoe / self.rhoa) * (uj / ua) ** 2
        delrho = self.rhoe - self.rhoa
        froude = 1.688
        if delrho > 0.0:
            froude = self.rhoa * ua**2 / 9.81 / delrho / diajet
        froude = min(froude, 1.688)

        lr = math.log(rj) if rj > 0 else -50.0
        if rj <= 0.036:
            aj, bj = 18.519 * rj, 0.4
        elif rj <= 10.0:
            aj = math.exp(0.2476 + 0.3016 * math.log(froude) + 0.24386 * lr)
            bj = 0.4
        elif rj <= 50.0:
            aj = math.exp(0.405465 + 0.131386 * lr + 0.054931 * lr * lr)
            bj = math.exp(-0.744691 - 0.074525 * lr)
        elif rj <= 600.0:
            aj = math.exp(-2.55104 + 1.49202 * lr - 0.097623 * lr * lr)
            bj = math.exp(-0.446718 - 0.150694 * lr)
        else:
            # beyond the correlation's range; the original extrapolates and says so
            aj = math.exp(1.44099 + 0.243045 * lr)
            bj = math.exp(-0.446718 - 0.150694 * lr)

        bji = 1.0 / bj
        zod = sod
        for _ in range(100):
            xod = (zod / aj) ** bji
            fff = xod**2 + zod**2 - sod**2
            fffp = xod**2 * 2.0 / bj / zod + 2.0 * zod
            zodn = zod - fff / fffp
            if abs((zodn - zod) / zod) <= 1.0e-5:
                break
            zod = zodn
        else:
            raise RuntimeError("SETJET trajectory iteration did not converge")

        xod = (zod / aj) ** bji
        dist = xod * diajet
        if xod > 0.0:
            theta = math.atan(bj * zod / xod)
        elif xod == 0.0:
            theta = PI / 2.0
        else:
            raise ValueError("SETJET produced a negative downwind distance")

        zj = zod * diajet + elejet
        uc = uj - ua * math.cos(theta)
        rk1, rk2 = self.rk[0], self.rk[1]
        sysz = erate / cc / (rk1 * ua * math.cos(theta) + rk2 * uc)

        y = np.zeros(6)
        y[J_CC] = cc
        y[J_SYSZ] = sysz
        y[J_THETA] = theta
        y[J_UC] = uc
        y[J_X] = max(dist, 1.0e-30)
        y[J_Z] = zj
        return y

    # -- the driver ---------------------------------------------------------

    def run(
        self, y0: np.ndarray, *, distmx: float, tol: float = 1.0e-4,
        smax: float = 1.0e5,
    ) -> JetResult:
        """Integrate the plume until it lands or thins out.

        Port of ``MODOUT``'s two termination branches. On touchdown the state
        is interpolated linearly in elevation back to ``z = 0`` and the
        concentration doubled, because the ground reflects the half of the
        plume that would have gone below it.
        """
        rows: list[list[float]] = []
        prev: list[float] | None = None
        result = JetResult()
        holder: dict = {}

        def fct(s, yy, dd, prmt):
            holder["st"] = self.derivatives(s, yy, dd)

        def outp(s, yy, dd, ihlf, ndim, prmt):
            nonlocal prev
            st = holder["st"]
            cur = [
                st.x, st.z, st.cc, st.sy, st.sz, st.theta, st.uc, st.ua,
                st.rho, st.temp, st.yc, s,
            ]

            if st.z <= 0.0 and prev is not None:  # noqa: SIM102
                # interpolate back to the ground
                f = st.z / (prev[1] - st.z)
                land = [c - f * (p - c) for p, c in zip(prev, cur)]
                cc = 2.0 * land[2]  # image reflection doubles the centreline
                mix = self.th.table.from_concentration(cc)
                land[1] = 0.0
                land[2] = cc
                land[10] = mix.yc
                rows.append(land)
                result.touchdown = True
                result.distance = land[0]
                result.concentration = cc
                result.halfwidth = self.k.delta * land[3]
                prmt.halt()
                return

            # MODOUT stores the interpolation values *before* adding the
            # ground image, then reports the imaged concentration. So the
            # touchdown interpolation uses the bare plume while the reported
            # centreline includes the reflection.
            prev = list(cur)
            cimage = st.cc * math.exp(-0.5 * (2.0 * st.z / st.sz) ** 2)
            imaged = self.th.table.from_concentration(st.cc + cimage)
            cur[2] = st.cc + cimage
            cur[8] = imaged.rho
            cur[9] = imaged.temp
            cur[10] = imaged.yc
            rows.append(cur)
            ym = imaged.yc
            if ym < self.yclow:
                # never reaches the ground at a concentration of interest;
                # DIST = 0 is JETPLU's signal that DEGADIS must not continue
                result.touchdown = False
                result.distance = 0.0
                prmt.halt()

        # For a vertical free jet ``x`` is identically zero, so it cannot set
        # an integration length scale.  Use the initial Gaussian width only
        # in that limit; all historical cross-wind runs retain the original
        # ``x/20`` starting step exactly.
        first_step_scale = y0[J_X]
        if abs(first_step_scale) <= 1.0e-20:
            first_step_scale = math.sqrt(abs(y0[J_SYSZ]))
        prmt = Control([
            0.0, smax, max(first_step_scale / 20.0, 1.0e-30), tol, distmx
        ] + [0.0] * 20)
        rkgst(fct, outp, prmt, y0, [1.0] * 6, ndim=6)
        result.rows = np.array(rows) if rows else np.zeros((0, 12))
        return result
