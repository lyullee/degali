"""Lift-off of a ground-based buoyant plume.

DEGADIS cannot represent this.  Its ground-level model has no vertical
momentum equation at all: gravity enters only as lateral slumping, and only
while the cloud is *denser* than air.  Once a cloud becomes lighter than
ambient, ``PSS`` sets ``dB/dx = 0`` and the cloud simply stops spreading and
disperses passively, pinned to the ground for ever.

That gap is the whole problem for liquid hydrogen.  Cold hydrogen vapour is
dense at 20 K, so a release starts as a dense cloud, but hydrogen's molecular
weight is 2 and the cloud becomes violently buoyant as it warms and dilutes.
A model that cannot lift the cloud off the ground will keep a flammable layer
at head height that in reality has risen away.

The same is true of the release this theory was developed for: anhydrous
hydrogen fluoride, which is heavier than air on release and becomes buoyant as
it reacts with atmospheric moisture.

Where this comes from
---------------------
The integral model of Slawson and co-workers, extended to a ground-truncated
cross-section by AEA Technology under the EC URAHFREP project (report
AEAT/NOIL/27328006/001, June 2001).  The equations are implemented here from
that description; no code was taken from anywhere.

The plume cross-section is a *lozenge*: a rectangle of fixed length
:math:`L_s` with semicircular ends of radius :math:`R`, truncated at ground
level.  That one shape covers a circular source (:math:`L_s = 0`), a line
source (:math:`L_s \\gg R`), and the wide shallow cloud a dense release leaves
behind when it turns buoyant -- which is exactly the shape a lift-off model
has to start from and the reason a simple round plume will not do.

Six quantities are integrated along the plume axis :math:`s`:

=================  =====================================================
:math:`m_g`        contaminant mass flux, conserved
:math:`m`          total mass flux, grows by entrainment
:math:`M_x`        horizontal momentum flux
:math:`M_z`        vertical momentum flux, driven by buoyancy
:math:`x, z`       trajectory
=================  =====================================================

with entrainment through the *free* perimeter only -- the part not in contact
with the ground:

.. math::
    E = \\rho_a L_{free}
        \\left[\\alpha |u - u_a\\cos\\theta| + \\beta |u_a \\sin\\theta|
        + \\gamma u_a\\right]

The three terms are shear along the plume, cross-flow across it, and ambient
turbulence.  The report's own assessment of this model is that it is "very
crude", particularly in representing ambient turbulence by :math:`\\gamma u_a`
and in neglecting velocity shear and ground effects.  It is kept as stated
rather than improved, because its virtue is that it was compared against wind
tunnel data in that form.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..core.constants import GG, PI

#: Entrainment coefficients from the URAHFREP report, taken in turn from the
#: free-plume models it builds on.
ALPHA = 0.1  #: shear along the plume axis
BETA = 0.6  #: cross-flow, normal to the axis
GAMMA = 0.1  #: ambient turbulence

#: Form-drag coefficient on the rising cross-section.
#:
#: The URAHFREP equations have no drag term.  Mack and co-workers (*Process
#: Safety and Environmental Protection* **176**, 2023) report the consequence
#: from their own work on EFFECTS: lift-off itself is predicted well, but
#: "during the rising phase of strongly buoyant plumes buoyancy was initially
#: over predicted resulting in too high plume trajectories", and "applying the
#: additional pressure drag then shows a good prediction of the plume
#: trajectory during rise".  The drag is "mainly oriented in vertical
#: direction for strongly buoyant plumes with high plume rise velocity
#: compared to ambient wind speed" -- which is exactly the regime an LH2 cloud
#: reaches.
#:
#: The value matches the form drag ``JETPLU`` already applies to its own
#: inclined cross-section.
CD = 0.2

#: Critical values of the bulk Richardson number, from Hall and Walker's wind
#: tunnel measurements as interpreted in section 4.1 of the report.  Briggs'
#: earlier estimates agree.  ``Ri*`` here is formed on the buoyancy flux, so
#: it is *positive* for a buoyant plume -- the opposite sign convention from
#: the dense-phase ``Ri*`` in :mod:`degali.core.atmosphere`.
#:
#: Hall and Walker state the same two thresholds in their own variable,
#: :math:`F/(W u^3)`: about 0.01 where the concentration maximum leaves the
#: ground and about 0.035 where the ground-level value has fallen to 10-20
#: per cent of the maximum, both measured 15 to 30 source lengths downwind
#: and with "considerable scatter". Fixing the constant at 200 from the first
#: of those maps the second to ``Ri* = 7``, against the 10 taken here from the
#: report's own reduction. Two readings of the same experiments agreeing to
#: thirty per cent is about what "considerable scatter" allows, and it is an
#: independent check rather than a fit.
RI_ONSET = 2.0  #: the concentration maximum leaves the ground
RI_SUBSTANTIAL = 10.0  #: ground-level concentration down to 10-20 % of maximum

#: Hall and Walker's own thresholds, in their variable ``F/(W u**3)``.
HW_FIRST_RISE = 0.01
HW_LIFTOFF = 0.035
RI_CLEAR = 70.0  #: ground-level concentration below 5 % of maximum

#: State vector indices.
I_MASS, I_MX, I_MZ, I_X, I_Z = range(5)


@dataclass
class LiftoffState:
    """The plume at one point along its trajectory."""

    s: float  #: arc length, m
    x: float
    z: float  #: centroid height, m
    theta: float  #: trajectory angle above horizontal, rad
    u: float  #: plume velocity, m/s
    radius: float  #: semicircular end radius, m
    area: float  #: cross-sectional area above ground, m**2
    concentration: float  #: contaminant mass fraction, averaged over the section
    density: float  #: kg/m**3
    ground_contact: float  #: length of cross-section touching the ground, m
    segment_length: float = 0.0  #: the lozenge's flat part, m
    delta: float = 2.15  #: section boundary, in sigma

    @property
    def sigma(self) -> float:
        """Gaussian width of the concentration profile, m."""
        return self.radius / self.delta

    @property
    def shape_factor(self) -> float:
        r"""Cross-sectional mean over peak, for a truncated Gaussian profile.

        The integral model conserves the contaminant flux, so what it carries
        is the *mean* concentration over the section.  A grab bottle, or any
        point measurement near the axis, sees the *peak*.  Comparing one
        against the other under-predicts by this factor, and for a circular
        section cut at 2.15 sigma it is 0.39 -- a factor of 2.6, which is
        very nearly the discrepancy measured against the NASA sample bottles
        before the profile was introduced.

        Integrating :math:`\exp(-r^2/2\sigma^2)` over the lozenge: the flat
        part contributes a one-dimensional Gaussian across its width, the two
        semicircular ends a two-dimensional one.
        """
        sigma = self.sigma
        if sigma <= 0.0:
            return 1.0
        d = self.delta
        rect = (
            self.segment_length * sigma * math.sqrt(2.0 * PI)
            * math.erf(d / math.sqrt(2.0))
        )
        ends = 2.0 * PI * sigma**2 * (1.0 - math.exp(-0.5 * d * d))
        area = 2.0 * self.segment_length * self.radius + PI * self.radius**2
        return (rect + ends) / area if area > 0.0 else 1.0

    @property
    def peak_concentration(self) -> float:
        """Contaminant mass fraction on the plume axis."""
        f = self.shape_factor
        return min(self.concentration / f, 1.0) if f > 0.0 else self.concentration

    def at_height(self, height: float) -> float:
        """Mass fraction at an elevation, on the plume's vertical centre line.

        Zero outside the section, which is where the model stops claiming to
        know anything.
        """
        offset = abs(height - self.z)
        if offset > self.delta * self.sigma:
            return 0.0
        return self.peak_concentration * math.exp(
            -0.5 * (offset / max(self.sigma, 1e-12)) ** 2
        )

    @property
    def airborne(self) -> bool:
        """Whether the cross-section has cleared the ground."""
        return self.ground_contact <= 0.0


@dataclass
class LiftoffResult:
    """A lift-off trajectory."""

    states: list[LiftoffState] = field(default_factory=list)
    #: Distance at which the cross-section leaves ground it was in contact
    #: with, m.  ``nan`` if it never touches down, or never leaves.  A plume
    #: released already clear of the ground has not lifted off.
    liftoff_distance: float = float("nan")
    reason: str = ""

    @property
    def lifts_off(self) -> bool:
        return math.isfinite(self.liftoff_distance)

    def as_array(self) -> np.ndarray:
        """``(n, 6)``: x, z, theta, u, concentration, density."""
        return np.array([
            [s.x, s.z, s.theta, s.u, s.concentration, s.density]
            for s in self.states
        ])

    def height_at(self, x: float) -> float:
        a = self.as_array()
        return float(np.interp(x, a[:, 0], a[:, 1]))

    def concentration_at(self, x: float) -> float:
        a = self.as_array()
        return float(np.interp(x, a[:, 0], a[:, 4]))


def richardson_liftoff(
    buoyancy_flux: float, wind: float, width: float
) -> float:
    r"""Bulk Richardson number governing lift-off.

    Formed on the buoyancy flux per unit width, :math:`Ri^* \propto
    F_W / u_a^3`, so it is positive for a buoyant plume.  The report's
    equivalence is ``F/u_a^3 W ~ 0.01`` at ``Ri* ~ 2``, which fixes the
    constant at 200.

    Parameters
    ----------
    buoyancy_flux
        :math:`F = g Q (\rho_a - \rho)/\rho_a`, m**4/s**3, with ``Q`` the
        volume flux.
    wind
        Ambient speed at the plume, m/s.
    width
        Source width, m.
    """
    if wind <= 0.0 or width <= 0.0:
        return float("inf")
    return 200.0 * buoyancy_flux / (wind**3 * width)


def hall_walker_parameter(
    buoyancy_flux: float, wind: float, width: float
) -> float:
    """Hall and Walker's own lift-off variable, :math:`F/(W u^3)`.

    The same physics as :func:`richardson_liftoff` without the factor of 200,
    so a result can be quoted against their published thresholds directly:
    :data:`HW_FIRST_RISE` and :data:`HW_LIFTOFF`.
    """
    if wind <= 0.0 or width <= 0.0:
        return float("inf")
    return buoyancy_flux / (wind**3 * width)


def liftoff_regime(ri: float) -> str:
    """Describe what Hall and Walker measured at this Richardson number."""
    if ri < RI_ONSET:
        return "no lift-off: the concentration maximum stays on the ground"
    if ri < RI_SUBSTANTIAL:
        return "lift-off beginning: the maximum has left the ground"
    if ri < RI_CLEAR:
        return "ground-level concentration down to 10-20 % of the maximum"
    return "plume clear of the ground: below 5 % of the maximum"


class LiftoffPlume:
    """The URAHFREP ground-truncated buoyant plume model.

    Parameters
    ----------
    rho_ambient
        kg/m**3.
    wind
        Ambient speed, m/s.  Taken as uniform: the model neglects shear, and
        the report says so.
    segment_length
        :math:`L_s`, the flat part of the lozenge.  Zero gives a round plume;
        large values give a line source.  For a cloud handed over from a dense
        phase this is the width the dense phase left it with.
    density_of
        ``concentration -> density``.  Supplying
        :meth:`degali.core.thermo.AdiabaticTable.from_mass_fraction` couples
        this to the same mixing line the rest of the model uses, so a cloud
        that becomes buoyant does so for thermodynamic reasons rather than by
        assumption.
    """

    def __init__(
        self,
        *,
        rho_ambient: float,
        wind: float,
        segment_length: float = 0.0,
        density_of=None,
        alpha: float = ALPHA,
        beta: float = BETA,
        gamma: float = GAMMA,
        drag: float = CD,
    ):
        self.rhoa = rho_ambient
        self.ua = wind
        self.ls = segment_length
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self.drag = drag
        self._density_of = density_of

    # -- geometry ----------------------------------------------------------

    def _geometry(self, area: float, z: float, theta: float) -> tuple[float, float, float]:
        """Solve the lozenge for ``(R, ground contact, free perimeter)``.

        The cross-section area is known from the fluxes; the radius that
        produces it, truncated at the ground, is not, so it is found by a
        bracketed search.  ``z1`` is the centre height measured in the plane
        of the cross-section, which is ``z / cos(theta)``.
        """
        z1 = z / max(math.cos(theta), 1e-6)

        def area_of(r: float) -> float:
            if r <= 0.0:
                return 0.0
            if z1 >= r:  # clear of the ground: the full lozenge
                return self.ls * 2.0 * r + PI * r * r
            if z1 <= -r:  # entirely below ground, which cannot happen
                return 0.0
            # circle of radius r centred z1 above the ground, area above it
            cap = r * r * (PI - math.acos(z1 / r)) + z1 * math.sqrt(
                max(r * r - z1 * z1, 0.0)
            )
            return self.ls * (r + z1) + cap

        lo, hi = 1e-9, max(10.0, abs(z1) * 10.0 + 1.0)
        while area_of(hi) < area and hi < 1e6:
            hi *= 2.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if area_of(mid) < area:
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)

        if z1 >= r:
            ground = 0.0
            free = self.ls * 2.0 + 2.0 * PI * r
        else:
            ground = self.ls + 2.0 * math.sqrt(max(r * r - z1 * z1, 0.0))
            free = self.ls + 2.0 * r * (PI - math.acos(max(min(z1 / r, 1.0), -1.0)))
        return r, ground, free

    # -- closure -----------------------------------------------------------

    #: Points used to integrate the density over the cross-section.
    PROFILE_POINTS = 33

    def _density(self, concentration: float) -> float:
        if self._density_of is None:
            return self.rhoa
        return float(self._density_of(concentration))

    def _section_density(self, mean: float, radius: float) -> float:
        r"""Area-averaged density over the cross-section, kg/m**3.

        The buoyancy force is :math:`g(\rho_a - \bar\rho)A` with
        :math:`\bar\rho` the *mean density*, and the mean density is not the
        density at the mean concentration.  On a mixing line as curved as
        hydrogen's the two differ substantially, and taking the second for the
        first is the same error as evaluating any non-linear function at an
        average argument.

        It matters here for a specific reason.  Giannissi and co-workers'
        CFD study of LH2 releases finds that condensation and freezing of
        *atmospheric humidity* dominates the cloud's buoyancy, and that the
        region where it happens is spatially extended -- a shell around the
        cold core where cold hydrogen meets moist air.  An integral model that
        evaluates the thermodynamics once, at the section mean, cannot see
        that shell.  Integrating the same mixing line over the concentration
        profile does, because the curvature the condensation puts into
        :math:`\rho(c)` is exactly what the averaging picks up.

        The profile is the truncated Gaussian of :class:`LiftoffState`, and
        the integral is over the lozenge: a one-dimensional average across the
        flat part, two-dimensional over the semicircular ends.
        """
        if self._density_of is None or radius <= 0.0:
            return self._density(mean)
        state = LiftoffState(
            s=0.0, x=0.0, z=0.0, theta=0.0, u=0.0, radius=radius, area=0.0,
            concentration=mean, density=0.0, ground_contact=0.0,
            segment_length=self.ls,
        )
        peak = state.peak_concentration
        sigma = state.sigma
        d = state.delta
        n = self.PROFILE_POINTS
        r = np.linspace(0.0, d * sigma, n)
        c = peak * np.exp(-0.5 * (r / max(sigma, 1e-12)) ** 2)
        rho = np.array([self._density(float(v)) for v in c])
        # weights: the flat part contributes a strip of length Ls at each |r|,
        # the ends an annulus of circumference 2 pi r
        w = 2.0 * self.ls + 2.0 * PI * r
        num = np.trapezoid(rho * w, r)
        den = np.trapezoid(w, r)
        return float(num / den) if den > 0.0 else self._density(mean)

    def _unpack(self, y: np.ndarray, mg: float) -> LiftoffState:
        mass, mx, mz, x, z = y
        theta = math.atan2(mz, mx)
        u = math.hypot(mx, mz) / max(mass, 1e-30)
        c = mg / max(mass, 1e-30)
        rho = self._density(c)
        area = mass / max(rho * u, 1e-30)
        r, ground, _free = self._geometry(area, z, theta)
        return LiftoffState(
            s=0.0, x=x, z=z, theta=theta, u=u, radius=r, area=area,
            concentration=c, density=rho, ground_contact=ground,
            segment_length=self.ls,
        )

    def _derivatives(self, s: float, y: np.ndarray, mg: float) -> np.ndarray:
        mass, mx, mz, _x, z = y
        theta = math.atan2(mz, mx)
        u = math.hypot(mx, mz) / max(mass, 1e-30)
        mean_c = mg / max(mass, 1e-30)
        rho = self._density(mean_c)
        area = mass / max(rho * u, 1e-30)
        _r, _ground, free = self._geometry(area, z, theta)
        # buoyancy uses the section-averaged density, not the density at the
        # section-averaged concentration
        rho_buoy = self._section_density(mean_c, _r)

        entrainment = self.rhoa * free * (
            self.alpha * abs(u - self.ua * math.cos(theta))
            + self.beta * abs(self.ua * math.sin(theta))
            + self.gamma * self.ua
        )
        d = np.zeros(5)
        d[I_MASS] = entrainment
        d[I_MX] = entrainment * self.ua
        # the only place buoyancy enters, and the equation DEGADIS lacks
        buoyancy = GG * (self.rhoa - rho_buoy) * area
        # form drag on the cross-section as it rises: the plume presents a
        # width Ls + 2R to its own vertical motion, and the resistance goes as
        # the square of the rise velocity. Without it a strongly buoyant plume
        # accelerates on buoyancy alone.
        w = mz / max(mass, 1e-30)
        width = self.ls + 2.0 * _r
        d[I_MZ] = buoyancy - self.drag * self.rhoa / 2.0 * w * abs(w) * width
        d[I_X] = math.cos(theta)
        d[I_Z] = math.sin(theta)
        return d

    # -- driver ------------------------------------------------------------

    def run(
        self,
        *,
        rate: float,
        concentration: float,
        velocity: float,
        height: float,
        theta: float = 0.0,
        max_distance: float = 1000.0,
        min_concentration: float = 1e-6,
        max_height: float = 500.0,
    ) -> LiftoffResult:
        """Integrate a plume from a ground-level source.

        Parameters
        ----------
        rate
            Contaminant mass rate, kg/s.
        concentration
            Contaminant mass fraction at the source.
        velocity
            Plume speed at the source, m/s.  For a cloud handed over from a
            dense phase this is the wind speed at its effective height.
        height
            Centroid height at the source, m.  Zero for a cloud sitting on the
            ground; the model then has the cross-section half buried, which is
            the truncation the geometry handles.
        """
        mg = rate
        mass = rate / max(concentration, 1e-12)
        y0 = np.array([
            mass,
            mass * velocity * math.cos(theta),
            mass * velocity * math.sin(theta),
            0.0,
            height,
        ])

        result = LiftoffResult()
        lifted = {"s": float("nan")}

        def stop_far(s, y, mg):
            return y[I_X] - max_distance
        def stop_high(s, y, mg):
            return y[I_Z] - max_height
        def stop_thin(s, y, mg):
            return mg / max(y[I_MASS], 1e-30) - min_concentration
        for e in (stop_far, stop_high, stop_thin):
            e.terminal = True
        stop_far.direction = 1.0
        stop_high.direction = 1.0
        stop_thin.direction = -1.0

        sol = solve_ivp(
            self._derivatives, (0.0, 10.0 * max_distance), y0,
            args=(mg,), events=(stop_far, stop_high, stop_thin),
            rtol=1e-6, atol=1e-9, dense_output=True, max_step=max_distance / 50.0,
        )
        if not sol.success:
            result.reason = f"integration failed: {sol.message}"
            return result

        touched = False
        for k in range(sol.y.shape[1]):
            st = self._unpack(sol.y[:, k], mg)
            st.s = float(sol.t[k])
            result.states.append(st)
            if not st.airborne:
                touched = True
            elif touched and not math.isfinite(lifted["s"]):
                # lift-off means leaving ground it was on, so a plume released
                # already clear of the ground has not lifted off
                lifted["s"] = st.x

        result.liftoff_distance = lifted["s"]
        if result.lifts_off:
            result.reason = f"cross-section cleared the ground at {lifted['s']:.1f} m"
        elif result.states and result.states[0].airborne:
            result.reason = "released clear of the ground; never in contact"
        elif sol.t_events[1].size:
            result.reason = "reached the height limit while still touching down"
        elif sol.t_events[2].size:
            result.reason = "diluted below the concentration floor without lifting off"
        else:
            result.reason = "reached the distance limit without lifting off"
        return result
