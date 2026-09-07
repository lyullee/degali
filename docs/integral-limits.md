# Applying an integral model to liquid hydrogen: the structural limits

Not a data-fitting exercise. This is what stops an integral dense-gas model
from describing an LH₂ release, which of those limits have been removed here,
and which are properties of the class of model rather than of any particular
one.

## 1. The dense-to-buoyant transition was a seam — removed

An LH₂ release passes continuously through neutral buoyancy: dense and
slumping near the source, buoyant and rising a few metres later. Nothing in
the physics switches there. The available models did:

| model | has | lacks |
|---|---|---|
| `DegadisClosure` | gravity slumping | any vertical momentum |
| `BuoyantClosure` | vertical momentum | perimeter entrainment; rise runs away |
| `LiftoffPlume` | perimeter entrainment | ground-layer entrainment |

Running them in sequence needs a handover height, which is **a free parameter
standing in for physics**, and at which neither model is valid — the
ground-layer entrainment is out of calibration above it and the perimeter law
does not apply below.

`addons.unified.UnifiedClosure` removes it. The cross-section always has a
**ground-contact fraction**

    f = L_ground / (L_ground + L_free)

which the geometry already computes: 1 for a cloud lying flat, 0 once it has
cleared. Entrainment is then

    E = f · E_layer + (1 − f) · E_perimeter

so a ground-hugging cloud entrains through its top as DEGADIS says, a risen
plume through its whole perimeter as URAHFREP says, and a cloud in between
does both in proportion to where it is. Buoyancy is scaled by the same
fraction, so a cloud resting on the ground cannot lift the part still on it.

Nothing is fitted. Both entrainment laws are the ones already validated in
their own regimes; the blend is geometry. One integration now carries a cloud
from `rho/rho_a = 1.30` on the ground to twenty metres up, with the regime
label moving from "on the ground" through "lifting" to "airborne" without a
threshold anywhere.

## 2. Buoyancy from the mean state — partly removed

The buoyancy force is `g(rho_a − rho_bar)A` with `rho_bar` the *mean density*,
and on hydrogen's mixing line that is not the density at the mean
concentration. `LiftoffPlume._section_density` now integrates the mixing line
over the concentration profile.

That matters for a specific reason. Giannissi and co-workers' CFD of these
releases finds that **condensation and freezing of atmospheric humidity
dominates LH₂ cloud buoyancy**, in a spatially extended shell where cold
hydrogen meets moist air. Evaluating the thermodynamics once at the section
mean cannot see that shell; integrating over the profile is closer.

It is only closer. Every annulus is still assumed to be at adiabatic
equilibrium, and a shell reaction is not that. Removing this properly means
carrying a condensation rate rather than an equilibrium state — a second
conserved quantity, not a reinterpretation of the first.

## 3. Cross-section shape — removed

The lozenge was a top hat: uniform concentration across a section whose radius
grows to tens of metres. A point measurement sees the peak. The peak-to-mean
ratio for a section cut at 2.15 sigma is 2.57, and the discrepancy measured
against the NASA grab bottles with a top-hat section was 2.45. The section now
carries a Gaussian profile.

## 4. Form drag — removed

The URAHFREP equations have none. Mack and co-workers report that integral
models over-predict the rise of strongly buoyant plumes and that pressure drag
fixes the trajectory. `JETPLU` already applied exactly that to its own
inclined section; the same term is now in the buoyant plume. It changes these
cases by 2 %, because they rise at about a metre per second and the resistance
is quadratic — but it bounds the behaviour of a faster one.

## What remains, and why it is not a coefficient

**Air condensation.** Below 90 K oxygen condenses and below 77 K nitrogen
does, removing mass from the gas phase and releasing latent heat. On the LH₂
mixing line that is above 81 mole per cent — near the source, exactly where
the dense phase lives. Modelling it means a mixing line with a second
condensing species and a composition that changes as it goes, not a fixed
line.

**Finite-rate condensation.** The mixing line assumes instantaneous
equilibrium across a 270 K temperature difference. Whether the condensation
keeps up with the mixing is an open question for LH₂ and would need a rate
law.

**One thermodynamic state per section.** This is the defining property of an
integral model. The profile integration above works around it for density; it
cannot work around it for a process whose *rate* varies across the section.
That is what a CFD code has and an integral model does not, and it is why
Giannissi's results resolve a shell that this cannot.

## Where that leaves the model

Removing the handover was the one that mattered most, because it was a free
parameter the answer depended on and it was of my own making. What remains is
the genuine boundary of the method: an integral model can carry a profile, but
not a reaction that lives in a shell within that profile.

For LH₂ specifically the consequence is bounded and stateable. The buoyancy
regime — grounded, lifting, airborne — comes out right, because it depends on
the integrated buoyancy and that survives averaging. Concentrations near the
source come out right, because the shell has not developed yet. What the shell
costs is accuracy in between, and that is where the residual sits.
