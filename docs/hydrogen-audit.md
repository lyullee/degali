# What liquid hydrogen breaks, and what it does not

DEGADIS was built for LNG and ammonia. Hydrogen is four hundred kelvin and a
factor of eight in molecular weight away from either. This is an audit of
every module against that, with the fixes and — as important — the things that
turned out to be fine.

Thermodynamics is CoolProp throughout when `backend="coolprop"`: water, air
and the contaminant all come from equations of state.

## Fixed

**The adiabatic mixing table outgrew its array bound.** `IGEN = 42` is a
Fortran dimension. The adaptive thinning keeps whatever nodes linear
interpolation needs, and hydrogen's mixing line spans 20 to 289 K — a factor
of fourteen against LNG's under three — so it needs 46. The bound is raised
rather than the tolerance loosened, because loosening it would degrade every
lookup silently.

**The equivalent-source search had a fixed temperature floor.** 120 K, chosen
so an ammonia answer near 205 K could not run away. Hydrogen's lands near
20 K, so the floor excluded the answer. No fixed value serves both; the
fluid's own triple point governs now.

**Hydrogen's critical temperature is below ambient.** 33.1 K against 289 K.
The closure condition — the vapour is exactly saturated at its own partial
pressure — has no meaning above the critical point, so the search is capped
there. A solution at the cap is the physical answer, not a failure: it says
the aerosol finished evaporating before the mixture warmed that far. Ammonia,
with a critical temperature of 405 K, never meets this.

**The CoolProp property grid started at 80 K.** Fine for methane at 112 K and
ammonia at 240 K. For hydrogen at 20 K the widening rule produced a lower
bound of **−30 K**: a sizeable part of a 4001-point grid spent on temperatures
that do not exist, and coarser interpolation over the ones that do. The grid
now starts from the contaminant's own boiling point and is floored at 4 K.

**`SETJET` assumes a vertical jet.** Its Kamotani and Greber trajectory is for
a jet issuing vertically into a crossflow, and 19 of the 24 HSL releases are
horizontal. The trajectory correlation is dropped for a directed release; the
zone-of-flow-development length is kept, because that is a property of the jet
rather than of its direction.

**The ground-level closure had no way to let a buoyant cloud rise.** Made
swappable; see [liftoff.md](liftoff.md).

## Checked and adequate

**`SERIES` overflow guard.** The flammable-mass integrals overflow above
`ln(cc/c) = 13.8`. Hydrogen's source-to-LFL ratio gives 6.6, and even
source-to-one-per-cent-of-LFL gives 11.2. There is margin, but less than for
LNG's 4.0, and a query far below the LFL would trip it.

**`PHIF` with a negative Richardson number.** The entrainment suppression is
written with `abs(Ri)` on the unstable branch, so a buoyant cloud is handled
rather than producing a domain error.

**Water condensation.** At 20–60 K the water vapour pressure is effectively
zero and the enthalpy already adds the heat of fusion below 273 K. Nothing
needed.

## Known limits, not fixed

**Air condenses and is not modelled.** The mixing line passes 90 K, oxygen's
boiling point, at 81 mole per cent hydrogen, and 77 K, nitrogen's, at 84 per
cent. Real LH₂ clouds do liquefy air, with oxygen condensing preferentially —
a recognised hazard in its own right. `Thermo` models water condensation only.
Below about 100 K the air heat capacity is clamped, so the mixture properties
there are approximate.

This does not change the buoyancy conclusion: it happens well inside the
buoyant range, and condensing air out of the gas phase makes what remains
lighter still.

**Ground heat transfer is extrapolated.** `SURFAC`'s natural-convection term
goes as `((rho/M)² ΔT)^(1/3)`, and hydrogen has both a small molecular weight
and a large temperature difference. For a 40 K cloud under a 289 K surface it
gives 76 kW/m² against 6 kW/m² for a comparable LNG cloud — a factor of
twelve, of which the molecular weight contributes twenty and the temperature
difference under two.

That magnitude is not obviously wrong: measured fluxes for LH₂ on concrete are
of that order. But the correlation was fitted for heavier gases and a smaller
temperature difference, and it is being used well outside where it was
established.

**Rainout is not modelled.** The equivalent-source calculation assumes all the
flashed liquid evaporates. The HSL trials were designed to study rainout
precisely because it happens, so this assumption is the weakest link in the
source term and should be checked against the report before results are
quoted.

## The correction the audit forced

An earlier reading of the mixing line said a hydrogen cloud is buoyant
everywhere below 99.9 mole per cent. That used the *pure vapour* line, which
is not what a release produces.

A flashing jet is much colder at a given concentration, because the liquid's
latent heat has to come from the air it entrains. Starting from the computed
equivalent source instead:

| source | source rho/rho_a | dense above |
|---|---|---|
| pure vapour | 1.094 at 20 K | 99.9 mol % |
| flashing jet, 1 barg | 1.134 at 41 K | 85.2 mol % |
| flashing jet, 5 barg | 1.029 at 41 K | 90.2 mol % |

The dense range is a hundred times wider than the vapour line suggested, and
DEGADIS's dense phase has real work to do near the source.

The conclusion survives, more precisely stated: the lower flammable limit at
4 mole per cent, stoichiometric at 29.5, and the *upper* flammable limit at 75
are all buoyant. The dense phase exists only above the UFL.
