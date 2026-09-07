# Pre-registration: argon phase completeness in the Raman jet

**Frozen before the first argon-enabled model result.** Do not edit above the
`RESULTS` line after a candidate result is known.

## Question

The accepted dry-air equilibrium model represents nitrogen and oxygen but
renormalises them after omitting atmospheric argon. Argon is only about
0.934 mol % of dry air, yet its 83.806 K triple point lies inside the measured
cryogenic jet range. Omitting it therefore removes both condensed volume and
phase enthalpy from the coldest part of the plume.

## Candidate fixed before running

1. Represent dry air as N2/O2/Ar mole fractions 0.78084/0.20946/0.00934,
   normalised over those three explicitly modelled components. Do not fit the
   composition to the Raman data.
2. Above the argon triple point, use CoolProp liquid-vapour saturation,
   density and latent heat. Below it, use a Clausius-Clapeyron solid-vapour
   extrapolation anchored to the CoolProp triple pressure and a fixed
   7.79 kJ/mol sublimation enthalpy. Hold condensed density at CoolProp's
   saturated-liquid triple-point value below the triple point; condensed
   volume is negligible at this atmospheric loading, and this avoids inventing
   a solid-density fit. These are explicit property-data approximations, not
   tunable closures.
3. Include component ideal-gas `h(T)` for argon on the same 14.1--400 K table
   as H2, N2, O2 and H2O. Count phase enthalpy exactly once.
4. Keep every accepted model choice and numerical setting unchanged:
   `scalar_peak` establishment, density/species profiles, zero humidity,
   81 radial points, 0.25 mm maximum step and `5e-8` relative tolerance.

NIST SRD 69 reports argon's molecular weight and phase data, including an
83.78--83.8 K triple temperature. CoolProp supplies the internally consistent
fluid values used by the calculation. The low-temperature sublimation value
is independently reported as 7.79 kJ/mol in the NIST thin-film measurement.

## Decision rule

Adopt argon in the recommended dry-air research configuration only if:

1. all boundary residuals are below `1e-8`, and species/energy drift remains
   below `2e-4` in every case;
2. all four Raman slopes and both radial coefficients still meet the existing
   acceptance limits under the corrected 369-point protocol;
3. every state is physical and a 121-point rerun changes every slope by less
   than 0.2%; and
4. no empirical coefficient is introduced.

The change is called materially relevant only if at least one printed slope
moves by 0.2% relative to the N2/O2 result. A smaller change may still be kept
as a composition-completeness correction if all four rules above pass, but it
must not be described as a validated error reduction.

---

## RESULTS

Interim result retained across an intentional laptop-move interruption: the
81-point original 549-point audit gives slopes `0.21565704`, `0.05969394`,
`0.02211815`, and `0.07536714`; all remain inside 25%. Relative to the
accepted N2/O2 model, the changes are +0.28%, +0.09%, -1.49%, and -0.56%.
Argon therefore exceeds the 0.2% materiality threshold, but does not move all
errors in a favourable direction: centreline-temperature error worsens while
temperature-width error improves.

Primary 369-point coverage, conservation aggregation and 121-point
convergence subsequently completed. At 81 points the primary slopes are
`0.20881427`, `0.05957851`, `0.02190436`, and `0.07667047`, with relative
errors -20.482%, -8.383%, -22.490%, and +23.902%; both radial coefficients
(`60.1982`, `34.9909`) remain inside the reported ranges. Maximum boundary,
species and energy errors are `1.59e-14`, `1.50e-5`, and `2.95e-5`.
Temperatures remain 50.95--230.74 K and all states are physical.

The 121-point slopes are `0.20887131`, `0.05957544`, `0.02191007`, and
`0.07666313`. Relative changes from 81 points are +0.0273%, -0.0052%,
+0.0261%, and -0.0096%, so the 0.2% convergence rule passes.

**Decision:** retain `equilibrium_argon_condensation=True` as an explicit
composition-completeness sensitivity, but do not promote it into the
recommended error-reduction configuration. It improves mass centre, mass
width and temperature width slightly, but worsens the centreline-temperature
error by about 1.41% relative to the accepted N2/O2 result. The candidate is
conservative and parameter-free but does not demonstrate a consistent error
reduction.

### Final-journal source correction

The revised primary errors are -24.64%, -15.72%, -22.49% and +23.90%.
The candidate remains provisionally 4/4 but still fails the pre-registered
consistent-improvement rule. See `hecht-panda-journal-benchmark-correction.md`.
