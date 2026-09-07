# Mixed N2/O2 liquid phase: bounded physical screen

2026-09-06, fixed before computing candidate flash or PRESLHY effect sizes.

## Identified assumption and scope

`axisymmetric_jet._condensed_air_state_exact.phase_at` sets each condensed
component's gas partial pressure to its PURE saturation pressure. It has no
shared N2/O2 liquid composition/activity. This is an independent-condensate
approximation, not the vapor-liquid equilibrium of a nitrogen/oxygen solution.

NBS Report3921 (Armstrong, Goldstein and Roberts,1955) measures this binary
solution at65,70,77.5K; see printed pp1,10-14,18-25 (PDF11,39-79). The report's
tables, equations1a/1b and definitions were visually checked in the original
scan. [Official original](https://nvlpubs.nist.gov/nistpubs/Legacy/RPT/nbsreport3921.pdf).
NIST's [Lemmon2000 air/mixture EOS description](https://www.nist.gov/publications/thermodynamic-properties-air-and-mixtures-nitrogen-argon-and-oxygen-60-2000-k-pressures)
also treats shared liquid and vapor phases, down to the composition-dependent
solidification boundary. It does NOT establish validity for H2-bearing solid
mixtures below that boundary. The existing model's phase assumptions are thus
worth screening; this is not yet an experimental dispersion improvement.

## First implementation

Add a separate, coefficient-free ideal-solution flash at63.151-120K:
`y_N P=x_N p_N_sat(T)`, `y_O P=x_O p_O_sat(T)`, `x_N+x_O=1`.
H2 is an insoluble gas in THIS LIMIT, not a measured solubility claim. Ideal
gas vapor, additive pure-liquid volumes and no excess mixing enthalpy are
explicit approximations. H2O ice remains a separate reservoir when present;
its saturation partial pressure and finite available mass must be respected.
Argon and dissolved H2 are not silently reclassified as N2/O2.

Solve the Rachford-Rice equation including the exactly non-condensing component.
No fitted flash quality, mixture coefficient, rate constant or trial correction.
Return both phase amounts and residuals; reject out-of-domain inputs. Do not
splice this liquid branch at63.151K into the old independent-solid closure.
Such a splice could create a nonphysical energy/density jump.

## Checks before actual-field diagnosis

- Exact component inventory, gas/liquid non-negativity, gas sum, Raoult equality,
  all-gas dew stability, pure-component and inert-only limits; extreme dilution.
- Rachford-Rice residual checked independently from the amount ledger.
- Nine pre-selected NBS measurements (runs11,14,23,24,33,37,40,48,52), spanning
  the three isotherms and liquid compositions, as an external property test.
  Use actual measured T,P,x,y; do not move them to nominal isotherms. Record
  ideal-solution error without fitting it away. Published activity correlation
  may be reproduced separately, labelled same-data correlation, not holdout.
- Original pressure-loss CONTROL seven fields replay before diagnostic reuse.
- Only local states with exact old T in63.151-120K enter the property screen;
  report excluded cold/warm states explicitly. Same T and bulk composition:
  compare condensed amounts, gas volume and enthalpy. Also solve temperature
  at fixed OLD enthalpy and composition if a root lies in the bounded domain.
  Fixed-state property changes are NOT new integrated field predictions.

Sample predetermined radial Gaussian shapes1,.9,.75,.5,.25,.1,.01 at the
handoff and all38 existing concentration stations, plus original41thermocouples.
Record shape/position coverage, maxima and source-trial breakdown; no new
dispersion score and no core/default promotion from this local diagnostic.

## Before any later full-model use

Need a thermodynamically consistent solid-liquid-vapor closure below63.151K,
gas non-ideality and H2 solubility bounds, phase-table value/derivative tests,
then conservative source/interface regeneration and actual downstream runs.
Adding a separate latent source on top of the phase enthalpy is prohibited.
