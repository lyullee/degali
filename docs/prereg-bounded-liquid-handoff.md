# Bounded liquid-only downstream handoff —2026-09-06

Use the original pressure-loss density CONTROL10/23 handoff states and5fluxes,
not Coriolis and not a newly fitted supply. This is a partial downstream-property
contrast; the upstream LH2/air formation model stays unchanged. No full-source
common-EOS claim. Preserve field/default/source files and previous failed gates.

Construct a separate wrapper around the old explicit N2/O2/H2O atmosphere.
Replace only N2/O2 air-condensation partition/volume/latent reference with the
completed ideal-liquid G. Retain legacy ideal-gas h interpolation and water
phase ledger; quantify numerical and residual reference approximations.
H2 is an insoluble ideal gas, without a claim of measured zero solubility.
Use a64K lower implementation guard, above the pure N2 triple vicinity,
NOT a newly fitted phase boundary. Do not continue below it, blend with pure
solids or extrapolate the liquid table. This excludes mixed-solid work from
the partial candidate; omitted mixture stability remains an explicit limitation.

The liquid table ends100K. Above100K retain only the gas/water branch, after
checking both pure ideal-liquid equilibrium K at100K exceed1 at the operating
pressure and the old saturation pressures remain>P in the warmer table.
No liquid properties above100K are extrapolated. Preserve water inventory.
For fixedP use Hermite interpolation of each liquid and gas standard G with
slope -s. Generate h by G-T*G', not an independent latent interpolation.
Interpolate liquid volume consistently using v,v'. Test versus exact potential
at non-grid states before any handoff result; retain current species mole basis.

Use a vectorized, inventory-checked binary/inert flash and monotone density
inversion only inside64--300K. The quadratic form of the3component
Rachford-Rice equation is a numerical acceleration, not a new physical fit.
Test it independently against the existing bracketing solver and all endpoint,
trace, scale and reservoir limits. Check gas-only parity, monotonicity,
forward/inverse residuals, water inventory and explicit domain failures.

Project the same target5fluxes at unchanged x,z and velocity-width exponent.
Keep existing1e-8 projection and1e-5 quadrature checks, geometry/temperature
interface criteria as originally specified. Record a candidate that conserves
fluxes but changes width/temperature beyond those criteria as failed, not
accepted. Check centre AND quadrature local temperatures. A valid fit alone
does not justify integration with unsupported local states. If an optimizer
trial hits the guard, retain counts/diagnostics without mistaking a rejected
optimizer probe for a physical trajectory crossing.

Only accepted boundaries may launch actual10/23 fields. Record domains and
conservation failures with last saved states, not fabricated coverage/scores.
If the partial candidate is blocked, identify exact missing mixed-phase states
and continue independent downstream residual analysis instead of repeating
the same property or synthetic-turbulence tests.
