# Sub-metre atmospheric-spreading consistency screen

2026-09-06. Registered after the matched-flow audit, before candidate field runs.

## Located inconsistency

The research crosswind path inherits this legacy JETPLU convention:

    sigma_ya = delta_y * x**beta_y
    sigma_za = delta_z * x**beta_z * exp(gamma_z * log(x)**2)

Both widths are evaluated for every x>0, including the present LH2 handoffs
at0.2--0.8m. However their along-path derivatives used in ambient entrainment
are set to zero at x<=1m. The actual width functions are not constant there.
This convention is present in legacy/Deg_for/JETPLU.FOR lines354--363; it is
not a newly introduced Python translation error.

At the frozen matched-flow diagnostic stations below1m, differentiating the
widths already in use produces a missing ambient mass-source contribution
of4.8--13.5% of the implemented total local mass source. These percentages
are local equation terms, not percentages of observed error or cumulative
entrained mass. No new field candidate had been run when this was found.

## Prospective candidate

For x>0 use the analytic along-path derivatives everywhere:

    d sigma_ya / ds = sigma_ya/x * beta_y * cos(theta)
    d sigma_za / ds = sigma_za/x * (beta_z + 2 gamma_z log(x)) * cos(theta)

Keep the existing atmospheric widths, all coefficients and plume equations.
The existing source formula then transports additional ambient material as

    dM = rho_a * U_a * I0 * d(sigma_ya*sigma_za)/ds
    dPx = U_a*dM; dPz = 0; dF_H2 = 0; dE = U_a**2/2*dM

where E is ambient-relative total enthalpy plus mean kinetic energy. Ambient
enthalpy is already zero in the phase-consistent convention. Do not add a
latent-heat or ambient-enthalpy source twice. Beyond1m the implementation
returns the exact historical geometry/source calculations. At x<=0 the new
option refuses use, since these empirical powers are not regularized there.

This is a consistency option within an inherited atmospheric spread model,
not proof that its empirical atmospheric correlation is experimentally valid
below1m. The physical applicability remains a limitation. The legacy Fortran
port, Stage1 defaults, source plane, phase model and coefficients are kept.

## Tests and comparison

Before field scoring, independently differentiate the width-product function
at0.2,0.5,1,2m; verify the derivative and the five source increments. Check
zero-spreading, zero-wind, x>1 exact parity, continuity at1m, original-object
immutability and invalid coordinates. All initial fluxes, temperature and
width must reproduce the passed original interface exactly.

The primary comparison uses the Stage1 pressure-loss/TC3/PT2 density CONTROL,
with identical initial states and all original settings. First run trials10
and23 to the same arc-length limits and fixed41 temperature sensors. Then
extend to the remaining five cases if the path is numerically valid, retaining
the exact38/17/41 population. Check conservation<=1e-5 and .02/.01m sensitivity
before any improvement/adoption claim. Report all three temperature statistics,
concentration, width and centre-height changes, including adverse changes.

The added ambient mixing may dilute hydrogen and warm the cloud. Its effect
on buoyancy, trajectory and off-axis sensors is not assumed monotone. No
coefficient or gate is adjusted based on these predictions. A failure remains
recorded; a physical correction is not equated with overall accuracy improvement.
