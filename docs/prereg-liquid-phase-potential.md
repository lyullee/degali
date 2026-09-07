# Common liquid/ideal-gas potential —2026-09-06

Recorded before candidate benchmarks and local effect sizes. Do not modify
frozen Stage1 or the completed N2/O2 solid modules and liquid-screen outputs.

## Bounded candidate

Use stable pure saturated-liquid HEOS states as references, never force an
unknown PT state to be liquid. Official CoolProp documentation warns that
imposing the wrong phase can return incorrect results:
[phase warning](https://coolprop.org/coolprop/HighLevelAPI.html#imposing-the-phase-optional).
Use its analytical saturation derivatives:
[API](https://coolprop.org/coolprop/LowLevelAPI.html#two-phase-and-saturation-derivatives).
The installed second saturation derivative accepts P,P, so transform those
derivatives to T analytically, not by an unreported finite-difference step.

For each species, construct G_l(T,P)=G_sat(T)+(P-P_sat(T))*v_sat(T).
This is the first-order pressure expansion, incompressible at fixed T, NOT an
exact high-pressure liquid EOS. All properties follow from this one G:
s=s_sat-dP*v_sat', h=h_sat+dP*(v_sat-T*v_sat'),
Cp=h_sat'-P_sat'*(v_sat-T*v_sat')-T*dP*v_sat''.
Use O2>=54.361K and N2>=63.151K, T<=100K,0<P<=200kPa. The100K cap is a
declared diagnostic range to avoid near-critical extrapolation, not a measured
validity boundary. Reference liquids below their pure triple points are absent.
Compare the pressure approximation to HEOS PT only in known stable-liquid
regions and record its error. Mathematical consistency is not physical accuracy.

Use HEOS ideal h/s/Cp at a common100kPa standard reference and universal
R=8.31446261815324J/mol/K in mixing and gas-pressure terms. Backend fitted
gas constants differ slightly; record them, do not mix different R values in
one ideal mixture. This convention leaves standard h/s unchanged and gives
one ideal-mixture chemical-potential and volume rule. No reference-state resets.

Add the ideal-liquid mixing term RT*sum(x_i ln x_i) to G. Liquid h,Cp,v are
linear sums under this zero-excess-enthalpy approximation. Derive K from the
liquid/ideal-gas chemical-potential difference at mechanical total pressure,
not merely from real-vapor saturation pressure. H2 is insoluble ideal gas;
this is a limiting hypothesis, not measured zero solubility. Reuse the tested
Rachford-Rice solver for gas/liquid amounts. Do not add the published NBS
activity factors to equilibrium without consistent excess h/s/G derivatives.
Mixed-solid stability and subtriple liquid remain explicitly unimplemented.

## Independent checks and effect audit

- Pure saturated G/h/s/v reproduction; analytic versus numerical T/P
  derivatives, Maxwell identities, pressure approximation in stable liquid.
- Match completed pure-solid triple G and measured fusion h/s.
- Ideal mixture extensive G derivatives versus component chemical potentials,
  Euler and Gibbs-Duhem identities; equilibrium component mu equality;
  independent constrained Gibbs minimization and inventory/scale/zero limits.
- Compare new K-based ideal-liquid bubble P/y to the already-transcribed
 9 NBS3921 rows without refitting. This is a property check, not plume validation.
- Reuse the hash-verified completed local-screen samples; quantify only
  in-domain changes, including all exclusions. Isolate the caloric-reference
  correction from equilibrium and volume changes. No replacement sensor scores.
- Preserve frozen134 files, source/code/input hashes and failed checks.

Next actionable gate: a bounded conservative phase closure in actual source
conditions. Do not splice mixed liquid to independent pure solids at63.151K;
obtain mixed-solid stability or leave that domain explicitly unsupported.
