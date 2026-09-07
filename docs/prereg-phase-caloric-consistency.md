# Phase caloric consistency and bounded gamma-O2 potential

2026-09-06, specified before the numerical candidate audit.

## Diagnosis

The frozen air closure uses piecewise linear log(P) vs 1/T solid pressures,
but one constant sublimation enthalpy per species. Its liquid latent heat is
real saturated vapor minus liquid enthalpy, subtracted from an IDEAL gas
enthalpy reference. Neither construction guarantees a common Gibbs potential.
The temperature grid also interpolates across triple points without explicit
phase fractions. Quantify these separately; do not hide a physical fusion
enthalpy with smoothing or add a second latent heat source.

## Sources inspected

- NBS Circular564 (1955), Table7-11 and7-11/b, printed360-361/PDF374-375:
  [official PDF](https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular564.pdf).
  The published solid nitrogen expression is log10(P/mmHg)=7.65894-359.093/T.
  The table includes62K,73.6mmHg, absent from the frozen transcription.
- Hans M. Roder, NBSIR77-859 (1977), sections6.1-6.5, printed13-14/PDF19-20:
  [official PDF](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir77-859.pdf).
  Original equations and Table3 were visually checked. The report provides
  gamma-solid O2 heat capacity, molar volume, fusion heat106.3cal/mol, and a
  thermodynamic pressure correlation. Its caloric/volume functions reproduce
  older measurements; this is not an independent new data set. It explicitly
  notes up to30% disagreement with older vapor pressure measurements.

## Bounded implementation, separate from the frozen core

Build a gamma-solid O2 Gibbs potential on43.801--54.361K,0<P<=200kPa only.
The upper endpoint is the current CoolProp triple temperature,0.002K above
the1977 value; record this reference adaptation. Do not extend gamma-solid
below its beta transition or use this as a N2/O2 solid-solution model.

At the current pure liquid triple point, set hs=hl-106.3*4.184J/mol and
ss=sl-(106.3*4.184)/Tt, with the SAME CoolProp reference. Integrate the
published heat capacity polynomial for h and s, initially approximating its
saturation-path heat capacity as isobaric Cp at Pt. Quantify the omitted
T*v'(T)*dP_sat/dT term against the reported0.5% caloric uncertainty; this
approximation must not silently be called exact experimental thermodynamics.

Use gs(T,P)=href(T)-T*sref(T)+(P-Pt)*v(T). Derive every returned property
from that one potential: s=sref-(P-Pt)*v', h=href+(P-Pt)*(v-T*v'),
Cp=Cref-T*(P-Pt)*v''. The solid is incompressible at fixedT in this limit.
The published volume polynomial is used as-is, not fitted to plume data.

Use the matching CoolProp ideal-gas h,s and its gas constant for the gas
chemical potential, with the pressure entropy transformed analytically.
Return the ideal-mixture equilibrium partial pressure from equality to gs
at the TOTAL mechanical pressure. For a pure vapor, solve totalP=partialP.
At the triple point this ideal-gas equilibrium pressure is the REAL vapor
fugacity, not exactly its saturation pressure. Explicitly distinguish them.
Do not copy a real-vapor latent heat into the ideal-gas ledger.

## Predetermined verification and reporting

- Finite inputs/domain rejection; h=g+Ts; -dg/dT=s; dg/dP=v; dh/dT=Cp;
  Maxwell ds/dP=-dv/dT. Positive Cp and volume throughout the domain.
- Exact triple-point Gibbs equality and specified fusion h and s differences.
- Pure vapor equilibrium and Clapeyron agreement using returned volumes and h.
- Compare against ALL12 Table3 enthalpy entries and the published pressure
  equation, labelled shared-source reproduction with a changed gas EOS.
- Audit ALL frozen solid pressure intervals, triple pressure mismatch, latent
  changes and pure liquid ideal-reference residual at fixed63.151,65,70,77.5,
  90,100,120K where valid. Quantify the artificial grid transition width.
- Hash inputs/code, run related tests, verify Stage1's134files unchanged.

No new integrated field score or default promotion. Next requires a bounded
nitrogen potential, mixed-solid stability/equilibrium, and conservative
source/interface regeneration; the pure gamma-O2 operator alone is not a
completed mixed air or LH2 dispersion model.
