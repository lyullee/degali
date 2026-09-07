# Beta-N2 caloric potential and reference check

2026-09-06, before computing candidate benchmarks/effect sizes.

## Data and bounded assumptions

Use Ziegler & Mullins(1963), Technical Report1, Project A-663, Georgia Tech,
prepared under NBS contract CST-7404. The original report and July2 corrections
were downloaded from the author's institution through its public repository
API. The UI download URL returned HTML, not a denied scientific file; the API
returned the PDF with MD5 matching the repository's171058caf8abd521e1e871c8e3bc0de7.

[Original UI link](https://repository.gatech.edu/bitstreams/4d6e76d0-4a55-4486-ab62-499c4013e7c4/download)
and [direct public PDF](https://repository.gatech.edu/server/api/core/bitstreams/4d6e76d0-4a55-4486-ab62-499c4013e7c4/content).

Printed7,10-13,55-56 and correction2 (PDF33,36-39,81-82,16) visually checked.
Use published beta-solid Cp polynomial TableIV, fusion172.3cal/mol TableIII,
and the explicitly approximate mean molar volume28.67cm3/mol. No new data fit.
Cp is based on Giauque & Clayton1933, whose publisher PDF is access-blocked;
this is reproduction of the1963 authors' caloric representation, not a new
independent measurement. Record solid thermal expansion and compressibility
as omitted by the mean-volume approximation. Do not infer precision from it.

Domain35.62--63.151K,0<P<=200kPa. Preserve phase identity: beta-solid only,
not alpha-N2 below35.62K or stable liquid above the triple point. The report
uses63.152K; anchor at current CoolProp63.151K and record the0.001K change.

## Candidate

As for the completed gamma-O2 potential, anchor h_s=h_l-L_f and
s_s=s_l-L_f/T at the current liquid triple point. Integrate the published
Cp(T) and Cp(T)/T. Construct G(T,P)=h_ref-T*s_ref+(P-Pt)*Vbar. Return h,s,Cp,
volume and gas equilibrium from the same potential. The constant-volume
approximation makes Cp and saturation heat capacity equal within this model.
Use matching HEOS ideal-gas enthalpy/entropy and its gas constant. Never mix
the real-vapor saturation latent heat with the ideal-gas enthalpy ledger.

## Predetermined checks

- Domain rejection, positive Cp/volume, triple Gibbs and fusion h/s equality,
  finite-difference Gibbs/Maxwell/Clapeyron identities, independent quadrature.
- Corrected TableXVI beta-solid integer rows36--63K (28rows), plus the report's
  lower transition endpoint35.62K. Record63.152K as outside the new domain,
  not silently shift a comparison datum. Use July2 corrected41--56K values.
- Pressure and latent comparison are shared-source reproduction with changed
  reference/gas EOS, not held-out data. Compare both the ideal-gas candidate
  and a real-vapor evaluation on the SAME solid potential when stable HEOS
  gas evaluation is available. Do not change the plume gas to this pure-gas EOS.
- Quantify old versus candidate enthalpy, Cp and density at59,60,61,62,63K,
  the range directly present in stored PRESLHY handoff/local states. These
  are pure material changes, not predicted sensor shifts.
- Run related tests and frozen134file validation; hash original documents,
  corrected data, implementation and audit tool. No default promotion.

Next: build the consistent pure-liquid reference and bounded mixed-phase
stability needed to couple these components. The state of mixed solids cannot
be assumed from pure-species triple points. Do not run invalid-domain fields
to create a premature new dispersion score.
