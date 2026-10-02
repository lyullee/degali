# Validation integration register

Date: 2026-09-17

## Why a register instead of one headline score

DEGALI now has validation evidence for the original dense-gas reconstruction,
cryogenic gas jets, outdoor LH2 releases, a cryogenic source transient, and
several mixed-phase and turbulent-energy mechanisms. They do not share a
source boundary, sensor geometry, averaging window, observable or uncertainty
model. Pooling their MG/VG/FAC2 values would therefore create a number with no
physical interpretation.

The executable register in `degali.validation.integration` records the claim
each evidence item can support and rejects automatic default promotion. It
ships no external observation: all third-party files remain ignored under
`reference/`.

Run it on any clean clone:

```powershell
.\.venv\Scripts\python.exe tools\audit_validation_portfolio.py
```

## Integrated status on 2026-09-17

| branch | evidence that can support it | integrated status | what it cannot establish |
|---|---|---|---|
| Legacy DEGADIS reconstruction | source-built Fortran oracle | reconstruction only | LH2 accuracy |
| Cryogenic gas axisymmetric jet | Sandia Raman mean profiles | provisional research screen | humidity-conditioned predictive accuracy |
| Fast single-velocity LH2 path | PRESLHY E3.5 and independent Spadeadam reported arcs | qualified field screen in declared conditions | universal accuracy, transient meander or phase-resolved behaviour |
| Cryo-compressed source | PRESLHY E3.1 80 K pressure histories | qualified source screen | a universal wall/pipe/two-phase source closure |
| Finite TKE transport | LES lower bound plus public PIV capability record | research only | experimental TKE, epsilon or length-scale closure |
| Finite-rate droplets / slip | coefficient-free source bounds | research only | particle size, phase inventory or gas-particle slip prediction |
| Equilibrium air-condensation limit | Li (2026) limiting comparison | rejected research comparison | a sub-triple N2/O2 phase closure |
| Time-resolved pool-source sequence | E3.4 source intervals and public FFI channel timing intake | quasi-steady source-to-plume boundary | outdoor arrival, puff storage, meander or turbulent dispersion |
| Cuboid/wall site geometry | exact contact screen plus public neutral Case-H intake | validity/invalidation screen | a quantitative LH2 obstacle-wake correction |

“Qualified field screen” is intentionally not a general accuracy claim and is
never a default-promotion decision. It requires an independent experimental,
uncalibrated prediction with the decisive state observed. The register makes
no such claim for a numerical LES, a manufactured conservation test, a
particle-tracer PIV display, or a model-form limiting comparison.

## Promotion rule

A new branch must separately provide an experimental result that is all of:

1. an uncalibrated prediction on the same declared branch;
2. independent of branch selection;
3. direct on the state that makes the new physics necessary; and
4. reported with its original sensor/source limitations.

Until then it remains research-only even when it satisfies conservation, PSD
or numerical-convergence gates. This is why finite TKE needs gas-velocity RMS
and epsilon/length data, while phase/slip transport needs particle-resolved
measurements. No model coefficient is selected from the thermal or
concentration residuals of the same campaign.

## Relationship to existing audits

The register does not replace the detailed calculations. It points to their
separate scopes: `validation.md` for reconstruction, `field-validation.md`
for non-LH2 field trials, `lh2-data-readiness.md` for data capability,
`public-turbulence-evidence-audit.md` for TKE provenance, and
`cold-mixed-phase-evidence.md` for the N2/O2 closure boundary.
