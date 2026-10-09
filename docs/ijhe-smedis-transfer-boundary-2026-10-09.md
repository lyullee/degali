# Non-LH2 SMEDIS/REDIPHEM transfer boundary — 2026-10-09

This note records a useful external-data finding without promoting it to the
LH2 validation claim.

## Candidate package found

The separate SLABx workspace contains 28 SMEDIS workbook files under
`파이썬 코드/validation-data/smedis/`, together with reduced CSV tables for
conditions, receptor observations and arc maxima. The package covers generic
dense-gas experiments, including ammonia, propane, R-12, SO2 and LNG/SF6/CO2
source records. The reduced tables expose source rate or equivalent-source
information, meteorology, receptor coordinates and concentration or dose
observations for several trials.

The same workspace also records an EEC control/fence pair (EEC360/EEC361)
and the obstacle-aware SMEDIS inventory. Those files are useful for testing a
no-fit observation operator and for rejecting a scalar wall-attenuation
factor. They do not provide a validated LH2 wake closure.

## Why it is not promoted into the IJHE LH2 score

The candidate package fails the LH2 promotion boundary for four independent
reasons:

1. none of the identified source records is an LH2 release;
2. the source thermodynamics and molecular-weight regime are different from
   the PRESLHY/FFI LH2 lanes;
3. the reduced observations are static means, doses or arc maxima rather than
   a common-clock LH2 concentration time series with channel calibration; and
4. obstacle cases retain a fence count or reduced geometry, not a complete
   three-dimensional site geometry and matched wake concentration field.

The data therefore remain a **non-LH2 transfer boundary**. They may support a
separate generic dense-gas benchmark in a future paper, but their records are
not pooled with E3.5, FFI/DNV, the Hecht–Panda boundary or the conditional
SLABx comparison. The original workbooks are third-party material and are not
redistributed by the DEGALI IJHE package.

## Reproducibility and decision

The repository already contains the generic reader and comparison boundary in
`src/degali/validation/smedis.py`, `src/degali/validation/compare.py` and
`docs/field-validation.md`. The local audit of this finding is documentary:
it records dataset identity and promotion status, not a new LH2 score. A
future promotion would require a declared generic-gas protocol, licence and
hash manifest, independently reproduced results, and a separate claim lane.

The repository field-evidence audit scanned the reduced candidate directory
and returned `status=withheld` with `promotion_allowed=false`. Its five
promotion channels (`source_boundary`, `weather`, `obstacle_geometry`,
`receptor_observations` and `common_clock`) were all `missing`; the machine
readable record is
`outputs/smedis-transfer-field-audit-2026-10-09.json`.
