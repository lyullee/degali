# FFI/DNV site-geometry boundary for the IJHE package

The FFI/DNV outdoor release report describes a test pad with two containers
stacked above one another, a plastic drum and an instrument box near the
release point. It also describes a horizontal-release orientation, a mast with
high and low wind sensors, and oxygen sensors whose depletion signal was
translated to hydrogen concentration. This is relevant context for the six-arc
FFI screen: the published campaign is not a featureless, obstacle-free field.

## Evidence provenance

- Public source: Aaneby, Gjesdal and Voie, *Large scale leakage of liquid
  hydrogen (LH2)*, FFI Report 20/03101, cited in the manuscript.
- Official public report page: https://www.ffi.no/en/publications-archive/large-scale-leakage-of-liquid-hydrogen-lh2-tests-related-to-bunkering-and-maritime-use-of-liquid-hydrogen
- Local report inspected:
  `<local-slabx-repository>/slabx-lh2/액화수소 검증/21-03101-잠금 해제됨.pdf`
- Local extracted text inspected:
  `<local-slabx-repository>/tmp/report_weather_instrument_pages.txt`
- SHA-256 of the local PDF:
  `D500B61DA23A44043648D71E06EE7D061E7B46671AF5ECC8FC370CF9C9A23A4A`
- SHA-256 of the local page extraction:
  `A8A6474104EA277ED6C3125AAC63C13DBB59CF615B8659E7CB2A9E4930D7CF94`

The extracted pages are an audit aid and are not redistributed in the public
snapshot. The report itself is the citable source.

## Modeling consequence

The current DEGALI FFI calculation does not resolve those containers, drum,
instrument box, mast, or any three-dimensional wake field. The FFI lane must
therefore remain a **conditional far-field/transport screen**, not obstacle
validation. The known pad geometry is recorded to prevent a reviewer from
interpreting the six-arc result as a clean free-field obstacle test. The
geometry may contribute to local residuals, but the public evidence does not
identify its effect separately from source-state, wind/time variability and
the unresolved transient/two-phase transport.

This note does not promote the FFI lane, change the headline metrics, or supply
a machine-readable measured geometry. A measured three-dimensional geometry
and common-clock receptor package remain required for an obstacle-resolved
operational claim.
