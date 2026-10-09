# Public open-channel H₂ source/sensor boundary audit

The USN/FFI open-channel dataset (DOI `10.23642/usn.26117989.v2`) contains a
logged source mass-flow channel and 29 hydrogen concentration channels per
experiment. It is useful for checking source/sensor timing and detector
threshold handling, but it is not an outdoor LH₂ plume experiment and does not
provide the weather, site-coordinate, obstacle or field common-clock package
needed for atmospheric LH₂ qualification.

## Local audit result

The read-only audit was run on 22 extracted CSV files:

- 22/22 files parsed successfully;
- 638 sensor channels were retained;
- 455 sensor channels crossed the declared 4 vol% threshold;
- threshold first-arrival times ranged from 36.89 s to 147.49 s across the
  archived runs;
- the source and sensor clocks are retained separately and are not implicitly
  treated as one synchronized clock;
- raw negative flow samples are reported, not clipped or silently converted to
  zero (the v2 audit counted 42,558 such samples).

The generated, hash-pinned local artifact is
`outputs/open-channel-h2-audit-2026-10-08-v2.json`. It is intentionally ignored
from the public snapshot because the raw dataset is third-party material.

Reproduce it when the licensed data are available:

```powershell
.venv\\Scripts\\python.exe tools\\audit_open_channel_h2_dataset.py `
  reference\\hydrogen\\open_channel_usn\\extracted `
  --threshold-percent 4 `
  --output outputs\\open-channel-h2-audit-2026-10-08.json
```

The audit schema is `degali.open-channel-h2-audit.v1`. Its status is
`boundary_only`, the serialized invariant is `promotion_allowed=false`, and the gate codes
include `not_lh2_pool_validation` and
`weather_geometry_common_clock_missing`. The integrated source mass uses the
recorded flow values as-is; negative samples are not clipped because clipping
would be an unregistered signal-processing choice.

## Use in the IJHE paper

This dataset can support a short supplementary boundary note on hydrogen
source-to-sensor timing and 4 vol% detector thresholds. It must not be pooled
with the PRESLHY or FFI/DNV concentration scores, used to fit a dispersion
coefficient, or described as validation of the atmospheric LH₂ model. The
appropriate claim is:

> The public open-channel H₂ data provide an independently timed source and
> detector-timing boundary, while their missing atmospheric geometry and
> meteorological evidence prevent promotion to an LH₂ dispersion-validation
> case.
