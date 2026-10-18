# TI021 3-D handoff packet audit — 2026-10-15

## Session objective

Repeat **cause search → physical repair → identical 600 s dose calculation →
conservation/path gate** until the TI021 positive-sensor dose MAPE is at or
below 10%. A candidate is not promoted from a packet audit alone.

## Evidence checked

The durable inner handoff is
`outputs/thorney-ti021-wide-inner-domain-600s-x450y640-centered-2026-10-13-v1/inner/ti021-inner-handoff.h5`.
It contains 4,801 frames over 600 s, with an 8 × 128 height/transverse plane
(5 m cells). Using the existing lazy reader and exact rectangular overlap
remap, with the declared ambient flux `rho_air * 3.4 m/s = 4.0885981633
kg m-2 s-1`, gave:

| check | result |
|---|---:|
| integrated positive species packet | 3,292.59968495 kg |
| integrated total-mass anomaly | 2,501.94815383 kg |
| peak species packet rate | 96.11726259 kg/s at 33.375 s |
| minimum remapped total flux | 1.1812033253 kg m-2 s-1 (600 s) |
| minimum remapped species flux | 0 kg m-2 s-1 (non-negative) |
| maximum overlap-remap residual | 1.94e-16 |

The packet is therefore admissible for a conservative interface **after** the
ambient carrier is restored. The raw total anomaly must not be passed to a
receiver as a standalone mass flux: it is signed and represents only the
departure from ambient flow.

## Promotion decision

**Not promoted to a MAPE candidate yet.** The integrated species packet
(3,292.6 kg) is essentially the same order as the already handed-off source
inventory. Injecting it into the current dense-layer run would double-count
the release unless the dense-layer initial state is replaced by a true
interior interface state. The next implementation must therefore be an
interior, time-resolved, height-to-layer packet exchange (or a clean
replacement of the current initial handoff), with a mass ledger proving that
the packet is not counted twice.

Current promoted best remains
`outputs/ti021-semidisk-backflow-reservoir-c005-r005-radius10-strict.json`:
MAPE **29.5441456118%**, species residual `1.00e-11 kg`, momentum residual
`9.31e-10 kg m/s`, runtime about 16 s. The 10% target gap is **19.5441
percentage points**.
