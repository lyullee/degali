# Continuous-barrier physics session goal — 2026-10-18

## Objective

Replace the animation-only scalar barrier path with DEGALI's coupled
variable-density momentum, pressure-projection and thermal/species solver.
Only render a wall case after the calculated field demonstrates the obstacle
mechanisms that were absent from the old path.

## Acceptance gates

1. Zero normal velocity and zero transported hydrogen through solid faces.
2. Conservative hydrogen ledger and converged variable-density continuity.
3. Windward deceleration and higher windward pressure than the lee region.
4. Acceleration and shear over the wall top.
5. A resolved lee reverse-flow region and negative vertical velocity
   (downwash), rather than a camera-induced appearance of descent.
6. The 120 s coarse-grid demonstration must solve in less than 30 s on the
   local workstation.
7. Concentration and arrival-time prediction remain conditional until the
   same path passes an external obstacle experiment/CFD comparison; internal
   mechanism checks alone do not authorize digital-twin safety use.

## Current evidence

The 120 s blocking-wall run in
`outputs/continuous-blocking-wall-3d-resolved-physics-2026-10-18` passes gates
1–6.  The authoritative numerical values and implementation hashes are in
`visualization_metadata.json`.  Gate 7 is intentionally still open.

