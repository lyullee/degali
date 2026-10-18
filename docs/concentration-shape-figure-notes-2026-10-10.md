# Concentration-shape comparison figures (2026-10-10)

The figure bundle in
`outputs/concentration-shape-comparisons-2026-10-10/` separates measured
spatial evidence, development diagnostics and independent validation.  All
concentrations are in vol% H₂.  Model residuals are in volume-percentage
points (vol%p), not percent relative error.

## Figure 1 — FFI Test 4 footprint

The three panels compare the measured 275 s concentration, DEGALI's
measured-wind replay and the retained SLABx-LH2 research path.  Each horizontal
position uses the maximum across the three instrument heights.  Filled
contours are linear interpolation only inside the sparse sensor convex hull;
they are a visual observation operator, not a CFD field.  The red line is the
4 vol% contour.  FFI Test 4 was used during development, source mappings differ
between the models, and the figure must not be described as an independent
model ranking.

## Figure 2 — HyTunnel measured plume shapes

The four panels use the same public GH₂ source condition: 4 mm nozzle and
nominal 1 g/s.  They show the published average matrices for still air,
1.5 m/s co-flow, 3.5 m/s co-flow and 3.5 m/s counter-flow.  The red/white line
is the measured 4 vol% contour.  Small negative instrument offsets are
preserved in the raw audit, while the colour display begins at physical zero.
This is GH₂ ventilated-jet component evidence and is not quantitative LH₂ field
validation.

## Figure 3 — USN Experiment 20 dynamic cloud

The four times use the exact public common clock and show the room-array side
view.  Sensors sharing a longitudinal coordinate and height are averaged
across lateral offsets solely for this two-dimensional view.  The green line
is 4 vol%.  This is confined GH₂ ventilation evidence for a dynamic
source-to-sensor operator, not an outdoor free-field DEGALI result.

## Figure 4 — vol% error dashboard

The upper-left panel reports MAE and RMSE in vol%p for all 30 FFI Test 4 sensor
rows.  The upper-right panel shows signed bias.  The lower-left panel gives
the raw 4 vol% confusion counts; only three rows are measured positives, so
100% sensitivity is not a broad confidence statement.  The lower-right panel
shows the independent corrected DEGALI residuals for the six FFI arc maxima.

Key values are also written to `concentration_error_summary.csv`, a readable
`concentration_error_report.md`, and the source-hashed `figure_manifest.json`.

At a nominal 4 vol% indication, the available declared specifications imply
an absolute width of ±3.00 vol%p for the Xensor XEN-5320 (±3% of a 0–100 vol%
full scale) and ±0.08 vol%p for the Dräger XXS H₂ HC (±2% of reading).  These
are specification-only bounds, not event calibration errors and not model
residual corrections.

Sensor-accuracy claims are campaign-specific. For the separate USN/FFI open-
channel obstacle data, Henriksen et al. report ±1% for the XEN-5320, and the
manufacturer's preliminary datasheet specifies 1% full-scale inaccuracy.
Keep that attribution distinct from the ±3% full-scale value cited for the
FFI Test 4 dataset above; neither statement substitutes for an event-specific
calibration certificate. If 1% full scale applies over 0–100 vol%, it is ±1
vol%-point (conditionally 50% of a 2 vol% reading and 25% of a 4 vol% reading),
not ±1% relative error at those concentrations.
