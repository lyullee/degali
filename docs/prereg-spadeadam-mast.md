# Pre-registration: Spadeadam ventilation-mast comparison

Written 2026-09-03 after recovering DNV GL Report 902696 source tables and
before running a mast-outlet dispersion model.

## Question

Can the existing jet equations represent the cold vertical hydrogen flow from
the ventilation mast without mistaking the supply release inside the TCS for
the atmospheric source?

This is a **separate source class** from both outdoor Spadeadam paths. It does
not enlarge the n=2 horizontal-jet or n=5 ground-impingement samples.

## Included tests

- Include tests 8--14, which have an intact ventilation mast and a DNV steady
  outlet interpretation.
- Exclude test 15: the mast had been destroyed by the Test 14 explosion.
- Use only external field sensors `OC_01`--`OC_24` and `OC_28`--`OC_30`.
- Exclude `OC_25`--`OC_27`, which are beside the TCS low-level vent, and
  `OC_31`--`OC_40`, which are inside the TCS.
- Drop rows already flagged `over_range`.

## Source, fixed before the run

Read the source from `reference/spadeadam/mast_conditions.csv`:

- 450 mm circular outlet;
- centre at east 5.205 m, north 0 m, elevation 11.5 m relative to the common
  TCS-floor origin;
- vertical upward release (`theta = pi/2`);
- test-specific outlet temperature, density, mass flow and averaging window;
- 100 vol % hydrogen, explicitly an assumption made by DNV GL;
- no flashing or equivalent-orifice expansion at the outlet.

Use the corrected plume entrainment coefficient and density-scaled
entrainment, because this is already a developed low-density gas source. Do
not apply ground contact at an 11.5 m outlet.

## Wind and coordinate transform

Use `wind_high` at a 10 m reference height, the closest measurement to the
outlet. Treat the tabulated direction as the meteorological direction **from**
which the wind blows.

Sensor bearings are clockwise from north. For each sensor,

    east  = radius * sin(bearing)
    north = radius * cos(bearing)

subtract the mast outlet `(east=5.205, north=0)`, rotate onto the downwind axis
opposite the reported wind-from direction, and retain only `x_downwind > 0`.
Evaluate the model at the resulting downwind, crosswind and absolute-height
coordinates. Do not maximise over an arc: the actual bearing is now known and
is part of the test.

## Observation and metric

DNV states that small field drift below about 0.5 vol % was not removed. Set
`0.5 vol %` as the field detection threshold.

The primary metric is a negative constraint:

- count sensors where the model predicts at least the 4 vol % LFL while the
  measured peak is below 0.5 vol %;
- report that count and the number of downwind external sensors;
- list each false-positive location rather than averaging it away.

Do not compute MG or VG across nondetects. Log statistics are undefined at
zero, and substituting an arbitrary half-detection value would make the result
depend on that choice.

Secondary outputs are the maximum predicted concentration at any measured
external location and the modelled centre height at each positive downwind
distance. The low-mast wind is a sensitivity run, not the primary result.

## Interpretation fixed in advance

- Zero false-positive flammable sensors is **consistency with a weak negative
  constraint**, not validation: the sparse sensors can miss a narrow plume.
- Any false-positive flammable sensor is a model/data discrepancy to inspect,
  but a single point is not sufficient to tune a coefficient.
- This comparison cannot validate near-source concentration because no sensor
  array crosses the mast plume at its outlet.
- No coefficient will be fitted on these seven tests.

## Results added after the run

Primary high-mast-wind result:

| test | downwind external sensors | observed >=0.5 vol % | false flammable nondetects | maximum predicted at a sensor (vol %) |
|---:|---:|---:|---:|---:|
| 8 | 9 | 0 | 0 | 5.13e-4 |
| 9 | 12 | 5 | 0 | 1.81e-8 |
| 10 | 0 | 0 | 0 | — |
| 11 | 21 | 8 | 0 | 7.07e-9 |
| 12 | 21 | 3 | 0 | 1.46e-8 |
| 13 | 24 | 0 | 0 | 2.41e-12 |
| 14 | 24 | 2 | 0 | 1.48e-7 |
| **total** | **111** | **18** | **0** | |

The primary negative constraint passes: the model does not create a flammable
prediction at any field nondetect. It is nevertheless not a successful
concentration comparison. All 18 detected sensors are predicted below 1e-6
vol %, and ten of those detections lie within 3.5 m crosswind of the calculated
plume axis. At the detected locations the modelled centre is at least 8 m
above the sensor.

The low-mast-wind sensitivity also produces zero false-flammable nondetects
and does not recover the positive signals. The result is therefore another
trajectory constraint: the existing vertical jet rises too far above the
5.24--6.94 m sensors, rather than merely missing them laterally. It is
consistent with a nonflammable external field but cannot reproduce the
measured low-concentration skirt.
