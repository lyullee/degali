"""NASA White Sands liquid hydrogen spills, transcribed from the paper.

Witcofski, R.D. and Chirivella, J.E. (1984), "Experimental and analytical
analyses of the mechanisms governing the dispersion of flammable clouds formed
by liquid hydrogen spills", *Int. J. Hydrogen Energy* **9**(5), 425-435.

Transcribed from the printed tables, not digitised from figures.
"""

#: Table 1: the seven spills. 5.7 m3 of liquid hydrogen except test 7, 2.8 m3.
#:
#: **Only tests 2, 4, 5 and 6 carry concentration data.** Tests 3 and 7 were
#: run at deliberately slow spill rates "to gain insight in the role of spill
#: rate", and Table 1 records test 3 as yielding motion picture only. Neither
#: appears in Table 3 or Table 4.
#:
#: That matters because the release *conditions* for tests 3 and 7 are
#: published -- they appear in Mack's Table 2 alongside the others -- and it
#: is easy to mistake conditions for measurements and believe the lift-off
#: sample can go from four to six. It cannot.
SPILLS = {
    #        spill_s  wind_ms       tamb_C  rh_%  concentration data
    1: dict(spill=60, wind=(2.7, 3.1), tamb=30, rh=18, data=False),
    2: dict(spill=40, wind=(1.3, 1.8), tamb=24, rh=49, data=True),
    3: dict(spill=85, wind=(4.5, 4.5), tamb=26, rh=27, data=False),
    4: dict(spill=33, wind=(3.1, 3.6), tamb=15, rh=43, data=True),
    5: dict(spill=24, wind=(6.3, 6.3), tamb=12, rh=43, data=True),
    6: dict(spill=35, wind=(2.2, 2.2), tamb=15, rh=29, data=True),
    7: dict(spill=240, wind=(3.1, 3.1), tamb=17, rh=29, data=False),
}

#: Table 3: maximum hydrogen concentration at the furthest tower row, vol %.
#:
#: Keyed by test, then height in metres, to ``(sample_bottles, sensors)``.
#: ``None`` where no bottle was mounted -- there were none below 9.4 m on the
#: last row.
#:
#: **The bottles are the accurate source.** The paper is explicit that sensor
#: readings above 8 % "are meaningless" and above 4-6 % questionable, so a
#: sensor value of 6.4 or 8.0 carries no information beyond "more than the
#: instrument can read".
#:
#: This is a **vertical profile at a fixed distance**, which is the quantity
#: that constrains a trajectory. The package previously used only Table 4,
#: a single height per test.
FAR_TOWER = {
    2: {1.0: (None, 0.2), 9.4: (0.0, 1.8), 18.6: (0.0, 6.4)},
    4: {1.0: (None, 0.2), 9.4: (4.2, 7.4), 18.6: (0.5, 4.6)},
    5: {1.0: (None, 7.5), 9.4: (29.2, 6.5), 18.6: (0.0, 0.0)},
    6: {1.0: (None, 2.7), 9.4: (18.7, 8.0), 18.6: (19.0, 4.5)},
}

#: Distance of the furthest tower row from the pond centre, m.
FAR_TOWER_DISTANCE = 33.8

#: Table 4: minimum height at which a flammable cloud was measured, m, with
#: the concentration there. The height is what the package already compares
#: against; the concentration is new here.
LOWEST_FLAMMABLE = {
    2: (18.3, 6.4),
    4: (6.4, 14.9),
    5: (0.3, 17.2),
    6: (3.4, 20.0),
}

#: Evaporation rates and conditions as used by Mack and co-workers, kg/s.
#: These are evaporation rates, not spill rates: 5.7 m3 over the Table 1 spill
#: time gives 10 to 17 kg/s, which is the liquid leaving the dewar rather than
#: the vapour entering the cloud.
RATES = {2: 9.23, 3: 4.23, 4: 10.29, 5: 9.95, 6: 9.48, 7: 1.66}

#: Diameter of the spill pond, m.
POND = 9.1
