# Public turbulence evidence audit for the finite-TKE branch

Date: 2026-09-10

## Decision

Do not wait for an author-supplied file before continuing the bounded
finite-TKE work. Public sources are sufficient to add two independent
constraints now:

1. Sandia's cryogenic-hydrogen Raman study supplies concentration RMS,
   intermittency and PDF-shape evidence for a future scalar-fluctuation
   diagnostic.
2. Mohammadpour et al.'s open institutional manuscript supplies a published
   LES contour of mean axial velocity and axial-velocity RMS for the Sandia
   5 bar, 50 K, 1 mm case. It can bound one normal Reynolds stress and hence
   the minimum realizable TKE.

Neither source identifies epsilon or a TKE dissipation time. The finite-TKE
transport remains opt-in and no default closure coefficient is promoted.

## Sources and provenance

### Sandia scalar fluctuation experiment

E. S. Hecht et al., "Concentration fluctuations and flammability of
cryo-compressed hydrogen and methane jets," *Fuel* 358 (2024) 130230,
doi:[10.1016/j.fuel.2023.130230](https://doi.org/10.1016/j.fuel.2023.130230).
The accepted manuscript is public through
[OSTI record 2311268](https://www.osti.gov/biblio/2311268).

The ignored local source has SHA-256
`0BFB93132FEDCDC15598B4A7BFA0024104A1E460C8DED0359DBD3CCE93EC947C`.
The paper describes a 1.25 mm, 3 bar, 51 K hydrogen jet with 334 retained
Raman images. It reports approximately 0.22--0.23 centreline relative
concentration RMS, a study value `K=0.049`, and the intermittency relation

```text
gamma = (K + 1) / (1 + (f_rms / f_mean)^2).
```

The centreline PDF is approximately Gaussian, while the radial shear layer
becomes bimodal. These observations matter for probabilistic flammability
and plume-edge occupancy. They do not determine velocity TKE or its
dissipation.

### Public LES velocity envelope

J. Mohammadpour et al., "An efficient dimensionality reduction approach for
modelling cryogenic hydrogen release," *International Journal of Hydrogen
Energy* 91 (2024) 649--658,
doi:[10.1016/j.ijhydene.2024.10.182](https://doi.org/10.1016/j.ijhydene.2024.10.182).
The author manuscript is publicly available from
[Macquarie University](https://research-management.mq.edu.au/ws/portalfiles/portal/399409510/394825139.pdf).

The ignored local source has SHA-256
`CE5677815FFAF43B54CDFEB7702C5E27CBBD58BCD1687900B44B2B8DA2D0CECF`.
It models the Sandia 5 bar, 50 K, 1 mm round-nozzle case with a 1-million-cell
LES using a dynamic kinetic-energy SGS model. Figure 3 shows mean and RMS
axial velocity in the first approximately 64 mm. The paper explicitly says
that Sandia's published report did not include velocity distributions.

Consequently this figure is **not** an experimental validation target. It is
used only as an independently published, physically plausible screening
envelope for the reduced transport state.

### Sandia PIV/Raman status

Sandia's public capability and project presentations document simultaneous
Raman concentration/temperature imaging and PIV. The 2019 flow-measurement
presentation includes mean axial-velocity comparisons and lists several PIV
cases, including 43 K/4 bar, 59 K/3 bar, 56 K/3.5 bar and 80 K/3 bar through
a 1 mm nozzle. The 2021 project presentation states that experiment data are
available on request. A public numeric archive containing the PIV arrays,
velocity RMS, Reynolds stresses, `k`, or epsilon was not located in OSTI,
DOE Data Explorer, or the Sandia publication records searched for this audit.

This absence is why the LES material is retained as a bound rather than
mistaken for a measurement.

### Follow-up public-archive search (2026-09-13)

A second search of OSTI, DOE data pages and Sandia's cryogenic-release
capability pages found confirmation that simultaneous Raman/PIV is used, but
no downloadable numeric PIV fields, velocity-RMS arrays, Reynolds stresses,
epsilon or particle-size archive. Sandia's public capability page describes
the diagnostics and their model-validation role; it does not offer a data
download. The public search therefore adds no calibration datum and does not
alter this branch's no-default decision.

## Registered figure reduction

The published Figure 3 was extracted from page 5 without editing. Printed
axis ticks and colour-bar ticks define a linear pixel-to-value mapping. Values
were sampled over a 5 x 5 pixel neighbourhood; shear-layer maxima were found
on the same image row and mapped to the corresponding radial coordinate in
the mean-velocity panel. Reporting precision is limited to about 10 m/s
because of rasterisation, colour quantisation and contour interpolation.

Approximate readings are:

| z (m) | centre mean U (m/s) | centre axial RMS (m/s) | shear-layer axial RMS (m/s) | mean U at RMS peak (m/s) | RMS/U | weakest k / (U^2/2) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.010 | 630 | 20 | 120 | 370 | 0.33 | 0.11 |
| 0.020 | 650 | 80 | 180 | 510 | 0.36 | 0.13 |
| 0.030 | 470 | 160 | 170 | 430 | 0.39 | 0.15 |
| 0.040 | 270 | 110 | 110 | 270 | 0.40 | 0.16 |
| 0.050 | 170 | 60 | 70 | 170 | 0.39 | 0.15 |
| 0.060 | 120 | 40 | 40 | 110 | 0.37 | 0.14 |

The non-monotone first two centreline velocities are consistent with the
under-expanded shock train described in the paper. These near-nozzle points
must not be treated as a smooth far-field decay law.

### Additional two-normal-RMS plot screen (2026-09-16)

Li et al., [*Flow and turbulence characteristics of underexpanded cryogenic
hydrogen jets*](https://hysafe.info/uploads/papers/2023/207.pdf), publish
centreline axial and radial velocity-RMS plots (their Figures 10 and 11) for
six **LES** cases. The 5 MPa, 50 K, 1.5 mm case is closest to the cold,
high-pressure source envelope considered here. At `z/d >= 90`, a direct
reading of the printed curves gives approximately

```text
u_rms / U_cl = 0.22--0.24
v_rms / U_cl = 0.18--0.21.
```

The paper separately states that the axial component exceeds the radial one;
the plot is therefore internally consistent. It also reports the far-field
radial-RMS comparison level `v_rms/U_cl = 0.23` from a different, small-scale
hydrogen experiment. Neither item supplies a matched covariance, azimuthal
normal stress, dissipation, or integral length. They are consequently not
measurements for calibration and not an initial condition for DEGALI.

They do give a second, coefficient-free realizability screen. For two
observed normal variances,

```text
k >= 0.5 * (u_rms^2 + v_rms^2),
2 k / U_cl^2 >= (u_rms/U_cl)^2 + (v_rms/U_cl)^2 ~= 0.08--0.10.
```

The omitted azimuthal variance and covariances have their smallest
realizable values in this *lower bound*; no isotropy is imposed. The new
`minimum_tke_from_two_normal_rms` and
`tke_two_normal_rms_realizability_margin` helpers implement only this reject
test. They do not ingest this digitisation, choose a coefficient, or change
the default prediction.

### Official Sandia PIV presentation re-audit (2026-09-17)

The public Sandia conference presentation, *Experimental validation of a
model for cryogenic hydrogen jet dispersion* (Hecht and Panda,
SAND2018-2834C; [OSTI 1503433](https://www.osti.gov/biblio/1503433)), was
read directly after the archive search. It confirms Raman imaging together
with particle PIV and shows mean-velocity comparison panels for the 1 mm
campaign. It does **not** tabulate or archive velocity RMS, Reynolds
covariances, sampling interval/count, integral scale, or epsilon. Its PIV
tracers are described as condensed entrained moisture particles, so even a
future RMS extracted from its display could not silently be treated as a
gas-phase TKE measurement without a tracer-response assessment.

Accordingly, the presentation adds no numeric closure input. It strengthens
the provenance boundary: only a declared gas-velocity observation may enter
either `FiniteTkeModalTransport` realizability gate. The new two-normal-RMS
gate is an optional, stricter algebraic check for a source that explicitly
provides both gas velocity RMS components; it is not a permission to convert
the Sandia mean-PIV panels or the LES figures into a calibrated TKE field.

## Physics admitted into DEGALI

For axial RMS `u_rms`, the known normal stress is

```text
R_ss = u_rms^2.
```

For the two axial/transverse shear covariances collected in vector `r`, a
positive-semidefinite Reynolds-stress tensor requires the sharp lower bound

```text
k_min = 0.5 * (R_ss + dot(r, r) / R_ss),  R_ss > 0.
```

If shear is not available, setting `r=0` gives only the weakest bound
`k_min=R_ss/2`. Relative to local mean axial kinetic energy, this becomes
`k_min/(U^2/2)=(u_rms/U)^2`, which produces the final column above.

`minimum_tke_from_axial_rms` implements this algebraic constraint without an
isotropy assumption and without clipping an impossible covariance. It does
not create an epsilon law, choose a normal-stress ratio, or convert the LES
figure into an initial condition automatically.

## Consequence for the next run

The bounded finite-TKE sensitivity can now reject trajectories whose local
TKE falls below the axial-RMS realizability envelope. It still must span a
declared dissipation-time range because no public matched epsilon profile was
found. A candidate may advance only if it:

1. remains above the registered TKE lower bound where the 5 bar/50 K case is
   compared;
2. satisfies covariance PSD and the existing energy ledger;
3. improves temperature and hydrogen metrics simultaneously at identical
   sensors; and
4. is labelled a sensitivity, not a validated closure, until independent
   velocity-fluctuation data are available.

No third-party source PDF, figure, or digitised row is distributed in the
repository.
