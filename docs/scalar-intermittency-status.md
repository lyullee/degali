# Scalar intermittency observation status

Date: 2026-09-12

## Reason for this branch

At matched PRESLHY stations, a steady mean prediction can lie between the
observed time-window minimum and median. That difference must not be hidden by
retuning a mean enthalpy or diffusivity. It can arise from unresolved scalar
fluctuations, plume-edge occupancy and sensor response, none of which is a
resolved steady-state transport coefficient.

The public Sandia Raman study by [Hecht et al.](https://doi.org/10.1016/j.fuel.2023.130230)
reports, for its 1.25 mm, 3 bar, 51 K hydrogen jet, centreline relative
concentration RMS of roughly 0.22--0.23 and the reported relation

```text
gamma = (K + 1) / (1 + (f_rms/f_mean)^2),  K = 0.049.
```

The centreline PDF was approximately Gaussian while the shear layer was
bimodal. This is the correct motivation for keeping mean, RMS and occupancy
separate; it is not evidence that the formula is universal across the PRESLHY
trials or all LH2 releases.

Using the paper's rounded `0.22` and `K=.049` literally gives `gamma` just
above one, while rounded `0.23` gives a value just below one. DEGALI retains
this small incompatibility rather than clipping it; the unrounded experimental
values are not publicly tabulated.

## Implementation

`validation.scalar_fluctuations` evaluates the published relation verbatim,
reports whether the resulting occupancy is a valid probability, and never
clips it. It also exposes the identity

```text
mean scalar = occupancy * conditional active-state mean.
```

The latter is labelled a two-state reporting envelope only. Neither function
adds random forcing, changes a plume mean, or estimates RMS from a DEGALI
solution. A physically impossible occupancy is an incompatibility result, not
a reason to tune `K`.

## Use gate

This diagnostic may be used only when a comparison supplies scalar mean and
RMS from the same spatial location and sampling window. It cannot turn
PRESLHY minimum/median temperature records into a velocity-TKE boundary, nor
can it justify modifying the thermal-width or source-strength model. A future
stochastic receptor extension needs an independently validated, location-aware
RMS/occupancy field.

## Public-data search update

An open 2026 LES/POD/BiLSTM paper by Mohammadpour et al.,
[*Innovative approaches for predicting cryogenic hydrogen behaviour*]
(https://doi.org/10.1016/j.ijheatfluidflow.2025.110025), uses 540 simulated
snapshots of the same 5 bar, 50 K Sandia case and presents instantaneous
hydrogen/temperature fields. Its numerical snapshots are stated to be
available on request, not in a public archive. The reported 2.85% flammable
threshold deviation is against its LES reconstruction target, while the paper
reports about 27% mean deviation of LES from the experimental 0.04-mole-
fraction contour. It therefore cannot be treated as public experimental
RMS/occupancy data or as a DEGALI calibration target.

The public search found no downloadable, location-resolved cryogenic-H2
RMS/PDF dataset. The closest fully downloadable scalar/velocity archives are
non-cryogenic reacting H2/He jet-flame datasets, whose thermochemistry,
density ratio and geometry do not justify transfer to the LH2 branch.
