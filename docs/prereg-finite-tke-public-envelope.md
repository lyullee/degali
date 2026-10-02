# Public finite-TKE end-member protocol

Date: 2026-09-12

## Purpose and scope

This protocol advances the optional finite-TKE branch without pretending that
an unmeasured cryogenic-hydrogen dissipation rate has been recovered. It
defines only a reproducible **rejection/sensitivity** study for the Sandia
5 bar, 50 K, 1 mm free jet. It does not alter the default DEGALI prediction,
fit a thermal result, or make the finite-TKE branch a validated closure.

## Public evidence

Saini et al., [*Effect of storage conditions on the characteristics of
cryogenic hydrogen jet dispersion*](https://doi.org/10.1016/j.ijhydene.2024.04.040),
solve cryogenic 3--7 bar, 50--70 K releases with a compressible RANS
`k-epsilon` formulation. Their stated conventional constants include
`C_mu=0.09`, `sigma_k=1`, and the standard production/dissipation transport
structure. This establishes nomenclature and a comparison model; it does
**not** validate those constants for DEGALI or supply an experimental
epsilon profile.

Mohammadpour et al., [*An efficient dimensionality reduction approach for
modelling cryogenic hydrogen release*](https://doi.org/10.1016/j.ijhydene.2024.10.182),
provide the only located open cryogenic 5 bar/50 K velocity-RMS contour. As
documented in [the public evidence audit](public-turbulence-evidence-audit.md),
it is LES output, not an experiment. It supplies a lower TKE screen only.

No public matched `epsilon`, integral-length-scale profile, or full Reynolds
stress tensor was located. Consequently no numeric end-member range is
encoded in software.

## Inputs that must be declared for every candidate

At every comparison location, a candidate must disclose all of:

1. `k` and the integral length `L` used for dissipation;
2. the dimensionless dissipation coefficient `C_e`;
3. eddy viscosity `nu_t` and the TKE transport ratio `sigma_k`; and
4. the resulting fields

```text
epsilon = C_e k^(3/2)/L
tau     = L/(C_e sqrt(k))
chi_k   = nu_t/sigma_k .
```

`degali.addons.tke_parameterization` implements these identities but has no
default values. The familiar conversion `C_e=C_mu^(3/4)` is supplied only as
algebra when a user explicitly elects a `k-epsilon` comparison closure. It
must not be read as a transferability claim across the shock-containing
cryogenic near field.

## Required filters and reporting

Each end-member run must pass all of the following before its thermal/species
residuals may be compared:

1. positive finite `L`, `k`, `C_e`, `nu_t`, and `sigma_k` everywhere;
2. `k` no lower than an observed-stress PSD bound at matched locations: either
   axial RMS plus its two shear covariances, or axial and one transverse RMS;
3. the existing finite-TKE covariance and total-energy ledgers; and
4. simultaneous, sensor-identical reporting of temperature and hydrogen
   concentration residuals.

Report the whole declared set of end members, including failures. A result
that improves one thermal quantity after choosing `tau` from that same
quantity is excluded. A result that fails a TKE lower bound or a conservation
ledger is excluded rather than clipped or reinitialised.

`tke_realizability_margin(k, u_rms, R_sperp)` and
`tke_two_normal_rms_realizability_margin(k, u_rms, v_rms)` are executable
forms of filter 2. Each returns the unmodified `k`, its sharp PSD lower bound
and their signed difference. A negative margin is a rejected candidate state;
neither helper raises `k`, completes an isotropic tensor or selects epsilon.

`FiniteTkeModalTransport` accepts exactly one explicit observation pair:
`axial_rms` plus `axial_shear_covariance`, or `axial_rms` plus
`transverse_rms`. Supplying an incomplete pair or both pairs is rejected;
omitting all three retains the pre-existing unobserved-TKE research path. The
two-normal path deliberately leaves the other transverse normal stress and all
covariances unmeasured. Every RMS input must also state
`velocity_rms_provenance='gas_velocity'`; the interface rejects a
particle-tracer declaration or an omitted declaration. This is a
realizability gate, not an RMS-to-epsilon conversion or a default LH2
turbulence model.

## Decision after the sensitivity

If one or more *predeclared* end members improves both observed fields while
passing the filters, it remains an opt-in sensitivity result. Promotion to a
physical closure still requires an independent velocity-fluctuation/
dissipation dataset. If none survive, the next physics target is not another
arbitrary TKE constant: it is pressure--strain/normal-stress transport or a
measured thermal/density-profile boundary.
