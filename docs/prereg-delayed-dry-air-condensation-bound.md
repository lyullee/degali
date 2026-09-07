# Pre-registration: delayed dry-air condensation bound

Date: 2026-09-05

## Why this is a distinct missing mechanism

The transported 1/10/100 um particle and stationary-condensate audits begin
after condensed N2/O2 exists. They bound particle carry, dropout and
re-evaporation; they do not test a finite delay in the onset of nucleation.
The trial-10 momentum budget now shows that the phase-equilibrium research
branch is already 64.3 K and weakly buoyant at 10D, whereas the corrected
JETPLU comparison remains 30.7 K and negatively buoyant. The resulting
candidate cumulative buoyancy is 5.695 N by 6 m versus 2.734 N.

Non-equilibrium cryogenic-nitrogen nozzle studies report a finite Wilson-point
delay, approximately 9--11 K supercooling, and separate sensitivity of droplet
growth/thermal relaxation after onset (Sun et al., *Cryogenics* 68, 2015,
DOI 10.1016/j.cryogenics.2015.01.010; Sun et al., *International Journal of
Multiphase Flow* 121, 2019, DOI 10.1016/j.ijmultiphaseflow.2019.103118).
Those nozzle data are not transferable kinetics for an atmospheric LH2 jet,
but they reject instantaneous equilibrium as the only physically admissible
limit. Sun et al. (*International Journal of Hydrogen Energy* 50, 2024,
DOI 10.1016/j.ijhydene.2023.06.201) also identify the opposing effects of
condensed mass and released latent heat in an LH2 release.

## Frozen bound

- Add `equilibrium_dry_air_condensation=False` as an explicit research-only
  metastable-gas bound. N2, O2 and Ar remain gaseous and release no phase
  latent heat. Ambient water retains the existing equilibrium ice treatment.
- Preserve total mass, component enthalpy, density state, H2 fraction,
  velocity profile, source, co-flow, entrainment, force, 10D handoff and all
  numerical tolerances. Do not replace the gas by vacuum or remove its heat
  capacity.
- Keep `True` as the general API default and the validated Raman/default
  configuration. The false value must emit a research warning and appear in
  result metadata.
- First evaluate trial 10 at 10D and 1.78/4/6 m. If the branch remains
  positive and passes the five-flux interface, run the frozen seven trials,
  42 concentration arcs and 17 vertical sections with Li enthalpy transport,
  corrected Houf velocity width and the 0.02 m step.

## Gates and decision fixed before calculation

1. Conservation gates remain `1e-8` at the interface, `1e-5` quadrature and
   `1e-6` downstream; all seven interfaces and 42/17 samples must survive.
2. The 10D candidate temperature, density, mass flux and buoyancy sign are
   reported before the field score.
3. A useful direction requires smaller internal centre MAE without worsening
   `abs(log MG)`, VG, FAC2 or internal width-ratio error relative to the
   corrected equilibrium candidate.
4. Complete suppression of dry-air condensation is a model-form bound, not a
   physical replacement. It is never promoted directly, even if it passes the
   field direction. A successful bound licenses a finite nucleation/growth
   state based on independent kinetics; it does not license fitting a delay
   length to PRESLHY.
5. If the bound fails or increases rise, nucleation delay is rejected as the
   next error-reduction route and the larger pre-10D mass entrainment becomes
   the next isolated target.

## Result

The coefficient-free bound was applied to trial 10 first, exactly as frozen.
The native near-field integration remains conservative: the largest of the
five normalized flux residuals is `3.51e-16`, energy-quadrature refinement is
`4.50e-6`, and the H2-width mismatch is 4.704%.  However, the projected
single-section state is 2.222 K colder than the transferred near-field centre
(`63.337` versus `65.559 K`).  This exceeds the unchanged 2 K interface gate,
so the interface is rejected and the seven-trial field run was not performed.

The force check points in the same direction.  The near-field section remains
weakly buoyant at `0.01792 N/m`; projection gives `0.01496 N/m`, a 16.50%
difference that also exceeds the earlier 10% diagnostic screen.  By comparison,
the equilibrium-air trial-10 near-field force was `0.01848 N/m`.  Completely
suppressing N2/O2/Ar condensation therefore reduces the 10D force only
slightly and changes the transferred centre temperature by only about 0.42 K.
It cannot explain the approximately 34 K difference between the independent
near field and the corrected JETPLU baseline at 10D.

Per the pre-registration, no tolerance was relaxed and no partial 7-trial or
field score was issued.  Delayed dry-air nucleation is rejected as the next
error-reduction route.  The off-default switch is retained only as an explicit
metastable-gas model-form bound; the next isolated target is the larger total
mass entrainment accumulated before 10D.

Machine-readable values are in
`reference/preslhy/delayed_dry_air_condensation_bound_2026-09-05.json`.
