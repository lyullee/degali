# PRESLHY finite-rate LH2 droplets: frozen protocol

Status: pre-registered after the measured-source ablations and before any
two-phase field score was generated.

## Why this route is next

The seven-trial measured-source candidate conserved mass, momentum and energy
at every accepted interface, but worsened all field promotion metrics.  Its
two ablations did not isolate an admissible single-phase change.  `flow_only`
passed four of seven interfaces and failed the 2 K centre-temperature gate in
trials 12, 23 and 24.  `nozzle_only` passed six; trial 23 reached its
momentum-balance handoff at 4.79 source diameters, before the axisymmetric jet
had entered the self-similar development zone.  Neither failed ablation is
eligible for field scoring or per-trial selection.

The source thermodynamics show the omitted state explicitly: decompression of
the measured nozzle states leaves most hydrogen liquid immediately after the
pressure flash.  The former source instead moved directly to an all-vapour
evaporation endpoint.  That shortcut transfers liquid mass into gas
concentration before the required heat has entered the jet.  The next
candidate therefore keeps liquid and vapour hydrogen separate in space.

## Frozen population and topology

The validation population remains trials 10, 11, 12, 22, 23, 24 and 25 with
the already frozen windows in `reference/preslhy/e35_reduced.json`.

- Nozzle trials 11, 12, 23 and 24 use the pressure-loss flow rate, corrected
  TC3 temperature and measured PT2 pressure in
  `reference/preslhy/measured_pipe_source_2026-09-05.json`.
- Open-pipe trials 10, 22 and 25 retain the preceding source construction
  until outlet quality or density is independently available.  Their PT2
  readings are effectively atmospheric, so temperature and pressure alone do
  not identify their two-phase state.
- No trial may be excluded, reassigned to a source mode or given a separate
  coefficient because of downstream agreement.

## Frozen post-flash source equations

At the measured pipe plane the area is `A_e`, velocity is
`u_e = m_dot / (rho_e A_e)`, and total axial momentum includes pressure
thrust.  The atmospheric plane is defined by

`u_f = u_e + (P_e - P_a) A_e / m_dot`.

Post-flash enthalpy is fixed by total specific-energy conservation,

`h_f = h_e + (u_e^2 - u_f^2) / 2`,

and atmospheric saturation enthalpies determine vapour quality.  The liquid
and vapour rates are `(1-Q_f)m_dot` and `Q_f m_dot`; they are not recombined
until finite-rate evaporation actually transfers mass.

The Yellow Book atomisation correlation is used only to initialize a
monodisperse sensitivity scale.  With post-flash radius `b_f`,

`Re_f = 2 b_f u_f / nu_L`,

`We_f = 2 b_f u_f^2 rho_L / sigma`.

The published intact-jet branch is used where its two stated conditions hold.
Otherwise `d = C_ds sigma / (u_f^2 rho_air)`.  The central value is fixed at
`C_ds = 15`; the only allowed sensitivity values are 10 and 20.  `C_ds` is
never fitted to a receptor, and a conclusion that changes across this range
must be reported as atomisation-uncertain rather than promoted.

## Required two-phase balances

Gas concentration uses vapour hydrogen only.  Total mixture mass, axial and
vertical momentum, total hydrogen, liquid hydrogen and total energy must each
close to relative residual below `1e-8` at the source and every numerical
step.  Evaporation transfers equal and opposite hydrogen mass and latent heat
between the liquid and gas ledgers.  Entrainment may add ambient gas and its
enthalpy but may not create hydrogen.

The first candidate is a no-rainout, common-velocity homogeneous-equilibrium
bound, consistent with the selected horizontal tests showing no release-time
rainout.  Droplet settling and velocity slip are calculated as diagnostics,
not silently discarded.  If particle relaxation or settling is not small on
the local residence-time scale, a second liquid momentum equation is required
before field scoring.

Gas and liquid temperatures remain separate whenever finite-rate heat
transfer is active.  An isolated droplet exposed directly to ambient
temperature may be computed only as a fastest-heating diagnostic; it cannot
be used as the production evaporation law because droplets inside the cold
jet see the local gas temperature and share its finite heat capacity.

## Gates fixed before field scoring

The candidate is admissible only if:

- all four measured nozzle sources have finite, physical phase rates and
  source mass, momentum and energy residuals below `1e-8`;
- `C_ds = 10, 15, 20` is reported without calibration;
- every retained droplet march keeps both phase masses non-negative and all
  six balances below `1e-8`;
- the spatial source reaches a valid crosswind handoff after the relevant jet
  development zone;
- settling and relaxation diagnostics support the common-velocity bound, or
  the model is upgraded to separate liquid momentum before scoring;
- all seven original trials pass their existing interface gates.

Only then is the unchanged paired field score run.  Promotion still requires
concentration MG closer to 1, lower VG and no lower FAC2, vertical width ratio
closer to 1, centre-height MAE no larger, and lower source-temperature median
absolute error and RMSE.  Failure leaves this as a documented physical bound
and does not change the validated default.

## Results after pre-registration

The measured post-flash states retain 98.08--100% liquid hydrogen in nozzle
trials 11, 12, 23 and 24.  At fixed `C_ds=15`, the Yellow Book shattered-jet
correlation gives 0.782--2.011 um droplets.  The published `C_ds=10--20`
sensitivity scales that diameter linearly but, as required, does not change
the collective equilibrium endpoint.  Isolated droplets exposed directly to
ambient air have fastest-transfer lifetimes of 0.35--2.39 us; velocity
relaxation is 0.13--0.89 us and terminal settling speed only
1.3--8.5 um/s.  These diagnostics support the common-velocity, no-rainout
*fast bound*, but not an ambient-temperature production evaporation law.

The independent PRESLHY D3.1/GASFLOW-MPI mass-transfer screen gives the same
regime without using the NASA-fitted Lee coefficient.  Equation 45,
`c = 6 Sh D/d^2`, is implemented directly.  To make `1/c` as long as the
0.178--0.593 ms post-flash transit would require an effective diffusivity of
`2.86e-10--5.68e-10 m2/s` at `C_ds=15`.  NBS Monograph 168 Table 25 reports
`0.8016e-6 m2/s` for H2 self-diffusion at 20.4 K and 0.1 MPa.  Although that
self-diffusion datum is a scale comparison rather than an exact spray
effective diffusivity, the required value is only 0.036--0.071% of it, and a
turbulent contribution would increase rather than decrease transfer.  The
corresponding reference relaxation time is 0.064--0.421 us.  A free phase-
delay coefficient is therefore not admitted: within the published droplet
range it would manufacture a slow phase rather than model one.

The same equation inverted for diameter makes the scale mismatch directly
testable.  At the NBS reference diffusivity, delaying transfer for the
collective-source transit requires 41.4--75.5 um droplets, whereas the
`C_ds=15` correlation gives 0.782--2.011 um.  The required diameter is
37.6--52.9 times larger.  Consequently, a finite-rate two-phase transport
branch is warranted only if a measured large-droplet tail on that scale is
provided; the current micrometre correlation cannot support it.

The collective enthalpy balance requires 0.835--0.868 kg ambient air per kg
H2 to finish evaporation.  Its constant-entrainment formation distance is
0.031--0.064 m.  All four measured nozzle endpoints close hydrogen mass,
total mass, momentum and total energy below `1e-12` relative.  Together with
the unchanged three open-pipe sources, all seven 10D interfaces pass: maximum
five-flux residual `7.39e-16`, H2-width mismatch 2.158% and centre-temperature
mismatch 0.680 K.

The first full field run exposed a numerical rather than physical failure in
trial 23: a single-start nonlinear flux inversion stopped at residual
`5.57e-6`.  A strict multi-start retry using factor-of-two neighbouring width
and velocity guesses, with the same equations and `1e-9` acceptance
tolerance, recovered the state and closed the downstream balance to
`3.83e-8`.  The other six solved trajectories were unchanged.  Disjoint
pooling of those solved shards is exact for the definitions of MG, VG, FAC2
and the arithmetic geometry means.

On all 38 common arcs, baseline MG/VG/FAC2 is
`1.0839 / 1.2344 / 0.9474`; the fast-bound candidate gives
`0.8207 / 1.1974 / 0.9211`.  Scatter improves but mean bias and FAC2 worsen.
Across 17 common vertical fits, width ratio improves from `1.0915` to
`0.9801`, while centre-height MAE worsens from `0.0542` to `0.0700 m`.
Trial 23 alone improves strongly in concentration
(`MG 1.5790 -> 1.0906`, `VG 1.2632 -> 1.0102`) but its centre MAE increases
from 0.0593 to 0.0735 m.

## Decision

The homogeneous-equilibrium fast bound is not promoted.  It demonstrates
that retaining the measured liquid inventory and phase-correct source density
repairs scatter and width, and it materially repairs trial 23 concentration,
but it does not reduce error across all required outputs.  The result also
narrows the remaining physics: micrometre droplets track the gas too closely
for gravitational slip to explain the campaign-scale trajectory residual,
while their collective evaporation remains limited by entrained-jet heat
capacity rather than isolated-droplet transfer.  A future two-temperature
finite-rate model must therefore be tested as a slower bound, not assumed to
be an automatic correction.

Machine-readable results:
`reference/preslhy/measured_pipe_droplet_source_2026-09-05.json`,
`reference/preslhy/measured_pipe_droplet_equilibrium_interface_2026-09-05.json`
and
`reference/preslhy/measured_pipe_droplet_equilibrium_field_complete_2026-09-05.json`.
