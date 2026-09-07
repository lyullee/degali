# Finite turbulent-energy modal transport: operator construction specification

2026-09-06. Implementation specification, NOT a completed real-trial model.
This follows the shear-covariance lower-bound diagnosis. It does not select
an actual inlet k, dissipation time, turbulent Prandtl/Schmidt number or a
circulation law. Those inputs must remain explicit, not inferred from a
desired concentration or temperature score.

## Variables and scope

Keep the existing enriched thermodynamic C/H fields, phase table and reduced
axial velocity/force conventions. Introduce Q=rho*k as a distinct positive
field with one log-center variable and the existing center-zero square modes.
Use Q=Q0 exp(-kappa_u*q+psi.pQ) as a numerical basis reference; the base width
is not an independently evolving gauge. Q is NOT included in the EOS thermal
enthalpy. The total axial energy flux is H-flux + mean axial KE + A integral(Q*u).

This is still a reduced mean-momentum model: no full Reynolds tensor, pressure
transport, buoyancy-density correlations or compressible k-energy EOS is being
claimed. Preserve the original model defaults and flag these approximations.
The existing PSD k lower bound becomes a necessary diagnostic, not a clipping
rule. A caller must provide nonnegative ambient k, positive transverse k mixing,
and a positive finite dissipation time field tau. D=A*Q/tau is a declared
one-equation relaxation assumption, not a calibrated or validated LH2 law.

## Coupled conservation rows to implement

Let the old reconstructed fluxes be FM,FP and shear production
P=-grad(u).(FP-u FM). They are affine in the current mean-state rates.

1. Retain mass, hydrogen and both momentum global rows. Replace the global
   energy Jacobian by JH+JmeanKE+JQ; its source adds k_ambient times the
   incoming ambient mass source to the existing reduced energy source.
2. Keep every hydrogen weak row. Every thermal weak row gets D instead of
   P. Do not retain immediate shear heating while also storing P in Q.
3. Add Q transport tested against [1,psi]. Its volume flux is
   k FM - A*rho*chi_k/(2q0)*grad(k); its source is P-D. Use the same moving
   transverse coordinates and quadrature factors as the other scalar rows.
4. The natural boundary Q flux is k_ambient FM. Check its agreement with the
   constitutive gradient flux independently at faces. This is an additional
   boundary requirement, not automatically satisfied by a chosen Q profile.
5. Since boundary Q is k_ambient FM, the thermal reservoir boundary remains
   .5*Ua² FM - (u FP-.5*u² FM). The combined boundary energy is then
   (.5*Ua²+k_ambient) FM, with no double-counted turbulent inflow.

For m square modes there are m+1 new Q rates and m+1 new weak equations.
The global total-energy and Q zeroth rows together must imply the thermal
zeroth budget. Derive/check this identity, including face terms. Do not omit
the old thermal or species rows to force a square system.

## Required verification before any physical trajectory

- Positive finite Q/tau/mixing and explicit ambient k; reject missing inputs.
- Independent actual-value Q advective derivatives on fixed fields.
- P cancels between mean KE and Q, D cancels between Q and thermal energy,
  and the natural face ledger sums to the declared ambient total energy.
- Complete matrix rank/conditioning, refinement and independent full-energy
  differences on coherent exact geometry, without recycling an old tangent.
- Independently check Q face-gradient agreement, k>=shear covariance bound,
  nonnegative production and all earlier scalar/momentum/physical gates.
- Recover the immediate-heating energy identity for a manufactured steady
  Q case with P=D and zero net Q transport. This is an algebraic limiting
  test, not proof that arbitrary real fields have that equilibrium limit.

First build and test the operator without fitting a real initial Q. A real
trial requires jointly consistent Q initial shape, energy/momentum source
matching, and dissipation/mixing inputs. A generic tau chosen for convenience
must not be promoted to a physical closure or measured-accuracy result.
