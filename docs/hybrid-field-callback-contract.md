# Coupled field callback contract

Each coarse and patch field callback used by
`HybridAMRExecutor.run_bidirectional` must perform these operations in order:

1. Accept `(state, elapsed_s, dt_s, incoming_mass_kg)`.
2. Apply `incoming_mass_kg` to a declared interface-cell mask with
   `apply_interface_mass_to_field`; do not add it to an arbitrary concentration
   cell or rescale the whole field.
3. Advance the local momentum/scalar/thermal operator for exactly `dt_s`.
4. Integrate the local patch/coarse interface face flux and return its outgoing
   mass for that same step.
5. Return `(new_state, outgoing_mass_kg)`.

The executor clips outgoing mass to the donor inventory and passes the
equal-and-opposite amount to the other callback on the next substep. A solver
that only returns a post-hoc residual or an independently completed field does
not satisfy this contract and must remain conditional.

For the resolved thermal path, the same shared-clock step also supplies a
conserved sensible-enthalpy callback.  The low-Mach corrector forwards the
species and enthalpy face fluxes on its final iteration, so a thermal patch is
not silently disconnected from the bridge.  Momentum packets use the
characteristic boundary mode when requested: the incoming carrier velocity is
reconstructed as `p_in/m_in`, relaxed onto the interface cells, and applied to
the normal staggered face.  The legacy additive impulse remains available only
for compatibility and is not the production-shaped FFI path.
