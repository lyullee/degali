# Pre-registration: arbitrary-wind obstacle and wall applicability screen

Date frozen: 2026-09-17, before using any obstacle concentration observation
to select a wake or bypass closure.

## Question and scope

DEGALI's independent-energy crosswind path has an unobstructed Gaussian
control volume.  The extension must make an arbitrary horizontal wind bearing
and simple facility geometry usable without representing a cuboid as a CFD
mesh or silently assigning it a dilution coefficient.

The new component accepts a global Cartesian trajectory and either:

- a solid global axis-aligned cuboid; or
- a wall normal to an explicit wind frame.  The wind angle is a **towards**
  bearing, in radians counter-clockwise from global `+x`; it is not a
  meteorological *from* bearing.

Its only model result is a geometric applicability decision.  A closed
line-segment/solid intersection is computed exactly with a slab test.  Face
contact and zero-thickness wall contact count as an intersection.

## Frozen physical rule

If the calculated **centre trajectory** enters a declared solid, the
unobstructed integral field is marked inapplicable downstream of first
contact.  The screen returns the entry/exit position and travelled arclength,
but it changes none of the source mass, H2 mass, momentum, enthalpy,
entrainment, width, temperature or concentration.

For a laterally unbounded transverse wall, centreline contact also returns the
wall-top height as the minimum geometric height a continuous centre path would
need to clear the wall.  This is not a predicted rise: a real cloud may
separate, recirculate, exchange heat, or have finite lateral bypasses outside
the idealisation.  For a finite wall no over/around route is selected.

## Explicit non-claims

This extension does **not** calculate:

- building wake velocity, recirculation, or pressure loss;
- scalar reflection, bypass split, wake dilution or concentration enhancement;
- wall heat transfer, frost/condensate deposition, or obstacle-induced phase
  change; or
- blockage effects on the near-field source.

Thus clearing the centreline does not demonstrate that the finite Gaussian
envelope clears a structure, and an intercepted line must not be reported as
an obstacle-corrected concentration prediction.  Those questions require
obstacle-resolved validation data and a separately declared transport closure.

## Required verification and evidence before a quantitative extension

Unit checks must show rigid horizontal rotation/translation covariance,
correct first contact, exact wall-plane contact, continuity of a multi-segment
encounter, and no mutation of the input plume state.  These are geometry
checks, not validation of obstacle dispersion.

A future wake model must be pre-registered and tested against cases with, at a
minimum, measured obstacle dimensions, wind vector at release height,
upstream/receptor H2 concentration, temperature, and enough time-resolved
locations to distinguish vertical overpass from lateral bypass.  No such
external obstacle dataset is used or distributed by this screen.
