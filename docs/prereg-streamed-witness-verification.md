# Streaming restart of independent witness verification

2026-09-06, after the vectorized fine audit was deliberately stopped when
the32GB laptop ran out of physical memory. Its completed coarse result and
partial file remain immutable. This is a computational implementation change,
not changed fields, equations, quadrature, witness, or scientific gates.

Use `stream_retained_operator` with8 angles per batch, accumulating the5
global and2m weak rows, all8 flux columns, RHS and actual moment Jacobian.
The hard boundary-moment equations are irrelevant to a frozen inequality
witness and are not solved. Fixed-witness full vector fluxes, production,
mass/heat boundaries and pointwise defects are evaluated in each batch.
No fields or amplitudes are reoptimized. Peak batch node count is recorded.

Repeat radial4/angular8 and8/16 with the SAME current phase angular/radial
splits and mixing events as the aborted audit. Compare streamed coarse
matrix, RHS and advective Jacobian with the retained dense coarse values
(relative/absolute denominator max(abs(reference),1), <=1e-8). Different
addition order is allowed, different integral definitions are not.

Then finish the original independent matrix/RHS refinement and two-sided
actual-moment FD tests with unchanged1e-8 retained-equation and1e-5
refinement/FD gates. Retain the coarse witness's already detected1.123e-5
weak failure; a successful restart does not relabel that candidate as passed.
Independent65 rays/radial16 remain required. A future refined witness would
need its own registration/result, not overwrite this fixed candidate.

Hash the stopped partial and its106 dependencies plus new implementation,
tests, tool and this specification before/after. Synthetic tests compare
the retained dense/streamed matrices and batch2/6 invariance. No default,
downstream or observation claim. Do not resume the memory-heavy original.
