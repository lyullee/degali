# Same-table root polishing and the energy finite difference

2026-09-06, after the unpolished common-grid pair also failed. Its errors
for h/2,h,2h,4h are6.066e-4,1.202e-4,1.220e-4,4.106e-5. Preserve that
negative result: local subtraction/accurate summation alone was insufficient.

The existing phase inverter stops at an enthalpy residual of
1e-10*max(abs(H),1). Differencing small perturbations of large conserved
energy may amplify root accuracy errors. Test this as a NUMERICAL hypothesis
without changing the thermodynamic table, phase branch or constitutive law.

New verification-only PolishedPhaseMassEnthalpyInverter first obtains the
original accepted root. Perform exactly3 Newton-polishing rounds on the
SAME exact bilinear H(rho,C) table and analytic root slope. Each round also
checks the adjacent representable density in the residual-reducing direction.
Keep a candidate only when its absolute enthalpy residual strictly improves;
retain the old root otherwise. Respect original physical bounds and property
oracle consistency checks. Report maximum original-normalized residual
before/after and improved-point counts. No new root-acceptance relaxation.

The paired actual-field finite difference still does NOT use analytic EOS
derivatives to manufacture a derivative; the slope is used only to solve
the same underlying state-value equation more accurately. Both perturbed
field views explicitly use the verification inverter. Binary64 values and
math.fsum remain; this is not arbitrary precision or a default replacement.

Repeat the same four distances0.5,1,2,4 times the frozen phase-consistent
candidate's original h. Compare all actual5+2m derivatives with its frozen
analytic Jacobian direction. Original h,h/2 and their difference must all
meet1e-5 for a same-step diagnostic pass. Larger distances remain sensitivity
only. A changed oracle, worse accepted root, branch violation or failed
original pair is a failure, not grounds for looser gates.

If successful, this establishes an independent value-difference verification
with better inverse precision; it does not re-label the unpolished failure,
prove stress closure, integrate downstream or improve an observed score.
Freeze the failed paired result and all dependencies plus new code/tests/tool/
specification; refuse overwrite. This is a separate follow-on experiment.
