# Architecture: what is 1960s and what is not, and why

The port is about 11 900 lines. Roughly a tenth of it deliberately keeps
Fortran idiom; the rest is ordinary modern Python. The split is not
accidental and it is not laziness, so it is worth setting out where the line
falls.

## Where the old shape is kept, and why

`core/rkgst.py` is a transcription of the `RKGST` Runge-Kutta-Gill routine,
callbacks and all: a derivative routine that writes into its output array in
place, an output routine that can terminate the run, and a single control
array passed to both.

**That is not a stylistic choice. The integrator is part of the answer.**

`RKGST` estimates its error by taking each step twice — once whole, once as
two halves — and bisecting until they agree. Where the step sequence goes is
what decides which points get recorded, and those points *are* the rows of
DEGADIS's own output. Substituting `scipy.integrate.solve_ivp` gives a more
accurate integration and a different model.

Two places in this repository show that concretely. Guessing `ERBND = 0.005`
instead of 0.0025 made the step sequence diverge at step three, and the ratio
between the two — 2^(−1/4) — is what identified the expansion exponent and led
to the `.ER1` reader. And `ALPH` fits the wind-profile exponent by driving
`ZBRENT` around an `RKGST` quadrature whose own error bound, 0.005, moves the
fitted exponent further than the root finder's tolerance does; integrating it
accurately shifts alpha by 2e-5 and everything downstream with it.

So the derivative routines that `RKGST` calls — `PSS`, `SSG`, `MODEL`, `OB`,
`SRC1` — keep their Fortran signatures too. They have to: the integrator
calls them.

## Where it is not kept

Everything else. `core/thermo.py` is a thousand lines with no raw index
access at all; the three add-ons have none; the readers, the validation
machinery, the drivers and the entry points are written the way any Python
would be:

- `@dataclass` for state — `CloudState`, `JetResult`, `FlashResult`,
  `Assessment`, `Profile`, `Snapshot`
- `Protocol` for the swappable pieces, so `DegadisClosure`,
  `BuoyantClosure` and `UnifiedClosure` are interchangeable and the downwind
  model does not know which it has
- properties for derived quantities — `shape_factor`, `peak_concentration`,
  `contact_fraction`, `mass_between`
- type hints throughout, numpy for the array work, `pyproject.toml`, a CLI,
  a CI matrix, 161 tests

## The control array

The one piece of Fortran shape that leaked into readable code was `PRMT`: a
single list carrying the integration bounds, step, tolerance and stop flag in
its first five slots, and beyond them whatever the caller's derivative routine
wants to stash. `prmt[9]` at a call site tells nobody anything.

It is now wrapped in `rkgst.Control`, which is *still a list* — it subclasses
`MutableSequence`, every index the Fortran uses works, and the ported
routines cannot tell the difference — but the five that matter have names,
and terminating a run is `prmt.halt()` rather than `prmt[4] = 1.0`.

The refactor changed no number anywhere; the 159 parity tests that existed
before it passed unchanged after it, and two more now pin the behaviour.

## The rule

Anything that decides *what number comes out* stays as DEGADIS wrote it, and
is tested against DEGADIS to 1e-12. Anything that decides *how the code reads*
is modern. Where those two conflict — as they do in the control array — the
wrapper preserves the first and fixes the second.
