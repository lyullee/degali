"""Comparing degali against field measurements.

Pairing a model with a field trial requires three decisions that are easy to
make implicitly and that change the answer more than most modelling choices
do.  This module makes each one an argument.

**Height.**  Trials instrument several elevations on the same mast.  A
ground-level prediction compared against whichever height read highest looks
good for the wrong reason, so the model's own vertical profile is evaluated at
each sensor's elevation,

.. math:: C(z) = C_c \\exp\\!\\left[-(z/S_z)^{1+\\alpha}\\right]

which matters most in the near field, where :math:`S_z` is a fraction of a
metre and the correction is a factor of several.

The default is to use *every* instrumented height.  Fixing one is available
but does two bad things: it makes the comparison circular, since the vertical
profile then does the work and never gets tested, and it silently discards
trials that happen not to instrument that elevation -- Thorney Island samples
at 0.4 m and FLADIS at 1.5 m, so asking for 1.0 m drops both series entirely
while reporting "no usable measurements".

**Averaging.**  Field data is sampled every second; a model reports a
concentration averaged over its own averaging time.  The two have to be put
on the same footing, and there are two ways to do it: average the
measurements up, or run the model with a short averaging time and compare
against short-time peaks.

The second is used by default, because the first is a trap.  The obvious
averaging time -- the release duration -- is far longer than the time a cloud
takes to pass one sensor, so averaging over it mixes the passage with the
empty record either side and understates the observation badly.  On Burro 9 it
moves the geometric variance from 1.4 to 10.8 and drops every point out of a
factor of two, which would be read as a model failure rather than as an
artefact of the reduction. Pass ``average_measurements=True`` together with an
``averaging_time`` short enough to be physical if the other convention is
wanted.

**Crosswind position.**  The model reports its centreline; the trial reports
whatever the sensors on that arc happened to see, and the plume does not
oblige by passing over one.  The arc maximum is used, which biases *in the
model's favour* -- if no sensor sat on the centreline the observation is low.
That direction is worth stating plainly, because it means a model that reads
high here reads higher still against a properly sampled centreline.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..run import run_steady
from .rediphem import Trial
from .statistics import Statistics, statistics
import inspect

from .trialcase import (
    SMEDIS_SOURCES,
    TrialCase,
    computed_source,
    jet_to_case,
    puff_to_case,
    to_case,
)


@dataclass
class Comparison:
    """One trial's predictions beside its measurements."""

    trial: Trial
    height: float | None  #: the single height compared, or None for all
    arcs: np.ndarray  #: downwind distance, m
    heights: np.ndarray | None = None  #: sensor elevation of each pair, m
    observed: np.ndarray = None  #: mole per cent
    predicted: np.ndarray = None  #: mole per cent
    case: TrialCase = None
    note: str = ""

    @property
    def ratio(self) -> np.ndarray:
        return self.predicted / self.observed

    def statistics(self, **kw) -> Statistics:
        return statistics(self.observed, self.predicted, **kw)

    def __str__(self) -> str:
        single = self.heights is None or len(set(self.heights)) == 1
        title = (
            f"{self.trial.series} {self.trial.name}"
            + (f" at z = {self.height:g} m" if single else " (all heights)")
        )
        cols = f"  {'arc (m)':>9}" + ("" if single else f" {'z (m)':>6}")
        cols += f" {'observed':>10} {'predicted':>10} {'ratio':>7}"
        rows = []
        for i, (x, o, p) in enumerate(zip(self.arcs, self.observed, self.predicted)):
            line = f"  {x:>9.1f}"
            if not single:
                line += f" {self.heights[i]:>6.1f}"
            rows.append(f"{line} {o:>10.3f} {p:>10.3f} {p / o:>7.2f}")
        return "\n".join([title, cols, *rows])

    def by_height(self) -> dict[float, Statistics]:
        """Statistics split by sensor elevation.

        Comparing one height and correcting the model down to it with its own
        vertical profile is circular. Splitting by height is what tests the
        profile: a model whose vertical structure is wrong shows a bias that
        marches with elevation even when the overall number looks fine.
        """
        if self.heights is None:
            return {}
        out = {}
        for z in sorted(set(self.heights)):
            mask = self.heights == z
            try:
                out[float(z)] = statistics(self.observed[mask], self.predicted[mask])
            except ValueError:
                pass
        return out


#: Which of the shared deck-building options ``puff_to_case`` accepts. An
#: instantaneous release has no jet reference height and no averaging over a
#: plume passage, so passing the full set raises rather than being ignored.
_PUFF_ARGUMENTS = set(
    inspect.signature(puff_to_case).parameters
)


#: The most recent :class:`TrialCase` that failed to run, per trial name, so
#: that the series driver can report what actually went wrong.
_LAST_PROBLEM: dict = {}


def _case_for(trial, **kw):
    """Build a deck by the route the trial's release type calls for.

    Used when reporting why a trial was skipped, so that the reason describes
    the deck that was actually run rather than a pool deck built for the
    occasion.
    """
    equivalent = SMEDIS_SOURCES.get((trial.series, trial.name))
    if equivalent is None and trial.release_type == "jet":
        try:
            equivalent = computed_source(trial)
        except Exception:
            equivalent = None
    if equivalent is not None:
        return jet_to_case(trial, equivalent)
    if trial.release_type == "puff":
        return puff_to_case(trial)
    return to_case(trial)


def compare(
    trial: Trial,
    *,
    height: float | None = None,
    backend: str | None = None,
    floor: float = 0.1,
    average_measurements: bool = False,
    surface_temperature: float | None = None,
    **case_kw,
) -> Comparison | None:
    """Run a trial and pair the result with its measurements.

    Returns ``None`` when the trial cannot be modelled -- the reason is on the
    :class:`~degali.validation.trialcase.TrialCase`, and callers that want
    to report skipped trials should build the case themselves.
    """
    equivalent = SMEDIS_SOURCES.get((trial.series, trial.name))
    if equivalent is None and trial.release_type == "jet":
        # no published source term: compute one, since the alternative is to
        # send a flashing jet through the pool source, which is meaningless
        equivalent = computed_source(trial)
    if equivalent is not None:
        tc = jet_to_case(trial, equivalent,
                         surface_temperature=surface_temperature, **case_kw)
    elif trial.release_type == "puff":
        # An instantaneous release is not a pool. Sending one through the
        # steady pool route turns a 14 m cylinder of Freon into a source of
        # 3970 kg/s over a 2050 m radius -- a number with no physical
        # meaning, which then fails a degeneracy check and is reported as
        # "no usable measurements". The route has to be chosen from the
        # release type, not defaulted to.
        accepted = {
            k: v for k, v in case_kw.items()
            if k in _PUFF_ARGUMENTS
        }
        tc = puff_to_case(
            trial, surface_temperature=surface_temperature, **accepted
        )
    else:
        tc = to_case(trial, surface_temperature=surface_temperature, **case_kw)
    if not tc.usable:
        return None

    averaging = tc.case.avtime if average_measurements else 0.0
    if height is None:
        raw = trial.arc_maxima_by_height(averaging=averaging)
    else:
        raw = {(x, height): c
               for x, c in trial.arc_maxima(height=height, averaging=averaging).items()}
    observed = {k: c for k, c in raw.items() if c > floor}
    if not observed:
        return None

    try:
        profile, source = run_steady(tc.case, backend=backend)
    except (RuntimeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        # Record why on the *trial case*, and hand it back: the caller used to
        # rebuild a pool deck to ask what went wrong, and a pool deck's verdict
        # on a puff has nothing to do with a van Ulden momentum balance that
        # did not converge. Thorney Island was being reported as "no usable
        # measurements" -- a true-sounding sentence for a false reason.
        tc.problems.append(f"model failed: {exc}")
        _LAST_PROBLEM[trial.name] = tc
        return None

    # a jet case reports distance from the release point; the equivalent
    # source sits downstream of it and OODIST already carries that offset
    x_model = profile.rows[:, 0]
    yc = profile.rows[:, 1]
    sz = profile.rows[:, 7]
    alpha1 = source.alpha + 1.0

    arcs, zs, obs, pred = [], [], [], []
    for x, z in sorted(observed):
        if x < x_model[0] or x > x_model[-1]:
            continue
        centre = float(np.interp(x, x_model, yc))
        depth = float(np.interp(x, x_model, sz))
        arcs.append(x)
        zs.append(z)
        obs.append(observed[(x, z)])
        pred.append(centre * np.exp(-((z / depth) ** alpha1)) * 100.0)

    if not arcs:
        return None
    return Comparison(
        trial=trial, height=height, arcs=np.array(arcs),
        heights=np.array(zs),
        observed=np.array(obs), predicted=np.array(pred), case=tc,
        note=f"{'block-averaged over %.0f s' % averaging if averaging else 'peak values'}",
    )


@dataclass
class SeriesResult:
    """Every trial in a series, compared."""

    series: str
    comparisons: list[Comparison] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)

    @property
    def observed(self) -> np.ndarray:
        return np.concatenate([c.observed for c in self.comparisons])

    @property
    def predicted(self) -> np.ndarray:
        return np.concatenate([c.predicted for c in self.comparisons])

    def statistics(self, **kw) -> Statistics:
        return statistics(self.observed, self.predicted, **kw)

    def by_height(self, **kw) -> dict[float, Statistics]:
        """Statistics split by sensor elevation, pooled over the series."""
        pairs: dict[float, list] = {}
        for c in self.comparisons:
            if c.heights is None:
                continue
            for z, o, p in zip(c.heights, c.observed, c.predicted):
                pairs.setdefault(float(z), []).append((o, p))
        out = {}
        for z, vals in sorted(pairs.items()):
            o, p = zip(*vals)
            try:
                out[z] = statistics(o, p, **kw)
            except ValueError:
                pass
        return out

    def by_trial(self, **kw) -> dict[str, Statistics]:
        out = {}
        for c in self.comparisons:
            try:
                out[c.trial.name] = c.statistics(**kw)
            except ValueError:
                pass
        return out


def compare_series(
    trials: list[Trial], series: str, **kw
) -> SeriesResult:
    """Compare every trial in one series.

    A field a trial omits is filled from its siblings where they agree:
    Desert Tortoise D2 records no anemometer height while D1, D3 and D4 all
    say 2 m, and dropping the trial over that would lose a quarter of the
    series to a transcription gap rather than to anything physical.
    """
    result = SeriesResult(series=series)
    members = [t for t in trials if t.series == series]
    heights = {t.wind_height for t in members}
    heights.discard(None)
    if len(heights) == 1:
        kw.setdefault("reference_height", heights.pop())
    for t in members:
        try:
            c = compare(t, **kw)
        except (RuntimeError, ValueError) as exc:
            # The model refused to run. Say so: asking `to_case` afterwards
            # rebuilds the deck through the *pool* route regardless of what
            # the trial is, and its verdict on a puff or a jet has nothing to
            # do with why the run stopped. Thorney Island failing the van
            # Ulden momentum balance was being reported as "no usable
            # measurements", which is a true-sounding sentence for a false
            # reason.
            result.skipped[t.name] = str(exc).split(";")[0].strip()
            continue
        if c is None:
            tc = _LAST_PROBLEM.pop(t.name, None) or _case_for(t, **kw)
            result.skipped[t.name] = (
                "; ".join(tc.problems) if tc.problems else "no usable measurements"
            )
        else:
            result.comparisons.append(c)
    return result
