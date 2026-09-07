"""Model evaluation statistics.

The measures used to judge dense-gas models come from Hanna, Chang and
Strimaitis (1993) and are the ones the EPA and FERC evaluations of DEGADIS
report.  All of them compare a set of predictions :math:`C_p` against
observations :math:`C_o`.

The geometric measures are the ones that matter here.  Dispersion
concentrations span orders of magnitude down a plume, so an arithmetic bias
is dominated by the near field and says almost nothing about the far field
where an exclusion zone is actually drawn.

.. math::
    MG = \\exp\\overline{(\\ln C_o - \\ln C_p)}, \\qquad
    VG = \\exp\\overline{(\\ln C_o - \\ln C_p)^2}

``MG`` is the geometric mean bias: 1 is unbiased, above 1 means the model
under-predicts, below 1 means it over-predicts.  ``VG`` is the geometric
variance, the scatter about that bias, and 1 is perfect.  ``FAC2`` is the
fraction of predictions within a factor of two of the observation, which is
the measure practitioners trust most because it is hard to game.

Note the sign convention: ``MG`` is defined observed-over-predicted, so a
model that reads high gives ``MG < 1``.  Getting this backwards is easy and
the published tables are not always explicit about it, so
:attr:`Statistics.reads` says which way round it came out in words.

Acceptance
----------
Hanna's criteria for an acceptable dense-gas model are ``0.7 < MG < 1.3``,
``VG < 1.6`` and ``FAC2 > 0.5``.  They are conventions rather than physics --
a model can fail them and still be the right tool, and pass them on a dataset
too easy to discriminate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

#: Hanna, Chang and Strimaitis (1993) acceptance criteria.
MG_RANGE = (0.7, 1.3)
VG_LIMIT = 1.6
FAC2_LIMIT = 0.5
FB_RANGE = (-0.3, 0.3)
NMSE_LIMIT = 4.0


@dataclass(frozen=True)
class Statistics:
    """Performance measures for one set of paired predictions."""

    n: int
    mg: float  #: geometric mean bias, observed over predicted
    vg: float  #: geometric variance
    fac2: float  #: fraction within a factor of two
    fb: float  #: fractional bias
    nmse: float  #: normalised mean square error
    mean_observed: float
    mean_predicted: float

    @property
    def reads(self) -> str:
        """Which way the bias runs, in words."""
        if self.mg > 1.05:
            return f"under-predicts by {self.mg:.2f}x"
        if self.mg < 0.95:
            return f"over-predicts by {1.0 / self.mg:.2f}x"
        return "essentially unbiased"

    @property
    def passes(self) -> dict[str, bool]:
        return {
            "MG": MG_RANGE[0] < self.mg < MG_RANGE[1],
            "VG": self.vg < VG_LIMIT,
            "FAC2": self.fac2 > FAC2_LIMIT,
            "FB": FB_RANGE[0] < self.fb < FB_RANGE[1],
            "NMSE": self.nmse < NMSE_LIMIT,
        }

    @property
    def acceptable(self) -> bool:
        """Whether the three measures Hanna leads with are all satisfied."""
        p = self.passes
        return p["MG"] and p["VG"] and p["FAC2"]

    def __str__(self) -> str:
        marks = self.passes
        flag = lambda k: "ok" if marks[k] else "--"
        return (
            f"n={self.n:<4d} MG={self.mg:6.3f} [{flag('MG')}] "
            f"VG={self.vg:6.3f} [{flag('VG')}] "
            f"FAC2={self.fac2:5.2f} [{flag('FAC2')}] "
            f"FB={self.fb:+6.3f} NMSE={self.nmse:6.3f}"
        )


def statistics(
    observed, predicted, *, floor: float | None = None
) -> Statistics:
    """Compute the evaluation measures for paired values.

    Parameters
    ----------
    floor
        Values at or below this are dropped.  The logarithmic measures are
        undefined at zero, and a sensor reading its own noise carries no
        information about the model, so a threshold is unavoidable; making it
        an explicit argument keeps it from being buried in the result.
    """
    obs = np.asarray(observed, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    if obs.shape != pred.shape:
        raise ValueError("observed and predicted must have the same shape")

    keep = np.isfinite(obs) & np.isfinite(pred)
    if floor is not None:
        keep &= (obs > floor) & (pred > floor)
    else:
        keep &= (obs > 0.0) & (pred > 0.0)
    obs, pred = obs[keep], pred[keep]
    if obs.size == 0:
        raise ValueError("no usable pairs after filtering")

    log_ratio = np.log(obs) - np.log(pred)
    ratio = pred / obs
    mean_o, mean_p = obs.mean(), pred.mean()

    return Statistics(
        n=int(obs.size),
        mg=float(np.exp(log_ratio.mean())),
        vg=float(np.exp((log_ratio**2).mean())),
        fac2=float(np.mean((ratio >= 0.5) & (ratio <= 2.0))),
        fb=float(2.0 * (mean_o - mean_p) / (mean_o + mean_p)),
        nmse=float(np.mean((obs - pred) ** 2) / (mean_o * mean_p)),
        mean_observed=float(mean_o),
        mean_predicted=float(mean_p),
    )


def table(results: dict[str, Statistics], *, title: str = "") -> str:
    """Format several statistics as a table."""
    lines = []
    if title:
        lines += [title, "=" * len(title)]
    lines.append(
        f"  {'set':<22} {'n':>4} {'MG':>7} {'VG':>7} {'FAC2':>6} "
        f"{'FB':>7} {'NMSE':>7}  verdict"
    )
    lines.append("  " + "-" * 78)
    for name, s in results.items():
        verdict = "acceptable" if s.acceptable else s.reads
        lines.append(
            f"  {name:<22} {s.n:>4d} {s.mg:>7.3f} {s.vg:>7.3f} "
            f"{s.fac2:>6.2f} {s.fb:>+7.3f} {s.nmse:>7.3f}  {verdict}"
        )
    lines.append("")
    lines.append(
        f"  criteria: {MG_RANGE[0]} < MG < {MG_RANGE[1]}, VG < {VG_LIMIT}, "
        f"FAC2 > {FAC2_LIMIT};  MG is observed/predicted"
    )
    return "\n".join(lines)
