"""One-call entry points for a whole simulation.

The stage modules mirror the original's six programs, which is the right shape
for validating a port but a poor shape for using one.  This module wires them
together:

===============================  ==========================================
:func:`run_steady`               a ground-level deck to a concentration profile
:func:`run_transient`            a ground-level deck to cloud snapshots
:func:`run_jet`                  a jet deck to a trajectory and touchdown
:func:`run_jet_to_ground`        a jet deck all the way to a ground profile
===============================  ==========================================

Each takes file paths and returns typed results.  The parameter files
(``.ER1``, ``.ER2``, ``.ER3``) are optional; the values EPA ships are used
when they are not supplied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .core.atmosphere import fit_alpha, psi
from .core.crfg import SourceVectors, build_source_vectors
from .core.downwind import Downwind
from .core.driver import DriverParameters, SourceRun
from .core.dose import DoseHistory, DoseRun, Receptor
from .core.jetplume import JetCoefficients, JetPlume, JetResult
from .core.numerics import gamma as _gamma
from .core.observer import release_times
from .core.steady import Profile, SteadyStateRun
from .core.thermo import (
    AdiabaticTable,
    AmbientConditions,
    GasProperties,
    LegacyBackend,
    Thermo,
    ThermoBackend,
)
from .core.timesort import Snapshot, TimeSort, sort_times
from .core.transient import TransientResult, TransientRun
from .io.bridge import bridge, read_in
from .io.inp import Case, read_inp
from .io.jetdeck import JetDeck, Touchdown, read_ino
from .io.params import (
    NumericalParameters,
    alph_settings,
    read_er1,
    read_er2,
    read_parameters,
)

#: Defaults matching the ``EXAMPLE.ER1`` / ``EXAMPLE.ER2`` files EPA ships.
_DEFAULT_ER1 = dict(
    stpin=0.01, erbnd=0.0025, wtrg=1.0, wttm=1.0, wtya=1.0, wtyc=1.0,
    wteb=1.0, wtmb=1.0, xli=0.05, xri=0.50, eps=1.0e-5, zlow=0.01,
    stpinz=-0.1, erbndz=0.005, srcoer=0.007, srcss=5.2, srccut=1.0e-5,
    ernobl=1.0005, noblpt=100, crfger=0.008, epsilon=0.59, ce=1.15,
    delrmn=0.0, szstp0=0.01, szerr=0.001, szsz0=0.01, ialpfl=1, alpco=0.2,
    iphifl=3, dellay=2.15, vua=1.3, vub=1.2, vue=20.0, vud=0.64,
    vudelta=0.20,
)
_DEFAULT_ER2 = dict(
    sy0er=0.0, erro=0.003, sz0er=1.0e-5, wtaio=1.0, wtqoo=1.0, wtszo=1.0,
    errp=0.003, smxp=80.0, wtszp=1.0, wtsyp=1.0, wtbep=1.0, wtdh=1.0,
    errg=0.003, smxg=120.0, ertdnf=5.0e-4, ertupf=5.0e-4, wtruh=1.0,
    wtdhg=1.0, stpo=0.05, stpp=0.05, odlp=0.06, odllp=80.0, stpg=60.0,
    odlg=0.045, odllg=80.0, nobs=30,
)


def _params(er1: str | Path | None):
    if er1 is None:
        raw = NumericalParameters(values=dict(_DEFAULT_ER1), labels=[])
        v = raw.values
        return DriverParameters(
            ce=v["ce"], epsilon=v["epsilon"], delrmn=v["delrmn"],
            dellay=v["dellay"], iphifl=int(v["iphifl"]), srccut=v["srccut"],
            srcss=v["srcss"], srcoer=v["srcoer"], vua=v["vua"], vub=v["vub"],
            vue=v["vue"], vud=v["vud"], stpin=v["stpin"], erbnd=v["erbnd"],
            wtrg=v["wtrg"], wttm=v["wttm"], wtyc=v["wtyc"], wtya=v["wtya"],
            wteb=v["wteb"], wtmb=v["wtmb"], ernobl=v["ernobl"],
            noblpt=int(v["noblpt"]),
        ), raw
    return read_er1(er1)


def _er2(path: str | Path | None) -> NumericalParameters:
    if path is None:
        return NumericalParameters(values=dict(_DEFAULT_ER2), labels=[])
    return read_er2(path)


@dataclass
class SourceResult:
    """What the source stage produces, and everything the next stage needs."""

    case: Case
    thermo: Thermo
    alpha: float
    gammaf: float
    ustar: float
    vectors: SourceVectors
    blanket: object  #: :class:`~degali.core.blanket.BlanketResult`
    params: DriverParameters
    raw: NumericalParameters

    def handoff(self) -> dict:
        """The scalars ``DEG2`` would have read from a ``.TR2`` file."""
        b = self.blanket
        return {
            "ess": b.ess, "outl": b.outl, "outb": b.outb, "outsz": b.outsz,
            "outcc": b.outcc, "swcl": b.swcl, "swal": b.swal, "senl": b.senl,
            "srhl": b.srhl, "yclow": self.case.yclow,
            "humid": self.case.ambient.humid, "oodist": self.case.oodist,
            "rhoe": self.thermo.table.rhoe, "rhoa": self.thermo.table.rhoa,
            "alpha": self.alpha, "gammaf": self.gammaf, "ustar": self.ustar,
            "deltay": self.case.stability.deltay,
            "betay": self.case.stability.betay,
            "aleph": b.aleph, "rmax": b.rmax, "tend": b.tend,
            "sigxco": self.case.stability.sigxco,
            "sigxp": self.case.stability.sigxp,
            "sigxmd": self.case.stability.sigxmd,
            "u0": self.case.u0, "z0": self.case.z0,
        }


def run_source(
    case: Case | str | Path,
    *,
    er1: str | Path | None = None,
    backend: ThermoBackend | str | None = None,
) -> SourceResult:
    """Run the source calculation: ``DEG1`` plus ``CRFG``."""
    if not isinstance(case, Case):
        # the deck is always *read* with the legacy correlations, because the
        # humidity reconciliation and the source enthalpy columns are part of
        # the deck's own definition; the backend choice applies to the
        # simulation that follows
        case = read_inp(case)
    params, raw = _params(er1)
    thermo = case.make_thermo(backend=backend)
    run = SourceRun(case, thermo, params)
    alpha = fit_alpha(
        case.u0, case.z0, case.zr, case.rml, ustar=run.ustar, legacy=True,
        **alph_settings(raw),
    )
    run.set_alpha(alpha)
    blanket = run.run()
    vectors = build_source_vectors(run.raw_records, crfger=raw["crfger"])
    return SourceResult(
        case=case, thermo=thermo, alpha=alpha, gammaf=run.gammaf,
        ustar=run.ustar, vectors=vectors, blanket=blanket, params=params,
        raw=raw,
    )


def _downwind(src: SourceResult, *, rebuild_table: bool) -> Downwind:
    """Build the downwind evaluator.

    ``DEG2S`` rebuilds the mixing table on the secondary-source composition
    before it starts, because a steady plume leaves the source with one fixed
    composition. ``DEG2`` does not: each observer collects its own mixture and
    ``SSSUP`` rebuilds the table per observer, so the table handed over stays
    the pure-contaminant one ``DEG1`` wrote.
    """
    h = src.handoff()
    th = src.thermo
    if rebuild_table and th.ambient.isofl != 1 and h["swcl"] > 0.0:
        th.ambient.humsrc = (
            1.0 - h["swcl"] - h["swal"] * (1.0 + th.ambient.humid)
        ) / h["swcl"]
        th.build_adiabatic_table(h["swcl"], h["swal"], h["senl"])
    return Downwind(
        src.case, th, src.params, alpha=src.alpha, gammaf=src.gammaf,
        ustar=src.ustar, deltay=h["deltay"], betay=h["betay"],
        rhoa=th.table.rhoa,
    )


def run_steady(
    case: Case | str | Path,
    *,
    er1: str | Path | None = None,
    er2: str | Path | None = None,
    backend: ThermoBackend | str | None = None,
) -> tuple[Profile, SourceResult]:
    """A steady ground-level release, from deck to concentration profile."""
    src = run_source(case, er1=er1, backend=backend)
    dw = _downwind(src, rebuild_table=True)
    profile = SteadyStateRun(
        dw, src.handoff(), _er2(er2), oodist=src.case.oodist
    ).run()
    return profile, src


@dataclass
class TransientOutput:
    """A transient run: observer profiles, snapshots, and the source behind them."""

    field: TransientResult
    snapshots: list[Snapshot]
    times: np.ndarray
    source: SourceResult
    sorter: TimeSort
    release_times: np.ndarray

    def dose(self, receptors, *, maxnt: int = 40) -> list[DoseHistory]:
        """Concentration time histories at fixed receptors (``DEG4``)."""
        h = self.source.handoff()
        dr = DoseRun(
            self.sorter, oodist=self.source.case.oodist,
            alpha1=self.source.alpha + 1.0,
        )
        return dr.run(
            receptors, self.release_times, rmax=h["rmax"], aleph=h["aleph"],
            maxnt=maxnt,
        )


def run_transient(
    case: Case | str | Path,
    *,
    er1: str | Path | None = None,
    er2: str | Path | None = None,
    times: np.ndarray | None = None,
    sigxfl: bool = True,
    backend: ThermoBackend | str | None = None,
) -> TransientOutput:
    """An unsteady ground-level release, from deck to cloud snapshots."""
    src = run_source(case, er1=er1, backend=backend)
    dw = _downwind(src, rebuild_table=False)
    h = src.handoff()
    er2p = _er2(er2)
    tr = TransientRun(dw, h, er2p, src.vectors, oodist=src.case.oodist)
    tr.hmrte = src.thermo.hmrte
    field = tr.run()
    t0s = release_times(tr.kin, src.vectors, h["tend"], tr.nobs)
    if times is None:
        times = sort_times(
            t0s, rmax=h["rmax"], aleph=h["aleph"], alpha1=src.alpha + 1.0
        )
    sorter = TimeSort(
        field, sigxco=h["sigxco"], sigxp=h["sigxp"], sigxmd=h["sigxmd"],
        sigxfl=sigxfl, alpha1=src.alpha + 1.0, gammaf=src.gammaf,
        table=src.thermo.table, gas=src.case.gas, ambient=src.thermo.ambient,
    )
    return TransientOutput(
        field=field, snapshots=sorter.run(times), times=np.asarray(times),
        source=src, sorter=sorter, release_times=t0s,
    )


def run_jet(
    deck: JetDeck | str | Path,
    *,
    backend: ThermoBackend | str | None = None,
) -> tuple[JetResult, JetDeck]:
    """A pressurised release, from ``.INO`` deck to trajectory and touchdown."""
    if not isinstance(deck, JetDeck):
        deck = read_ino(deck)
    rhoa = float(deck.den[0, 2])
    rhoe = float(deck.den[-1, 2])
    gas = GasProperties(
        mw=deck.gasmw, temp=deck.gastem, rho=rhoe, cpk=deck.gascpk,
        cpp=deck.gascpp, ulc=deck.gasulc, llc=deck.gasllc, zzc=deck.gaszzc,
    )
    ambient = AmbientConditions(
        tamb=deck.tamb, pamb=deck.pamb, humid=deck.humid, tsurf=deck.tsurf,
        isofl=deck.isofl, ihtfl=0, iwtfl=0,
    )
    if isinstance(backend, str):
        from .core.fluids import resolve
        from .core.thermo import make_backend

        if backend == "coolprop":
            gas.coolprop_name = resolve(None, deck.gasmw).fluid
            backend = make_backend("coolprop", contaminant=gas.coolprop_name)
        else:
            backend = make_backend(backend)
    th = Thermo(gas=gas, ambient=ambient, backend=backend or LegacyBackend())
    # the jet deck supplies its own density table, so no mixing line is built
    th.table = AdiabaticTable(
        yc=deck.den[:, 0], cc=deck.den[:, 1], rho=deck.den[:, 2],
        h=deck.den[:, 3], t=deck.den[:, 4], humid=deck.humid, humsrc=0.0,
        gasmw=deck.gasmw,
    )
    jp = JetPlume(
        th, coefficients=JetCoefficients(alfa1=deck.alfa1, alfa2=deck.alfa2),
        u0=deck.u0, z0=deck.z0, zr=deck.zr, rml=deck.rml, ustar=deck.ustar,
        rhoa=rhoa, rhoe=rhoe, deltay=deck.deltay, betay=deck.betay,
        deltaz=deck.deltaz, betaz=deck.betaz, gammaz=deck.gammaz,
        yclow=deck.yclow,
    )
    ua = deck.ustar / 0.35 * (
        math.log((deck.elejet + deck.zr) / deck.zr) - psi(deck.elejet, deck.rml)
    )
    y0 = jp.initial_conditions(
        erate=deck.erate, diajet=deck.diajet, elejet=deck.elejet, ua=ua
    )
    return jp.run(y0, distmx=deck.distmx), deck


def run_jet_to_ground(
    ino: str | Path,
    jet_in: str | Path,
    *,
    er1: str | Path | None = None,
    er2: str | Path | None = None,
    backend: ThermoBackend | str | None = None,
) -> tuple[Profile, JetResult, SourceResult]:
    """A pressurised release all the way to a ground-level profile.

    Runs the jet, bridges the touchdown into a ground-level case, and carries
    it downwind. Raises :class:`ValueError` if the plume never lands, which is
    itself the answer: there is no ground-level hazard.
    """
    jet, _deck = run_jet(ino, backend=backend)
    td = Touchdown(jet.distance, jet.concentration, jet.halfwidth)
    case = bridge(td, read_in(jet_in))
    profile, src = run_steady(case, er1=er1, er2=er2, backend=backend)
    return profile, jet, src


__all__ = [
    "SourceResult", "TransientOutput", "Receptor",
    "run_source", "run_steady", "run_transient", "run_jet",
    "run_jet_to_ground",
]
