"""Regression of the Python port against the original Fortran.

Every assertion here compares against full double-precision state pulled out
of a live run of DEGADIS 2.1, not against its printed output.  ``RTOL`` is set
at 1e-12: tight enough that a genuine algebraic difference fails, loose enough
that operation ordering does not.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from degali.core.atmosphere import (
    absolute_humidity,
    ambient_density,
    fit_alpha,
    friction_velocity,
    psi,
    stability_defaults,
)
from degali.core.entrainment import phi, phi_hat, surface_exchange
from degali.core.numerics import gamma
from degali.core.thermo import (
    AmbientConditions,
    GasProperties,
    LegacyBackend,
    Thermo,
)
from degali.validation.reference import Reference, TEST_CASES

RTOL = 1.0e-12
REFERENCE_ROOT = Path(__file__).resolve().parents[1] / "reference"


def _rel(got: float, ref: float) -> float:
    return abs(got - ref) / max(abs(ref), 1.0e-300)


@pytest.fixture(scope="session")
def reference() -> Reference:
    ref = Reference(REFERENCE_ROOT)
    if os.name == "nt" and not (ref.bin / "probe.exe").exists():
        pytest.skip(
            "the packaged Fortran oracle and state probes are Linux binaries; "
            "run parity tests under Linux/WSL or provide a native Windows build"
        )
    if not (ref.bin / "deg1").exists():
        ref.build()
    return ref


@pytest.fixture(scope="session")
def b9(reference, tmp_path_factory):
    """Burro 9 steady state, with both probe dumps."""
    work = tmp_path_factory.mktemp("b9")
    run = reference.run("b9", work)
    import subprocess

    subprocess.run([str(reference.bin / "probe2"), "b9"], cwd=work, capture_output=True)
    p2 = json.loads((work / "probe2.json").read_text())
    return run, p2


# ==========================================================================
# End to end: the Fortran we build must still reproduce EPA's golden output
# ==========================================================================


@pytest.mark.parametrize("case", sorted(TEST_CASES))
def test_fortran_reproduces_epa_golden(reference, tmp_path_factory, case):
    """The ported Fortran must match EPA's 2012 listings apart from timestamps.

    This guards the *oracle*.  If a portability patch ever changes a number,
    every downstream comparison in this file becomes meaningless, so it is
    checked first.
    """
    import re

    run = reference.run(case, tmp_path_factory.mktemp(case))
    stamp = re.compile(r"\d{1,2}-[A-Z]{3}-\d{4}")

    def norm(text: str) -> list[str]:
        return [
            line.rstrip()
            for line in text.replace("\r", "").replace("\x1a", "").splitlines()
        ]

    golden, produced = norm(run.golden), norm(run.listing)
    assert len(produced) >= len(golden) - 1

    mismatched = [
        (a, b)
        for a, b in zip(golden, produced)
        if a != b and not (stamp.search(a) or stamp.search(b))
    ]
    # b9t drifts in the 7th significant figure on 8 of 1086 lines; the EPA
    # readme calls that out explicitly as expected between compilers.
    assert len(mismatched) <= (8 if case == "b9t" else 0), mismatched[:5]


# ==========================================================================
# atmosphere
# ==========================================================================


def test_stability_defaults_table():
    """ATMDEF is pure data; check one row per class and the averaging-time
    floor that class F applies below 4.6 s."""
    assert stability_defaults(0.1, "D", 600.0).rml == 0.0  # neutral sentinel
    assert stability_defaults(0.1, "C", 600.0).deltay == pytest.approx(0.210)
    assert stability_defaults(0.1, "F", 1.0).deltay == stability_defaults(
        0.1, "F", 4.6
    ).deltay
    for cls, idx in zip("ABCDEF", range(1, 7)):
        assert stability_defaults(0.1, cls, 600.0) == stability_defaults(
            0.1, idx, 600.0
        )


def test_psi_neutral_is_zero():
    assert psi(10.0, 0.0) == 0.0
    assert psi(10.0, 100.0) < 0.0  # stable
    assert psi(10.0, -100.0) > 0.0  # unstable


def test_friction_velocity(b9):
    run, _ = b9
    got = friction_velocity(run["u0"], run["z0"], run["zr"], run["rml"])
    assert _rel(got, run["ustar"]) < RTOL


def test_alpha_legacy_matches_fortran(b9):
    """The RKGST+ZBRENT path must reproduce ALPH to round-off."""
    run, _ = b9
    got = fit_alpha(
        run["u0"], run["z0"], run["zr"], run["rml"], ustar=run["ustar"], legacy=True
    )
    assert _rel(got, run["alpha"]) < 1.0e-11


def test_alpha_modern_within_original_tolerance(b9):
    """The accurate path may differ, but only inside the original's own EPS."""
    run, _ = b9
    got = fit_alpha(
        run["u0"], run["z0"], run["zr"], run["rml"], ustar=run["ustar"], legacy=False
    )
    assert abs(got - run["alpha"]) < 1.0e-5


def test_gammaf(b9):
    run, _ = b9
    assert _rel(gamma(1.0 / (run["alpha"] + 1.0)), run["gammaf"]) < RTOL


def test_humidity_and_ambient_density(b9):
    run, _ = b9
    be = LegacyBackend()
    humid, relhum = absolute_humidity(
        run["tamb"], run["pamb"], be.water_vapour_pressure, relhum=13.0
    )
    assert _rel(humid, run["humid"]) < RTOL
    assert relhum == pytest.approx(13.0)
    assert _rel(ambient_density(run["tamb"], run["pamb"], run["humid"]), run["rhoa"]) < RTOL


# ==========================================================================
# thermo
# ==========================================================================


def _thermo(run) -> Thermo:
    gas = GasProperties(
        mw=run["gasmw"],
        temp=run["gastem"],
        rho=run["gasrho"],
        cpk=run["gascpk"],
        cpp=run["gascpp"],
        ulc=run["gasulc"],
        llc=run["gasllc"],
        zzc=run["gaszzc"],
    )
    amb = AmbientConditions(
        tamb=run["tamb"],
        pamb=run["pamb"],
        humid=run["humid"],
        tsurf=run["tsurf"],
        isofl=run["isofl"],
        ihtfl=run["ihtfl"],
        iwtfl=run["iwtfl"],
    )
    return Thermo(gas=gas, ambient=amb, backend=LegacyBackend())


def test_reference_enthalpies(b9):
    run, _ = b9
    th = _thermo(run)
    hmrte, harte, hwrte = th.reference_enthalpies()
    assert _rel(hmrte, run["hmrte"]) < RTOL
    assert harte == run["harte"] == 0.0
    assert _rel(hwrte, run["hwrte"]) < RTOL


def test_adiabatic_table_matches_setden(b9):
    """SETDEN's adaptive node placement must be reproduced exactly.

    Every later lookup interpolates on these nodes, so a single displaced node
    would shift the whole downwind solution.
    """
    run, _ = b9
    th = _thermo(run)
    th.reference_enthalpies()
    table = th.build_adiabatic_table(1.0, 0.0, run["hmrte"])
    ref = np.array(run["den"])
    got = table.as_array()

    assert got.shape == ref.shape, f"node count {got.shape[0]} != {ref.shape[0]}"
    for j, name in enumerate(["yc", "cc", "rho", "h", "T"]):
        dev = np.abs(got[:, j] - ref[:, j]) / np.maximum(np.abs(ref[:, j]), 1e-300)
        assert dev.max() < RTOL, f"{name}: max rel {dev.max():.3e}"


def test_setden_single_precision_grid_quirk(b9):
    """The dilution grid uses FLOAT(), i.e. single precision.

    Dropping that quirk shifts every node by about 1e-6 relative -- small, but
    far outside the tolerance the rest of the port holds to.
    """
    run, _ = b9
    th = _thermo(run)
    th.reference_enthalpies()
    faithful = th.build_adiabatic_table(1.0, 0.0, run["hmrte"]).as_array()
    exact = th.build_adiabatic_table(
        1.0, 0.0, run["hmrte"], exact_grid=True
    ).as_array()
    dev = np.abs(faithful[:, 0] - exact[:, 0]) / np.maximum(exact[:, 0], 1e-30)
    assert dev.max() > 1e-7


def test_adiabatic_lookups_are_self_consistent(b9):
    """from_concentration and from_mass_fraction must agree on the nodes."""
    run, _ = b9
    th = _thermo(run)
    th.reference_enthalpies()
    table = th.build_adiabatic_table(1.0, 0.0, run["hmrte"])
    for i in range(1, table.n - 1):
        a = table.from_concentration(float(table.cc[i]))
        assert _rel(a.rho, float(table.rho[i])) < 1e-10
        assert _rel(a.temp, float(table.t[i])) < 1e-8


# ==========================================================================
# entrainment
# ==========================================================================


def test_phif_all_prescriptions(b9):
    """All five IPHIFL branches over the sign and magnitude range of Ri."""
    _, p2 = b9
    assert len(p2["phif"]) == 180
    for ifl, ri, rit, ref in p2["phif"]:
        assert _rel(phi(ri, rit, int(ifl)), ref) < RTOL


def test_phif_rejects_bad_flag():
    with pytest.raises(ValueError):
        phi(1.0, 0.0, iphifl=9)


def test_phihat(b9):
    run, p2 = b9
    kw = dict(
        rhoa=run["rhoa"],
        alpha=run["alpha"],
        ustar=run["ustar"],
        u0=run["u0"],
        z0=run["z0"],
        gammaf=run["gammaf"],
        dellay=run["dellay"],
    )
    for rho, fetch, ref in p2["phihat"]:
        assert _rel(phi_hat(rho, fetch, **kw), ref) < RTOL


def test_phihat_neutral_cloud():
    kw = dict(
        rhoa=1.2, alpha=0.1, ustar=0.2, u0=5.0, z0=10.0, gammaf=1.06, dellay=2.15
    )
    assert phi_hat(1.2, 100.0, **kw) == 0.88
    assert phi_hat(1.0, 100.0, **kw) == 0.88


def test_gseries_clamp_costs_accuracy(b9):
    """GSERIES clamps |z| at 0.999 and stops at 7e-5 relative.

    That is a real approximation, not round-off: it moves phi_hat by up to
    ~2e-4.  Documented here so the legacy default is a deliberate choice.
    """
    run, p2 = b9
    kw = dict(
        rhoa=run["rhoa"],
        alpha=run["alpha"],
        ustar=run["ustar"],
        u0=run["u0"],
        z0=run["z0"],
        gammaf=run["gammaf"],
        dellay=run["dellay"],
    )
    dev = max(
        _rel(phi_hat(rho, fetch, legacy=False, **kw), ref)
        for rho, fetch, ref in p2["phihat"]
    )
    assert 1e-6 < dev < 1e-3


def test_surfac_all_prescriptions(b9):
    """Heat and water fluxes for IHTFL in {-1, 1, 2, 3}."""
    run, p2 = b9
    be = LegacyBackend()
    common = dict(
        tsurf=run["tsurf"],
        pamb=run["pamb"],
        zr=run["zr"],
        u0=run["u0"],
        z0=run["z0"],
        alpha=run["alpha"],
        ustar=run["ustar"],
        isofl=run["isofl"],
        iwtfl=run["iwtfl"],
        wtco=0.0,
        vapour_pressure=be.water_vapour_pressure,
    )
    assert len(p2["surfac"]) == 48
    for ihtfl, temp, height, w_ref, q_ref in p2["surfac"]:
        w, q = surface_exchange(
            temp, height, 1.2, 25.0, 1200.0, 0.005, ihtfl=int(ihtfl), htco=0.0, **common
        )
        assert _rel(w, w_ref) < RTOL
        assert _rel(q, q_ref) < RTOL


def test_surfac_no_transfer_when_surface_colder():
    be = LegacyBackend()
    w, q = surface_exchange(
        320.0, 1.0, 1.2, 25.0, 1200.0, 0.005,
        tsurf=300.0, pamb=1.0, zr=0.01, u0=5.0, z0=10.0, alpha=0.1,
        ustar=0.2, isofl=0, ihtfl=1, iwtfl=1, htco=0.0, wtco=0.0,
        vapour_pressure=be.water_vapour_pressure,
    )
    assert (w, q) == (0.0, 0.0)


def test_surfac_disabled_paths():
    be = LegacyBackend()
    kw = dict(
        tsurf=310.0, pamb=1.0, zr=0.01, u0=5.0, z0=10.0, alpha=0.1, ustar=0.2,
        htco=0.0, wtco=0.0, vapour_pressure=be.water_vapour_pressure,
    )
    assert surface_exchange(300.0, 1.0, 1.2, 25.0, 1200.0, 0.005,
                            isofl=1, ihtfl=1, iwtfl=1, **kw) == (0.0, 0.0)
    assert surface_exchange(300.0, 1.0, 1.2, 25.0, 1200.0, 0.005,
                            isofl=0, ihtfl=0, iwtfl=1, **kw) == (0.0, 0.0)
    assert surface_exchange(300.0, 0.0, 1.2, 25.0, 1200.0, 0.005,
                            isofl=0, ihtfl=1, iwtfl=1, **kw) == (0.0, 0.0)


# ==========================================================================
# numerics
# ==========================================================================


def test_rkgst_quadrature():
    from degali.core.rkgst import quadrature

    assert quadrature(lambda x: x * x, 0.0, 3.0, 0.01, 1e-8) == pytest.approx(9.0, rel=1e-6)
    assert quadrature(math.sin, 0.0, math.pi, 0.01, 1e-9) == pytest.approx(2.0, rel=1e-6)


def test_zbrent_precision_floor():
    """ZBRENT cannot converge below its hard-wired EPS = 3e-8.

    The convergence test is ``tol1 = 2*EPS*|b| + tol/2``, so asking for
    ``tol = 1e-14`` still stops at a bracket half-width near ``6e-8 |b|``.
    Callers that need more accuracy must use :func:`brentq`.
    """
    from scipy.optimize import brentq as sp

    from degali.core.numerics import RootBracketError, zbrent

    f = lambda x: x**3 - 2 * x - 5
    truth = sp(f, 1.0, 3.0)
    assert zbrent(f, 1.0, 3.0, 1e-14) == pytest.approx(truth, abs=1e-7)
    assert zbrent(f, 1.0, 3.0, 1e-14) != pytest.approx(truth, abs=1e-10)
    with pytest.raises(RootBracketError):
        zbrent(lambda x: x * x + 1, -1.0, 1.0, 1e-6)


# ==========================================================================
# input decks
# ==========================================================================


@pytest.fixture(scope="session")
def b9_case(b9):
    from degali.io.inp import read_inp

    run, _ = b9
    return read_inp(run.workdir / "b9.inp")


def test_inp_reader_scalars(b9, b9_case):
    """IO.FOR derives more than it reads: stability defaults, the humidity
    reconciliation, the CPP==0 shorthand and the source enthalpy/density
    columns are all computed on load."""
    run, _ = b9
    c = b9_case
    for name, got in [
        ("u0", c.u0), ("z0", c.z0), ("zr", c.zr), ("rml", c.rml),
        ("tamb", c.ambient.tamb), ("pamb", c.ambient.pamb),
        ("humid", c.ambient.humid), ("tsurf", c.ambient.tsurf),
        ("gasmw", c.gas.mw), ("gastem", c.gas.temp), ("gasrho", c.gas.rho),
        ("gascpk", c.gas.cpk), ("gascpp", c.gas.cpp),
        ("gasulc", c.gas.ulc), ("gasllc", c.gas.llc), ("gaszzc", c.gas.zzc),
        ("yclow", c.yclow), ("gmass", c.gmass0),
    ]:
        assert _rel(got, run[name]) < RTOL, name
    assert (c.istab, c.ambient.isofl, c.ambient.ihtfl, c.ambient.iwtfl) == (
        run["istab"], run["isofl"], run["ihtfl"], run["iwtfl"]
    )
    assert c.steady_state is True
    assert c.instantaneous is False


def test_inp_reader_source_table(b9, b9_case):
    run, _ = b9
    s = b9_case.source
    ref = np.array(run["source"])[: len(s)]
    got = np.column_stack(
        [s.time, s.rate, s.radius, s.wc, s.temp, s.fracv, s.enthalpy, s.rho]
    )
    dev = np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)
    assert dev.max() < RTOL


def test_inp_round_trip(b9_case, tmp_path):
    from degali.io.inp import read_inp, write_inp

    write_inp(b9_case, tmp_path / "rt.inp")
    back = read_inp(tmp_path / "rt.inp")
    assert _rel(back.u0, b9_case.u0) < RTOL
    assert _rel(back.yclow, b9_case.yclow) < RTOL
    assert _rel(back.gas.cpk, b9_case.gas.cpk) < 1e-13
    assert back.steady_state == b9_case.steady_state
    assert np.allclose(back.source.rate, b9_case.source.rate, rtol=1e-13)


# ==========================================================================
# secondary source blanket
# ==========================================================================


@pytest.fixture(scope="session")
def b9_probe3(b9):
    import subprocess

    run, _ = b9
    exe = REFERENCE_ROOT / "fortran" / "bin" / "probe3"
    subprocess.run([str(exe), "b9"], cwd=run.workdir, capture_output=True)
    return json.loads((run.workdir / "probe3.json").read_text())


def _blanket(run, p3, case):
    """Blanket configured from the Fortran's own constants.

    ``alpha`` and ``ustar`` are injected rather than refitted: they carry a
    3e-13 residual from ZBRENT's hard EPS floor, and several blanket
    derivatives are differences of nearly equal quantities that amplify it by
    three orders of magnitude. Injecting them tests the blanket algebra alone.
    """
    from degali.core.blanket import Blanket, BlanketParameters

    bl = Blanket(
        case,
        case.make_thermo(),
        BlanketParameters(
            ce=p3["ce"], epsilon=run["epsilon"], delrmn=p3["delrmn"],
            dellay=run["dellay"], srccut=p3["srccut"],
            vua=p3["vua"], vub=p3["vub"], vue=p3["vuc"], vud=p3["vud"],
        ),
    )
    bl.ustar = run["ustar"]
    bl.alpha = run["alpha"]
    bl.gammaf = run["gammaf"]
    bl.rhoa = bl.th.table.rhoa
    return bl


def _state(case, y1, y2):
    from degali.core.blanket import I_E, I_M, I_MA, I_MC, I_R

    pwcp = case.source.wc_at(0.0)
    pwap = (1.0 - pwcp) / (1.0 + case.ambient.humid)
    y = np.zeros(6)
    y[I_R], y[I_M] = y1, y2
    y[I_MC] = y2 * pwcp
    y[I_MA] = y2 * pwap
    y[I_E] = y2 * case.source.enthalpy_at(0.0)
    return y


def test_blanket_initial_conditions(b9, b9_case, b9_probe3):
    from degali.core.constants import PI

    run, _ = b9
    c, p3 = b9_case, b9_probe3
    pwcp = c.source.wc_at(0.0)
    r0 = c.source.radius_at(0.0)
    m0 = max(c.gmass0 / pwcp, PI * r0**2 * 1.1 * p3["srccut"] * c.source.rho_at(0.0))
    assert _rel(r0, p3["y1_0"]) < RTOL
    assert _rel(m0, p3["y2_0"]) < RTOL
    assert _rel(m0 * pwcp, p3["y3_0"]) < RTOL
    assert _rel(m0 * c.source.enthalpy_at(0.0), p3["y5_0"]) < RTOL


def test_blanket_derivatives_gravity_slumping(b9, b9_case, b9_probe3):
    """SRC1 over 30 states: 6 times x 5 mass and radius scalings."""
    from degali.core.blanket import I_R

    run, _ = b9
    bl = _blanket(run, b9_probe3, b9_case)
    assert len(b9_probe3["src1"]) == 30
    for row in b9_probe3["src1"]:
        t, y1, y2 = row[0], row[1], row[2]
        y = _state(b9_case, y1, y2)
        d, st = bl.derivatives(t, y, {"vuflag": False, "rmax": b9_probe3["rmax0"]})
        got = [
            d[0], d[1], d[2], d[3], d[4],
            st.qstar, st.sz, st.height, st.mixture.rho, st.richardson,
            st.mixture.yc, st.mixture.ya, d[I_R],
            st.mixture.wc, st.mixture.wa, st.mixture.enthalpy,
            st.mixture.temp, st.rate, y[I_R],
        ]
        for g, r in zip(got, row[3:]):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-11


def test_blanket_derivatives_momentum_balance(b9, b9_case, b9_probe3):
    """The van Ulden branch over 90 states.

    No EPA test case reaches this code: all five have GMASS0 = 0, so the
    height-to-diameter ratio is zero and VUFLAG never gets set. The probe
    forces it on, because it is the most intricate part of SRC1 -- a bracketed
    iteration on the frontal velocity with head and tail depths solved
    alongside it.
    """
    from degali.core.blanket import I_P, I_R

    run, _ = b9
    bl = _blanket(run, b9_probe3, b9_case)
    assert len(b9_probe3["src1_vu"]) == 90
    for row in b9_probe3["src1_vu"]:
        t, y1, y2, y6, hh_prev = row[:5]
        y = _state(b9_case, y1, y2)
        y[I_P] = y6
        d, st = bl.derivatives(
            t, y,
            {"vuflag": True, "rmax": b9_probe3["rmax0"], "vel": 0.5,
             "ht_last": 1.0, "hh_last": hh_prev},
        )
        got = [d[I_R], d[I_P], st.height, st.mixture.rho, st.richardson,
               d[I_R], st.ht, st.hh]
        for g, r in zip(got, row[5:]):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-11


def test_blanket_clamps_radius_to_pool(b9, b9_case, b9_probe3):
    """SRC1 mutates the state vector: when the blanket is smaller than the
    pool feeding it, Y(iR) is overwritten with RADP + 1e-10 and the growth
    rate becomes the pool's own. The mutation is visible to the integrator."""
    from degali.core.blanket import I_R

    run, _ = b9
    bl = _blanket(run, b9_probe3, b9_case)
    radp = b9_case.source.radius_at(0.0)
    y = _state(b9_case, radp * 0.5, b9_probe3["y2_0"])
    d, _ = bl.derivatives(0.0, y, {"vuflag": False, "rmax": radp})
    assert y[I_R] == pytest.approx(radp + bl.p.zero, rel=1e-15)
    assert d[I_R] == 0.0  # B9's pool radius is constant in time


# ==========================================================================
# numerical parameter files
# ==========================================================================


def test_er1_reader(b9, b9_probe3):
    """The .ER1 file is read positionally with FORMAT(10X, G10.4).

    The labels are decorative -- ESTRT1 never parses them -- so the reader
    must not depend on them either.
    """
    from degali.io.params import read_er1

    run, _ = b9
    params, raw = read_er1(run.workdir / "b9.er1")
    p3 = b9_probe3
    assert params.erbnd == 0.0025  # not the 0.005 the ALPHI bound uses
    for name, ref in [
        ("epsilon", run["epsilon"]), ("dellay", run["dellay"]),
        ("ce", p3["ce"]), ("delrmn", p3["delrmn"]), ("srccut", p3["srccut"]),
        ("vua", p3["vua"]), ("vub", p3["vub"]), ("vue", p3["vuc"]),
        ("vud", p3["vud"]),
    ]:
        assert getattr(params, name) == pytest.approx(ref, rel=1e-15), name
    assert params.iphifl == run["iphifl"]
    assert isinstance(params.noblpt, int) and params.noblpt == 100
    assert raw.labels[1].upper() == "ERBND"


# ==========================================================================
# the RKGST step sequence
# ==========================================================================


@pytest.fixture(scope="session")
def b9_probe4(b9):
    import subprocess

    run, _ = b9
    exe = REFERENCE_ROOT / "fortran" / "bin" / "probe4"
    subprocess.run([str(exe), "b9"], cwd=run.workdir, capture_output=True)
    return json.loads((run.workdir / "probe4.json").read_text())


def test_rkgst_step_sequence_matches(b9, b9_case, b9_probe3, b9_probe4):
    """RKGST driving the real SRC1: every step, not just the answer.

    The step-size controller is what decides where output is recorded, so
    reproducing the trajectory is not enough -- the *sequence* has to match.
    It does, to 1e-12, for the first 200 output points, which covers the whole
    of B9's blanket phase (it ends at t = 7.66 s, around step 170).

    Beyond that the two diverge, and that is expected rather than a defect:
    the adaptive controller compares a computed error against a fixed bound,
    and once round-off puts the two sides of that comparison on opposite sides
    of the bound, the step histories separate. The solutions stay close (the
    radius agrees to 1e-5 at step 1000) but the recorded points no longer line
    up, so a step-by-step comparison stops being meaningful.
    """
    from degali.core.blanket import Blanket, I_R
    from degali.core.constants import PI
    from degali.core.rkgst import rkgst
    from degali.io.params import read_er1

    run, _ = b9
    params, _raw = read_er1(run.workdir / "b9.er1")
    bl = Blanket(b9_case, b9_case.make_thermo(), params)
    bl.ustar, bl.alpha = run["ustar"], run["alpha"]
    bl.gammaf, bl.rhoa = run["gammaf"], bl.th.table.rhoa

    src = b9_case.source
    pwcp, rhop = src.wc_at(0.0), src.rho_at(0.0)
    y = np.zeros(6)
    y[I_R] = src.radius_at(0.0)
    y[1] = max(b9_case.gmass0 / pwcp, PI * y[I_R] ** 2 * 1.1 * params.srccut * rhop)
    y[2] = y[1] * pwcp
    y[3] = y[1] * (1.0 - pwcp) / (1.0 + b9_case.ambient.humid)
    y[4] = y[1] * src.enthalpy_at(0.0)
    ctx = {"vuflag": False, "rmax": float(y[I_R]), "vel": 0.0,
           "ht_last": 0.0, "hh_last": 0.0, "wc_last": pwcp}

    steps: list[list[float]] = []
    holder: dict = {}

    def fct(x, yy, dd, prmt):
        d, st = bl.derivatives(x, yy, ctx)
        dd[:] = d
        holder["st"] = st

    def outp(x, yy, dd, ihlf, ndim, prmt):
        st = holder["st"]
        steps.append([x, ihlf] + list(yy) + [st.height, st.mixture.rho, dd[I_R]])
        if len(steps) >= 220:
            prmt[4] = 1.0

    prmt = [0.0, 6.023e23, params.stpin, params.erbnd, params.stpmax] + [0.0] * 20
    weights = [params.wtrg, params.wttm, params.wtyc,
               params.wtya, params.wteb, params.wtmb]
    rkgst(fct, outp, prmt, y, weights, ndim=6)

    ref = np.array(b9_probe4["steps"][:200])
    got = np.array(steps[:200])
    assert got.shape == ref.shape
    assert (got[:, 1] == ref[:, 1]).all(), "IHLF sequence diverged"
    dev = np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)
    # 6e-10 by step 200: the state itself tracks to 1e-13, but dR/dt is a
    # difference of near-equal fluxes and amplifies that by ~10^3.
    assert dev.max() < 1e-9, f"max rel dev {dev.max():.3e}"
    assert np.abs(got[:, 2] - ref[:, 2]).max() / ref[:, 2].max() < 1e-12


# ==========================================================================
# the full DEG1 source calculation
# ==========================================================================


@pytest.fixture(scope="session")
def b9_source(b9, b9_case):
    from degali.core.crfg import build_source_vectors
    from degali.core.driver import SourceRun
    from degali.io.params import read_er1

    run, _ = b9
    params, raw = read_er1(run.workdir / "b9.er1")
    sr = SourceRun(b9_case, b9_case.make_thermo(), params)
    sr.ustar = run["ustar"]
    sr.blanket.ustar = run["ustar"]
    sr.set_alpha(run["alpha"])
    result = sr.run()
    vectors = build_source_vectors(sr.raw_records, crfger=raw["crfger"])
    return sr, result, vectors


def _scl_table(run) -> np.ndarray:
    text = run.golden.replace("\r", "")
    block = text.split("CALCULATED SOURCE PARAMETERS")[1].split("Source strength")[0]
    rows = [
        [float(x) for x in line.split()]
        for line in block.splitlines()
        if len(line.split()) == 9 and line.split()[0][0].isdigit()
    ]
    return np.array(rows)


def test_source_parameter_table(b9, b9_source):
    """The printed source table, after both thinning stages.

    SRC1O records 9 points; CRFG keeps 8. Comparing against the listing tests
    the integration, the output controller and the thinning together.
    """
    run, _ = b9
    _sr, _result, vectors = b9_source
    ref = _scl_table(run)
    got = vectors.listing
    assert got.shape == ref.shape, f"{got.shape} vs {ref.shape}"
    dev = np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)
    # the listing carries six significant figures
    assert dev.max() < 1e-5, f"max rel dev {dev.max():.3e}"


def _tr2_sections(run):
    lines = run.workdir.joinpath("b9.tr2").read_text(errors="replace").replace(
        "\r", ""
    ).split("\n")
    i = 4
    n = int(lines[i]); i += 1 + n
    n = int(lines[i]); i += 1 + n
    n = int(lines[i]); i += 1
    gen3 = np.array([[float(x) for x in lines[i + k].split()] for k in range(n)])
    return gen3, lines[i + n:]


def test_gen3_vectors_match_handoff(b9, b9_source):
    """The /GEN3/ vectors as written to the .TR2 handoff file.

    These six functions of time are the entire interface between the source
    model and the downwind model, so this is the check that matters most.
    """
    run, _ = b9
    _sr, _result, v = b9_source
    gen3, _tail = _tr2_sections(run)
    got = np.column_stack(
        [v.time, v.radg, v.qstr, v.srcden, v.srcwc, v.srcwa, v.srcen]
    )
    assert got.shape == gen3.shape
    dev = np.abs(got - gen3) / np.maximum(np.abs(gen3), 1e-300)
    assert dev.max() < 1e-5, f"max rel dev {dev.max():.3e}"  # 7 sig figs in .TR2


def test_secondary_source_summary(b9, b9_source):
    """RM, SZM, EMAX, RMAX, TSC1, ALEPH and the secondary source state."""
    run, _ = b9
    sr, result, _v = b9_source
    _gen3, tail = _tr2_sections(run)

    def find(width, first_value, rel=1e-2):
        for j, line in enumerate(tail):
            parts = line.split()
            if len(parts) != width:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if abs(vals[0] - first_value) <= rel * abs(first_value):
                return j, vals
        raise AssertionError(f"no {width}-field row starting near {first_value}")

    j, vals = find(5, result.rm)
    for got, ref in zip(
        (result.rm, result.szm, result.emax, result.rmax, sr.tsc1), vals
    ):
        assert _rel(got, ref) < 1e-6
    aleph, tend = [float(x) for x in tail[j + 1].split()]
    assert _rel(result.aleph, aleph) < 1e-6
    assert _rel(result.tend, tend) < 1e-6

    j, vals = find(4, result.outcc)
    for got, ref in zip(
        (result.outcc, result.outsz, result.outb, result.outl), vals
    ):
        assert _rel(got, ref) < 1e-6
    for got, ref in zip(
        (result.swcl, result.swal, result.senl, result.srhl),
        [float(x) for x in tail[j + 1].split()],
    ):
        assert _rel(got, ref) < 1e-6


def test_steady_exit_recomputes_qstar_from_primary_rate(b9_source):
    """The steady-state exit reuses the local Qstar.

    Label 122 overwrites it with PRMT(21)/(pi R^2) -- and PRMT(21) is ERTE,
    the *primary* source rate, not the take-up -- then falls through to the
    terminal write. So the last row of the table is not the integrated
    take-up. Reading PRMT(21) as the take-up instead puts that one row 3e-4
    out while every other value stays exact, which is how this was found.
    """
    sr, _result, v = b9_source
    from degali.core.constants import PI

    last = v.listing[-1]
    expected = sr.case.source.rate_at(last[0]) / PI / last[1] ** 2
    assert last[3] == pytest.approx(expected, rel=1e-12)


# ==========================================================================
# SZF: the blanket-free source layer
# ==========================================================================


def test_szf(b9, b9_case, b9_probe3):
    """SZF integrates the layer along the fetch instead of in time."""
    from degali.core.szf import sigma_z_over_source
    from degali.io.params import read_er1

    run, _ = b9
    params, raw = read_er1(run.workdir / "b9.er1")
    th = b9_case.make_thermo()
    kw = dict(
        table=th.table, alpha=run["alpha"], gammaf=run["gammaf"],
        u0=b9_case.u0, z0=b9_case.z0, ustar=run["ustar"], rhoa=th.table.rhoa,
        dellay=params.dellay, iphifl=params.iphifl,
        szstp0=raw["szstp0"], szerr=raw["szerr"],
    )
    wcp = b9_case.source.wc_at(0.0)
    assert len(b9_probe3["szf"]) == 16
    for q, length, sz, cclay, wclay, rholay in b9_probe3["szf"]:
        st = sigma_z_over_source(q, length, wcp, **kw)
        for got, ref in ((st.sz, sz), (st.cclay, cclay),
                         (st.wclay, wclay), (st.rholay, rholay)):
            assert _rel(got, ref) < RTOL


def test_adiabat_treats_wa_as_input(b9_case):
    """ADIABAT(ifl=1) reads the dry-air fraction it is passed.

    It writes ``wa`` only when ``wc`` is out of range, so a caller that leaves
    the argument uninitialised silently selects different table panels.
    SZLOCAL is such a caller: under /noauto its WALAY is static, zero, and
    never written, so SZF runs its lookups as though there were no air in the
    mixture at all. Reproducing that is what took SZF from 5e-3 to 4e-15.
    """
    th = b9_case.make_thermo()
    consistent = th.table.from_mass_fraction(0.1)
    as_fortran = th.table.from_mass_fraction(0.1, wa=0.0)
    assert _rel(as_fortran.rho, consistent.rho) > 1e-3
    assert as_fortran.wa == 0.0


# ==========================================================================
# downwind dispersion
# ==========================================================================


@pytest.fixture(scope="session")
def b9_probe5(b9):
    import subprocess

    run, _ = b9
    exe = REFERENCE_ROOT / "fortran" / "bin" / "probe5"
    subprocess.run([str(exe), "b9"], cwd=run.workdir, capture_output=True)
    return json.loads((run.workdir / "probe5.json").read_text())


@pytest.fixture(scope="session")
def b9_downwind(b9, b9_case, b9_probe5):
    """A Downwind configured from the .TR2 handoff, as DEG2S configures itself."""
    from degali.core.downwind import Downwind
    from degali.io.params import read_er1
    from degali.io.tr2 import read_tr2

    run, _ = b9
    handoff = read_tr2(run.workdir / "b9.tr2")
    params, _raw = read_er1(run.workdir / "b9.er1")
    th = b9_case.make_thermo()
    th.ambient.humsrc = (
        1.0 - handoff["swcl"] - handoff["swal"] * (1.0 + handoff["humid"])
    ) / handoff["swcl"]
    th.build_adiabatic_table(handoff["swcl"], handoff["swal"], handoff["senl"])
    dw = Downwind(
        b9_case, th, params,
        alpha=handoff["alpha"], gammaf=handoff["gammaf"], ustar=handoff["ustar"],
        deltay=handoff["deltay"], betay=handoff["betay"], rhoa=handoff["rhoa"],
    )
    return dw, handoff


def test_tr2_reader(b9, b9_downwind, b9_source):
    """The .TR2 handoff must read back what DEG1 wrote."""
    _dw, h = b9_downwind
    _sr, result, vectors = b9_source
    assert h["check4"] is True
    assert h.gen3.shape == (len(vectors), 7)
    assert _rel(h["outcc"], result.outcc) < 1e-6
    assert _rel(h["aleph"], result.aleph) < 1e-6
    assert h["gasnam"] == "LNG"


def test_pss_derivatives(b9_probe5, b9_downwind):
    """PSS over 48 states: 4 distances x 4 mass fluxes x 3 added-heat levels."""
    dw, _h = b9_downwind
    assert len(b9_probe5["pss"]) == 48
    for row in b9_probe5["pss"]:
        y = np.zeros(6)
        y[0:4] = row[1:5]
        dery = np.zeros(6)
        dw.sz_seed = b9_probe5["sz0"]
        st = dw.pss(row[0], y, dery, b9_probe5["erate"])
        got = list(dery) + [
            st.cc, st.bb, st.yc, st.rho, st.temp, st.gamma, st.rholay, st.sz
        ]
        for g, r in zip(got, row[5:]):
            # the flammable-mass derivatives amplify gamma, itself a small
            # difference of two densities near 1.07
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-6


def test_ssg_derivatives(b9_probe5, b9_downwind):
    """SSG over 48 states, with a virtual origin of 25 m."""
    dw, _h = b9_downwind
    assert len(b9_probe5["ssg"]) == 48
    for row in b9_probe5["ssg"]:
        y = np.zeros(6)
        y[0:2] = row[1:3]
        dery = np.zeros(6)
        dw.sz_seed = b9_probe5["sz0"]
        st = dw.ssg(row[0], y, dery, b9_probe5["erate"], 25.0)
        got = list(dery[:4]) + [
            st.cc, st.yc, st.rho, st.temp, st.gamma, st.rholay, st.sz
        ]
        for g, r in zip(got, row[3:]):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-5


def test_deg2s_rebuilds_the_mixing_table(b9_probe5, b9_downwind):
    """DEG2S calls SETDEN again on the secondary-source composition.

    The table DEG1 wrote to .TR2 is for the pure release material; the cloud
    leaving the source has already entrained air, so the mixing line differs.
    """
    dw, _h = b9_downwind
    ref = np.array(b9_probe5["den2"])
    got = dw.th.table.as_array()
    assert got.shape == ref.shape
    dev = np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)
    assert dev.max() < 1e-9


def test_gaminc_is_unregularised():
    """GAMINC undoes Numerical Recipes' normalisation on its last line.

    Treating it as scipy's regularised gammainc makes the flammable-mass
    derivative low by exactly Gamma(1/(1+alpha)) -- 6.5 per cent here -- and
    nothing else in PSS changes, which is what made it hard to spot.
    """
    from scipy.special import gamma as sp_gamma
    from scipy.special import gammainc as sp_gammainc

    from degali.core.numerics import gaminc

    a, x = 0.9039, 3.9758
    assert gaminc(a, x) == pytest.approx(sp_gammainc(a, x) * sp_gamma(a), rel=1e-13)
    assert gaminc(a, x) > 1.0  # the regularised form cannot exceed 1
    assert gaminc(a, 0.0) == 0.0


def test_series_overflow_guard():
    from degali.core.downwind import series

    assert series(1.0, 1.1063) > 0.0
    with pytest.raises(OverflowError):
        series(14.0, 1.1063)


# ==========================================================================
# the full steady-state downwind run
# ==========================================================================


@pytest.fixture(scope="session")
def b9_profile(b9, b9_downwind):
    from degali.core.steady import SteadyStateRun
    from degali.io.params import read_er2

    run, _ = b9
    dw, h = b9_downwind
    er2 = read_er2(run.workdir / "b9.er2")
    return SteadyStateRun(dw, h, er2, oodist=h["oodist"]).run(), er2


def test_er2_reader(b9):
    """.ER2 is read positionally too, and its fields feed three consumers."""
    from degali.io.params import read_er2

    run, _ = b9
    p = read_er2(run.workdir / "b9.er2")
    assert p["sy0er"] == 0.0
    assert (p["errp"], p["smxp"], p["stpp"]) == (0.003, 80.0, 0.05)
    assert (p["odlp"], p["odllp"]) == (0.06, 80.0)
    assert (p["errg"], p["smxg"], p["stpg"]) == (0.003, 120.0, 60.0)
    assert isinstance(p["nobs"], int) and p["nobs"] == 30


def test_deg2s_initial_conditions(b9_downwind, b9_probe5):
    """The material balance over the secondary source is closed first."""
    from degali.core.steady import SteadyStateRun
    from degali.io.params import read_er2

    dw, h = b9_downwind
    er2 = read_er2(REFERENCE_ROOT.parent / "x")  if False else None
    run = SteadyStateRun(dw, h, {"sy0er": b9_probe5["sy0er"]}, oodist=h["oodist"])
    y, _rl, _bb, sz0, cc, rholay = run._initial()
    assert _rel(sz0, b9_probe5["sz0"]) < RTOL
    assert _rel(cc, b9_probe5["cc0"]) < RTOL
    assert _rel(rholay, b9_probe5["rholay0"]) < 1e-9
    assert _rel(y[0], b9_probe5["y1_0"]) < 1e-9
    assert y[1] == b9_probe5["y2_0"]
    assert _rel(y[2], b9_probe5["y3_0"]) < RTOL


def _reference_profile(run) -> np.ndarray:
    text = run.golden.replace("\r", "")
    block = text.split("(m)                (kg/m**3)")[1].split("___")[0]
    rows = []
    for line in block.splitlines():
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            rows.append([float(x) for x in parts[:9]])
        except ValueError:
            pass
    return np.array(rows)


def test_downwind_profile(b9, b9_profile):
    """The concentration profile, end to end, against the .SR3 listing.

    The listing carries three significant figures, and the two runs place
    their output points at slightly different distances (see
    ``test_added_heat_integral_is_tolerance_limited``), so degali is
    interpolated onto the reference distances rather than compared row by row.
    """
    run, _ = b9
    profile, _er2 = b9_profile
    ref = _reference_profile(run)
    assert len(profile.rows) == len(ref)

    x = profile.rows[:, 0]
    # 3 % on concentration: the listing rounds to three figures and the
    # added-heat drift shifts the profile slightly in the mid field
    for col, name, tol in ((1, "yc", 0.03), (2, "cc", 0.03), (7, "sz", 0.02)):
        got = np.interp(ref[:, 0], x, profile.rows[:, col])
        dev = np.abs(got - ref[:, col]) / ref[:, col]
        assert dev.max() < tol, f"{name}: max rel dev {dev.max():.3e}"


def test_flammable_mass(b9_profile):
    """Mass above the LFL and between the two levels of concern.

    These are the numbers a dispersion study actually reports, and they
    integrate the whole model: source, both downwind stages and the
    concentration profile shape.
    """
    profile, _er2 = b9_profile
    assert profile.mass_above_lfl == pytest.approx(6740.2, rel=0.005)
    assert profile.mass_between == pytest.approx(3521.4, rel=0.005)


def test_lfl_distance(b9, b9_profile):
    """Distance to 5 mol %, the quantity 49 CFR 193.2059 turns on."""
    run, _ = b9
    profile, _er2 = b9_profile
    ref = _reference_profile(run)
    d = profile.distance_to(0.05)
    i = int(np.argmax(ref[:, 1] <= 0.05))
    assert ref[i - 1, 0] < d < ref[i, 0]


def test_added_heat_integral_is_tolerance_limited(b9, b9_downwind, b9_probe5):
    """Known limit: the added-heat state diverges from the Fortran by ~2 %.

    ``ADDHEAT`` resolves the heated temperature with ``ZBRENT`` at a 1e-3 K
    tolerance, and then forms the heat capacity as ``dh/(temp - amt)``. Early
    in the dense phase ``temp - amt`` is only a few hundredths of a kelvin, so
    that quotient inherits a percent-level uncertainty, which propagates into
    the ground heat flux and hence back into ``dh``.

    Everything else in the step tracks exactly -- distance, the bisection
    count, and the other three states all agree at the first output point --
    which is what identifies the added-heat integral as the sole source. The
    consequence is that output points land at different distances, so profile
    comparisons interpolate rather than align rows.
    """
    from degali.core.downwind import Downwind
    from degali.core.rkgst import rkgst
    from degali.core.steady import SteadyStateRun
    from degali.io.params import read_er2

    run, _ = b9
    dw, h = b9_downwind
    er2 = read_er2(run.workdir / "b9.er2")
    sr = SteadyStateRun(dw, h, er2, oodist=h["oodist"])
    y, rl, _bb, sz0, _cc, _rholay = sr._initial()
    dw.sz_seed = sz0

    steps: list[list[float]] = []
    holder: dict = {}

    def fct(x, yy, dd, prmt):
        holder["st"] = dw.pss(x, yy, dd, sr.erate)

    def outp(x, yy, dd, ihlf, ndim, prmt):
        dw.sz_seed = holder["st"].sz
        steps.append([x, ihlf] + list(yy))
        if len(steps) >= 2:
            prmt[4] = 1.0

    prmt = [rl / 2.0, 6.023e13, er2["stpp"], er2["errp"], er2["smxp"]] + [0.0] * 20
    rkgst(
        fct, outp, prmt, y,
        [er2["wtszp"], er2["wtsyp"], er2["wtbep"], er2["wtdh"], 1.0, 1.0], ndim=6,
    )

    ref = b9_probe5["dense_steps"][:2]
    assert len(steps) == 2
    for i in range(2):
        assert steps[i][0] == pytest.approx(ref[i][0], rel=1e-12)  # distance
        assert steps[i][1] == ref[i][1]  # bisection count
        for j in (2, 3, 4):  # rho*u*H, sy^2, beff
            assert steps[i][j] == pytest.approx(ref[i][j], rel=1e-6)
    # the added-heat state is the one that drifts
    assert steps[1][5] == pytest.approx(ref[1][5], rel=0.02)
    assert steps[1][5] != pytest.approx(ref[1][5], rel=1e-4)


# ==========================================================================
# transient observers
# ==========================================================================


@pytest.fixture(scope="session")
def b9t(reference, tmp_path_factory):
    """Burro 9 as a transient release, with the observer probe."""
    import subprocess

    work = tmp_path_factory.mktemp("b9t")
    run = reference.run("b9t", work)
    subprocess.run(
        [str(reference.bin / "probe6"), "b9t"], cwd=work, capture_output=True
    )
    return run, json.loads((work / "probe6.json").read_text())


@pytest.fixture(scope="session")
def b9t_observer(b9t):
    from degali.core.crfg import SourceVectors
    from degali.core.observer import Observer, ObserverKinematics
    from degali.io.inp import read_inp
    from degali.io.params import read_er1
    from degali.io.tr2 import read_tr2

    run, p6 = b9t
    h = read_tr2(run.workdir / "b9t.tr2")
    params, _raw = read_er1(run.workdir / "b9t.er1")
    case = read_inp(run.workdir / "b9t.inp")
    th = case.make_thermo()
    th.build_adiabatic_table(1.0, 0.0, p6["hmrte"])
    vectors = SourceVectors.from_handoff(h.gen3)
    kin = ObserverKinematics(alpha=h["alpha"], aleph=h["aleph"], rmax=h["rmax"])
    ob = Observer(
        th, vectors, kin, params=params, gammaf=h["gammaf"], ustar=h["ustar"],
        rhoa=h["rhoa"], u0=h["u0"], z0=h["z0"], humid=h["humid"],
    )
    return ob, kin, vectors, h


def test_transient_case_is_not_steady_state(b9t):
    """B9T is the same spill as B9 but declared transient."""
    from degali.io.inp import read_inp

    run, _p6 = b9t
    case = read_inp(run.workdir / "b9t.inp")
    assert case.steady_state is False


def test_observer_release_window(b9t, b9t_observer):
    """SSSUP spreads NOBS observers across the window in which one can be
    released and still cross the source.

    The earliest release is not simply the one that arrives when the source
    appears: for a case that grows a blanket, a later, larger radius can imply
    an even earlier release, so every recorded radius has to be checked.
    """
    from degali.core.observer import release_times

    _run, p6 = b9t
    _ob, kin, vectors, h = b9t_observer
    t0s = release_times(kin, vectors, h["tend"], p6["nobs"])
    assert len(t0s) == p6["nobs"] == 30
    assert t0s[0] == pytest.approx(p6["t01"], rel=1e-15)
    assert (t0s[1] - t0s[0]) == pytest.approx(p6["dtob"], rel=1e-12)
    assert vectors.radg_at(h["tend"]) == pytest.approx(p6["xend"], rel=1e-15)


def test_observer_kinematics_and_crossings(b9t, b9t_observer):
    """All 30 observers: release time, both edge crossings, positions, speeds.

    An observer that cannot reach an edge is not given a root to solve -- it
    is handed the endpoint instead. Calling TUPF unconditionally makes the
    Fortran abort in LIMIT, which is how the guards were found.
    """
    from degali.core.observer import (
        crossing_downwind,
        crossing_upwind,
        release_times,
    )

    _run, p6 = b9t
    _ob, kin, vectors, h = b9t_observer
    t0s = release_times(kin, vectors, h["tend"], p6["nobs"])
    xend = vectors.radg_at(h["tend"])
    r0 = vectors.radg_at(0.0)

    assert len(p6["obs"]) == 30
    for t0, row in zip(t0s, p6["obs"]):
        pdn = not (xend > kin.position(h["tend"], t0))
        pup = not (t0 <= 0.0 and kin.position(0.0, t0) > -r0)
        tup = (
            crossing_upwind(t0, kin, vectors.radg_at, tol=p6["ertupf"])
            if pup else 0.0
        )
        tdown = (
            crossing_downwind(t0, kin, vectors.radg_at, tol=p6["ertdnf"])
            if pdn else h["tend"]
        )
        got = [
            t0, tup, tdown,
            kin.position(tup, t0), kin.position(tdown, t0),
            kin.velocity(tup, t0), kin.velocity(tdown, t0),
        ]
        for g, r in zip(got, row):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-9


def test_ob_derivatives(b9t, b9t_observer):
    """OB over 60 states: 5 times x 3 half-widths x 4 collected masses."""
    import numpy as np

    _run, p6 = b9t
    ob, kin, vectors, _h = b9t_observer
    assert len(p6["ob"]) == 60
    for row in p6["ob"]:
        t = row[0]
        y = np.array(row[1:6])
        dery = np.zeros(5)
        st = ob.derivatives(t, y, dery, p6["ob_t0"], p6["ob_rlen"])
        outside = abs(kin.position(t, p6["ob_t0"])) - vectors.radg_at(t)
        if st.outside:
            got = list(dery) + [0.0] * 5 + [outside]
            ref = list(row[6:11]) + [0.0] * 5 + [row[16]]
        else:
            got = list(dery) + [
                st.cclay, st.wclay, st.walay, st.enthlay, st.rholay, outside
            ]
            ref = row[6:]
        for g, r in zip(got, ref):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-12


def test_obout_keeps_the_last_valid_layer(b9t_observer):
    """An observer outside the source contributes nothing and must not
    overwrite the layer state, because that state seeds the downwind run."""
    import numpy as np

    ob, _kin, _vectors, _h = b9t_observer
    from degali.core.observer import ObserverState

    good = ObserverState(1.0, 0.5, 0.4, -1.0, 1.2, outside=False)
    ob.record(good)
    ob.record(ObserverState(0.0, 0.0, 0.0, 0.0, 0.0, outside=True))
    assert ob.last is good


@pytest.fixture(scope="session")
def b9t_transient(b9t):
    from degali.core.crfg import SourceVectors
    from degali.core.downwind import Downwind
    from degali.core.transient import TransientRun
    from degali.io.inp import read_inp
    from degali.io.params import read_er1, read_er2
    from degali.io.tr2 import read_tr2

    run, p6 = b9t
    h = read_tr2(run.workdir / "b9t.tr2")
    params, _raw = read_er1(run.workdir / "b9t.er1")
    er2 = read_er2(run.workdir / "b9t.er2")
    case = read_inp(run.workdir / "b9t.inp")
    th = case.make_thermo()
    th.build_adiabatic_table(1.0, 0.0, p6["hmrte"])
    vectors = SourceVectors.from_handoff(h.gen3)
    dw = Downwind(
        case, th, params, alpha=h["alpha"], gammaf=h["gammaf"], ustar=h["ustar"],
        deltay=h["deltay"], betay=h["betay"], rhoa=h["rhoa"],
    )
    tr = TransientRun(dw, h, er2, vectors, oodist=h["oodist"])
    tr.hmrte = p6["hmrte"]
    return tr, h


def test_plume_seeding_per_observer(b9t, b9t_transient):
    """Each observer's OB integration and the plume it seeds.

    Fourteen quantities per observer: the five integrated states, the layer
    composition OBOUT preserved, and the derived source strength, half-width,
    concentration and sigma_z0.
    """
    from degali.core.observer import release_times

    _run, p6 = b9t
    tr, h = b9t_transient
    ref = np.array(p6["seed"])
    assert len(ref) == tr.nobs == 30

    t0s = release_times(tr.kin, tr.v, h["tend"], tr.nobs)
    for i, t0 in enumerate(t0s):
        tup, tdown, _pup, _pdn = tr._crossings(t0)
        xup = tr.kin.position(tup, t0)
        xdown = tr.kin.position(tdown, t0)
        layer, y, rlen = tr._collect(t0, tup, tdown, xup, xdown)
        erate, b, sz0, _sy0er, cc, _rholay = tr._seed_plume(layer, y, rlen, tdown)
        got = [
            y[0], y[1], y[2], y[3], y[4],
            layer.cclay, layer.wclay, layer.walay, layer.enthlay, layer.rholay,
            cc, y[2] / (b * rlen), erate, sz0,
        ]
        for g, r in zip(got, ref[i][1:]):
            assert (_rel(g, r) if r != 0 else abs(g)) < 1e-8


def test_rhoe_is_a_fixed_scalar_not_the_table(b9t_transient):
    """RHOE lives in /parm/ and is read once from the handoff.

    SSSUP rebuilds the mixing table for every observer, so reading the clamp
    from the table's last row instead would let it drift: the second observer
    would be clamped against the first observer's extrapolated centreline
    rather than the release density. That mis-set sigma_z0 by 13 %.
    """
    tr, h = b9t_transient
    assert tr.rhoe == h["rhoe"]
    assert tr.rhoe == pytest.approx(1.68448, rel=1e-5)


def test_transient_field(b9t_transient):
    """Every observer produces a profile, and the field is monotone downwind."""
    tr, _h = b9t_transient
    result = tr.run()
    assert len(result.observers) == 30
    assert all(len(o.rows) > 0 for o in result.observers)
    assert len(result.rows) > 1000

    # each row carries a time as well as a distance -- that is what makes the
    # set of steady profiles into a transient field
    for o in result.observers:
        assert np.all(np.diff(o.rows[:, 0]) >= 0)  # distance increases
        assert np.all(np.diff(o.rows[:, 1]) >= 0)  # so does arrival time
        assert o.rows[0, 1] >= o.tup - 1e-6

    peak_near, _t = result.max_concentration_at(100.0)
    peak_far, _t = result.max_concentration_at(800.0)
    assert peak_near > peak_far > 0.0


# ==========================================================================
# DEG3: cloud snapshots at fixed times
# ==========================================================================


@pytest.fixture(scope="session")
def b9t_snapshots(b9t, b9t_transient):
    from degali.core.observer import release_times
    from degali.core.timesort import TimeSort, sort_times

    _run, p6 = b9t
    tr, h = b9t_transient
    result = tr.run()
    t0s = release_times(tr.kin, tr.v, h["tend"], tr.nobs)
    times = sort_times(
        t0s, rmax=h["rmax"], aleph=h["aleph"], alpha1=tr.dw.alpha1
    )
    sorter = TimeSort(
        result, sigxco=h["sigxco"], sigxp=h["sigxp"], sigxmd=h["sigxmd"],
        sigxfl=True, alpha1=tr.dw.alpha1, gammaf=h["gammaf"],
        table=tr.dw.th.table, gas=tr.dw.case.gas, ambient=tr.dw.th.ambient,
    )
    return sorter.run(times), times


def _reference_snapshots(run) -> dict:
    import re

    text = run.golden.replace("\r", "")
    blocks = re.split(
        r"Time after beginning of spill\s+([0-9.E+]+)\s+sec", text
    )
    out = {}
    for i in range(1, len(blocks), 2):
        t = float(blocks[i])
        body = blocks[i + 1].split("For the ULC")[0]
        rows = []
        for line in body.splitlines():
            parts = line.split()
            if len(parts) >= 9 and parts[0][0].isdigit():
                try:
                    rows.append([float(x) for x in parts[:9]])
                except ValueError:
                    pass
        m = re.search(r"above the LLC is:\s*([0-9.E+-]+)", blocks[i + 1])
        out[t] = (np.array(rows), float(m.group(1)) if m else 0.0)
    return out


def test_sort_times_match(b9t, b9t_snapshots):
    """GETTS picks the snapshot instants without any user input.

    The window runs from the first observer reaching 2*Rmax to the last
    reaching 6*Rmax, with the interval rounded to whole seconds and floored at
    five. For B9T that gives 11 s spacing starting at 2 s; the reference drops
    the first because fewer than three observers had arrived.
    """
    run, _p6 = b9t
    _snaps, times = b9t_snapshots
    ref = _reference_snapshots(run)
    assert list(times[1:]) == sorted(ref)
    assert times[1] - times[0] == pytest.approx(11.0)


def test_snapshot_structure(b9t, b9t_snapshots):
    """Every snapshot has exactly the observers the Fortran found.

    The row count is a strict check: it says the same observers were in range
    at the same instants, which depends on the release times, both edge
    crossings, the whole downwind integration and the arrival-time mapping all
    being right.
    """
    run, _p6 = b9t
    snaps, _times = b9t_snapshots
    ref = _reference_snapshots(run)
    assert len(snaps) == len(ref) == 19
    for s in snaps:
        rows, _mass = ref[s.time]
        assert len(s) == len(rows), f"t={s.time}: {len(s)} vs {len(rows)}"


def test_snapshot_extent(b9t, b9t_snapshots):
    """How far the cloud has travelled at each instant."""
    run, _p6 = b9t
    snaps, _times = b9t_snapshots
    ref = _reference_snapshots(run)
    for s in snaps:
        rows, _mass = ref[s.time]
        # the listing rounds distance to three significant figures
        got = s.column("dist")[-1]
        assert got == pytest.approx(rows[-1, 0], rel=3e-3)


def test_snapshot_peak_concentration(b9t, b9t_snapshots):
    """Peak mole fraction in each snapshot.

    Tight early and looser late: by 200 s the cloud has been over warm ground
    for minutes, so the added-heat drift documented in
    ``test_added_heat_integral_is_tolerance_limited`` has had time to
    accumulate. It reaches about 20 % at the last snapshot, where the
    concentration is 2 mol % -- an order of magnitude below the LFL.
    """
    run, _p6 = b9t
    snaps, _times = b9t_snapshots
    ref = _reference_snapshots(run)
    for s in snaps:
        rows, _mass = ref[s.time]
        got = s.column("yc").max()
        tol = 0.05 if s.time <= 100.0 else 0.25
        assert got == pytest.approx(rows[:, 1].max(), rel=tol), f"t={s.time}"


def test_flammable_mass_history(b9t, b9t_snapshots):
    """Mass above the lower level of concern, snapshot by snapshot.

    This peaks and then decays as the cloud dilutes, and the peak is what a
    vapour-cloud-explosion assessment turns on. It is also the single most
    integrated number the model produces: source, both downwind stages, the
    observer machinery, the time sort and the along-wind correction all feed
    into it.
    """
    run, _p6 = b9t
    snaps, _times = b9t_snapshots
    ref = _reference_snapshots(run)

    got = np.array([s.mass_above_llc for s in snaps])
    want = np.array([ref[s.time][1] for s in snaps])
    assert got.argmax() == want.argmax()  # same snapshot holds the peak
    assert got.max() == pytest.approx(want.max(), rel=0.05)

    significant = want > 0.1 * want.max()
    dev = np.abs(got[significant] - want[significant]) / want[significant]
    assert dev.max() < 0.08
    # tighter where the mass matters: within a factor of two of the peak
    peak = want > 0.5 * want.max()
    assert (np.abs(got[peak] - want[peak]) / want[peak]).max() < 0.035


def test_along_wind_correction_can_be_disabled(b9t_transient, b9t_snapshots):
    """SIGXFL off leaves the centreline concentration untouched.

    The correction convolves each snapshot with a Gaussian of scale
    sigma_x = a x^p, because a set of quasi-steady observer plumes carries no
    along-wind spreading of its own. The scale uses the distance travelled
    *from the source* -- SORTS sets TABLE(21) once and never updates it -- not
    the spacing between output points. Reading it as the spacing makes the
    correction never fire, since consecutive points are far closer together
    than the 130 m cut-off.
    """
    from degali.core.timesort import TimeSort

    tr, h = b9t_transient
    snaps, times = b9t_snapshots
    result = tr.run()
    plain = TimeSort(
        result, sigxco=h["sigxco"], sigxp=h["sigxp"], sigxmd=h["sigxmd"],
        sigxfl=False, alpha1=tr.dw.alpha1, gammaf=h["gammaf"],
        table=tr.dw.th.table, gas=tr.dw.case.gas, ambient=tr.dw.th.ambient,
    ).run(times)
    assert len(plain) == len(snaps)
    for s in plain:
        assert np.allclose(s.column("ccstr"), s.column("cc"))
    # and it does something once the plume has travelled past SIGXMD
    late = snaps[-1]
    assert not np.allclose(late.column("ccstr"), late.column("cc"), rtol=1e-6)


# ==========================================================================
# DEG4: concentration time history at a receptor
# ==========================================================================


def test_receptor_times(b9t, b9t_transient):
    """GETTD brackets the window in which the cloud passes one receptor.

    Unlike GETTS, which spans the whole cloud lifetime, this runs from the
    first observer's arrival at that distance to the last one's.
    """
    from degali.core.dose import receptor_times
    from degali.core.observer import release_times

    _run, _p6 = b9t
    tr, h = b9t_transient
    t0s = release_times(tr.kin, tr.v, h["tend"], tr.nobs)
    times = receptor_times(
        t0s, 400.0, oodist=h["oodist"], rmax=h["rmax"], aleph=h["aleph"],
        alpha1=tr.dw.alpha1,
    )
    assert len(times) >= 4
    assert np.all(np.diff(times) > 0)
    # the window must bracket the arrivals it was built from
    def ts(t0):
        return t0 + (400.0 + h["rmax"]) ** (1.0 / tr.dw.alpha1) / h["aleph"]

    # GETTD rounds the start down to a whole second and caps the count at
    # MAXNT, so the window can stop short of the last observer's arrival
    assert times[0] <= ts(t0s[0]) + 1.0
    assert times[-1] > times[0]
    assert len(times) <= 40


def test_dose_history(b9t_transient, b9t_snapshots):
    """Concentration against time at three downwind receptors.

    The peak at each distance must agree with what the snapshots give, since
    the two views are assembled from the same observer profiles by different
    routes -- DEG3 interpolates in time at fixed instants, DEG4 interpolates
    in distance at fixed receptors.
    """
    from degali.core.dose import DoseRun, Receptor
    from degali.core.observer import release_times
    from degali.core.timesort import TimeSort

    tr, h = b9t_transient
    result = tr.run()
    t0s = release_times(tr.kin, tr.v, h["tend"], tr.nobs)
    sorter = TimeSort(
        result, sigxco=h["sigxco"], sigxp=h["sigxp"], sigxmd=h["sigxmd"],
        sigxfl=True, alpha1=tr.dw.alpha1, gammaf=h["gammaf"],
        table=tr.dw.th.table, gas=tr.dw.case.gas, ambient=tr.dw.th.ambient,
    )
    dr = DoseRun(sorter, oodist=h["oodist"], alpha1=tr.dw.alpha1)
    receptors = [Receptor(x=200.0), Receptor(x=400.0, offsets=[(0.0, 1.0), (40.0, 1.0)])]
    histories = dr.run(receptors, t0s, rmax=h["rmax"], aleph=h["aleph"])

    assert len(histories) == 2
    for hist in histories:
        assert len(hist.rows) > 5
        assert np.all(np.diff(hist.column("time")) > 0)
        peak, t_peak = hist.peak
        assert 0.0 < peak < 1.0
        assert hist.dose() > 0.0
        # nearer receptors see more
    assert histories[0].peak[0] > histories[1].peak[0]

    # off-axis is weaker than the centreline, and further off is weaker still
    off = histories[1].offaxis
    assert off.shape[1] == 2
    assert off[:, 0].max() > off[:, 1].max()
    assert off[:, 0].max() <= histories[1].peak[0] * 1.001


def test_dose_matches_snapshot_peak(b9t_transient, b9t_snapshots):
    """DEG3 and DEG4 must agree where they overlap."""
    from degali.core.dose import DoseRun, Receptor
    from degali.core.observer import release_times
    from degali.core.timesort import TimeSort

    tr, h = b9t_transient
    result = tr.run()
    t0s = release_times(tr.kin, tr.v, h["tend"], tr.nobs)
    sorter = TimeSort(
        result, sigxco=h["sigxco"], sigxp=h["sigxp"], sigxmd=h["sigxmd"],
        sigxfl=True, alpha1=tr.dw.alpha1, gammaf=h["gammaf"],
        table=tr.dw.th.table, gas=tr.dw.case.gas, ambient=tr.dw.th.ambient,
    )
    dr = DoseRun(sorter, oodist=h["oodist"], alpha1=tr.dw.alpha1)
    (hist,) = dr.run([Receptor(x=400.0)], t0s, rmax=h["rmax"], aleph=h["aleph"])
    from_dose, _t = hist.peak
    from_field, _t2 = result.max_concentration_at(400.0)
    assert from_dose == pytest.approx(from_field, rel=0.02)


# ==========================================================================
# the jet/plume model
# ==========================================================================


@pytest.fixture(scope="session", params=["ex1", "ex2", "ex3"])
def jet_case(reference, tmp_path_factory, request):
    """A jet case, run in Fortran, with the deck read back in Python."""
    from degali.io.jetdeck import read_ind, read_ino

    case = request.param
    run = reference.run(case, tmp_path_factory.mktemp(case))
    deck = read_ino(run.workdir / f"{case}.ino")
    vel = np.array([
        [float(x) for x in line.split()]
        for line in (run.workdir / f"{case}.vel").read_text(
            errors="replace"
        ).replace("\r", "").split("\n")
        if len(line.split()) == 8
    ])
    ind_path = run.workdir / f"{case}.ind"
    ind = read_ind(ind_path) if ind_path.exists() else None
    return case, deck, vel, ind


def _build_jet(deck):
    import math

    from degali.core.atmosphere import psi
    from degali.core.jetplume import JetCoefficients, JetPlume
    from degali.core.thermo import (
        AdiabaticTable,
        AmbientConditions,
        GasProperties,
        LegacyBackend,
        Thermo,
    )

    rhoa = float(deck.den[0, 2])
    rhoe = float(deck.den[-1, 2])
    gas = GasProperties(
        mw=deck.gasmw, temp=deck.gastem, rho=rhoe, cpk=deck.gascpk,
        cpp=deck.gascpp, ulc=deck.gasulc, llc=deck.gasllc, zzc=deck.gaszzc,
    )
    amb = AmbientConditions(
        tamb=deck.tamb, pamb=deck.pamb, humid=deck.humid, tsurf=deck.tsurf,
        isofl=deck.isofl, ihtfl=0, iwtfl=0,
    )
    th = Thermo(gas=gas, ambient=amb, backend=LegacyBackend())
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
    return jp, y0


def test_jet_initial_conditions(jet_case):
    """SETJET skips the zone of flow development.

    A jet leaves the orifice with a top-hat velocity profile the similarity
    model cannot start from, so SETJET uses the Kamotani and Greber
    wind-tunnel correlations to place the plume where the profile has become
    Gaussian. The first row of the .VEL file is that starting state.
    """
    _case, deck, vel, _ind = jet_case
    _jp, y0 = _build_jet(deck)
    first = vel[0]
    assert y0[4] == pytest.approx(first[0], rel=1e-5, abs=1e-15)  # x
    assert y0[5] == pytest.approx(first[1], rel=1e-5)  # z
    assert y0[0] == pytest.approx(first[2], rel=1e-5)  # cc
    assert y0[3] == pytest.approx(first[5], rel=1e-5)  # uc
    assert y0[2] == pytest.approx(first[7], rel=1e-5)  # theta


def test_jet_trajectory(jet_case):
    """The whole trajectory, against the full-precision .VEL dump.

    Four coupled balances -- contaminant mass, total mass, and momentum in x
    and z -- solved as a 4x4 linear system at every step, with the elliptical
    cross-section split by a root find and the wind averaged over it.
    """
    case, deck, vel, _ind = jet_case
    jp, y0 = _build_jet(deck)
    result = jp.run(y0, distmx=deck.distmx)
    got = result.rows
    assert len(got) > 20

    inside = (vel[:, 0] >= got[0, 0]) & (vel[:, 0] <= got[-1, 0])
    ref = vel[inside]
    assert len(ref) > 10
    # ex2's first few metres are extremely stiff -- the excess velocity falls
    # from 430 to 15 m/s in ten metres -- so the two step sequences separate
    # there and only reconverge downstream. ex3 runs for three kilometres and
    # agrees to 4e-4 except in the last few metres before touchdown.
    tolerances = {"ex1": 1e-4, "ex2": 3e-2, "ex3": 2e-3}
    tol = tolerances[case]
    for col, name in ((1, "z"), (2, "cc"), (3, "sy"), (4, "sz")):
        interp = np.interp(ref[:, 0], got[:, 0], got[:, col])
        # the final row of a landing case sits at z = 0 exactly, where a
        # relative comparison is meaningless; drop it and check the touchdown
        # distance separately in test_jet_touchdown
        # elevation goes to zero at touchdown, so the last few metres of a
        # landing plume have no useful relative scale; the landing point
        # itself is checked in test_jet_touchdown
        scale = np.abs(ref[:, col])
        if name == "z":
            scale = np.maximum(scale, 1.0)
        dev = np.abs(interp - ref[:, col]) / scale
        assert dev.max() < tol, f"{case} {name}: {dev.max():.3e}"


def test_jet_touchdown(jet_case):
    """Where the plume lands -- the handoff to the ground-level model.

    EX1 never lands, and JETPLU signals that by writing a distance of zero so
    DEGADIS does not continue; EX2 and EX3 do, and the three numbers written
    to the .IND file are what DEGBRIDG turns into a ground-level deck.
    """
    case, deck, _vel, ind = jet_case
    jp, y0 = _build_jet(deck)
    result = jp.run(y0, distmx=deck.distmx)

    if case == "ex1":
        assert not result.touchdown
        assert result.distance == 0.0
        return

    assert result.touchdown
    assert ind is not None and ind.lands
    assert result.distance == pytest.approx(ind.distance, rel=1e-4)
    assert result.concentration == pytest.approx(ind.concentration, rel=1e-4)
    assert result.halfwidth == pytest.approx(ind.halfwidth, rel=1e-4)


def test_ellipse_against_closed_forms():
    """ELLIPS reduces to the circle and the half-ellipse correctly."""
    from degali.core.constants import PI
    from degali.core.jetplume import ellipse

    perim, area = ellipse(2.0, 2.0, 2.0)
    assert area == pytest.approx(PI * 4.0)
    assert perim == pytest.approx(2.0 * PI * 2.0)

    perim, area = ellipse(3.0, 1.0, 0.0)
    assert area == pytest.approx(PI * 3.0 / 2.0)
    full, _ = ellipse(3.0, 1.0, 1.0)
    assert perim == pytest.approx(full / 2.0)


def test_ino_reader(jet_case):
    """The .INO deck is positional, like every other DEGADIS file."""
    case, deck, _vel, _ind = jet_case
    assert deck.erate > 0.0
    assert deck.diajet > 0.0
    assert 1 <= deck.istab <= 6
    assert deck.den.shape[1] == 5
    # a continuous release integrates to the level of concern, a finite one
    # to a quarter of it
    expected = deck.gasllc / 4.0 if deck.tend > 0.0 else deck.gasllc
    assert deck.yclow == expected


# ==========================================================================
# DEGBRIDG: jet touchdown to ground-level case
# ==========================================================================


@pytest.fixture(scope="session", params=["ex2", "ex3"])
def bridged_case(reference, tmp_path_factory, request):
    """A landing jet case, bridged in Python, beside the Fortran's own deck."""
    from degali.io.bridge import bridge, read_in
    from degali.io.inp import read_inp
    from degali.io.jetdeck import read_ind

    case = request.param
    run = reference.run(case, tmp_path_factory.mktemp(case))
    deck = read_in(run.workdir / f"{case}.in")
    td = read_ind(run.workdir / f"{case}.ind")
    return case, bridge(td, deck), read_inp(run.workdir / f"{case}.inp"), run


def test_bridge_scalars(bridged_case):
    """Every scalar DEGBRIDG writes into the ground-level deck.

    The touchdown distance becomes OODIST, an offset added to every reported
    distance so the two stages share one coordinate; the release temperature
    becomes the jet temperature; the density is the release density, not the
    ambient one.
    """
    _case, got, ref, _run = bridged_case
    pairs = [
        ("u0", got.u0, ref.u0), ("z0", got.z0, ref.z0), ("zr", got.zr, ref.zr),
        ("oodist", got.oodist, ref.oodist), ("avtime", got.avtime, ref.avtime),
        ("rml", got.rml, ref.rml),
        ("tamb", got.ambient.tamb, ref.ambient.tamb),
        ("pamb", got.ambient.pamb, ref.ambient.pamb),
        ("humid", got.ambient.humid, ref.ambient.humid),
        ("tsurf", got.ambient.tsurf, ref.ambient.tsurf),
        ("gasmw", got.gas.mw, ref.gas.mw),
        ("gastem", got.gas.temp, ref.gas.temp),
        ("gasrho", got.gas.rho, ref.gas.rho),
        ("cpk", got.gas.cpk, ref.gas.cpk), ("cpp", got.gas.cpp, ref.gas.cpp),
        ("yclow", got.yclow, ref.yclow), ("gmass0", got.gmass0, ref.gmass0),
    ]
    for name, g, r in pairs:
        assert (_rel(g, r) if r != 0 else abs(g)) < 1e-6, name
    assert (got.istab, got.ambient.isofl, got.ambient.ihtfl) == (
        ref.istab, ref.ambient.isofl, ref.ambient.ihtfl
    )
    assert got.steady_state is ref.steady_state is True


def test_bridge_source_table(bridged_case):
    """The equivalent area source.

    A circular source of the plume's landed half-width, releasing the same
    mass rate but already diluted to the concentration the plume had when it
    landed -- which is the whole point of the bridge. Restarting from pure
    contaminant would over-predict everything downwind.
    """
    _case, got, ref, _run = bridged_case
    a = np.column_stack([
        got.source.time, got.source.rate, got.source.radius, got.source.wc,
        got.source.temp, got.source.fracv, got.source.enthalpy, got.source.rho,
    ])
    b = np.column_stack([
        ref.source.time, ref.source.rate, ref.source.radius, ref.source.wc,
        ref.source.temp, ref.source.fracv, ref.source.enthalpy, ref.source.rho,
    ])
    assert a.shape == b.shape == (4, 8)
    dev = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
    assert dev.max() < 1e-6
    # the source composition must be dilute, not pure contaminant
    assert 0.0 < got.source.wc[0] < 0.01


def test_bridge_density_table(bridged_case):
    """The jet deck's NDEN flag selects the density treatment; the bridge
    carries the choice across unchanged."""
    _case, got, ref, _run = bridged_case
    assert got.density_table is not None
    assert got.density_table.shape == ref.density_table.shape
    assert np.allclose(got.density_table, ref.density_table, rtol=1e-12)


def test_bridge_refuses_a_plume_that_never_lands(reference, tmp_path_factory):
    """EX1's plume stays aloft, and JETPLU signals that with a zero distance."""
    from degali.io.bridge import bridge, read_in
    from degali.io.jetdeck import Touchdown

    run = reference.run("ex1", tmp_path_factory.mktemp("ex1b"))
    deck = read_in(run.workdir / "ex1.in")
    with pytest.raises(ValueError, match="never reached the ground"):
        bridge(Touchdown(0.0, 0.0, 0.0), deck)


def test_jet_to_ground_pipeline(bridged_case):
    """EX2 and EX3 end to end in Python: jet, bridge, source, downwind.

    Nothing in this test touches the Fortran except to read the reference
    profile at the end. It is the first point at which the port replaces the
    original for a pressurised release.
    """
    import math

    from degali.core.atmosphere import fit_alpha, psi
    from degali.core.crfg import build_source_vectors
    from degali.core.downwind import Downwind
    from degali.core.driver import SourceRun
    from degali.core.numerics import gamma as gamma_fn
    from degali.core.steady import SteadyStateRun
    from degali.io.bridge import bridge
    from degali.io.jetdeck import Touchdown, read_ino
    from degali.io.params import alph_settings, read_er1, read_er2

    case, _got, _ref, run = bridged_case
    w = run.workdir

    # 1. the jet
    deck = read_ino(w / f"{case}.ino")
    jp, y0 = _build_jet(deck)
    jet = jp.run(y0, distmx=deck.distmx)
    assert jet.touchdown

    # 2. the bridge
    from degali.io.bridge import read_in

    td = Touchdown(jet.distance, jet.concentration, jet.halfwidth)
    c = bridge(td, read_in(w / f"{case}.in"))

    # 3. the source
    params, raw = read_er1(w / f"{case}.er1")
    er2 = read_er2(w / f"{case}.er2")
    th = c.make_thermo()
    sr = SourceRun(c, th, params)
    alpha = fit_alpha(
        c.u0, c.z0, c.zr, c.rml, ustar=sr.ustar, legacy=True,
        **alph_settings(raw),
    )
    sr.set_alpha(alpha)
    result = sr.run()
    vectors = build_source_vectors(sr.raw_records, crfger=raw["crfger"])
    assert len(vectors) >= 2

    # 4. downwind
    th.ambient.humsrc = (
        1.0 - result.swcl - result.swal * (1.0 + c.ambient.humid)
    ) / result.swcl
    if c.ambient.isofl != 1:
        th.build_adiabatic_table(result.swcl, result.swal, result.senl)
    dw = Downwind(
        c, th, params, alpha=alpha, gammaf=gamma_fn(1.0 / (alpha + 1.0)),
        ustar=sr.ustar, deltay=c.stability.deltay, betay=c.stability.betay,
        rhoa=th.table.rhoa,
    )
    handoff = {
        "ess": result.ess, "outl": result.outl, "outb": result.outb,
        "outsz": result.outsz, "outcc": result.outcc, "swcl": result.swcl,
        "swal": result.swal, "senl": result.senl, "srhl": result.srhl,
        "yclow": c.yclow, "rhoe": th.table.rhoe,
    }
    profile = SteadyStateRun(dw, handoff, er2, oodist=c.oodist).run()
    assert len(profile.rows) > 20

    # against the Fortran's own profile for the same case
    text = run.golden.replace("\r", "")
    block = text.split("(m)                (kg/m**3)")[1].split("___")[0]
    rows = []
    for line in block.splitlines():
        parts = line.split()
        if len(parts) >= 9 and parts[0][0].isdigit():
            try:
                rows.append([float(x) for x in parts[:9]])
            except ValueError:
                pass
    ref_rows = np.array(rows)
    assert len(ref_rows) > 20

    inside = (ref_rows[:, 0] >= profile.rows[0, 0]) & (
        ref_rows[:, 0] <= profile.rows[-1, 0]
    )
    assert inside.sum() > len(ref_rows) // 2
    interp = np.interp(
        ref_rows[inside, 0], profile.rows[:, 0], profile.rows[:, 1]
    )
    dev = np.abs(interp - ref_rows[inside, 1]) / ref_rows[inside, 1]
    assert dev.max() < 0.03, f"{case}: max rel dev {dev.max():.3e}"


# ==========================================================================
# the high-level API and the command line
# ==========================================================================


def test_run_steady_matches_the_staged_path(b9, b9_profile):
    """The one-call API must give the same answer as driving the stages.

    Not to round-off: the staged fixture seeds the downwind model from the
    ``.TR2`` handoff the Fortran wrote, which carries seven significant
    figures, while ``run_steady`` passes the source result straight through at
    full precision. The one-call path is the more accurate of the two.
    """
    from degali.run import run_steady

    run, _p6 = b9
    staged, _er2 = b9_profile
    profile, src = run_steady(
        run.workdir / "b9.inp",
        er1=run.workdir / "b9.er1",
        er2=run.workdir / "b9.er2",
    )
    assert len(profile.rows) == len(staged.rows)
    assert profile.mass_above_lfl == pytest.approx(
        staged.mass_above_lfl, rel=1e-5
    )
    assert src.alpha == pytest.approx(run["alpha"], rel=1e-11)


def test_run_steady_uses_shipped_defaults(b9):
    """Omitting the parameter files must fall back to EPA's own values."""
    from degali.run import run_steady

    run, _p6 = b9
    with_files, _s1 = run_steady(
        run.workdir / "b9.inp",
        er1=run.workdir / "b9.er1",
        er2=run.workdir / "b9.er2",
    )
    defaults, _s2 = run_steady(run.workdir / "b9.inp")
    assert defaults.mass_above_lfl == pytest.approx(
        with_files.mass_above_lfl, rel=1e-9
    )


def test_run_transient_and_dose(b9t):
    """A whole transient run, then receptor histories from the same field."""
    from degali.run import Receptor, run_transient

    run, _p6 = b9t
    out = run_transient(
        run.workdir / "b9t.inp",
        er1=run.workdir / "b9t.er1",
        er2=run.workdir / "b9t.er2",
    )
    assert len(out.field.observers) == 30
    assert len(out.snapshots) == 19

    histories = out.dose([Receptor(x=200.0), Receptor(x=800.0)])
    assert len(histories) == 2
    near, far = (h.peak[0] for h in histories)
    assert near > far > 0.0
    # the peak at a distance must agree with the snapshot view of it
    from_field, _t = out.field.max_concentration_at(800.0)
    assert far == pytest.approx(from_field, rel=0.02)


def test_run_jet_and_bridge(reference, tmp_path_factory):
    """run_jet on a plume that stays aloft, run_jet_to_ground on one that lands."""
    from degali.run import run_jet, run_jet_to_ground

    aloft = reference.run("ex1", tmp_path_factory.mktemp("ex1r"))
    jet, deck = run_jet(aloft.workdir / "ex1.ino")
    assert not jet.touchdown
    assert deck.erate > 0.0

    lands = reference.run("ex2", tmp_path_factory.mktemp("ex2r"))
    profile, jet2, src = run_jet_to_ground(
        lands.workdir / "ex2.ino", lands.workdir / "ex2.in",
        er1=lands.workdir / "ex2.er1", er2=lands.workdir / "ex2.er2",
    )
    assert jet2.touchdown
    # every reported distance is offset by the touchdown point so the two
    # stages share one coordinate
    assert src.case.oodist == pytest.approx(jet2.distance, rel=1e-12)
    assert profile.rows[0, 0] > jet2.distance


def test_run_jet_to_ground_refuses_a_plume_that_stays_aloft(
    reference, tmp_path_factory
):
    from degali.run import run_jet_to_ground

    run = reference.run("ex1", tmp_path_factory.mktemp("ex1s"))
    with pytest.raises(ValueError, match="never reached the ground"):
        run_jet_to_ground(run.workdir / "ex1.ino", run.workdir / "ex1.in")


@pytest.mark.parametrize(
    "argv_key",
    ["steady", "transient", "dose", "jet_aloft", "jet_bridged"],
)
def test_cli(b9, b9t, reference, tmp_path_factory, capsys, argv_key):
    """Each subcommand runs and prints something specific to it."""
    from degali.cli import main

    b9_run, _ = b9
    b9t_run, _ = b9t
    ex1 = reference.run("ex1", tmp_path_factory.mktemp("ex1c"))
    ex2 = reference.run("ex2", tmp_path_factory.mktemp("ex2c"))

    argvs = {
        "steady": (
            ["steady", str(b9_run.workdir / "b9.inp"),
             "--er1", str(b9_run.workdir / "b9.er1"),
             "--er2", str(b9_run.workdir / "b9.er2")],
            "lower level of concern",
        ),
        "transient": (
            ["transient", str(b9t_run.workdir / "b9t.inp"),
             "--er1", str(b9t_run.workdir / "b9t.er1"),
             "--er2", str(b9t_run.workdir / "b9t.er2"),
             "--snapshot", "60", "--snapshot", "120"],
            "peak flammable mass",
        ),
        "dose": (
            ["dose", str(b9t_run.workdir / "b9t.inp"),
             "--er1", str(b9t_run.workdir / "b9t.er1"),
             "--er2", str(b9t_run.workdir / "b9t.er2"), "--at", "400"],
            "dose (mol frac.s)",
        ),
        "jet_aloft": (
            ["jet", str(ex1.workdir / "ex1.ino")],
            "never reaches the ground",
        ),
        "jet_bridged": (
            ["jet", str(ex2.workdir / "ex2.ino"),
             "--bridge", str(ex2.workdir / "ex2.in"),
             "--er1", str(ex2.workdir / "ex2.er1"),
             "--er2", str(ex2.workdir / "ex2.er2")],
            "bridged to a",
        ),
    }
    argv, expected = argvs[argv_key]
    assert main(argv) == 0
    assert expected in capsys.readouterr().out


def test_cli_reports_errors_without_a_traceback(reference, tmp_path_factory, capsys):
    from degali.cli import main

    run = reference.run("ex1", tmp_path_factory.mktemp("ex1e"))
    code = main([
        "jet", str(run.workdir / "ex1.ino"),
        "--bridge", str(run.workdir / "ex1.in"),
    ])
    assert code == 1
    assert "never reached the ground" in capsys.readouterr().err


# ==========================================================================
# the CoolProp backend and the modern numerics
# ==========================================================================


def test_fluid_resolver():
    """Deck labels and molecular weights map onto CoolProp fluids.

    The weight match is deliberately narrow: naming the wrong fluid is worse
    than naming none, because the 1989 correlations are at least
    approximations of the right substance.
    """
    from degali.core.fluids import resolve

    assert resolve("LNG", 16.04).fluid == "Methane"
    assert resolve("LNG", 16.04).how == "alias"
    assert resolve("nh3", 17.03).fluid == "Ammonia"
    assert resolve("Methane", None).fluid == "Methane"  # already canonical

    by_weight = resolve("ZZZ", 44.01)
    assert by_weight.fluid == "CarbonDioxide"
    assert by_weight.how == "molecular weight"

    unmatched = resolve("ZZZ", 123.4)
    assert unmatched.fluid is None
    assert not unmatched
    assert resolve(None, None).fluid is None


def test_legacy_backend_is_the_default(b9, b9_case):
    """Asking for nothing must give the 1989 behaviour, unchanged."""
    from degali.core.thermo import LegacyBackend

    run, _p6 = b9
    th = b9_case.make_thermo()
    assert isinstance(th.backend, LegacyBackend)
    assert th.legacy_numerics is True
    th.reference_enthalpies()
    assert _rel(th.hmrte, run["hmrte"]) < RTOL


def test_coolprop_backend_selection(b9_case):
    """Selecting CoolProp by name resolves the contaminant and switches the
    numerics with it: someone asking for real properties has not also asked
    for the original's root-finder tolerances."""
    pytest.importorskip("CoolProp")
    from degali.core.thermo import CoolPropBackend

    th = b9_case.make_thermo(backend="coolprop")
    assert isinstance(th.backend, CoolPropBackend)
    assert th.legacy_numerics is False
    assert b9_case.gas.coolprop_name == "Methane"
    # the contaminant enthalpy moves: the deck's two fitted heat-capacity
    # constants are a poor fit for methane near saturation
    legacy = b9_case.make_thermo()
    legacy.reference_enthalpies()
    th.reference_enthalpies()
    assert abs(th.hmrte - legacy.hmrte) / abs(legacy.hmrte) > 0.01


def test_modern_numerics_can_be_forced_on_the_legacy_backend(b9_case):
    """The property correlations and the numerics are separate choices."""
    th = b9_case.make_thermo(legacy_numerics=False)
    from degali.core.thermo import LegacyBackend

    assert isinstance(th.backend, LegacyBackend)
    assert th.legacy_numerics is False


def test_heat_capacity_is_continuous_through_condensation(b9_case):
    """Mixture enthalpy has a slope kink where water starts to condense.

    A derivative jumps there; a secant over a finite interval does not,
    because enthalpy itself is continuous. Narrowing the interval towards a
    derivative recovers the jump, and the downwind integrator then bisects to
    a standstill -- which is how the interval width was chosen.
    """
    th = b9_case.make_thermo()
    th.reference_enthalpies()
    th.build_adiabatic_table(1.0, 0.0, th.hmrte)
    wc, wa = 0.02, 0.96
    temps = np.linspace(250.0, 300.0, 400)
    cp = np.array([th.heat_capacity(wc, wa, t) for t in temps])
    assert np.all(cp > 0.0)
    # no step changes: successive values differ by far less than cp itself
    jumps = np.abs(np.diff(cp)) / cp[:-1]
    assert jumps.max() < 0.05


def test_addheat_conditioning_is_fixed_by_accurate_inversion(b9_case):
    """The one place the legacy path loses accuracy.

    ADDHEAT forms cp as dh/(T - T_adiabatic) while resolving T to only 1e-3 K.
    Early in a dense plume that denominator is a few hundredths of a kelvin,
    so the quotient carries percent-level noise. Resolving T accurately makes
    the same secant well conditioned, and the result stops depending on where
    the root finder happened to stop.
    """
    legacy = b9_case.make_thermo()
    legacy.reference_enthalpies()
    legacy.build_adiabatic_table(1.0, 0.0, legacy.hmrte)
    modern = b9_case.make_thermo(legacy_numerics=False)
    modern.reference_enthalpies()
    modern.table = legacy.table

    cc = float(legacy.table.cc[3])
    # a small heat addition: the regime where the quotient is ill conditioned
    small = [legacy.add_heat(cc, dh).cp for dh in (40.0, 40.5, 41.0)]
    smooth = [modern.add_heat(cc, dh).cp for dh in (40.0, 40.5, 41.0)]
    # the modern heat capacities vary smoothly with the heat added
    assert max(smooth) / min(smooth) < 1.02
    # and all values remain physical
    assert all(c > 0.0 for c in small + smooth)


@pytest.mark.slow
def test_coolprop_runs_every_case_type(b9, b9t, reference, tmp_path_factory):
    """The CoolProp path has to survive the whole model, not just the tables.

    Structure must be preserved -- the same observers in range at the same
    instants -- while the numbers move by a few per cent.
    """
    pytest.importorskip("CoolProp")
    from degali.run import run_jet_to_ground, run_steady, run_transient

    b9_run, _ = b9
    steady, src = run_steady(
        b9_run.workdir / "b9.inp",
        er1=b9_run.workdir / "b9.er1", er2=b9_run.workdir / "b9.er2",
        backend="coolprop",
    )
    legacy, _ = run_steady(
        b9_run.workdir / "b9.inp",
        er1=b9_run.workdir / "b9.er1", er2=b9_run.workdir / "b9.er2",
    )
    assert src.case.gas.coolprop_name == "Methane"
    for got, want in (
        (steady.distance_to(0.05), legacy.distance_to(0.05)),
        (steady.mass_above_lfl, legacy.mass_above_lfl),
    ):
        assert 0.9 < got / want < 1.1  # moves, but not by much

    b9t_run, _ = b9t
    out = run_transient(
        b9t_run.workdir / "b9t.inp",
        er1=b9t_run.workdir / "b9t.er1", er2=b9t_run.workdir / "b9t.er2",
        backend="coolprop",
    )
    assert len(out.snapshots) == 19
    assert all(len(o.rows) > 10 for o in out.field.observers)


def test_jet_decks_ignore_the_backend(reference, tmp_path_factory):
    """A jet deck supplies its own density table.

    ``ISOFL = 1`` means no mixing line is built, so the backend has nothing to
    change. Identical results are the correct outcome here, not a sign the
    selection was dropped.
    """
    pytest.importorskip("CoolProp")
    from degali.run import run_jet

    run = reference.run("ex2", tmp_path_factory.mktemp("ex2cp"))
    legacy, deck = run_jet(run.workdir / "ex2.ino")
    modern, _ = run_jet(run.workdir / "ex2.ino", backend="coolprop")
    assert deck.isofl == 1
    assert modern.distance == pytest.approx(legacy.distance, rel=1e-12)
    assert modern.concentration == pytest.approx(legacy.concentration, rel=1e-12)


def test_cli_backend_option(b9, capsys):
    """--backend switches the property model and says which one it used."""
    pytest.importorskip("CoolProp")
    from degali.cli import main

    run, _p6 = b9
    base = [str(run.workdir / "b9.inp"),
            "--er1", str(run.workdir / "b9.er1"),
            "--er2", str(run.workdir / "b9.er2")]

    assert main(["steady", *base]) == 0
    legacy_out = capsys.readouterr().out
    assert "properties: legacy" in legacy_out

    assert main(["steady", *base, "--backend", "coolprop"]) == 0
    modern_out = capsys.readouterr().out
    assert "properties: CoolProp (Methane)" in modern_out
    assert legacy_out != modern_out


def test_coolprop_grid_survives_a_two_phase_gap():
    """A cryogenic release sits below its own saturation line at ambient
    pressure, so part of the tabulated span has no single-phase data.

    The nearest valid value is held across the gap. Abandoning the grid
    instead -- which one failed point used to do -- silently dropped the whole
    run back to the 1989 correlation while still reporting success.
    """
    pytest.importorskip("CoolProp")
    from degali.core.thermo import CoolPropBackend

    b = CoolPropBackend("Methane")
    # 111.7 K at 0.94 atm is below methane's saturation temperature
    cp = b.cp_contaminant_eos(111.7, 0.94)
    assert cp is not None and 1500.0 < cp < 3000.0
    rho = b.rho_contaminant_eos(111.7, 0.94)
    assert rho == pytest.approx(1.7029, rel=1e-3)
    # and the real-gas density differs from the deck's ideal-gas value
    assert abs(rho - 1.68448) / 1.68448 > 0.005


def test_coolprop_grid_is_built_once():
    """Growing the grid outward as queries arrive costs hundreds of rebuilds,
    because the model walks its temperature range gradually."""
    pytest.importorskip("CoolProp")
    from degali.core.thermo import CoolPropBackend

    b = CoolPropBackend()
    for t in np.linspace(120.0, 380.0, 200):
        b.cp_air(float(t))
    lo, hi, xs, _ys = b._grids["cpa"]
    assert (lo, hi) == CoolPropBackend.SPAN
    assert len(xs) == CoolPropBackend.POINTS


# ==========================================================================
# field-trial validation machinery
# ==========================================================================

# The REDIPHEM database is not redistributable, so these tests run only when
# $REDIPHEM_ROOT points at an unpacked copy.
from degali.validation.rediphem import ENV_VAR, default_root

REDIPHEM_ROOT = default_root()
needs_rediphem = pytest.mark.skipif(
    REDIPHEM_ROOT is None, reason=f"set {ENV_VAR} to run the field-trial tests"
)


def test_statistics_sign_convention():
    """MG is observed over predicted, so a model reading high gives MG < 1.

    Getting this backwards is easy and the published tables are not always
    explicit, so the convention is pinned here.
    """
    from degali.validation.statistics import statistics

    obs = np.array([1.0, 2.0, 4.0])
    high = statistics(obs, obs * 2.0)
    assert high.mg == pytest.approx(0.5)
    assert "over-predicts" in high.reads
    low = statistics(obs, obs / 2.0)
    assert low.mg == pytest.approx(2.0)
    assert "under-predicts" in low.reads

    perfect = statistics(obs, obs)
    assert perfect.mg == pytest.approx(1.0)
    assert perfect.vg == pytest.approx(1.0)
    assert perfect.fac2 == 1.0
    assert perfect.acceptable


def test_statistics_needs_positive_pairs():
    from degali.validation.statistics import statistics

    with pytest.raises(ValueError):
        statistics([1.0, 2.0], [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        statistics([0.0, 0.0], [1.0, 2.0])


def test_running_peak_beats_block_averaging():
    """A cloud passes one sensor in far less time than a release lasts.

    Fixed blocks straddle the arrival and dilute the peak with the empty
    record either side; a running mean does not. On Burro 9 the difference
    moves the geometric variance from 1.4 to 10.8, which would read as a model
    failure rather than as an artefact of the reduction.
    """
    from degali.validation.rediphem import _running_peak

    series = np.zeros(300)
    series[140:170] = 10.0  # a 30 s passage in a 300 s record
    assert _running_peak(series, 1) == pytest.approx(10.0)
    assert _running_peak(series, 30) == pytest.approx(10.0)
    # a window far longer than the passage necessarily dilutes it
    assert _running_peak(series, 120) == pytest.approx(2.5)
    # gaps are skipped rather than counted as zero
    gappy = series.copy()
    gappy[150:155] = np.nan
    assert _running_peak(gappy, 30) == pytest.approx(10.0)


@needs_rediphem
def test_rediphem_reader():
    """DATA.DBF is raw float32, not the dBase file its extension claims."""
    from degali.validation.rediphem import load

    trials = load(REDIPHEM_ROOT)
    assert len(trials) > 300
    b9 = next(t for t in trials if t.series == "BURRO" and t.name == "B9")
    assert b9.substance == "lng"
    assert b9.release_type == "pool"
    assert b9.rate == pytest.approx(136.0)
    assert b9.stability == 4  # D

    ids, time, values = b9.data()
    assert ids.size == 216
    assert values.shape == (time.size, ids.size)
    assert np.median(np.diff(time)) == pytest.approx(1.0)
    assert np.isnan(values).any()  # -1234 became nan


@needs_rediphem
def test_channel_types_are_read_not_assumed():
    """Concentration channel numbers differ between series.

    Burro uses 21 and 22, Eagle 21 and 24. Hardcoding them drops whole series
    silently -- Eagle produced no comparisons at all until this was fixed.
    """
    from degali.validation.rediphem import load

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    assert trials[("BURRO", "B9")].channel_types() == {21, 22}
    assert trials[("EAGLE", "E2")].channel_types() == {21, 24}
    for key in (("BURRO", "B9"), ("EAGLE", "E2"), ("TORTOISE", "D1")):
        assert len(trials[key].arc_maxima()) > 0


@needs_rediphem
def test_trial_to_case_records_its_assumptions():
    """A trial specification is not a model deck; the gap must be explicit."""
    from degali.validation.rediphem import load
    from degali.validation.trialcase import to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    tc = to_case(trials[("BURRO", "B9")])
    assert tc.usable
    c = tc.case
    assert c.u0 == pytest.approx(5.7)
    assert c.z0 == pytest.approx(2.0)
    assert c.zr == pytest.approx(0.00025)
    assert c.istab == 4
    assert c.source.rate[0] == pytest.approx(136.0)
    assert c.source.radius[0] == pytest.approx(29.0)  # 58 m recorded diameter
    assert c.gas.name == "LNG"
    for expected in ("surface temperature", "release temperature", "source radius"):
        assert expected in tc.assumptions

    # a substance with no property data is reported, not guessed at
    unknown = trials[("TNO", "TUV01")]
    assert not to_case(unknown).usable


@needs_rediphem
@pytest.mark.slow
def test_burro_evaluation():
    """The LNG pool spills, end to end against measurement.

    This is the check that the whole chain -- reader, deck builder, model,
    statistics -- produces a defensible number rather than merely running.
    DEGADIS reads about 20 % high on Burro with a factor-of-two hit rate above
    a half, which is in line with the published evaluations.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load

    trials = load(REDIPHEM_ROOT)
    result = compare_series(trials, "BURRO", height=1.0, averaging_time=18.4)
    assert len(result.comparisons) == 8
    assert not result.skipped

    s = result.statistics()
    assert s.n > 50
    assert 0.6 < s.mg < 1.0     # reads high, but not by much
    assert s.fac2 > 0.5         # Hanna's factor-of-two criterion


@needs_rediphem
def test_equivalent_source_route_for_flashing_jets():
    """A flashing jet cannot go through the pool source.

    A pressurised release of a liquefied gas flashes, throws an aerosol and
    entrains air while the droplets evaporate. None of that is in DEGADIS.
    Forcing Desert Tortoise through the pool source assumes a 293 m pool for
    an 81 kg/s ammonia jet; the SMEDIS equivalent source hands the model the
    plume state where the liquid has gone instead.
    """
    from degali.validation.compare import compare
    from degali.validation.rediphem import load
    from degali.validation.trialcase import SMEDIS_SOURCES, to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    d1 = trials[("TORTOISE", "D1")]

    # the pool route: a nonsense source, and it says so
    pool = to_case(d1)
    assert pool.case.source.radius[0] > 200.0
    assert "pressurised jet" in pool.assumptions["source"]

    # the equivalent-source route is picked automatically where one exists
    comparison = compare(d1, height=1.0, averaging_time=18.4)
    assert comparison is not None
    assert "SMEDIS equivalent source" in comparison.case.assumptions["source"]
    assert "momentum" in comparison.case.assumptions
    assert comparison.case.case.oodist == pytest.approx(51.0)
    assert comparison.case.case.source.radius[0] == pytest.approx(6.40)
    # and it brings the prediction inside a factor of two
    assert 0.5 <= comparison.ratio[0] <= 2.0

    assert ("TORTOISE", "D3") not in SMEDIS_SOURCES  # falls back to the pool


@needs_rediphem
def test_missing_fields_are_filled_from_sibling_trials():
    """Desert Tortoise D2 records no anemometer height; D1, D3 and D4 all say
    2 m. Dropping the trial over that loses a quarter of the series to a
    transcription gap rather than to anything physical."""
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load
    from degali.validation.trialcase import to_case

    trials = load(REDIPHEM_ROOT)
    d2 = next(t for t in trials if t.series == "TORTOISE" and t.name == "D2")
    assert "reference height missing" in to_case(d2).problems

    result = compare_series(trials, "TORTOISE", height=1.0, averaging_time=18.4)
    assert not result.skipped
    assert len(result.comparisons) == 4


@needs_rediphem
def test_model_declines_to_predict_inside_its_own_source():
    """Desert Tortoise instruments an arc at 100 m; with the equivalent source
    at 51 m the secondary source is 103 m long, so the first prediction is at
    103 m. Having no prediction there is right -- the arc is inside the source
    -- and the harness drops the pair rather than extrapolating."""
    from degali.validation.compare import compare
    from degali.validation.rediphem import load

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    d1 = trials[("TORTOISE", "D1")]
    assert set(d1.arc_maxima(height=1.0)) == {100.0, 800.0}
    comparison = compare(d1, height=1.0, averaging_time=18.4)
    assert list(comparison.arcs) == [800.0]


@needs_rediphem
def test_computed_source_reproduces_the_published_one():
    """The flash calculation is checked against source terms it did not see.

    SMEDIS published equivalent sources for Desert Tortoise: 13 mole per cent
    at 205 K, half-width 6.4 m. Computing one from the storage conditions and
    an energy balance gives 13.6 per cent at 204.9 K and 5.3 m. That agreement
    is what licenses applying the same calculation to FLADIS, where no source
    term was issued.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.flashing import equivalent_source, plume_half_width

    result = equivalent_source(
        "Ammonia", storage_temperature=294.65, ambient_temperature=302.45,
        ambient_pressure=0.909e5, molecular_weight=17.03,
    )
    assert result.temperature == pytest.approx(205.0, abs=1.0)
    assert result.molar_percent == pytest.approx(13.0, abs=1.0)
    assert plume_half_width(81.0, result, 7.5) == pytest.approx(6.4, rel=0.25)
    # the equivalent-source temperature is far below the pure boiling point,
    # because the plume is dilute and the partial pressure is a fraction of
    # ambient
    assert result.temperature < 240.0


def test_equivalent_source_spans_ammonia_to_hydrogen():
    """One calculation has to cover substances four hundred kelvin apart.

    Ammonia's critical temperature is 405 K, above ambient, so its whole
    liquid range is available and the equivalent source lands near 205 K.
    Hydrogen's is 33 K, far *below* ambient: there is no saturation line
    above it, so the search is capped at the critical point and a solution
    there means the aerosol evaporated before the mixture warmed that far.

    A fixed temperature floor cannot serve both -- any value that excludes
    nonsense for ammonia excludes the answer for hydrogen -- so the fluid's
    own triple and critical points set the bracket.
    """
    pytest.importorskip("CoolProp")
    from CoolProp.CoolProp import PropsSI

    from degali.validation.flashing import equivalent_source

    ammonia = equivalent_source(
        "Ammonia", storage_temperature=294.65, ambient_temperature=302.45,
        ambient_pressure=0.909e5, molecular_weight=17.03,
    )
    assert ammonia.temperature == pytest.approx(205.0, abs=1.0)

    boil = PropsSI("T", "P", 6.013e5, "Q", 0, "Hydrogen")
    hydrogen = equivalent_source(
        "Hydrogen", storage_temperature=boil, ambient_temperature=288.6,
        ambient_pressure=101325.0, molecular_weight=2.016,
    )
    assert hydrogen.temperature < PropsSI("Tcrit", "Hydrogen")
    assert 15.0 < hydrogen.temperature < 33.0
    assert 80.0 < hydrogen.molar_percent < 99.0
    # the flash is a much larger fraction than for ammonia, and the mixture
    # it leaves is denser than air despite being hydrogen
    assert hydrogen.flash_fraction > 0.1
    assert hydrogen.density > 1.2


def test_phase_safe_hydrogen_source_skips_the_condensed_bulk_air_zone():
    """The Gaussian jet must not start with impossible 20 K gaseous air.

    The optional plug-flow boundary conserves component enthalpy until the
    N2/O2/Ar mixture reaches its first all-gas state.  The temperature is not
    fitted: oxygen's partial pressure closes on its saturation pressure.
    """
    pytest.importorskip("CoolProp")
    from CoolProp.CoolProp import PropsSI

    from degali.validation.flashing import equivalent_source

    stored = PropsSI("T", "P", 6.013e5, "Q", 0, "Hydrogen")
    old = equivalent_source(
        "Hydrogen", storage_temperature=stored,
        ambient_temperature=288.6, ambient_pressure=101325.0,
        molecular_weight=2.016,
    )
    safe = equivalent_source(
        "Hydrogen", storage_temperature=stored,
        ambient_temperature=288.6, ambient_pressure=101325.0,
        molecular_weight=2.016, bulk_air_phase_safe=True,
    )

    assert old.temperature < 21.0
    assert 67.0 < safe.temperature < 70.0
    assert safe.temperature > PropsSI("Tcrit", "Hydrogen")
    assert safe.air_ratio > 3.0 * old.air_ratio
    assert 0.19 < safe.mass_fraction < 0.22
    assert safe.bulk_air_phase_safe
    assert safe.limiting_air_species == "Oxygen"
    assert safe.max_air_saturation_ratio == pytest.approx(1.0, abs=1.0e-9)


def test_phase_safe_source_is_not_silently_applied_to_other_fluids():
    pytest.importorskip("CoolProp")
    from degali.validation.flashing import equivalent_source

    with pytest.raises(ValueError, match="implemented for Hydrogen only"):
        equivalent_source(
            "Ammonia", storage_temperature=294.65,
            ambient_temperature=302.45, ambient_pressure=0.909e5,
            molecular_weight=17.03, bulk_air_phase_safe=True,
        )


def test_equivalent_source_refuses_a_gas_with_no_liquid_range():
    pytest.importorskip("CoolProp")
    from degali.validation.flashing import equivalent_source

    with pytest.raises(ValueError, match="sub-critical liquid range"):
        equivalent_source(
            "Helium", storage_temperature=4.5, ambient_temperature=293.0,
            ambient_pressure=1.0e5, molecular_weight=4.003,
            minimum_temperature=100.0,
        )


@needs_rediphem
def test_thorney_island_hits_a_real_model_limit():
    """The instantaneous route reaches DEGADIS's own failure, not the port's.

    A Thorney Island release is a 14 m cylinder with a height-to-diameter
    ratio near one, which turns on the van Ulden momentum balance -- a branch
    no EPA test case exercises. The balance does not converge for that shape.
    The original Fortran stops in the same place with "STOP SRC1 velocity
    loop", so this is a limit of the model.
    """
    from degali.run import run_transient
    from degali.validation.rediphem import load
    from degali.validation.trialcase import puff_to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    tc = puff_to_case(trials[("TI", "TI08")])
    assert tc.usable
    assert tc.case.instantaneous is True
    assert tc.case.gmass0 > 0.0
    assert np.all(tc.case.source.rate == 0.0)  # a puff has no flux
    # the cylinder held a mixture, not pure contaminant
    assert 0.4 < tc.case.source.wc[0] < 0.6

    with pytest.raises(RuntimeError, match="DEGADIS 2.1 stops here too"):
        run_transient(tc.case)


def test_degenerate_secondary_source_is_reported():
    """A source model that produces nothing must say so, not divide by zero."""
    from degali.core.steady import SteadyStateRun

    class _Stub:
        alpha1 = 1.1
        p = None

        class case:
            u0, z0 = 5.0, 10.0

    run = SteadyStateRun.__new__(SteadyStateRun)
    run.dw = _Stub()
    run.h = {"outl": 0.0, "outb": 0.0, "outsz": 1.0, "outcc": 0.0}
    run.er2 = {"sy0er": 0.0}
    run.erate = 1.0
    with pytest.raises(ValueError, match="degenerate secondary source"):
        run._initial()


@needs_rediphem
def test_channel_numbering_is_never_guessed():
    """A series with no CHANDEF.DAT yields nothing, and that is correct.

    Channel type numbers are per series -- Burro uses 21 and 22, Eagle 21 and
    24. FLADIS ships no definitions, and a plausible fallback list selects
    channels reading 302 and 23.5 at 20 m downwind, which are not
    concentrations. That produced a clean-looking FLADIS result that meant
    nothing, so the reader now refuses rather than guesses.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load

    trials = load(REDIPHEM_ROOT)
    fladis = next(t for t in trials if t.series == "FLADIS")
    assert not (fladis.path.parent / "CHANDEF.DAT").exists()
    assert fladis.channel_types() == set()
    assert fladis.arc_maxima() == {}

    result = compare_series(trials, "FLADIS", averaging_time=18.4)
    assert not result.comparisons
    assert all("no usable measurements" in v for v in result.skipped.values())


@needs_rediphem
@pytest.mark.slow
def test_whole_cloud_comparison_fails_where_one_height_passes():
    """Reporting only the lowest height flatters the model.

    DEGADIS gets the ground-level centreline roughly right and the vertical
    distribution badly wrong, so a single-height statistic is carried by the
    one elevation where the two errors cancel. Pooling every instrumented
    height is a different quantity, not merely a noisier one.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load
    from degali.validation.statistics import statistics

    trials = load(REDIPHEM_ROOT)
    result = compare_series(trials, "BURRO", averaging_time=18.4)
    lowest = min(z for c in result.comparisons for z in c.heights)
    obs = np.concatenate([c.observed[c.heights == lowest] for c in result.comparisons])
    pred = np.concatenate([c.predicted[c.heights == lowest] for c in result.comparisons])

    single = statistics(obs, pred)
    pooled = result.statistics()
    assert 0.6 < single.mg < 1.0
    assert pooled.mg > 10.0
    assert pooled.vg > single.vg * 1e6


@needs_rediphem
def test_release_clock_aligns_with_the_data():
    """The release time is three lines, only the first of which is labelled.

    Reading the labelled line alone gives the hour and silently loses 37
    minutes, which puts every arrival time out by that much.
    """
    from degali.validation.rediphem import load

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    for key in (("BURRO", "B9"), ("TORTOISE", "D1")):
        t = trials[key]
        _ids, time, _v = t.data()
        assert abs(t.start_time - time[0]) < 5.0, key


@needs_rediphem
def test_mast_profiles_do_not_mix_positions():
    """A vertical profile has to come from one mast.

    Taking the arc maximum at each height separately mixes masts, and can put
    the highest reading at 8 m and the lowest at 3 m simply because different
    sensors were nearest the plume centre. On Burro 9 that artefact makes the
    cloud look like it grows with height.
    """
    from degali.validation.rediphem import load

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    b9 = trials[("BURRO", "B9")]

    masts = b9.mast_profiles()
    assert masts
    for (x, _y), levels in masts.items():
        assert x > 0.0
        assert len(levels) >= 2

    # the arc-maximum view can be non-monotonic in height; the mast view is
    # the one that means something
    by_height = b9.arc_maxima_by_height()
    at_57 = {z: c for (x, z), c in by_height.items() if x == 57.0}
    assert at_57[8.0] > at_57[1.0]  # the artefact
    mast_57 = next(v for (x, _y), v in masts.items() if x == 57.0)
    assert mast_57[8.0] > mast_57[1.0]  # here it is real, and worth reporting


@needs_rediphem
@pytest.mark.slow
def test_vertical_profile_is_too_steep():
    """DEGADIS confines the cloud far more tightly than Burro measured.

    At the same mast, 50 to 140 m downwind, the measured concentration at 3 m
    is 0.57 to 1.05 of the value at 1 m; the model gives 0.02 to 0.18. At 8 m
    the model is four orders of magnitude down while the measurements are
    within a factor of a few.

    This is a property of the model, not of the port, and it is why a
    single-height comparison flatters it: with all the contaminant in a thin
    layer the 1 m concentration comes out roughly right while the vertical
    distribution is wrong.
    """
    from degali.run import run_steady
    from degali.validation.rediphem import load
    from degali.validation.trialcase import to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    b9 = trials[("BURRO", "B9")]
    profile, source = run_steady(to_case(b9, averaging_time=18.4).case)
    alpha1 = source.alpha + 1.0
    x_model, sz = profile.rows[:, 0], profile.rows[:, 7]

    checked = 0
    for (x, _y), levels in b9.mast_profiles().items():
        if 1.0 not in levels or 3.0 not in levels:
            continue
        if not (x_model[0] <= x <= 200.0) or levels[1.0] < 1.0:
            continue
        depth = float(np.interp(x, x_model, sz))
        observed = levels[3.0] / levels[1.0]
        modelled = np.exp(-((3.0 / depth) ** alpha1)) / np.exp(
            -((1.0 / depth) ** alpha1)
        )
        assert observed > 0.4
        assert modelled < 0.25
        checked += 1
    assert checked >= 3


@needs_rediphem
@pytest.mark.slow
def test_burro_bias_is_robust_to_the_assumptions():
    """The over-prediction is not an artefact of the deck-building choices.

    Surface temperature, ground exchange and averaging time all move MG by
    less than 0.15 -- much less than the bias itself.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load

    trials = load(REDIPHEM_ROOT)
    variants = {
        "base": {},
        "no exchange": dict(heat_transfer=0, water_transfer=0),
        "long averaging": dict(averaging_time=60.0),
    }
    values = []
    for kw in variants.values():
        options = {"averaging_time": 18.4} | kw
        values.append(
            compare_series(trials, "BURRO", height=1.0, **options).statistics().mg
        )
    assert max(values) - min(values) < 0.2
    assert all(0.7 < v < 1.0 for v in values)


@needs_rediphem
def test_crosswind_sampling_limits_what_can_be_tested():
    """Burro cannot resolve lateral spread; Desert Tortoise can.

    Burro instrumented at most two crosswind positions per arc, and they are a
    symmetric pair. With the plume off-centre -- at 49 m on B7 one reads 15
    mole per cent and its mirror reads 0.09 -- two points cannot constrain a
    width, and they also mean the arc maximum under-samples the true
    centreline. Desert Tortoise has four to seven positions.
    """
    import collections

    from degali.validation.rediphem import load

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}

    def positions(trial):
        kinds = trial.channel_types()
        arcs = collections.defaultdict(set)
        for x, y, _z, kind in trial.sensors().values():
            if kind in kinds and x > 0.0:
                arcs[x].add(round(y, 1))
        return {x: len(ys) for x, ys in arcs.items()}

    assert max(positions(trials[("BURRO", "B9")]).values()) <= 2
    assert max(positions(trials[("TORTOISE", "D3")]).values()) >= 4

    # and the off-centre plume that makes two points useless
    b7 = trials[("BURRO", "B7")].crosswind_profiles(height=1.0, minimum_sensors=2)
    pair = b7[49.0]
    assert max(pair.values()) / (min(pair.values()) + 1e-9) > 50.0


@needs_rediphem
@pytest.mark.slow
def test_lateral_spread_is_too_wide():
    """The complement of the vertical finding.

    Desert Tortoise resolves the crosswind profile well enough to measure a
    half-width. The model's is consistently larger -- by a factor of one to
    three -- so DEGADIS spreads the cloud too far sideways while confining it
    too tightly in the vertical. Both errors move the ground-level centreline
    concentration in opposite directions, which is part of why the
    single-height arc-maximum statistics look better than the cloud does.
    """
    from degali.core.constants import SQRTPI
    from degali.core.steady import PROFILE_COLUMNS
    from degali.run import run_steady
    from degali.validation.rediphem import load
    from degali.validation.trialcase import (
        SMEDIS_SOURCES,
        computed_source,
        jet_to_case,
    )

    i_b = PROFILE_COLUMNS.index("b")
    i_sy = PROFILE_COLUMNS.index("sy")
    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}

    ratios = []
    for name in ("D1", "D2", "D3", "D4"):
        trial = trials[("TORTOISE", name)]
        source = SMEDIS_SOURCES.get(("TORTOISE", name)) or computed_source(trial)
        tc = jet_to_case(trial, source, averaging_time=18.4, reference_height=2.0)
        if not tc.usable:
            continue
        profile, _src = run_steady(tc.case)
        x_model = profile.rows[:, 0]
        for x, ys in trial.crosswind_profiles(
            height=1.0, minimum_sensors=3
        ).items():
            if not (x_model[0] <= x <= x_model[-1]):
                continue
            y = np.array(sorted(ys))
            c = np.array([ys[k] for k in y])
            if c.max() < 0.2:
                continue
            above = y[c >= c.max() / np.e]
            if len(above) < 2:
                continue
            observed = (above.max() - above.min()) / 2.0
            width = float(np.interp(x, x_model, profile.rows[:, i_b])) + (
                SQRTPI / 2.0 * float(np.interp(x, x_model, profile.rows[:, i_sy]))
            )
            ratios.append(width / observed)

    assert len(ratios) >= 4
    assert min(ratios) > 1.0          # always wider than measured
    assert 1.2 < float(np.median(ratios)) < 3.5


# ==========================================================================
# add-on: buoyant lift-off
# ==========================================================================


def test_degadis_has_no_vertical_momentum_at_ground_level():
    """The gap the add-on exists to fill, asserted against the Fortran itself.

    `PSS` and `SSG` use the density difference for exactly two things: lateral
    gravity spreading, guarded so it only acts while the cloud is *denser*
    than air, and the density-concentration slope. There is no z-momentum
    equation, so a cloud that becomes lighter than ambient stops spreading and
    stays on the ground.

    `JETPLU` does have buoyancy -- its vertical momentum balance carries a
    `-RK1*gg*gamma*ccsysz` term -- which is why EX1's plume never touches
    down. The gap is specifically in the ground-level model.
    """
    fortran = REFERENCE_ROOT / "fortran"
    ground = (fortran / "PSS.for").read_text() + (fortran / "SSG.for").read_text()
    # gravity appears only in the common block declaration, never in an equation
    equations = [
        line for line in ground.splitlines()
        if "gg" in line.lower() and "common" not in line.lower()
        and not line.lstrip().startswith(("c", "C", "."))
        and "/parm/" not in line
    ]
    assert not equations, equations

    jet = (fortran / "JETPLU.for").read_text()
    assert "gg*gamma*ccsysz" in jet.replace(" ", "")


def test_liftoff_rises_only_when_buoyant():
    """A buoyant plume leaves the ground; a neutral one does not."""
    from degali.addons import LiftoffPlume

    rhoa = 1.2

    def light(c):
        return 1.0 / (c / 0.5 + (1.0 - c) / rhoa)

    buoyant = LiftoffPlume(rho_ambient=rhoa, wind=2.0, density_of=light).run(
        rate=1.0, concentration=0.5, velocity=2.0, height=0.0
    )
    assert buoyant.lifts_off
    assert buoyant.as_array()[-1, 1] > 10.0

    neutral = LiftoffPlume(
        rho_ambient=rhoa, wind=2.0, density_of=lambda c: rhoa
    ).run(rate=1.0, concentration=0.5, velocity=2.0, height=0.0,
          max_distance=200.0)
    assert not neutral.lifts_off
    assert neutral.as_array()[-1, 1] < 0.5   # no rise at all

    # a plume released already clear of the ground has not lifted off
    already_up = LiftoffPlume(rho_ambient=rhoa, wind=2.0, density_of=light).run(
        rate=1.0, concentration=0.5, velocity=2.0, height=20.0,
        max_distance=100.0,
    )
    assert not already_up.lifts_off
    assert "never in contact" in already_up.reason

    dense = LiftoffPlume(
        rho_ambient=rhoa, wind=2.0,
        density_of=lambda c: 1.0 / (c / 3.0 + (1.0 - c) / rhoa),
    ).run(rate=1.0, concentration=0.5, velocity=2.0, height=0.0,
          max_distance=200.0)
    assert not dense.lifts_off  # sinks instead


def test_liftoff_conserves_contaminant():
    """Contaminant mass flux is a conserved quantity in this model."""
    from degali.addons import LiftoffPlume

    rhoa = 1.2
    r = LiftoffPlume(
        rho_ambient=rhoa, wind=3.0,
        density_of=lambda c: 1.0 / (c / 0.5 + (1.0 - c) / rhoa),
    ).run(rate=2.0, concentration=0.4, velocity=3.0, height=0.0,
          max_distance=300.0)
    a = r.as_array()
    # concentration falls monotonically as air is entrained
    assert np.all(np.diff(a[:, 4]) <= 1e-12)
    assert a[0, 4] == pytest.approx(0.4, rel=1e-6)
    assert a[-1, 4] < 0.05


def test_liftoff_geometry_covers_line_and_round_sources():
    """The lozenge is chosen so one shape spans round, line, and the wide
    shallow cloud a dense phase leaves behind when it turns buoyant."""
    from degali.addons import LiftoffPlume
    from degali.core.constants import PI

    p = LiftoffPlume(rho_ambient=1.2, wind=2.0, segment_length=0.0)
    # clear of the ground, Ls = 0: a circle
    r, ground, free = p._geometry(PI * 4.0, z=10.0, theta=0.0)
    assert r == pytest.approx(2.0, rel=1e-4)
    assert ground == 0.0
    assert free == pytest.approx(2.0 * PI * 2.0, rel=1e-4)

    # centre on the ground: half the circle, and the diameter in contact
    r, ground, free = p._geometry(PI * 4.0 / 2.0, z=0.0, theta=0.0)
    assert r == pytest.approx(2.0, rel=1e-3)
    assert ground == pytest.approx(4.0, rel=1e-3)

    wide = LiftoffPlume(rho_ambient=1.2, wind=2.0, segment_length=50.0)
    r_wide, ground_wide, _f = wide._geometry(PI * 4.0, z=0.0, theta=0.0)
    assert r_wide < 0.3          # same area spread over a long segment
    assert ground_wide > 50.0    # and almost all of it touching the ground


def test_liftoff_richardson_thresholds():
    """Hall and Walker's measured thresholds, as the report reports them."""
    from degali.addons import RI_ONSET, liftoff_regime, richardson_liftoff

    # the report's equivalence: F/(u^3 W) ~ 0.01 corresponds to Ri* ~ 2
    assert richardson_liftoff(0.01 * 8.0 * 1.0, 2.0, 1.0) == pytest.approx(
        RI_ONSET, rel=1e-9
    )
    assert "no lift-off" in liftoff_regime(1.0)
    assert "beginning" in liftoff_regime(5.0)
    assert "10-20" in liftoff_regime(20.0)
    assert "clear of the ground" in liftoff_regime(100.0)


# ==========================================================================
# swappable ground-level closures
# ==========================================================================


def test_default_closure_is_degadis_exactly(b9, b9_downwind):
    """Making the closure pluggable must not change the default path."""
    from degali.core.closures import DegadisClosure

    dw, _h = b9_downwind
    assert isinstance(dw.closure, DegadisClosure)
    assert dw.n_closure == 0
    assert dw.closure.coefficient == pytest.approx(dw.consts.c_spread)


def test_closures_agree_while_the_cloud_is_dense():
    """A buoyant closure must reduce to DEGADIS wherever DEGADIS applies.

    Swapping it in for an LNG case has to give the same slumping, or the
    substitution would be changing answers it has no business changing.
    """
    import numpy as np

    from degali.addons import BuoyantClosure
    from degali.core.closures import CloudState, DegadisClosure

    dense = CloudState(
        dist=50.0, beff=20.0, sz=1.5, heff=1.4, rho=1.35, rhoa=1.0718,
        cc=0.05, temp=280.0, wind=4.0, ustar=0.2, own=np.zeros(2),
        mass_flux=2.0, mass_flux_rate=0.05,
    )
    original = DegadisClosure(1.64, 0.106, 8.0)
    buoyant = BuoyantClosure(1.64, 0.106, 8.0)
    assert buoyant.lateral_rate(dense) == pytest.approx(
        original.lateral_rate(dense), rel=1e-12
    )
    # and a dense cloud is pushed down, not up
    assert buoyant.derivatives(dense)[0] < 0.0


def test_buoyant_closure_lifts_what_degadis_freezes():
    """Hydrogen at its lower flammable limit is lighter than air.

    DEGADIS sets the spreading rate to zero there and has no vertical
    momentum equation, so the cloud stays on the ground for ever. The
    replacement gives it an upward acceleration instead.
    """
    import numpy as np

    from degali.addons import BuoyantClosure
    from degali.core.closures import CloudState, DegadisClosure

    light = CloudState(
        dist=50.0, beff=20.0, sz=1.5, heff=1.4, rho=1.199, rhoa=1.2178,
        cc=0.05, temp=282.0, wind=4.0, ustar=0.2, own=np.zeros(2),
        mass_flux=2.0, mass_flux_rate=0.05,
    )
    assert DegadisClosure(1.64, 0.106, 8.0).lateral_rate(light) == 0.0
    assert BuoyantClosure(1.64, 0.106, 8.0).derivatives(light)[0] > 0.0


def test_rise_is_bounded_by_entrainment():
    """Momentum, not velocity, is the integrated quantity.

    The air a cloud entrains carries no vertical momentum, so it dilutes what
    the cloud has. Integrating dw/dx directly omits that and lets a hydrogen
    cloud accelerate past 20 m/s, an order of magnitude faster than a buoyant
    plume rises.
    """
    import numpy as np

    from degali.addons import BuoyantClosure
    from degali.core.closures import CloudState

    cl = BuoyantClosure(1.64, 0.106, 8.0)
    strong = dict(
        dist=50.0, beff=5.0, sz=1.0, heff=1.0, rho=0.9, rhoa=1.2178,
        cc=0.05, temp=150.0, wind=2.0, ustar=0.15,
    )
    at_rest = cl.derivatives(CloudState(own=np.array([0.0, 0.5]),
                                        mass_flux=2.0, mass_flux_rate=0.5,
                                        **strong))[0]
    moving = cl.derivatives(CloudState(own=np.array([5.0, 0.5]),
                                       mass_flux=2.0, mass_flux_rate=0.5,
                                       **strong))[0]
    assert moving < at_rest      # the dilution term opposes the rise
    assert at_rest > 0.0


def test_closure_hands_over_rather_than_extrapolating():
    """The downwind model's entrainment is calibrated for a ground-hugging
    cloud. Once the cloud has left the ground it entrains through its whole
    perimeter instead, so the closure stops and the plume model takes over."""
    import numpy as np

    from degali.addons import BuoyantClosure
    from degali.core.closures import CloudState

    cl = BuoyantClosure(1.64, 0.106, 8.0, handover_height=2.0)
    assert not np.isfinite(cl.lifted_off)
    high = CloudState(
        dist=120.0, beff=5.0, sz=1.0, heff=1.0, rho=0.9, rhoa=1.2178,
        cc=0.05, temp=150.0, wind=2.0, ustar=0.15,
        own=np.array([1.0, 5.0]), mass_flux=2.0, mass_flux_rate=0.5,
    )
    assert np.all(cl.derivatives(high) == 0.0)
    assert cl.lifted_off == pytest.approx(120.0)


def test_hydrogen_is_buoyant_across_the_whole_flammable_range():
    """Why the closure has to be swappable at all.

    On the adiabatic mixing line for liquid hydrogen into ambient air, the
    cloud is denser than air only above about 99.9 mole per cent. At the
    lower flammable limit of 4 mole per cent it is lighter, and the minimum
    density ratio sits near 72 mole per cent, inside the flammable range.
    Every concentration a flammability assessment cares about is buoyant.
    """
    pytest.importorskip("CoolProp")
    import numpy as np
    from CoolProp.CoolProp import PropsSI

    from degali.core.atmosphere import absolute_humidity
    from degali.core.thermo import (
        AmbientConditions,
        CoolPropBackend,
        GasProperties,
        Thermo,
    )

    backend = CoolPropBackend("Hydrogen")
    tamb = 288.6
    humid, _rh = absolute_humidity(
        tamb, 1.0, backend.water_vapour_pressure, relhum=65.0
    )
    tboil = PropsSI("T", "P", 101325, "Q", 0, "Hydrogen")
    gas = GasProperties(
        name="LH2", mw=2.016, temp=tboil,
        rho=PropsSI("D", "T", tboil, "Q", 1, "Hydrogen"),
        coolprop_name="Hydrogen",
    )
    ambient = AmbientConditions(
        tamb=tamb, pamb=1.0, humid=humid, tsurf=tamb, ihtfl=1, iwtfl=1
    )
    th = Thermo(gas=gas, ambient=ambient, backend=backend, legacy_numerics=False)
    th.reference_enthalpies()
    table = th.build_adiabatic_table(1.0, 0.0, th.hmrte)

    a = table.as_array()
    yc, rho = a[:, 0], a[:, 2]
    rhoa = table.rhoa

    # dense only at the very top of the mixing line
    assert yc[rho > rhoa].min() > 0.99
    # buoyant at both flammable limits and at stoichiometric
    for fraction in (0.04, 0.295, 0.75):
        assert float(np.interp(fraction, yc, rho)) < rhoa
    # and the strongest buoyancy is inside the flammable range
    assert 0.5 < yc[int(np.argmin(rho))] < 0.8


def test_property_grid_adapts_to_the_substance():
    """A grid chosen for LNG puts hydrogen's below absolute zero.

    The default span starts at 80 K, which is fine for methane at 112 K and
    ammonia at 240 K. Hydrogen boils at 20 K, and simply widening by the
    margin gives a lower bound of -30 K: a sizeable part of a 4001-point grid
    spent on temperatures that do not exist, and coarser interpolation over
    the ones that do.
    """
    pytest.importorskip("CoolProp")
    from degali.core.thermo import CoolPropBackend

    assert CoolPropBackend("Methane").SPAN == (80.0, 400.0)
    hydrogen = CoolPropBackend("Hydrogen")
    assert hydrogen.SPAN[0] >= CoolPropBackend.FLOOR
    assert hydrogen.SPAN[0] < 25.0

    # and walking the full range costs almost no rebuilds
    rebuilds = 0
    original = CoolPropBackend._interp

    def counted(self, key, temp, compute, pressure=0.0):
        nonlocal rebuilds
        grid = self._grids.get(key)
        if grid is None or not (grid[0] <= temp <= grid[1]):
            rebuilds += 1
        return original(self, key, temp, compute, pressure)

    CoolPropBackend._interp = counted
    try:
        for t in np.linspace(288.0, 20.4, 300):
            hydrogen.cp_air(float(t))
    finally:
        CoolPropBackend._interp = original
    assert rebuilds <= 2
    assert hydrogen._grids["cpa"][0] >= CoolPropBackend.FLOOR


def test_mixing_table_bound_is_a_fortran_limit_not_a_physical_one():
    """IGEN = 42 is an array dimension. Liquid hydrogen needs more nodes.

    The mixing line spans 20 to 289 K, a factor of fourteen, against LNG's
    factor of under three, so the adaptive thinning keeps more points. The
    bound is raised rather than the tolerance loosened, because loosening it
    would degrade every lookup silently.
    """
    pytest.importorskip("CoolProp")
    from CoolProp.CoolProp import PropsSI

    from degali.core.constants import IGEN
    from degali.core.thermo import (
        AmbientConditions,
        CoolPropBackend,
        GasProperties,
        Thermo,
    )

    tboil = PropsSI("T", "P", 101325, "Q", 0, "Hydrogen")
    th = Thermo(
        gas=GasProperties(name="LH2", mw=2.016, temp=tboil,
                          rho=PropsSI("D", "T", tboil, "Q", 1, "Hydrogen"),
                          coolprop_name="Hydrogen"),
        ambient=AmbientConditions(tamb=288.6, pamb=1.0, humid=0.007,
                                  tsurf=288.6, ihtfl=1, iwtfl=1),
        backend=CoolPropBackend("Hydrogen"), legacy_numerics=False,
    )
    th.reference_enthalpies()
    assert th.build_adiabatic_table(1.0, 0.0, th.hmrte).n > IGEN
    with pytest.raises(RuntimeError, match="overflow"):
        th.build_adiabatic_table(1.0, 0.0, th.hmrte, maxrows=IGEN)


def test_lozenge_profile_shape_factor():
    """The cross-section carries a profile, not a top hat.

    An integral model conserves the contaminant flux, so what it carries is the
    *mean* concentration over the section; a grab bottle sees the peak.
    Comparing one against the other under-predicts by the peak-to-mean ratio.

    For a circular section cut at 2.15 sigma that ratio is 2.57, and the
    discrepancy measured against the NASA sample bottles with a top-hat
    section was 2.45 — which is what identified the shape as the cause rather
    than any coefficient.
    """
    import math

    from degali.addons.liftoff import LiftoffState

    def state(ls, r):
        return LiftoffState(
            s=0.0, x=0.0, z=10.0, theta=0.0, u=2.0, radius=r, area=0.0,
            concentration=0.01, density=1.0, ground_contact=0.0,
            segment_length=ls,
        )

    circular = state(0.0, 10.0)
    assert 1.0 / circular.shape_factor == pytest.approx(2.57, rel=0.02)
    # a long flat section is closer to one-dimensional and less peaked
    line = state(50.0, 10.0)
    assert 1.0 / line.shape_factor < 1.0 / circular.shape_factor

    # the profile is peaked on the axis and vanishes at the boundary
    s = circular
    assert s.at_height(s.z) == pytest.approx(s.peak_concentration)
    assert s.at_height(s.z + 0.99 * s.radius) < 0.2 * s.peak_concentration
    assert s.at_height(s.z + 1.01 * s.radius) == 0.0
    # and it conserves nothing it should not: the peak exceeds the mean
    assert s.peak_concentration > s.concentration


# ==========================================================================
# the liquid hydrogen entry point
# ==========================================================================


@pytest.mark.slow
def test_lh2_assess_reproduces_the_nasa_regimes():
    """One call has to give what the staged path gave.

    The four NASA spills, asked at the 33.8 m tower row: the model must place
    the cloud in the right band for all four, which is the best-established
    output it has.
    """
    pytest.importorskip("CoolProp")
    from degali.lh2 import assess

    trials = [
        # rate kg/s, wind, T, RH, measured lowest flammable height
        (9.23, 1.55, 297.15, 49.0, 18.3, "aloft"),
        (10.29, 3.35, 288.15, 43.0, 6.4, "low"),
        (9.95, 6.30, 285.15, 43.0, 0.3, "grounded"),
        (9.48, 2.20, 288.15, 29.0, 3.4, "low"),
    ]
    errors = []
    for rate, wind, tamb, rh, measured, band in trials:
        r = assess(
            rate=rate, pool_diameter=9.1, wind=wind, ambient_temperature=tamb,
            relative_humidity=rh, max_distance=60.0, at_distance=33.8,
        )
        assert r.regime == band, (rate, wind, r.regime, band)
        errors.append(r.lowest_flammable_height - measured)

    rms = float(np.sqrt(np.mean(np.square(errors))))
    assert rms < 4.0     # against clouds ten to fifteen metres deep
    assert abs(np.mean(errors)) < 2.0   # and no systematic offset


@pytest.mark.slow
def test_lh2_assess_handles_a_pressurised_release():
    """A jet goes through the flashing source, a pool does not."""
    pytest.importorskip("CoolProp")
    from degali.lh2 import assess

    jet = assess(
        rate=0.285, orifice=0.0254, storage_pressure=6.0, wind=2.46,
        height=0.5, ambient_temperature=288.6,
    )
    assert jet.notes["source"] == "flashing jet"
    assert jet.notes["model path"] == "JetPlume"
    assert float(jet.notes["source mass fraction"]) < 1.0
    assert jet.distance_to_lfl > jet.distance_to_stoichiometric > 0.0

    pool = assess(rate=0.285, pool_diameter=2.0, wind=2.46)
    assert pool.notes["source"] == "pool or vapour"
    assert pool.notes["model path"] == "LiftoffPlume"
    assert float(pool.notes["source mass fraction"]) == pytest.approx(1.0)


@pytest.mark.slow
def test_lh2_assess_says_when_it_is_extrapolating():
    """Outside the range it was checked in, it still runs and says so.

    The range is per release type. A jet and a pool have been checked by
    different campaigns over different conditions, and treating them as one
    range reported a jet asked about 30 m as inside a limit set by a pool
    experiment 33.8 m from its pond.
    """
    pytest.importorskip("CoolProp")
    from degali.lh2 import VALIDATED, assess

    inside = assess(rate=9.5, pool_diameter=9.1, wind=2.5, max_distance=30.0)
    assert not any("rate" in w or "wind" in w for w in inside.warnings)

    outside = assess(rate=500.0, pool_diameter=20.0, wind=0.2, max_distance=30.0)
    assert any("rate" in w for w in outside.warnings)
    assert any("wind" in w for w in outside.warnings)
    assert VALIDATED["pool"]["rate"][1] < 500.0

    # the same rate is far outside the jet range and inside the pool one
    assert VALIDATED["jet"]["rate"][1] < 9.5 < VALIDATED["pool"]["rate"][1]

    # a one-point "range" says that rather than reading as a near miss
    edge = assess(rate=9.0, pool_diameter=9.1, wind=2.5, max_distance=30.0)
    assert any("not a series" in w for w in edge.warnings), edge.warnings


def test_lh2_assess_needs_a_source_geometry():
    pytest.importorskip("CoolProp")
    from degali.lh2 import assess

    with pytest.raises(ValueError, match="exactly one"):
        assess(rate=1.0, wind=2.0)
    with pytest.raises(ValueError, match="exactly one"):
        assess(rate=1.0, wind=2.0, height=0.5,
               orifice=0.0254, pool_diameter=1.0)
    with pytest.raises(ValueError, match="positive release elevation"):
        assess(rate=1.0, wind=2.0, orifice=0.0254)
    with pytest.raises(ValueError, match="at_distance"):
        assess(rate=1.0, wind=2.0, pool_diameter=1.0,
               max_distance=10.0, at_distance=20.0)


@pytest.mark.slow
def test_lh2_report_carries_its_own_confidence():
    """Each number is quoted with how far it has been checked, because they
    are not equally well established."""
    pytest.importorskip("CoolProp")
    from degali.lh2 import assess

    pool = assess(
        rate=9.5, pool_diameter=9.1, wind=2.5, max_distance=30.0
    ).report()
    assert "4 of 4" in pool            # the regime
    assert "RMS 2.9 m" in pool         # the height
    assert "no direct pool-concentration validation" in pool
    assert "mixing line" in pool       # the buoyancy crossover

    # A pool result must not borrow grade-B concentration evidence from a jet.
    assert "grade A" in pool and "grade C" in pool and "grade B" not in pool

    jet = assess(
        rate=0.26, orifice=0.0254, storage_pressure=6.0,
        wind=2.5, height=0.5, max_distance=6.0, at_distance=4.0,
    ).report()
    assert "MG 1.047" in jet and "n=62" in jet and "grade B" in jet


def test_unified_closure_has_no_handover_threshold():
    """One formulation from dense ground cloud to risen plume.

    The staged approach needed a handover height — a free parameter standing
    in for physics, at which neither model was valid: the ground-layer
    entrainment is out of calibration above it and the perimeter law does not
    apply below. The contact fraction replaces it with geometry.
    """
    import numpy as np

    from degali.addons import UnifiedClosure
    from degali.core.closures import CloudState

    cl = UnifiedClosure(
        coefficient=1.64, alpha_wind=0.106, z0=8.0, segment_length=9.1
    )
    assert not hasattr(cl, "handover_height")

    # the fraction is continuous and monotone as the cloud lifts
    fractions = [cl.contact_fraction(3.0, z) for z in (-3.0, -1.0, 0.0, 1.5, 2.9, 3.0)]
    assert fractions[0] == 1.0          # entirely below ground
    assert fractions[-1] == 0.0         # just cleared
    assert all(a >= b for a, b in zip(fractions, fractions[1:]))

    def state(z, radius, rho):
        return CloudState(
            dist=10.0, beff=9.1, sz=radius, heff=radius, rho=rho, rhoa=1.218,
            cc=0.05, temp=200.0, wind=2.0, ustar=0.15,
            own=np.array([0.0, z, radius]), mass_flux=5.0, mass_flux_rate=0.3,
        )

    # a dense cloud on the ground slumps and does not rise
    dense = state(0.0, 1.0, 1.30 * 1.218)
    assert cl.lateral_rate(dense) > 0.0
    assert cl.derivatives(dense)[0] <= 0.0

    # a buoyant cloud that has cleared rises and no longer slumps
    light = state(6.0, 1.0, 0.85 * 1.218)
    assert cl.lateral_rate(light) == 0.0
    assert cl.derivatives(light)[0] > 0.0


def test_unified_entrainment_blends_between_the_two_laws():
    """Ground-layer entrainment where the cloud touches, perimeter where it
    does not, in proportion — not one or the other."""
    import numpy as np

    from degali.addons import UnifiedClosure
    from degali.core.closures import CloudState

    cl = UnifiedClosure(
        coefficient=1.64, alpha_wind=0.106, z0=8.0, segment_length=9.1
    )

    def growth(z):
        s = CloudState(
            dist=10.0, beff=9.1, sz=2.0, heff=2.0, rho=1.0, rhoa=1.218,
            cc=0.05, temp=200.0, wind=2.0, ustar=0.15,
            own=np.array([1.0, z, 2.0]), mass_flux=5.0, mass_flux_rate=0.3,
        )
        return cl.derivatives(s)[2]

    on_ground = growth(-2.0)
    airborne = growth(5.0)
    between = growth(1.0)
    # the free-perimeter law entrains faster than the thin ground layer
    assert airborne > on_ground
    assert on_ground < between < airborne


@needs_rediphem
def test_unified_closure_runs_a_whole_transition():
    """A cloud that starts dense and ends buoyant, in one integration."""
    import numpy as np

    from degali.addons import UnifiedClosure
    from degali.core.closures import CloudState

    cl = UnifiedClosure(
        coefficient=1.64, alpha_wind=0.106, z0=8.0, segment_length=9.1
    )
    own = cl.initial(
        CloudState(dist=0.0, beff=9.1, sz=1.0, heff=1.0, rho=1.5, rhoa=1.218,
                   cc=0.1, temp=100.0, wind=2.0, ustar=0.15)
    )
    mass, x, dx = 5.0, 0.0, 0.5
    heights = []
    for ratio in (1.30, 1.05, 0.95, 0.88, 0.92):
        for _ in range(20):
            s = CloudState(
                dist=max(x, 1.0), beff=9.1, sz=own[2], heff=own[2],
                rho=ratio * 1.218, rhoa=1.218, cc=0.05, temp=200.0,
                wind=2.0, ustar=0.15, own=own, mass_flux=mass,
                mass_flux_rate=0.3,
            )
            own = own + cl.derivatives(s) * dx
            own[1] = max(own[1], 0.0)
            mass += 0.3 * dx
            x += dx
        heights.append(float(own[1]))

    assert heights[0] == pytest.approx(0.0, abs=1e-9)  # dense: stays down
    assert heights[-1] > 5.0                            # buoyant: has risen
    assert all(a <= b + 1e-9 for a, b in zip(heights, heights[1:]))
    assert own[2] > 1.0                                 # and has grown


def test_control_is_the_prmt_array_with_names():
    """Named access to PRMT, without changing what it is.

    ``RKGST`` takes one array carrying the integration bounds, step,
    tolerance and stop flag, and beyond them whatever the derivative routine
    wants to stash. Reading ``prmt[9]`` at a call site tells nobody anything,
    so the array is wrapped -- but it is still the same list, so every index
    the Fortran uses works and the numerics are untouched.
    """
    from degali.core.rkgst import Control

    c = Control([0.0, 10.0, 0.05, 0.003, 80.0] + [0.0] * 17)
    assert (c.lower, c.upper, c.step, c.tolerance) == (0.0, 10.0, 0.05, 0.003)
    assert len(c) == 22

    # halting is what PRMT(5) = 1 does, and it is visible both ways
    assert c.stop == 80.0
    c.halt()
    assert c.stop == 1.0 and c[4] == 1.0

    # the caller's own slots are untouched by any of this
    c[9] = 1.64
    assert c[9] == 1.64
    assert "caller slots" in repr(c)

    # and it is a sequence, so the ported routines cannot tell the difference
    assert list(c)[:5] == [0.0, 10.0, 0.05, 0.003, 1.0]


def test_no_magic_stop_index_left_in_the_drivers():
    """Every driver halts by name rather than by index."""
    core = REFERENCE_ROOT.parent / "src" / "degali" / "core"
    for name in ("steady.py", "transient.py", "jetplume.py", "timesort.py",
                 "driver.py"):
        text = (core / name).read_text()
        assert "prmt[4] = 1.0" not in text, name
        if "prmt" in text:
            assert "Control(" in text or "prmt.halt()" in text, name


# ==========================================================================
# PRESLHY E3.5: the raw dataset
# ==========================================================================

# The dataset is not redistributable (DOI 10.35097/1481, KIT open data), so
# it is located from the environment rather than from a path on the machine
# this was written on. E35_ROOT holds the 24 workbooks; the report supplies
# the sensor positions, which the workbooks do not carry.
#
#     export DEGALI_E35_ROOT=/path/to/10.35097-1481/data/dataset
#     export DEGALI_E35_REPORT=/path/to/PRESLHY_D3_6_...pdf
#
# The reduced CSVs (conditions and far-field arcs) are derived products of
# this work and travel with the repository, so they need no variable.
E35_ROOT = Path(os.environ.get("DEGALI_E35_ROOT", "e35"))
E35_REPORT = Path(os.environ.get("DEGALI_E35_REPORT", "e35-report.pdf"))
REDUCED = REFERENCE_ROOT.parent / "reference" / "preslhy"
needs_e35 = pytest.mark.skipif(
    not (E35_ROOT.exists() and E35_REPORT.exists()),
    reason="set DEGALI_E35_ROOT and DEGALI_E35_REPORT to run these",
)
needs_reduced = pytest.mark.skipif(
    not REDUCED.exists(), reason="the reduced PRESLHY tables are not present"
)


@needs_e35
def test_flow_window_is_the_longest_run_not_the_outermost_crossings():
    """One spike must not stretch the window across the whole record.

    On trial 2 the Coriolis meter oscillates between -18.8 and +30.5 g/s.
    Taking the first and last crossing of half the peak gives a window in
    which only six per cent of the samples are flowing, and a mean mass flow
    of *minus* 0.5 g/s -- which then drives the source term.
    """
    from degali.validation.preslhy import read_trial, sensor_positions

    positions = sensor_positions(E35_REPORT)
    trial = read_trial(E35_ROOT / "trial_2_11-09-2019alldata.xlsx", positions)
    assert trial.flow_mean > 0.0
    assert trial.flow_mean > 0.3 * trial.flow_peak
    first, last = trial.window
    assert last - first >= 5


@needs_e35
def test_sensor_serials_match_the_report_table():
    """The workbook columns carry the serials Table A3 gives coordinates for.

    Without that match the raw data is a pile of unlabelled time series: the
    concentrations are in the workbook and the positions are only in the PDF.
    """
    from degali.validation.preslhy import read_trial, sensor_positions

    positions = sensor_positions(E35_REPORT)
    assert len(positions) >= 25
    trial = read_trial(E35_ROOT / "trial_10_13-09-2019alldata.xlsx", positions)
    assert len(trial.readings) >= 20
    # the near-field array: a few metres downwind, low, and narrow
    assert all(0.3 <= r.x <= 6.5 for r in trial.readings)
    assert all(0.0 <= r.z <= 1.0 for r in trial.readings)
    assert {round(abs(r.y), 1) for r in trial.readings} <= {0.0, 1.0}


@needs_e35
def test_the_measurements_are_uncensored():
    """The point of the raw data.

    The far-field Dräger devices range 0 to 4 vol %, and hydrogen's lower
    flammable limit is 4 %, so every reading that matters is at the ceiling.
    The near-field Xensor sensors range 0 to 100 and are not.
    """
    from degali.validation.preslhy import read_trial, sensor_positions

    positions = sensor_positions(E35_REPORT)
    trial = read_trial(E35_ROOT / "trial_10_13-09-2019alldata.xlsx", positions)
    peaks = [r.peak for r in trial.readings]
    assert max(peaks) > 50.0        # far above the Dräger ceiling
    assert max(peaks) <= 100.0      # and over-range readings are dropped


@needs_e35
def test_reader_covers_the_whole_campaign():
    from degali.validation.preslhy import load

    trials = load(E35_ROOT, E35_REPORT)
    assert len(trials) == 24
    assert all(t.flow_mean > 0.0 for t in trials)
    assert sum(len(t.readings) for t in trials) > 400


def test_presets_are_two_distinct_models():
    """The original and the hydrogen configuration are named, not flagged.

    They are different models. A caller chooses one; neither applies silently.
    """
    from degali.presets import DEGADIS_21, LIQUID_HYDROGEN, PRESETS

    assert DEGADIS_21.closure == "degadis"
    assert DEGADIS_21.legacy_numerics and DEGADIS_21.backend == "legacy"
    assert LIQUID_HYDROGEN.closure == "buoyant"
    assert not LIQUID_HYDROGEN.legacy_numerics
    assert set(PRESETS) == {"DEGADIS 2.1", "liquid hydrogen"}

    # the original carries no applicability limits: it is the reference
    assert DEGADIS_21.applies(rate=1e6, wind=0.01) == []
    # the hydrogen one does, and the momentum ratio is among them
    assert "velocity_ratio" in LIQUID_HYDROGEN.applicability.limits
    assert LIQUID_HYDROGEN.applies(velocity_ratio=1.2)

    for preset in (DEGADIS_21, LIQUID_HYDROGEN):
        text = preset.describe()
        assert preset.name in text and "checked" in text


@pytest.mark.slow
def test_lh2_warns_when_the_wind_steers_the_plume():
    """A steady jet model does not describe a wind-steered release.

    On the PRESLHY 1 barg trials the exit velocity is 4 to 13 m/s against a
    1.5 to 4 m/s wind, and nominally identical releases gave arc maxima of
    83 % and 4 %: the plume went wherever the wind pointed. Those trials
    scatter over a factor of seventy; the momentum-dominated ones agree to
    within sixteen per cent between trials.
    """
    pytest.importorskip("CoolProp")
    from degali.lh2 import MOMENTUM_RATIO, assess

    steered = assess(
        rate=0.030, orifice=0.0254, storage_pressure=2.0, wind=3.0, height=0.5
    )
    assert any("wind steers" in w for w in steered.warnings)

    driven = assess(
        rate=0.285, orifice=0.0254, storage_pressure=6.0, wind=2.5, height=0.5
    )
    assert not any("wind steers" in w for w in driven.warnings)
    assert MOMENTUM_RATIO == 10.0


def test_hall_walker_thresholds_cross_check_the_richardson_scale():
    """Two readings of the same wind-tunnel experiments, and they agree.

    Hall and Walker state their thresholds as F/(W u**3): about 0.01 where the
    concentration maximum leaves the ground, about 0.035 where the ground
    value has dropped to 10-20 per cent of the maximum. Fixing the constant
    at 200 from the first maps the second to Ri* = 7, against the 10 taken
    independently from the URAHFREP report's own reduction.

    That is a check rather than a fit: the constant is set by one threshold
    and the other one then has to land somewhere sensible on its own.
    """
    from degali.addons.liftoff import (
        HW_FIRST_RISE,
        HW_LIFTOFF,
        RI_ONSET,
        RI_SUBSTANTIAL,
        hall_walker_parameter,
        richardson_liftoff,
    )

    flux, wind, width = 0.5, 2.0, 9.1
    assert richardson_liftoff(flux, wind, width) == pytest.approx(
        200.0 * hall_walker_parameter(flux, wind, width), rel=1e-12
    )

    # the constant is anchored on the first threshold
    assert RI_ONSET / HW_FIRST_RISE == pytest.approx(200.0)
    # and the second lands within thirty per cent of the independent value
    implied = 200.0 * HW_LIFTOFF
    assert abs(implied - RI_SUBSTANTIAL) / RI_SUBSTANTIAL < 0.35


@needs_e35
@pytest.mark.slow
def test_sensor_heights_are_measured_from_the_release_axis():
    """Table A3's z is above the nozzle, not above the ground.

    The report's Figure 4 is captioned "near-field array configured for
    releases at 1.5 m height": the array was rebuilt for each release height.
    The measurements say the same thing -- at 0.35 m downwind the 0.5 m sensor
    reads 80.8 vol % for a 0.5 m release and 87.9 for a 1.5 m one, and a jet
    cannot fall a metre in thirty-five centimetres.

    Reading it as a height above ground puts the 1.5 m releases' sensors a
    metre below their own axis, and the model then predicts ten to the minus
    nineteen where 43 to 98 vol % was measured.
    """
    import csv

    from degali.validation.preslhy import load

    with open(REDUCED / "conditions.csv") as fh:
        heights = {
            int(row["trial"]): float(row["release_height_m"])
            for row in csv.DictReader(fh)
        }
    trials = load(E35_ROOT, E35_REPORT, heights=heights)
    by_number = {t.number: t for t in trials}

    low = next(t for t in trials if heights.get(t.number) == 0.5 and t.readings)
    high = next(t for t in trials if heights.get(t.number) == 1.5 and t.readings)
    # the axis-relative coordinate is the same for both
    assert min(r.z_axis for r in low.readings) == pytest.approx(
        min(r.z_axis for r in high.readings)
    )
    # the ground-relative one differs by the release height
    assert min(r.z for r in high.readings) - min(r.z for r in low.readings) == (
        pytest.approx(1.0)
    )
    assert by_number  # the mapping covered the campaign


# ==========================================================================
# attribution: the three findings, pinned to line numbers
# ==========================================================================


def test_adiabat_omits_wa_on_one_path_only():
    """`ADIABAT` computes `wa` for four of its five `ifl` values.

    A defect claim needs the source and a line number, not a plausible
    reading. Here it is: `TPROP.for` assigns `wa` on the branches guarded for
    `ifl` of 0, -1, -2 and 2, and on no branch for 1. `SZF` calls it that way
    and passes a name that appears nowhere else in that file.
    """
    import re

    lines = (REFERENCE_ROOT / "fortran" / "TPROP.for").read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if "subroutine adiabat" in l)

    assigned, current = set(), None
    for line in lines[start:]:
        bare = line.replace(" ", "").replace("\t", "")
        guard = re.match(r"^(?:\d+)?if\(ifl\.ne\.(-?\d+)\)goto", bare)
        if guard:
            current = int(guard.group(1))
        elif bare.startswith("wa=") and current is not None:
            assigned.add(current)
        elif bare.startswith("end") and assigned:
            break
    assert assigned == {0, -1, -2, 2}, assigned
    assert 1 not in assigned

    szf = (REFERENCE_ROOT / "fortran" / "SZF.for").read_text()
    assert "adiabat(1,wclay,walay" in szf.replace(" ", "")
    assert szf.count("walay") == 1  # no declaration, no assignment


def test_pss_and_ssg_differ_by_one_letter():
    """The steady and transient routines write the layer temperature into
    different variables, and every other argument to the call matches.

    Both files carry two `adiabat(0, ...)` calls; the layer one is identified
    by its `yclay` argument.
    """
    def layer_call(name):
        for line in (REFERENCE_ROOT / "fortran" / name).read_text().splitlines():
            bare = line.replace(" ", "").replace("\t", "")
            if "adiabat(0," in bare and "yclay" in bare:
                return bare
        raise AssertionError(f"no layer adiabat call in {name}")

    a, b = layer_call("PSS.for"), layer_call("SSG.for")
    assert a != b
    assert a.replace("temlay", "T") == b.replace("temlam", "T")


def test_gaminc_is_unregularised_on_purpose():
    """Not a defect: the comment above the line says what it is doing.

    A parallel reimplementation of SLAB reported six defects in the original
    and, on obtaining the Fortran, found only one was: three were artefacts of
    a JavaScript transcription treated as the original, two were deliberate
    sentinels misread as errors. The same care applies here.
    """
    text = (REFERENCE_ROOT / "fortran" / "INCGAMMA.for").read_text()
    lines = text.splitlines()
    i = next(i for i, line in enumerate(lines) if "gaminc =" in line)
    preceding = " ".join(lines[max(0, i - 4) : i]).lower()
    assert "multiply" in preceding and "gamma(alpha)" in preceding
    assert "gln" in lines[i]  # the log-gamma is added back deliberately


@needs_rediphem
@pytest.mark.slow
def test_epa_1991_degadis_score_is_an_evaluation_height_artefact():
    """EPA's published DEGADIS bias on Burro is reproducible, and explicable.

    EPA-450/4-90-018 Table 5-7 reports an average fractional bias of −1.07 for
    DEGADIS on Burro — the worst of the four models compared, and about a
    factor of three high. Section 4.3 of the same report says measurements
    were taken from 1 m samplers and compared against model predictions *at
    ground level*.

    DEGADIS's vertical profile is too steep, which this work documents
    separately, so concentration piles up at the ground. Evaluating this port
    at ground level against the same 1 m measurements reproduces the published
    bias; evaluating it at the height the sensors were actually at gives an
    unbiased result.
    """
    from degali.run import run_steady
    from degali.validation.rediphem import load
    from degali.validation.statistics import statistics
    from degali.validation.trialcase import to_case

    trials = [t for t in load(REDIPHEM_ROOT) if t.series == "BURRO"]
    bias = {}
    for height in (0.0, 1.0):
        obs, pred = [], []
        for trial in trials:
            try:
                profile, source = run_steady(
                    to_case(trial, averaging_time=18.4).case
                )
            except Exception:
                continue
            alpha1 = source.alpha + 1.0
            x, centre, depth = (
                profile.rows[:, 0], profile.rows[:, 1], profile.rows[:, 7]
            )
            for arc, value in trial.arc_maxima(
                height=1.0, averaging=18.4
            ).items():
                if value <= 0.05 or not (x[0] <= arc <= x[-1]):
                    continue
                sz = float(np.interp(arc, x, depth))
                model = (
                    float(np.interp(arc, x, centre))
                    * np.exp(-((height / sz) ** alpha1)) * 100.0
                )
                if model > 0.05:
                    obs.append(value)
                    pred.append(model)
        bias[height] = statistics(obs, pred).fb

    # at ground level, EPA's number; at the sensor height, unbiased
    assert bias[0.0] < -1.0                       # EPA reported -1.07
    assert abs(bias[1.0]) < 0.4
    assert bias[0.0] < bias[1.0] - 0.8


# ==========================================================================
# SMEDIS: the model-evaluation spreadsheets
# ==========================================================================

# Also not redistributable. One directory holding the batch folders:
#
#     export DEGALI_SMEDIS_ROOT=/path/to/smedis
#
SMEDIS_ROOTS = [
    Path(os.environ.get("DEGALI_SMEDIS_ROOT", "smedis")) / name
    for name in ("batch1", "batch2", "batch3")
]
needs_smedis = pytest.mark.skipif(
    not any(p.exists() for p in SMEDIS_ROOTS),
    reason="set DEGALI_SMEDIS_ROOT to run these",
)


@needs_smedis
def test_smedis_reader_covers_the_exercise():
    from degali.validation.smedis import load

    trials = [t for root in SMEDIS_ROOTS if root.exists() for t in load(root)]
    assert len(trials) >= 25
    assert sum(len(t.sensors) for t in trials) > 1000
    names = {t.dataset for t in trials}
    assert {"FLADIS", "Desert Tortoise", "Thorney Island"} <= names


@needs_smedis
def test_smedis_recovers_fladis():
    """FLADIS is readable here and not in REDIPHEM.

    REDIPHEM stores concentrations as numbered channels whose meaning is per
    series, defined in a file FLADIS does not ship. Guessing the numbering
    selected channels reading 302 and 23.5 at 20 m -- not concentrations --
    and produced a clean-looking result that meant nothing. In these files
    every sensor is a row carrying its own position and value.
    """
    from degali.validation.rediphem import ENV_VAR, default_root
    from degali.validation.smedis import load

    trials = [t for root in SMEDIS_ROOTS if root.exists() for t in load(root)]
    fladis = [t for t in trials if t.dataset == "FLADIS"]
    assert len(fladis) == 3
    for trial in fladis:
        assert trial.substance == "ammonia"
        assert len(trial.sensors) > 40
        assert max(s.x for s in trial.sensors) > 200.0
        assert not trial.suspect_units

    if default_root() is not None:
        from degali.validation.rediphem import load as load_rediphem

        for t in load_rediphem():
            if t.series == "FLADIS":
                assert t.channel_types() == set(), ENV_VAR


@needs_smedis
def test_smedis_records_the_wind_direction_spread():
    """The meander, as a number rather than an assumption.

    Every fixed-array comparison in this package has had to reason about how
    far the plume wandered. FLADIS recorded under two degrees, so a
    point-to-point comparison there is meaningful; the PRESLHY hydrogen
    trials wandered 33 to 49 degrees, and there it is not.
    """
    from degali.validation.smedis import load

    trials = [t for root in SMEDIS_ROOTS if root.exists() for t in load(root)]
    fladis = {t.trial: t for t in trials if t.dataset == "FLADIS"}
    spreads = [
        t.wind_direction_spread for t in fladis.values()
        if not np.isnan(t.wind_direction_spread)
    ]
    assert spreads and max(spreads) < 10.0


@needs_smedis
def test_smedis_coordinates_are_made_downwind_distances():
    """Thorney Island uses site grid coordinates, and says so.

    The release sits at (400, 200) there, so a sensor at x = 100 is 300 m
    *upwind*, not 100 m downwind. FLADIS and Desert Tortoise put the release
    at the origin. Reading all of them the same way would silently mis-place
    an entire series.
    """
    from degali.validation.smedis import load

    trials = {
        (t.dataset, t.trial): t
        for root in SMEDIS_ROOTS if root.exists() for t in load(root)
    }
    ti = trials[("Thorney Island", "Trial 008")]
    assert ti.origin == (400.0, 200.0)
    assert min(s.x for s in ti.sensors) < 0.0   # genuinely upwind sensors

    dt = trials[("Desert Tortoise", "DT1")]
    assert dt.origin == (0.0, 0.0)
    assert min(s.x for s in dt.sensors) > 0.0


@needs_smedis
def test_smedis_flags_concentrations_that_are_not_percentages():
    """Thorney Island peaks at 2060 under a `mean_C(%)` header.

    Something other than a volume percentage is being reported, so it is
    flagged rather than silently used or silently scaled.
    """
    from degali.validation.smedis import load

    trials = {
        (t.dataset, t.trial): t
        for root in SMEDIS_ROOTS if root.exists() for t in load(root)
    }
    assert trials[("Thorney Island", "Trial 021")].suspect_units
    assert trials[("Prairie Grass", "PG8")].suspect_units
    assert not trials[("FLADIS", "Trial009")].suspect_units
    assert not trials[("Desert Tortoise", "DT1")].suspect_units


@needs_e35
@pytest.mark.slow
def test_the_near_field_cannot_test_transient_behaviour():
    """Worth knowing before anyone tries.

    The transient path -- the observers, the time sort, the receptor
    histories -- has been compared against Burro, where the cloud takes a
    hundred seconds to cross four hundred metres. It cannot be tested on the
    PRESLHY near field: the array spans 0.35 to 6 m from a jet leaving at 16
    to 660 m/s, so the travel time is under a second everywhere and the
    sampling interval is 0.3 s. There is nothing to resolve.

    The measured arrival times bear that out -- a median of about zero, which
    is the right answer and an uninformative one -- and their tails are
    contaminated by hydrogen left over from earlier releases, giving arrivals
    of minus a hundred seconds.
    """
    import csv

    from degali.validation.preslhy import load

    with open(REDUCED / "conditions.csv") as fh:
        rows = list(csv.DictReader(fh))
    heights = {int(r["trial"]): float(r["release_height_m"]) for r in rows}
    winds = {int(r["trial"]): float(r["wind_mean_ms"]) for r in rows}

    trials = load(E35_ROOT, E35_REPORT, heights=heights)
    arrivals, advection = [], []
    for trial in trials:
        wind = winds.get(trial.number)
        if not wind:
            continue
        for timing in trial.timings:
            if timing.peak < 1.0:
                continue
            arrivals.append(timing.arrival)
            advection.append(timing.x / wind)

    assert len(arrivals) > 200
    # advection is of order a second: nothing for a transient model to predict
    assert float(np.median(advection)) < 3.0
    # the median arrival is a few seconds, which is the right answer and an
    # uninformative one: any model of this release predicts the same
    assert abs(float(np.median(arrivals))) < 15.0
    # while the spread is dominated by contamination, not by physics
    assert float(np.percentile(arrivals, 10)) < -50.0


@needs_e35
@pytest.mark.slow
def test_the_momentum_criterion_catches_the_vertical_releases_too():
    """The criterion was set from the horizontal trials; it generalises.

    Three of the four non-horizontal releases with usable flow leave at three
    to five times the wind speed and are wind-steered by the same test that
    separated the horizontal ones. They fail against the model exactly as
    that criterion predicts, so their failure is not new information about the
    model. The fourth is a downward jet from half a metre, which reaches the
    ground inside its own zone of flow development.
    """
    import csv
    import math

    from CoolProp.CoolProp import PropsSI

    from degali.lh2 import MOMENTUM_RATIO
    from degali.validation.flashing import equivalent_source
    from degali.validation.preslhy import load

    with open(REDUCED / "conditions.csv") as fh:
        conditions = {r["trial"]: r for r in csv.DictReader(fh)}
    heights = {
        int(k): float(v["release_height_m"]) for k, v in conditions.items()
    }
    trials = {
        t.number: t
        for t in load(E35_ROOT, E35_REPORT, heights=heights)
    }

    ratios = {}
    for number, trial in trials.items():
        record = conditions.get(str(number))
        if not record or record["orientation"] == "horizontal":
            continue
        rate = trial.flow_mean / 1000.0
        if rate <= 0.005:
            continue
        diameter = float(record["orifice_mm"]) / 1000.0
        stored = PropsSI(
            "T", "P", (float(record["tanker_barg"]) + 1.013) * 1e5, "Q", 0,
            "Hydrogen",
        )
        flash = equivalent_source(
            "Hydrogen", storage_temperature=stored,
            ambient_temperature=float(record["T_mean_C"]) + 273.15,
            ambient_pressure=101325.0, molecular_weight=2.016,
        )
        area = math.pi * diameter * diameter / 4.0
        exit_speed = rate / (flash.orifice_density * area)
        ratios[number] = exit_speed / float(record["wind_mean_ms"])

    assert len(ratios) >= 4
    steered = [n for n, r in ratios.items() if r < MOMENTUM_RATIO]
    assert len(steered) >= 3


@needs_e35
def test_the_far_field_array_cannot_be_placed():
    """Why the 686 far-field readings stay unused.

    The far-field stands are fixed and the plume goes wherever the wind
    points, so pairing a stand with a model position needs the wind direction
    through the release. The dataset has one: a domestic weather station
    reporting a 16-point compass string at five-minute intervals.

    A release lasts two to four minutes, so **every trial contains exactly one
    direction sample**, quantised to 22.5 degrees, and some of those are
    "---". At 14 m a 22.5 degree uncertainty is +/- 5.4 m of lateral
    displacement against stands spaced 2.5 m apart: the array cannot be placed
    relative to the plume even in principle.

    This is a property of the measurement programme, not of the model, and it
    is recorded so that the far-field data is not mistaken for something the
    model failed at.
    """
    import math

    import openpyxl

    path = E35_ROOT / "trial_10_13-09-2019alldata.xlsx"
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    weather = list(book["LocalWeather"].iter_rows(values_only=True))
    flow = list(book["Flowmeter"].iter_rows(values_only=True))

    header = [str(h) for h in weather[0]]
    direction = header.index("Loc_WindDirection")
    start, end = str(flow[1][0]), str(flow[-1][0])
    inside = [r for r in weather[1:] if start <= str(r[0]) <= end]

    assert len(inside) == 1                      # one sample per release
    value = str(inside[0][direction])
    assert value in {"N", "NE", "E", "SE", "S", "SW", "W", "NW",
                     "NNE", "ENE", "ESE", "SSE", "SSW", "WSW", "WNW", "NNW",
                     "---"}                       # a compass string, not degrees

    # 22.5 degrees at the far-field radius, against 2.5 m stand spacing
    lateral = 14.0 * math.tan(math.radians(22.5 / 2))
    assert lateral > 2.5


def test_equation_of_state_fallbacks_are_counted_not_silenced():
    """A caller asking for the equation of state must find out if it failed.

    Three lookups in `CoolPropBackend` fell back to a 1989 correlation or to
    `None` when CoolProp raised, and said nothing. A run could therefore use
    the legacy correlation for one property over one temperature range while
    reporting itself as `coolprop`.

    That is the same class of defect as a misplaced bare `raise` — the answer
    changes and the error does not surface — and it has already cost this
    project once, when a two-phase gap was being filled with the nearest valid
    value. The behaviour is unchanged; it is now counted.
    """
    pytest.importorskip("CoolProp")
    from degali.core.thermo import CoolPropBackend

    backend = CoolPropBackend("Hydrogen")
    backend.latent_heat(250.0)
    backend.cp_contaminant_eos(200.0, 1.0)
    assert backend.fallback_report() == ""   # nothing fell back
    assert backend.fallbacks == {}

    # force the failure the counter exists for
    broken = CoolPropBackend("Hydrogen")

    def unavailable(*_a, **_k):
        raise RuntimeError("no such fluid")

    broken._props = unavailable
    broken._grids.clear()
    assert broken.cp_contaminant_eos(200.0, 1.0) is None
    assert broken.rho_contaminant_eos(200.0, 1.0) is None
    assert broken._latent(250.0) > 0.0        # falls back to a constant

    report = broken.fallback_report()
    assert "fell back" in report
    assert set(broken.fallbacks) == {
        "cp of contaminant", "density of contaminant", "latent heat"
    }


@needs_e35
def test_the_far_field_position_is_fitted_rather_than_looked_up():
    """The wind record cannot place the plume, and does not need to.

    The far-field stands sit on two arcs, at 10 and 14 m, spaced 2 and 2.5 m
    apart in the crosswind direction — so adjacent stands are 10 to 12 degrees
    apart as seen from the release. Placing the plume among them needs the
    wind direction to better than that.

    The raw dataset records wind direction as a compass point at 22.5 degree
    resolution, once every five minutes. Most release windows contain a single
    sample. So the resolution is twice the stand spacing and there is no
    within-trial spread at all: the plume cannot be located, and a
    point-to-point comparison there would be measuring which way the wind
    happened to be pointing.

    Five stands across an arc constrain a Gaussian, so the centreline
    concentration and the lateral spread come out of the measurements
    directly, with the plume's position as a nuisance parameter. What is
    compared is the magnitude and the width, which is what the model predicts
    anyway.
    """
    import math
    import re
    from collections import Counter

    import openpyxl

    stands = {
        8: (10.0, 0.0), 7: (10.0, 2.0), 9: (10.0, -2.0),
        6: (10.0, 4.0), 10: (10.0, -4.0),
        3: (14.0, 0.0), 2: (14.0, 2.5), 4: (14.0, -2.5),
        1: (14.0, 5.0), 5: (14.0, -5.0),
    }
    angles = sorted(
        math.degrees(math.asin(y / r)) for r, y in stands.values()
    )
    spacings = [b - a for a, b in zip(angles, angles[1:]) if b - a > 1e-9]
    assert min(spacings) < 12.0        # stands are ~10-12 degrees apart

    # and the wind record is coarser than that, with one sample per release
    path = E35_ROOT / "trial_10_13-09-2019alldata.xlsx"
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    weather = list(book["LocalWeather"].iter_rows(values_only=True))
    column = next(
        i for i, h in enumerate(weather[0]) if h and "WindDirection" in str(h)
    )
    values = [str(r[column]).strip() for r in weather[1:] if r[column]]
    compass = Counter(v for v in values if re.fullmatch(r"[NSEW]{1,3}", v))
    assert compass                       # compass points, not degrees
    assert len(compass) <= 16            # 22.5 degree resolution at best

    flow = list(book["Flowmeter"].iter_rows(values_only=True))
    first, last = str(flow[1][0]), str(flow[-1][0])
    in_window = [
        r for r in weather[1:] if r[0] and first <= str(r[0]) <= last
    ]
    assert len(in_window) <= 3           # no meander estimate within a release

    # but the arc fits recover what the wind record cannot supply
    import csv

    from degali.validation.preslhy import fit_arcs

    with open(REDUCED / "farfield.csv") as fh:
        fits = fit_arcs(list(csv.DictReader(fh)))
    good = [f for f in fits if f.well_constrained]
    assert len(good) >= 15
    assert all(0.2 < f.sigma_y < 40.0 / math.sqrt(2.0) for f in good)
    assert all(abs(f.offset) < 4.0 for f in good)

    # the fitted spread is several times the modelled plume width, and the
    # difference implies a wind-direction wander of order fifteen degrees --
    # derived from the measurements rather than assumed
    widths = [f.sigma_y for f in good if f.radius == 14.0]
    assert widths and 2.0 < float(np.median(widths)) < 6.0


@needs_reduced
@pytest.mark.slow
def test_far_field_extends_the_validated_range_and_reverses_the_bias():
    """Fitting the arcs takes the comparison from 6 m out to 14 m.

    Near the source the model reads about a third high; at 10 to 14 m it reads
    several times low. A single coefficient cannot do both. What does is the
    trajectory: the modelled plume climbs from 1.4 m at 6 m downwind to 3.6 m
    at 14 m, so sensors at 0.5 to 2.5 m are inside it near the source and
    underneath it further out.
    """
    import csv

    from degali.validation.preslhy import fit_arcs

    with open(REDUCED / "farfield.csv") as fh:
        fits = [f for f in fit_arcs(list(csv.DictReader(fh))) if f.well_constrained]

    inner = [f.centreline for f in fits if f.radius == 10.0]
    outer = [f.centreline for f in fits if f.radius == 14.0]
    assert inner and outer
    # the measured cloud is still substantial at 14 m
    assert float(np.median(outer)) > 0.5
    # and thins with distance, as it must
    assert float(np.median(inner)) > float(np.median(outer))


@needs_reduced
def test_non_horizontal_releases_order_as_the_physics_requires():
    """One arc fit each, but the ordering is the informative part.

    A downward jet impinges and stays compact; an upward one spreads and
    dilutes; a horizontal one falls between. That ordering is a grade C
    result -- a direction on one fit per orientation -- and is recorded as
    such rather than as a performance number.
    """
    import csv

    from degali.validation.preslhy import fit_arcs

    with open(REDUCED / "conditions.csv") as fh:
        orientation = {int(r["trial"]): r["orientation"] for r in csv.DictReader(fh)}
    with open(REDUCED / "farfield.csv") as fh:
        fits = [f for f in fit_arcs(list(csv.DictReader(fh))) if f.well_constrained]

    by_orientation = {}
    for fit in fits:
        by_orientation.setdefault(orientation.get(fit.trial, ""), []).append(fit)

    down = by_orientation.get("vertical_down")
    up = by_orientation.get("vertical_up")
    assert down and up
    # downward: narrow and concentrated; upward: wide and dilute
    assert down[0].sigma_y < up[0].sigma_y
    assert down[0].centreline > up[0].centreline


@needs_e35
@pytest.mark.slow
def test_the_measured_plume_does_not_rise():
    """The trajectory, measured directly rather than inferred.

    Four or five heights at each distance constrain a Gaussian, so the plume's
    centre and its vertical spread come out of the measurements. Over six
    metres the measured centre stays within a tenth of a metre of the release
    height; the model lifts it by a metre.

    That gap is about 2.1 times the measured sigma_z, so the predicted
    value at sensor height is cut by exp(-0.5*2.1**2) — a factor of ten — which is
    what the far-field comparison reports as the model reading several times
    low. Near the source the sensors are still inside the plume and the same
    model reads a third high. One defect, two apparent signs.
    """
    import csv

    from degali.validation.preslhy import fit_vertical, load

    with open(REDUCED / "conditions.csv") as fh:
        records = {r["trial"]: r for r in csv.DictReader(fh)}
    heights = {int(k): float(v["release_height_m"]) for k, v in records.items()}
    trials = load(E35_ROOT, E35_REPORT, heights=heights)

    rise, spread = [], []
    for trial in trials:
        record = records.get(str(trial.number))
        if not record or record["orientation"] != "horizontal":
            continue
        release = float(record["release_height_m"])
        for fit in fit_vertical(trial):
            if fit.well_constrained and fit.x >= 3.0:
                rise.append(fit.centre - release)
                spread.append((fit.x, fit.sigma_z))

    assert len(rise) >= 10
    # the plume stays at the release height out to six metres
    assert abs(float(np.median(rise))) < 0.25
    # and its standard deviation grows to about 0.64 m by then
    far = [s for x, s in spread if x >= 5.0]
    if far:
        assert 0.4 < float(np.median(far)) < 1.0


@needs_e35
def test_vertical_fits_separate_two_faults_a_concentration_test_mixes():
    """Why local quantities are worth measuring.

    A concentration comparison gives MG 0.74 near the source and 3.6 in the
    far field -- a signal no single coefficient produces. Fitting the plume
    centre and the spread separately shows why: the spread is about 1.4 times
    narrower than measured everywhere, which raises the centreline value, and
    the trajectory climbs when it should not, which removes the sensors from
    the plume further out. The two pull opposite ways and cancel near the
    source.
    """
    import csv

    from degali.validation.preslhy import fit_vertical, load

    with open(REDUCED / "conditions.csv") as fh:
        records = {r["trial"]: r for r in csv.DictReader(fh)}
    heights = {int(k): float(v["release_height_m"]) for k, v in records.items()}
    trials = load(E35_ROOT, E35_REPORT, heights=heights)

    fits = [
        f for t in trials for f in fit_vertical(t)
        if f.well_constrained and records.get(str(t.number), {}).get(
            "orientation"
        ) == "horizontal"
    ]
    assert len(fits) >= 25
    # the spread grows monotonically with distance, as a jet's must
    near = [f.sigma_z for f in fits if f.x < 1.5]
    mid = [f.sigma_z for f in fits if 1.5 <= f.x < 3.0]
    far = [f.sigma_z for f in fits if f.x >= 3.0]
    assert near and mid and far
    assert float(np.median(near)) < float(np.median(mid)) < float(np.median(far))


def test_preslhy_fit_width_is_a_gaussian_standard_deviation():
    """The fitted width and JETPLU use the same Gaussian convention.

    The old fit used ``exp(-(x/w)**2)`` and stored ``w`` as ``sigma``.  Since
    ``w=sqrt(2)*sigma``, that made a correct model appear 29.3% too narrow.
    """
    from degali.validation.preslhy import _gaussian

    sigma = 0.4
    assert _gaussian(sigma, 1.0, 0.0, sigma) == pytest.approx(
        math.exp(-0.5)
    )
    assert _gaussian(math.sqrt(2.0) * sigma, 1.0, 0.0, sigma) == \
        pytest.approx(math.exp(-1.0))


# ==========================================================================
# code audit: patterns that hide errors
# ==========================================================================


def test_no_bare_raise_outside_an_except_block():
    """A `raise` with no active exception raises `RuntimeError` instead.

    A parallel reimplementation of SLAB found exactly this: a bare `raise`
    that had drifted below the `except` that gave it meaning, so a genuine
    property-lookup failure surfaced as `RuntimeError: No active exception to
    reraise` — the original error, and its message, gone.
    """
    import ast

    package = REFERENCE_ROOT.parent / "src" / "degali"
    offenders = []
    for path in sorted(package.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))

        class Walk(ast.NodeVisitor):
            def __init__(self):
                self.depth = 0

            def visit_ExceptHandler(self, node):
                self.depth += 1
                self.generic_visit(node)
                self.depth -= 1

            def visit_Raise(self, node):
                if node.exc is None and self.depth == 0:
                    offenders.append(f"{path.name}:{node.lineno}")
                self.generic_visit(node)

        Walk().visit(tree)
    assert not offenders, offenders


def test_no_exception_is_swallowed_without_narrowing_it():
    """`except Exception: pass` turns a wrong answer into a plausible one.

    Every place this package catches broadly, it either re-raises with a
    message or falls back for a reason stated at the site. The one that used
    to swallow silently -- the vapour enthalpy lookup -- now distinguishes a
    missing saturation line, which is expected above the critical point, from
    an unknown fluid, which is not.
    """
    package = REFERENCE_ROOT.parent / "src" / "degali"
    offenders = []
    for path in sorted(package.rglob("*.py")):
        lines = path.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines[:-1]):
            if line.strip() in ("except Exception:", "except:"):
                following = lines[i + 1].strip()
                if following == "pass":
                    offenders.append(f"{path.name}:{i + 1}")
    assert not offenders, offenders


def test_unknown_fluid_is_an_error_not_a_number():
    pytest.importorskip("CoolProp")
    from CoolProp.CoolProp import PropsSI

    from degali.validation.flashing import _vapour_enthalpy

    with pytest.raises(ValueError, match="critical temperature"):
        _vapour_enthalpy("NotAFluid", 200.0, 101325.0, PropsSI)

    # and the legitimate cases still work, on both sides of the critical point
    assert _vapour_enthalpy("Hydrogen", 20.0, 101325.0, PropsSI) > 0.0
    assert _vapour_enthalpy("Hydrogen", 300.0, 101325.0, PropsSI) > 0.0
    assert _vapour_enthalpy("Ammonia", 205.0, 90900.0, PropsSI) > 0.0


def test_the_package_carries_no_third_party_data():
    """Nothing redistributable ships, and the tests say so by skipping.

    REDIPHEM, the SMEDIS spreadsheets and the PRESLHY workbooks are all
    someone else's to distribute. The package contains none of them, and
    without them the suite still passes -- the field tests skip.
    """
    package = REFERENCE_ROOT.parent / "src" / "degali"
    for pattern in ("*.xls", "*.xlsx", "*.DBF", "*.dbf", "*.csv"):
        assert not list(package.rglob(pattern)), pattern


def test_no_absolute_paths_in_the_package():
    """A path from the machine it was written on is not a default.

    This covers the tests as well as the package. It did not, and the tests
    had accumulated eleven of them -- ``/home/claude/e35``, ``/home/claude/exp``
    and eight ``/mnt/project`` CSVs. The suite still passed on a clone, because
    those tests skipped; it just skipped for the wrong reason, and nobody else
    could have run them. A guard that exempts the place the paths actually
    are is not a guard.
    """
    root = REFERENCE_ROOT.parent
    offenders = []
    for base in (root / "src" / "degali", root / "tests"):
        for path in sorted(base.rglob("*.py")):
            for i, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1
            ):
                if line.lstrip().startswith("#"):
                    continue  # prose about a path is not a path
                # assembled rather than written out, so this line is not
                # itself a hit
                if any('"' + "/" + p + "/" in line for p in ("home", "mnt")):
                    offenders.append(f"{path.name}:{i}")
    assert not offenders, offenders


def test_the_citation_file_is_complete_enough_to_archive():
    """`authors` is required by CFF 1.2.0, and this file had none.

    Zenodo rejects it, which is the kind of thing found at the moment of
    publication rather than before it. The placeholder is allowed while the
    version is a dev one and blocked once it is not, so this fails at the
    release rather than silently letting a placeholder be minted a DOI.
    """
    yaml = pytest.importorskip("yaml")

    root = REFERENCE_ROOT.parent
    meta = yaml.safe_load((root / "CITATION.cff").read_text())
    for key in ("cff-version", "message", "title", "authors", "version"):
        assert key in meta, key
    assert meta["authors"], "CFF 1.2.0 requires at least one author"

    placeholder = any(
        "REPLACE" in str(v) for a in meta["authors"] for v in a.values()
    )
    if "dev" not in str(meta["version"]):
        assert not placeholder, "fill in the authors before releasing"
        assert meta.get("doi"), "reserve the DOI before archiving"
        assert meta.get("date-released")

    # every reference carries enough to find it again
    for ref in meta.get("references", []):
        assert ref.get("title") and (ref.get("year") or ref.get("doi")), ref


def test_no_retired_figure_survives_in_the_prose():
    """Documentation drifts; this is what stops it drifting silently.

    A review found one quantity reported four different ways across the
    README, the handover, a test assertion and the program's own output, and
    the near-field statistic quoted from two different generations of the
    result. Both were true once. `evidence.RETIRED` lists the superseded
    figures, and any of them reappearing in prose fails here.
    """
    from degali.evidence import RETIRED

    root = REFERENCE_ROOT.parent
    offenders = []
    for path in sorted(root.rglob("*.md")):
        if "reference" in path.parts:
            continue
        # Mention is not use: a sentence that records *why* a figure was
        # retired has to contain it. Those are exempt, and are the only
        # exemption -- which is why they have to say so in words.
        #
        # The check is per *paragraph* rather than per line. Markdown wraps,
        # so the marker word and the retired figure routinely land on
        # different lines of the same sentence, and a line-granular check
        # rejected its own changelog entry.
        excuses = ("earlier version", "an earlier draft", "was stale",
                   "superseded", "had drifted", "against an actual")
        line_no, paragraph = 1, []
        for i, line in enumerate(
            path.read_text(encoding="utf-8").splitlines() + [""], 1
        ):
            if line.strip():
                if not paragraph:
                    line_no = i
                paragraph.append(line)
                continue
            text = "\n".join(paragraph)
            paragraph = []
            if not text or any(e in text.lower() for e in excuses):
                continue
            for retired, instead in RETIRED.items():
                if retired in text:
                    offenders.append(
                        f"{path.name}:{line_no}: '{retired}' -> {instead}"
                    )
    assert not offenders, offenders


def test_the_reported_evidence_is_the_measured_evidence():
    """The numbers printed to a user must be the ones the suite reproduces.

    `evidence.py` exists because those two had come apart: `lh2.py` printed
    `MG 0.74, VG 1.22, n=53` from a hard-coded f-string while every document
    and the preset carried `MG 0.738, VG 1.41, n=69`.
    """
    from degali.evidence import (
        JET_CORRECTIONS,
        LIFTOFF_HEIGHT,
        NEAR_FIELD_CONCENTRATION,
        NEAR_FIELD_CORRECTED,
        NEUTRAL_BUOYANCY,
    )
    from degali.presets import LIQUID_HYDROGEN

    evidence = LIQUID_HYDROGEN.applicability.evidence
    assert "1.121" in evidence and "69 historical arcs" in evidence
    assert "1.047" in evidence and "62 arcs" in evidence
    assert NEAR_FIELD_CONCENTRATION.value == "MG 1.121"
    assert NEAR_FIELD_CONCENTRATION.n == 69
    assert NEAR_FIELD_CORRECTED.value == "MG 1.047"
    assert NEAR_FIELD_CORRECTED.n == 62
    assert JET_CORRECTIONS.value == "sigma_z ratio 0.901 -> 1.033"
    assert JET_CORRECTIONS.n == 23

    # the preset's limits are the jet ones, and they are the evidence module's
    from degali.evidence import RANGE

    assert LIQUID_HYDROGEN.applicability.limits["distance"] == RANGE["jet"]["distance"]

    # grades, so an n=4 result cannot be written as firmly as an n=69 one
    assert LIFTOFF_HEIGHT.grade == "C" and LIFTOFF_HEIGHT.n == 4
    assert NEAR_FIELD_CONCENTRATION.grade == "B"
    assert NEUTRAL_BUOYANCY.grade == "A" and NEUTRAL_BUOYANCY.n is None


@pytest.mark.slow
def test_the_scope_warning_fires_on_the_answer_not_the_integration_limit():
    """The regression this was written for.

    `assess` checked `min(max_distance, ...)` -- the integration limit, which
    defaults to 100 m -- rather than the distance any reported answer relies
    on. Every call taking the default warned that 100 m was out of range while
    every number it returned was inside it, which trains a reader to ignore
    the warnings.
    """
    pytest.importorskip("CoolProp")
    from degali.lh2 import assess

    # a pool asked at 30 m, inside the NASA tower row at 33.8 m: no warning
    pool = assess(rate=9.2, pool_diameter=9.1, wind=1.6, at_distance=30.0)
    assert not any("distance" in w for w in pool.warnings), pool.warnings

    # the same pool asked beyond the tower row: warned
    far = assess(rate=9.2, pool_diameter=9.1, wind=1.6,
                 max_distance=200.0, at_distance=120.0)
    assert any("distance 120" in w for w in far.warnings), far.warnings

    # a jet is held to the jet range, not the pool one: 30 m is far outside
    # the 6 m array even though it is inside the pool's 33.8 m
    jet = assess(rate=0.26, orifice=0.0254, storage_pressure=6.0,
                 wind=2.5, height=0.5, at_distance=30.0)
    assert any("jet range" in w for w in jet.warnings), jet.warnings

    # and the integration limit alone never raises one
    quiet = assess(rate=9.2, pool_diameter=9.1, wind=1.6,
                   max_distance=1000.0, at_distance=10.0)
    assert not any("distance" in w for w in quiet.warnings), quiet.warnings


@needs_rediphem
def test_the_source_route_follows_the_release_type():
    """An instantaneous release must not go through the pool route.

    Thorney Island is a 14 m cylinder of Freon-air released all at once.
    Sending it through the steady pool route produces a source of 3970 kg/s
    over a 2050 m radius -- a number with no physical meaning -- which then
    fails a degeneracy check and was reported as "no usable measurements".

    That is the worst kind of failure: a true-sounding message for a false
    reason. The route now follows `release_type`, and the trial reaches the
    real obstacle instead.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load
    from degali.validation.trialcase import puff_to_case, to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    ti = trials[("TI", "TI08")]
    assert ti.release_type == "puff"

    # the pool route on a puff builds a case that looks fine and is not: the
    # source rate and radius are wrong by orders of magnitude, and the failure
    # only surfaces when it runs
    pooled = to_case(ti, averaging_time=18.4)
    assert pooled.usable            # nothing in the case build objects
    source = pooled.case.source
    assert float(np.max(source.rate)) > 1000.0      # kg/s, for a 14 m cylinder
    assert float(np.max(source.radius)) > 500.0     # m

    # the puff route builds a case that is at least physical
    puffed = puff_to_case(ti, averaging_time=18.4)
    assert puffed.usable
    assert ti.arc_maxima(averaging=18.4)     # and there *are* measurements

    # and the reported reason for TI08 is now the model limit, not a claim
    # about the data. TI21 is skipped earlier, for a genuinely missing
    # stability class.
    result = compare_series(load(REDIPHEM_ROOT), "TI", averaging_time=18.4)
    assert "stability" in result.skipped["TI21"].lower()
    reason = result.skipped["TI08"].lower()
    assert "velocity" in reason or "van ulden" in reason, reason
    assert "no usable measurements" not in reason


@needs_rediphem
def test_a_skip_reason_describes_the_deck_that_was_run():
    """The reason must come from the failing case, not a fresh pool deck.

    `compare` records why a run stopped on the trial case and returns None.
    The series driver used to throw that away and rebuild a *pool* deck to ask
    what was wrong -- so a puff that failed the van Ulden momentum balance was
    reported as "no usable measurements", which is a true-sounding sentence
    for a false reason, and the trial had eighteen arcs of measurements.

    Where the message is still "no usable measurements" it is now accurate:
    the Lathen trials it names have every concentration channel at or upwind
    of the release.
    """
    from degali.validation.compare import compare_series
    from degali.validation.rediphem import load

    trials = load(REDIPHEM_ROOT)
    result = compare_series(trials, "TI", averaging_time=18.4)
    assert result.skipped["TI08"].startswith("model failed")

    lathen = compare_series(trials, "LATHEN", averaging_time=18.4)
    empty = [n for n, why in lathen.skipped.items() if "no usable" in why]
    assert empty
    by_name = {t.name: t for t in trials if t.series == "LATHEN"}
    for name in empty[:3]:
        trial = by_name[name]
        downwind = [
            info for info in trial.sensors().values()
            if info[3] in trial.channel_types() and info[0] > 0.0
        ]
        assert not downwind, name


@needs_rediphem
@pytest.mark.slow
def test_the_van_ulden_limit_is_not_a_resolution_artefact():
    """A parallel reimplementation of SLAB found a failure that looked
    structural and was marginal instability: at the default grid the first
    step killed the signal that triggers source expansion, and a sixteen-fold
    refinement fixed it with no code change.

    The same diagnosis applied here does not reproduce that. The frontal
    velocity iteration succeeds thirty-three times and fails once, so the root
    is lost part-way rather than absent from the start -- which is what
    marginal instability looks like. But refining the integrator's initial
    step by a factor of a hundred, and tightening its error bound tenfold,
    leaves the failure in the same place. The original Fortran stops there
    too.
    """
    import copy

    import degali.core.driver as driver
    from degali.run import run_source
    from degali.validation.rediphem import load
    from degali.validation.trialcase import puff_to_case

    trials = {(t.series, t.name): t for t in load(REDIPHEM_ROOT)}
    case = puff_to_case(trials[("TI", "TI08")], averaging_time=18.4).case

    original = driver.DriverParameters.__init__
    failures = 0
    for step in (0.5, 0.05, 0.005):
        def patched(self, *a, _step=step, **k):
            original(self, *a, **k)
            self.stpin = _step

        driver.DriverParameters.__init__ = patched
        try:
            run_source(copy.deepcopy(case))
        except RuntimeError as exc:
            assert "velocity" in str(exc)
            failures += 1
        except Exception:
            pass
        finally:
            driver.DriverParameters.__init__ = original

    assert failures == 3      # every refinement fails in the same way


# ==========================================================================
# the five adopted jet-path corrections
# ==========================================================================


def test_expanded_source_start_is_self_consistent():
    """The historical expanded-plane candidate remains reconstructible.

    Taking the density from after the flash while keeping the concentration
    and the area at the orifice is the error EPA's 1991 evaluation records the
    SLAB developer objecting to. It was made here before it was noticed:
    pushing the expanded source through the orifice area gives 457 m/s for a
    25.4 mm liquid hydrogen release.
    """
    pytest.importorskip("CoolProp")
    import math

    from CoolProp.CoolProp import PropsSI

    from degali.core.jetplume import JetPlume
    from degali.validation.flashing import equivalent_source

    stored = PropsSI("T", "P", 6.013e5, "Q", 0, "Hydrogen")
    flash = equivalent_source(
        "Hydrogen", storage_temperature=stored, ambient_temperature=288.6,
        ambient_pressure=101325.0, molecular_weight=2.016,
    )
    rate, orifice = 0.26, 0.0254
    density, fraction, diameter = JetPlume.expanded_source_start(
        mass_flow=rate, orifice_diameter=orifice,
        orifice_density=flash.orifice_density,
        expanded_density=flash.density, expanded_fraction=flash.mass_fraction,
        entrainment_momentum_conserving=False,
    )
    # the expanded plane is wider, less dense and less concentrated
    assert diameter > orifice
    assert density < flash.orifice_density
    assert fraction < 1.0

    # and the velocity it implies is physical, not the 457 m/s the orifice
    # area would give for the same state
    area = math.pi * diameter * diameter / 4.0
    velocity = rate / (fraction * density * area)
    assert 80.0 < velocity < 250.0
    wrong = rate / (flash.mass_fraction * flash.density
                    * math.pi * orifice**2 / 4.0)
    assert wrong > 400.0

    # This density scaling is retained only to reproduce the historical
    # candidate.  It does not conserve the total momentum after stationary
    # air is entrained; the following invariant test covers the adopted path.
    u_orifice = rate / (flash.orifice_density * math.pi * orifice**2 / 4.0)
    assert velocity == pytest.approx(
        u_orifice * math.sqrt(flash.orifice_density / density), rel=1e-6
    )
    assert velocity > u_orifice     # lighter, so faster at the same flux


def test_flashing_source_survives_the_first_thermodynamic_lookup():
    """The mixed flash plane must not turn back into nearly pure hydrogen.

    This is the source-table defect that was invisible while the initial
    density, fraction and area were checked only against each other.  The ODE
    immediately looked the concentration up on a different, pure-H2 mixing
    line and deleted almost all of the air already entrained by the flash.
    """
    pytest.importorskip("CoolProp")

    from degali.core.jetplume import J_CC
    from degali.validation.nearfield import hydrogen_jet

    jp, y0 = hydrogen_jet(
        rate=0.833, diameter=0.0254, wind=2.5, height=0.5,
        ambient_temperature=277.15, relative_humidity=90.0,
        storage_pressure_barg=2.53, wind_reference_height=10.0,
        corrections=True,
    )
    state = jp.th.table.from_concentration(float(y0[J_CC]))

    assert state.wc == pytest.approx(0.4000011710, rel=2e-6)
    assert state.rho == pytest.approx(3.050196834, rel=2e-6)
    assert state.yc == pytest.approx(0.9052186854, rel=2e-6)
    assert state.temp == pytest.approx(20.03846730, abs=2e-5)
    assert float(y0[J_CC]) / state.rho == pytest.approx(state.wc, rel=1e-12)


def test_corrected_flash_retains_conserved_axisymmetric_source_plane():
    """The coupled model must receive the exact corrected source plane."""
    pytest.importorskip("CoolProp")

    from degali.validation.nearfield import hydrogen_jet

    rate = 0.099944
    plume, _initial = hydrogen_jet(
        rate=rate, diameter=0.006, wind=2.0, height=0.5,
        ambient_temperature=289.4, relative_humidity=60.0,
        storage_pressure_barg=5.0, wind_reference_height=1.5,
        corrections=True,
    )
    source = plume.axisymmetric_source

    assert source.theta == 0.0
    assert source.x == 0.0
    assert source.y == 0.5
    assert source.temperature == pytest.approx(plume.th.table.t[-1])
    assert source.density == pytest.approx(plume.th.table.rhoe)
    assert source.fuel_mass_flow == pytest.approx(rate, rel=1e-12)
    assert source.mass_flow == pytest.approx(rate / source.mass_fraction)
    assert source.velocity - plume.local_source_wind > 0.0


@pytest.mark.slow
def test_phase_safe_source_reaches_the_jet_table_without_reverting():
    """The optional 68 K boundary must survive the same table handoff."""
    pytest.importorskip("CoolProp")

    from degali.core.jetplume import J_CC
    from degali.validation.nearfield import hydrogen_jet

    kwargs = dict(
        rate=0.099944, diameter=0.006, wind=2.0, height=0.5,
        ambient_temperature=289.4, relative_humidity=60.0,
        storage_pressure_barg=5.0, wind_reference_height=1.5,
        corrections=True,
    )
    ordinary, _ = hydrogen_jet(**kwargs)
    safe, y0 = hydrogen_jet(**kwargs, bulk_air_phase_safe_source=True)
    state = safe.th.table.from_concentration(float(y0[J_CC]))

    assert ordinary.th.table.t[-1] < 21.0
    assert 67.0 < safe.th.table.t[-1] < 70.0
    assert 67.0 < state.temp < 70.0
    assert 0.19 < state.wc < 0.22
    assert float(y0[J_CC]) / state.rho == pytest.approx(state.wc, rel=1e-12)


def test_flash_entrainment_can_conserve_total_momentum():
    """Station 2--3 entrainment cannot manufacture mixture momentum."""
    import math

    from degali.core.jetplume import JetPlume

    rate, diameter = 0.26, 0.0254
    rho_in, rho_out, fraction = 5.30, 3.05, 0.40
    out_rho, out_fraction, out_diameter = JetPlume.expanded_source_start(
        mass_flow=rate, orifice_diameter=diameter,
        orifice_density=rho_in, expanded_density=rho_out,
        expanded_fraction=fraction,
        entrainment_momentum_conserving=True,
    )
    area_in = math.pi * diameter**2 / 4.0
    area_out = math.pi * out_diameter**2 / 4.0
    u_in = rate / (rho_in * area_in)
    u_out = rate / (out_fraction * out_rho * area_out)

    assert u_out == pytest.approx(fraction * u_in, rel=1e-12)
    assert (rate / fraction) * u_out == pytest.approx(rate * u_in, rel=1e-12)
    assert out_rho == rho_out and out_fraction == fraction


def test_flash_gas_branch_is_scoped_to_the_corrected_jet():
    """Avoid liquid EOS states without changing independent pool validation."""
    pytest.importorskip("CoolProp")

    from degali.core.thermo import CoolPropBackend

    historical = CoolPropBackend("Hydrogen")
    jet = CoolPropBackend("Hydrogen", force_contaminant_gas=True)
    assert historical.force_contaminant_gas is False
    assert jet.force_contaminant_gas is True
    assert jet.rho_contaminant_eos(20.04, 1.0) < 2.0


def test_density_scaled_entrainment_has_the_right_sign():
    """Ricou and Spalding: entrainment scales *inversely* with sqrt(rho_jet).

    Panda and Hecht restate it for cryogenic hydrogen -- "scaling inversely by
    the square root of the density of the jet at the nozzle". Applied the
    other way round it moved the modelled vertical spread from 0.70 of the
    measured value to 0.65; applied correctly, to 0.76.
    """
    from degali.core.jetplume import JetCoefficients

    plain = JetCoefficients()
    scaled = JetCoefficients(density_scaled_entrainment=True)

    # a jet lighter than ambient entrains more, not less
    assert scaled.shear_coefficient(0.0, density_ratio=0.5) > plain.alfa1
    # and a dense one entrains less
    assert scaled.shear_coefficient(0.0, density_ratio=4.0) < plain.alfa1
    # off, the coefficient is untouched
    assert plain.shear_coefficient(0.0, density_ratio=0.5) == plain.alfa1


def test_the_corrections_are_all_off_by_default():
    """Nothing changes unless a caller asks. Every parity test depends on it."""
    import inspect

    from degali.core.jetplume import JetCoefficients
    from degali.validation.nearfield import hydrogen_jet

    c = JetCoefficients()
    assert c.alfa1 == 0.057                    # the DATA statement value
    assert not c.density_scaled_entrainment
    assert not c.plume_transition
    assert not c.vertical_shear
    assert not c.ground_layer_entrainment
    assert c.momentum_entrainment_beta == 0.0
    assert c.rise_drag == 0.0
    assert (
        inspect.signature(hydrogen_jet)
        .parameters["bulk_air_phase_safe_source"]
        .default
        is False
    )


def test_the_lh2_preset_records_what_it_changes():
    """A caller reading the preset must see the corrections and their source."""
    from degali.presets import DEGADIS_21, LIQUID_HYDROGEN

    assert DEGADIS_21.jet_corrections == ()
    assert len(LIQUID_HYDROGEN.jet_corrections) == 5
    text = LIQUID_HYDROGEN.describe()
    for source in ("EPA-450/4-90-018", "Ricou", "Papanicolaou & List 1988"):
        assert source in text


def test_primary_source_entrainment_constants_are_not_fischers_values():
    """Keep the primary-source correction from silently reverting.

    Papanicolaou and List (1988) measured 0.0545, 0.0875 and 0.716. The
    frequently attributed 0.0533, 0.0833 and 0.557 are Fischer et al.'s
    earlier proposed values, which the paper quotes only for comparison.
    """
    from degali.core.jetplume import JetCoefficients

    assert JetCoefficients.ALPHA_JET == 0.0545
    assert JetCoefficients.ALPHA_PLUME == 0.0875
    assert JetCoefficients.RI_PLUME == 0.716


@pytest.mark.slow
def test_hydrogen_jet_defaults_to_the_measurement_but_can_rebuild_the_dump():
    """The model advances while the historical reconstruction stays exact."""
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import hydrogen_jet

    kwargs = dict(
        rate=0.189425, diameter=0.0254, wind=2.4667, height=0.5,
        ambient_temperature=15.4333 + 273.15, relative_humidity=58.0,
        storage_pressure_barg=5.0, wind_reference_height=1.5,
        corrections=True,
    )
    measured, _ = hydrogen_jet(**kwargs)
    historical, _ = hydrogen_jet(**kwargs, alfa1=0.0833)
    assert measured.k.alfa1 == 0.0875
    assert historical.k.alfa1 == 0.0833


# ==========================================================================
# the liquid hydrogen jet: reproducing the configuration, not remembering it
# ==========================================================================


@pytest.mark.slow
def test_the_hydrogen_jet_configuration_is_in_the_package():
    """It was not, and that is why none of the LH2 results had a test.

    Every liquid hydrogen jet number in this project -- the near-field
    statistic, the vertical spread ratio, the trajectory correction -- was
    produced by assembling the thermodynamics, the boundary layer, the flash
    and the coefficients by hand in a working session. None of that assembly
    was in the repository, so none of the results was reproducible from it,
    and an audit found that **not one statistical claim in this package was
    computed by the suite**: they were typed into documents and into an
    f-string.

    `validation.nearfield.hydrogen_jet` is that configuration. This test only
    checks it runs and behaves like a jet; the numbers it produces are pinned
    separately, and where they disagree with what was reported that
    disagreement is itself pinned.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import Trajectory, hydrogen_jet

    jp, y0 = hydrogen_jet(
        rate=0.2853, diameter=0.0254, wind=2.467, height=0.5,
        ambient_temperature=288.6, relative_humidity=58.0,
        storage_pressure_barg=5.0, corrections=False,
    )
    traj = Trajectory(jp.th.table, jp.run(y0, distmx=40.0).rows)
    assert traj.ok

    near, far = traj.at(0.5), traj.at(6.0)
    assert near is not None and far is not None
    # a jet: it spreads, it dilutes, and being buoyant it rises
    assert far.sz > near.sz and far.sy > near.sy
    assert far.cc < near.cc
    assert far.z > near.z
    # and the centreline is near-pure at the nozzle
    assert traj.concentration_at(0.35, 0.0, 0.5) > 90.0
    # the ground image is included, so a sensor on the ground is not zero
    assert traj.concentration_at(6.0, 0.0, 0.0) > 0.0
    assert traj.temperature_at(0.35, 0.0, 0.5) < 100.0
    assert traj.temperature_at(6.0, 0.0, 0.0) < jp.th.ambient.tamb


@pytest.mark.slow
def test_the_corrections_widen_the_section_towards_the_measurement():
    """The vertical spread deficit closes, and this now runs rather than
    being recalled.

    Measured `sigma_z` is 0.34/sqrt(2) m at 2 m downwind: 0.34 m was the
    e-folding width formerly mislabeled as sigma. The corrected single trial
    reaches about 0.91 of the measured standard deviation.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import Trajectory, hydrogen_jet

    def sigma_z(corrections):
        jp, y0 = hydrogen_jet(
            rate=0.2853, diameter=0.0254, wind=2.467, height=0.5,
            ambient_temperature=288.6, relative_humidity=58.0,
            storage_pressure_barg=5.0, corrections=corrections,
        )
        traj = Trajectory(jp.th.table, jp.run(y0, distmx=40.0).rows)
        return traj.at(2.0).sz

    shipped, corrected = sigma_z(False), sigma_z(True)
    measured = 0.34 / math.sqrt(2.0)

    assert shipped < measured          # too narrow as shipped
    assert corrected > shipped         # the corrections widen it
    assert abs(corrected - measured) < abs(shipped - measured)   # towards it
    assert corrected / measured == pytest.approx(0.912, abs=0.03)


# `test_the_reported_trajectory_correction_does_not_reproduce` lived here.
#
# It pinned a bound -- the as-shipped rise stays under 1.0 m across every
# candidate setting -- from when the configuration was still being
# reconstructed and the published 1.07 m looked unreachable. It was written to
# be deleted only when the disagreement was resolved rather than forgotten,
# and it has been: with the validated configuration the figure reproduces at
# 1.00 m over the band it was quoted for.
#
# `test_the_published_trajectory_pair_compares_two_aggregations` replaces it
# and carries the finding that survived, which is not that the number was
# wrong but that it was paired with one taken over a different subset.


def test_the_arc_maximum_is_a_tail_value_and_the_tail_is_sensitive():
    """Why two runs that agree on the centreline disagree by a factor of sixty
    at a sensor, and why the near-field statistic is so sensitive to the
    configuration.

    Rebuilding the saved prediction table gives centreline concentrations
    within 10 to 20 per cent of it at every arc. Off the centreline the two
    diverge by up to a factor of a hundred, and the mechanism is here: at
    0.79 m downwind the two runs differ in `sigma_z` by twenty per cent, and
    the reading at a sensor two and a half sigma out differs by a factor of
    four.

    That matters because **the arc maximum becomes a tail value** as soon as
    the plume centre climbs above the array, which it has by 6 m. The
    statistic is then carried by exactly the quantity the two runs disagree
    about most, which is why a modest difference in the vertical spread moves
    it so far.
    """
    pytest.importorskip("CoolProp")
    import math

    from degali.validation.nearfield import Trajectory, hydrogen_jet

    jp, y0 = hydrogen_jet(
        rate=0.2853, diameter=0.0254, wind=2.467, height=0.5,
        ambient_temperature=288.6, relative_humidity=58.0,
        storage_pressure_barg=5.0, corrections=False,
    )
    traj = Trajectory(jp.th.table, jp.run(y0, distmx=40.0).rows)
    state, table = traj.at(0.79), jp.th.table

    def reading(sigma, z):
        direct = math.exp(-0.5 * ((z - state.z) / sigma) ** 2)
        image = math.exp(-0.5 * ((z + state.z) / sigma) ** 2)
        return 100.0 * table.from_concentration(state.cc * (direct + image)).yc

    # two spreads twenty per cent apart, the size of difference that
    # separated candidate configurations while this was being reconstructed
    narrow, wide = state.sz, state.sz * 1.2

    # on the centreline they are indistinguishable
    assert reading(narrow, state.z) == pytest.approx(
        reading(wide, state.z), rel=1e-9
    )
    # two and a half sigma out they differ by a factor of several
    assert reading(wide, 0.2) / reading(narrow, 0.2) > 3.0

    # and by 6 m the plume centre is above the topmost sensor, so the arc
    # maximum is a tail value rather than a centreline one
    assert traj.at(6.0).z > 1.0


# ==========================================================================
# the LH2 statistics, computed instead of quoted
# ==========================================================================

E35_REDUCED = REFERENCE_ROOT / "preslhy" / "e35_reduced.json"
needs_e35_reduced = pytest.mark.skipif(
    not E35_REDUCED.exists(), reason="the reduced PRESLHY table is not present"
)


@pytest.mark.slow
@needs_e35_reduced
def test_the_exclusion_chain_reproduces_the_published_provenance():
    """24 workbooks to 9 trials to 69 arcs, by criteria fixed in advance.

    `docs/lh2-results.md` publishes this chain so a reader can trace the
    statistic back to the dataset. Until now nothing checked that the code
    produced it. It does, exactly: seven releases are not horizontal, eight
    more leave the nozzle at under ten times the wind speed, and the nine that
    remain give 69 arcs.

    The wind-steered exclusion is the one that matters, and it is a statement
    about where a steady jet model applies rather than a filter on the result:
    the criterion is the exit-velocity ratio, computed from the flash and the
    measured flow, never from agreement.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import from_reduced

    field = from_reduced(E35_REDUCED, corrections=False)

    not_horizontal = [k for k, v in field.excluded.items() if "horizontal" in v]
    steered = [k for k, v in field.excluded.items() if "wind-steered" in v]
    assert len(not_horizontal) == 7
    assert len(steered) == 8
    assert len(field.trials) == 9
    assert len(field.pairs) == 69


@needs_e35_reduced
@pytest.mark.slow
def test_the_near_field_statistic_recomputed():
    """Recompute the headline LH2 result downstream of source establishment.

    The mass-, enthalpy- and momentum-consistent source gives sensor-based
    MG/VG/FAC2 of 1.047/1.425/0.84.  All three pass the usual screening
    limits. Seven 0.35/0.53 m arcs are upstream of its established plane, so
    this comparison starts at the next instrumented arc, 0.79 m.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import common_arcs, from_reduced

    shipped, corrected = common_arcs(
        from_reduced(E35_REDUCED, corrections=False),
        from_reduced(E35_REDUCED, corrections=True),
    )
    assert len(shipped.pairs) == len(corrected.pairs) == 62

    a, b = shipped.statistics(), corrected.statistics()

    # At-sensor, which is the convention Burro was already using and which
    # the LH2 near field is now brought into line with.
    assert a.mg == pytest.approx(1.132, abs=0.03)
    assert b.mg == pytest.approx(1.047, abs=0.03)
    lo, hi = shipped.interval()
    assert lo < 1.0 < hi

    # Under the centreline convention the baseline reproduces the published
    # figures to two per cent -- MG 0.738, CI [0.591, 0.918], VG 1.41,
    # FAC2 0.83 -- which is how the configuration was validated. Kept as a
    # check on the rebuild, not as the reported statistic.
    centreline = common_arcs(
        from_reduced(E35_REDUCED, corrections=False, convention="centreline"),
        from_reduced(E35_REDUCED, corrections=True, convention="centreline"),
    )[0]
    c = centreline.statistics()
    assert c.mg == pytest.approx(0.722, abs=0.02)
    assert c.vg == pytest.approx(1.44, abs=0.05)
    assert c.fac2 == pytest.approx(0.82, abs=0.02)
    assert centreline.interval() == pytest.approx((0.574, 0.903), abs=0.03)

    # and the VG the convention costs is the tail sensitivity surfacing
    assert a.vg > 10.0 * c.vg


@needs_e35_reduced
@pytest.mark.slow
def test_the_comparison_convention_moves_the_headline_number():
    """Which convention is used is a decision, not a detail.

    Against the model centreline the as-shipped MG is 0.855; against the model
    evaluated at each sensor and maximised over the arc -- the same operation
    the measurement gets -- it is 1.077. A quarter of the answer.

    The centreline convention also gives a far better VG, and that is not the
    model being better: it is the comparison being smoother. Once the plume
    centre climbs above the array the arc maximum is a tail value, and tails
    are where the two configurations disagree most, so the sensor convention
    inherits that scatter and the centreline convention hides it.

    Neither is wrong. Mixing them is, and mixing them is what produced the
    0.785-against-1.228 puzzle in the figures.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import from_reduced

    centreline = from_reduced(
        E35_REDUCED, corrections=False, convention="centreline"
    ).statistics()
    sensor = from_reduced(
        E35_REDUCED, corrections=False, convention="sensor"
    ).statistics()

    assert centreline.mg == pytest.approx(0.729, abs=0.02)
    assert sensor.mg == pytest.approx(1.132, abs=0.03)
    # the centreline convention is smoother, not better
    assert centreline.vg < 1.5 < sensor.vg


@needs_e35_reduced
@pytest.mark.slow
def test_the_corrected_run_reproduces_the_reference_trajectory():
    """The acceptance test for the whole rebuild.

    The liquid hydrogen jet configuration was never committed, so it was
    reconstructed, and for a while the reconstruction disagreed with the
    published results badly enough that nothing could be trusted. The session
    holding the archives supplied a full parameter dump for one run -- trial
    10, corrections on -- and this checks the rebuild against it.

    Four things had to be corrected to get here, and each was worth a factor:

    ``distmx`` **is the integration step, not the limit.** The rebuild passed
    40.0 there and 40.0 is the *arclength* limit, ``smax``. Integrating the
    whole plume in one stride produced a trajectory wrong by a factor of two
    that still looked entirely plausible -- monotone, smooth, right order of
    magnitude. Nothing about the output said it was under-resolved.

    ``rhoe`` **stays the saturated-vapour density** even with the expanded
    source on. Only ``rho_exit`` at the starting plane takes the expanded
    value: the expansion changes the state at the plane, not the property of
    the substance.

    The averaging time is **60 s**, not 18.4, and the roughness **0.001 m**,
    not 0.01. Neither was written down; both were recovered by matching
    ``deltay`` and ``ustar`` against the dump.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import (
        REACH, STEP, Trajectory, hydrogen_jet,
    )

    jp, y0 = hydrogen_jet(
        rate=0.189425, diameter=0.0254, wind=2.4667, height=0.5,
        ambient_temperature=15.4333 + 273.15, relative_humidity=58.0,
        storage_pressure_barg=5.0, wind_reference_height=1.5,
        corrections=True, alfa1=0.0833, source_table_consistency=False,
        source_momentum_consistency=False,
    )

    # the boundary layer and the plume constants, against the dump
    assert jp.ustar == pytest.approx(0.11804, abs=1e-5)
    assert jp.rhoa == pytest.approx(1.21839, abs=1e-5)
    assert jp.rhoe == pytest.approx(1.33217, abs=1e-5)
    assert jp.deltay == pytest.approx(0.08581, abs=1e-5)
    assert jp.deltaz == pytest.approx(0.04134, abs=1e-5)
    assert jp.betaz == pytest.approx(1.17370, abs=1e-5)
    assert jp.gammaz == pytest.approx(-0.031600, abs=1e-6)
    assert jp.yclow == 1.0e-5

    # the starting state
    assert float(y0[0]) == pytest.approx(1.124108, rel=1e-5)
    assert float(y0[5]) == pytest.approx(0.5, abs=1e-9)

    traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)

    # and the trajectory itself, to the precision the dump was given to
    for x, z, sigma_z, sigma_y in (
        (1.01, 0.502, 0.1553, 0.1728),
        (2.09, 0.542, 0.3019, 0.3311),
        (4.08, 0.787, 0.5027, 0.5515),
        (6.04, 1.158, 0.6909, 0.7549),
    ):
        state = traj.at(x)
        assert state.z == pytest.approx(z, abs=0.003), x
        assert state.sz == pytest.approx(sigma_z, abs=0.002), x
        assert state.sy == pytest.approx(sigma_y, abs=0.002), x


@needs_e35_reduced
@pytest.mark.slow
def test_both_configurations_match_their_reference_dumps():
    """Trial 10, as shipped and corrected, against parameter dumps of both.

    The baseline source plane is the one thing that could not be inferred and
    it is worth stating, because it is what the expanded-source correction
    exists to fix: **the density is the two-phase value after the flash, while
    the concentration and the area are the orifice's.** One plane for one
    quantity and a different plane for the others. The correction makes the
    three consistent.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import (
        REACH, STEP, Trajectory, hydrogen_jet,
    )

    def trajectory(corrections):
        jp, y0 = hydrogen_jet(
            rate=0.189425, diameter=0.0254, wind=2.4667, height=0.5,
            ambient_temperature=15.4333 + 273.15, relative_humidity=58.0,
            storage_pressure_barg=5.0, wind_reference_height=1.5,
            corrections=corrections, alfa1=0.0833,
            source_table_consistency=False, source_momentum_consistency=False,
        )
        return y0, Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)

    # as shipped: orifice density 5.34689, concentration 1.0, orifice diameter
    y0, shipped = trajectory(False)
    assert float(y0[0]) == pytest.approx(5.346894, rel=1e-6)
    assert float(y0[5]) == pytest.approx(0.5, abs=1e-9)
    for x, z, sigma_z in ((1.01, 0.502, 0.0855), (2.04, 0.534, 0.1702),
                          (4.02, 0.847, 0.3479), (6.05, 1.377, 0.5542)):
        state = shipped.at(x)
        assert state.z == pytest.approx(z, abs=0.003), x
        assert state.sz == pytest.approx(sigma_z, abs=0.002), x

    # corrected: the expanded plane, 2.54663 and a 0.04602 m equivalent bore
    y0, corrected = trajectory(True)
    assert float(y0[0]) == pytest.approx(1.124108, rel=1e-5)
    for x, z, sigma_z in ((1.01, 0.502, 0.1553), (2.09, 0.542, 0.3019),
                          (4.08, 0.787, 0.5027), (6.04, 1.158, 0.6909)):
        state = corrected.at(x)
        assert state.z == pytest.approx(z, abs=0.003), x
        assert state.sz == pytest.approx(sigma_z, abs=0.002), x


@needs_e35_reduced
@pytest.mark.slow
def test_the_corrections_improve_the_trajectory_by_a_quarter():
    """`1.07 -> 0.19 m` is withdrawn at both ends.

    Neither figure comes from a configuration that can now be reconstructed,
    and they were taken over different subsets besides -- 5-7 m over all
    release heights against 3 m and beyond for the 0.5 m releases alone, four
    fits against seven.

    From the two dumped configurations, on trial 10 at 6 m, the pair is
    **0.877 -> 0.658 m**: the corrections take about a quarter off the rise,
    not four fifths. Over the 5-7 m band across trials it is 1.00 -> 0.85 m,
    the same fraction.

    The fault is untouched by this and is the largest open defect in the
    liquid hydrogen path. The measurement puts the plume *below* the nozzle at
    that distance and the model puts it a metre above; a quarter off a metre
    is not a fix.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import (
        REACH, STEP, Trajectory, hydrogen_jet, vertical,
    )

    def rise_at_six(corrections):
        jp, y0 = hydrogen_jet(
            rate=0.189425, diameter=0.0254, wind=2.4667, height=0.5,
            ambient_temperature=15.4333 + 273.15, relative_humidity=58.0,
            storage_pressure_barg=5.0, wind_reference_height=1.5,
            corrections=corrections, alfa1=0.0833,
            source_table_consistency=False, source_momentum_consistency=False,
        )
        traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)
        return traj.at(6.05).z - 0.5

    a, b = rise_at_six(False), rise_at_six(True)
    assert a == pytest.approx(0.877, abs=0.005)
    assert b == pytest.approx(0.658, abs=0.006)
    improvement = (a - b) / a
    assert 0.15 < improvement < 0.35, improvement

    # the same fraction over the band, across trials
    shipped = vertical(
        E35_REDUCED, corrections=False, source_table_consistency=False,
        source_momentum_consistency=False,
    )
    corrected = {
        (r["trial"], r["x"]): r
        for r in vertical(
            E35_REDUCED, corrections=True, source_table_consistency=False,
            source_momentum_consistency=False,
        )
    }
    far = [r for r in shipped if 5.0 <= r["x"] <= 7.0]
    med = lambda key: float(np.median([key(r) for r in far]))
    band_a = med(lambda r: r["modelled_centre"] - r["release_height"])
    band_b = med(
        lambda r: corrected[(r["trial"], r["x"])]["modelled_centre"]
        - r["release_height"]
    )
    measured = med(lambda r: r["measured_centre"] - r["release_height"])
    assert band_a == pytest.approx(1.00, abs=0.06)
    assert band_b == pytest.approx(0.85, abs=0.06)
    assert 0.10 < (band_a - band_b) / band_a < 0.30

    # and the fault: the measurement puts the plume below the nozzle here
    assert measured < 0.0
    assert band_b - measured > 0.9


@needs_e35_reduced
@pytest.mark.slow
def test_the_vertical_spread_correction_reproduces():
    """The Gaussian convention correction removes a false width deficit.

    Three things had to be identified, none of them written down:

    **The basis is 23 fits, not 42.** The concentration statistic's momentum
    filter *is* applied to the spread ratio, even though it is not applied to
    the fits used for the trajectory. Two different populations in the same
    table.

    **The statistic is the mean of per-fit ratios**, not the ratio of medians.
    The two remain separately asserted.

    **The measurement is the fit, the model is evaluated at the fit's own
    distance.** Both configurations are validated against parameter dumps to
    four significant figures, so nothing here rests on a reconstruction.

    The former fit used ``exp(-(z/w)**2)`` and called ``w`` sigma, while
    JETPLU uses ``exp(-0.5*(z/sigma_z)**2)``.  Since ``w=sqrt(2)*sigma_z``,
    the old 0.64 and 0.731 ratios were low by exactly ``sqrt(2)``.  With one
    standard-deviation convention the historical model is 0.901 and the
    mass-, enthalpy- and momentum-consistent source is 1.033.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import vertical

    shipped = vertical(E35_REDUCED, corrections=False, momentum_filter=True)
    corrected = {
        (r["trial"], r["x"]): r
        for r in vertical(E35_REDUCED, corrections=True, momentum_filter=True)
    }
    assert len(shipped) == 23, "the spread ratio uses the filtered basis"

    def ratio(table=None):
        return float(np.mean([
            (r if table is None else table[(r["trial"], r["x"])])[
                "modelled_sigma_z"
            ] / r["measured_sigma_z"]
            for r in shipped
        ]))

    a, b = ratio(), ratio(corrected)
    assert a == pytest.approx(0.901, abs=0.02)
    assert b == pytest.approx(1.033, abs=0.02)
    assert a < b
    assert abs(b - 1.0) < abs(a - 1.0)
    assert b < 1.08

    # the form of the statistic is worth this much, so it gets asserted too
    medians = float(np.median([r["modelled_sigma_z"] for r in shipped])) / (
        float(np.median([r["measured_sigma_z"] for r in shipped]))
    )
    assert abs(medians - a) > 0.08, "ratio-of-medians is a different number"


@needs_e35_reduced
def test_the_published_fit_counts_identify_the_basis():
    """42 and 18 are what pinned down which fits the figures used.

    Neither count was written down as a *criterion*, only as a number, and the
    two plausible readings -- every horizontal trial, or only the nine
    momentum-driven ones -- give different answers for the vertical spread.
    The counts settle it.
    """
    import json

    data = json.loads(E35_REDUCED.read_text())
    assert data["gaussian_width_definition"].startswith("standard deviation")

    horizontal = [
        f
        for t in data["trials"] if t["orientation"] == "horizontal"
        for f in t["vertical_fits"] if f["well_constrained"]
    ]
    assert len(horizontal) == 42

    arcs = data["arc_fits"]
    assert len(arcs) == 65
    assert len([f for f in arcs if f["well_constrained"]]) == 18


# `test_the_rediphem_reduction_is_marked_incomplete` lived here. It held a
# flag on a reduction that carried 840 arc maxima and none of the conditions
# a prediction needs, so that a file which looked usable and was not could not
# quietly cost the next reader a session. The complete reduction arrived and
# the incomplete one is gone; `test_the_complete_reduction_builds_runnable_cases`
# is what replaced it.


def test_a_reduced_rediphem_trial_serves_the_readers_unchanged():
    """A reduction is a different transport, not a different data model.

    `from_reduced` rebuilds the specification as the same `specs` dict the
    archive parser produces, so every existing accessor works: `to_case` does
    not need to know where the trial came from. Round-tripping a synthetic
    trial checks that, without needing the archive.

    The one thing a reduction cannot do is re-average: it carries maxima at
    the averaging time recorded in the file, so asking for a mast profile or a
    time series raises rather than quietly returning something else.
    """
    import json
    import tempfile

    from degali.validation.rediphem import from_reduced

    record = {
        "source": "synthetic", "note": "round-trip check",
        "averaging_time_s": 18.4, "complete": True,
        "trials": [{
            "series": "BURRO", "name": "B9", "substance": "LNG",
            "release_type": "pool", "stability": 5,
            "conditions": {
                "pool diameter": 22.0, "release rate": 130.0,
                "release duration": 79.0, "site average windspeed": 5.7,
                "reference height for wind": 1.0, "surface roughness": 0.0002,
                "ambient temperature": 308.6, "ambient pressure": 0.94,
                "relative humidity": 7.2,
            },
            "qualifiers": {"exit temperature": "est"},
            "arcs": {"z=1.0": {"57.0": 12.5, "140.0": 4.1, "400.0": 1.2}},
        }],
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(record, fh)
        path = fh.name

    trial, = from_reduced(path)
    assert trial.series == "BURRO" and trial.name == "B9"
    assert trial.substance == "LNG"
    assert trial.release_type == "pool"
    assert trial.rate == pytest.approx(130.0)
    assert trial.duration == pytest.approx(79.0)
    assert trial.wind_height == pytest.approx(1.0)
    assert trial.stability == 5
    assert trial.value("pool diameter") == pytest.approx(22.0)
    assert trial.qualifiers["exit temperature"] == "est"

    arcs = trial.arc_maxima(height=1.0)
    assert arcs == {57.0: 12.5, 140.0: 4.1, 400.0: 1.2}
    assert trial.arc_maxima(height=8.0) == {}

    # what a reduction cannot answer, it refuses rather than approximates
    with pytest.raises(ValueError, match="give one"):
        trial.arc_maxima()
    with pytest.raises(NotImplementedError):
        trial.time_series()


REDIPHEM_FULL = REFERENCE_ROOT / "rediphem_reduce_full.json"
needs_rediphem_full = pytest.mark.skipif(
    not REDIPHEM_FULL.exists(),
    reason="the complete REDIPHEM reduction is not present",
)


@needs_rediphem_full
def test_the_complete_reduction_builds_runnable_cases():
    """313 trials with every `to_case` input, read back as cases.

    The first reduction carried 840 arc maxima and no conditions, so no
    statistic could be computed from it. This one carries the built `Case`
    rather than the specification, which is the right thing to store: the
    specification needs `to_case` to interpret it, and `to_case` makes
    judgements -- which source route a release type takes, what to assume for
    a missing pool diameter -- that should be made once and recorded.

    Two fields are *not* in the reduction and must not be left at their
    dataclass defaults, which is how this was found:

    **The levels of concern.** With `llc` at its default the flammable-mass
    integral evaluates `log(cc / clow)` at 28 and `SERIES` raises its overflow
    guard, on every trial. They are properties of the substance, so they come
    from `trialcase.SUBSTANCES`, and an unknown substance raises rather than
    silently taking a default.

    **The dispersion coefficients.** They follow from roughness, stability
    class and averaging time, so they are recomputed. The Monin-Obukhov length
    is stored, because the archive measures it and it does not follow from the
    class.
    """
    from degali.validation.rediphem import cases_from_reduced

    cases = cases_from_reduced(REDIPHEM_FULL)
    assert len(cases) == 313
    assert sum(1 for c in cases if c.usable) == 221

    burro = [c for c in cases if c.trial.series == "BURRO"]
    assert len(burro) == 8 and all(c.usable for c in burro)

    b9, = [c for c in burro if c.trial.name == "B9"]
    assert b9.case.gas.llc == 0.05 and b9.case.gas.ulc == 0.15
    assert b9.case.stability is not None
    assert b9.case.stability.rml == pytest.approx(-140.0)
    assert b9.case.source.rate[0] == pytest.approx(136.0)
    assert b9.case.source.radius[0] == pytest.approx(29.0)


@needs_rediphem_full
@pytest.mark.slow
def test_burro_pairs_the_model_at_the_sensor_height():
    """What the 61 pairs were, and how close that gets.

    The published statistic does not compare the model's ground-level
    centreline. It evaluates DEGADIS's own vertical profile at the sensor
    height:

        pred = centre * exp(-((z / sigma_z) ** (alpha + 1))) * 100

    At 1 m that is worth a factor of eight in the near field, far more than a
    Gaussian would give, because the exponent exceeds one and `sigma_z` is
    small close in. **The profile being too steep -- a defect this project
    documents separately -- is what makes the height matter so much here.**

    Two further details fix the sample at 61: arcs are *not* merged, and
    observations at or below 0.1 vol % are dropped. That cut removes twelve
    readings which are off-centreline sensors the plume missed -- 0.002 and
    0.004 vol % against neighbours reading 5 -- and they are what a variance
    statistic cannot survive: keeping them gives VG 1574.

    ===============  ===  ======  =======  =====
    pairing           n   MG      VG       FAC2
    ===============  ===  ======  =======  =====
    ground centreline 73  0.081   115000   0.04
    at 1 m, all obs   73  0.248   1574     0.45
    at 1 m, obs > 0.1 61  0.620   3.32     0.54
    published         61  0.811   2.53     0.56
    ===============  ===  ======  =======  =====

    The sample size and FAC2 match; MG is thirty per cent low and VG somewhat
    high. Whatever remains is smaller than every step taken to get here, and
    it is recorded rather than tuned away: neither interpolation order
    (0.620 against 0.625) nor the thermodynamic backend (CoolProp gives 0.576
    and n = 63) accounts for it.
    """
    pytest.importorskip("CoolProp")
    import math

    from degali.run import run_steady
    from degali.validation.rediphem import cases_from_reduced
    from degali.validation.statistics import statistics

    cases = [
        c for c in cases_from_reduced(REDIPHEM_FULL)
        if c.trial.series == "BURRO" and c.usable
    ]
    observed, predicted = [], []
    for case in cases:
        profile, source = run_steady(case.case, backend="legacy")
        rows = profile.rows
        x, centre, depth = rows[:, 0], rows[:, 1], rows[:, 7]
        alpha1 = source.alpha + 1.0
        for d, v in sorted(case.trial.arc_maxima(height=1.0).items()):
            if v <= 0.1 or not (x[0] <= d <= x[-1]):
                continue
            sigma_z = float(np.interp(d, x, depth))
            p = float(np.interp(d, x, centre))
            p *= math.exp(-((1.0 / sigma_z) ** alpha1)) * 100.0
            if p > 0.0:
                observed.append(v)
                predicted.append(p)

    assert len(observed) == 61, "the published sample size"
    s = statistics(observed, predicted)
    assert s.fac2 == pytest.approx(0.56, abs=0.03)
    assert s.mg == pytest.approx(0.620, abs=0.02)
    assert s.vg == pytest.approx(3.32, abs=0.15)

    # the height is what does it: on the ground centreline the same pairs
    # give an order of magnitude more
    ground = []
    for case in cases:
        profile, _ = run_steady(case.case, backend="legacy")
        rows = profile.rows
        x, centre = rows[:, 0], rows[:, 1]
        for d, v in sorted(case.trial.arc_maxima(height=1.0).items()):
            if v <= 0.1 or not (x[0] <= d <= x[-1]):
                continue
            ground.append(float(np.interp(d, x, centre)) * 100.0)
    on_ground = statistics(observed, ground)
    assert on_ground.mg == pytest.approx(0.212, abs=0.02)
    assert on_ground.fac2 < 0.10, "and it is not a near miss"
    assert s.mg / on_ground.mg > 2.5


def test_a_scope_check_on_a_name_that_is_not_a_limit_is_an_error():
    """A misspelled keyword used to be silently ignored.

    `check_range` looked each name up with `.get` and skipped a miss, so a
    typo turned a scope check into no check at all with nothing saying so.
    Raising instead immediately found a live instance: `lh2.assess` passes
    `diameter=pool_diameter` on both paths, and the jet range has no diameter
    limit, so a jet was being checked on four criteria while the caller
    believed it was five.

    `None` still means "not supplied", so a caller can pass everything it has
    without testing each one first. That is the whole distinction: an absent
    value is a fact about the release, an unknown name is a mistake in the
    code.
    """
    from degali.evidence import RANGE, check_range

    assert check_range("jet", diameter=None, wind=2.0) == []
    with pytest.raises(TypeError, match="no jet limit"):
        check_range("jet", diamter=9.1)
    with pytest.raises(TypeError, match="no jet limit"):
        check_range("jet", diameter=9.1)

    # a pool does have one, and a jet does not
    assert "diameter" in RANGE["pool"] and "diameter" not in RANGE["jet"]
    assert check_range("pool", diameter=9.1) == []


def test_the_liquid_hydrogen_machinery_does_not_perturb_degadis():
    """Importing the LH2 path must leave the reproduction claim untouched.

    Every switch that the liquid hydrogen work added is off by default, and
    the modules that add them do not mutate anything at import. Both halves
    matter: an option defaulting to on would change the ported model, and a
    module that patched a constant on import would change it only for
    programs that happened to import it -- which is worse, because the test
    suite would still pass whenever it ran the two in isolation.
    """
    from degali.core.jetplume import JetCoefficients, JetPlume
    from degali.io.inp import read_inp
    from degali.run import run_steady

    deck = REFERENCE_ROOT / "testcases" / "B9.INP"
    before = run_steady(read_inp(deck), backend="legacy")[0].rows.copy()

    import degali.lh2  # noqa: F401
    import degali.presets  # noqa: F401
    import degali.validation.nearfield  # noqa: F401

    after = run_steady(read_inp(deck), backend="legacy")[0].rows
    assert np.array_equal(before, after), "importing LH2 changed a DEGADIS run"

    # and every option is off, at its published DEGADIS 2.1 value
    c = JetCoefficients()
    assert (c.alfa1, c.alfa2, c.cd, c.sc, c.delta) == (0.057, 0.5, 0.2, 1.42, 2.15)
    assert not c.plume_transition
    assert not c.density_scaled_entrainment
    assert not c.vertical_shear
    assert not c.ground_layer_entrainment
    assert c.momentum_entrainment_beta == 0.0
    assert c.rise_drag == 0.0

    import inspect

    defaults = {
        n: p.default
        for n, p in inspect.signature(JetPlume.__init__).parameters.items()
        if isinstance(p.default, bool) or p.default == 0.0
    }
    for name in ("ground_effect", "non_boussinesq", "spread_floor"):
        assert defaults[name] is False, name
    for name in ("liquid_fraction", "evaporation_ratio"):
        assert defaults[name] == 0.0, name


@needs_e35_reduced
@pytest.mark.slow
@pytest.mark.parametrize("energy_transport", ["total", "enthalpy"])
def test_independent_energy_state_recovers_all_preslhy_interfaces(
    energy_transport,
):
    """The fifth state removes the two single-mixing-line rejections."""
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import (
        independent_energy_interfaces_from_reduced,
    )

    result = independent_energy_interfaces_from_reduced(
        E35_REDUCED, energy_transport=energy_transport
    )
    assert result.selected_trials == [10, 11, 12, 22, 23, 24, 25]
    assert result.all_interfaces_accepted
    assert not result.failures
    assert set(result.interfaces) == set(result.selected_trials)
    assert max(
        max(interface.relative_residuals.values())
        for interface in result.interfaces.values()
    ) < 1.0e-8
    assert max(
        interface.energy_quadrature_residual
        for interface in result.interfaces.values()
    ) < 1.0e-5
    assert max(
        interface.halfwidth_residual
        for interface in result.interfaces.values()
    ) < 0.05
    assert max(
        interface.temperature_residual
        for interface in result.interfaces.values()
    ) < 2.0
    assert result.interfaces[10].accepted
    assert result.interfaces[25].accepted
    assert all(
        interface.energy_transport == energy_transport
        for interface in result.interfaces.values()
    )
    assert all(
        interface.houf_width_mapping == "velocity"
        for interface in result.interfaces.values()
    )


def test_independent_energy_reduced_entry_points_expose_source_rate():
    """The raw and reduced PRESLHY routes must be able to use one source."""
    import inspect

    from degali.validation.nearfield import (
        independent_energy_from_reduced,
        independent_energy_interfaces_from_reduced,
    )

    for function in (
        independent_energy_from_reduced,
        independent_energy_interfaces_from_reduced,
    ):
        parameter = inspect.signature(function).parameters["source"]
        assert parameter.default == "flow_mean_gs"
        temperature = inspect.signature(function).parameters[
            "liquid_temperature_mode"
        ]
        assert temperature.default == "tank_saturation"
        measured = inspect.signature(function).parameters[
            "measured_source_table"
        ]
        assert measured.default is None
        spreading = inspect.signature(function).parameters[
            "fit_velocity_spreading"
        ]
        assert spreading.default is False
        mode = inspect.signature(function).parameters["measured_source_mode"]
        assert mode.default == "full"
        equilibrium = inspect.signature(function).parameters[
            "measured_lh2_equilibrium_bound"
        ]
        assert equilibrium.default is False

    field_filter = inspect.signature(independent_energy_from_reduced).parameters[
        "trial_filter"
    ]
    assert field_filter.default is None

    from degali.validation.nearfield import INDEPENDENT_ENERGY_TRIALS

    assert INDEPENDENT_ENERGY_TRIALS == (10, 11, 12, 22, 23, 24, 25)


def test_measured_pipe_source_is_topology_aware():
    from degali.validation.nearfield import (
        _measured_pipe_source_rows,
        _trial_source_inputs,
    )

    rows = _measured_pipe_source_rows(
        REFERENCE_ROOT / "preslhy" / "measured_pipe_source_2026-09-05.json"
    )
    rate, temperature, pressure = _trial_source_inputs(
        {"trial": 23, "flow_mean_gs": 121.098},
        source="flow_mean_gs",
        liquid_temperature_mode="tank_saturation",
        measured_rows=rows,
    )
    assert rate == pytest.approx(0.2371533942133789)
    assert temperature == pytest.approx(20.22341379150174)
    assert pressure == pytest.approx(2.249465067105)

    flow_rate, flow_temperature, flow_pressure = _trial_source_inputs(
        {"trial": 23, "flow_mean_gs": 121.098},
        source="flow_mean_gs",
        liquid_temperature_mode="tank_saturation",
        measured_rows=rows,
        measured_source_mode="flow_only",
    )
    assert flow_rate == pytest.approx(0.2371533942133789)
    assert flow_temperature is None
    assert flow_pressure is None

    nozzle_rate, nozzle_temperature, nozzle_pressure = _trial_source_inputs(
        {"trial": 23, "flow_mean_gs": 121.098},
        source="flow_mean_gs",
        liquid_temperature_mode="tank_saturation",
        measured_rows=rows,
        measured_source_mode="nozzle_only",
    )
    assert nozzle_rate == pytest.approx(0.121098)
    assert nozzle_temperature == pytest.approx(20.22341379150174)
    assert nozzle_pressure == pytest.approx(2.249465067105)

    rate, temperature, pressure = _trial_source_inputs(
        {"trial": 10, "flow_mean_gs": 189.425},
        source="flow_mean_gs",
        liquid_temperature_mode="tank_saturation",
        measured_rows=rows,
    )
    assert rate == pytest.approx(0.2843402449528671)
    assert temperature is None
    assert pressure is None


def test_momentum_dominance_distance_keeps_density_and_area():
    from types import SimpleNamespace

    from degali.validation.nearfield import _momentum_dominance_distance

    source = SimpleNamespace(density=2.0, velocity=3.0, area=4.0)
    distance = _momentum_dominance_distance(
        source, ambient_density=1.2, local_wind=5.0
    )
    assert distance == pytest.approx(math.sqrt(72.0 / 30.0))


def test_preslhy_liquid_temperature_is_independent_of_driving_pressure():
    """A pressure-driven subcooled liquid is not a 6 bar saturated liquid."""
    pytest.importorskip("CoolProp")
    from CoolProp.CoolProp import PropsSI

    from degali.validation.nearfield import _liquid_source_temperature

    trial = {"tanker_barg": 5.0}
    boiling = _liquid_source_temperature(trial, "ambient_boiling")
    tank_saturation = PropsSI(
        "T", "P", (trial["tanker_barg"] + 1.013) * 1.0e5,
        "Q", 0, "Hydrogen",
    )
    assert boiling == pytest.approx(20.36890353912106)
    assert tank_saturation == pytest.approx(28.267382228947373)
    assert _liquid_source_temperature(trial, "tank_saturation") is None
    with pytest.raises(ValueError, match="liquid temperature mode"):
        _liquid_source_temperature(trial, "fitted")


@needs_e35_reduced
@pytest.mark.slow
def test_one_trial_carries_the_far_arc_variance():
    """The consistent source removes the catastrophic far-arc variance.

    On the old pure-H2 table, trial 20 drove the 3--7 m band to VG 50.7 and
    had trial-only VG above 100.  Retaining the flash-entrained air and its
    momentum balance reduces campaign VG to 2.12; excluding trial 20 gives
    1.78.
    It remains the weakest low-wind trial, but no longer makes the campaign
    statistic singular.  The velocity-ratio applicability check is retained
    below because it still cannot identify this edge case.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import common_arcs, from_reduced
    from degali.validation.statistics import statistics

    _shipped, corrected = common_arcs(
        from_reduced(E35_REDUCED, corrections=False),
        from_reduced(E35_REDUCED, corrections=True),
    )
    far = [p for p in corrected.pairs if 3.0 <= p.x < 7.0]
    assert len(far) == 15

    def stats(pairs):
        return statistics(
            [p.observed for p in pairs], [p.predicted for p in pairs]
        )

    with_it = stats(far)
    without = stats([p for p in far if p.trial != 20])
    assert with_it.mg == pytest.approx(0.904, rel=0.05)
    assert with_it.vg == pytest.approx(2.115, rel=0.05)
    assert without.mg == pytest.approx(0.736, rel=0.05)
    assert without.vg == pytest.approx(1.775, abs=0.08)
    assert len([p for p in far if p.trial != 20]) == 13

    # and it is one trial, not a band: every other trial is tight
    for trial in sorted({p.trial for p in corrected.pairs} - {20}):
        s = stats([p for p in corrected.pairs if p.trial == trial])
        assert s.vg < 2.1, (trial, s.vg)
    s20 = stats([p for p in corrected.pairs if p.trial == 20])
    assert s20.vg < 25.0
    assert with_it.vg < 3.0

    # the filter cannot separate it: the lowest exit velocity in the campaign
    # passes a ratio test because the wind is smaller still
    import json

    from degali.validation.nearfield import MOMENTUM_RATIO, _exit_ratio

    data = {t["trial"]: t for t in json.loads(E35_REDUCED.read_text())["trials"]}
    speeds = {}
    for trial in (10, 11, 12, 20, 21, 22, 23, 24, 25):
        d = data[trial]
        ratio = _exit_ratio(d, d["flow_mean_gs"] / 1000.0)
        speeds[trial] = ratio * d["wind_ms"]
        assert ratio > MOMENTUM_RATIO, trial

    assert speeds[20] == min(speeds.values())
    assert speeds[11] / speeds[20] > 30.0


# ==========================================================================
# the buoyant trajectory at facility distances
# ==========================================================================

SPADEADAM = REFERENCE_ROOT / "spadeadam"
needs_spadeadam = pytest.mark.skipif(
    not (SPADEADAM / "sensors.csv").exists(),
    reason="the Spadeadam tables are not present",
)


def _spadeadam(test, *, ground_effect, roughness=0.001, **coefficients):
    """Run one Spadeadam horizontal release."""
    import csv
    import dataclasses

    from degali.validation.nearfield import (
        REACH, STEP, Trajectory, hydrogen_jet,
    )

    with open(SPADEADAM / "conditions.csv") as fh:
        c = {int(r["test"]): r for r in csv.DictReader(fh)}[test]
    wind = (float(c["wind_high_ms"]) + float(c["wind_low_ms"])) / 2.0
    jp, y0 = hydrogen_jet(
        rate=float(c["flow_kgs"]), diameter=float(c["orifice_mm"]) / 1000.0,
        wind=wind, height=0.5, ambient_temperature=277.15,
        relative_humidity=90.0, roughness=roughness,
        storage_pressure_barg=float(c["P04_barg"]),
        wind_reference_height=10.0, corrections=True,
        ground_effect=ground_effect,
    )
    if coefficients:
        jp.k = dataclasses.replace(jp.k, **coefficients)
    return Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=250.0).rows)


def _gradient(traj, x):
    """`c(1.8 m) / c(0.1 m)`, the metric fixed in the pre-registration."""
    low = traj.concentration_at(x, 0.0, 0.1)
    high = traj.concentration_at(x, 0.0, 1.8)
    return high / low if low > 0 else float("nan")


@needs_spadeadam
@pytest.mark.slow
def test_consistent_source_does_not_need_forced_ground_contact():
    """The corrected source grounds test 4 without permanently attaching it.

    At 30 m the free-to-detach solution has centre 1.43 m and a downward
    vertical concentration gradient, already matching the observed grounded
    character.  Forcing contact lowers it to 0.66 m but raises the farther
    concentrations.  The earlier decision to force contact was therefore an
    artefact of the inconsistent source table.
    """
    pytest.importorskip("CoolProp")

    off, on = _spadeadam(4, ground_effect=False), _spadeadam(4, ground_effect=True)

    assert off.at(30.0).z == pytest.approx(1.425, abs=0.06)
    assert on.at(30.0).z == pytest.approx(0.657, abs=0.06)

    assert _gradient(off, 30.0) == pytest.approx(0.833, abs=0.02)
    assert _gradient(on, 30.0) == pytest.approx(0.708, abs=0.02)
    assert _gradient(off, 30.0) < 1.0

    # and at every arc, not only the first
    for x in (30.0, 50.0, 100.0):
        assert _gradient(on, x) < _gradient(off, x), x
        assert _gradient(on, x) < 1.0, x

    # The free-to-detach run is inside the measured mean-to-peak band and is
    # closer on the farther arcs, where permanent attachment overpredicts.
    assert 8.4 < off.concentration_at(30.0, 0.0, 1.0) < 17.2
    assert abs(off.concentration_at(50.0, 0.0, 1.0) - 6.4) < \
        abs(on.concentration_at(50.0, 0.0, 1.0) - 6.4)


@needs_spadeadam
@pytest.mark.slow
def test_ground_contact_is_not_adopted_after_source_consistency():
    """It barely affects PRESLHY and still fails the independent regime test.

    Concentration statistics change negligibly, but the measured/modelled
    vertical-width ratio moves from 1.033 to about 1.01. That small geometric
    improvement does not override its wrong Spadeadam detachment behaviour,
    so the option stays experimental.
    """
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import common_arcs, from_reduced, vertical

    def near(ground):
        a, b = common_arcs(
            from_reduced(E35_REDUCED, corrections=False, ground_effect=ground),
            from_reduced(E35_REDUCED, corrections=True, ground_effect=ground),
        )
        return a.statistics(), b.statistics()

    (a0, b0), (a1, b1) = near(False), near(True)
    assert a1.mg == pytest.approx(a0.mg, abs=0.01)
    assert b1.mg == pytest.approx(b0.mg, abs=0.02)
    assert a1.vg == pytest.approx(a0.vg, rel=0.02)
    assert a1.fac2 == pytest.approx(a0.fac2, abs=0.01)

    def spread(ground):
        rows = vertical(
            E35_REDUCED, corrections=True, momentum_filter=True,
            ground_effect=ground,
        )
        return float(np.mean([
            r["modelled_sigma_z"] / r["measured_sigma_z"] for r in rows
        ]))

    assert spread(True) == pytest.approx(spread(False), abs=0.04)
    assert 0.98 < spread(True) < 1.04
    assert abs(spread(True) - 1.0) < abs(spread(False) - 1.0)


@needs_spadeadam
@pytest.mark.slow
def test_three_predicted_mechanisms_were_falsified():
    """Pre-registered, tried, and none of them is the fault.

    **M1, a lift-off criterion.** The model has none and every model in AEA's
    survey does, so it looked like the obvious gap. It is not: the bulk
    Richardson number runs 187 to 363 over the whole trajectory against
    Briggs' threshold of 20 to 30. **The criterion agrees the plume should
    lift.** Adding it changes the answer in the fourth decimal.

    **M2, shape-dependent pressure drag.** At a coefficient of 2.0 it moves
    the 30 m plume centre by about seven per cent.

    **M3, vertical velocity in the shear entrainment.** Mack's other change.
    Less than one per cent.

    The consistent source table resolves the large discrepancy before any of
    these optional mechanisms is applied; even the largest moves the 30 m
    height by less than eight per cent.
    """
    pytest.importorskip("CoolProp")

    base = _spadeadam(6, ground_effect=False)
    z0 = base.at(30.0).z

    # M1: the criterion fires nowhere, because Ri* is far above the threshold
    m1 = _spadeadam(6, ground_effect=False, liftoff_richardson=30.0)
    assert m1.at(30.0).z == pytest.approx(z0, abs=0.02)

    # M2: seven per cent at an implausible coefficient
    m2 = _spadeadam(6, ground_effect=False, shape_drag=2.0)
    assert 0.90 < m2.at(30.0).z / z0 < 0.95

    # M3: less than one per cent
    m3 = _spadeadam(6, ground_effect=False, vertical_shear=True)
    assert 0.98 < m3.at(30.0).z / z0 < 1.0

    # The consistent source has already brought the gradient below one; none
    # of these optional mechanisms materially changes the trajectory.
    for traj in (m1, m2, m3):
        assert abs(traj.at(30.0).z / z0 - 1.0) < 0.08


@needs_spadeadam
@pytest.mark.slow
def test_the_low_wind_release_detaches_while_test_four_stays_low():
    """The consistent source reproduces the observed regime split.

    The two releases have the same geometry and mass rate.  At the mean mast
    winds, high-wind test 4 remains below the 1.8 m array top at 30 m, while
    low-wind test 6 is above it and reaches 16.6 m at 100 m.  The experimental
    report and an independent EFFECTS analysis describe the same grounded /
    lifted split.  Permanent ground contact suppresses that result.
    """
    pytest.importorskip("CoolProp")

    four = _spadeadam(4, ground_effect=False)
    six = _spadeadam(6, ground_effect=False)
    held = _spadeadam(6, ground_effect=True)

    assert four.at(30.0).z < 1.8
    assert six.at(30.0).z > 1.8
    assert six.at(100.0).z > 10.0
    assert six.at(100.0).z > 1.5 * held.at(100.0).z
    assert four.at(100.0).z < six.at(100.0).z / 3.0


@needs_spadeadam
@pytest.mark.slow
def test_the_release_is_not_buoyancy_conserving():
    """Why the lift-off literature does not settle the low-wind case.

    Every lift-off criterion and correlation in AEA Technology's survey --
    Briggs' `Lp`, the Hall and Walker distance correlations, Hanna's
    exponential factor -- is derived from wind-tunnel plumes whose buoyancy
    flux is conserved. AEA says twice what that costs: "lift-off parameters
    based on non-dimensional fluxes are of limited use for non-buoyancy
    conserving flows", and "the lift-off distance correlations are unlikely to
    be valid for non-buoyancy conserving flows".

    A liquid hydrogen release is exactly such a flow, and by a wide margin.
    The thermodynamically consistent expanded source is about 3.05 kg/m3. It
    remains denser than air through three metres, crosses to buoyant near
    3.2 m, peaks near ten metres, then loses buoyancy to dilution.

    So the buoyancy flux at the source is negative, its maximum is somewhere
    downstream, and no single flux characterises the release. That is not a
    modelling choice; it is the substance. It is why source-only lift-off
    correlations are not used here.

    """
    pytest.importorskip("CoolProp")
    import numpy as np

    from degali.validation.nearfield import REACH, STEP, hydrogen_jet

    jp, y0 = hydrogen_jet(
        rate=0.833, diameter=0.0254, wind=2.5, height=0.5,
        ambient_temperature=277.15, relative_humidity=90.0,
        storage_pressure_barg=2.53, wind_reference_height=10.0,
        corrections=True, ground_effect=False,
    )
    # the source is denser than air, on both readings of "the source"
    assert jp.rhoe > jp.rhoa
    rows = jp.run(y0, distmx=STEP, smax=250.0).rows
    x, rho = rows[:, 0], rows[:, 8]

    deficit = jp.rhoa - rho
    assert deficit.min() < 0.0 or x[0] < 0.5   # dense at or before the start

    # The conserved source begins at about 0.82 m, stays dense through 3 m,
    # becomes buoyant between 3 and 5 m, peaks downstream, then decays.
    def at(distance):
        return float(np.interp(distance, x, deficit))

    assert at(1.0) < at(1.5) < at(3.0) < 0.0
    assert at(5.0) > 0.0
    assert at(10.0) > at(5.0) > at(30.0) > 0.0

    # which is what makes a single source buoyancy flux meaningless here
    assert at(1.0) == pytest.approx(-1.825, abs=0.04)
    assert at(30.0) == pytest.approx(0.0266, abs=0.006)


@needs_spadeadam
@pytest.mark.slow
def test_the_far_field_comparison_aggregated():
    """The 30 to 100 m comparison, through the packaged reader.

    Six arcs from the two horizontal releases, arc maximum against the model
    at the same heights:

    ===============  ===  ======  =====  =====  ==========
    ground_effect     n   MG      VG     FAC2   within 2x
    ===============  ===  ======  =====  =====  ==========
    off               6   1.245   1.37   0.83   5 of 6
    on                6   0.359   3.93   0.33   1 of 6
    ===============  ===  ======  =====  =====  ==========

    With the consistent source, free detachment has much smaller bias and
    variance. Five of six arcs improve; only the test-6 30 m peak reading
    favours the attached solution.

    Six arcs from two releases is a thin sample and this is not a validation.
    What it establishes is a direction and a size at a distance nothing in
    this package had reached.

    All seven outdoor releases carry signal on the 30 m arc. Five are
    downward releases, which impinge and spread as a ground-level source that
    this path does not model. Tests 8--15 are a separate closed-room and
    ventilation-mast campaign; their supply-nozzle conditions are not the
    atmospheric source conditions at the mast outlet. Both exclusions are
    physical properties of the experiment, and the reader returns all
    fifteen so the arithmetic remains visible.
    """
    pytest.importorskip("CoolProp")
    from degali.validation import spadeadam as sp
    from degali.validation.statistics import statistics

    trials = sp.load(SPADEADAM)
    assert len(trials) == 15
    assert [t.test for t in trials if t.outdoor] == list(range(1, 8))
    assert [t.test for t in trials if t.closed_room] == list(range(8, 16))
    assert [t.test for t in trials if t.downward_outdoor] == [1, 2, 3, 5, 7]
    horizontal = [t for t in trials if t.horizontal]
    assert [t.test for t in horizontal] == [4, 6]

    # the campaign arithmetic, before any model runs
    with_signal = [t for t in trials if t.outdoor and t.arc(30.0)
                   and max(t.arc(30.0).values()) >= 0.5]
    assert len(with_signal) == 7
    assert len([t for t in with_signal if t.downward_outdoor]) == 5

    def aggregate(ground_effect):
        arcs = [
            a for t in horizontal
            for a in sp.compare(
                t, corrections=True, ground_effect=ground_effect,
                roughness=0.001, wind_reference_height=10.0,
            )
        ]
        return arcs, statistics(
            [a.observed_max for a in arcs], [a.predicted_max for a in arcs]
        )

    off_arcs, off = aggregate(False)
    on_arcs, on = aggregate(True)
    assert len(off_arcs) == len(on_arcs) == 6

    assert off.mg == pytest.approx(1.245, rel=0.03)
    assert on.mg == pytest.approx(0.359, rel=0.03)
    assert off.vg == pytest.approx(1.373, rel=0.05)
    assert on.vg == pytest.approx(3.927, rel=0.05)
    assert off.fac2 == pytest.approx(0.83, abs=0.02)
    assert on.fac2 == pytest.approx(0.33, abs=0.02)

    # Free detachment passes bias and FAC2 and is close on VG; permanent
    # ground contact fails the bias and variance criteria badly.
    assert 0.7 < off.mg < 1.3
    assert off.fac2 > 0.5
    assert off.vg < 2.0
    assert on.vg > 3.0

    # Five arcs improve; only the test-6 30 m peak favours attachment.
    improved = 0
    for a, b in zip(off_arcs, on_arcs):
        assert a.radius == b.radius and a.test == b.test
        improved += abs(a.observed_max / a.predicted_max - 1.0) < \
            abs(b.observed_max / b.predicted_max - 1.0)
    assert improved == 5


@needs_spadeadam
def test_a_downward_release_is_refused_rather_than_guessed():
    """Five outdoor releases point at the ground.

    A downward jet impinges within centimetres and spreads as a ground-level
    source. That is a different source term from the horizontal jet path, and
    the tempting thing is to run it anyway with the release height as a
    stand-in, which would produce numbers that look like a five-fold larger
    sample and mean nothing.

    `spadeadam.compare` raises instead. The eight closed-room/ventilation-mast
    tests are refused by both atmospheric source paths until mast-outlet
    conditions are available. All trials stay loaded and classified, so the
    exclusions are visible rather than hidden in a filter.
    """
    from degali.validation import spadeadam as sp

    trials = {t.test: t for t in sp.load(SPADEADAM)}
    assert sum(t.downward_outdoor for t in trials.values()) == 5
    assert sum(t.closed_room for t in trials.values()) == 8
    assert trials[8].mast_source is not None
    assert trials[8].mast_source.velocity == pytest.approx(6.4)
    assert trials[8].mast_source.temperature == pytest.approx(125.15)
    assert trials[14].mast_source.mass_flow == pytest.approx(0.364)
    assert trials[15].mast_source is None  # stack lost in the Test 14 explosion
    for n in range(8, 15):
        source = trials[n].mast_source
        assert source is not None
        assert source.density * source.volume_flow == pytest.approx(
            source.mass_flow, rel=0.015
        )

    down = trials[5]
    assert down.downward_outdoor
    assert down.rate == pytest.approx(0.739)
    # it carries real signal, which is what makes refusing it a choice
    assert max(down.arc(30.0).values()) > 5.0
    with pytest.raises(ValueError, match="downward"):
        sp.compare(down, corrections=True)

    mast = trials[15]
    with pytest.raises(ValueError, match="closed-room/ventilation-mast"):
        sp.compare(mast, corrections=True)
    with pytest.raises(ValueError, match="closed-room/ventilation-mast"):
        sp.compare_downward(mast)


@needs_spadeadam
@pytest.mark.slow
def test_the_ventilation_mast_is_a_separate_negative_constraint():
    """The recovered mast source agrees with nondetects but over-rises.

    This comparison was pre-registered in ``docs/prereg-spadeadam-mast.md``.
    It uses the actual vertical 450 mm outlet for tests 8--14, wind-aligned
    sensor coordinates, and no log statistic across field nondetects.

    No sensor below the 0.5 vol % drift threshold is falsely predicted above
    the 4 vol % LFL. That is only a weak negative constraint. The more useful
    result is independent corroboration of the trajectory fault: 18 downwind
    sensors did register at least 0.5 vol %, ten of them within 3.5 m of the
    plume axis, while the model puts its centre more than eight metres above
    every one and predicts effectively zero at them.
    """
    pytest.importorskip("CoolProp")
    from degali.validation import spadeadam as sp

    trials = {t.test: t for t in sp.load(SPADEADAM)}
    rows = [
        row for n in range(8, 15)
        for row in sp.compare_mast(trials[n])
    ]

    assert len(rows) == 111
    assert sum(row.false_flammable for row in rows) == 0

    detected = [row for row in rows if not row.censored]
    assert len(detected) == 18
    near_axis = [row for row in detected if abs(row.y) <= 3.5]
    assert len(near_axis) == 10
    assert max(row.predicted for row in detected) < 1.0e-6
    assert min(row.centre - row.z for row in detected) > 8.0

    with pytest.raises(ValueError, match="not a ventilation-mast"):
        sp.compare_mast(trials[4])
    with pytest.raises(ValueError, match="no intact"):
        sp.compare_mast(trials[15])


@needs_spadeadam
@pytest.mark.slow
def test_the_downward_releases_show_the_same_fault_independently():
    """Five more tests, through a different plume model, same fault.

    The downward releases cannot be run as jets: at 0.32 m with a 25.4 mm
    orifice the fluid reaches the ground inside the jet's own development
    length, and `initial_conditions_directed` refuses them for that reason.
    The fluid arrives essentially undiluted, so what is left is a ground-level
    source of unknown footprint.

    **The footprint was measured before a value was picked for it.** Sweeping
    it eighty-fold, 0.05 m to 4 m, moves the 30 m concentration by seven per
    cent and the plume centre by half a metre. The unknown does not matter,
    which is what makes these five outdoor tests usable at all.

    They go through `addons.LiftoffPlume` -- the URAHFREP ground-truncated
    buoyant plume -- and not through `JetPlume`. So the two halves of the
    Spadeadam comparison are independent code paths:

    ==========================================  ========  =====  ====  =====
    comparison                                   n         MG     VG    FAC2
    ==========================================  ========  =====  ====  =====
    horizontal, JetPlume, ground contact off     6 arcs    1.25   1.37  0.83
    downward, LiftoffPlume                       5 tests   4.79  12.32  0.00
    ==========================================  ========  =====  ====  =====

    **Two buoyant plume models, written by different people for different
    purposes, underpredict concentration on the same campaign, but the
    conserved JetPlume is much closer.** The downward impingement path remains
    strongly biased and cannot be used as a quantitative validation of the
    horizontal source correction.

    It is worth noting what the URAHFREP model is: AEA Technology's own, built
    for this problem, and reported by them as over-predicting rise. It does so
    here too.
    """
    pytest.importorskip("CoolProp")
    from degali.validation import spadeadam as sp
    from degali.validation.statistics import statistics

    trials = {t.test: t for t in sp.load(SPADEADAM)}

    # the footprint does not matter
    spread = [
        sp.compare_downward(trials[5], footprint=d)[1]
        for d in (0.05, 0.5, 4.0)
    ]
    assert max(spread) / min(spread) < 1.15

    observed, predicted, centres = [], [], []
    for n in sp.OUTDOOR_DOWNWARD:
        out = sp.compare_downward(trials[n])
        assert out is not None, n
        observed.append(out[0])
        predicted.append(out[1])
        centres.append(out[2])

    s = statistics(observed, predicted)
    assert s.n == 5
    assert s.mg == pytest.approx(4.79, rel=0.05)
    assert s.vg == pytest.approx(12.32, rel=0.05)
    assert s.fac2 == pytest.approx(0.0, abs=0.01)

    # and the plume is put metres up on every one of them
    assert min(centres) > 2.0
    assert max(centres) > 5.0

    # same bias direction, but the unresolved impingement path is much worse
    jet_mg = 1.245
    assert 3.0 < s.mg / jet_mg < 4.5

    # a horizontal test is refused by this path, as a downward one is by the
    # other: neither is quietly run through the wrong source term
    with pytest.raises(ValueError, match="horizontal"):
        sp.compare_downward(trials[4])


@pytest.mark.slow
def test_the_nasa_vertical_profile_shows_the_same_over_rise():
    """A third campaign, and the fault scales with wind the same way.

    The package used one number per NASA test: the minimum height at which a
    flammable cloud was measured, Witcofski's Table 4. **Table 3 is a vertical
    profile** -- the maximum concentration at the furthest tower row, at 1,
    9.4 and 18.6 m, for the same four tests -- and it had never been used. It
    is the quantity that constrains a trajectory, and it is transcribed from
    the printed table rather than digitised from a figure.

    Against the height at which the measurement peaks:

    ======  ======  ==============  =============  =====
    test    wind     measured peak   model centre   ratio
    ======  ======  ==============  =============  =====
    5       6.30     9.4 m           6.2 m          0.66
    4       3.35     9.4 m           14.8 m         1.57
    6       2.20     ~14 m           23.9 m         1.71
    2       1.55     >= 18.6 m       40.3 m         > 2.2
    ======  ======  ==============  =============  =====

    **The model over-predicts the height, and monotonically worse as the wind
    falls.** That is the Spadeadam finding on a different campaign, four years
    of instrumentation apart, at ten times the mass flow, from a 9.1 m pond
    instead of a 25 mm orifice, through the ground-level buoyant path rather
    than the jet.

    Three campaigns, two plume models, one fault.

    Test 5 is the only one the model puts too low, and it is the highest wind
    at 6.3 m/s -- twice any other. The paper notes the same thing from the
    data side: its cloud travelled downstream fastest and gave the highest far
    concentration, 29.2 % at 9.4 m.
    """
    pytest.importorskip("CoolProp")
    import numpy as np

    from degali.lh2 import assess
    from degali.validation import witcofski as w

    # tests 3 and 7 have published conditions and no concentration data
    assert set(w.FAR_TOWER) == {2, 4, 5, 6}
    assert not w.SPILLS[3]["data"] and not w.SPILLS[7]["data"]
    assert 3 in w.RATES and 7 in w.RATES   # the trap: conditions exist

    heights = {}
    for test, profile in w.FAR_TOWER.items():
        spill = w.SPILLS[test]
        result = assess(
            rate=w.RATES[test], pool_diameter=w.POND,
            wind=sum(spill["wind"]) / 2.0,
            ambient_temperature=spill["tamb"] + 273.15,
            relative_humidity=spill["rh"], max_distance=120.0,
            at_distance=w.FAR_TOWER_DISTANCE,
        )
        traj = result.trajectory
        assert len(traj) and traj[-1, 0] >= w.FAR_TOWER_DISTANCE, test
        heights[test] = float(
            np.interp(w.FAR_TOWER_DISTANCE, traj[:, 0], traj[:, 1])
        )

    assert heights[5] == pytest.approx(6.2, abs=0.4)
    assert heights[4] == pytest.approx(14.8, abs=0.6)
    assert heights[6] == pytest.approx(23.9, abs=1.0)
    assert heights[2] == pytest.approx(40.3, abs=1.5)

    # monotone in wind: the weaker the wind, the higher the model puts it
    winds = {t: sum(w.SPILLS[t]["wind"]) / 2.0 for t in heights}
    order = sorted(heights, key=lambda t: winds[t])
    assert [heights[t] for t in order] == sorted(
        (heights[t] for t in order), reverse=True
    )

    # and against the height where the bottles peak, it is high on three of
    # four, with the exception being the fastest wind
    peak_at = {}
    for test, profile in w.FAR_TOWER.items():
        bottles = {z: v[0] for z, v in profile.items() if v[0] is not None}
        peak_at[test] = max(bottles, key=bottles.get) if max(bottles.values()) > 0 else 18.6
    assert peak_at[4] == 9.4 and peak_at[5] == 9.4 and peak_at[6] == 18.6
    assert heights[5] < peak_at[5]
    for test in (2, 4, 6):
        assert heights[test] > peak_at[test], test


@needs_e35_reduced
@needs_spadeadam
@pytest.mark.slow
def test_buoyant_entrainment_is_the_only_mechanism_that_touches_wind_dependence():
    """The one thing AEA recommended and nobody tried, and it works.

    Having found that suppressing rise through the vertical velocity spoils
    the dilution, AEA Technology wrote one sentence about what to do instead:

        Modifying the entrainment formula to depend upon a buoyant velocity
        scale rather than vertical component of velocity may circumvent such
        problems

    `sqrt(g' H)` does not scale with the wind, and every other mechanism tried
    here does. That is why it is the only one that changes the *shape* of the
    error rather than its size: the model over-predicts the plume height by
    0.66 at 6.3 m/s and by more than 2.2 at 1.55 m/s, and a correction that
    multiplies the rise by a constant cannot fix both.

    ============  =========  =========  ==========  =========  =========
    coefficient    near MG    near VG    near FAC2   sigma_z    far z6/z4
    ============  =========  =========  ==========  =========  =========
    0              1.132      19.15      0.79        1.37       2.46
    2              0.950      2.76       0.83        1.53       2.22
    4              0.938      1.73       0.85        1.81       1.88
    8              1.107      1.75       0.82        2.47       1.47
    ============  =========  =========  ==========  =========  =========

    **AEA's coupling prediction is falsified on their own escape route.** The
    pre-registration said the near-field dilution would degrade unless the
    buoyant scale was used, and set a five per cent tolerance on MG. It does
    not degrade: MG moves *towards* unity and the variance collapses by a
    factor of seven, on data the corrections were never tuned against.

    **It is not adopted.** The coefficient has no source, and the
    pre-registration says a mechanism that works is not adopted until it does.
    What it costs is visible in the same table: after correcting the Gaussian
    width convention, this pushes an already over-wide 1.37 to 1.53 at a
    coefficient of 2 and 1.81 at 4.
    """
    pytest.importorskip("CoolProp")
    import dataclasses

    import numpy as np

    import degali.validation.nearfield as nf
    from degali.validation import spadeadam as sp
    from degali.validation.nearfield import (
        REACH, STEP, Trajectory, common_arcs, from_reduced, hydrogen_jet,
        vertical,
    )

    original = nf.hydrogen_jet

    def near(coefficient):
        def patched(*a, **kw):
            kw["source_table_consistency"] = False
            kw["source_momentum_consistency"] = False
            jp, y0 = original(*a, **kw)
            jp.k = dataclasses.replace(
                jp.k, buoyant_entrainment=coefficient
            )
            return jp, y0

        nf.hydrogen_jet = patched
        try:
            shipped, _ = common_arcs(
                from_reduced(E35_REDUCED, corrections=False),
                from_reduced(E35_REDUCED, corrections=True),
            )
            rows = vertical(
                E35_REDUCED, corrections=True, momentum_filter=True
            )
            spread = float(np.mean([
                r["modelled_sigma_z"] / r["measured_sigma_z"] for r in rows
            ]))
            return shipped.statistics(), spread
        finally:
            nf.hydrogen_jet = original

    trials = {t.test: t for t in sp.load(SPADEADAM)}

    def far(test, coefficient):
        t = trials[test]
        jp, y0 = hydrogen_jet(
            rate=t.rate, diameter=t.orifice, wind=t.wind_low, height=0.50,
            ambient_temperature=277.15, relative_humidity=90.0,
            roughness=0.001, storage_pressure_barg=t.line_pressure,
            wind_reference_height=10.0, corrections=True, ground_effect=True,
            source_table_consistency=False,
            source_momentum_consistency=False,
        )
        jp.k = dataclasses.replace(jp.k, buoyant_entrainment=coefficient)
        traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=250.0).rows)
        return traj.at(30.0).z

    base, base_spread = near(0.0)
    two, two_spread = near(2.0)
    four, four_spread = near(4.0)

    # the variance collapses; the bias moves towards unity
    assert base.vg == pytest.approx(19.15, rel=0.05)
    assert two.vg == pytest.approx(2.76, rel=0.08)
    assert four.vg == pytest.approx(1.73, rel=0.08)
    assert abs(two.mg - 1.0) < abs(base.mg - 1.0)
    assert two.fac2 > base.fac2

    # AEA's coupling prediction, on their own escape route: falsified
    assert two.mg > 0.9, "the dilution does not degrade"

    # the cost is the quantity the corrections were tuned on
    assert base_spread == pytest.approx(1.37, abs=0.04)
    assert two_spread == pytest.approx(1.53, abs=0.06)
    assert four_spread > 1.7

    # and the wind dependence compresses, which nothing else did
    ratios = {c: far(6, c) / far(4, c) for c in (0.0, 2.0, 4.0, 8.0)}
    assert ratios[0.0] == pytest.approx(2.46, rel=0.05)
    assert ratios[8.0] < 1.6
    assert list(ratios.values()) == sorted(ratios.values(), reverse=True)

    # default off, because the coefficient has no source
    from degali.core.jetplume import JetCoefficients

    assert JetCoefficients().buoyant_entrainment == 0.0
