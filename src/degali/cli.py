"""Command-line interface.

The original ships six executables that pass files between them; ``degali``
runs the equivalent stages in one process and prints a summary:

.. code-block:: console

    degali steady    B9.INP --er1 EXAMPLE.ER1 --er2 EXAMPLE.ER2
    degali transient B9T.INP --snapshot 60 --snapshot 120
    degali jet       EX1.INO
    degali jet       EX2.INO --bridge EX2.IN
    degali dose      B9T.INP --at 200 --at 400 --at 800
    degali steady    B9.INP --backend coolprop
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .run import (
    Receptor,
    run_jet,
    run_jet_to_ground,
    run_steady,
    run_transient,
)


def _profile_table(rows: np.ndarray, limit: int = 20) -> str:
    head = f"{'dist (m)':>10} {'mole frac':>11} {'kg/m3':>11} {'T (K)':>8} {'Sz (m)':>8} {'Sy (m)':>8}"
    step = max(1, len(rows) // limit)
    lines = [head, "-" * len(head)]
    for r in rows[::step]:
        lines.append(
            f"{r[0]:10.4g} {r[2]:11.4g} {r[3]:11.4g} {r[5]:8.1f} {r[7]:8.3g} {r[8]:8.3g}"
        )
    return "\n".join(lines)


def _cmd_steady(args) -> int:
    profile, src = run_steady(
        args.deck, er1=args.er1, er2=args.er2, backend=args.backend
    )
    gas = src.case.gas
    fluid = gas.coolprop_name
    props = f"CoolProp ({fluid})" if fluid else src.thermo.backend.name
    print(f"Steady release: {gas.name}, {src.blanket.ess:.4g} kg/s "
          f"[properties: {props}]")
    print(f"  wind-profile exponent alpha = {src.alpha:.5f}")
    print(f"  secondary source: {src.blanket.outl:.3g} m long, "
          f"{src.blanket.outb:.3g} m half-width")
    print(f"  profile: {len(profile.rows)} points, "
          f"{profile.n_dense} in the dense phase, "
          f"transition at {profile.transition:.4g} m")
    for name, level in (("upper", gas.ulc), ("lower", gas.llc)):
        d = profile.distance_to(level)
        shown = "not reached" if np.isnan(d) else f"{d:.4g} m"
        print(f"  distance to the {name} level of concern "
              f"({level * 100:.4g} mol %): {shown}")
    print(f"  mass above the lower level of concern: "
          f"{profile.mass_above_lfl:.5g} kg")
    print(f"  mass between the two levels: {profile.mass_between:.5g} kg")
    print()
    print(_profile_table(profile.rows))
    return 0


def _cmd_transient(args) -> int:
    times = np.array(args.snapshot, dtype=float) if args.snapshot else None
    out = run_transient(
        args.deck, er1=args.er1, er2=args.er2, times=times,
        backend=args.backend,
    )
    gas = out.source.case.gas
    print(f"Transient release: {gas.name}, "
          f"{len(out.field.observers)} observers, {len(out.snapshots)} snapshots")
    head = (f"{'t (s)':>8} {'points':>7} {'extent (m)':>11} "
            f"{'peak mol frac':>14} {'mass > LLC (kg)':>16}")
    print()
    print(head)
    print("-" * len(head))
    for s in out.snapshots:
        print(f"{s.time:8.4g} {len(s):7d} {s.column('dist')[-1]:11.4g} "
              f"{s.column('yc').max():14.5g} {s.mass_above_llc:16.5g}")
    peak = max(out.snapshots, key=lambda s: s.mass_above_llc, default=None)
    if peak is not None:
        print(f"\npeak flammable mass {peak.mass_above_llc:.5g} kg at "
              f"t = {peak.time:.4g} s")
    return 0


def _cmd_dose(args) -> int:
    out = run_transient(
        args.deck, er1=args.er1, er2=args.er2, backend=args.backend
    )
    receptors = [Receptor(x=x) for x in args.at]
    histories = out.dose(receptors)
    print(f"Concentration histories at {len(histories)} receptors")
    head = (f"{'x (m)':>9} {'points':>7} {'peak mol frac':>14} "
            f"{'at t (s)':>9} {'dose (mol frac.s)':>18}")
    print()
    print(head)
    print("-" * len(head))
    for h in histories:
        peak, t = h.peak
        print(f"{h.receptor.x:9.4g} {len(h.rows):7d} {peak:14.5g} "
              f"{t:9.4g} {h.dose():18.5g}")
    return 0


def _cmd_jet(args) -> int:
    if args.bridge is None:
        jet, deck = run_jet(args.deck, backend=args.backend)
        print(f"Jet release: {deck.erate:.4g} kg/s through a "
              f"{deck.diajet:.3g} m orifice at {deck.elejet:.3g} m")
        if not jet.touchdown:
            print("  the plume never reaches the ground at a concentration "
                  "of interest")
            print(f"  final elevation {jet.rows[-1, 1]:.4g} m at "
                  f"{jet.rows[-1, 0]:.4g} m downwind")
            return 0
        print(f"  touchdown at {jet.distance:.5g} m, "
              f"{jet.concentration:.5g} kg/m3, "
              f"half-width {jet.halfwidth:.5g} m")
        print("  pass --bridge with the .IN deck to continue downwind")
        return 0

    profile, jet, src = run_jet_to_ground(
        args.deck, args.bridge, er1=args.er1, er2=args.er2,
        backend=args.backend,
    )
    gas = src.case.gas
    print(f"Jet release: touchdown at {jet.distance:.5g} m")
    print(f"  bridged to a {jet.halfwidth:.4g} m source, "
          f"diluted to {src.case.source.wc[0]:.4g} mass fraction")
    d = profile.distance_to(gas.llc)
    shown = "not reached" if np.isnan(d) else f"{d:.4g} m"
    print(f"  distance to the lower level of concern "
          f"({gas.llc * 100:.4g} mol %): {shown}")
    print()
    print(_profile_table(profile.rows))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="degali",
        description=(
            "DEGALI - Dense Gas Dispersion for Liquid Hydrogen, built on a "
            "verified DEGADIS 2.1 reimplementation"
        ),
    )
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, deck_help):
        sp.add_argument("deck", type=Path, help=deck_help)
        sp.add_argument("--er1", type=Path, help="source-model parameter file")
        sp.add_argument("--er2", type=Path, help="downwind parameter file")
        sp.add_argument(
            "--backend", choices=("legacy", "coolprop"), default="legacy",
            help="property model: 'legacy' reproduces DEGADIS 2.1 exactly, "
                 "'coolprop' uses equations of state and solves the "
                 "thermodynamic inversions accurately (default: legacy)",
        )

    s = sub.add_parser("steady", help="a steady ground-level release")
    common(s, "the .INP input deck")
    s.set_defaults(func=_cmd_steady)

    t = sub.add_parser("transient", help="an unsteady ground-level release")
    common(t, "the .INP input deck")
    t.add_argument(
        "--snapshot", type=float, action="append",
        help="time to report the cloud at, in seconds; repeatable. "
             "Defaults to the window DEGADIS chooses itself.",
    )
    t.set_defaults(func=_cmd_transient)

    d = sub.add_parser("dose", help="concentration history at fixed receptors")
    common(d, "the .INP input deck")
    d.add_argument(
        "--at", type=float, action="append", required=True,
        help="downwind distance of a receptor, in metres; repeatable",
    )
    d.set_defaults(func=_cmd_dose)

    j = sub.add_parser("jet", help="a pressurised release")
    common(j, "the .INO jet deck")
    j.add_argument(
        "--bridge", type=Path,
        help="the .IN jet deck; supplying it continues past touchdown "
             "into the ground-level model",
    )
    j.set_defaults(func=_cmd_jet)

    args = p.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, RuntimeError) as exc:
        print(f"degali: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
