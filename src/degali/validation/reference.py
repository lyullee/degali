"""Reference oracle: run the original DEGADIS 2.1 Fortran and read its state.

``degali`` is validated against the original code, not against printed
output.  Printed output carries only 5-6 significant figures, which is not
enough to distinguish a correct port from one that is merely close.  The
:mod:`reference` package therefore builds the Fortran from source, runs it, and
extracts full double-precision internal state through a purpose-built probe.

Layout::

    reference/
        fortran/        patched DEGADIS 2.1 sources + build.sh
        testcases/      B9, B9T, EX1, EX2, EX3 inputs and EPA golden output
        probe.f         the state-dump program linked against the objects

The patches applied to the Fortran are enumerated in
:data:`PORTABILITY_PATCHES` and are all portability fixes; none changes the
numerics.  Two of the five test cases reproduce EPA's 2012 golden output
byte for byte apart from timestamps, and the other three differ only in the
seventh significant figure, which the EPA readme explicitly calls normal.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

#: Portability changes made to the 1989/2012 Fortran so it builds with
#: gfortran on a case-sensitive filesystem.  Recorded here so the port is
#: auditable; none of them alters a computed value.
PORTABILITY_PATCHES = {
    "include-case": (
        "INCLUDE statements name 'DEG1.prm' while the files are DEG1.PRM. "
        "VMS and DOS ignore case; Linux does not."
    ),
    "eof-marker": "DEG1.for ends with a DOS ^Z (0x1A) byte.",
    "getarg": (
        "GETARG(n, buf, status) is a DEC 3-argument extension; gfortran takes "
        "two arguments. Replaced with GETARG + LEN_TRIM."
    ),
    "empty-arglist": (
        "FUNCTION TUPF / FUNCTION TDNF declared without '()', which Fortran 90 "
        "and later reject."
    ),
    "ms-extensions": (
        "GETTIM, GETDAT and SYSTEMQQ are Microsoft Fortran 5.0 extensions. "
        "Reimplemented over DATE_AND_TIME and SYSTEM in MSSHIM.for."
    ),
    "extension-case": (
        "File-name extensions appear as both '.ER2' and '.er2' in different "
        "units; normalised to lower case throughout."
    ),
    "d-lines": (
        "Debug lines flagged by 'd' in column 1 (a VAX convention); compiled "
        "with -fd-lines-as-comments."
    ),
    "estrt1-equivalence": (
        "COMMON /ERROR/ interleaves 18 REAL*8, one INTEGER*4 (NOBLpt) and two "
        "more REAL*8. ESTRT1 overlays read buffers onto it with EQUIVALENCE, "
        "which forces gfortran to pad that common differently from every other "
        "unit; CRFGER and EPSILON were silently read as zero, disabling air "
        "entrainment into the source blanket. Intel's /align:dcommons masked "
        "this. The overlay is replaced by explicit assignment; values are "
        "unchanged."
    ),
}

#: Test cases shipped by EPA, with the program sequence each one runs.
TEST_CASES = {
    "b9": ("steady-state, ground-level plume (Burro 9)", ["deg1", "deg2s"]),
    "b9t": ("transient, ground-level plume (Burro 9)", ["deg1", "deg2", "deg3"]),
    "ex1": ("vertical jet, no touchdown (MIC)", ["jetplu", "degbridg"]),
    "ex2": (
        "vertical jet with touchdown (ammonia)",
        ["jetplu", "degbridg", "deg1", "deg2s"],
    ),
    "ex3": (
        "vertical jet with touchdown, simplified density",
        ["jetplu", "degbridg", "deg1", "deg2s"],
    ),
}


@dataclass
class ReferenceRun:
    """A completed run of the original Fortran."""

    case: str
    workdir: Path
    state: dict

    def __getitem__(self, key: str):
        return self.state[key]

    @property
    def listing(self) -> str:
        """Concatenated ``.lis`` output, as the batch files assemble it."""
        parts = []
        for ext in ("out", "scl", "sr3"):
            p = self.workdir / f"{self.case}.{ext}"
            if p.exists():
                parts.append(p.read_text(errors="replace"))
        return "".join(parts)

    @property
    def golden(self) -> str:
        """EPA's 2012 golden listing for this case."""
        return (self.workdir / f"{self.case}.ref.lis").read_text(errors="replace")


class Reference:
    """Driver for the ported Fortran reference implementation."""

    def __init__(self, root: Path | str):
        self.root = Path(root)
        self.bin = self.root / "fortran" / "bin"
        self.testcases = self.root / "testcases"

    def build(self) -> None:
        subprocess.run(
            ["bash", "build.sh"], cwd=self.root / "fortran", check=True,
            capture_output=True,
        )

    def run(self, case: str, workdir: Path | str) -> ReferenceRun:
        """Run ``case`` in ``workdir`` and return its full-precision state."""
        case = case.lower()
        if case not in TEST_CASES:
            raise KeyError(f"unknown case {case!r}; have {sorted(TEST_CASES)}")
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)

        for src in self.testcases.iterdir():
            if src.stem.lower() == case:
                dst = work / src.name.lower()
                dst.write_bytes(src.read_bytes().replace(b"\r\n", b"\n"))
        for ext in ("er1", "er2", "er3"):
            src = self.testcases / f"EXAMPLE.{ext.upper()}"
            if src.exists():
                (work / f"{case}.{ext}").write_bytes(
                    src.read_bytes().replace(b"\r\n", b"\n")
                )
        golden = work / f"{case}.lis"
        if golden.exists():
            shutil.move(str(golden), str(work / f"{case}.ref.lis"))

        for prog in TEST_CASES[case][1]:
            subprocess.run(
                [str(self.bin / prog), case], cwd=work, capture_output=True
            )

        state = {}
        probe = self.bin / "probe"
        if probe.exists() and (work / f"{case}.inp").exists():
            subprocess.run([str(probe), case], cwd=work, capture_output=True)
            pj = work / "probe.json"
            if pj.exists():
                state = json.loads(pj.read_text())

        return ReferenceRun(case=case, workdir=work, state=state)
