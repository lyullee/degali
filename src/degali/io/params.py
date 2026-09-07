"""Reading the ``.ER1`` / ``.ER2`` / ``.ER3`` numerical parameter files.

These carry every tolerance, weight and closure constant DEGADIS uses.  They
sit outside the input deck so that a user can retune convergence without
touching the physical case, and the EPA distribution ships one canonical set
(``EXAMPLE.ER1`` and friends) that all the test cases use.

The format is positional and column-sensitive at the same time.  ``ESTRT1``
reads with ``FORMAT(10X, G10.4)``: the label occupies columns 1-10, the value
columns 11-20, and everything past column 20 is a comment.  Lines beginning
with ``!`` are skipped.  The *labels are never parsed* -- values are assigned
purely in order, so a deck with the right numbers in the wrong order is read
without complaint.

This reader keeps that ordering as the source of truth but also returns the
labels, so a caller can check them if it wants to.  :func:`read_er1` returns a
:class:`~degali.core.driver.DriverParameters`, ready to hand to
:class:`~degali.core.driver.SourceRun`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..core.driver import DriverParameters

#: Field order in ``.ER1``, exactly as ``ESTRT1`` consumes it.  The names are
#: the Fortran variables, not the labels in the file (which differ in a few
#: places: ``delrhomin`` is ``DELRMN``, ``vuc`` is ``VUE``).
ER1_FIELDS = (
    "stpin", "erbnd", "wtrg", "wttm", "wtya", "wtyc", "wteb", "wtmb",
    "xli", "xri", "eps", "zlow", "stpinz", "erbndz",
    "srcoer", "srcss", "srccut", "ernobl", "noblpt",
    "crfger", "epsilon",
    "ce", "delrmn",
    "szstp0", "szerr", "szsz0",
    "ialpfl", "alpco",
    "iphifl", "dellay",
    "vua", "vub", "vue", "vud", "vudelta",
)

#: Field order in ``.ER2``, as ``ESTRT2`` consumes it.  The first block lands
#: in ``/ERRORP/``, the second in ``/STP/``, and the last value in ``/CNOBS/``.
ER2_FIELDS = (
    "sy0er", "erro", "sz0er", "wtaio", "wtqoo", "wtszo",
    "errp", "smxp", "wtszp", "wtsyp", "wtbep", "wtdh",
    "errg", "smxg", "ertdnf", "ertupf", "wtruh", "wtdhg",
    "stpo", "stpp", "odlp", "odllp", "stpg", "odlg", "odllg",
    "nobs",
)

#: Fields that ``ESTRT1``/``ESTRT2`` round to an integer with ``NINT``.
_INTEGER_FIELDS = {"noblpt", "ialpfl", "iphifl", "nobs"}


@dataclass
class NumericalParameters:
    """Raw contents of a parameter file, keyed by Fortran variable name."""

    values: dict[str, float]
    labels: list[str]
    path: Path | None = None

    def __getitem__(self, key: str) -> float:
        return self.values[key]

    def get(self, key: str, default=None):
        return self.values.get(key, default)


def _data_lines(path: str | Path) -> list[str]:
    text = Path(path).read_text(errors="replace").replace("\r\n", "\n")
    return [ln for ln in text.split("\n") if ln.strip() and not ln.startswith("!")]


def read_parameters(path: str | Path, fields: tuple[str, ...]) -> NumericalParameters:
    """Read a parameter file positionally.

    ``ESTRT1`` uses ``FORMAT(10X, G10.4)``, so only columns 11-20 hold the
    value.  Lines are consumed in order and matched to ``fields`` by position;
    the label text is captured but never used to decide anything, because the
    Fortran does not use it either.
    """
    lines = _data_lines(path)
    if len(lines) < len(fields):
        raise ValueError(
            f"{path}: expected at least {len(fields)} data lines, found {len(lines)}"
        )

    values: dict[str, float] = {}
    labels: list[str] = []
    for name, line in zip(fields, lines):
        padded = line.ljust(20)
        labels.append(padded[:10].strip())
        try:
            value = float(padded[10:20])
        except ValueError as exc:
            raise ValueError(
                f"{path}: could not read {name!r} from columns 11-20 of {line!r}"
            ) from exc
        values[name] = round(value) if name in _INTEGER_FIELDS else value

    return NumericalParameters(values=values, labels=labels, path=Path(path))


def read_er1(path: str | Path) -> tuple[DriverParameters, NumericalParameters]:
    """Read ``.ER1`` and build the driver parameters from it.

    Returns both the typed parameter object and the raw file contents, since a
    few consumers outside the source model (``ALPH``, ``SZF``, ``CRFG``) need
    fields the driver does not carry.
    """
    raw = read_parameters(path, ER1_FIELDS)
    v = raw.values
    params = DriverParameters(
        # blanket closure
        ce=v["ce"],
        epsilon=v["epsilon"],
        delrmn=v["delrmn"],
        dellay=v["dellay"],
        iphifl=int(v["iphifl"]),
        srccut=v["srccut"],
        srcss=v["srcss"],
        srcoer=v["srcoer"],
        vua=v["vua"],
        vub=v["vub"],
        vue=v["vue"],
        vud=v["vud"],
        # integration
        stpin=v["stpin"],
        erbnd=v["erbnd"],
        wtrg=v["wtrg"],
        wttm=v["wttm"],
        wtyc=v["wtyc"],
        wtya=v["wtya"],
        wteb=v["wteb"],
        wtmb=v["wtmb"],
        ernobl=v["ernobl"],
        noblpt=int(v["noblpt"]),
    )
    return params, raw


def read_er2(path: str | Path) -> NumericalParameters:
    """Read ``.ER2``: the tolerances and output criteria for ``DEG2``/``DEG2S``.

    Unlike ``.ER1`` this one is consumed directly rather than folded into a
    typed object, because its fields belong to three different consumers (the
    dense-phase integrator, the Gaussian integrator, and the transient
    observer machinery) that are configured separately.
    """
    return read_parameters(path, ER2_FIELDS)


def alph_settings(raw: NumericalParameters) -> dict:
    """Keyword arguments for :func:`degali.core.atmosphere.fit_alpha`."""
    v = raw.values
    return dict(
        zlow=v["zlow"], xli=v["xli"], xri=v["xri"], eps=v["eps"],
        stpinz=v["stpinz"], erbndz=v["erbndz"],
        ialpfl=int(v["ialpfl"]), alpco=v["alpco"],
    )
