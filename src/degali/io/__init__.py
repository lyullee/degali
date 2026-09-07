"""Reading and writing DEGADIS input and transfer files.

Every DEGADIS file is free-format and positional: values carry no labels and
are identified purely by the order the Fortran reads them in. So each reader
here is a transcription of one ``READ`` sequence.

=========  =================================================================
``.INP``   ground-level input deck (``IO.FOR``)
``.ER1``   source-model tolerances and closure constants (``ESTRT1``)
``.ER2``   downwind and observer tolerances (``ESTRT2``)
``.TR2``   the DEG1 to DEG2 handoff (``TRANS``/``STRT2``)
``.INO``   jet/plume input deck (``JETPLUIN``)
``.IND``   jet touchdown, the JETPLU to DEGBRIDG handoff
=========  =================================================================
"""

from .inp import Case, SourceTable, read_inp, write_inp
from .jetdeck import JetDeck, Touchdown, read_ind, read_ino, write_ind
from .params import (
    ER1_FIELDS,
    ER2_FIELDS,
    NumericalParameters,
    alph_settings,
    read_er1,
    read_er2,
    read_parameters,
)
from .tr2 import Handoff, read_tr2

__all__ = [
    "Case", "SourceTable", "read_inp", "write_inp",
    "NumericalParameters", "ER1_FIELDS", "ER2_FIELDS",
    "read_er1", "read_er2", "read_parameters", "alph_settings",
    "Handoff", "read_tr2",
    "JetDeck", "Touchdown", "read_ino", "read_ind", "write_ind",
]
