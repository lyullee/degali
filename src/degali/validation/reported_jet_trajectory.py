"""Opt-in receptor adapter for JetPlume's already-imaged reported rows.

Not for independent-energy states or touchdown output. The historical
Trajectory remains unchanged so frozen field results are reproducible.
"""
import numpy as np

from .nearfield import Trajectory


class ReportedJetTrajectory(Trajectory):
    """Recover the bare Gaussian at stored nodes, then image exactly once."""

    def __init__(self, table, reported_rows):
        reported = np.asarray(reported_rows, dtype=float)
        if reported.ndim != 2 or reported.shape[1] != 12 or len(reported) < 3:
            raise ValueError('expected at least three 12-column JetPlume reported rows')
        if not np.all(np.isfinite(reported)):
            raise ValueError('nonfinite reported rows')
        if np.any(reported[:, 1] <= 0):
            raise ValueError('touchdown/nonpositive elevation requires a separate adapter')
        if np.any(reported[:, 2] < 0) or np.any(reported[:, 3:5] <= 0):
            raise ValueError('invalid concentration or Gaussian widths')
        if np.any(np.diff(reported[:, 0]) <= 0):
            raise ValueError('reported downwind positions must be strictly increasing')
        self.reported_rows = reported.copy()
        bare = reported[:, :5].copy()
        bare[:, 2] /= 1.0 + np.exp(-2.0 * (reported[:, 1]/reported[:, 4])**2)
        super().__init__(table, bare)

    def _require_receptor(self, x, y, z):
        if not np.all(np.isfinite([x, y, z])) or z < 0:
            raise ValueError('receptor must be finite and above/on the ground')
        if self.at(x) is None:
            raise ValueError('receptor is outside the integrated domain')

    def concentration_at(self, x, y, z):
        self._require_receptor(x, y, z)
        return super().concentration_at(x, y, z)

    def temperature_at(self, x, y, z):
        self._require_receptor(x, y, z)
        return super().temperature_at(x, y, z)
