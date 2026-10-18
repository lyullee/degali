import json
from pathlib import Path

from tools.audit_bos_obstacle_benchmark import audit


ROOT = Path(__file__).resolve().parents[1]


def test_bos_benchmark_is_fail_closed_and_shape_registered():
    descriptor = ROOT / "reference" / "obstacle_wake" / "bos_obstacle_concentration_benchmark.json"
    result = audit(descriptor)
    assert result["status"] == "withheld_literature_aggregate_only"
    assert result["promotion_allowed"] is False
    assert result["reported_inversion_error_percent"] == 3.81
    assert result["shape_signatures"]["flat_plate_with_sidewalls"] == "W_radial_profile"
