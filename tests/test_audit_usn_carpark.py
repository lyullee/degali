import json
from pathlib import Path
import runpy

import pytest


ROOT = Path("reference/datasets/usn-carpark-14405903/raw")


@pytest.mark.skipif(
    not (ROOT / "Exp19.mat").is_file(),
    reason="requires the locally controlled public USN car-park data",
)
def test_public_usn_carpark_audit_reproduces_paper_table3(tmp_path):
    tool = runpy.run_path("tools/audit_usn_carpark.py")
    record = tool["audit_usn_carpark"](ROOT)

    assert record["status"] == "component_evidence_ready"
    assert record["reproduction"]["paper_table3_operator_reproduced"]
    assert record["reproduction"]["maximum_first_arrival_difference_s"] < 0.5
    assert record["reproduction"]["maximum_flammable_duration_difference_s"] < 6.0
    assert record["promotion_allowed"] is False
    assert all(event["published_md5_verified"] for event in record["events"])

    output = tmp_path / "usn-audit.json"
    tool["write_audit"](record, output)
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded["events"][0]["threshold_operator"]["sensor_count"] == 24
    with pytest.raises(FileExistsError):
        tool["write_audit"](record, output)
