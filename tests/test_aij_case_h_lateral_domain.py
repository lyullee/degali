from __future__ import annotations

import json
from pathlib import Path

from tools.compare_aij_case_h_lateral_domain_tracer import compare


ROOT = Path(__file__).resolve().parents[1]


def test_lateral_domain_tracer_comparison_uses_common_point_ids(tmp_path: Path) -> None:
    result = compare(
        ROOT / "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1",
        ROOT / "outputs/aij-case-h-tracer-fv-momentum-sgs-wide05-12s-dx2p5cm-dt1ms-2026-10-14-v1",
        tmp_path,
    )
    assert result["common_point_count"] == 582
    assert result["calibration_used"] is False
    assert result["promotion_allowed"] is False
    assert result["rmse_reduction_fraction"] > 0.20
    record = json.loads(
        (tmp_path / "aij-case-h-lateral-domain-tracer-comparison.json").read_text(
            encoding="utf-8"
        )
    )
    assert record["wide_common_point_metrics"]["pearson_r"] < record["base_common_point_metrics"]["pearson_r"]
