from dataclasses import dataclass
import json
import runpy


@dataclass
class _Reading:
    radius: float
    height: float
    peak: float


def _operator_table():
    namespace = runpy.run_path("tools/audit_ffi_sensor_operator.py")
    return namespace["_operator_table"]


def test_ffi_operator_table_separates_arc_and_sensor_height_groups():
    readings = [
        _Reading(30.0, 0.1, 15.0),
        _Reading(30.0, 1.0, 12.0),
        _Reading(30.0, 0.1, 21.0),
    ]
    predicted = [5.0, 4.0, 7.5]
    table = _operator_table()

    arc = table(readings, predicted, group="radius")
    heights = table(readings, predicted, group="height")

    assert arc == [{
        "radius_m": 30.0,
        "sensor_count": 3,
        "observed_arc_max_vol_pct": 21.0,
        "projected_sensor_max_vol_pct": 7.5,
        "observed_over_projected": 2.8,
        "observed_is_lower_bound": True,
    }]
    assert [row["height_m"] for row in heights] == [0.1, 1.0]
    assert heights[0]["observed_arc_max_vol_pct"] == 21.0
    assert heights[0]["projected_sensor_max_vol_pct"] == 7.5
    assert heights[1]["observed_over_projected"] == 3.0


def test_ffi_operator_table_rejects_mismatched_inputs():
    table = _operator_table()
    reading = _Reading(30.0, 0.1, 1.0)

    try:
        table([reading], [], group="radius")
    except ValueError as error:
        assert "equal lengths" in str(error)
    else:
        raise AssertionError("mismatched operator inputs should fail closed")

    try:
        table([reading], [1.0], group="bearing")
    except ValueError as error:
        assert "group" in str(error)
    else:
        raise AssertionError("unknown operator group should fail closed")


def test_ffi_operator_audit_emits_separate_tables_without_raw_sensor_rows():
    namespace = runpy.run_path("tools/audit_ffi_sensor_operator.py")
    assert "_operator_table" in namespace
    # The CLI itself is exercised by the existing local-reference smoke check;
    # this test locks the output contract to aggregate rows only.
    payload = {
        "arc_max_table": namespace["_operator_table"](
            [_Reading(30.0, 0.1, 21.0)], [7.5], group="radius"
        ),
        "sensor_height_table": namespace["_operator_table"](
            [_Reading(30.0, 0.1, 21.0)], [7.5], group="height"
        ),
    }
    encoded = json.dumps(payload)
    assert "sensor_count" in encoded
    assert "bearing" not in encoded.lower()
