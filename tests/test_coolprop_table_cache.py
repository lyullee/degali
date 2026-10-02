import pytest


pytest.importorskip("CoolProp")

from degali.core.thermo import CoolPropBackend


@pytest.fixture(autouse=True)
def isolated_shared_tables(monkeypatch):
    CoolPropBackend.clear_shared_tables()
    monkeypatch.setattr(CoolPropBackend, "POINTS", 11)
    yield
    CoolPropBackend.clear_shared_tables()


def test_compatible_backends_share_one_immutable_property_table():
    calls = []

    def properties(_output, _input1, temperature, _input2, _pressure, _fluid):
        calls.append(float(temperature))
        return 1000.0 + float(temperature)

    first = CoolPropBackend("Hydrogen", span=(100.0, 200.0))
    second = CoolPropBackend("Hydrogen", span=(100.0, 200.0))
    first._props = second._props = properties

    assert first.cp_air(150.0) == pytest.approx(1150.0)
    assert len(calls) == CoolPropBackend.POINTS

    assert second.cp_air(160.0) == pytest.approx(1160.0)
    assert len(calls) == CoolPropBackend.POINTS
    assert second._grids["cpa"] is first._grids["cpa"]
    assert not second._grids["cpa"][2].flags.writeable
    assert not second._grids["cpa"][3].flags.writeable
    assert CoolPropBackend.shared_table_count() == 1


def test_contaminant_tables_are_separated_by_pressure():
    calls = []

    def properties(_output, _input1, temperature, _input2, pressure, _fluid):
        calls.append((float(temperature), float(pressure)))
        return float(temperature) + float(pressure) * 1.0e-6

    backend = CoolPropBackend(
        "Hydrogen", span=(100.0, 200.0), force_contaminant_gas=True
    )
    backend._props = properties

    at_one_atmosphere = backend.cp_contaminant_eos(150.0, 1.0)
    assert len(calls) == CoolPropBackend.POINTS

    at_two_atmospheres = backend.cp_contaminant_eos(150.0, 2.0)
    assert len(calls) == 2 * CoolPropBackend.POINTS
    assert at_two_atmospheres > at_one_atmosphere
    assert CoolPropBackend.shared_table_count() == 2
