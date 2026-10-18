from tools.compare_thorney_ti021_geometry_representation_prefixes import _scope_warning


def test_scope_warning_reports_actual_comparison_window():
    warning = _scope_warning(18.0)

    assert "same 18-s computed window" in warning
    assert "first-second" not in warning
    assert "not an experimental validation" in warning
