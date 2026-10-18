from __future__ import annotations

import ast
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from tools.score_ti021_dense_layer_hybrid import _handoff_vertical_geometry


SCRIPT = Path(__file__).resolve().parents[1] / "tools/score_ti021_dense_layer_hybrid.py"


def test_finite_prism_multilayer_initialization_uses_current_handoff_height() -> None:
    state = SimpleNamespace(height_m=6.6)

    gaussian, domain_top = _handoff_vertical_geometry(
        "finite_prism", state, puff=None
    )

    assert gaussian == pytest.approx((3.3, 6.6 / math.sqrt(12.0)))
    assert domain_top == pytest.approx(6.6)


def test_gaussian_multilayer_initialization_uses_current_handoff_moments() -> None:
    state = SimpleNamespace(
        centre_position_m=(1.0, 2.0, 4.0), vertical_half_width_m=3.0
    )
    puff = SimpleNamespace(config=SimpleNamespace(width_cutoff_sigma=1.5))

    gaussian, domain_top = _handoff_vertical_geometry("gaussian", state, puff)

    assert gaussian == pytest.approx((4.0, 2.0))
    assert domain_top is None


def test_cli_forwards_safety_and_anisotropic_physics_options() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "run"
    ]
    assert len(calls) == 1
    forwarded = {keyword.arg for keyword in calls[0].keywords}
    assert {
        "puff_top_ambient_turbulence_coefficient",
        "puff_side_ambient_turbulence_coefficient",
        "horizontal_diffusion_respect_walls",
        "cavity_core_entrainment_reduction",
        "cavity_reservoir_capture_mode",
        "cavity_reservoir_vertical_capture_mode",
        "multilayer_entrainment_placement",
    } <= forwarded


def test_cli_has_no_hidden_process_global_option_channel() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert "_CLI_" not in source
