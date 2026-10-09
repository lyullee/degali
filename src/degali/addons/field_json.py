"""Strict JSON helpers for auditable field artifacts.

JSON object keys are part of the input contract. A normal ``json.loads``
silently keeps the last duplicate key, which can change a schema field,
fingerprint, or applicability flag without leaving a parse error. Field
artifact readers therefore use this module before mapping/schema validation.
"""

from __future__ import annotations

import json
from typing import Any


def _strict_object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    values: dict[str, Any] = {}
    normalized: dict[str, str] = {}
    for key, value in pairs:
        if not isinstance(key, str):
            raise ValueError("JSON object keys must be strings")
        normalized_key = key.strip().casefold()
        previous = normalized.get(normalized_key)
        if previous is not None:
            raise ValueError(
                "duplicate or normalization-colliding JSON object key: "
                f"{previous!r} and {key!r}"
            )
        normalized[normalized_key] = key
        values[key] = value
    return values


def strict_json_loads(value: str) -> Any:
    """Decode JSON while refusing duplicate or whitespace/case-colliding keys."""
    return json.loads(value, object_pairs_hook=_strict_object_pairs)


__all__ = ["strict_json_loads"]
