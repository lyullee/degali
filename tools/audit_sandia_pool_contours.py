"""Screen local Sandia LH2 pool contour reductions without redistributing data."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.validation.sandia_pool import (
    read_pool_contour_lower_bounds,
    screen_pool_contours,
    screen_report,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("observations", type=Path,
                        help="ignored local JSON contour reduction")
    parser.add_argument("predictions", type=Path,
                        help="local JSON object: contour identifier -> centreline reach [m]")
    parser.add_argument("--output", type=Path, help="optional local JSON report")
    args = parser.parse_args()
    observations = read_pool_contour_lower_bounds(args.observations)
    predictions = json.loads(args.predictions.read_text(encoding="utf-8"))
    if not isinstance(predictions, dict):
        raise ValueError("predictions must be a JSON object")
    report = screen_report(screen_pool_contours(observations, predictions))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
