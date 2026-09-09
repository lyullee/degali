"""Confirm a smooth phase lookup against the existing exact property oracle."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator, RectBivariateSpline, RegularGridInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_reservoir_thermal_segments import actual_source_model
from degali.addons.enthalpy_profile import PhaseMassEnthalpyInverter


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


class HermiteSurface:
    def __init__(self, x, y, values):
        self.x, self.y = np.asarray(x), np.asarray(y)
        values = np.asarray(values)
        fx = PchipInterpolator(self.x, values, axis=0).derivative()(self.x)
        fy = PchipInterpolator(self.y, values, axis=1).derivative()(self.y)
        fxy_x = PchipInterpolator(self.y, fx, axis=1).derivative()(self.y)
        fxy_y = PchipInterpolator(self.x, fy, axis=0).derivative()(self.x)
        fxy = .5*(fxy_x+fxy_y)
        dx, dy = np.diff(self.x)[:, None], np.diff(self.y)[None, :]
        p = np.empty((len(self.x)-1, len(self.y)-1, 4, 4))
        p[:, :, 0, 0], p[:, :, 0, 1] = values[:-1, :-1], values[:-1, 1:]
        p[:, :, 1, 0], p[:, :, 1, 1] = values[1:, :-1], values[1:, 1:]
        p[:, :, 0, 2], p[:, :, 0, 3] = fy[:-1, :-1]*dy, fy[:-1, 1:]*dy
        p[:, :, 1, 2], p[:, :, 1, 3] = fy[1:, :-1]*dy, fy[1:, 1:]*dy
        p[:, :, 2, 0], p[:, :, 2, 1] = fx[:-1, :-1]*dx, fx[:-1, 1:]*dx
        p[:, :, 3, 0], p[:, :, 3, 1] = fx[1:, :-1]*dx, fx[1:, 1:]*dx
        scale = dx*dy
        p[:, :, 2, 2], p[:, :, 2, 3] = fxy[:-1, :-1]*scale, fxy[:-1, 1:]*scale
        p[:, :, 3, 2], p[:, :, 3, 3] = fxy[1:, :-1]*scale, fxy[1:, 1:]*scale
        matrix = np.array([
            [1., 0., 0., 0.], [0., 0., 1., 0.],
            [-3., 3., -2., -1.], [2., -2., 1., 1.],
        ])
        self.coefficients = np.einsum("ab,ijbc,cd->ijad", matrix, p, matrix.T)

    def evaluate(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        shape = x.shape
        x, y = x.ravel(), y.ravel()
        i = np.clip(np.searchsorted(self.x, x, side="right")-1, 0, len(self.x)-2)
        j = np.clip(np.searchsorted(self.y, y, side="right")-1, 0, len(self.y)-2)
        dx, dy = self.x[i+1]-self.x[i], self.y[j+1]-self.y[j]
        tx, ty = (x-self.x[i])/dx, (y-self.y[j])/dy
        bx = np.column_stack([np.ones_like(tx), tx, tx*tx, tx*tx*tx])
        by = np.column_stack([np.ones_like(ty), ty, ty*ty, ty*ty*ty])
        dbx = np.column_stack([np.zeros_like(tx), np.ones_like(tx), 2*tx, 3*tx*tx])/dx[:, None]
        dby = np.column_stack([np.zeros_like(ty), np.ones_like(ty), 2*ty, 3*ty*ty])/dy[:, None]
        coefficients = self.coefficients[i, j]
        value = np.einsum("ni,nij,nj->n", bx, coefficients, by)
        derivative_x = np.einsum("ni,nij,nj->n", dbx, coefficients, by)
        derivative_y = np.einsum("ni,nij,nj->n", bx, coefficients, dby)
        return value.reshape(shape), derivative_x.reshape(shape), derivative_y.reshape(shape)


def errors(actual, estimate, floor):
    scaled = abs(estimate-actual)/np.maximum(abs(actual), floor)
    return {
        "maximum": float(max(scaled)),
        "p99": float(np.quantile(scaled, .99)),
        "rms": float(np.sqrt(np.mean(scaled*scaled))),
    }


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite smooth phase lookup evidence")
    dependencies = [
        ROOT / "docs/prereg-smooth-phase-hermite-lookup.md",
        ROOT / "src/degali/addons/axisymmetric_jet.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    reduced = next(
        row for row in json.loads((REF / "e35_reduced.json").read_text(encoding="utf-8"))["trials"]
        if row["trial"] == 24
    )
    measured = {
        row["trial"]: row
        for row in json.loads((
            REF / "measured_pipe_source_2026-09-05.json"
        ).read_text(encoding="utf-8"))["trials"]
    }
    _, thermodynamics, _ = actual_source_model(reduced, measured)
    inverse = PhaseMassEnthalpyInverter(thermodynamics)
    ideal_grid, fraction_grid, temperature_linear, enthalpy_linear = (
        thermodynamics._condensed_lookup
    )
    temperature_values = temperature_linear.values
    enthalpy_values = enthalpy_linear.values
    temperature_spline = RectBivariateSpline(
        ideal_grid, fraction_grid, temperature_values, kx=3, ky=3, s=0.,
    )
    enthalpy_spline = RectBivariateSpline(
        ideal_grid, fraction_grid, enthalpy_values, kx=3, ky=3, s=0.,
    )
    temperature_pchip = RegularGridInterpolator(
        (ideal_grid, fraction_grid), temperature_values,
        method="pchip", bounds_error=True,
    )
    enthalpy_pchip = RegularGridInterpolator(
        (ideal_grid, fraction_grid), enthalpy_values,
        method="pchip", bounds_error=True,
    )
    temperature_hermite = HermiteSurface(
        ideal_grid, fraction_grid, temperature_values,
    )
    enthalpy_hermite = HermiteSurface(
        ideal_grid, fraction_grid, enthalpy_values,
    )

    rng = np.random.default_rng(240911)
    count = 8000
    ideal = rng.uniform(ideal_grid[0], ideal_grid[-1], count)
    fraction = rng.uniform(fraction_grid[0], fraction_grid[-1], count)
    molecular_weight = 1. / (
        fraction/thermodynamics.fuel_molecular_weight
        +(1.-fraction)/thermodynamics._humid_ambient_molecular_weight
    )
    density = thermodynamics.ambient_pressure*molecular_weight/(8.31446261815324*ideal)
    exact_temperature, exact_enthalpy = thermodynamics._condensed_air_state_exact(
        density, fraction,
    )
    points = np.column_stack([ideal, fraction])
    linear_temperature = temperature_linear(points)
    linear_enthalpy = enthalpy_linear(points)
    spline_temperature = temperature_spline.ev(ideal, fraction)
    spline_enthalpy = enthalpy_spline.ev(ideal, fraction)
    pchip_temperature = temperature_pchip(points)
    pchip_enthalpy = enthalpy_pchip(points)
    hermite_temperature = temperature_hermite.evaluate(ideal, fraction)[0]
    hermite_enthalpy = enthalpy_hermite.evaluate(ideal, fraction)[0]

    dense_ideal = np.linspace(ideal_grid[0], ideal_grid[-1], 801)
    dense_fraction = np.linspace(fraction_grid[0], fraction_grid[-1], 641)
    mesh_t, mesh_y = np.meshgrid(dense_ideal, dense_fraction, indexing="ij")
    flat_t, flat_y = mesh_t.ravel(), mesh_y.ravel()
    h_t = enthalpy_spline.ev(flat_t, flat_y, dx=1, dy=0)
    h_y = enthalpy_spline.ev(flat_t, flat_y, dx=0, dy=1)
    rho = thermodynamics.ambient_pressure*(1./(
        flat_y/thermodynamics.fuel_molecular_weight
        +(1.-flat_y)/thermodynamics._humid_ambient_molecular_weight
    ))/(8.31446261815324*flat_t)
    c = rho*flat_y
    slope = -h_t*flat_t/(rho+inverse.k*c)-h_y*flat_y/rho

    monotone_count = 30000
    check_t = rng.uniform(ideal_grid[2], ideal_grid[-3], monotone_count)
    check_y = rng.uniform(fraction_grid[2], fraction_grid[-3], monotone_count)
    check_mw = 1. / (
        check_y/thermodynamics.fuel_molecular_weight
        +(1.-check_y)/thermodynamics._humid_ambient_molecular_weight
    )
    check_rho = thermodynamics.ambient_pressure*check_mw/(8.31446261815324*check_t)
    check_c = check_rho*check_y
    epsilon = 1e-5
    rho_plus, rho_minus = check_rho*np.exp(epsilon), check_rho*np.exp(-epsilon)
    y_plus, y_minus = check_c/rho_plus, check_c/rho_minus
    t_plus = inverse.a/(rho_plus+inverse.k*check_c)
    t_minus = inverse.a/(rho_minus+inverse.k*check_c)
    pchip_plus = enthalpy_pchip(np.column_stack([t_plus, y_plus]))
    pchip_minus = enthalpy_pchip(np.column_stack([t_minus, y_minus]))
    pchip_slope = (pchip_plus-pchip_minus)/(rho_plus-rho_minus)
    _, hermite_h_t, hermite_h_y = enthalpy_hermite.evaluate(flat_t, flat_y)
    hermite_slope = (
        -hermite_h_t*flat_t/(rho+inverse.k*c)-hermite_h_y*flat_y/rho
    )

    node_t, node_y = np.meshgrid(ideal_grid, fraction_grid, indexing="ij")
    output = {
        "completed": False,
        "random_seed": 240911,
        "samples": count,
        "temperature_linear_absolute_K": errors(
            exact_temperature, linear_temperature, 1.,
        ),
        "temperature_spline_absolute_K": errors(
            exact_temperature, spline_temperature, 1.,
        ),
        "enthalpy_linear_relative": errors(exact_enthalpy, linear_enthalpy, 1.),
        "enthalpy_spline_relative": errors(exact_enthalpy, spline_enthalpy, 1.),
        "temperature_pchip_absolute_K": errors(
            exact_temperature, pchip_temperature, 1.,
        ),
        "enthalpy_pchip_relative": errors(exact_enthalpy, pchip_enthalpy, 1.),
        "temperature_hermite_absolute_K": errors(
            exact_temperature, hermite_temperature, 1.,
        ),
        "enthalpy_hermite_relative": errors(exact_enthalpy, hermite_enthalpy, 1.),
        "temperature_node_max_abs_K": float(max(abs(
            temperature_spline.ev(node_t.ravel(), node_y.ravel())
            -temperature_values.ravel()
        ))),
        "enthalpy_node_max_abs": float(max(abs(
            enthalpy_spline.ev(node_t.ravel(), node_y.ravel())
            -enthalpy_values.ravel()
        ))),
        "temperature_pchip_node_max_abs_K": float(max(abs(
            temperature_pchip(np.column_stack([node_t.ravel(), node_y.ravel()]))
            -temperature_values.ravel()
        ))),
        "enthalpy_pchip_node_max_abs": float(max(abs(
            enthalpy_pchip(np.column_stack([node_t.ravel(), node_y.ravel()]))
            -enthalpy_values.ravel()
        ))),
        "temperature_hermite_node_max_abs_K": float(max(abs(
            temperature_hermite.evaluate(node_t.ravel(), node_y.ravel())[0]
            -temperature_values.ravel()
        ))),
        "enthalpy_hermite_node_max_abs": float(max(abs(
            enthalpy_hermite.evaluate(node_t.ravel(), node_y.ravel())[0]
            -enthalpy_values.ravel()
        ))),
        "spline_density_slope_min": float(min(slope)),
        "spline_density_slope_max": float(max(slope)),
        "spline_nonnegative_density_slope_points": int(np.count_nonzero(slope >= 0.)),
        "spline_dense_points": int(len(slope)),
        "pchip_density_slope_min": float(min(pchip_slope)),
        "pchip_density_slope_max": float(max(pchip_slope)),
        "pchip_nonnegative_density_slope_points": int(np.count_nonzero(pchip_slope >= 0.)),
        "pchip_slope_points": monotone_count,
        "hermite_density_slope_min": float(min(hermite_slope)),
        "hermite_density_slope_max": float(max(hermite_slope)),
        "hermite_nonnegative_density_slope_points": int(np.count_nonzero(hermite_slope >= 0.)),
        "hermite_slope_points": int(len(hermite_slope)),
    }
    output["confirmation_passed"] = bool(
        output["temperature_hermite_absolute_K"]["rms"]
            < output["temperature_linear_absolute_K"]["rms"]
        and output["temperature_hermite_absolute_K"]["p99"]
            < output["temperature_linear_absolute_K"]["p99"]
        and output["enthalpy_hermite_relative"]["rms"]
            < output["enthalpy_linear_relative"]["rms"]
        and output["enthalpy_hermite_relative"]["p99"]
            < output["enthalpy_linear_relative"]["p99"]
        and output["temperature_hermite_node_max_abs_K"] <= 1e-10
        and output["enthalpy_hermite_node_max_abs"] <= 1e-6
        and output["hermite_nonnegative_density_slope_points"] == 0
        and np.isfinite(output["hermite_density_slope_min"])
        and np.isfinite(output["hermite_density_slope_max"])
        and output["hermite_density_slope_max"] < 0.
    )
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("a smooth lookup screening input changed during execution")
    output["hashes"] = after
    output["completed"] = True
    output["finished_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "completed": True,
        "confirmation_passed": output["confirmation_passed"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
