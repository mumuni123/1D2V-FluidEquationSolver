#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import math
import os
import re
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import config as cfg
import normalization as norm


def _resolve_default_paths(script_dir: str) -> Tuple[str, str]:
    project_root = os.path.dirname(script_dir)
    output_dir = os.path.join(project_root, "output")
    results_dir = os.path.join(project_root, "results")
    return output_dir, results_dir


def _load_csv(path: str) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError("empty file: {0}".format(path))
    return data


def _extract_step(path: str) -> int:
    m = re.search(r"state_(\d+)\.csv$", os.path.basename(path))
    if m is None:
        return -1
    return int(m.group(1))


def _list_state_files(output_dir: str) -> List[str]:
    files = glob.glob(os.path.join(output_dir, "state_*.csv"))
    files = [f for f in files if _extract_step(f) >= 0]
    files = sorted(files, key=_extract_step)
    if not files:
        raise FileNotFoundError("No state_*.csv found in: {0}".format(output_dir))
    return files


def _parse_positions_um(positions_text: str) -> List[float]:
    out: List[float] = []
    for token in positions_text.split(","):
        s = token.strip()
        if not s:
            continue
        out.append(float(s))
    if not out:
        raise ValueError("No valid positions provided")
    return out


def _nearest_indices(z_um: np.ndarray, pos_um: Sequence[float]) -> Tuple[List[int], List[float]]:
    idx: List[int] = []
    actual: List[float] = []
    for p in pos_um:
        i = int(np.argmin(np.abs(z_um - p)))
        idx.append(i)
        actual.append(float(z_um[i]))
    return idx, actual


def _build_time_axis(output_dir: str, files: Sequence[str], n0: float) -> Tuple[np.ndarray, str]:
    diag_path = os.path.join(output_dir, "diagnostics.csv")
    if os.path.isfile(diag_path):
        diag = _load_csv(diag_path)
        if diag.size == len(files):
            try:
                return norm.diagnostics_time_axis(diag, False, n0)
            except KeyError:
                pass
        else:
            try:
                diag_t, diag_label = norm.diagnostics_time_axis(diag, False, n0)
            except KeyError:
                diag_t = np.asarray([], dtype=np.float64)
                diag_label = ""

            diag_t = np.asarray(diag_t, dtype=np.float64)
            diag_t = diag_t[np.isfinite(diag_t)]
            if diag_t.size >= 2:
                diag_dt = np.diff(diag_t)
                median_diag_dt = float(np.median(diag_dt))
                steps = np.asarray([_extract_step(f) for f in files], dtype=np.float64)
                unique_steps = np.unique(steps[steps >= 0])
                if (
                    median_diag_dt > 0.0
                    and np.isfinite(median_diag_dt)
                    and unique_steps.size >= 2
                    and diag_label
                ):
                    step_stride = float(np.median(np.diff(unique_steps)))
                    if step_stride > 0.0 and np.isfinite(step_stride):
                        return (steps - steps[0]) / step_stride * median_diag_dt + diag_t[0], diag_label

    steps = np.asarray([_extract_step(f) for f in files], dtype=np.float64)
    if np.all(steps >= 0):
        return steps, "step"

    return np.arange(len(files), dtype=np.float64), "snapshot index"


def _safe_line(ax: plt.Axes, x: np.ndarray, y: np.ndarray, label: str) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) == 0:
        return
    ax.plot(x[mask], y[mask], label=label)


def _robust_ylim(ax: plt.Axes, y_all: np.ndarray) -> None:
    arr = y_all[np.isfinite(y_all)]
    if arr.size < 10:
        return
    lo = np.percentile(arr, 1.0)
    hi = np.percentile(arr, 99.0)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return
    pad = 0.08 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)


def _var_ylabel(var_name: str, normalized: bool) -> str:
    if normalized:
        ratios = {
            "Ex": "Ex/E0",
            "Ez": "Ez/E0",
            "By": "By/B0",
            "ne": "ne/n0",
            "vx": "vx/c",
            "vz": "vz/c",
            "Pe": "Pe/P0",
        }
        return ratios.get(var_name, var_name)

    units = {
        "Ex": "V/m",
        "Ez": "V/m",
        "By": "T",
        "ne": "m^-3",
        "vx": "m/s",
        "vz": "m/s",
        "Pe": "Pa",
    }
    unit = units.get(var_name, "")
    if unit:
        return "{0} ({1})".format(var_name, unit)
    return var_name


def _nearest_value_at_position(
    state_phys: Dict[str, np.ndarray],
    state_values: Dict[str, np.ndarray],
    var_name: str,
    position_um: float,
) -> Tuple[float, float]:
    z_um = np.asarray(state_phys["z"], dtype=np.float64) * 1.0e6
    values = np.asarray(state_values[var_name], dtype=np.float64)
    mask = np.isfinite(z_um) & np.isfinite(values)
    if np.count_nonzero(mask) == 0:
        return np.nan, np.nan

    z_valid = z_um[mask]
    v_valid = values[mask]
    i = int(np.argmin(np.abs(z_valid - position_um)))
    return float(v_valid[i]), float(z_valid[i])


def run(
    output_dir: str,
    results_dir: str,
    positions_um: Sequence[float],
    variables: Sequence[str],
    normalize: bool,
    n0: float,
) -> str:
    if not os.path.isdir(output_dir):
        raise FileNotFoundError("Output directory not found: {0}".format(output_dir))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    files = _list_state_files(output_dir)
    first = _load_csv(files[0])

    _, is_norm = norm.convert_state_columns(first, normalize, n0)
    s0_phys, _ = norm.convert_state_columns(first, False, n0)
    z_um = s0_phys["z"] * 1.0e6
    _, actual_pos = _nearest_indices(z_um, positions_um)

    time_x, x_label = _build_time_axis(output_dir, files, n0)

    for v in variables:
        if v not in ["Ex", "Ez", "By", "ne", "vx", "vz", "Pe"]:
            raise KeyError("Variable '{0}' is not supported".format(v))

    series: Dict[str, np.ndarray] = {}
    for v in variables:
        series[v] = np.full((len(positions_um), len(files)), np.nan, dtype=np.float64)

    for j, f in enumerate(files):
        s = _load_csv(f)
        sv, _ = norm.convert_state_columns(s, is_norm, n0, cfg.VELOCITY_DENSITY_MIN_RATIO)
        s_phys, _ = norm.convert_state_columns(s, False, n0)
        for vi, pos_um in enumerate(positions_um):
            for v in variables:
                value, actual = _nearest_value_at_position(s_phys, sv, v, pos_um)
                series[v][vi, j] = value
                if np.isfinite(actual):
                    actual_pos[vi] = actual

    n_var = len(variables)
    n_cols = 2
    n_rows = int(math.ceil(float(n_var) / float(n_cols)))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(13, 4.2 * n_rows))
    axes_flat = np.atleast_1d(axes).ravel()

    for k, v in enumerate(variables):
        ax = axes_flat[k]
        y_collect: List[np.ndarray] = []
        for vi, pos in enumerate(actual_pos):
            y = series[v][vi, :]
            y_collect.append(y)
            _safe_line(ax, time_x, y, "z={0:.3f} um".format(pos))
        ax.set_title("{0} at selected positions".format(v))
        ax.set_xlabel(x_label)
        ax.set_ylabel(_var_ylabel(v, is_norm))
        ax.grid(True)
        ax.legend(loc="best", fontsize=8)
        _robust_ylim(ax, np.concatenate(y_collect))

        # Always show the full simulated time range (e.g. up to 120 fs),
        # even if data after some point is NaN and therefore not drawable.
        finite_t = time_x[np.isfinite(time_x)]
        if finite_t.size > 1:
            ax.set_xlim(float(np.min(finite_t)), float(np.max(finite_t)))

    for k in range(n_var, n_rows * n_cols):
        axes_flat[k].axis("off")

    fig.tight_layout()
    out_path = os.path.join(results_dir, "time_evolution_at_positions.png")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    print("Selected positions (requested -> actual):")
    for req, act in zip(positions_um, actual_pos):
        print("- {0:.4f} um -> {1:.4f} um".format(req, act))

    return out_path


def parse_args() -> argparse.Namespace:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output, default_results = _resolve_default_paths(script_dir)

    parser = argparse.ArgumentParser(
        description="Plot full-time evolution at user-specified positions."
    )
    parser.add_argument("--output", default=cfg.OUTPUT_DIR or default_output, help="Path to output directory")
    parser.add_argument("--results", default=cfg.RESULTS_DIR or default_results, help="Path to results directory")
    parser.add_argument(
        "--positions-um",
        default=cfg.positions_to_text(cfg.TIME_POSITIONS_UM),
        help="Comma-separated positions in um, e.g. 0.5,2.0,4.0",
    )
    parser.add_argument(
        "--variables",
        default=cfg.variables_to_text(cfg.TIME_VARIABLES),
        help="Comma-separated variables to plot",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Enable normalized plotting variables",
    )
    parser.add_argument(
        "--n0",
        type=float,
        default=cfg.INITIAL_DENSITY_M3,
        help="Initial density n0 (m^-3) used for normalization references",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    positions_um = _parse_positions_um(args.positions_um)
    variables = [v.strip() for v in args.variables.split(",") if v.strip()]
    normalize = cfg.NORMALIZE_OUTPUT or args.normalize

    out_path = run(args.output, args.results, positions_um, variables, normalize, args.n0)
    print("Saved: {0}".format(out_path))


if __name__ == "__main__":
    main()
