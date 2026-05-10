#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import os
import re
from typing import List, Tuple

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


def _list_state_files(output_dir: str) -> List[str]:
    files = glob.glob(os.path.join(output_dir, "state_*.csv"))
    files = [f for f in files if _extract_step(f) >= 0]
    files = sorted(files, key=_extract_step)
    if not files:
        raise FileNotFoundError("No state_*.csv found in: {0}".format(output_dir))
    return files


def _extract_step(path: str) -> int:
    m = re.search(r"state_(\d+)\.csv$", os.path.basename(path))
    if m is None:
        return -1
    return int(m.group(1))


def _finite_ratio(arr: np.ndarray) -> float:
    x = np.asarray(arr, dtype=np.float64)
    if x.size == 0:
        return 0.0
    return float(np.isfinite(x).sum()) / float(x.size)


def _is_state_usable(data: np.ndarray, min_ratio: float = 0.10) -> bool:
    has_phys = all(k in data.dtype.names for k in ["Ex", "Ez", "By", "ne", "vx", "vz", "Pe"])
    has_tilde = all(
        k in data.dtype.names
        for k in ["Ex_tilde", "Ez_tilde", "By_tilde", "ne_tilde", "vx_tilde", "vz_tilde", "Pe_tilde"]
    )
    if has_phys:
        needed = ["Ex", "Ez", "By", "ne", "vx", "vz", "Pe"]
    elif has_tilde:
        needed = ["Ex_tilde", "Ez_tilde", "By_tilde", "ne_tilde", "vx_tilde", "vz_tilde", "Pe_tilde"]
    else:
        return False

    for key in needed:
        if _finite_ratio(data[key]) < min_ratio:
            return False
    return True


def _find_nearest_valid_state(files: List[str], target_path: str) -> str:
    if target_path not in files:
        return target_path

    idx = files.index(target_path)
    order: List[int] = [idx]
    for d in range(1, len(files)):
        left = idx - d
        right = idx + d
        if left >= 0:
            order.append(left)
        if right < len(files):
            order.append(right)

    for i in order:
        p = files[i]
        try:
            s = _load_csv(p)
        except Exception:
            continue
        if _is_state_usable(s):
            return p

    return target_path


def _resolve_state_path(output_dir: str, state_file: str) -> str:
    if state_file.lower() == "latest":
        return _list_state_files(output_dir)[-1]

    if os.path.isabs(state_file):
        if not os.path.isfile(state_file):
            raise FileNotFoundError("State file not found: {0}".format(state_file))
        return state_file

    cand1 = os.path.join(output_dir, state_file)
    if os.path.isfile(cand1):
        return cand1

    if os.path.isfile(state_file):
        return state_file

    raise FileNotFoundError("State file not found: {0}".format(state_file))


def _safe_line(ax: plt.Axes, x: np.ndarray, y: np.ndarray, label: str = "") -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) == 0:
        return
    if label:
        ax.plot(x[mask], y[mask], label=label)
    else:
        ax.plot(x[mask], y[mask])


def _robust_ylim(ax: plt.Axes, y: np.ndarray) -> None:
    arr = y[np.isfinite(y)]
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
    else:
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


def run(output_dir: str, results_dir: str, state_file: str, normalize: bool, n0: float) -> str:
    if not os.path.isdir(output_dir):
        raise FileNotFoundError("Output directory not found: {0}".format(output_dir))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    state_files = _list_state_files(output_dir)
    state_path = _resolve_state_path(output_dir, state_file)
    state_path = _find_nearest_valid_state(state_files, state_path)
    data = _load_csv(state_path)

    if not _is_state_usable(data):
        raise ValueError(
            "Selected state file is invalid (too many NaN/Inf): {0}".format(state_path)
        )

    s, is_norm = norm.convert_state_columns(data, normalize, n0, cfg.VELOCITY_DENSITY_MIN_RATIO)
    s_phys, _ = norm.convert_state_columns(data, False, n0)
    z_axis = s_phys["z"] * 1.0e6
    x_label = "z (um)"

    fig, axes = plt.subplots(3, 2, figsize=(12, 10))

    _safe_line(axes[0, 0], z_axis, s["Ex"])
    axes[0, 0].set_title("Ex")
    axes[0, 0].set_xlabel(x_label)
    axes[0, 0].set_ylabel(_var_ylabel("Ex", is_norm))
    axes[0, 0].grid(True)
    _robust_ylim(axes[0, 0], s["Ex"])

    _safe_line(axes[0, 1], z_axis, s["Ez"])
    axes[0, 1].set_title("Ez")
    axes[0, 1].set_xlabel(x_label)
    axes[0, 1].set_ylabel(_var_ylabel("Ez", is_norm))
    axes[0, 1].grid(True)
    _robust_ylim(axes[0, 1], s["Ez"])

    _safe_line(axes[1, 0], z_axis, s["By"])
    axes[1, 0].set_title("By")
    axes[1, 0].set_xlabel(x_label)
    axes[1, 0].set_ylabel(_var_ylabel("By", is_norm))
    axes[1, 0].grid(True)
    _robust_ylim(axes[1, 0], s["By"])

    _safe_line(axes[1, 1], z_axis, s["ne"])
    axes[1, 1].set_title("ne")
    axes[1, 1].set_xlabel(x_label)
    axes[1, 1].set_ylabel(_var_ylabel("ne", is_norm))
    axes[1, 1].grid(True)
    _robust_ylim(axes[1, 1], s["ne"])

    vx = s["vx"]
    vz = s["vz"]
    _safe_line(axes[2, 0], z_axis, vx, "vx")
    _safe_line(axes[2, 0], z_axis, vz, "vz")
    axes[2, 0].set_title("Velocity")
    axes[2, 0].set_xlabel(x_label)
    axes[2, 0].set_ylabel(_var_ylabel("vx", is_norm))
    axes[2, 0].grid(True)
    handles, labels = axes[2, 0].get_legend_handles_labels()
    if handles:
        axes[2, 0].legend(loc="best")
    _robust_ylim(axes[2, 0], np.concatenate([vx, vz]))

    _safe_line(axes[2, 1], z_axis, s["Pe"])
    axes[2, 1].set_title("Pe")
    axes[2, 1].set_xlabel(x_label)
    axes[2, 1].set_ylabel(_var_ylabel("Pe", is_norm))
    axes[2, 1].grid(True)
    _robust_ylim(axes[2, 1], s["Pe"])

    fig.tight_layout()

    stem = os.path.splitext(os.path.basename(state_path))[0]
    out_path = os.path.join(results_dir, "spatial_profile_{0}.png".format(stem))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print("Using state file: {0}".format(state_path))
    return out_path


def parse_args() -> argparse.Namespace:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output, default_results = _resolve_default_paths(script_dir)

    parser = argparse.ArgumentParser(
        description="Plot full-position profile for a specified state file."
    )
    parser.add_argument("--output", default=cfg.OUTPUT_DIR or default_output, help="Path to output directory")
    parser.add_argument("--results", default=cfg.RESULTS_DIR or default_results, help="Path to results directory")
    parser.add_argument(
        "--state-file",
        default=cfg.SPACE_STATE_FILE,
        help="State file name/path, or 'latest'",
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
    normalize = cfg.NORMALIZE_OUTPUT or args.normalize
    out_path = run(args.output, args.results, args.state_file, normalize, args.n0)
    print("Saved: {0}".format(out_path))


if __name__ == "__main__":
    main()
