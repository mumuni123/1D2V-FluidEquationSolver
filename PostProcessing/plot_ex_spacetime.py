#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import config as cfg
import normalization as norm
from plot_ex_time_fft_at_position import _build_time_axis_with_inference
from plot_time_evolution_at_positions import (
    _list_state_files,
    _load_csv,
    _resolve_default_paths,
    _var_ylabel,
)


def _edge_axis(x: np.ndarray) -> np.ndarray:
    if x.size < 2:
        raise ValueError("Need at least 2 grid points")

    edges = np.empty(x.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def _finite_limits(arr: np.ndarray, percentile: float) -> Tuple[float, float]:
    x = np.asarray(arr, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("No finite Ex values to plot")

    if percentile <= 0.0 or percentile >= 100.0:
        vmax = float(np.max(np.abs(x)))
    else:
        vmax = float(np.percentile(np.abs(x), percentile))

    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = float(np.max(np.abs(x)))
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = 1.0
    return -vmax, vmax


def _range_mask(axis: np.ndarray, lo: float | None, hi: float | None) -> np.ndarray:
    mask = np.isfinite(axis)
    if lo is not None:
        mask &= axis >= lo
    if hi is not None:
        mask &= axis <= hi
    if np.count_nonzero(mask) == 0:
        raise ValueError("Selected axis range contains no data")
    return mask


def _range_tag(name: str, lo: float | None, hi: float | None) -> str:
    if lo is None and hi is None:
        return ""
    lo_text = "min" if lo is None else "{0:g}".format(lo)
    hi_text = "max" if hi is None else "{0:g}".format(hi)
    return "_{0}_{1}_to_{2}".format(name, lo_text, hi_text)


def run(
    output_dir: str,
    results_dir: str,
    normalize: bool,
    n0: float,
    z_min_um: float | None,
    z_max_um: float | None,
    time_min: float | None,
    time_max: float | None,
    vlim_percentile: float,
    cmap: str,
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

    time_x, time_label = _build_time_axis_with_inference(output_dir, files, n0)

    z_mask = _range_mask(z_um, z_min_um, z_max_um)
    t_mask = _range_mask(time_x, time_min, time_max)
    selected_files = [f for f, keep in zip(files, t_mask) if keep]
    selected_time = time_x[t_mask]
    selected_z = z_um[z_mask]

    ex_zt = np.full((len(selected_files), np.count_nonzero(z_mask)), np.nan, dtype=np.float64)
    for j, f in enumerate(selected_files):
        s = _load_csv(f)
        sv, _ = norm.convert_state_columns(s, is_norm, n0)
        s_phys, _ = norm.convert_state_columns(s, False, n0)
        z_frame = np.asarray(s_phys["z"], dtype=np.float64) * 1.0e6
        ex_frame = np.asarray(sv["Ex"], dtype=np.float64)
        mask = np.isfinite(z_frame) & np.isfinite(ex_frame)
        if np.count_nonzero(mask) < 2:
            continue

        z_valid = z_frame[mask]
        ex_valid = ex_frame[mask]
        order = np.argsort(z_valid)
        z_valid = z_valid[order]
        ex_valid = ex_valid[order]
        ex_zt[j, :] = np.interp(selected_z, z_valid, ex_valid, left=np.nan, right=np.nan)

    vmin, vmax = _finite_limits(ex_zt, vlim_percentile)
    z_edges = _edge_axis(selected_z)
    t_edges = _edge_axis(selected_time)

    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    mesh = ax.pcolormesh(
        z_edges,
        t_edges,
        ex_zt,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(_var_ylabel("Ex", is_norm))

    ax.set_title("Ex spacetime map")
    ax.set_xlabel("z (um)")
    ax.set_ylabel(time_label)
    ax.set_xlim(float(np.min(selected_z)), float(np.max(selected_z)))
    ax.set_ylim(float(np.min(selected_time)), float(np.max(selected_time)))

    fig.tight_layout()
    name = "ex_spacetime"
    name += _range_tag("z", z_min_um, z_max_um)
    name += _range_tag("t", time_min, time_max)
    out_path = os.path.join(results_dir, "{0}.png".format(name))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    print("z range: {0:.4g} to {1:.4g} um".format(float(np.min(selected_z)), float(np.max(selected_z))))
    print("time range: {0:.4g} to {1:.4g} ({2})".format(
        float(np.min(selected_time)),
        float(np.max(selected_time)),
        time_label,
    ))
    return out_path


def parse_args() -> argparse.Namespace:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output, default_results = _resolve_default_paths(script_dir)

    parser = argparse.ArgumentParser(description="Plot Ex(z,t) spacetime map.")
    parser.add_argument("--output", default=cfg.OUTPUT_DIR or default_output, help="Path to output directory")
    parser.add_argument("--results", default=cfg.RESULTS_DIR or default_results, help="Path to results directory")
    parser.add_argument("--z-min-um", type=float, default=None, help="Minimum z in um")
    parser.add_argument("--z-max-um", type=float, default=None, help="Maximum z in um")
    parser.add_argument("--time-min", type=float, default=None, help="Minimum time in the plotted time-axis unit")
    parser.add_argument("--time-max", type=float, default=None, help="Maximum time in the plotted time-axis unit")
    parser.add_argument(
        "--vlim-percentile",
        type=float,
        default=99.0,
        help="Symmetric color limit percentile of |Ex|; use 100 for full range",
    )
    parser.add_argument("--cmap", default="RdBu_r", help="Matplotlib colormap name")
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
    out_path = run(
        args.output,
        args.results,
        normalize,
        args.n0,
        args.z_min_um,
        args.z_max_um,
        args.time_min,
        args.time_max,
        args.vlim_percentile,
        args.cmap,
    )
    print("Saved: {0}".format(out_path))


if __name__ == "__main__":
    main()
