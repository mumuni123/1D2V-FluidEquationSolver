#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import config as cfg
import normalization as norm
from plot_time_evolution_at_positions import (
    _build_time_axis,
    _list_state_files,
    _load_csv,
    _nearest_indices,
    _resolve_default_paths,
    _robust_ylim,
    _safe_line,
    _var_ylabel,
)


def _frequency_axis(time_x: np.ndarray, x_label: str) -> Tuple[np.ndarray, str]:
    dt = float(np.median(np.diff(time_x)))
    if dt <= 0.0 or not np.isfinite(dt):
        raise ValueError("Time axis must be strictly increasing for FFT")

    if x_label == "time (fs)":
        freq = np.fft.rfftfreq(time_x.size, d=dt * 1.0e-15) / 1.0e12
        return freq, "frequency (THz)"

    if x_label == "time":
        freq = np.fft.rfftfreq(time_x.size, d=dt)
        return freq, "frequency (1/normalized time)"

    freq = np.fft.rfftfreq(time_x.size, d=dt)
    return freq, "frequency (1/{0})".format(x_label)


def _extract_step(path: str) -> int:
    m = re.search(r"state_(\d+)\.csv$", os.path.basename(path))
    if m is None:
        return -1
    return int(m.group(1))


def _build_time_axis_with_inference(output_dir: str, files: list[str], n0: float) -> Tuple[np.ndarray, str]:
    time_x, x_label = _build_time_axis(output_dir, files, n0)
    if x_label != "step":
        return time_x, x_label

    diag_path = os.path.join(output_dir, "diagnostics.csv")
    if not os.path.isfile(diag_path):
        return time_x, x_label

    diag = _load_csv(diag_path)
    try:
        diag_t, diag_label = norm.diagnostics_time_axis(diag, False, n0)
    except KeyError:
        return time_x, x_label

    diag_t = np.asarray(diag_t, dtype=np.float64)
    diag_t = diag_t[np.isfinite(diag_t)]
    if diag_t.size < 2:
        return time_x, x_label

    diag_dt = np.diff(diag_t)
    median_diag_dt = float(np.median(diag_dt))
    if median_diag_dt <= 0.0 or not np.isfinite(median_diag_dt):
        return time_x, x_label

    steps = np.asarray([_extract_step(f) for f in files], dtype=np.float64)
    if np.all(steps >= 0):
        unique_steps = np.unique(steps)
        if unique_steps.size >= 2:
            step_stride = float(np.median(np.diff(unique_steps)))
            if step_stride > 0.0 and np.isfinite(step_stride):
                return (steps - steps[0]) / step_stride * median_diag_dt + diag_t[0], diag_label

    return np.arange(len(files), dtype=np.float64) * median_diag_dt + diag_t[0], diag_label


def _uniform_series(time_x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(time_x) & np.isfinite(y)
    if np.count_nonzero(mask) < 4:
        raise ValueError("Need at least 4 finite samples for FFT")

    t = np.asarray(time_x[mask], dtype=np.float64)
    v = np.asarray(y[mask], dtype=np.float64)
    order = np.argsort(t)
    t = t[order]
    v = v[order]

    unique_t, unique_idx = np.unique(t, return_index=True)
    t = unique_t
    v = v[unique_idx]
    if t.size < 4:
        raise ValueError("Need at least 4 unique time samples for FFT")

    dt = np.diff(t)
    median_dt = float(np.median(dt))
    if median_dt <= 0.0 or not np.isfinite(median_dt):
        raise ValueError("Invalid time spacing for FFT")

    rel_jitter = np.max(np.abs(dt - median_dt)) / median_dt
    if rel_jitter <= 1.0e-6:
        return t, v

    t_uniform = t[0] + np.arange(t.size, dtype=np.float64) * median_dt
    return t_uniform, np.interp(t_uniform, t, v)


def _fft_amplitude(time_x: np.ndarray, ex: np.ndarray, x_label: str) -> Tuple[np.ndarray, np.ndarray, str]:
    t, y = _uniform_series(time_x, ex)
    y = y - np.mean(y)
    if y.size > 1:
        y = y * np.hanning(y.size)

    spectrum = np.fft.rfft(y)
    amp = 2.0 * np.abs(spectrum) / float(y.size)
    if amp.size:
        amp[0] *= 0.5

    freq, freq_label = _frequency_axis(t, x_label)
    return freq, amp, freq_label


def run(
    output_dir: str,
    results_dir: str,
    position_um: float,
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
    pos_idx, actual_pos = _nearest_indices(z_um, [position_um])
    i_cell = pos_idx[0]
    actual_um = actual_pos[0]

    time_x, x_label = _build_time_axis_with_inference(output_dir, files, n0)
    ex = np.full(len(files), np.nan, dtype=np.float64)

    for j, f in enumerate(files):
        s = _load_csv(f)
        sv, _ = norm.convert_state_columns(s, is_norm, n0)
        arr = sv["Ex"]
        if 0 <= i_cell < arr.size:
            ex[j] = arr[i_cell]

    freq, amp, freq_label = _fft_amplitude(time_x, ex, x_label)

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.0))

    _safe_line(axes[0], time_x, ex, "z={0:.3f} um".format(actual_um))
    axes[0].set_title("Ex time evolution at z={0:.3f} um".format(actual_um))
    axes[0].set_xlabel(x_label)
    axes[0].set_ylabel(_var_ylabel("Ex", is_norm))
    axes[0].grid(True)
    axes[0].legend(loc="best", fontsize=8)
    _robust_ylim(axes[0], ex)

    finite_t = time_x[np.isfinite(time_x)]
    if finite_t.size > 1:
        axes[0].set_xlim(float(np.min(finite_t)), float(np.max(finite_t)))

    _safe_line(axes[1], freq, amp, "FFT amplitude")
    axes[1].set_title("FFT spectrum of Ex")
    axes[1].set_xlabel(freq_label)
    axes[1].set_ylabel("amplitude")
    axes[1].grid(True)
    if freq.size > 1:
        axes[1].set_xlim(float(freq[0]), float(freq[-1]))

    fig.tight_layout()
    out_path = os.path.join(results_dir, "ex_time_fft_at_{0:g}um.png".format(position_um))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    print("Selected position (requested -> actual):")
    print("- {0:.4f} um -> {1:.4f} um".format(position_um, actual_um))

    return out_path


def parse_args() -> argparse.Namespace:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output, default_results = _resolve_default_paths(script_dir)

    parser = argparse.ArgumentParser(
        description="Plot Ex time evolution at one position and its FFT spectrum."
    )
    parser.add_argument("--output", default=cfg.OUTPUT_DIR or default_output, help="Path to output directory")
    parser.add_argument("--results", default=cfg.RESULTS_DIR or default_results, help="Path to results directory")
    parser.add_argument(
        "--position-um",
        type=float,
        default=cfg.TIME_POSITIONS_UM[0],
        help="Position in um, e.g. 23.0",
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
    out_path = run(args.output, args.results, args.position_um, normalize, args.n0)
    print("Saved: {0}".format(out_path))


if __name__ == "__main__":
    main()
