#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from typing import List, Sequence, Tuple

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
    _nearest_value_at_position,
    _nearest_indices,
    _parse_positions_um,
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


def _default_fft_xlim(freq: np.ndarray, freq_label: str, n0: float) -> Tuple[float, float]:
    if freq.size == 0:
        return 0.0, 1.0

    f_max = float(np.nanmax(freq))
    if not np.isfinite(f_max) or f_max <= 0.0:
        return 0.0, 1.0

    refs = norm.build_refs(n0)
    f_pe_hz = refs["omega_pe"] / (2.0 * np.pi)
    if freq_label == "frequency (THz)":
        limit = 10.0 * f_pe_hz / 1.0e12
    elif freq_label == "frequency (1/normalized time)":
        limit = 10.0 / (2.0 * np.pi)
    else:
        limit = f_max

    if not np.isfinite(limit) or limit <= 0.0:
        limit = f_max
    return 0.0, min(float(limit), f_max)


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
    positions_um: Sequence[float],
    normalize: bool,
    n0: float,
    fft_max: float | None = None,
) -> str:
    if len(positions_um) == 0:
        raise ValueError("No FFT positions provided")

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

    time_x, x_label = _build_time_axis_with_inference(output_dir, files, n0)
    ex = np.full((len(positions_um), len(files)), np.nan, dtype=np.float64)

    for j, f in enumerate(files):
        s = _load_csv(f)
        sv, _ = norm.convert_state_columns(s, is_norm, n0)
        s_phys, _ = norm.convert_state_columns(s, False, n0)
        for vi, pos_um in enumerate(positions_um):
            value, actual = _nearest_value_at_position(s_phys, sv, "Ex", pos_um)
            ex[vi, j] = value
            if np.isfinite(actual):
                actual_pos[vi] = actual

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.0))

    y_collect: List[np.ndarray] = []
    for vi, pos in enumerate(actual_pos):
        y = ex[vi, :]
        y_collect.append(y)
        _safe_line(axes[0], time_x, y, "z={0:.3f} um".format(pos))

    axes[0].set_title("Ex time evolution at selected positions")
    axes[0].set_xlabel(x_label)
    axes[0].set_ylabel(_var_ylabel("Ex", is_norm))
    axes[0].grid(True)
    axes[0].legend(loc="best", fontsize=8)
    _robust_ylim(axes[0], np.concatenate(y_collect))

    finite_t = time_x[np.isfinite(time_x)]
    if finite_t.size > 1:
        axes[0].set_xlim(float(np.min(finite_t)), float(np.max(finite_t)))

    freq = np.asarray([], dtype=np.float64)
    freq_label = "frequency"
    for vi, pos in enumerate(actual_pos):
        freq, amp, freq_label = _fft_amplitude(time_x, ex[vi, :], x_label)
        _safe_line(axes[1], freq, amp, "z={0:.3f} um".format(pos))

    axes[1].set_title("FFT spectrum of Ex at selected positions")
    axes[1].set_xlabel(freq_label)
    axes[1].set_ylabel("amplitude")
    axes[1].grid(True)
    axes[1].legend(loc="best", fontsize=8)
    if freq.size > 1:
        if fft_max is None:
            axes[1].set_xlim(*_default_fft_xlim(freq, freq_label, n0))
        else:
            axes[1].set_xlim(float(freq[0]), min(float(fft_max), float(freq[-1])))

    fig.tight_layout()
    out_path = os.path.join(results_dir, "ex_time_fft_at_positions.png")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    print("Selected FFT positions (requested -> actual):")
    for req, act in zip(positions_um, actual_pos):
        print("- {0:.4f} um -> {1:.4f} um".format(req, act))

    return out_path


def parse_args() -> argparse.Namespace:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output, default_results = _resolve_default_paths(script_dir)

    parser = argparse.ArgumentParser(
        description="Plot Ex time evolution and FFT spectra at selected positions."
    )
    parser.add_argument("--output", default=cfg.OUTPUT_DIR or default_output, help="Path to output directory")
    parser.add_argument("--results", default=cfg.RESULTS_DIR or default_results, help="Path to results directory")
    parser.add_argument(
        "--positions-um",
        default=cfg.positions_to_text(cfg.FFT_POSITIONS_UM),
        help="Comma-separated FFT positions in um, e.g. 0.5,2.0,4.0",
    )
    parser.add_argument(
        "--position-um",
        type=float,
        default=None,
        help="Single FFT position in um. Overrides --positions-um when provided.",
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
    parser.add_argument(
        "--fft-max",
        type=float,
        default=None,
        help="Maximum FFT x-axis value in the displayed frequency unit; default is 3*f_pe",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    normalize = cfg.NORMALIZE_OUTPUT or args.normalize
    positions_um = [args.position_um] if args.position_um is not None else _parse_positions_um(args.positions_um)
    out_path = run(args.output, args.results, positions_um, normalize, args.n0, args.fft_max)
    print("Saved: {0}".format(out_path))


if __name__ == "__main__":
    main()
