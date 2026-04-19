#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os

import config as cfg
from plot_spatial_profile_for_state import run as run_space
from plot_time_evolution_at_positions import run as run_time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified entry for post-processing plots."
    )
    parser.add_argument(
        "mode",
        choices=["time", "space"],
        nargs="?",
        default=cfg.DEFAULT_MODE,
        help="time: full-time evolution at specified positions; space: full-position profile for one state file",
    )
    parser.add_argument("--output", default=None, help="Path to output directory")
    parser.add_argument("--results", default=None, help="Path to results directory")
    parser.add_argument(
        "--positions-um",
        default=cfg.positions_to_text(cfg.TIME_POSITIONS_UM),
        help="Used in mode=time, comma-separated positions in um",
    )
    parser.add_argument(
        "--variables",
        default=cfg.variables_to_text(cfg.TIME_VARIABLES),
        help="Used in mode=time, comma-separated variable names",
    )
    parser.add_argument(
        "--state-file",
        default=cfg.SPACE_STATE_FILE,
        help="Used in mode=space, state file name/path or 'latest'",
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
    script_dir = os.path.dirname(os.path.abspath(__file__))

    if args.mode == "time":
        from plot_time_evolution_at_positions import _resolve_default_paths as p1

        default_output, default_results = p1(script_dir)
        output_dir = args.output or cfg.OUTPUT_DIR or default_output
        results_dir = args.results or cfg.RESULTS_DIR or default_results

        positions_um = [float(x.strip()) for x in args.positions_um.split(",") if x.strip()]
        variables = [x.strip() for x in args.variables.split(",") if x.strip()]

        normalize = cfg.NORMALIZE_OUTPUT or args.normalize
        out_path = run_time(output_dir, results_dir, positions_um, variables, normalize, args.n0)
        print("Saved: {0}".format(out_path))
        return

    if args.mode == "space":
        from plot_spatial_profile_for_state import _resolve_default_paths as p2

        default_output, default_results = p2(script_dir)
        output_dir = args.output or cfg.OUTPUT_DIR or default_output
        results_dir = args.results or cfg.RESULTS_DIR or default_results

        normalize = cfg.NORMALIZE_OUTPUT or args.normalize
        out_path = run_space(output_dir, results_dir, args.state_file, normalize, args.n0)
        print("Saved: {0}".format(out_path))
        return


if __name__ == "__main__":
    main()
