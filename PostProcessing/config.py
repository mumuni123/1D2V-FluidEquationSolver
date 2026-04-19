from __future__ import annotations

import os
from typing import List

# Base paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")

# Unified entry default mode: "time" or "space"
DEFAULT_MODE = "time"

# Output normalization switch for plotting.
# False: use physical SI units from CSV.
# True: convert to normalized variables using initial density n0.
NORMALIZE_OUTPUT = True

# Initial density n0 used for normalization reference.
INITIAL_DENSITY_M3 = 3.0e26

# Time-evolution defaults
TIME_POSITIONS_UM: List[float] = [4.0]
TIME_VARIABLES: List[str] = ["Ex", "Ez", "By", "ne", "vx", "vz", "Pe"]

# Spatial-profile defaults
SPACE_STATE_FILE = "state_1798800.csv"


def positions_to_text(values: List[float]) -> str:
    return ",".join("{0:g}".format(v) for v in values)


def variables_to_text(values: List[str]) -> str:
    return ",".join(values)
