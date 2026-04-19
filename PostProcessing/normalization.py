from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np

QE = 1.602176634e-19
ME = 9.1093837015e-31
EPS0 = 8.8541878128e-12
C = 299792458.0


def build_refs(n0: float) -> Dict[str, float]:
    omega_pe = math.sqrt(n0 * QE * QE / (ME * EPS0))
    return {
        "n0": n0,
        "omega_pe": omega_pe,
        "v_scale": C,
        "e_scale": ME * C * omega_pe / QE,
        "b_scale": ME * omega_pe / QE,
        "p_scale": ME * C * C * n0,
    }


def _col(data: np.ndarray, name: str) -> np.ndarray:
    if name not in data.dtype.names:
        raise KeyError("Column '{0}' is missing".format(name))
    return np.asarray(data[name], dtype=np.float64)


def convert_state_columns(data: np.ndarray, normalize: bool, n0: float) -> Tuple[Dict[str, np.ndarray], bool]:
    ref = build_refs(n0)
    names = set(data.dtype.names)

    has_phys = {"z", "Ex", "Ez", "By", "ne", "vx", "vz", "Pe"}.issubset(names)
    has_tilde = {
        "z_tilde",
        "Ex_tilde",
        "Ez_tilde",
        "By_tilde",
        "ne_tilde",
        "vx_tilde",
        "vz_tilde",
        "Pe_tilde",
    }.issubset(names)

    if normalize:
        if has_tilde:
            out = {
                "z": _col(data, "z_tilde"),
                "Ex": _col(data, "Ex_tilde"),
                "Ez": _col(data, "Ez_tilde"),
                "By": _col(data, "By_tilde"),
                "ne": _col(data, "ne_tilde"),
                "vx": _col(data, "vx_tilde"),
                "vz": _col(data, "vz_tilde"),
                "Pe": _col(data, "Pe_tilde"),
            }
            return out, True

        if not has_phys:
            raise KeyError("State file lacks required columns for normalization")

        out = {
            "z": _col(data, "z") * ref["omega_pe"] / ref["v_scale"],
            "Ex": _col(data, "Ex") / ref["e_scale"],
            "Ez": _col(data, "Ez") / ref["e_scale"],
            "By": _col(data, "By") / ref["b_scale"],
            "ne": _col(data, "ne") / ref["n0"],
            "vx": _col(data, "vx") / ref["v_scale"],
            "vz": _col(data, "vz") / ref["v_scale"],
            "Pe": _col(data, "Pe") / ref["p_scale"],
        }
        return out, True

    if has_phys:
        out = {
            "z": _col(data, "z"),
            "Ex": _col(data, "Ex"),
            "Ez": _col(data, "Ez"),
            "By": _col(data, "By"),
            "ne": _col(data, "ne"),
            "vx": _col(data, "vx"),
            "vz": _col(data, "vz"),
            "Pe": _col(data, "Pe"),
        }
        return out, False

    if has_tilde:
        out = {
            "z": _col(data, "z_tilde") * ref["v_scale"] / ref["omega_pe"],
            "Ex": _col(data, "Ex_tilde") * ref["e_scale"],
            "Ez": _col(data, "Ez_tilde") * ref["e_scale"],
            "By": _col(data, "By_tilde") * ref["b_scale"],
            "ne": _col(data, "ne_tilde") * ref["n0"],
            "vx": _col(data, "vx_tilde") * ref["v_scale"],
            "vz": _col(data, "vz_tilde") * ref["v_scale"],
            "Pe": _col(data, "Pe_tilde") * ref["p_scale"],
        }
        return out, False

    raise KeyError("State file lacks required columns")


def diagnostics_time_axis(diag: np.ndarray, normalize: bool, n0: float) -> Tuple[np.ndarray, str]:
    names = set(diag.dtype.names)
    ref = build_refs(n0)

    if normalize:
        if "time_tilde" in names:
            return np.asarray(diag["time_tilde"], dtype=np.float64), "time"
        if "time_fs" in names:
            t_tilde = np.asarray(diag["time_fs"], dtype=np.float64) * 1.0e-15 * ref["omega_pe"]
            return t_tilde, "time"
    else:
        if "time_fs" in names:
            return np.asarray(diag["time_fs"], dtype=np.float64), "time (fs)"
        if "time_tilde" in names:
            t_fs = np.asarray(diag["time_tilde"], dtype=np.float64) / ref["omega_pe"] * 1.0e15
            return t_fs, "time (fs)"

    raise KeyError("No recognized time column in diagnostics.csv")