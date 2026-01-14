# src/dbmle/stata_bridge.py

from __future__ import annotations
from typing import Any, Dict, Tuple, List, Optional

import dbmle.core as core
from dbmle.core import DBMLEResult

try:
    # Only available inside Stata's embedded Python
    from sfi import Scalar, Macro, Matrix
except Exception:
    Scalar = Macro = Matrix = None


# Use Python NaN as "missing" fill for Stata matrices
MISSING = float("nan")


THETA_TO_NAME = {
    "theta11": "always",
    "theta10": "complier",
    "theta01": "defier",
    "theta00": "never",
}


def _require_stata() -> None:
    if Scalar is None or Macro is None or Matrix is None:
        raise RuntimeError(
            "Stata bridge requires Stata's embedded Python (module 'sfi' not found). "
            "Run this from within Stata using a python: ... end block."
        )


def _rname(prefix: str, name: str) -> str:
    return f"r({prefix}{name})"


def _safe_int(x: Any) -> int:
    return int(x)


def _set_intervals_matrix(
    *,
    prefix: str,
    base_name: str,
    intervals: Optional[List[Tuple[Any, Any]]],
    store_k: bool = True,
) -> None:
    _require_stata()

    mat = _rname(prefix, base_name)
    ksc = _rname(prefix, f"{base_name}_k")

    if not intervals:
        Matrix.create(mat, 1, 2, MISSING)
        if store_k:
            Scalar.setValue(ksc, 0.0)
        return

    k = len(intervals)
    Matrix.create(mat, k, 2, MISSING)

    for i, (lo, hi) in enumerate(intervals):
        Matrix.storeAt(mat, i, 0, float(_safe_int(lo)))
        Matrix.storeAt(mat, i, 1, float(_safe_int(hi)))

    if store_k:
        Scalar.setValue(ksc, float(k))


def _set_four_interval_mats_from_obj(prefix: str, suffix: str, obj: Dict[str, Any]) -> None:
    for theta in ("theta11", "theta10", "theta01", "theta00"):
        name = THETA_TO_NAME[theta]
        intervals = obj.get(f"{theta}_intervals")
        _set_intervals_matrix(prefix=prefix, base_name=f"{name}_{suffix}", intervals=intervals)


def _set_four_interval_mats_from_intervals_dict(prefix: str, suffix: str, intervals_dict: Dict[str, Any]) -> None:
    for theta in ("theta11", "theta10", "theta01", "theta00"):
        name = THETA_TO_NAME[theta]
        intervals = intervals_dict.get(theta)
        _set_intervals_matrix(prefix=prefix, base_name=f"{name}_{suffix}", intervals=intervals)


def _set_four_simple_interval_mats(prefix: str, suffix: str, bounds: Dict[str, Tuple[Any, Any]]) -> None:
    for theta in ("theta11", "theta10", "theta01", "theta00"):
        name = THETA_TO_NAME[theta]
        lo, hi = bounds[theta]
        _set_intervals_matrix(prefix=prefix, base_name=f"{name}_{suffix}", intervals=[(lo, hi)])


def set_r_from_result(res: Dict[str, Any], *, prefix: str = "") -> None:
    _require_stata()

    inp = res["summary"]["inputs"]  # expects n,m,c
    outmode = res["meta"]["output"]

    mle = res.get("mle", {})
    mle_list: List[Tuple[int, int, int, int]] = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("Result missing mle.mle_list; cannot populate r(mle_list).")

    Scalar.setValue(_rname(prefix, "n"), float(inp["n"]))
    Scalar.setValue(_rname(prefix, "intervention"), float(inp["m"]))
    Scalar.setValue(_rname(prefix, "control"), float(inp["c"]))

    t11, t10, t01, t00 = mle_list[0]
    Scalar.setValue(_rname(prefix, "always_mle"), float(t11))
    Scalar.setValue(_rname(prefix, "complier_mle"), float(t10))
    Scalar.setValue(_rname(prefix, "defier_mle"), float(t01))
    Scalar.setValue(_rname(prefix, "never_mle"), float(t00))

    k = len(mle_list)
    
    Matrix.create(_rname(prefix, "mle"), 4, k, MISSING)
    
    try:
        Matrix.setRowNames(_rname(prefix, "mle"), ["always", "complier", "defier", "never"])
    except Exception:
        pass
    
    try:
        Matrix.setColNames(_rname(prefix, "mle"), [f"{i+1}" for i in range(k)])
    except Exception:
        pass
    
    for j, (a, c_, d, n_) in enumerate(mle_list):  # j = 0..k-1
        Matrix.storeAt(_rname(prefix, "mle"), 0, j, float(a))   # always
        Matrix.storeAt(_rname(prefix, "mle"), 1, j, float(c_))  # complier
        Matrix.storeAt(_rname(prefix, "mle"), 2, j, float(d))   # defier
        Matrix.storeAt(_rname(prefix, "mle"), 3, j, float(n_))  # never

    num_mles = len(mle_list)
    Scalar.setValue(_rname(prefix, "num_mles"), float(num_mles))

    # Store as r(always_scs), r(complier_scs), ... matrices
    if outmode in ("basic", "auxiliary"):
        g = res.get("global_95_scs")
        if isinstance(g, dict):
            if isinstance(g.get("intervals"), dict):
                _set_four_interval_mats_from_intervals_dict(prefix, "scs", g["intervals"])
            else:
                # Create "empty" matrices with k=0 so Stata code doesn't break.
                for name in ("always", "complier", "defier", "never"):
                    _set_intervals_matrix(prefix=prefix, base_name=f"{name}_scs", intervals=None)

    if outmode == "auxiliary":
        supports = res.get("supports", {})

        # Largest Possible Support (simple contiguous bounds)
        lps = supports.get("largest_possible_support")
        if isinstance(lps, dict) and all(k in lps for k in ("theta11", "theta10", "theta01", "theta00")):
            _set_four_simple_interval_mats(prefix, "lps", lps)

        # Estimated Fréchet bounds: stored as *_intervals
        efb = supports.get("estimated_frechet_bounds")
        if isinstance(efb, dict):
            _set_four_interval_mats_from_obj(prefix, "frechet", efb)
        else:
            for name in ("always", "complier", "defier", "never"):
                _set_intervals_matrix(prefix=prefix, base_name=f"{name}_frechet", intervals=None)

        # 95% SCS within estimated Fréchet set: stored as *_intervals
        fscs = res.get("frechet_95_scs")
        if isinstance(fscs, dict):
            _set_four_interval_mats_from_obj(prefix, "frechet_scs", fscs)
        else:
            for name in ("always", "complier", "defier", "never"):
                _set_intervals_matrix(prefix=prefix, base_name=f"{name}_frechet_scs", intervals=None)

        # Optional diagnostics string: keep as local (string)
        diag = res.get("meta", {}).get("diagnostics_str")
        if isinstance(diag, str):
            Macro.setLocal(f"{prefix}diagnostics", diag)

    # ---- approx-mode exports ----
    if outmode == "approx":
        approx_meta = res.get("meta", {}).get("approx")
        if isinstance(approx_meta, str):
            Macro.setLocal(f"{prefix}approx", approx_meta)

    # Printable report: keep as local (string)
    if "report" in res and isinstance(res["report"], str):
        Macro.setLocal(f"{prefix}report", res["report"])


def dbmle_to_r(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    *,
    output: str = "basic",
    level: float = 0.95,
    show_progress: bool = False,
    prefix: str = "",
    return_result: bool = False,
) -> Optional[DBMLEResult]:
    """
    Compute dbmle() and populate Stata outputs.
    """
    res = core.dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    set_r_from_result(res, prefix=prefix)

    report = res.get("report")
    if isinstance(report, str) and report.strip():
        print(report)

    if return_result:
        return res
    return None
