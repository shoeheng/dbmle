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


def _require_stata() -> None:
    if Scalar is None or Macro is None or Matrix is None:
        raise RuntimeError(
            "Stata bridge requires Stata's embedded Python (module 'sfi' not found). "
            "Run this from within Stata using a python: ... end block."
        )


def _rname(prefix: str, name: str) -> str:
    """Return-name helper for Stata r() objects."""
    return f"r({prefix}{name})"


def _safe_int(x: Any) -> int:
    """Best-effort conversion to Python int (handles numpy ints)."""
    return int(x)


def _set_intervals_matrix(
    *,
    prefix: str,
    base_name: str,
    intervals: Optional[List[Tuple[Any, Any]]],
    store_k: bool = True,
) -> None:
    """
    Store a union of integer intervals as an r() matrix with 2 columns [lo hi].

      intervals = [(lo,hi), (lo,hi), ...] -> r(<prefix><base_name>) is k x 2
      also store r(<prefix><base_name>_k) = k (number of pieces)

    If intervals is empty/None:
      create a 1x2 matrix filled with missing and set _k = 0
      (avoids 0-row matrix quirks in some Stata contexts)
    """
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
    """
    obj must have keys like:
      theta11_intervals, theta10_intervals, theta01_intervals, theta00_intervals
    Each value is a list of (lo,hi) segments (possibly empty).
    """
    _set_intervals_matrix(prefix=prefix, base_name=f"theta11_{suffix}", intervals=obj.get("theta11_intervals"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta10_{suffix}", intervals=obj.get("theta10_intervals"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta01_{suffix}", intervals=obj.get("theta01_intervals"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta00_{suffix}", intervals=obj.get("theta00_intervals"))


def _set_four_interval_mats_from_intervals_dict(prefix: str, suffix: str, intervals_dict: Dict[str, Any]) -> None:
    """
    intervals_dict must have keys:
      theta11, theta10, theta01, theta00
    Each value is a list of (lo,hi) segments (possibly empty).
    """
    _set_intervals_matrix(prefix=prefix, base_name=f"theta11_{suffix}", intervals=intervals_dict.get("theta11"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta10_{suffix}", intervals=intervals_dict.get("theta10"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta01_{suffix}", intervals=intervals_dict.get("theta01"))
    _set_intervals_matrix(prefix=prefix, base_name=f"theta00_{suffix}", intervals=intervals_dict.get("theta00"))


def _set_four_simple_interval_mats(prefix: str, suffix: str, bounds: Dict[str, Tuple[Any, Any]]) -> None:
    """
    bounds must have keys:
      theta11: (lo,hi), theta10: (lo,hi), theta01: (lo,hi), theta00: (lo,hi)

    Stored as 1x2 matrices in r().
    """
    for key in ("theta11", "theta10", "theta01", "theta00"):
        lo, hi = bounds[key]
        _set_intervals_matrix(prefix=prefix, base_name=f"{key}_{suffix}", intervals=[(lo, hi)])


def set_r_from_result(res: Dict[str, Any], *, prefix: str = "") -> None:
    """
    Write results into Stata r().

    - Numeric outputs -> r() scalars / matrices
    - Interval unions -> r() matrices (k x 2), with companion scalar r(<name>_k)

    prefix: optional prefix for names, e.g. prefix="dbmle_" -> r(dbmle_n), etc.
    """
    _require_stata()

    # ---- required structure ----
    inp = res["summary"]["inputs"]  # expects n,m,c
    outmode = res["meta"]["output"]

    mle = res.get("mle", {})
    mle_list: List[Tuple[int, int, int, int]] = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("Result missing mle.mle_list; cannot populate r(mle_list).")

    # ---- scalars always available ----
    Scalar.setValue(_rname(prefix, "n"), float(inp["n"]))
    Scalar.setValue(_rname(prefix, "m"), float(inp["m"]))
    Scalar.setValue(_rname(prefix, "c"), float(inp["c"]))

    # First MLE as scalars
    t11, t10, t01, t00 = mle_list[0]
    Scalar.setValue(_rname(prefix, "theta11_mle"), float(t11))
    Scalar.setValue(_rname(prefix, "theta10_mle"), float(t10))
    Scalar.setValue(_rname(prefix, "theta01_mle"), float(t01))
    Scalar.setValue(_rname(prefix, "theta00_mle"), float(t00))

    # All MLEs (ties) as matrix kx4
    k = len(mle_list)
    Matrix.create(_rname(prefix, "mle_list"), k, 4, MISSING)
    try:
        Matrix.setColNames(_rname(prefix, "mle_list"), ["always", "complier", "defier", "never"])
    except Exception:
        pass

    for i, (a, c_, d, n_) in enumerate(mle_list):
        Matrix.storeAt(_rname(prefix, "mle_list"), i, 0, float(a))
        Matrix.storeAt(_rname(prefix, "mle_list"), i, 1, float(c_))
        Matrix.storeAt(_rname(prefix, "mle_list"), i, 2, float(d))
        Matrix.storeAt(_rname(prefix, "mle_list"), i, 3, float(n_))

    num_mles = len(mle_list)
    Scalar.setValue(_rname(prefix, "num_mles"), float(num_mles))

    # ---- global SCS unions (only in exhaustive modes) ----
    # Store as r(theta**_scs) matrices (k x 2), not locals.
    if outmode in ("basic", "auxiliary"):
        g = res.get("global_95_scs")
        if isinstance(g, dict):
            # Preferred: structured intervals dict from core
            if isinstance(g.get("intervals"), dict):
                _set_four_interval_mats_from_intervals_dict(prefix, "scs", g["intervals"])
            else:
                # If intervals aren't provided, we can't reliably reconstruct pieces from union strings.
                # Create "empty" matrices with k=0 so Stata code doesn't break.
                _set_intervals_matrix(prefix=prefix, base_name="theta11_scs", intervals=None)
                _set_intervals_matrix(prefix=prefix, base_name="theta10_scs", intervals=None)
                _set_intervals_matrix(prefix=prefix, base_name="theta01_scs", intervals=None)
                _set_intervals_matrix(prefix=prefix, base_name="theta00_scs", intervals=None)

    # ---- auxiliary-only exports ----
    if outmode == "auxiliary":
        supports = res.get("supports", {})

        # Largest Possible Support (simple contiguous bounds)
        lps = supports.get("largest_possible_support")
        if isinstance(lps, dict) and all(k in lps for k in ("theta11", "theta10", "theta01", "theta00")):
            _set_four_simple_interval_mats(prefix, "lps", lps)
        else:
            # ensure presence (optional)
            pass

        # Estimated Fréchet bounds: stored as *_intervals
        efb = supports.get("estimated_frechet_bounds")
        if isinstance(efb, dict):
            _set_four_interval_mats_from_obj(prefix, "frechet", efb)
        else:
            _set_intervals_matrix(prefix=prefix, base_name="theta11_frechet", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta10_frechet", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta01_frechet", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta00_frechet", intervals=None)

        # 95% SCS within estimated Fréchet set: stored as *_intervals
        fscs = res.get("frechet_95_scs")
        if isinstance(fscs, dict):
            _set_four_interval_mats_from_obj(prefix, "frechet_scs", fscs)
        else:
            _set_intervals_matrix(prefix=prefix, base_name="theta11_frechet_scs", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta10_frechet_scs", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta01_frechet_scs", intervals=None)
            _set_intervals_matrix(prefix=prefix, base_name="theta00_frechet_scs", intervals=None)

        # Optional diagnostics string: keep as local (string)
        diag = res.get("meta", {}).get("diagnostics_str")
        if isinstance(diag, str):
            Macro.setLocal(f"{prefix}diagnostics", diag)

    # ---- approx-mode exports (best-effort) ----
    if outmode == "approx":
        approx_meta = res.get("meta", {}).get("approx")
        if isinstance(approx_meta, str):
            Macro.setLocal(f"{prefix}approx", approx_meta)

    # Printable report: keep as local (string)
    if "report" in res and isinstance(res["report"], str):
        Macro.setLocal(f"{prefix}report", res["report"])

    if num_mles > 1:
        print(
            f"[dbmle] NOTE: {num_mles} tied MLEs detected. "
            f"See matrix {_rname(prefix, 'mle_list')} for all solutions."
        )


def dbmle_to_r(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    *,
    output: str = "basic",
    level: float = 0.95,
    show_progress: bool = True,
    prefix: str = "",
    return_result: bool = False,
) -> Optional[DBMLEResult]:
    """
    Compute dbmle() and populate Stata outputs.

    - r() scalars: n, m, c, theta??_mle, num_mles
    - r() matrix:  mle_list  (k x 4)

    Interval unions stored as r() matrices (k x 2) with companion scalar *_k:
      - r(theta11_scs), r(theta11_scs_k)  (and theta10/theta01/theta00)
      - auxiliary mode also:
          r(theta??_lps), r(theta??_lps_k)                 (1x2, k=1)
          r(theta??_frechet), r(theta??_frechet_k)
          r(theta??_frechet_scs), r(theta??_frechet_scs_k)

    Returns the DBMLEResult iff return_result=True, else None.
    """
    res = core.dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    set_r_from_result(res, prefix=prefix)

    if return_result:
        return res
    return None
