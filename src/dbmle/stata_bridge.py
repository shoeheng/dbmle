# src/dbmle/stata_bridge.py

from __future__ import annotations
from typing import Any, Dict, Tuple, List

import dbmle.core as core
from dbmle.core import DBMLEResult

try:
    # Only available inside Stata's embedded Python
    from sfi import Scalar, Macro, Matrix
except Exception:
    Scalar = Macro = Matrix = None


def _require_stata() -> None:
    if Scalar is None or Macro is None or Matrix is None:
        raise RuntimeError(
            "Stata bridge requires Stata's embedded Python (module 'sfi' not found). "
            "Run this from within Stata using a python: ... end block."
        )


def _set_union_locals(prefix: str, base_name: str, union_str: Dict[str, str]) -> None:
    """
    Store union-of-interval strings as Stata locals like:
      `<prefix>theta11_<base_name>' etc.
    """
    Macro.setLocal(f"{prefix}theta11_{base_name}", union_str["theta11"])
    Macro.setLocal(f"{prefix}theta10_{base_name}", union_str["theta10"])
    Macro.setLocal(f"{prefix}theta01_{base_name}", union_str["theta01"])
    Macro.setLocal(f"{prefix}theta00_{base_name}", union_str["theta00"])

def _intervals_to_union_str(intervals) -> str:
    """
    Convert a list of [lo, hi] (or (lo, hi)) into a union string:
      [[1,2],[5,7]] -> "[1,2] U [5,7]"
    Accepts intervals=None, empty list, numpy ints, etc.
    """
    if not intervals:
        return ""

    parts = []
    for seg in intervals:
        lo, hi = seg
        lo_i, hi_i = int(lo), int(hi)
        parts.append(f"[{lo_i},{hi_i}]")
    return " U ".join(parts)

def _set_union_locals_from_intervals(prefix: str, base_name: str, obj: Dict[str, Any]) -> None:
    """
    obj must have keys like theta11_intervals, theta10_intervals, ...
    """
    Macro.setLocal(f"{prefix}theta11_{base_name}", _intervals_to_union_str(obj.get("theta11_intervals")))
    Macro.setLocal(f"{prefix}theta10_{base_name}", _intervals_to_union_str(obj.get("theta10_intervals")))
    Macro.setLocal(f"{prefix}theta01_{base_name}", _intervals_to_union_str(obj.get("theta01_intervals")))
    Macro.setLocal(f"{prefix}theta00_{base_name}", _intervals_to_union_str(obj.get("theta00_intervals")))

def _set_interval_locals(prefix: str, base_name: str, intervals) -> None:
    """
    Store simple [lo,hi] integer intervals as Stata locals like:
      `<prefix>theta11_<base_name>' = "[lo,hi]"
    """
    for key in ("theta11", "theta10", "theta01", "theta00"):
        lo, hi = intervals[key]
        lo_i = int(lo)
        hi_i = int(hi)

        name = f"{prefix}{key}_{base_name}"
        value = f"[{lo_i},{hi_i}]"
        value = value.replace("\n", " ").replace("\r", " ")
        Macro.setLocal(name, value)


def set_r_from_result(res: Dict[str, Any], *, prefix: str = "") -> None:
    """
    Write results into Stata.

    - Numeric outputs -> r() scalars / matrices
    - Set-valued outputs (interval unions) -> locals (strings)

    prefix: optional prefix for names, e.g. prefix="dbmle_" -> r(dbmle_n), etc.
    """
    _require_stata()

    # ---- required structure ----
    inp = res["summary"]["inputs"]  # expects xI1,xI0,xC1,xC0,n,m,c
    outmode = res["meta"]["output"]

    mle = res.get("mle", {})
    mle_list: List[Tuple[int, int, int, int]] = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("Result missing mle.mle_list; cannot populate r(mle_list).")

    def rname(name: str) -> str:
        return f"r({prefix}{name})"

    # ---- scalars always available ----
    Scalar.setValue(rname("n"), float(inp["n"]))
    Scalar.setValue(rname("m"), float(inp["m"]))
    Scalar.setValue(rname("c"), float(inp["c"]))

    # First MLE as scalars
    t11, t10, t01, t00 = mle_list[0]
    Scalar.setValue(rname("theta11_mle"), float(t11))
    Scalar.setValue(rname("theta10_mle"), float(t10))
    Scalar.setValue(rname("theta01_mle"), float(t01))
    Scalar.setValue(rname("theta00_mle"), float(t00))

    # All MLEs (ties) as matrix kx4 (0-based indexing for sfi.Matrix.storeAt)
    k = len(mle_list)
    Matrix.create(rname("mle_list"), k, 4, 0)

    # Try to set column names inside the matrix itself (avoids Stata parser quirks later).
    # If the installed sfi.Matrix doesn't support setColNames, this silently does nothing.
    try:
        Matrix.setColNames(rname("mle_list"), ["always", "complier", "defier", "never"])
    except Exception:
        pass

    for i, (a, c_, d, n_) in enumerate(mle_list):  # i = 0..k-1
        Matrix.storeAt(rname("mle_list"), i, 0, float(a))
        Matrix.storeAt(rname("mle_list"), i, 1, float(c_))
        Matrix.storeAt(rname("mle_list"), i, 2, float(d))
        Matrix.storeAt(rname("mle_list"), i, 3, float(n_))

    # ---- global SCS unions (only in exhaustive modes) ----
    if outmode in ("basic", "auxiliary"):
        g = res.get("global_95_scs")
        if g and "union_str" in g:
            _set_union_locals(prefix, "scs", g["union_str"])

    # ---- auxiliary-only exports (correct key paths) ----
    if outmode == "auxiliary":
        supports = res.get("supports", {})

        # Largest Possible Support
        lps = supports.get("largest_possible_support")
        if lps and all(k in lps for k in ("theta11", "theta10", "theta01", "theta00")):
            _set_interval_locals(prefix, "lps", lps)

        # Estimated Fréchet bounds: stored as *_intervals in supports["estimated_frechet_bounds"]
        efb = supports.get("estimated_frechet_bounds")
        if isinstance(efb, dict):
            _set_union_locals_from_intervals(prefix, "frechet", efb)
        
        # 95% SCS within estimated Fréchet set: stored as *_intervals in res["frechet_95_scs"]
        fscs = res.get("frechet_95_scs")
        if isinstance(fscs, dict):
            _set_union_locals_from_intervals(prefix, "frechet_scs", fscs)

        # Optional: store additional diagnostics if present
        diag = res.get("meta", {}).get("diagnostics_str")
        if isinstance(diag, str):
            Macro.setLocal(f"{prefix}diagnostics", diag)

    # ---- approx-mode exports (best-effort) ----
    if outmode == "approx":
        approx_meta = res.get("meta", {}).get("approx")
        if isinstance(approx_meta, str):
            Macro.setLocal(f"{prefix}approx", approx_meta)

    # Store full printable report as a local if present
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
    show_progress: bool = True,
    prefix: str = "",
) -> DBMLEResult:
    """
    Compute dbmle() and populate Stata outputs.

    - r() scalars: n, m, c, theta??_mle
    - r() matrix:  mle_list  (k x 4)
    - locals:      <prefix>theta??_scs
                   plus auxiliary locals when output="auxiliary":
                     <prefix>theta??_lps
                     <prefix>theta??_frechet
                     <prefix>theta??_frechet_scs

    Returns the DBMLEResult as well.
    """
    res = core.dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )
    set_r_from_result(res, prefix=prefix)
    return res
