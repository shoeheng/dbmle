# src/dbmle/stata_eclass.py
"""
Stata e()-class bridge for dbmle.

Call from Stata:

    python:
        from dbmle import dbmle_to_eclass
        dbmle_to_eclass(50, 11, 23, 31, output="basic", level=0.95, show_progress=False)
    end
    ereturn list
    matrix list e(b)
    matrix list e(V)

Design goals:
- Post ONLY:
    * MLE joint counts (and all tied MLEs, if any)
    * Global (design-based) smallest credible sets at the requested level
      (NOT the Frechet-conditional SCS)
- If output="approx": post only MLEs (no credible sets)
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple, Optional

from dbmle.core import dbmle

Theta = Tuple[int, int, int, int]


def _require_stata_modules():
    """
    Import Stata Python integration modules.

    Returns (stata, Matrix) where:
      - stata.run(cmd) executes a Stata command
      - Matrix.store(name, list_of_lists) stores a Stata matrix
    """
    try:
        import stata  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Could not import the 'stata' Python module. "
            "This function must be called from within Stata's Python integration."
        ) from e

    try:
        from sfi import Matrix  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Could not import 'sfi.Matrix'. "
            "This function must be called from within Stata (version with Python integration)."
        ) from e

    return stata, Matrix


def _matrix_store(
    Matrix,
    name: str,
    data: Sequence[Sequence[float]],
    rownames: Optional[List[str]] = None,
    colnames: Optional[List[str]] = None,
) -> None:
    """Store a matrix in Stata via sfi.Matrix, and optionally attach names."""
    Matrix.store(name, data)

    # Some Stata versions expose setRowNames/setColNames; guard for safety.
    if rownames is not None and hasattr(Matrix, "setRowNames"):
        Matrix.setRowNames(name, rownames)
    if colnames is not None and hasattr(Matrix, "setColNames"):
        Matrix.setColNames(name, colnames)


def _intervals_to_matrix(intervals: List[Tuple[int, int]]) -> List[List[float]]:
    """Convert [(lo,hi), ...] to [[lo,hi], ...] for Stata matrix storage."""
    return [[float(a), float(b)] for a, b in intervals]


def dbmle_to_eclass(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    *,
    output: str = "basic",
    level: float = 0.95,
    show_progress: bool = True,
    prefix: str = "dbmle",
) -> Dict[str, Any]:
    """
    Run dbmle(...) and post selected outputs into Stata's e() returns.

    Posted items:
      e(b)     : 1x4 row vector of the FIRST MLE joint count (A,C,D,N)
      e(V)     : 4x4 matrix of missing values (no variance defined here)
      e(b_all) : (#ties)x4 matrix of ALL tied MLEs (1 row if unique)

    If output != "approx", also:
      e(scs_A), e(scs_C), e(scs_D), e(scs_N)
          : matrices with 2 columns [lo hi] (one row per disjoint interval)
      e(scs_A_str), ... : locals with printable union strings

    Additionally:
      e(N) scalar, e(level) scalar, e(output) local, e(cmd) local

    Returns the underlying DBMLEResult (dict-like) for convenience.
    """
    stata, Matrix = _require_stata_modules()

    out_mode = (output or "basic").strip().lower()
    if out_mode not in {"basic", "auxiliary", "approx"}:
        raise ValueError("output must be one of {'basic','auxiliary','approx'}.")

    res = dbmle(
        xI1, xI0, xC1, xC0,
        output=out_mode,
        level=level,
        show_progress=show_progress,
    )

    # ---- Extract MLEs ----
    mle_list: List[Theta] = list(res["mle"]["mle_list"])
    if len(mle_list) == 0:
        raise RuntimeError("dbmle returned no MLEs; this should not happen.")

    # Always-takers, Compliers, Defiers, Never-takers
    colnames = ["A", "C", "D", "N"]

    # e(b): first MLE only (1x4)
    b_first = [list(map(float, mle_list[0]))]
    b_all = [list(map(float, t)) for t in mle_list]

    b_name = f"{prefix}_b"
    ball_name = f"{prefix}_b_all"
    V_name = f"{prefix}_V"

    _matrix_store(Matrix, b_name, b_first, rownames=["b"], colnames=colnames)
    _matrix_store(
        Matrix,
        ball_name,
        b_all,
        rownames=[f"mle{j+1}" for j in range(len(b_all))],
        colnames=colnames,
    )

    # No variance defined -> missing 4x4
    V = [[float("nan")] * 4 for _ in range(4)]
    _matrix_store(Matrix, V_name, V, rownames=colnames, colnames=colnames)

    # ---- Post into e() ----
    stata.run("ereturn clear")
    stata.run(f"ereturn matrix b = {b_name}")
    stata.run(f"ereturn matrix V = {V_name}")
    stata.run(f"ereturn matrix b_all = {ball_name}")

    # Basic scalars/locals
    n = int(res["summary"]["inputs"]["n"])
    stata.run(f"ereturn scalar N = {n}")
    stata.run(f"ereturn scalar level = {float(level)}")
    stata.run(f'ereturn local output "{out_mode}"')
    stata.run('ereturn local cmd "dbmle"')
    stata.run('ereturn local cmdline "dbmle_to_eclass (via Python)"')

    # ---- Credible sets (global SCS only) ----
    if out_mode != "approx":
        scs = res.get("global_95_scs", {})
        intervals = scs.get("intervals", {})
        union_str = scs.get("union_str", {})

        key_map = {"theta11": "A", "theta10": "C", "theta01": "D", "theta00": "N"}

        for key, lbl in key_map.items():
            ivs: List[Tuple[int, int]] = list(intervals.get(key, []))
            mat = _intervals_to_matrix(ivs)

            scs_name = f"{prefix}_scs_{lbl}"
            if len(mat) > 0:
                _matrix_store(Matrix, scs_name, mat, colnames=["lo", "hi"])
                stata.run(f"ereturn matrix scs_{lbl} = {scs_name}")

            s = union_str.get(key, "EMPTY")
            stata.run(f'ereturn local scs_{lbl}_str "{s}"')

    return res
