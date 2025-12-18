# src/dbmle/stata_eclass.py
"""
Bridge dbmle() -> Stata e() returns.

This is meant to be called ONLY from inside Stata's:

    python:
        ...
    end

It posts:
  - e(b):     first MLE (A,C,D,N) as 1x4
  - e(b_all): all tied MLEs as (#ties)x4
  - e(V):     missing 4x4 (no variance defined)

If output != "approx", also posts global smallest credible sets (NOT Frechet-conditional):
  - e(scs_A), e(scs_C), e(scs_D), e(scs_N): kx2 matrices [lo hi]
  - e(scs_A_str) etc: locals with printable union strings
"""

from __future__ import annotations
from typing import Any, Dict, List, Tuple, Optional, Sequence

from dbmle.core import dbmle

Theta = Tuple[int, int, int, int]


def _require_stata():
    """
    Require Stata embedded Python integration (python: ... end).
    """
    try:
        from sfi import Matrix, SFIToolkit  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Stata integration modules not found. "
            "This function must be run inside Stata's `python:` block."
        ) from e
    return Matrix, SFIToolkit


def _store_matrix(
    Matrix,
    name: str,
    data: Sequence[Sequence[float]],
    rownames: Optional[List[str]] = None,
    colnames: Optional[List[str]] = None,
) -> None:
    """
    Store a matrix into Stata, optionally naming rows/cols when supported.
    """
    Matrix.store(name, data)
    if rownames is not None and hasattr(Matrix, "setRowNames"):
        Matrix.setRowNames(name, rownames)
    if colnames is not None and hasattr(Matrix, "setColNames"):
        Matrix.setColNames(name, colnames)


def _intervals_to_mat(intervals: List[Tuple[int, int]]) -> List[List[float]]:
    return [[float(lo), float(hi)] for (lo, hi) in intervals]


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
    Run dbmle(...) and post selected outputs into Stata e().

    Only posts:
      - MLE counts (with ties)
      - global smallest credible sets (NOT conditional on Frechet)
    If output='approx': only MLEs.
    """
    Matrix, SFIToolkit = _require_stata()

    out_mode = (output or "basic").strip().lower()
    if out_mode not in {"basic", "auxiliary", "approx"}:
        raise ValueError("output must be one of {'basic','auxiliary','approx'}.")

    res = dbmle(
        xI1, xI0, xC1, xC0,
        output=out_mode,
        level=level,
        show_progress=show_progress,
    )

    # ---- Extract tied MLEs (A,C,D,N = theta11,theta10,theta01,theta00) ----
    mle_list: List[Theta] = list(res["mle"]["mle_list"])
    if not mle_list:
        raise RuntimeError("No MLEs returned by dbmle().")

    colnames = ["A", "C", "D", "N"]

    b_first = [list(map(float, mle_list[0]))]                 # 1x4
    b_all   = [list(map(float, t)) for t in mle_list]         # (#ties)x4
    V       = [[float("nan")] * 4 for _ in range(4)]          # 4x4 missing

    # Store matrices in Stata (temporary names)
    b_name    = f"{prefix}__b"
    ball_name = f"{prefix}__b_all"
    V_name    = f"{prefix}__V"

    _store_matrix(Matrix, b_name, b_first, rownames=["b"], colnames=colnames)
    _store_matrix(
        Matrix,
        ball_name,
        b_all,
        rownames=[f"mle{j+1}" for j in range(len(b_all))],
        colnames=colnames,
    )
    _store_matrix(Matrix, V_name, V, rownames=colnames, colnames=colnames)

    # ---- Post into e() ----
    run = SFIToolkit.stata
    run("ereturn clear")
    run(f"ereturn matrix b = {b_name}")
    run(f"ereturn matrix V = {V_name}")
    run(f"ereturn matrix b_all = {ball_name}")

    # Basic scalars / locals
    n = int(res["summary"]["inputs"]["n"])
    run(f"ereturn scalar N = {n}")
    run(f"ereturn scalar level = {float(level)}")
    run(f'ereturn local output "{out_mode}"')
    run('ereturn local cmd "dbmle"')
    run('ereturn local cmdline "dbmle_to_eclass (via Stata python)"')

    # ---- Global smallest credible sets (ONLY if not approx) ----
    if out_mode != "approx":
        # These are the GLOBAL SCS (not Frechet-conditional).
        scs = res.get("global_95_scs", {})
        intervals = scs.get("intervals", {})
        union_str = scs.get("union_str", {})

        # Map internal keys to A,C,D,N
        key_map = {"theta11": "A", "theta10": "C", "theta01": "D", "theta00": "N"}

        for key, lbl in key_map.items():
            ivs: List[Tuple[int, int]] = list(intervals.get(key, []))
            mat = _intervals_to_mat(ivs)

            # Only post matrix if non-empty
            if mat:
                scs_name = f"{prefix}__scs_{lbl}"
                _store_matrix(Matrix, scs_name, mat, colnames=["lo", "hi"])
                run(f"ereturn matrix scs_{lbl} = {scs_name}")

            s = str(union_str.get(key, "EMPTY"))
            run(f'ereturn local scs_{lbl}_str "{s}"')

    return res
