# src/dbmle/stata_eclass.py
# ---------------------------------------------------------------------
# Pure-Python Stata eclass bridge for dbmle (no shipped .ado).
#
# Key constraint:
#   Stata forbids setting e() unless you're inside an eclass program.
#
# Robust approach:
#   - Python computes dbmle()
#   - Python creates temp matrices + globals
#   - Python ensures an in-memory eclass program exists by writing a
#     temporary .do file and `do`-ing it (runtime only; not shipped)
#   - Python calls the eclass program to populate e()
#
# Behavior:
#   - Always populates: e(b), e(V), e(N), e(level), e(k_mle), e(cmd), e(title)
#   - If output == "approx": DOES NOT populate SCS objects in e()
#   - Else: populates e(scs_minmax) and e(scs_theta11/10/01/00)
#   - Optional snapshots with prefix (outside e()) for looping
# ---------------------------------------------------------------------

from __future__ import annotations

import os
import re
import tempfile
from typing import Any, Dict, List, Optional, Tuple

from dbmle.core import dbmle as _dbmle


# ----------------------------
# Stata execution helpers
# ----------------------------
def _require_stata() -> None:
    try:
        import sfi  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "Run inside Stata's embedded Python, e.g.\n"
            "    python:\n"
            "        from dbmle import dbmle_to_eclass\n"
            "        dbmle_to_eclass(...)\n"
            "    end"
        ) from e


def _stata(cmd: str) -> None:
    from sfi import SFIToolkit
    # SFIToolkit.stata expects a single valid Stata command
    SFIToolkit.stata(cmd)


def _matrix_drop(name: str) -> None:
    _stata(f"capture matrix drop {name}")


def _scalar_drop(name: str) -> None:
    _stata(f"capture scalar drop {name}")


def _validate_prefix(prefix: str) -> str:
    if prefix is None:
        return ""
    prefix = str(prefix).strip()
    if prefix == "":
        return ""
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", prefix):
        raise ValueError(
            f"Invalid prefix '{prefix}'. Use letters/digits/underscore; start with letter/underscore."
        )
    return prefix


def _global_set(name: str, value: str) -> None:
    """
    Set a Stata global macro safely (bypasses Stata command parsing).
    """
    from sfi import Macro
    v = "" if value is None else str(value)
    v = v.replace("\r", " ").replace("\n", " ")
    Macro.setGlobal(name, v)


# ----------------------------
# Credible set formatting
# ----------------------------
def _intervals_to_union_str(intervals: Optional[List[List[int]]]) -> str:
    if not intervals:
        return ""
    return " U ".join(f"[{int(lo)},{int(hi)}]" for lo, hi in intervals)


def _intervals_minmax(intervals: Optional[List[List[int]]]) -> Tuple[Optional[int], Optional[int]]:
    if not intervals:
        return None, None
    lo = min(int(seg[0]) for seg in intervals)
    hi = max(int(seg[1]) for seg in intervals)
    return lo, hi


# ----------------------------
# eclass program loader (runtime-only .do)
# ----------------------------
_PROG_LOADED_FLAG_GLOBAL = "_DBMLE_PY_SETE_LOADED"


def _ensure_sete_program() -> None:
    """
    Ensure the Stata eclass program _dbmle_py_sete exists.

    Robust method: write a tiny .do file and `do` it.
    This avoids all the "1>>>" program-entry-mode weirdness.
    """
    from sfi import Macro

    if Macro.getGlobal(_PROG_LOADED_FLAG_GLOBAL) == "1":
        return

    # Minimal eclass setter program: reads temp matrices/globals created by Python
    do_text = r"""
capture program drop _dbmle_py_sete
program define _dbmle_py_sete, eclass
    version 17
    // Required
    syntax , BMAT(name) VMAT(name) N(real) LEVEL(real) KMLE(integer) ///
            [ SCSMAT(name) CMD(string) TITLE(string) ]

    ereturn clear
    ereturn post `bmat' `vmat'

    ereturn scalar N     = `n'
    ereturn scalar level = `level'
    ereturn scalar k_mle = `kmle'

    if ("`cmd'"=="")   ereturn local cmd   "dbmle"
    else              ereturn local cmd   "`cmd'"

    if ("`title'"=="") ereturn local title "Design-based MLE counts (Always/Complier/Defier/Never)"
    else               ereturn local title "`title'"

    ereturn local properties "b V"

    if ("`scsmat'" != "") {
        ereturn matrix scs_minmax = `scsmat'
        ereturn local scs_theta11 "${_DBMLE_SCS_THETA11}"
        ereturn local scs_theta10 "${_DBMLE_SCS_THETA10}"
        ereturn local scs_theta01 "${_DBMLE_SCS_THETA01}"
        ereturn local scs_theta00 "${_DBMLE_SCS_THETA00}"
    }
end
"""

    # Write & run the do-file
    fd, path = tempfile.mkstemp(prefix="dbmle_sete_", suffix=".do")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(do_text)
        # Use quotes around path (may include spaces)
        _stata(f'do "{path}"')
        Macro.setGlobal(_PROG_LOADED_FLAG_GLOBAL, "1")
    finally:
        try:
            os.remove(path)
        except OSError:
            pass


# ----------------------------
# Public API
# ----------------------------
def dbmle_to_eclass(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    *,
    output: str = "basic",
    level: float = 0.95,
    show_progress: bool = True,
    prefix: str = "",
    store_snapshots: bool = True,
    cmd_name: str = "dbmle",
    title: str = "Design-based MLE counts (Always/Complier/Defier/Never)",
    return_result: bool = False,
) -> Optional[Dict[str, Any]]:
    """
    Compute dbmle in Python; populate Stata e() using an eclass program.

    - output="approx": ONLY MLE fields in e()
    - otherwise: also populate credible set fields in e() if available
    """
    _require_stata()
    from sfi import Matrix

    pfx = _validate_prefix(prefix)
    populate_scs = (output != "approx")

    # Compute in Python
    res = _dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    # Pull MLE
    mle = res.get("mle", {})
    mle_list = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("dbmle result missing mle.mle_list; cannot populate e().")

    a, c_, d, n_ = mle_list[0]
    k_mle = len(mle_list)

    inp = res.get("summary", {}).get("inputs", {})
    n_obs = float(inp.get("n", float("nan")))

    # Temp matrix names (drop every run)
    tmp_b = "__dbmle_b"
    tmp_V = "__dbmle_V"
    tmp_scs = "__dbmle_scs"
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)
    _matrix_drop(tmp_scs)

    # Build b and V
    Matrix.create(tmp_b, 1, 4, 0)
    try:
        Matrix.setColNames(tmp_b, ["always", "complier", "defier", "never"])
    except Exception:
        pass
    Matrix.storeAt(tmp_b, 0, 0, float(a))
    Matrix.storeAt(tmp_b, 0, 1, float(c_))
    Matrix.storeAt(tmp_b, 0, 2, float(d))
    Matrix.storeAt(tmp_b, 0, 3, float(n_))

    Matrix.create(tmp_V, 4, 4, 0)  # placeholder zeros

    # SCS prep (only if requested)
    scsmat_clause = ""
    union_cache: Dict[str, str] = {}
    if populate_scs:
        g = res.get("global_95_scs", None)
        if isinstance(g, dict):
            intervals_obj = g.get("intervals", None)
            union_str_obj = g.get("union_str", None)

            def get_union(key: str) -> str:
                if isinstance(union_str_obj, dict) and union_str_obj.get(key):
                    return str(union_str_obj[key])
                if isinstance(intervals_obj, dict):
                    return _intervals_to_union_str(intervals_obj.get(key))
                return ""

            union_cache = {
                "theta11": get_union("theta11"),
                "theta10": get_union("theta10"),
                "theta01": get_union("theta01"),
                "theta00": get_union("theta00"),
            }

            # globals consumed by the eclass program
            _global_set("_DBMLE_SCS_THETA11", union_cache["theta11"])
            _global_set("_DBMLE_SCS_THETA10", union_cache["theta10"])
            _global_set("_DBMLE_SCS_THETA01", union_cache["theta01"])
            _global_set("_DBMLE_SCS_THETA00", union_cache["theta00"])

            # numeric envelope matrix (min/max per type)
            if isinstance(intervals_obj, dict):
                Matrix.create(tmp_scs, 4, 2, 0)
                try:
                    Matrix.setRowNames(tmp_scs, ["always", "complier", "defier", "never"])
                    Matrix.setColNames(tmp_scs, ["min", "max"])
                except Exception:
                    pass
                keys = ["theta11", "theta10", "theta01", "theta00"]
                for r, key in enumerate(keys):
                    lo, hi = _intervals_minmax(intervals_obj.get(key))
                    if lo is None or hi is None:
                        Matrix.storeAt(tmp_scs, r, 0, float("nan"))
                        Matrix.storeAt(tmp_scs, r, 1, float("nan"))
                    else:
                        Matrix.storeAt(tmp_scs, r, 0, float(lo))
                        Matrix.storeAt(tmp_scs, r, 1, float(hi))
                scsmat_clause = f" scsmat({tmp_scs})"

    # Ensure eclass program exists, then call it
    _ensure_sete_program()

    # Note: cmd() and title() are parsed by Stata; keep them simple.
    # We'll avoid complicated quoting: strip newlines and double quotes.
    cmd_clean = (cmd_name or "dbmle").replace("\r", " ").replace("\n", " ").replace('"', "''")
    title_clean = (title or "").replace("\r", " ").replace("\n", " ").replace('"', "''")

    call = (
        f'_dbmle_py_sete, bmat({tmp_b}) vmat({tmp_V}) '
        f'n({n_obs}) level({float(level)}) kmle({int(k_mle)})'
        f'{scsmat_clause} cmd("{cmd_clean}") title("{title_clean}")'
    )
    _stata(call)

    # Snapshots (outside e()) if prefix requested
    if store_snapshots and pfx:
        _matrix_drop(f"{pfx}_b")
        _matrix_drop(f"{pfx}_V")
        Matrix.create(f"{pfx}_b", 1, 4, 0)
        try:
            Matrix.setColNames(f"{pfx}_b", ["always", "complier", "defier", "never"])
        except Exception:
            pass
        Matrix.storeAt(f"{pfx}_b", 0, 0, float(a))
        Matrix.storeAt(f"{pfx}_b", 0, 1, float(c_))
        Matrix.storeAt(f"{pfx}_b", 0, 2, float(d))
        Matrix.storeAt(f"{pfx}_b", 0, 3, float(n_))
        Matrix.create(f"{pfx}_V", 4, 4, 0)

        _scalar_drop(f"{pfx}_N")
        _scalar_drop(f"{pfx}_level")
        _scalar_drop(f"{pfx}_k_mle")
        _stata(f"scalar {pfx}_N = {n_obs}")
        _stata(f"scalar {pfx}_level = {float(level)}")
        _stata(f"scalar {pfx}_k_mle = {int(k_mle)}")

        if populate_scs and scsmat_clause:
            _matrix_drop(f"{pfx}_scs_minmax")
            Matrix.create(f"{pfx}_scs_minmax", 4, 2, 0)
            for r in range(4):
                for c in range(2):
                    Matrix.storeAt(f"{pfx}_scs_minmax", r, c, Matrix.getAt(tmp_scs, r, c))

            _global_set(f"{pfx}_scs_theta11", union_cache.get("theta11", ""))
            _global_set(f"{pfx}_scs_theta10", union_cache.get("theta10", ""))
            _global_set(f"{pfx}_scs_theta01", union_cache.get("theta01", ""))
            _global_set(f"{pfx}_scs_theta00", union_cache.get("theta00", ""))

    # Clean up temp matrices
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)
    _matrix_drop(tmp_scs)

    return res if return_result else None
