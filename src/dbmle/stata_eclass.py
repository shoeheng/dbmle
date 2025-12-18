# src/dbmle/stata_eclass.py
# ---------------------------------------------------------------------
# Pure-Python Stata eclass bridge for dbmle (no .ado files shipped).
#
# How it works (important Stata constraint):
# - Stata only allows setting e() from an eclass program.
# - So Python:
#     1) computes dbmle()
#     2) writes temporary Stata matrices/globals
#     3) defines an in-memory eclass program (once per call; safe)
#     4) calls that eclass program to populate e()
#
# Behavior:
# - Always populates e(b), e(V), e(N), e(level), e(k_mle), plus e(cmd), e(title), e(properties).
# - If output == "approx": DOES NOT populate credible set objects in e().
# - Else (basic/auxiliary): populates e(scs_minmax) and e(scs_theta11/10/01/00).
# - Optional prefix snapshots (non-e()) for looping without overwriting:
#     matrices: <prefix>_b, <prefix>_V, <prefix>_scs_minmax
#     scalars : scalar <prefix>_N, <prefix>_level, <prefix>_k_mle
#     globals : ${<prefix>_scs_theta11}, ... (exact unions as strings)
# ---------------------------------------------------------------------

from __future__ import annotations

import re
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
    SFIToolkit.stata(cmd)


def _matrix_drop(name: str) -> None:
    _stata(f"capture matrix drop {name}")


def _stata_compound_quote(s: str) -> str:
    """
    Return a Stata compound-quoted string literal: `"...""'
    We also strip CR/LF to prevent Stata 'invalid syntax' failures in command parsing.
    """
    s = "" if s is None else str(s)
    s = s.replace("\r", " ").replace("\n", " ")
    s = s.replace('"', '""')
    return f'`"{s}"\''


def _ereturn_local(key: str, value: str) -> None:
    lit = _stata_compound_quote(value)
    _stata(f"ereturn local {key} {lit}")


def _global_set(name: str, value: str) -> None:
    """
    Set a Stata global macro safely (no Stata command parsing / quoting).
    """
    from sfi import Macro
    v = "" if value is None else str(value)
    v = v.replace("\r", " ").replace("\n", " ")
    Macro.setGlobal(name, v)


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
# In-memory eclass setter program
# ----------------------------
def _ensure_sete_program() -> None:
    """
    Define a tiny Stata eclass program in memory (no .ado).
    Safe to call repeatedly.
    """
    _stata("capture program drop _dbmle_py_sete")
    _stata(r"""
program define _dbmle_py_sete, eclass
    version 17
    // Required: bmat, vmat, N, level, kmle
    // Optional: scsmat, and globals holding union strings
    syntax , BMAT(name) VMAT(name) N(real) LEVEL(real) KMLE(integer) ///
            [ SCSMAT(name) CMD(string) TITLE(string) ]

    ereturn clear
    // Register b and V
    ereturn post `bmat' `vmat'

    // Standard scalars
    ereturn scalar N     = `n'
    ereturn scalar level = `level'
    ereturn scalar k_mle = `kmle'

    // Metadata
    if ("`cmd'"=="")   ereturn local cmd   "dbmle"
    else              ereturn local cmd   "`cmd'"
    if ("`title'"=="") ereturn local title "Design-based MLE counts (Always/Complier/Defier/Never)"
    else               ereturn local title "`title'"
    ereturn local properties "b V"

    // Credible set if provided
    if ("`scsmat'" != "") {
        ereturn matrix scs_minmax = `scsmat'
        ereturn local scs_theta11 "${_DBMLE_SCS_THETA11}"
        ereturn local scs_theta10 "${_DBMLE_SCS_THETA10}"
        ereturn local scs_theta01 "${_DBMLE_SCS_THETA01}"
        ereturn local scs_theta00 "${_DBMLE_SCS_THETA00}"
    }
end
""")


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
    Compute dbmle in Python; populate Stata e() using an in-memory eclass program.

    Rules:
    - output="approx"  -> ONLY MLE is populated (no SCS in e()).
    - output!="approx" -> populate SCS in e() if available.

    Note:
    - Stata has only one active e() at a time; each call overwrites e().
    - Use `estimates store name` in Stata to keep multiple e() results.
    - prefix snapshots are an extra convenience for loops (non-e()).
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

    # Temporary Stata matrix names (drop to avoid conflicts)
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

    # Prepare SCS (if requested and available)
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

            # These globals are read by the eclass Stata program
            _global_set("_DBMLE_SCS_THETA11", union_cache["theta11"])
            _global_set("_DBMLE_SCS_THETA10", union_cache["theta10"])
            _global_set("_DBMLE_SCS_THETA01", union_cache["theta01"])
            _global_set("_DBMLE_SCS_THETA00", union_cache["theta00"])

            # Numeric envelope matrix for easy Stata use
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

    # Define and call the eclass program
    _ensure_sete_program()
    cmd_lit = _stata_compound_quote(cmd_name)
    title_lit = _stata_compound_quote(title)

    _stata(
        f"_dbmle_py_sete, bmat({tmp_b}) vmat({tmp_V}) "
        f"n({n_obs}) level({float(level)}) kmle({int(k_mle)})"
        f"{scsmat_clause} cmd({cmd_lit}) title({title_lit})"
    )

    # Optional snapshots (non-e())
    if store_snapshots and pfx:
        # matrices
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

        # scalars
        _stata(f"capture scalar drop {pfx}_N")
        _stata(f"capture scalar drop {pfx}_level")
        _stata(f"capture scalar drop {pfx}_k_mle")
        _stata(f"scalar {pfx}_N = {n_obs}")
        _stata(f"scalar {pfx}_level = {float(level)}")
        _stata(f"scalar {pfx}_k_mle = {int(k_mle)}")

        # SCS snapshots only if included
        if populate_scs and scsmat_clause:
            _matrix_drop(f"{pfx}_scs_minmax")
            Matrix.create(f"{pfx}_scs_minmax", 4, 2, 0)
            for r in range(4):
                for c in range(2):
                    Matrix.storeAt(f"{pfx}_scs_minmax", r, c, Matrix.getAt(tmp_scs, r, c))

            # exact unions as globals
            _global_set(f"{pfx}_scs_theta11", union_cache.get("theta11", ""))
            _global_set(f"{pfx}_scs_theta10", union_cache.get("theta10", ""))
            _global_set(f"{pfx}_scs_theta01", union_cache.get("theta01", ""))
            _global_set(f"{pfx}_scs_theta00", union_cache.get("theta00", ""))

    # Clean up temp matrices
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)
    _matrix_drop(tmp_scs)

    return res if return_result else None
