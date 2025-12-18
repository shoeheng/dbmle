# src/dbmle/stata_eclass.py

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from dbmle.core import dbmle as _dbmle


def _require_stata() -> None:
    try:
        import sfi  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "Run inside Stata's embedded Python (within a `python:` ... `end` block)."
        ) from e


def _stata(cmd: str) -> None:
    from sfi import SFIToolkit
    SFIToolkit.stata(cmd)


def _matrix_drop(name: str) -> None:
    _stata(f"capture matrix drop {name}")


def _global_set(name: str, value: str) -> None:
    v = str(value).replace('"', '""')
    _stata(f'global {name} "{v}"')


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
    // Optional: scsmat (4x2), and globals holding union strings
    syntax , BMAT(name) VMAT(name) N(real) LEVEL(real) KMLE(integer) ///
            [ SCSMAT(name) CMD(string) TITLE(string) ]

    ereturn clear
    // register b and V
    ereturn post `bmat' `vmat'

    // standard fields
    ereturn scalar N     = `n'
    ereturn scalar level = `level'
    ereturn scalar k_mle = `kmle'

    // metadata
    if ("`cmd'"=="")   ereturn local cmd   "dbmle"
    else              ereturn local cmd   "`cmd'"
    if ("`title'"=="") ereturn local title "Design-based MLE counts (Always/Complier/Defier/Never)"
    else               ereturn local title "`title'"
    ereturn local properties "b V"

    // credible set if provided
    if ("`scsmat'" != "") {
        ereturn matrix scs_minmax = `scsmat'
        ereturn local scs_theta11 "${_DBMLE_SCS_THETA11}"
        ereturn local scs_theta10 "${_DBMLE_SCS_THETA10}"
        ereturn local scs_theta01 "${_DBMLE_SCS_THETA01}"
        ereturn local scs_theta00 "${_DBMLE_SCS_THETA00}"
    }
end
""")


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
    Compute dbmle in Python; populate Stata e() via an in-memory eclass program.

    output="approx"  -> only MLE in e(b) (no SCS fields)
    output!= "approx"-> also populate SCS fields if present

    prefix snapshots (optional) create:
      matrices: <prefix>_b, <prefix>_V, <prefix>_scs_minmax (if applicable)
      scalars : scalar <prefix>_N, <prefix>_level, <prefix>_k_mle
      globals : ${<prefix>_scs_theta11}, ... (if applicable)
    """
    _require_stata()
    from sfi import Matrix

    pfx = _validate_prefix(prefix)
    populate_scs = (output != "approx")

    # Run dbmle
    res = _dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    mle = res.get("mle", {})
    mle_list = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("dbmle result missing mle.mle_list; cannot populate e().")

    a, c_, d, n_ = mle_list[0]
    k_mle = len(mle_list)

    inp = res.get("summary", {}).get("inputs", {})
    n_obs = float(inp.get("n", float("nan")))

    # Temp matrix names (drop every run to avoid conflicts)
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

    Matrix.create(tmp_V, 4, 4, 0)  # placeholder

    # Prepare SCS globals + matrix if needed
    scsmat_clause = ""
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

            # These globals are read by the eclass Stata program
            _global_set("_DBMLE_SCS_THETA11", get_union("theta11"))
            _global_set("_DBMLE_SCS_THETA10", get_union("theta10"))
            _global_set("_DBMLE_SCS_THETA01", get_union("theta01"))
            _global_set("_DBMLE_SCS_THETA00", get_union("theta00"))

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

    # Define eclass program and invoke it (this is what allows setting e())
    _ensure_sete_program()

    cmd_esc = cmd_name.replace('"', '""')
    title_esc = title.replace('"', '""')

    _stata(
        f'_dbmle_py_sete, bmat({tmp_b}) vmat({tmp_V}) '
        f'n({n_obs}) level({float(level)}) kmle({int(k_mle)})'
        f'{scsmat_clause} cmd("{cmd_esc}") title("{title_esc}")'
    )

    # Optional snapshots (outside e())
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

        _stata(f"capture scalar drop {pfx}_N")
        _stata(f"capture scalar drop {pfx}_level")
        _stata(f"capture scalar drop {pfx}_k_mle")
        _stata(f"scalar {pfx}_N = {n_obs}")
        _stata(f"scalar {pfx}_level = {float(level)}")
        _stata(f"scalar {pfx}_k_mle = {int(k_mle)}")

        if populate_scs and scsmat_clause:
            _matrix_drop(f"{pfx}_scs_minmax")
            Matrix.create(f"{pfx}_scs_minmax", 4, 2, 0)
            for r in range(4):
                for c in range(2):
                    Matrix.storeAt(f"{pfx}_scs_minmax", r, c, Matrix.getAt(tmp_scs, r, c))

            _global_set(f"{pfx}_scs_theta11", _stata('display "${_DBMLE_SCS_THETA11}"') or "")
            # Above is awkward; better: just re-set from Python strings:
            # (so do that instead)
            g = res.get("global_95_scs", {})
            intervals_obj = g.get("intervals", {})
            union_str_obj = g.get("union_str", {})

            def get_union(key: str) -> str:
                if isinstance(union_str_obj, dict) and union_str_obj.get(key):
                    return str(union_str_obj[key])
                if isinstance(intervals_obj, dict):
                    return _intervals_to_union_str(intervals_obj.get(key))
                return ""

            _global_set(f"{pfx}_scs_theta11", get_union("theta11"))
            _global_set(f"{pfx}_scs_theta10", get_union("theta10"))
            _global_set(f"{pfx}_scs_theta01", get_union("theta01"))
            _global_set(f"{pfx}_scs_theta00", get_union("theta00"))

    # Clean temp matrices
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)
    _matrix_drop(tmp_scs)

    return res if return_result else None
