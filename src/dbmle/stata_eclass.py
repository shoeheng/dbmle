# src/dbmle/stata_eclass.py

from __future__ import annotations
from typing import Any, Dict, List, Tuple, Optional
from dbmle.core import dbmle as _dbmle


def _require_stata() -> None:
    try:
        import sfi  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "Run inside Stata's embedded Python (within a `python:` ... `end` block)."
        ) from e


def _intervals_to_union_str(intervals: Optional[List[List[int]]]) -> str:
    if not intervals:
        return ""
    parts = [f"[{int(lo)},{int(hi)}]" for lo, hi in intervals]
    return " U ".join(parts)


def _intervals_minmax(intervals: Optional[List[List[int]]]) -> Tuple[Optional[int], Optional[int]]:
    if not intervals:
        return None, None
    lo = min(int(seg[0]) for seg in intervals)
    hi = max(int(seg[1]) for seg in intervals)
    return lo, hi


def _stata_missing() -> float:
    return float("nan")


def _ereturn_local(key: str, value: str) -> None:
    from sfi import SFIToolkit
    v = value.replace('"', '""')
    SFIToolkit.stata(f'ereturn local {key} "{v}"')
    
def _matrix_drop(name: str) -> None:
    from sfi import SFIToolkit
    SFIToolkit.stata(f"capture matrix drop {name}")


def _hard_ereturn_clear() -> None:
    """
    Clear e() and also drop estimation matrices that can remain in the matrix namespace.
    This prevents `ereturn post` from throwing 'name conflict' on repeated calls.
    """
    from sfi import SFIToolkit
    # Drop common e() matrices first (they can persist as matrices even after ereturn clear)
    SFIToolkit.stata("capture matrix drop e(b)")
    SFIToolkit.stata("capture matrix drop e(V)")
    SFIToolkit.stata("capture matrix drop e(mle_list)")
    SFIToolkit.stata("capture matrix drop e(scs_minmax)")
    # Now clear e()
    SFIToolkit.stata("ereturn clear")


def _global_macro_set(name: str, value: str) -> None:
    """
    Store a persistent string as a GLOBAL macro so it survives outside the python: block
    and can be uniquely named with a prefix in loops.
    """
    from sfi import SFIToolkit
    v = value.replace('"', '""')
    SFIToolkit.stata(f'global {name} "{v}"')


def _ereturn_clear() -> None:
    from sfi import SFIToolkit
    SFIToolkit.stata("ereturn clear")


def _ereturn_post(b_name: str, V_name: str) -> None:
    from sfi import SFIToolkit
    SFIToolkit.stata(f"ereturn post {b_name} {V_name}")


def _validate_prefix(prefix: str) -> str:
    """
    Prefix is used to create Stata object names like <prefix>_b, <prefix>_N, etc.
    Keep it conservative: letters/digits/underscore only, must start with letter/underscore.
    """
    import re
    if prefix is None:
        return ""
    prefix = str(prefix).strip()
    if prefix == "":
        return ""
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", prefix):
        raise ValueError(
            f"Invalid prefix '{prefix}'. Use only letters/digits/underscore and start with a letter or underscore."
        )
    return prefix


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
    Run dbmle() and populate Stata e() results (pure Python, no .ado).

    e() always receives the latest estimation:
      - e(b): 1x4 MLE counts [always, complier, defier, never]
      - e(V): 4x4 zero placeholder
      - e(N), e(level), e(k_mle)
      - e(mle_list) if multiple tied MLEs
      - e(scs_minmax) and e(scs_theta**) ONLY if output != "approx"

    Optional snapshots (to avoid overwriting when looping):
      - matrices: <prefix>_b, <prefix>_V, <prefix>_mle_list, <prefix>_scs_minmax
      - scalars : <prefix>_N, <prefix>_level, <prefix>_k_mle
      - globals : ${<prefix>_scs_theta11}, ${...theta10}, ${...theta01}, ${...theta00}

    Notes:
      - prefix does NOT create multiple e() namespaces (Stata doesn’t support that).
      - For native Stata workflows, prefer `estimates store <name>` after each run.
    """
    _require_stata()
    from sfi import Matrix, Scalar

    pfx = _validate_prefix(prefix)
    populate_scs = (output != "approx")

    # Compute in Python
    res = _dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    # Reset e()
    _hard_ereturn_clear()

    # Pull MLE list
    mle = res.get("mle", {})
    mle_list: List[Tuple[int, int, int, int]] = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("dbmle result missing mle.mle_list; cannot populate e().")

    a, c_, d, n_ = mle_list[0]

    # Build temporary matrices, then ereturn post to create e(b), e(V)
    _matrix_drop("DBMLE_b")
    _matrix_drop("DBMLE_V")
    
    Matrix.create("DBMLE_b", 1, 4, 0)
    try:
        Matrix.setColNames("DBMLE_b", ["always", "complier", "defier", "never"])
    except Exception:
        pass
    Matrix.storeAt("DBMLE_b", 0, 0, float(a))
    Matrix.storeAt("DBMLE_b", 0, 1, float(c_))
    Matrix.storeAt("DBMLE_b", 0, 2, float(d))
    Matrix.storeAt("DBMLE_b", 0, 3, float(n_))

    Matrix.create("DBMLE_V", 4, 4, 0)  # placeholder zeros

    _ereturn_post("DBMLE_b", "DBMLE_V")

    # Scalars
    inp = res.get("summary", {}).get("inputs", {})
    n_obs = inp.get("n", None)
    if n_obs is not None:
        Scalar.setValue("e(N)", float(n_obs))
    Scalar.setValue("e(level)", float(level))
    Scalar.setValue("e(k_mle)", float(len(mle_list)))

    # Locals (macros)
    _ereturn_local("cmd", cmd_name)
    _ereturn_local("title", title)
    _ereturn_local("properties", "b V")

    # e(mle_list) if ties
    if len(mle_list) > 1:
        Matrix.create("e(mle_list)", len(mle_list), 4, 0)
        try:
            Matrix.setColNames("e(mle_list)", ["always", "complier", "defier", "never"])
        except Exception:
            pass
        for i, (aa, cc, dd, nn) in enumerate(mle_list):
            Matrix.storeAt("e(mle_list)", i, 0, float(aa))
            Matrix.storeAt("e(mle_list)", i, 1, float(cc))
            Matrix.storeAt("e(mle_list)", i, 2, float(dd))
            Matrix.storeAt("e(mle_list)", i, 3, float(nn))

    # SCS only if NOT approx
    if populate_scs:
        g = res.get("global_95_scs", None)
        if g and isinstance(g, dict):
            intervals_obj = g.get("intervals", None)   # theta11/theta10/theta01/theta00 -> segments
            union_str_obj = g.get("union_str", None)   # optional preformatted strings

            # Exact unions as e()-locals (don’t lose disjointness)
            def get_union(key: str) -> str:
                if isinstance(union_str_obj, dict) and key in union_str_obj and union_str_obj[key]:
                    return str(union_str_obj[key])
                if isinstance(intervals_obj, dict):
                    return _intervals_to_union_str(intervals_obj.get(key))
                return ""

            _ereturn_local("scs_theta11", get_union("theta11"))
            _ereturn_local("scs_theta10", get_union("theta10"))
            _ereturn_local("scs_theta01", get_union("theta01"))
            _ereturn_local("scs_theta00", get_union("theta00"))

            # Min/max envelope matrix for easy tabulation
            if isinstance(intervals_obj, dict):
                Matrix.create("e(scs_minmax)", 4, 2, 0)
                try:
                    Matrix.setRowNames("e(scs_minmax)", ["always", "complier", "defier", "never"])
                    Matrix.setColNames("e(scs_minmax)", ["min", "max"])
                except Exception:
                    pass

                keys = ["theta11", "theta10", "theta01", "theta00"]
                for r, key in enumerate(keys):
                    lo, hi = _intervals_minmax(intervals_obj.get(key))
                    if lo is None or hi is None:
                        Matrix.storeAt("e(scs_minmax)", r, 0, _stata_missing())
                        Matrix.storeAt("e(scs_minmax)", r, 1, _stata_missing())
                    else:
                        Matrix.storeAt("e(scs_minmax)", r, 0, float(lo))
                        Matrix.storeAt("e(scs_minmax)", r, 1, float(hi))

    # Optional snapshots (prefix avoids overwriting snapshots across loop iterations)
    if store_snapshots and pfx:
        # Matrices
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

        if len(mle_list) > 1:
            Matrix.create(f"{pfx}_mle_list", len(mle_list), 4, 0)
            try:
                Matrix.setColNames(f"{pfx}_mle_list", ["always", "complier", "defier", "never"])
            except Exception:
                pass
            for i, (aa, cc, dd, nn) in enumerate(mle_list):
                Matrix.storeAt(f"{pfx}_mle_list", i, 0, float(aa))
                Matrix.storeAt(f"{pfx}_mle_list", i, 1, float(cc))
                Matrix.storeAt(f"{pfx}_mle_list", i, 2, float(dd))
                Matrix.storeAt(f"{pfx}_mle_list", i, 3, float(nn))

        # Scalars
        if n_obs is not None:
            Scalar.setValue(f"{pfx}_N", float(n_obs))
        Scalar.setValue(f"{pfx}_level", float(level))
        Scalar.setValue(f"{pfx}_k_mle", float(len(mle_list)))

        # SCS snapshots only if available (i.e., not approx)
        if populate_scs:
            # Copy the e()-locals into GLOBAL macros keyed by prefix
            from sfi import SFIToolkit
            # pull e()-locals by asking Stata to expand them
            def _expand_e_local(name: str) -> str:
                # Stata will print the expanded macro; capture via `return` isn't easy here,
                # so we just rebuild the same unions again from res (more robust).
                return ""

            g = res.get("global_95_scs", None)
            if g and isinstance(g, dict):
                intervals_obj = g.get("intervals", None)
                union_str_obj = g.get("union_str", None)

                def get_union(key: str) -> str:
                    if isinstance(union_str_obj, dict) and key in union_str_obj and union_str_obj[key]:
                        return str(union_str_obj[key])
                    if isinstance(intervals_obj, dict):
                        return _intervals_to_union_str(intervals_obj.get(key))
                    return ""

                _global_macro_set(f"{pfx}_scs_theta11", get_union("theta11"))
                _global_macro_set(f"{pfx}_scs_theta10", get_union("theta10"))
                _global_macro_set(f"{pfx}_scs_theta01", get_union("theta01"))
                _global_macro_set(f"{pfx}_scs_theta00", get_union("theta00"))

                if isinstance(intervals_obj, dict):
                    Matrix.create(f"{pfx}_scs_minmax", 4, 2, 0)
                    try:
                        Matrix.setRowNames(f"{pfx}_scs_minmax", ["always", "complier", "defier", "never"])
                        Matrix.setColNames(f"{pfx}_scs_minmax", ["min", "max"])
                    except Exception:
                        pass
                    keys = ["theta11", "theta10", "theta01", "theta00"]
                    for r, key in enumerate(keys):
                        lo, hi = _intervals_minmax(intervals_obj.get(key))
                        if lo is None or hi is None:
                            Matrix.storeAt(f"{pfx}_scs_minmax", r, 0, _stata_missing())
                            Matrix.storeAt(f"{pfx}_scs_minmax", r, 1, _stata_missing())
                        else:
                            Matrix.storeAt(f"{pfx}_scs_minmax", r, 0, float(lo))
                            Matrix.storeAt(f"{pfx}_scs_minmax", r, 1, float(hi))

    return res if return_result else None
