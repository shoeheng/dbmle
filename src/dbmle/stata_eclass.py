# src/dbmle/stata_eclass.py

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from dbmle.core import dbmle as _dbmle


# ----------------------------
# small Stata helpers
# ----------------------------
def _require_stata() -> None:
    try:
        import sfi  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "Run this inside Stata's embedded Python:\n"
            "    python:\n"
            "        ...\n"
            "    end"
        ) from e


def _stata(cmd: str) -> None:
    from sfi import SFIToolkit
    SFIToolkit.stata(cmd)


def _matrix_drop(name: str) -> None:
    # Works for normal matrices; for e(b) style names use quotes in capture if needed
    _stata(f"capture matrix drop {name}")


def _ereturn_clear_hard() -> None:
    """
    Clears e() AND drops common e() matrices that sometimes linger as matrices.
    """
    _stata("capture ereturn clear")
    # Defensive: drop both named and e()-namespaced matrices
    _stata("capture matrix drop e(b)")
    _stata("capture matrix drop e(V)")
    _stata("capture matrix drop e(mle_list)")
    _stata("capture matrix drop e(scs_minmax)")


def _ereturn_local(key: str, value: str) -> None:
    v = value.replace('"', '""')
    _stata(f'ereturn local {key} "{v}"')


def _ereturn_scalar(key: str, value: float) -> None:
    # key should be bare name; Stata stores as e(<key>)
    _stata(f"ereturn scalar {key} = {value}")


def _ereturn_matrix(key: str, source_matrix_name: str) -> None:
    """
    Store matrix into e(): ereturn matrix <key> = <source>
    Example: key="b" source="__dbmle_b123"
    """
    _stata(f"ereturn matrix {key} = {source_matrix_name}")


def _global_set(name: str, value: str) -> None:
    v = value.replace('"', '""')
    _stata(f'global {name} "{v}"')


def _validate_prefix(prefix: str) -> str:
    if prefix is None:
        return ""
    prefix = str(prefix).strip()
    if prefix == "":
        return ""
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", prefix):
        raise ValueError(
            f"Invalid prefix '{prefix}'. Use letters/digits/underscore; start with a letter/underscore."
        )
    return prefix


# ----------------------------
# credible set formatting
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
    Pure-Python Stata bridge: writes dbmle results into e().

    - Always populates:
        e(b) : 1x4 MLE counts [always, complier, defier, never]
        e(V) : 4x4 zero placeholder
        e(N), e(level), e(k_mle)
        e(cmd), e(title), e(properties)

    - Populates SCS only when output != "approx":
        e(scs_minmax) as 4x2 (min,max envelope)
        e(scs_theta11/10/01/00) as e()-locals (exact unions as strings)

    - Optional snapshots when prefix is given:
        matrices: <prefix>_b, <prefix>_V, <prefix>_scs_minmax (if applicable)
        scalars : <prefix>_N, <prefix>_level, <prefix>_k_mle
        globals : ${<prefix>_scs_theta11}, ... (if applicable)

    Note: e() always holds only the *latest* estimation (Stata design).
    """
    _require_stata()
    from sfi import Matrix  # only used to build numeric matrices safely

    pfx = _validate_prefix(prefix)
    populate_scs = (output != "approx")

    # Run dbmle in Python
    res = _dbmle(
        xI1, xI0, xC1, xC0,
        output=output,
        level=level,
        show_progress=show_progress,
    )

    # Hard clear e()
    _ereturn_clear_hard()

    # Pull MLE list
    mle = res.get("mle", {})
    mle_list = mle.get("mle_list", [])
    if not mle_list:
        raise RuntimeError("dbmle result missing mle.mle_list; cannot populate e().")

    a, c_, d, n_ = mle_list[0]
    k_mle = len(mle_list)

    # Build temporary matrices with unique-ish names and drop first
    tmp_b = "__dbmle_b"
    tmp_V = "__dbmle_V"
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)

    # temp b
    Matrix.create(tmp_b, 1, 4, 0)
    try:
        Matrix.setColNames(tmp_b, ["always", "complier", "defier", "never"])
    except Exception:
        pass
    Matrix.storeAt(tmp_b, 0, 0, float(a))
    Matrix.storeAt(tmp_b, 0, 1, float(c_))
    Matrix.storeAt(tmp_b, 0, 2, float(d))
    Matrix.storeAt(tmp_b, 0, 3, float(n_))

    # temp V (placeholder zeros)
    Matrix.create(tmp_V, 4, 4, 0)

    # Store into e() WITHOUT ereturn post
    _ereturn_matrix("b", tmp_b)
    _ereturn_matrix("V", tmp_V)

    # Scalars
    inp = res.get("summary", {}).get("inputs", {})
    n_obs = inp.get("n", None)
    if n_obs is not None:
        _ereturn_scalar("N", float(n_obs))
    _ereturn_scalar("level", float(level))
    _ereturn_scalar("k_mle", float(k_mle))

    # Locals (macros)
    _ereturn_local("cmd", cmd_name)
    _ereturn_local("title", title)
    _ereturn_local("properties", "b V")

    # Optional: store all tied MLEs as e(mle_list)
    if k_mle > 1:
        _matrix_drop("e(mle_list)")
        Matrix.create("e(mle_list)", k_mle, 4, 0)
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
        if isinstance(g, dict):
            intervals_obj = g.get("intervals", None)   # dict: theta11/theta10/theta01/theta00 -> segments
            union_str_obj = g.get("union_str", None)   # optional dict of strings

            def get_union(key: str) -> str:
                if isinstance(union_str_obj, dict) and union_str_obj.get(key):
                    return str(union_str_obj[key])
                if isinstance(intervals_obj, dict):
                    return _intervals_to_union_str(intervals_obj.get(key))
                return ""

            # exact unions as e()-locals
            _ereturn_local("scs_theta11", get_union("theta11"))
            _ereturn_local("scs_theta10", get_union("theta10"))
            _ereturn_local("scs_theta01", get_union("theta01"))
            _ereturn_local("scs_theta00", get_union("theta00"))

            # numeric envelope matrix
            if isinstance(intervals_obj, dict):
                _matrix_drop("e(scs_minmax)")
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
                        Matrix.storeAt("e(scs_minmax)", r, 0, float("nan"))
                        Matrix.storeAt("e(scs_minmax)", r, 1, float("nan"))
                    else:
                        Matrix.storeAt("e(scs_minmax)", r, 0, float(lo))
                        Matrix.storeAt("e(scs_minmax)", r, 1, float(hi))

    # Snapshots (non-e()) if prefix requested
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

        # scalars (as Stata scalars, not e())
        if n_obs is not None:
            _stata(f"scalar {pfx}_N = {float(n_obs)}")
        _stata(f"scalar {pfx}_level = {float(level)}")
        _stata(f"scalar {pfx}_k_mle = {float(k_mle)}")

        if populate_scs and isinstance(res.get("global_95_scs", None), dict):
            g = res["global_95_scs"]
            intervals_obj = g.get("intervals", None)
            union_str_obj = g.get("union_str", None)

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

            if isinstance(intervals_obj, dict):
                _matrix_drop(f"{pfx}_scs_minmax")
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
                        Matrix.storeAt(f"{pfx}_scs_minmax", r, 0, float("nan"))
                        Matrix.storeAt(f"{pfx}_scs_minmax", r, 1, float("nan"))
                    else:
                        Matrix.storeAt(f"{pfx}_scs_minmax", r, 0, float(lo))
                        Matrix.storeAt(f"{pfx}_scs_minmax", r, 1, float(hi))

    # Cleanup temp matrices (optional)
    _matrix_drop(tmp_b)
    _matrix_drop(tmp_V)

    return res if return_result else None
