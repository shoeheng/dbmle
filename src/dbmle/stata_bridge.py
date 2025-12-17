from __future__ import annotations
from typing import Any, Dict, List, Tuple, Optional

import dbmle.core as core

try:
    from sfi import Scalar, Macro, Matrix
except Exception:
    Scalar = Macro = Matrix = None


def _require_stata():
    if Scalar is None or Macro is None or Matrix is None:
        raise RuntimeError(
            "Stata bridge requires Stata's embedded Python (module 'sfi' not found). "
            "Run this from within Stata using a python: ... end block."
        )


def set_r_from_result(res: Dict[str, Any], *, prefix: str = "") -> None:
    """
    Write key results into Stata r().

    prefix: optional prefix for names, e.g. prefix="dbmle_" -> r(dbmle_n), etc.
    """
    _require_stata()

    inp = res["summary"]["inputs"]          # expects xI1,xI0,xC1,xC0,n,m,c
    mle_list = res["mle"]["mle_list"]       # list of (t11,t10,t01,t00)
    outmode = res["meta"]["output"]

    def rname(name: str) -> str:
        return f"r({prefix}{name})"

    # scalars
    Scalar.setValue(rname("n"), float(inp["n"]))
    Scalar.setValue(rname("m"), float(inp["m"]))
    Scalar.setValue(rname("c"), float(inp["c"]))

    # first MLE as scalars
    t11, t10, t01, t00 = mle_list[0]
    Scalar.setValue(rname("theta11_mle"), float(t11))
    Scalar.setValue(rname("theta10_mle"), float(t10))
    Scalar.setValue(rname("theta01_mle"), float(t01))
    Scalar.setValue(rname("theta00_mle"), float(t00))

    # all MLEs (ties) as matrix kx4
    k = len(mle_list)
    Matrix.create(f"r({prefix}mle_list)", k, 4, 0)
    
    for i, (a, c_, d, n_) in enumerate(mle_list):   # i = 0..k-1
        Matrix.storeAt(f"r({prefix}mle_list)", i, 0, float(a))
        Matrix.storeAt(f"r({prefix}mle_list)", i, 1, float(c_))
        Matrix.storeAt(f"r({prefix}mle_list)", i, 2, float(d))
        Matrix.storeAt(f"r({prefix}mle_list)", i, 3, float(n_))


    # locals (strings) for unions etc., only in exhaustive modes
    if outmode in ("basic", "auxiliary"):
        union = res["global_95_scs"]["union_str"]
        Macro.setLocal(f"{prefix}theta11_scs", union["theta11"])
        Macro.setLocal(f"{prefix}theta10_scs", union["theta10"])
        Macro.setLocal(f"{prefix}theta01_scs", union["theta01"])
        Macro.setLocal(f"{prefix}theta00_scs", union["theta00"])


    # store full printable report as a local
    if "report" in res:
        Macro.setLocal(f"{prefix}report", res["report"])



def dbmle_to_r(
    xI1: int, xI0: int, xC1: int, xC0: int,
    *, output: str = "basic", level: float = 0.95,
    show_progress: bool = True, prefix: str = ""
) -> DBMLEResult:
    """
    Compute dbmle() and immediately populate Stata's r().
    Returns the usual DBMLEResult as well.
    """
    res = core.dbmle(xI1, xI0, xC1, xC0, output=output, level=level, show_progress=show_progress)
    set_r_from_result(res, prefix=prefix)
    return res
