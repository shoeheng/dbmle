# src/dbmle/core.py
# ---------------------------------------------------------------------
# DBMLE: Design-Based Maximum Likelihood + optional supporting stats
#
# inputs (in this order): xI1, xI0, xC1, xC0
#
# method:
#   "approx"     -> fast approximation of MLE
#   "exhaustive" -> exhaustive grid search
#
# auxiliary:
#   False        -> MLE and (if exhaustive) smallest credible set
#   True         -> All auxiliary statistics
#
# Output: two text blocks:
#   1. "Standard Statistics"
#   2. "Design-Based Maximum Likelihood Estimates"
#      or "Design-Based Maximum Likelihood Estimates and Auxiliary Statistics"
#
# Requires: tqdm   (pip install tqdm)
# ---------------------------------------------------------------------

import math
import warnings
from typing import Dict, Tuple, List, Any, Optional
from tqdm import tqdm

# A "theta" here is the joint count of (always, complier, defier, never)
# in that order: (theta11, theta10, theta01, theta00)
Theta = Tuple[int, int, int, int]


# =====================================================================
# combinatorics + likelihood
# =====================================================================
def _logC_scalar(n: int, k: int) -> float:
    """log of "n choose k" that is safe for out-of-range k."""
    if 0 <= k <= n:
        return (
            math.lgamma(n + 1)
            - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)
        )
    return -math.inf


def _loglik(
    n: int,
    m: int,
    xI1: int,
    xC1: int,
    t11: int,
    t10: int,
    t01: int,
    t00: int,
) -> float:
    """
    Log-likelihood of seeing (xI1, xC1) when the true joint counts are
    (t11, t10, t01, t00) and we randomized m units into intervention.
    """
    # must be a partition of the total population
    if t11 + t10 + t01 + t00 != n:
        return -math.inf

    c = n - m
    xI0 = m - xI1
    xC0 = c - xC1

    # basic feasibility: do we have enough would-not-take types to match zeros?
    if t01 + t00 < xI0 or t10 + t00 < xC0:
        return -math.inf

    # compute inner loop bounds for how many always-takers show up in I=1
    lo = max(
        0,
        xI1 - t10,
        t11 - xC1,
        t11 + t01 + xI1 - m - xC1,
    )
    hi = min(
        t11,
        xI1,
        t11 + t01 - xC1,
        c - t10 - xC1 + xI1,
    )
    if lo > hi:
        return -math.inf

    # common log binomial term (randomization)
    log_choose = _logC_scalar(n, m)

    # we sum over all consistent allocations into the intervention arm
    logL = -math.inf
    for i in range(lo, hi + 1):
        I10 = xI1 - i
        I01 = t11 + t01 - xC1 - i
        I00 = m + xC1 + i - t11 - t01 - xI1
        if (
            I10 < 0
            or I10 > t10
            or I01 < 0
            or I01 > t01
            or I00 < 0
            or I00 > t00
        ):
            continue

        term = (
            _logC_scalar(t11, i)
            + _logC_scalar(t10, I10)
            + _logC_scalar(t01, I01)
            + _logC_scalar(t00, I00)
        )

        # log-sum-exp, but done by hand to avoid underflow
        if logL < term:
            logL, term = term, logL
        if term > -math.inf:
            logL = logL + math.log1p(math.exp(term - logL))

    return logL - log_choose


# =====================================================================
# validation
# =====================================================================
def _validate_inputs(xI1: int, xI0: int, xC1: int, xC0: int) -> None:
    """Check basic shape and nonnegativity of the raw data counts."""
    for name, v in (("xI1", xI1), ("xI0", xI0), ("xC1", xC1), ("xC0", xC0)):
        if not isinstance(v, int) or v < 0:
            raise ValueError(f"{name} must be a nonnegative int.")
    m = xI1 + xI0
    c = xC1 + xC0
    if m <= 0 or c <= 0:
        raise ValueError("Intervention and control sizes must both be positive.")


# =====================================================================
# Frechet (continuous, data-implied)
# theta11+theta10 = (n/m) * xI1
# theta11+theta01 = (n/(n-m)) * xC1
# =====================================================================
def _frechet_continuous_from_data(
    n: int, m: int, xI1: int, xC1: int
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Frechet bounds, treating the line-sum constraints as continuous.
    Returns lower/upper for each of the 4 types.
    """
    c = n - m
    A = (n / m) * xI1
    B = (n / c) * xC1

    U11 = min(A, B)
    L11 = max(0.0, A + B - n)

    U10 = min(A, n - B)
    L10 = max(0.0, A - B)

    U01 = min(B, n - A)
    L01 = max(0.0, B - A)

    U00 = min(n - A, n - B)
    L00 = max(0.0, n - A - B)

    return L11, U11, L10, U10, L01, U01, L00, U00


def _corners_float(
    n: int, m: int, xI1: int, xC1: int
) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
    """
    Pick the two opposite "corners" of the continuous Frechet region that
    are relevant for the observed data. These get integerized later.
    """
    L11, U11, L10, U10, L01, U01, L00, U00 = _frechet_continuous_from_data(
        n, m, xI1, xC1
    )
    c = n - m
    if xI1 > m / 2 and xC1 > c / 2:
        conj = (U11, L10, L01, U00)
        opp = (L11, U10, U01, L00)
    elif xI1 < m / 2 and xC1 > c / 2:
        conj = (L11, U10, U01, L00)
        opp = (U11, L10, L01, U00)
    elif xI1 < m / 2 and xC1 < c / 2:
        conj = (U11, L10, L01, U00)
        opp = (L11, U10, U01, L00)
    else:
        # xI1 >= m/2 and xC1 <= c/2
        conj = (L11, U10, U01, L00)
        opp = (U11, L10, L01, U00)

    return conj, opp


def _integerizations_sum_n(n: int, corner: Tuple[float, ...]) -> List[Theta]:
    """
    Take a continuous corner (4 floats), round up/down in all 16 ways, and
    keep the integer versions that sum exactly to n.
    """
    floors = [math.floor(v) for v in corner]
    ceils = [math.ceil(v) for v in corner]
    out: List[Theta] = []
    for mask in range(16):
        t = [ceils[i] if (mask >> i) & 1 else floors[i] for i in range(4)]
        if sum(t) == n and min(t) >= 0:
            out.append(tuple(t))  # type: ignore[arg-type]
    return sorted(set(out))


def _best_frechet_corner_point(
    n: int, m: int, xI1: int, xC1: int
) -> Theta:
    """
    Among the integerizations of the two main Frechet corners, pick the one
    that maximizes the likelihood.
    """
    conj, opp = _corners_float(n, m, xI1, xC1)
    cands = _integerizations_sum_n(n, conj) + _integerizations_sum_n(n, opp)
    if not cands:
        # fallback: round the conjugate corner and fix the last entry
        near = tuple(int(round(v)) for v in conj)
        if sum(near) != n:
            near = (near[0], near[1], near[2], n - near[0] - near[1] - near[2])
        return near  # type: ignore[return-value]

    scored = [(t, _loglik(n, m, xI1, xC1, *t)) for t in cands]
    best_ll = max(ll for _, ll in scored)
    bests = [t for t, ll in scored if abs(ll - best_ll) <= 1e-13]
    # deterministic tie-breaking
    return sorted(bests)[0]


# =====================================================================
# interval helpers
# =====================================================================
def _disjoint_intervals(vals: List[int]) -> List[Tuple[int, int]]:
    """
    Take a list like [1,2,3,7,8] and return [(1,3), (7,8)].
    """
    if not vals:
        return []
    vals = sorted(set(vals))
    start = prev = vals[0]
    out: List[Tuple[int, int]] = []
    for v in vals[1:]:
        if v == prev + 1:
            prev = v
            continue
        out.append((start, prev))
        start = prev = v
    out.append((start, prev))
    return out


def _format_union(intervals: List[Tuple[int, int]]) -> str:
    """
    Format a union of integer intervals.
    """
    if not intervals:
        return "EMPTY"
    parts = []
    for a, b in intervals:
        if a == b:
            parts.append("{" + str(a) + "}")
        else:
            parts.append("[" + str(a) + "," + str(b) + "]")
    return " U ".join(parts)


# =====================================================================
# Frechet-marginal 95% SCS (returns also per-interval denoms)
# =====================================================================
def _frechet_marginal_scs(
    n: int, m: int, xI1: int, xC1: int, level: float = 0.95
) -> Dict[str, Any]:
    """
    Build a 1D credible set along the Frechet line implied by the
    data-implied row/column sums. We weight candidate joints by likelihood,
    normalize, and then take the top mass until we hit "level".
    """
    theta_star = _best_frechet_corner_point(n, m, xI1, xC1)
    t11s, t10s, t01s, t00s = theta_star
    R = t11s + t10s
    C = t11s + t01s

    t11_lo = max(0, R + C - n)
    t11_hi = min(R, C)

    thetas: List[Theta] = []
    lls: List[float] = []
    for t11 in range(t11_lo, t11_hi + 1):
        t10 = R - t11
        t01 = C - t11
        t00 = n - R - C + t11
        if min(t10, t01, t00) < 0:
            continue
        ll = _loglik(n, m, xI1, xC1, t11, t10, t01, t00)
        if ll == -math.inf:
            continue
        thetas.append((t11, t10, t01, t00))
        lls.append(ll)

    if not thetas:
        return {
            "theta11_intervals": [],
            "theta10_intervals": [],
            "theta01_intervals": [],
            "theta00_intervals": [],
            "theta11_denoms": [],
            "theta10_denoms": [],
            "theta01_denoms": [],
            "theta00_denoms": [],
        }

    max_ll = max(lls)
    weights = [math.exp(ll - max_ll) for ll in lls]
    Z = sum(weights)
    post = [w / Z for w in weights]

    order = sorted(range(len(thetas)), key=lambda i: post[i], reverse=True)
    sorted_thetas = [thetas[i] for i in order]
    sorted_post = [post[i] for i in order]

    cum = 0.0
    cut = 0
    for cut, p in enumerate(sorted_post):
        cum += p
        if cum >= level:
            break

    cred = sorted_thetas[: cut + 1]

    vals11 = sorted({t[0] for t in cred})
    vals10 = sorted({t[1] for t in cred})
    vals01 = sorted({t[2] for t in cred})
    vals00 = sorted({t[3] for t in cred})

    int11 = _disjoint_intervals(vals11)
    int10 = _disjoint_intervals(vals10)
    int01 = _disjoint_intervals(vals01)
    int00 = _disjoint_intervals(vals00)

    # for this Frechet-based SCS, the denominator is always n
    denoms11 = [n] * len(int11)
    denoms10 = [n] * len(int10)
    denoms01 = [n] * len(int01)
    denoms00 = [n] * len(int00)

    return {
        "theta11_intervals": int11,
        "theta10_intervals": int10,
        "theta01_intervals": int01,
        "theta00_intervals": int00,
        "theta11_denoms": denoms11,
        "theta10_denoms": denoms10,
        "theta01_denoms": denoms01,
        "theta00_denoms": denoms00,
    }


# =====================================================================
# fast MLE
# =====================================================================
def _atnt_tuple(n: int, m: int, xI1: int, xC1: int) -> Theta:
    """All always-takers and never-takers only."""
    return (xI1 + xC1, 0, 0, n - xI1 - xC1)


def _cd_tuple(n: int, m: int, xI1: int, xC1: int) -> Theta:
    """All compliers and defiers only."""
    c = n - m
    return (0, xI1 + c - xC1, m - xI1 + xC1, 0)


def _delta_from_n(n: int) -> int:
    """Base cube half-width for the local likelihood search."""
    return int(math.ceil(max(8.0, 0.21 * (n ** 0.58))))


def _delta_expand_from_n(n: int) -> int:
    """Slightly bigger cube for a one-shot expansion if we hit an edge."""
    base = _delta_from_n(n)
    return int(math.ceil(max(base, 0.30 * (n ** 0.58))))


def _cube_mle(
    n: int,
    m: int,
    xI1: int,
    xC1: int,
    center: Theta,
    delta: int,
    tol: float,
) -> Tuple[List[Theta], float]:
    """
    Search a 4D cube centered at "center" with side length 2*delta+1,
    keeping only those that sum to n and have nonnegative counts.
    """
    t11c, t10c, t01c, t00c = center

    def rng(c0: int):
        lo = max(0, c0 - delta)
        hi = min(n, c0 + delta)
        return range(lo, hi + 1)

    best_ll = -math.inf
    best_set: set[Theta] = set()
    for t11 in rng(t11c):
        for t10 in rng(t10c):
            for t01 in rng(t01c):
                t00 = n - t11 - t10 - t01
                if t00 < 0:
                    continue
                ll = _loglik(n, m, xI1, xC1, t11, t10, t01, t00)
                if ll > best_ll + tol:
                    best_ll = ll
                    best_set = {(t11, t10, t01, t00)}
                elif abs(ll - best_ll) <= tol:
                    best_set.add((t11, t10, t01, t00))
    return sorted(best_set), best_ll


def fast_mle(
    n: int,
    m: int,
    xI1: int,
    xC1: int,
    delta: Optional[int] = None,
    tol: float = 1e-13,
    allow_one_shot_expand: bool = True,
) -> Tuple[List[Theta], float, Dict[str, Any]]:
    """
    Start from data-implied Frechet corners and 2 simple joints
    (AT-NT and C-D). If best is Frechet, complete a grid search around it. If a search hits the
    edge, we may do one expanded pass.
    """
    if not (0 < m < n) or not (0 <= xI1 <= m) or not (0 <= xC1 <= (n - m)):
        raise ValueError("inputs out of range")

    delta_base = _delta_from_n(n) if delta is None else int(delta)

    conj_f, opp_f = _corners_float(n, m, xI1, xC1)
    cands1 = _integerizations_sum_n(n, conj_f)
    cands2 = _integerizations_sum_n(n, opp_f)

    theta_atnt = _atnt_tuple(n, m, xI1, xC1)
    theta_cd = _cd_tuple(n, m, xI1, xC1)

    # Score all initialization candidates
    initial: List[Tuple[Theta, float]] = []
    for t in cands1:
        initial.append((t, _loglik(n, m, xI1, xC1, *t)))
    for t in cands2:
        initial.append((t, _loglik(n, m, xI1, xC1, *t)))
    ll_atnt = _loglik(n, m, xI1, xC1, *theta_atnt)
    ll_cd   = _loglik(n, m, xI1, xC1, *theta_cd)
    initial.append((theta_atnt, ll_atnt))
    initial.append((theta_cd, ll_cd))

    init_best = max(ll for _, ll in initial)

    # ---- NEW: Early-exit if a two-type candidate is among the best seeds ----
    two_type_mles: List[Theta] = []
    if abs(ll_atnt - init_best) <= tol:
        two_type_mles.append(theta_atnt)
    if abs(ll_cd - init_best) <= tol:
        two_type_mles.append(theta_cd)

    if two_type_mles:
        # Return immediately without any local cube search
        return (
            sorted(set(two_type_mles)),
            init_best,
            {
                "delta_used": 0,                # no local search performed
                "hit_edge": False,
                "expanded": False,
                "skipped_local_search": True,
                "reason": "two-type seed attains initial maximum",
            },
        )
    # ---- END NEW ----------------------------------------------------------------

    # Otherwise proceed with the usual local cube search
    seeds = sorted({t for (t, ll) in initial if ll >= init_best - tol})

    global_best_ll = -math.inf
    global_set: set[Theta] = set()
    hit_edge = False

    for center in seeds:
        mles_loc, ll_loc = _cube_mle(n, m, xI1, xC1, center, delta_base, tol)
        if ll_loc > global_best_ll + tol:
            global_best_ll = ll_loc
            global_set = set(mles_loc)
        elif abs(ll_loc - global_best_ll) <= tol:
            global_set.update(mles_loc)

        t11c, t10c, t01c, t00c = center
        for t11, t10, t01, t00 in mles_loc:
            if (
                abs(t11 - t11c) == delta_base
                or abs(t10 - t10c) == delta_base
                or abs(t01 - t01c) == delta_base
                or abs(t00 - t00c) == delta_base
            ):
                hit_edge = True

    expanded = False
    if allow_one_shot_expand and hit_edge:
        expanded = True
        delta2 = _delta_expand_from_n(n)
        for center in seeds:
            mles_loc, ll_loc = _cube_mle(n, m, xI1, xC1, center, delta2, tol)
            if ll_loc > global_best_ll + tol:
                global_best_ll = ll_loc
                global_set = set(mles_loc)
            elif abs(ll_loc - global_best_ll) <= tol:
                global_set.update(mles_loc)

    return (
        sorted(global_set),
        global_best_ll,
        {
            "delta_used": delta_base,
            "hit_edge": hit_edge,
            "expanded": expanded,
            "skipped_local_search": False,
        },
    )

# =====================================================================
# exhaustive grid
# =====================================================================
def _exhaustive_grid(
    n: int,
    m: int,
    xI1: int,
    xC1: int,
    level: float = 0.95,
    show_progress: bool = True,
) -> Dict[str, Any]:
    """
    Try every possible joint count
    that sums to n and is consistent with the observed 2x2 table.
    This is the exact method and will be slow for large n.
    """
    c = n - m
    xI0_req = m - xI1
    xC0_req = c - xC1

    # upper bound on total theta points to scan
    total_upper = 0
    for t11 in range(0, n + 1):
        rem1 = n - t11
        for t10 in range(0, rem1 + 1):
            rem2 = rem1 - t10
            total_upper += (rem2 + 1)

    pbar = tqdm(
        total=total_upper,
        desc="Enumerating Joint Distributions",
        unit="Joint Distribution",
    ) if show_progress else None

    thetas: List[Theta] = []
    lls: List[float] = []

    for t11 in range(0, n + 1):
        rem1 = n - t11
        for t10 in range(0, rem1 + 1):
            rem2 = rem1 - t10
            for t01 in range(0, rem2 + 1):
                t00 = rem2 - t01
                if pbar is not None:
                    pbar.update(1)

                # feasibility checks
                if t01 + t00 < xI0_req:
                    continue
                if t10 + t00 < xC0_req:
                    continue

                ll = _loglik(n, m, xI1, xC1, t11, t10, t01, t00)
                if ll == -math.inf:
                    continue
                thetas.append((t11, t10, t01, t00))
                lls.append(ll)

    if pbar is not None:
        pbar.close()

    enumerated = len(thetas)
    if enumerated == 0:
        return {
            "enumerated": 0,
            "mles": [],
            "max_loglik": float("-inf"),
            "credible_thetas": [],
            "credible_posteriors": [],
            "coverage": 0.0,
            "intervals": {
                "theta11": [],
                "theta10": [],
                "theta01": [],
                "theta00": [],
            },
            "union_str": {
                "theta11": "EMPTY",
                "theta10": "EMPTY",
                "theta01": "EMPTY",
                "theta00": "EMPTY",
            },
        }

    max_ll = max(lls)
    tol = 1e-13
    mles = [thetas[i] for i, ll in enumerate(lls) if abs(ll - max_ll) <= tol]

    # posterior weights over the feasible theta's
    weights = [math.exp(ll - max_ll) for ll in lls]
    Z = sum(weights)
    post = [w / Z for w in weights]

    # sort by posterior mass
    order = sorted(range(enumerated), key=lambda i: post[i], reverse=True)
    sorted_thetas = [thetas[i] for i in order]
    sorted_post = [post[i] for i in order]

    # take top draws until we reach the desired level
    cum = 0.0
    cutoff_idx = 0
    for cutoff_idx, p in enumerate(sorted_post):
        cum += p
        if cum >= level:
            break

    credible_thetas = sorted_thetas[: cutoff_idx + 1]
    credible_post = sorted_post[: cutoff_idx + 1]

    # marginal intervals
    vals11 = sorted({t[0] for t in credible_thetas})
    vals10 = sorted({t[1] for t in credible_thetas})
    vals01 = sorted({t[2] for t in credible_thetas})
    vals00 = sorted({t[3] for t in credible_thetas})

    ints = {
        "theta11": _disjoint_intervals(vals11),
        "theta10": _disjoint_intervals(vals10),
        "theta01": _disjoint_intervals(vals01),
        "theta00": _disjoint_intervals(vals00),
    }

    return {
        "enumerated": enumerated,
        "mles": sorted(set(mles)),
        "max_loglik": max_ll,
        "credible_thetas": credible_thetas,
        "credible_posteriors": credible_post,
        "coverage": cum,
        "intervals": ints,
        "union_str": {
            "theta11": _format_union(ints["theta11"]),
            "theta10": _format_union(ints["theta10"]),
            "theta01": _format_union(ints["theta01"]),
            "theta00": _format_union(ints["theta00"]),
        },
    }


# =====================================================================
# largest possible support
# =====================================================================
def _largest_possible_support(
    n: int, m: int, xI1: int, xC1: int
) -> Dict[str, Tuple[int, int]]:
    """
    For each type, return the broadest possible interval of counts that
    could ever be compatible with this observed 2x2 table.
    """
    c = n - m
    return {
        "theta11": (0, xI1 + xC1),
        "theta10": (0, xI1 + (c - xC1)),
        "theta01": (0, (m - xI1) + xC1),
        "theta00": (0, (m - xI1) + (c - xC1)),
    }


# =====================================================================
# estimated Frechet bounds (discrete)
# =====================================================================
def _estimated_frechet_bounds_discrete(
    n: int, m: int, xI1: int, xC1: int
) -> Dict[str, Any]:
    """
    Integer version of Frechet bounds: round the data-implied margins,
    then walk along that line and collect all feasible integer joints.
    """
    c = n - m
    A = (n / m) * xI1
    B = (n / c) * xC1
    R_est = int(round(A))
    C_est = int(round(B))

    t11_lo = max(0, R_est + C_est - n)
    t11_hi = min(R_est, C_est)
    thetas: List[Theta] = []
    for t11 in range(t11_lo, t11_hi + 1):
        t10 = R_est - t11
        t01 = C_est - t11
        t00 = n - R_est - C_est + t11
        if min(t10, t01, t00) < 0:
            continue
        thetas.append((t11, t10, t01, t00))

    # if rounding made it infeasible, fall back to the continuous-based corner
    if not thetas:
        theta_star = _best_frechet_corner_point(n, m, xI1, xC1)
        t11s, t10s, t01s, t00s = theta_star
        R_est = t11s + t10s
        C_est = t11s + t01s
        t11_lo = max(0, R_est + C_est - n)
        t11_hi = min(R_est, C_est)
        thetas = []
        for t11 in range(t11_lo, t11_hi + 1):
            t10 = R_est - t11
            t01 = C_est - t11
            t00 = n - R_est - C_est + t11
            if min(t10, t01, t00) < 0:
                continue
            thetas.append((t11, t10, t01, t00))

    vals11 = sorted({t[0] for t in thetas})
    vals10 = sorted({t[1] for t in thetas})
    vals01 = sorted({t[2] for t in thetas})
    vals00 = sorted({t[3] for t in thetas})

    return {
        "theta11_intervals": _disjoint_intervals(vals11),
        "theta10_intervals": _disjoint_intervals(vals10),
        "theta01_intervals": _disjoint_intervals(vals01),
        "theta00_intervals": _disjoint_intervals(vals00),
    }


# =====================================================================
# Fisher's exact 2x2
# =====================================================================
def _fisher_exact_2x2(xI1: int, xI0: int, xC1: int, xC0: int) -> float:
    """
    Right-tail Fisher exact test for a 2x2 with margins fixed.
    """
    a = xI1
    b = xI0
    c = xC1
    d = xC0
    m1 = a + b
    m2 = c + d
    n1 = a + c
    N = m1 + m2

    def hypergeom_prob(a_):
        return math.exp(
            _logC_scalar(m1, a_)
            + _logC_scalar(m2, n1 - a_)
            - _logC_scalar(N, n1)
        )

    p_obs = hypergeom_prob(a)

    lo = max(0, n1 - m2)
    hi = min(n1, m1)

    p_val = 0.0
    for aa in range(lo, hi + 1):
        p = hypergeom_prob(aa)
        if p <= p_obs + 1e-13:
            p_val += p
    return min(p_val, 1.0)


# =====================================================================
# formatting helpers
# =====================================================================
def _make_standard_stats_table(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    n: int,
    m: int,
    c: int,
) -> str:
    """
    Construct the basic table of sample size, takeup in each arm,
    difference in takeup, CI, and Fisher p-value.
    """
    # average effect = difference in takeup rates across arms
    diff_pct = (xI1 / m - xC1 / c) * 100

    # CI (difference in takeup rates)
    pI = xI1 / m
    pC = xC1 / c
    diff = pI - pC
    se = math.sqrt(pI * (1 - pI) / m + pC * (1 - pC) / c) if m > 0 and c > 0 else 0.0
    z = 1.96
    ci_lo = diff - z * se
    ci_hi = diff + z * se

    p_fisher = _fisher_exact_2x2(xI1, xI0, xC1, xC0)

    lines = []
    lines.append(" ")
    lines.append("Standard Statistics")
    lines.append("----------------------------------------------")
    lines.append(
        "Average Effect              "
        f"{xI1}/{m} - {xC1}/{c} = {diff_pct:.2f}%"
    )
    lines.append(
        "95% Confidence Interval     "
        f"[{ci_lo*100:.2f}%, {ci_hi*100:.2f}%]"
    )
    lines.append(
        "Fisher's Exact Test p-value " f"{p_fisher:.4g}"
    )
    lines.append(
        "Intervention Takeup Rate    "
        f"{xI1}/{m} = {(xI1/m)*100:.2f}%"
    )
    lines.append(
        "Control Takeup Rate         "
        f"{xC1}/{c} = {(xC1/c)*100:.2f}%"
    )
    lines.append(f"Sample Size                  {n}")
    return "\n".join(lines)


# =====================================================================
# table maker
# =====================================================================
def _make_mle_table(
    n: int,
    mle_list: List[Theta],
    auxiliary: bool,
    method: str,
    largest_support: Dict[str, Tuple[int, int]],
    global_scs: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    est_frechet: Optional[Dict[str, Any]] = None,
    frechet_scs: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Turn the MLE(s) and all optional supporting sets into a readable
    text block.
    """
    if mle_list:
        t11, t10, t01, t00 = mle_list[0]
    else:
        t11 = t10 = t01 = t00 = 0

    if not auxiliary:
        title = "Christy and Kowalski Design-Based Maximum Likelihood Estimates"
    else:
        title = "Christy and Kowalski Design-Based Maximum Likelihood Estimates and Auxiliary Statistics"

    def _get_global_scs_for(key: str) -> Optional[List[Tuple[int, int]]]:
        if global_scs is None:
            return None
        return global_scs.get(key)

    def _get_est_frechet_for(key: str) -> Optional[List[Tuple[int, int]]]:
        if est_frechet is None:
            return None
        return est_frechet.get(key + "_intervals")

    def _get_frechet_scs_for(key: str) -> Optional[List[Tuple[int, int]]]:
        if frechet_scs is None:
            return None
        return frechet_scs.get(key + "_intervals")

    def _get_frechet_denoms_for(key: str) -> Optional[List[int]]:
        if frechet_scs is None:
            return None
        return frechet_scs.get(key + "_denoms")

    def _fmt_single_interval_frac(a: int, b: int, denom: int) -> str:
        if a == b:
            pct = a / denom * 100
            return f"{{{a}}}/{denom} = {pct:.2f}%"
        else:
            pct_a = a / denom * 100
            pct_b = b / denom * 100
            return f"[{a},{b}]/{denom} = [{pct_a:.2f}%, {pct_b:.2f}%]"

    def _fmt_union_same_denom(intervals: List[Tuple[int, int]], denom: int) -> str:
        if not intervals:
            return "EMPTY"
        return " U ".join(_fmt_single_interval_frac(a, b, denom) for (a, b) in intervals)

    def _fmt_union_with_denoms(intervals: List[Tuple[int, int]], denoms: List[int]) -> str:
        if not intervals:
            return "EMPTY"
        left_parts = []
        right_parts = []
        for (a, b), d in zip(intervals, denoms):
            if a == b:
                p = a / d * 100
                left_parts.append(f"{{{a}}}/{d}")
                right_parts.append(f"{p:.2f}%")
            else:
                pa = a / d * 100
                pb = b / d * 100
                left_parts.append(f"[{a},{b}]/{d}")
                right_parts.append(f"[{pa:.2f}%, {pb:.2f}%]")
        left_side = " U ".join(left_parts)
        right_side = " U ".join(right_parts)
        return f"{left_side} = {right_side}"

    def _type_block(
        label: str,
        mle_count: int,
        largest_iv: Tuple[int, int],
        g_scs: Optional[List[Tuple[int, int]]],
        est_fr: Optional[List[Tuple[int, int]]],
        fr_scs: Optional[List[Tuple[int, int]]],
        fr_denoms: Optional[List[int]],
    ) -> List[str]:
        blk: List[str] = []
        blk.append(label)
        blk.append(f"  MLE: {mle_count}/{n} = {(mle_count/n)*100:.2f}%")

        # simplest thing and fast method → stop
        if not auxiliary and method == "approx":
            return blk

        # exhaustive but no auxiliaries → show at least global SCS
        if not auxiliary and method == "exhaustive":
            if g_scs is not None:
                blk.append("  95% Smallest Credible Set: " + _fmt_union_same_denom(g_scs, n))
            return blk

        # auxiliary=True → show everything
        if g_scs is not None:
            blk.append("  95% Smallest Credible Set: " + _fmt_union_same_denom(g_scs, n))

        blk.append("  Largest Possible Support: " + _fmt_union_same_denom([largest_iv], n))

        if est_fr is not None:
            blk.append("  Estimated Frechet Bounds: " + _fmt_union_same_denom(est_fr, n))

        if fr_scs is not None:
            if fr_denoms is not None and len(fr_denoms) == len(fr_scs):
                blk.append(
                    "  95% SCS within Est. Frechet: "
                    + _fmt_union_with_denoms(fr_scs, fr_denoms)
                )
            else:
                blk.append(
                    "  95% SCS within Est. Frechet: "
                    + _fmt_union_same_denom(fr_scs, n)
                )
        return blk

    lines: List[str] = []
    lines.append("")
    lines.append(title)
    lines.append("------------------------------------------------------------------------------------")

    lines += _type_block(
        "Always takers",
        t11,
        largest_support["theta11"],
        _get_global_scs_for("theta11"),
        _get_est_frechet_for("theta11"),
        _get_frechet_scs_for("theta11"),
        _get_frechet_denoms_for("theta11"),
    )
    lines.append("")

    lines += _type_block(
        "Compliers",
        t10,
        largest_support["theta10"],
        _get_global_scs_for("theta10"),
        _get_est_frechet_for("theta10"),
        _get_frechet_scs_for("theta10"),
        _get_frechet_denoms_for("theta10"),
    )
    lines.append("")

    lines += _type_block(
        "Defiers",
        t01,
        largest_support["theta01"],
        _get_global_scs_for("theta01"),
        _get_est_frechet_for("theta01"),
        _get_frechet_scs_for("theta01"),
        _get_frechet_denoms_for("theta01"),
    )
    lines.append("")

    lines += _type_block(
        "Never takers",
        t00,
        largest_support["theta00"],
        _get_global_scs_for("theta00"),
        _get_est_frechet_for("theta00"),
        _get_frechet_scs_for("theta00"),
        _get_frechet_denoms_for("theta00"),
    )

    return "\n".join(lines)


# =====================================================================
# result object
# =====================================================================
class DBMLEResult(dict):
    """
    Wrapper around a dict that carries both the structured output
    and a printable report.
    """
    def report(self) -> str:
        # keep the key name for compatibility
        return self.get("report", repr(self))

# ── helpers: parse a single value to binary 0/1 if possible ─────────────────────
def _try_parse_binary(v) -> Optional[int]:
    """
    Return 0 or 1 if v represents exactly binary {0,1}; otherwise return None.
    Accepts: bool, int in {0,1}, float in {0.0,1.0}, str in {"0","1"} (with whitespace).
    """
    if isinstance(v, bool):
        return int(v)
    if isinstance(v, int):
        return v if v in (0, 1) else None
    if isinstance(v, float):
        if not math.isfinite(v):
            return None
        return int(v) if v in (0.0, 1.0) else None
    if isinstance(v, str):
        s = v.strip()
        if s in ("0", "1"):
            return int(s)
        return None
    return None  # unsupported type


def _sanitize_ZD(
    Z,
    D,
    invalid_policy: str = "drop",
    warn_on_invalid: bool = True,
    warn_limit: int = 10,
) -> Tuple[int, int, int, int, Dict[str, Any]]:
    """
    Validate+aggregate (Z, D) to 2x2 counts with flexible invalid handling.

    Returns
    -------
    xI1, xI0, xC1, xC0, stats

    stats includes:
        total, used, dropped, coerced_Z, coerced_D,
        invalid_indices (up to warn_limit),
        policy
    """
    if Z is None or D is None:
        raise ValueError("Z and D must be sequences of equal length.")
    if len(Z) != len(D):
        raise ValueError(f"Z and D must have the same length (got {len(Z)} vs {len(D)}).")
    if invalid_policy not in {"drop", "coerce-0", "coerce-1", "raise"}:
        raise ValueError(
            "invalid_policy must be one of {'drop','coerce-0','coerce-1','raise'}."
        )

    default_val = None
    if invalid_policy == "coerce-0":
        default_val = 0
    elif invalid_policy == "coerce-1":
        default_val = 1

    xI1 = xI0 = xC1 = xC0 = 0
    total = len(Z)
    used = 0
    dropped = 0
    coerced_Z = 0
    coerced_D = 0
    bad_idxs: List[int] = []

    for i, (z_raw, d_raw) in enumerate(zip(Z, D)):
        z = _try_parse_binary(z_raw)
        d = _try_parse_binary(d_raw)

        if (z is None) or (d is None):
            if invalid_policy == "raise":
                which = []
                if z is None: which.append(f"Z[{i}]={z_raw!r}")
                if d is None: which.append(f"D[{i}]={d_raw!r}")
                raise ValueError(f"Invalid {', '.join(which)}; expected 0/1.")
            elif invalid_policy == "drop":
                dropped += 1
                if len(bad_idxs) < warn_limit:
                    bad_idxs.append(i)
                continue
            else:
                # coerce to default_val (0 or 1)
                if z is None:
                    z = default_val  # type: ignore[assignment]
                    coerced_Z += 1
                if d is None:
                    d = default_val  # type: ignore[assignment]
                    coerced_D += 1

        # aggregate
        if z == 1:
            if d == 1: xI1 += 1
            else:      xI0 += 1
        else:
            if d == 1: xC1 += 1
            else:      xC0 += 1
        used += 1

    # Warn once with a concise summary
    if warn_on_invalid and (dropped > 0 or coerced_Z > 0 or coerced_D > 0):
        details = []
        if dropped:
            details.append(f"dropped={dropped}")
        if coerced_Z:
            details.append(f"coerced_Z={coerced_Z}")
        if coerced_D:
            details.append(f"coerced_D={coerced_D}")
        where = f"; first invalid indices: {bad_idxs}" if bad_idxs else ""
        warnings.warn(
            f"dbmle_from_ZD: invalid entries handled by policy='{invalid_policy}' "
            f"({', '.join(details)}; total={total}, used={used}){where}",
            UserWarning,
            stacklevel=2,
        )

    # sanity: both arms must be non-empty after handling
    I = xI1 + xI0
    C = xC1 + xC0
    if I == 0 or C == 0:
        raise ValueError(
            f"After applying invalid_policy='{invalid_policy}', one arm is empty: I={I}, C={C}. "
            "Ensure Z (after cleaning) contains at least one 1 and one 0."
        )

    stats = {
        "total": total,
        "used": used,
        "dropped": dropped,
        "coerced_Z": coerced_Z,
        "coerced_D": coerced_D,
        "invalid_indices": bad_idxs,
        "policy": invalid_policy,
    }
    return xI1, xI0, xC1, xC0, stats


# =====================================================================
# main command
# =====================================================================
def dbmle(
    xI1: int,
    xI0: int,
    xC1: int,
    xC0: int,
    method: str = "approx",
    auxiliary: bool = False,
    level: float = 0.95,
    show_progress: bool = True,
) -> DBMLEResult:
    """
    Pick method (approx vs exhaustive) and depth of
    statistics (auxiliary = False vs True), and return a structured result
    plus a printable report string.
    """
    _validate_inputs(xI1, xI0, xC1, xC0)
    m = xI1 + xI0
    c = xC1 + xC0
    n = m + c

    out: Dict[str, Any] = {
        "summary": {
            "inputs": {
                "xI1": xI1,
                "xI0": xI0,
                "xC1": xC1,
                "xC0": xC0,
                "n": n,
                "m": m,
                "c": c,
            }
        },
        "meta": {
            "method": method,
            "auxiliary": auxiliary,
            "level": level,
        },
    }

    # Warning when approx + auxiliary=True
    if method.lower() == "approx" and auxiliary:
        warnings.warn(
            "Note: When 'auxiliary=True' is specified with method='approx', "
            "the function automatically performs an exhaustive grid search to "
            "compute the true MLE and 95% credible set. The MLE is exact as well.",
            UserWarning,
            stacklevel=2,
        )

    largest_support = _largest_possible_support(n, m, xI1, xC1)
    est_frechet = _estimated_frechet_bounds_discrete(n, m, xI1, xC1)

    # ================================================================
    # 1) EXHAUSTIVE
    # ================================================================
    if method.lower() == "exhaustive":
        grid = _exhaustive_grid(n, m, xI1, xC1, level=level, show_progress=show_progress)
        mles = grid["mles"]
        global_ints = grid["intervals"]

        if auxiliary:
            fre_scs = _frechet_marginal_scs(n, m, xI1, xC1, level=level)
        else:
            fre_scs = None

        std_tbl = _make_standard_stats_table(xI1, xI0, xC1, xC0, n, m, c)

        # If there are tied MLEs, print one block per MLE (even when auxiliary=False)
        if len(mles) > 1:
            blocks: List[str] = []
            for idx, mle_theta in enumerate(mles, start=1):
                this_tbl = _make_mle_table(
                    n,
                    [mle_theta],   # force THIS MLE to be first/only in the table
                    auxiliary,
                    "exhaustive",
                    largest_support,
                    global_scs=global_ints,                   # show SCS if method=exhaustive
                    est_frechet=est_frechet if auxiliary else None,
                    frechet_scs=fre_scs if auxiliary else None,
                )
                blocks.append(f"(tied MLE #{idx})\n{this_tbl}")
            mle_tbl = "\n\n".join(blocks)
        else:
            mle_tbl = _make_mle_table(
                n,
                mles,
                auxiliary,
                "exhaustive",
                largest_support,
                global_scs=global_ints,
                est_frechet=est_frechet if auxiliary else None,
                frechet_scs=fre_scs if auxiliary else None,
            )

        out.update(
            {
                "mle": {"mle_list": mles, "max_loglik": grid["max_loglik"]},
                "supports": {
                    "largest_possible_support": largest_support,
                    "estimated_frechet_bounds": est_frechet,
                },
                "global_95_scs": {
                    "intervals": global_ints,
                    "union_str": grid["union_str"],
                },
                "report": "\n" + std_tbl + "\n\n" + mle_tbl,
            }
        )
        if fre_scs is not None:
            out["frechet_95_scs"] = fre_scs

        return DBMLEResult(out)

    # ================================================================
    # 2) APPROX + NO AUXILIARY (old "proposed")
    # ================================================================
    if method.lower() == "approx" and not auxiliary:
        mles_fast, ll_fast, meta_fast = fast_mle(
            n, m, xI1, xC1, allow_one_shot_expand=True
        )
        out["mle"] = {"mle_list": mles_fast, "max_loglik": ll_fast}
        out["meta"]["fast_mle"] = meta_fast

        std_tbl = _make_standard_stats_table(xI1, xI0, xC1, xC0, n, m, c)

        # If there are tied MLEs, print one block per MLE even without auxiliaries
        if len(mles_fast) > 1:
            blocks2: List[str] = []
            for idx, mle_theta in enumerate(mles_fast, start=1):
                this_tbl = _make_mle_table(
                    n,
                    [mle_theta],   # single MLE per block
                    False,         # auxiliary=False
                    "approx",
                    largest_support,
                )
                blocks2.append(f"(tied MLE #{idx})\n{this_tbl}")
            mle_tbl = "\n\n".join(blocks2)
        else:
            mle_tbl = _make_mle_table(
                n,
                mles_fast,
                False,
                "approx",
                largest_support,
            )

        out["report"] = "\n" + std_tbl + "\n\n" + mle_tbl
        out["supports"] = {
            "largest_possible_support": largest_support,
            "estimated_frechet_bounds": est_frechet,
        }
        return DBMLEResult(out)


    # ================================================================
    # 3) APPROX + AUXILIARY (same behavior as exhaustive + auxiliary)
    # ================================================================
    grid = _exhaustive_grid(n, m, xI1, xC1, level=level, show_progress=show_progress)
    mles_exact = grid["mles"]
    global_ints = grid["intervals"]
    fre_scs = _frechet_marginal_scs(n, m, xI1, xC1, level=level)

    std_tbl = _make_standard_stats_table(xI1, xI0, xC1, xC0, n, m, c)

    if len(mles_exact) > 1:
        aux_blocks2: List[str] = []
        for idx, mle_theta in enumerate(mles_exact, start=1):
            this_tbl = _make_mle_table(
                n,
                [mle_theta],
                True,
                "approx",         # label: user asked approx
                largest_support,
                global_scs=global_ints,
                est_frechet=est_frechet,
                frechet_scs=fre_scs,
            )
            aux_blocks2.append(f"(tied MLE #{idx})\n{this_tbl}")
        mle_tbl = "\n\n".join(aux_blocks2)
    else:
        mle_tbl = _make_mle_table(
            n,
            mles_exact,
            True,
            "approx",                   # label remains "approx" because user requested it
            largest_support,
            global_scs=global_ints,
            est_frechet=est_frechet,
            frechet_scs=fre_scs,
        )

    out.update(
        {
            "mle": {"mle_list": mles_exact, "max_loglik": grid["max_loglik"]},
            "supports": {
                "largest_possible_support": largest_support,
                "estimated_frechet_bounds": est_frechet,
            },
            "global_95_scs": {
                "intervals": global_ints,
                "union_str": grid["union_str"],
            },
            "frechet_95_scs": fre_scs,
            "report": "\n" + std_tbl + "\n\n" + mle_tbl,
        }
    )
    return DBMLEResult(out)


def dbmle_from_ZD(
    Z: List[int],
    D: List[int],
    method: str = "approx",
    auxiliary: bool = False,
    level: float = 0.95,
    show_progress: bool = True,
    invalid_policy: str = "drop",     # "drop" | "coerce-0" | "coerce-1" | "raise"
    warn_on_invalid: bool = True,
) -> DBMLEResult:
    """
    Takes individual-level data (Z, D) and converts to counts, with robust handling of
    invalid/missing entries.

    Parameters
    ----------
    Z : sequence of 0/1 (bool/int/float/str accepted if exactly {0,1})
        Randomized assignment (1 = intervention, 0 = control).
    D : sequence of 0/1 (bool/int/float/str accepted if exactly {0,1})
        Observed take-up (1 = took up treatment, 0 = did not).
    method : {"approx","exhaustive"}, default "approx"
    auxiliary : bool, default False
    level : float, default 0.95
    show_progress : bool, default True
    invalid_policy : {"drop","coerce-0","coerce-1","raise"}, default "drop"
        - "drop": ignore any record with invalid Z or D; warn with counts.
        - "coerce-0": coerce invalid Z/D values to 0; warn with counts.
        - "coerce-1": coerce invalid Z/D values to 1; warn with counts.
        - "raise": strict mode; raise ValueError on the first invalid entry.
    warn_on_invalid : bool, default True
        Emit a single summary warning if any invalid entries were dropped/coerced.

    Raises
    ------
    ValueError
        - If Z and D lengths differ.
        - If invalid_policy is unknown.
        - If, after cleaning, either arm is empty (no Z==1 or no Z==0).

    Returns
    -------
    DBMLEResult
    """
    xI1, xI0, xC1, xC0, stats = _sanitize_ZD(
        Z, D, invalid_policy=invalid_policy, warn_on_invalid=warn_on_invalid
    )

    result = dbmle(
        xI1,
        xI0,
        xC1,
        xC0,
        method=method,
        auxiliary=auxiliary,
        level=level,
        show_progress=show_progress,
    )

    # attach cleaning stats for transparency/debugging
    meta = result.setdefault("meta", {})
    meta["from_ZD"] = stats
    return result
