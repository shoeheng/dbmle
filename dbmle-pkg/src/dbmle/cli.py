# src/dbmle/cli.py

import argparse

# Import the public API
from dbmle.core import dbmle


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Design-based MLE for always-takers, compliers, defiers, and never-takers "
            "from 2x2 counts (xI1, xI0, xC1, xC0). "
            "For Z/D individual-level input, call dbmle_from_ZD in Python."
        )
    )

    # Required 2x2 counts
    parser.add_argument("--xI1", type=int, required=True, help="Intervention: took up treatment")
    parser.add_argument("--xI0", type=int, required=True, help="Intervention: did not take up")
    parser.add_argument("--xC1", type=int, required=True, help="Control: took up treatment")
    parser.add_argument("--xC0", type=int, required=True, help="Control: did not take up")

    # Method (fast vs exact)
    parser.add_argument(
        "--method",
        choices=["approx", "exhaustive"],
        default="approx",
        help="MLE method: fast ('approx') or exact grid search ('exhaustive').",
    )

    # Auxiliary stats toggle (flag only; do NOT pass True/False after it)
    parser.add_argument(
        "--auxiliary",
        action="store_true",
        help="If set, include auxiliary statistics (credible sets, bounds).",
    )

    # Credible set level
    parser.add_argument(
        "--level",
        type=float,
        default=0.95,
        help="Credible set confidence level (default: 0.95).",
    )

    # Progress bar toggle
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Hide progress bar (relevant only for exhaustive mode).",
    )

    args = parser.parse_args()

    # Run the estimator
    res = dbmle(
        args.xI1,
        args.xI0,
        args.xC1,
        args.xC0,
        method=args.method,
        auxiliary=args.auxiliary,
        level=args.level,
        show_progress=not args.no_progress,
    )

    # Print the human-readable report
    print(res.report())


if __name__ == "__main__":
    main()
