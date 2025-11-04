# src/counting_defiers/cli.py

import argparse
from .core import counting_defiers_command


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Design-based MLE for always-takers, compliers, defiers, and never-takers "
            "from randomized encouragement (Z,D) data."
        )
    )

    # Define arguments
    parser.add_argument("--xI1", type=int, required=True, help="Intervention: took up treatment")
    parser.add_argument("--xI0", type=int, required=True, help="Intervention: did not take up")
    parser.add_argument("--xC1", type=int, required=True, help="Control: took up treatment")
    parser.add_argument("--xC0", type=int, required=True, help="Control: did not take up")

    parser.add_argument(
        "--method",
        choices=["approx", "exhaustive"],
        default="approx",
        help="MLE method: fast ('approx') or exact grid search ('exhaustive').",
    )

    parser.add_argument(
        "--which-stats",
        choices=["proposed", "auxiliary"],
        default="proposed",
        help="Level of detail: 'proposed' (MLE only) or 'auxiliary' (full credible sets).",
    )

    parser.add_argument(
        "--level",
        type=float,
        default=0.95,
        help="Credible set confidence level (default: 0.95).",
    )

    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Hide progress bar (relevant only for exhaustive mode).",
    )

    args = parser.parse_args()

    # Run the estimator
    res = counting_defiers_command(
        args.xI1,
        args.xI0,
        args.xC1,
        args.xC0,
        method=args.method,
        which_stats=args.which_stats,
        level=args.level,
        show_progress=not args.no_progress,
    )

    # Print the report 
    print(res.report())


if __name__ == "__main__":
    main()
