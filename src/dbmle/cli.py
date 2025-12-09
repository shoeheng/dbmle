# src/dbmle/cli.py

import argparse

from dbmle.core import dbmle


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Design-based MLE for always takers, compliers, defiers, and never takers "
            "from aggregate counts (xI1, xI0, xC1, xC0).\n\n"
            "Tip: For individual-level Z/D input, call `dbmle_from_ZD` from Python."
        )
    )

    parser.add_argument("--xI1", type=int, required=True, help="Takeup in Intervention")
    parser.add_argument("--xI0", type=int, required=True, help="No Takeup in Intervention")
    parser.add_argument("--xC1", type=int, required=True, help="Takeup in Control")
    parser.add_argument("--xC0", type=int, required=True, help="No Takeup in Control")

    parser.add_argument(
        "--output",
        choices=["basic", "auxiliary", "approx"],  
        default="basic",
        help=(
            "What to compute and display:\n"
            "  - basic     : exhaustive grid; Standard Stats + MLE(s) + global 95%% SCS\n"
            "  - auxiliary : exhaustive grid; adds more stats from auxiliary table\n"
            "  - approx    : fast approximate MLE only (no auxiliaries)\n"
        ),
    )

    parser.add_argument(
        "--level",
        type=float,
        default=0.95,
        help="Credible-set level (default: 0.95).",
    )

    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Hide progress bar (relevant only for exhaustive modes: 'basic'/'auxiliary').",
    )

    args = parser.parse_args()

    res = dbmle(
        args.xI1,
        args.xI0,
        args.xC1,
        args.xC0,
        output=args.output,
        level=args.level,
        show_progress=not args.no_progress,
    )

    print(res.report())


if __name__ == "__main__":
    main()



