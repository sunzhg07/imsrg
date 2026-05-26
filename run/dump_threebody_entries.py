#!/usr/bin/env /usr/bin/python3
"""Dump three-body diagram entries for the first N g1==g2 cases."""

import argparse
import contextlib
import io
from pathlib import Path

import test_diagram_comm as tdc


def main():
    parser = argparse.ArgumentParser(
        description="Write first N g1==g2 three-body entry cases to a CSV-like file."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output file.",
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=100,
        help="Number of g1==g2 cases to dump.",
    )
    args = parser.parse_args()

    with contextlib.redirect_stdout(io.StringIO()):
        ms2, hs, eom, configs, configs2 = tdc.setup_problem()

    rows = [
        "pair_i,pair_j,a,b,c,d,e,f,g,jab_seed,jde_seed,entry_jab,entry_jde,twoJ,val"
    ]

    checked = 0
    for pair_i, pair_j, twoJ in tdc.iter_pair_jobs_from_orbits(
        ms2,
        configs2,
        list(range(len(configs2))),
        list(range(len(configs2))),
    ):
        _, _, d, g1, a, b, j0 = configs2[pair_i]
        _, _, c, g2, f, e, j2 = configs2[pair_j]
        if g1 != g2:
            continue

        checked += 1
        entries = list(eom.ThreeBody_Diagram_Entries(a, b, c, d, e, f, g1, int(j0), int(j2)))
        if not entries:
            rows.append(
                f"{pair_i},{pair_j},{a},{b},{c},{d},{e},{f},{g1},{int(j0)},{int(j2)},-1,-1,{twoJ},0.0"
            )
        else:
            for entry_jab, entry_jde, entry_twoJ, value in entries:
                rows.append(
                    f"{pair_i},{pair_j},{a},{b},{c},{d},{e},{f},{g1},{int(j0)},{int(j2)},"
                    f"{int(entry_jab)},{int(entry_jde)},{int(entry_twoJ)},{float(value):.15g}"
                )

        if checked >= args.max_cases:
            break

    output_path = Path(args.output)
    output_path.write_text("\n".join(rows) + "\n")
    print(f"wrote {output_path}")
    print(f"cases={checked}")


if __name__ == "__main__":
    main()