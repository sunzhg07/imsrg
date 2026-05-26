#!/usr/bin/env /usr/bin/python3
"""Compare two dumps produced by dump_threebody_entries.py."""

import argparse
import csv
from collections import defaultdict


def load_rows(path):
    grouped = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            case_key = (
                int(row["pair_i"]),
                int(row["pair_j"]),
                int(row["a"]),
                int(row["b"]),
                int(row["c"]),
                int(row["d"]),
                int(row["e"]),
                int(row["f"]),
                int(row["g"]),
                int(row["jab_seed"]),
                int(row["jde_seed"]),
            )
            entry_key = (
                int(row["entry_jab"]),
                int(row["entry_jde"]),
                int(row["twoJ"]),
            )
            grouped[case_key].append((entry_key, float(row["val"])))
    return grouped


def normalize(entries):
    return sorted(entries, key=lambda item: item[0])


def main():
    parser = argparse.ArgumentParser(description="Compare two three-body entry dump files.")
    parser.add_argument("left")
    parser.add_argument("right")
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--max-mismatches", type=int, default=20)
    args = parser.parse_args()

    left = load_rows(args.left)
    right = load_rows(args.right)

    mismatches = 0
    all_cases = sorted(set(left) | set(right))
    for case in all_cases:
      left_entries = {k: v for k, v in normalize(left.get(case, []))}
      right_entries = {k: v for k, v in normalize(right.get(case, []))}
      if set(left_entries) != set(right_entries):
          mismatches += 1
          print(f"CASE {case} entry-keys differ")
          print(f"  left only: {sorted(set(left_entries) - set(right_entries))}")
          print(f"  right only: {sorted(set(right_entries) - set(left_entries))}")
      else:
          bad = []
          for key in sorted(left_entries):
              if abs(left_entries[key] - right_entries[key]) > args.atol:
                  bad.append((key, left_entries[key], right_entries[key]))
          if bad:
              mismatches += 1
              print(f"CASE {case} value mismatch")
              for key, lval, rval in bad[:8]:
                  print(f"  {key}: left={lval:.15g} right={rval:.15g}")

      if mismatches >= args.max_mismatches:
          break

    print(f"checked_cases={len(all_cases)}")
    print(f"mismatched_cases={mismatches}")


if __name__ == "__main__":
    main()