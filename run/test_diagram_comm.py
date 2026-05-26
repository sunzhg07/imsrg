#!/usr/bin/env /usr/bin/python3
"""Compare threebody_diagram_comm against comm223ss for one pair or many pairs."""

import argparse
import sys
import numpy as np

sys.path.insert(0, "/Users/wolf/work/imsrg/build")
from pyIMSRG import *


def load_configs(path, min_len):
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < min_len:
                continue
            rows.append(tuple([parts[0]] + [int(x) for x in parts[1:min_len]]))
    return rows


def setup_problem():
    emax = 3
    ms2 = ModelSpace(emax, "He4", "p-shell")
    ms2.SetHbarOmega(16)

    rank_j, parity, rank_Tz, particle_rank = 0, 0, 0, 2
    hs = Operator(ms2, rank_j, parity, rank_Tz, particle_rank)

    eom = EOM(hs, "he8.ref", 0, 0, 0)
    eom.ConstructConfigs()
    eom.ConstructNormMatrix()

    configs = load_configs("cfs", 6)
    configs2 = load_configs("cfs2", 7)

    return ms2, hs, eom, configs, configs2


def allowed_twoJ_values(ms2, p, q, r, pair_J):
    op = ms2.GetOrbit(p)
    oq = ms2.GetOrbit(q)
    or_ = ms2.GetOrbit(r)

    pair_twoJ = 2 * int(pair_J)
    pair_min = abs(op.j2 - oq.j2)
    pair_max = op.j2 + oq.j2
    if pair_twoJ < pair_min or pair_twoJ > pair_max:
        return []

    total_min = abs(pair_twoJ - or_.j2)
    total_max = pair_twoJ + or_.j2

    parity = (op.l + oq.l + or_.l) % 2
    two_tz = op.tz2 + oq.tz2 + or_.tz2

    twoJ_values = []
    for twoJ in range(total_min, total_max + 1):
        channel_index = ms2.GetThreeBodyChannelIndex(twoJ, parity, two_tz)
        if channel_index < 0:
            continue
        twoJ_values.append(twoJ)

    return twoJ_values


def build_seed_operators(hs, configs, configs2, pair_i, pair_j):
    cfsi = configs[pair_i]
    cfsj = configs[pair_j]
    _, _, d, g1, a, b, j0 = configs2[pair_i]
    _, _, c, g2, f, e, j2 = configs2[pair_j]
    _, _, c0_i, c1_i, c2_i, _ = cfsi
    _, _, c0_j, c1_j, c2_j, _ = cfsj

    if g1 != g2:
        raise RuntimeError(f"Expected g1 == g2, got g1={g1}, g2={g2}")

    x_test = 0.0 * hs
    x_test.SetHermitian()
    y_test = 0.0 * hs
    y_test.SetAntiHermitian()

    x_test.TwoBody.SetTBME_chij(c2_i, c2_i, c0_i, c1_i, 1.0)
    y_test.TwoBody.SetTBME_chij(c2_j, c2_j, c0_j, c1_j, 1.0)

    return {
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "e": e,
        "f": f,
        "g": g1,
        "j0": int(j0),
        "j2": int(j2),
        "x_test": x_test,
        "y_test": y_test,
        "x_seed": (c2_i, c0_i, c1_i),
        "y_seed": (c2_j, c0_j, c1_j),
    }


def compare_pair(ms2, hs, eom, configs, configs2, pair_i, pair_j, twoJ, atol=1e-10):
    seed = build_seed_operators(hs, configs, configs2, pair_i, pair_j)
    rank_j, parity, rank_Tz = 0, 0, 0

    z_test = Operator(ms2, rank_j, rank_Tz, parity, 3)
    z_test.SetHermitian()
    z_test.ThreeBody.SetMode("pn")

    Commutator.comm223ss(seed["x_test"], seed["y_test"], z_test)

    mat, J1min, J2min = eom.threebody_diagram_comm(
        seed["a"], seed["b"], seed["c"],
        seed["d"], seed["e"], seed["f"],
        twoJ,
        seed["a"], seed["b"], seed["d"], seed["g"], seed["j0"],
        seed["c"], seed["f"], seed["e"], seed["j2"],
        1, -1,
    )

    mismatches = []
    nonzero_reference = []

    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            J1 = J1min + row
            J2 = J2min + col
            diag_val = mat[row, col]
            z_direct = z_test.ThreeBody.GetME_pn(
                J1, J2, twoJ, seed["a"], seed["b"], seed["c"], seed["d"], seed["e"], seed["f"]
            )
            z_swapped = z_test.ThreeBody.GetME_pn(
                J2, J1, twoJ, seed["a"], seed["b"], seed["c"], seed["d"], seed["e"], seed["f"]
            )

            if abs(z_direct) > atol:
                nonzero_reference.append((J1, J2, z_direct))

            if abs(diag_val) < atol and abs(z_direct) < atol and abs(z_swapped) < atol:
                continue

            if np.isclose(diag_val, z_direct, atol=atol):
                match = "direct"
            elif np.isclose(diag_val, z_swapped, atol=atol):
                match = "swapped"
            elif abs(diag_val) < atol and abs(z_direct) > atol:
                match = "missing"
            elif abs(diag_val) > atol and abs(z_direct) < atol:
                match = "extra"
            else:
                match = "other"

            if match != "direct":
                mismatches.append({
                    "row": row,
                    "col": col,
                    "J1": J1,
                    "J2": J2,
                    "diagram": diag_val,
                    "z_direct": z_direct,
                    "z_swapped": z_swapped,
                    "match": match,
                })

    return {
        "pair_i": pair_i,
        "pair_j": pair_j,
        "twoJ": twoJ,
        "bra": (seed["a"], seed["b"], seed["c"]),
        "ket": (seed["d"], seed["e"], seed["f"]),
        "g": seed["g"],
        "j0": seed["j0"],
        "j2": seed["j2"],
        "x_seed": seed["x_seed"],
        "y_seed": seed["y_seed"],
        "mat": mat,
        "J1min": J1min,
        "J2min": J2min,
        "mismatches": mismatches,
        "nonzero_reference": nonzero_reference,
    }


def print_pair_details(result):
    pair_i = result["pair_i"]
    pair_j = result["pair_j"]
    twoJ = result["twoJ"]
    a, b, c = result["bra"]
    d, e, f = result["ket"]
    x_seed = result["x_seed"]
    y_seed = result["y_seed"]
    mat = result["mat"]

    print(f"pair_i={pair_i} pair_j={pair_j} twoJ={twoJ}")
    print(f"bra=({a},{b},{c}) ket=({d},{e},{f}) g={result['g']} j0={result['j0']} j2={result['j2']}")
    print(f"X seed: channel={x_seed[0]} local=({x_seed[1]},{x_seed[2]})")
    print(f"Y seed: channel={y_seed[0]} local=({y_seed[1]},{y_seed[2]})")
    print(f"J1min={result['J1min']} J2min={result['J2min']} shape={mat.shape}")
    print()
    print("Diagram matrix")
    print(mat)
    print()

    header = (
        f"{'row':>3} {'col':>3} {'J1':>3} {'J2':>3} "
        f"{'diagram':>14} {'Z(J1,J2)':>14} {'Z(J2,J1)':>14} {'match':>10}"
    )
    print(header)
    print("-" * len(header))

    for item in result["mismatches"]:
        print(
            f"{item['row']:3d} {item['col']:3d} {item['J1']:3d} {item['J2']:3d} "
            f"{item['diagram']:14.8f} {item['z_direct']:14.8f} {item['z_swapped']:14.8f} {item['match']:>10}"
        )

    if not result["mismatches"]:
        for row in range(mat.shape[0]):
            for col in range(mat.shape[1]):
                J1 = result["J1min"] + row
                J2 = result["J2min"] + col
                value = mat[row, col]
                if abs(value) > 1e-10:
                    print(f"{row:3d} {col:3d} {J1:3d} {J2:3d} {value:14.8f} {value:14.8f} {0.0:14.8f} {'direct':>10}")

    print()
    print("Nonzero reference entries Z(J1,J2)")
    for J1, J2, z_val in result["nonzero_reference"]:
        print(f"  Z({J1},{J2}) = {z_val: .8f}")


def iter_pair_jobs(configs2, pair_i_values, pair_j_values):
    for pair_i in pair_i_values:
        for pair_j in pair_j_values:
            twoJ = int(configs2[pair_j][-1])
            yield pair_i, pair_j, twoJ


def iter_pair_jobs_from_orbits(ms2, configs2, pair_i_values, pair_j_values):
    for pair_i in pair_i_values:
        _, _, d, g1, a, b, j0 = configs2[pair_i]
        bra_twoJ_values = set(allowed_twoJ_values(ms2, a, b, d, int(j0)))
        for pair_j in pair_j_values:
            _, _, c, g2, f, e, j2 = configs2[pair_j]
            if g1 != g2:
                continue
            ket_twoJ_values = set(allowed_twoJ_values(ms2, c, f, e, int(j2)))
            twoJ_values = sorted(bra_twoJ_values & ket_twoJ_values)
            for twoJ in twoJ_values:
                yield pair_i, pair_j, twoJ


def build_pair_range(limit, start, stop, size):
    lo = 0 if start is None else start
    hi = size if stop is None else min(stop, size)
    values = list(range(lo, hi))
    if limit is not None:
        values = values[:limit]
    return values


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-i", type=int, default=0)
    parser.add_argument("--pair-j", type=int, default=6)
    parser.add_argument("--twoJ", type=int, default=None)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--limit-i", type=int, default=None)
    parser.add_argument("--limit-j", type=int, default=None)
    parser.add_argument("--start-i", type=int, default=None)
    parser.add_argument("--stop-i", type=int, default=None)
    parser.add_argument("--start-j", type=int, default=None)
    parser.add_argument("--stop-j", type=int, default=None)
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--max-failures", type=int, default=20)
    parser.add_argument("--show-passing", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    ms2, hs, eom, configs, configs2 = setup_problem()

    if not args.sweep:
        twoJ = args.twoJ
        if twoJ is None:
            _, _, d, _, a, b, j0 = configs2[args.pair_i]
            _, _, c, _, f, e, j2 = configs2[args.pair_j]
            bra_twoJ_values = set(allowed_twoJ_values(ms2, a, b, d, int(j0)))
            ket_twoJ_values = set(allowed_twoJ_values(ms2, c, f, e, int(j2)))
            choices = sorted(bra_twoJ_values & ket_twoJ_values)
            if not choices:
                raise RuntimeError(f"No valid twoJ found for pair_i={args.pair_i}, pair_j={args.pair_j}")
            twoJ = choices[0]
        result = compare_pair(ms2, hs, eom, configs, configs2, args.pair_i, args.pair_j, twoJ)
        print_pair_details(result)
        return

    pair_i_values = build_pair_range(args.limit_i, args.start_i, args.stop_i, len(configs))
    pair_j_values = build_pair_range(args.limit_j, args.start_j, args.stop_j, len(configs))

    total = 0
    passing = 0
    failures = []

    for pair_i, pair_j, twoJ in iter_pair_jobs_from_orbits(ms2, configs2, pair_i_values, pair_j_values):
        if args.max_cases is not None and total >= args.max_cases:
            break
        total += 1
        try:
            result = compare_pair(ms2, hs, eom, configs, configs2, pair_i, pair_j, twoJ)
        except Exception as exc:
            failures.append({
                "pair_i": pair_i,
                "pair_j": pair_j,
                "twoJ": twoJ,
                "error": str(exc),
            })
            if len(failures) >= args.max_failures:
                break
            continue

        if result["mismatches"]:
            failures.append(result)
            if len(failures) >= args.max_failures:
                break
        else:
            passing += 1
            if args.show_passing:
                print(
                    f"PASS pair_i={pair_i} pair_j={pair_j} twoJ={twoJ} "
                    f"bra={result['bra']} ket={result['ket']}"
                )

    print(f"checked={total} passing={passing} failing={len(failures)}")

    for item in failures:
        if "error" in item:
            print(
                f"ERROR pair_i={item['pair_i']} pair_j={item['pair_j']} twoJ={item['twoJ']} "
                f"message={item['error']}"
            )
            continue

        print()
        print_pair_details(item)


if __name__ == "__main__":
    main()
