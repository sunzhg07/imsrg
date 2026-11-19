#!/usr/bin/env python3
"""Convert kshell_threebody transit TDTD output (rank 0, same state) to the
OSM `.ref` TRBTD format read by EOM::ReadOsmRdm.

  python3 tdtd2ref.py <interaction.snt> <transit_log> <out.ref> [il ir [twoJ_state]]

Conversion (kshell_threebody README, "TDTD output format"):
  OSM = TDTD * sqrt(2 J(ab)c + 1) * sqrt(g_abc * g_def),  g = 1, 2, 6
  for 0, 2, 3 equal orbits; Jab, Jde written doubled.
For a state with J != 0 pass twoJ_state: the rank-0 reduced density is
divided by sqrt(2J+1), as EOM::ReadKshellRdm does for OBTD/TBTD, so the
values match a J_file = 0 header (M-averaged scalar density).
Header: J_file = 0, orbit table from the .snt (KSHELL order: protons first).
Only the 3-body block is written; OBTD/TBTD counts are 0.
"""
import math
import sys


def read_orbits(snt):
    lines = [l.split('!')[0].split() for l in open(snt)]
    lines = [l for l in lines if l]
    npo, nno = int(lines[0][0]), int(lines[0][1])
    orbs = []
    for l in lines[1:1 + npo + nno]:
        idx, n, lq, j2, tz = map(int, l[:5])
        orbs.append((idx, n, lq, j2, tz))
    return orbs


def g(a, b, c):
    s = {a, b, c}
    return {3: 1, 2: 2, 1: 6}[len(s)]


def main():
    snt, log, out = sys.argv[1:4]
    il, ir = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (1, 1)
    twoJ = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    orbs = read_orbits(snt)
    entries = []
    for line in open(log):
        if not line.startswith('TDTD:'):
            continue
        orb_part, j_part, st_part, val = line[5:].split(':')
        a, b, c, d, e, f = map(int, orb_part.split())
        jab, jde, j3b, j3k, rank = map(int, j_part.split())
        i, j = map(int, st_part.split())
        if (i, j) != (il, ir) or rank != 0:
            continue
        if j3b != j3k:
            raise ValueError("rank 0 with J3b != J3k: " + line)
        if not (a <= b <= c and d <= e <= f):
            raise ValueError("non-canonical triple: " + line)
        v = float(val) * math.sqrt(j3b + 1.0) * math.sqrt(g(a, b, c) * g(d, e, f))
        v /= math.sqrt(twoJ + 1.0)
        entries.append((a, b, c, d, e, f, 2 * jab, 2 * jde, j3b, v))
    with open(out, 'w') as fo:
        fo.write(".0\n%d\n" % len(orbs))
        for o in orbs:
            fo.write("%d %d %d %d %d\n" % o)
        fo.write("0\n0\n%d\n" % len(entries))
        for t in entries:
            fo.write("TRBTD %d %d %d %d %d %d %d %d %d %.12e\n" % t)
    print("wrote %d TRBTD entries to %s" % (len(entries), out))


if __name__ == '__main__':
    main()
