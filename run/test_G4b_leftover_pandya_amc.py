#!/usr/bin/env python3
"""Γ^{IV_b} leftover W1: AMC CC-GEMM chain ≡ AMC five-6j gold.

Pure angular lock (no ModelSpace): random χ^{J2J3}(ai;bk), Ω^{J4J5}(jb;la)
tables on one (i,j,k,l,a,b) set of j-values.

  gold  learn/amc_tts/factored_GIV/output/G4b_W1_leftover.tex   (five 6j + j0)
  new   learn/amc_tts/factored_GIV/output/G4b_leftover_pandya.tex
          eq1 barChi (Pandya (ik)|(ab))   eq2 barO (Pandya (jl)|(ab))
          eq3 pure GEMM over (ab, J_ab)   eq4 inverse 6j

AMC `--collect-ninejs` drops the (-1)^{2x} of the 3x6j→9j identity
(YutsisGraph.sign is never emitted).  x = j0 is half-integer here, so
  sum_j0 (2j0+1) 6j 6j 6j = - NineJ(...)
and each Pandya below carries that minus explicitly.

  PYTHONPATH=build python3 run/test_G4b_leftover_pandya_amc.py [lambda=2] [ntrial=6]
"""

from __future__ import annotations

import math
import random
import sys

from pyIMSRG import NineJ, SixJ

lam = int(sys.argv[1]) if len(sys.argv) > 1 else 2
ntrial = int(sys.argv[2]) if len(sys.argv) > 2 else 6


def hat(J):
    return math.sqrt(2 * J + 1)


def iph(n):
    return 1.0 if int(round(n)) % 2 == 0 else -1.0


def tri(a, b, c):
    return abs(a - b) <= c <= a + b


def jrange(ja, jb):
    return range(int(round(abs(ja - jb))), int(round(ja + jb)) + 1)


def gold_W1(js, chi, om, J0):
    """AMC G4b_W1_leftover.tex (five 6j, dummy j0)."""
    ji, jj, jk, jl, ja, jb = js
    tot = 0.0
    for (J2, J3), c in chi.items():
        for (J4, J5), o in om.items():
            for J6 in jrange(ji, jk):
                if not tri(jl, jj, J6):
                    continue
                s5 = SixJ(jk, jl, J0, jj, ji, J6)
                if abs(s5) < 1e-14:
                    continue
                acc = 0.0
                j0 = 0.5
                while j0 <= ja + lam + 1e-9:
                    if tri(ja, lam, j0):
                        acc += (
                            (2 * j0 + 1)
                            * SixJ(J3, lam, J2, ja, ji, j0)
                            * SixJ(J4, lam, J5, ja, jl, j0)
                            * SixJ(ji, jk, J6, jb, j0, J3)
                            * SixJ(jl, jj, J6, jb, j0, J4)
                        )
                    j0 += 1.0
                if abs(acc) < 1e-16:
                    continue
                tot += (
                    iph(ja + jb + lam)
                    * hat(J2) * hat(J3) * hat(J4) * hat(J5)
                    * (2 * J6 + 1) / hat(lam)
                    * acc * s5 * c * o
                )
    return -iph(J0 + jj + jk) * tot


def bar_chi(js, chi, J_ik, J_ab):
    """eq1 as a 9j (AMC minus x collect-ninejs minus = +)."""
    ji, jj, jk, jl, ja, jb = js
    tot = 0.0
    for (J2, J3), c in chi.items():
        nj = NineJ(lam, J_ik, J_ab, J3, jk, jb, J2, ji, ja)
        if abs(nj) < 1e-14:
            continue
        tot += hat(J2) * hat(J3) * nj * c
    return iph(J_ik + J_ab + lam) * hat(J_ik) * hat(J_ab) * tot


def bar_om(js, om, J_jl, J_ab):
    """eq2 as a 9j (AMC plus x collect-ninejs minus = -)."""
    ji, jj, jk, jl, ja, jb = js
    tot = 0.0
    for (J4, J5), o in om.items():
        nj = NineJ(lam, J_jl, J_ab, J5, jl, ja, J4, jj, jb)
        if abs(nj) < 1e-14:
            continue
        tot += iph(J4 + J5) * hat(J4) * hat(J5) * nj * o
    return -iph(J_jl + jj + jl + lam) * hat(J_jl) * hat(J_ab) * tot


def new_W1(js, chi, om, J0):
    """eq3 GEMM then eq4 inverse."""
    ji, jj, jk, jl, ja, jb = js
    tot = 0.0
    for J2 in jrange(ji, jk):
        if not tri(jj, jl, J2):
            continue
        s = SixJ(jk, jl, J0, jj, ji, J2)
        if abs(s) < 1e-14:
            continue
        barW = 0.0
        for J_ab in jrange(ja, jb):
            if not tri(J2, J_ab, lam):
                continue
            barW += bar_chi(js, chi, J2, J_ab) * bar_om(js, om, J2, J_ab)
        barW /= hat(J2) * hat(lam)
        tot += iph(J2) * hat(J2) * s * barW
    return iph(J0 + jj + jk) * tot


rng = random.Random(7)
print(f"lambda={lam}  ntrial={ntrial}   gold(5x6j) vs CC-GEMM chain")
n_ok = n_bad = 0
for t in range(ntrial):
    js = [0.5 + rng.randint(0, 3) for _ in range(6)]
    ji, jj, jk, jl, ja, jb = js
    chi = {}
    for J2 in jrange(ja, ji):
        for J3 in jrange(jb, jk):
            if tri(J2, J3, lam):
                chi[(J2, J3)] = rng.uniform(-1, 1)
    om = {}
    for J4 in jrange(jj, jb):
        for J5 in jrange(jl, ja):
            if tri(J4, J5, lam):
                om[(J4, J5)] = rng.uniform(-1, 1)
    if not chi or not om:
        continue
    for J0 in jrange(ji, jj):
        if not tri(jk, jl, J0):
            continue
        g = gold_W1(js, chi, om, J0)
        n = new_W1(js, chi, om, J0)
        if abs(g) < 1e-10 and abs(n) < 1e-10:
            continue
        r = n / g if abs(g) > 1e-12 else float("nan")
        ok = abs(n - g) < 1e-9 * max(1.0, abs(g))
        n_ok += ok
        n_bad += not ok
        if not ok or t == 0:
            print(
                f"  j=({ji},{jj},{jk},{jl};{ja},{jb}) J0={J0}  "
                f"gold={g:12.5e} new={n:12.5e} r={r:.6f} "
                f"{'OK' if ok else 'FAIL'}"
            )

print(f"\n{n_ok} OK  {n_bad} FAIL")
print("PASS" if n_bad == 0 and n_ok > 0 else "FAIL")
sys.exit(0 if (n_bad == 0 and n_ok > 0) else 1)
