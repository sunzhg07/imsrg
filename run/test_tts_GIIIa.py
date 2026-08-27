#!/usr/bin/env python3
"""Γ^{III_a} ethS extract vs ladder gold.

ethS: AMC Path B χ^η (Pandya→DGEMM→inv) then Chi_AS×Γ ordinary-channel ladder.
Gold: test_GIIIa_ladder_mscheme.py (m ≡ from_chi ≡ Chi_AS DGEMM)

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIIIa.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-6
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
if not Eta.IsReduced():
    Eta.MakeReduced()
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def hat(J: float) -> float:
    return math.sqrt(2 * J + 1)


def jo(a: int) -> float:
    return ms.GetOrbit(a).j2 * 0.5


def j2i(a: int) -> int:
    return ms.GetOrbit(a).j2


def occ(a: int) -> float:
    return ms.GetOrbit(a).occ


def nbar(a: int) -> float:
    return 1.0 - occ(a)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def w_eta(a, b, k) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def tbme(J1, J2, a, b, c, d) -> float:
    return Eta.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


def gam_J(J0, a, b, c, d) -> float:
    return Gamma.TwoBody.GetTBME_J(J0, J0, a, b, c, d)


def chi_eta_red(i, j, k, l, J0) -> float:
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not (tri(ji, jj, J0) and tri(jk, jl, J0)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for b in orbits:
            w = w_eta(a, b, k)
            if abs(w) < 1e-12:
                continue
            ja = jo(a)
            for J2 in range(abs(j2i(i) - j2i(a)) // 2, (j2i(i) + j2i(a)) // 2 + 1):
                for J3 in range(
                    abs(j2i(b) - j2i(l)) // 2, (j2i(b) + j2i(l)) // 2 + 1
                ):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = tbme(J2, J3, i, a, b, l)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(b) - j2i(j)) // 2, (j2i(b) + j2i(j)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(k) - j2i(a)) // 2, (j2i(k) + j2i(a)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = tbme(J4, J5, b, j, k, a)
                            if abs(o2) < 1e-16:
                                continue
                            j0max = int(max(J2, J3, J4, J5, ji, jk, ja, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                s1 = SixJ(J3, lam, J2, ja, ji, j0)
                                s2 = SixJ(J4, lam, J5, ja, jk, j0)
                                n9 = NineJ(jl, jo(b), J3, jk, J4, j0, J0, jj, ji)
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                tot += (
                                    w
                                    * phase(J2 + J4 + lam)
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (2 * j0 + 1)
                                    / hat(lam)
                                    * s1
                                    * s2
                                    * n9
                                    * o1
                                    * o2
                                )
    return -phase(ji + jk) * hat(J0) * tot


print(f"emax={emax} λ={lam} seed={seed}")
print("Building χ^η_red (Path A ≡ Path B) ...")
t0 = time.time()
chi_red: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
                        continue
                    v = chi_eta_red(i, j, k, l, J0)
                    if abs(v) > 1e-16:
                        chi_red[(i, j, k, l, J0)] = v
print(f"  nonzero={len(chi_red)}  ({time.time() - t0:.2f}s)")


def matmul(A, B):
    n, m, p = len(A), len(A[0]), len(B[0])
    C = [[0.0] * p for _ in range(n)]
    for i in range(n):
        for k in range(m):
            aik = A[i][k]
            if abs(aik) < 1e-16:
                continue
            for j in range(p):
                C[i][j] += aik * B[k][j]
    return C


def transpose(A):
    return [list(row) for row in zip(*A)]


z_dgemm: dict = {}
for J0 in range(0, max_J + 1):
    pairs = [(p, q) for p in orbits for q in orbits if tri(jo(p), jo(q), J0)]
    if not pairs:
        continue
    n = len(pairs)
    idx = {pq: i for i, pq in enumerate(pairs)}
    Chi = [[0.0] * n for _ in range(n)]
    Gam = [[0.0] * n for _ in range(n)]
    for ia, (i, j) in enumerate(pairs):
        for ib, (a, b) in enumerate(pairs):
            Chi[ia][ib] = chi_red.get((i, j, a, b, J0), 0.0)
            Gam[ia][ib] = gam_J(J0, i, j, a, b)
    Chi_AS = [[0.0] * n for _ in range(n)]
    for ia, (i, j) in enumerate(pairs):
        pij = phase(J0 + jo(i) + jo(j))
        ji_idx = idx.get((j, i))
        for ib in range(n):
            if ji_idx is None:
                Chi_AS[ia][ib] = Chi[ia][ib]
            else:
                Chi_AS[ia][ib] = Chi[ia][ib] - pij * Chi[ji_idx][ib]
    T1 = matmul(Chi_AS, Gam)
    T2 = matmul(Gam, transpose(Chi_AS))
    for ia, (i, j) in enumerate(pairs):
        for ib, (k, l) in enumerate(pairs):
            v = -(T1[ia][ib] + T2[ia][ib])
            if abs(v) > 1e-16:
                z_dgemm[(i, j, k, l, J0)] = v

print("ethS GIIIa (AMC Path B χ^η → Chi_AS×Γ) ...")
t1 = time.time()
Z_fac = Operator(ms, 0, 0, 0, 2)
Z_fac.SetHermitian()
Commutator.FactorizedDoubleCommutator_eths.comm223_232_GIIIa(Eta, Gamma, Z_fac)
print(f"  ethS wall {time.time() - t1:.2f}s")

rats = Counter()
n = 0
max_ad = 0.0
for ch in range(ms.GetNumberTwoBodyChannels()):
    tbc = ms.GetTwoBodyChannel(ch)
    J0 = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        i, j = tbc.GetKet(ib).p, tbc.GetKet(ib).q
        for ik in range(nk):
            k, l = tbc.GetKet(ik).p, tbc.GetKet(ik).q
            z_red_eth = Z_fac.TwoBody.GetTBME_J(J0, J0, i, j, k, l) * hat(J0)
            z_red_g = z_dgemm.get((i, j, k, l, J0), 0.0)
            if abs(z_red_eth) < 1e-12 and abs(z_red_g) < 1e-12:
                continue
            n += 1
            max_ad = max(max_ad, abs(z_red_eth - z_red_g))
            if abs(z_red_g) > 1e-8:
                rats[round(z_red_eth / z_red_g, 6)] += 1

print(
    f"ethS ≡ ladder DGEMM: n={n}  max|Δ|={max_ad:.3e}  "
    f"eth/gold={rats.most_common(3)}"
)
print(f"Z_fac 2b norm = {Z_fac.TwoBodyNorm():.8e}")
ok = max_ad < tol
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
