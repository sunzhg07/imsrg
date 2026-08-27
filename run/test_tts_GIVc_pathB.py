#!/usr/bin/env python3
"""Bench ethS GIVc Path B (Pandya→DGEMM→inv) ≡ ring fold (≡ m).

λ≠0 only. Gold: tts_ring Path A + fermionic AS. ethS uses IMSRG
tensor Pandya (adcb map) → mid-J DGEMM → corrected inv → AS.

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIVc_pathB.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
seed = 11
tol = 1e-6

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1)


def occ(a):
    return ms.GetOrbit(a).occ


def nbar(a):
    return 1.0 - occ(a)


def jo(a):
    return ms.GetOrbit(a).j2 * 0.5


def j2i(a):
    return ms.GetOrbit(a).j2


def hat(J):
    return math.sqrt(2 * J + 1)


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def chi_amc_J(J0, J1, i, j, k, l) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w_l = nbar(a) * nbar(b) * occ(l) + occ(a) * occ(b) * nbar(l)
            w_j = nbar(a) * nbar(b) * occ(j) + occ(a) * occ(b) * nbar(j)
            if abs(w_l) > 1e-12 and tri(jo(a), jo(b), J0):
                sm += w_l * (
                    Gamma.TwoBody.GetTBME_J(J0, J0, i, j, a, b)
                    * Eta.TwoBody.GetTBME_J(J0, J1, a, b, k, l)
                )
            if abs(w_j) > 1e-12 and tri(jo(a), jo(b), J1):
                sm += w_j * (
                    Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b)
                    * Gamma.TwoBody.GetTBME_J(J1, J1, a, b, k, l)
                )
    return sm


print("Building χ^λ cache ...")
t0 = time.time()
chiJ: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        v = chi_amc_J(J0, J1, i, j, k, l)
                        if abs(v) > 1e-16:
                            chiJ[(J0, J1, i, j, k, l)] = v
print(f"  nonzero={len(chiJ)}  ({time.time()-t0:.2f}s)")


def ring_X(p, q, s, r, J0) -> float:
    jp, jq, js, jr = jo(p), jo(q), jo(s), jo(r)
    if not (tri(jp, jq, J0) and tri(js, jr, J0)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for b in orbits:
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(p) - j2i(b)) // 2, (j2i(p) + j2i(b)) // 2 + 1):
                for J3 in range(abs(j2i(a) - j2i(r)) // 2, (j2i(a) + j2i(r)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = chiJ.get((J2, J3, p, b, a, r), 0.0)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(a) - j2i(q)) // 2, (j2i(a) + j2i(q)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(s) - j2i(b)) // 2, (j2i(s) + j2i(b)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = Eta.TwoBody.GetTBME_J(J4, J5, a, q, s, b)
                            if abs(o2) < 1e-16:
                                continue
                            j0max = int(max(J2, J3, J4, J5, jp, js, jb, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                s1 = SixJ(J3, lam, J2, jb, jp, j0)
                                s2 = SixJ(J4, lam, J5, jb, js, j0)
                                n9 = NineJ(jr, ja, J3, js, J4, j0, J0, jq, jp)
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                tot += (
                                    phase(J2 + J4 + lam)
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (2 * j0 + 1)
                                    * hat_lam_inv
                                    * s1
                                    * s2
                                    * n9
                                    * o1
                                    * o2
                                )
    return -phase(jp + js) * hat(J0) * tot


def amc_fold_red(i, j, k, l, J0) -> float:
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    x = ring_X(i, j, k, l, J0)
    xkl = ring_X(i, j, l, k, J0)
    xij = ring_X(j, i, k, l, J0)
    xijkl = ring_X(j, i, l, k, J0)
    pkl = phase(jo(k) + jo(l) - J0)
    pij = phase(jo(i) + jo(j) - J0)
    return 0.5 * (x - pkl * xkl - pij * xij + pij * pkl * xijkl)


def fresh_Z():
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    return Z


print("ethS Path B ...")
t0 = time.time()
cm = Commutator.FactorizedDoubleCommutator_eths
Z_B = fresh_Z()
cm.comm223_232_GIVc(Eta, Gamma, Z_B)
print(f"  ‖Z_B‖={Z_B.TwoBodyNorm():.6g}  ({time.time()-t0:.2f}s)")

print("Compare GetTBME_J vs ring fold unreduced (channel upper triangle) ...")
t0 = time.time()
n = 0
max_abs = 0.0
max_rel = 0.0
sum_sq = 0.0
nch = ms.GetNumberTwoBodyChannels()
for ch in range(nch):
    tbc = ms.GetTwoBodyChannel(ch)
    J0 = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        bra = tbc.GetKet(ib)
        i, j = bra.p, bra.q
        for ik in range(ib, nk):
            ket = tbc.GetKet(ik)
            k, l = ket.p, ket.q
            aU = amc_fold_red(i, j, k, l, J0) / hat(J0)
            bU = Z_B.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
            if abs(aU) < 1e-10 and abs(bU) < 1e-10:
                continue
            n += 1
            d = abs(aU - bU)
            max_abs = max(max_abs, d)
            den = max(abs(aU), abs(bU), 1e-30)
            max_rel = max(max_rel, d / den)
            sum_sq += d * d
print(
    f"  n={n}  max|Δ|={max_abs:.3e}  max|rel|={max_rel:.3e}"
    f"  rms={math.sqrt(sum_sq / max(n, 1)):.3e}  ({time.time()-t0:.2f}s)"
)

ok = max_abs < tol
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
