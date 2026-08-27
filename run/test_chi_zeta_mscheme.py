#!/usr/bin/env python3
"""χ^ζ (tensor 1b rank λ): m-scheme ≡ AMC J (reduced tensor).

AMC: learn/amc_tts/factored_GII/input/chi_zeta.txt
  χ_ij^λ = 1/2 (−1)^{jj+λ} Σ_{abc J0 J1} w (−1)^{J0+ja}
           Ĵ0 Ĵ1 {λ J1 J0; ja ji jj} Γ_{aibc}^{J0} Ω_{bcaj}^{J0 J1 λ}

m: χ_ij(mi,mj) = 1/2 Σ_abc w Γ_aibc(m) Ω_bcaj(m)   (Γ×Ω→χ^λ, no extra CG)

Γ: unreduced preferred (MakeNotReduced if needed). Ω: reduced tensor.

Usage:
  PYTHONPATH=build python3 run/test_chi_zeta_mscheme.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-7
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω: reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ: reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}"
)


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


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def m_range(o):
    return range(-j2i(o), j2i(o) + 1, 2)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def chi_zeta_m(i, mi, j, mj) -> float:
    sm = 0.0
    for a in orbits:
        na, nba = occ(a), nbar(a)
        for ma in m_range(a):
            for b in orbits:
                nb, nbb = occ(b), nbar(b)
                for mb in m_range(b):
                    for c in orbits:
                        nc, nbc = occ(c), nbar(c)
                        w = na * nb * nbc + nba * nbb * nc
                        if abs(w) < 1e-12:
                            continue
                        for mc in m_range(c):
                            g = ut.GetMschemeMatrixElement_2b(
                                Gamma, a, ma, i, mi, b, mb, c, mc
                            )
                            if abs(g) < 1e-16:
                                continue
                            o = ut.GetMschemeMatrixElement_2b(
                                Eta, b, mb, c, mc, a, ma, j, mj
                            )
                            sm += w * g * o
    return 0.5 * sm


def chi_zeta_J(i, j) -> float:
    """AMC reduced χ^ζ RME (unreduced Γ)."""
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    max_J = max(j2i(o) for o in orbits)
    sm = 0.0
    for a in orbits:
        ja = jo(a)
        for b in orbits:
            for c in orbits:
                w = occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)
                if abs(w) < 1e-12:
                    continue
                for J0 in range(0, max_J + 1):
                    if not (tri(ja, ji, J0) and tri(jo(b), jo(c), J0)):
                        continue
                    g = Gamma.TwoBody.GetTBME_J(J0, J0, a, i, b, c)
                    if abs(g) < 1e-16:
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam) or not tri(ja, jj, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, b, c, a, j)
                        if abs(o) < 1e-16:
                            continue
                        six = SixJ(lam, J1, J0, ja, ji, jj)
                        if abs(six) < 1e-16:
                            continue
                        # (−1)^{jj+λ+J0+ja} — integer form like χ^β
                        ph = iphase((j2i(j) + j2i(a)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * g * o
    return 0.5 * sm


print("\nBuilding reduced χ^ζ_J ...")
t0 = time.time()
Chi = Operator(ms, lam, 0, 0, 1)
Chi.SetNonHermitian()
nch = 0
for i in orbits:
    for j in orbits:
        if not tri(jo(i), jo(j), lam):
            continue
        Chi.SetOneBody(i, j, chi_zeta_J(i, j))
        nch += 1
print(
    f"  channels={nch} reduced={Chi.IsReduced()} "
    f"||1b||={Chi.OneBodyNorm():.6e}  ({time.time()-t0:.2f}s)"
)

print("\nCompare m_lit vs GetMschemeMatrixElement_1b(Chi_J) ...")
t0 = time.time()
max_abs = 0.0
n = 0
rats = Counter()
worst = None
for i in orbits:
    for j in orbits:
        if not tri(jo(i), jo(j), lam):
            continue
        for mi in m_range(i):
            for mj in m_range(j):
                if abs(mi - mj) > 2 * lam:
                    continue
                zm = chi_zeta_m(i, mi, j, mj)
                zj = ut.GetMschemeMatrixElement_1b(Chi, i, mi, j, mj)
                if abs(zm) < 1e-14 and abs(zj) < 1e-14:
                    continue
                n += 1
                err = abs(zm - zj)
                if err > max_abs:
                    max_abs = err
                    worst = (i, j, mi, mj, zm, zj, err)
                if abs(zj) > 1e-12:
                    rats[round(zm / zj, 6)] += 1

print(
    f"  n={n}  max|m−J_m|={max_abs:.3e}  m/J_m={rats.most_common(3)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst:
    i, j, mi, mj, zm, zj, err = worst
    print(
        f"  worst ij=({i},{j}) m=({mi},{mj}): "
        f"m={zm:.6e} J_m={zj:.6e} Δ={err:.3e}"
    )

ok = max_abs < tol
print("\nPASS — χ^ζ m ≡ AMC J" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
