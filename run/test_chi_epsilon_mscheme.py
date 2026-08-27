#!/usr/bin/env python3
"""χ^ε (scalar 1b): m-scheme ≡ AMC J (unreduced ĵ^{-2}).

AMC: learn/amc_tts/factored_GI/input/chi_epsilon.txt
  χ_ij^0 = 1/2 δ_{jj,ji} ĵ_i^{-2} Σ_{abc J0 J1 λ}
           w (−1)^{J0+J1+λ} λ̂^{-1} Ω_{ciab}^{J0 J1 λ} Ω_{abcj}^{J1 J0 λ}

m (always unreduced) + [Ω^λ × Ω^λ]^(0):
  χ_ij(mi,mj) = 1/2 Σ w Ω_ciab(μ) Ω_abcj(−μ) CG(λ μ, λ −μ; 0 0)

Usage:
  PYTHONPATH=build python3 run/test_chi_epsilon_mscheme.py [emax=1] [lambda=2]
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
orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω: reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}"
)


def occ(a):
    return ms.GetOrbit(a).occ


def nbar(a):
    return 1.0 - occ(a)


def j2i(a):
    return ms.GetOrbit(a).j2


def m_range(o):
    return range(-j2i(o), j2i(o) + 1, 2)


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def hat_lam_inv():
    return 1.0 / math.sqrt(2 * lam + 1)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def cg0(mu: int) -> float:
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


def chi_eps_m(i, mi, j, mj) -> float:
    """m gold: analyze χ^ε + [Ω^λ×Ω^λ]^(0) (same CG pattern as f^I)."""
    if mi != mj or j2i(i) != j2i(j):
        return 0.0
    sm = 0.0
    for a in orbits:
        na, nba = occ(a), nbar(a)
        for ma in m_range(a):
            for b in orbits:
                nb, nbb = occ(b), nbar(b)
                for mb in m_range(b):
                    for c in orbits:
                        nc, nbc = occ(c), nbar(c)
                        w = nba * nbb * nc + na * nb * nbc
                        if abs(w) < 1e-12:
                            continue
                        for mc in m_range(c):
                            o1 = ut.GetMschemeMatrixElement_2b(
                                Eta, c, mc, i, mi, a, ma, b, mb
                            )
                            if abs(o1) < 1e-16:
                                continue
                            # m indices are 2m; μ = ΔM on first Ω
                            mu = (mc + mi - ma - mb) // 2
                            cg = cg0(mu)
                            if abs(cg) < 1e-16:
                                continue
                            o2 = ut.GetMschemeMatrixElement_2b(
                                Eta, a, ma, b, mb, c, mc, j, mj
                            )
                            if abs(o2) < 1e-16:
                                continue
                            sm += w * cg * o1 * o2
    return 0.5 * sm


def chi_eps_J(i, j) -> float:
    """AMC unreduced χ^ε (ĵ_i^{-2})."""
    if j2i(i) != j2i(j):
        return 0.0
    max_J = max(j2i(o) for o in orbits)
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = nbar(a) * nbar(b) * occ(c) + occ(a) * occ(b) * nbar(c)
                if abs(w) < 1e-12:
                    continue
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        o1 = Eta.TwoBody.GetTBME_J(J0, J1, c, i, a, b)
                        o2 = Eta.TwoBody.GetTBME_J(J1, J0, a, b, c, j)
                        if abs(o1 * o2) < 1e-16:
                            continue
                        sm += (
                            w
                            * iphase(J0 + J1 + lam)
                            * hat_lam_inv()
                            * o1
                            * o2
                        )
    return 0.5 * sm / (j2i(i) + 1.0)


print("\nBuilding unreduced χ^ε_J (AMC) ...")
t0 = time.time()
Chi = Operator(ms, 0, 0, 0, 1)
Chi.SetHermitian()
# scalar 1b default not reduced
assert not Chi.IsReduced()
for i in orbits:
    for j in orbits:
        if j2i(i) != j2i(j):
            continue
        Chi.SetOneBody(i, j, chi_eps_J(i, j))
print(
    f"  ||1b||={Chi.OneBodyNorm():.6e}  ({time.time()-t0:.2f}s)"
)

print("\nCompare m_lit vs GetMschemeMatrixElement_1b(Chi) ...")
t0 = time.time()
max_abs = 0.0
n = 0
rats = Counter()
worst = None
for i in orbits:
    for j in orbits:
        if j2i(i) != j2i(j):
            continue
        for mi in m_range(i):
            zm = chi_eps_m(i, mi, j, mi)
            zj = ut.GetMschemeMatrixElement_1b(Chi, i, mi, j, mi)
            if abs(zm) < 1e-14 and abs(zj) < 1e-14:
                continue
            n += 1
            err = abs(zm - zj)
            if err > max_abs:
                max_abs = err
                worst = (i, j, mi, zm, zj, err)
            if abs(zj) > 1e-12:
                rats[round(zm / zj, 6)] += 1

print(
    f"  n={n}  max|m−J_m|={max_abs:.3e}  m/J_m={rats.most_common(3)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst:
    i, j, mi, zm, zj, err = worst
    print(
        f"  worst ij=({i},{j}) m={mi}: m={zm:.6e} J_m={zj:.6e} Δ={err:.3e}"
    )

ok = max_abs < tol
print("\nPASS — χ^ε m ≡ AMC J (unreduced)" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
