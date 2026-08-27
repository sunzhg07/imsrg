#!/usr/bin/env python3
"""χ^λ (T×S / S×T → tensor): m-scheme ≡ AMC.

AMC: learn/amc_tts/factored_GIV/input/chi_lambda.txt
  Term1: χ^{J0 J1 λ} = Σ_ab w_l Γ^{J0}_{ijab} Ω^{J0 J1 λ}_{abkl}
  Term2: χ^{J0 J1 λ} = Σ_ab w_j Ω^{J0 J1 λ}_{ijab} Γ^{J1}_{abkl}
  (unreduced Γ; no 6j). Γ reduce=true → extra Ĵ0^{-1}/Ĵ1^{-1}.

m (physical T×S→T — no [T×T]^(0) CG):
  χ(m) = Σ_ab [w_l Γ_ijab(m) Ω_abkl(m) + w_j Ω_ijab(m) Γ_abkl(m)]

Compare via WE unpack of AMC χ_J (same formula as UnitTest::GetMschemeMatrixElement_2b):
  m = Σ_{J0 J1} CG(J1 M1; λ μ; J0 M0) CG_ab CG_cd / Ĵ0 × χ^{J0 J1}

Usage:
  PYTHONPATH=build python3 -B run/test_chi_lambda_mscheme.py [emax=1] [lambda=2]
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
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)  # Ω tensor
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)  # Γ scalar
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}"
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


def m_range(o):
    return range(-j2i(o), j2i(o) + 1, 2)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


max_J = max(j2i(o) for o in orbits)


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """Physical m: ΓΩ + ΩΓ (T×S→T)."""
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            for b in orbits:
                for mb in m_range(b):
                    w_l = (
                        nbar(a) * nbar(b) * occ(l)
                        + occ(a) * occ(b) * nbar(l)
                    )
                    w_j = (
                        nbar(a) * nbar(b) * occ(j)
                        + occ(a) * occ(b) * nbar(j)
                    )
                    if abs(w_l) > 1e-12:
                        g = ut.GetMschemeMatrixElement_2b(
                            Gamma, i, mi, j, mj, a, ma, b, mb
                        )
                        if abs(g) > 1e-16:
                            o = ut.GetMschemeMatrixElement_2b(
                                Eta, a, ma, b, mb, k, mk, l, ml
                            )
                            sm += w_l * g * o
                    if abs(w_j) > 1e-12:
                        o = ut.GetMschemeMatrixElement_2b(
                            Eta, i, mi, j, mj, a, ma, b, mb
                        )
                        if abs(o) > 1e-16:
                            g = ut.GetMschemeMatrixElement_2b(
                                Gamma, a, ma, b, mb, k, mk, l, ml
                            )
                            sm += w_j * o * g
    return sm


def chi_amc_J(J0, J1, i, j, k, l) -> float:
    """AMC chi_lambda.tex (unreduced Γ): Term1+Term2."""
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
                g = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, a, b)
                o = Eta.TwoBody.GetTBME_J(J0, J1, a, b, k, l)
                sm += w_l * g * o
            if abs(w_j) > 1e-12 and tri(jo(a), jo(b), J1):
                o = Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b)
                g = Gamma.TwoBody.GetTBME_J(J1, J1, a, b, k, l)
                sm += w_j * o * g
    return sm


# Cache χ_J
print("Building χ^λ_J cache from AMC ...")
t0 = time.time()
chiJ: dict = {}
nfill = 0
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        v = chi_amc_J(J0, J1, i, j, k, l)
                        if abs(v) < 1e-16:
                            continue
                        chiJ[(J0, J1, i, j, k, l)] = v
                        nfill += 1
print(f"  nonzero={nfill}  ({time.time()-t0:.2f}s)")


def chi_J_to_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """WE unpack of reduced tensor χ (UnitTest formula)."""
    if i == j and mi == mj:
        return 0.0
    if k == l and mk == ml:
        return 0.0
    if abs(mi + mj - mk - ml) > 2 * lam:
        return 0.0
    M0 = (mi + mj) // 2
    M1 = (mk + ml) // 2
    mu = M0 - M1
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0:
            continue
        if i == j and J0 % 2 > 0:
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        for J1 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if abs(M1) > J1:
                continue
            if k == l and J1 % 2 > 0:
                continue
            if not tri(J0, J1, lam):
                continue
            ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            cj = CG(J1, M1, lam, mu, J0, M0)
            if abs(cj) < 1e-15:
                continue
            v = chiJ.get((J0, J1, i, j, k, l), 0.0)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


print("\nCompare m_gold vs WE(AMC χ_J) ...")
t0 = time.time()
rats = Counter()
n = 0
max_abs = 0.0
worst = None

for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                if len({i, j, k, l}) < 3:
                    continue
                for mi in m_range(i):
                    for mj in m_range(j):
                        for mk in m_range(k):
                            for ml in m_range(l):
                                if abs(mi + mj - mk - ml) > 2 * lam:
                                    continue
                                m = chi_m(i, mi, j, mj, k, mk, l, ml)
                                jme = chi_J_to_m(i, mi, j, mj, k, mk, l, ml)
                                if abs(m) < 1e-8 and abs(jme) < 1e-8:
                                    continue
                                n += 1
                                d = abs(m - jme)
                                if d > max_abs:
                                    max_abs = d
                                    worst = (
                                        i,
                                        j,
                                        k,
                                        l,
                                        (mi, mj, mk, ml),
                                        m,
                                        jme,
                                    )
                                if abs(jme) > 1e-8:
                                    rats[round(m / jme, 4)] += 1

print(
    f"  n={n}  max|m−J|={max_abs:.3e}  m/J={rats.most_common(6)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst is not None:
    i, j, k, l, ms_, m, jme = worst
    print(
        f"  worst ({i},{j},{k},{l}) m={ms_}: m={m:.6e} J={jme:.6e} "
        f"r={m/jme if abs(jme)>1e-14 else float('nan'):g}"
    )

ok = max_abs < tol and (not rats or rats.most_common(1)[0][0] == 1.0)
print("\nPASS — χ^λ m ≡ AMC" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
