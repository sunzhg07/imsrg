#!/usr/bin/env python3
"""Bare χ^ι: m-scheme ↔ AMC direct (no (1−P), no RC / second Ω).

Trusted normal product (η-analog; Factorized bar_CHI_V twin):
  χ^ι_ijkl = Σ_ab (n̄_a n_b n̄_k + n_a n̄_b n_k) Ω_bjka Γ_iabl

Arxiv / older analyze Ω_bika Γ_iabl is **wrong** (spectator j; Yutsis fails).
It only matches this form on fold index patterns (χ_aibk, χ_akbi).

Packaging (same as χ^κ):
  Ω  : WE-reduced tensor for **all** λ (incl. 0)
  Γ  : unreduced scalar
  χ  : AMC tensor LHS as printed
  m  : physical

Phases: combined-integer (χ^β / κ rule).

AMC: learn/amc_tts/factored_GIV/output/chi_iota_analyze.tex
  (Ω·Γ order → κ-like print with mid j0)

Usage:
  PYTHONPATH=build python3 -B run/test_chi_iota_m_vs_amc.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import random
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-5
seed = 11
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 40

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
if lam != 0:
    ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()
if not Eta.IsReduced():
    Eta.MakeReduced()

orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
print(
    f"emax={emax} λ={lam} seed={seed} nsamp={nsamp}\n"
    f"  Ω reduced={Eta.IsReduced()} Jrank={Eta.GetJRank()} ||2b||={Eta.TwoBodyNorm():.4g}"
    f"  [tensor/WE path for all λ]\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  AMC χ^ι: Ω_bjka Γ_iabl (η-analog; arxiv Ω_bika discarded)"
)


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


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


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def w_occ(a: int, b: int, k: int) -> float:
    """occ_AbarBC on (a,b,k): n̄_a n_b n̄_k + n_a n̄_b n_k."""
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    w = w_occ(a, b, k)
                    if abs(w) < 1e-12:
                        continue
                    g = ut.GetMschemeMatrixElement_2b(
                        Gamma, i, mi, a, ma, b, mb, l, ml
                    )
                    if abs(g) < 1e-16:
                        continue
                    o = ut.GetMschemeMatrixElement_2b(
                        Eta, b, mb, j, mj, k, mk, a, ma
                    )
                    sm += w * g * o
    return sm


def chi_amc(J0: int, J1: int, i: int, j: int, k: int, l: int) -> float:
    """chi_iota_analyze.tex (Ω·Γ order) with combined-integer phases.

    −(−1)^{J0+j_i}(−1)^{J2+j_b} → −iphase(J0+J2+(j2_i+j2_b)//2)
    """
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    tot = 0.0
    for a in orbits:
        for b in orbits:
            w = w_occ(a, b, k)
            if abs(w) < 1e-12:
                continue
            ja, jb = jo(a), jo(b)
            # Ω_bjka^{J2 J3 λ}
            for J2 in range(
                abs(j2i(b) - j2i(j)) // 2, (j2i(b) + j2i(j)) // 2 + 1
            ):
                for J3 in range(
                    abs(j2i(k) - j2i(a)) // 2, (j2i(k) + j2i(a)) // 2 + 1
                ):
                    if not tri(J2, J3, lam):
                        continue
                    om = Eta.TwoBody.GetTBME_J(J2, J3, b, j, k, a)
                    if abs(om) < 1e-16:
                        continue
                    # Γ_iabl^{J4 J4}
                    J4min = max(abs(j2i(i) - j2i(a)), abs(j2i(b) - j2i(l))) // 2
                    J4max = min(j2i(i) + j2i(a), j2i(b) + j2i(l)) // 2
                    for J4 in range(J4min, J4max + 1):
                        gam = Gamma.TwoBody.GetTBME_J(J4, J4, i, a, b, l)
                        if abs(gam) < 1e-16:
                            continue
                        for J5 in range(0, max_J + 2):
                            j0max = int(max(J2, J3, J0, J1, jj, jb, lam) + 3)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                six = (
                                    SixJ(J3, lam, J2, jj, jb, j0)
                                    * SixJ(J1, lam, J0, jj, ji, j0)
                                    * SixJ(jb, ja, J5, jk, j0, J3)
                                    * SixJ(jl, ji, J5, ja, jb, J4)
                                    * SixJ(jk, jl, J1, ji, j0, J5)
                                )
                                if abs(six) < 1e-16:
                                    continue
                                # Ĵ2 Ĵ3 Ĵ4² Ĵ5² ĵ0²
                                hats = (
                                    hat(J2)
                                    * hat(J3)
                                    * (2 * J4 + 1)
                                    * (2 * J5 + 1)
                                    * (2 * j0 + 1)
                                )
                                ph = -iphase(
                                    J0 + J2 + (j2i(i) + j2i(b)) // 2
                                )
                                tot += ph * hats * six * w * om * gam
    return hat(J0) * hat(J1) * tot


print("Caching χ_AMC as printed ...")
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
                        if not (
                            tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)
                        ):
                            continue
                        v = chi_amc(J0, J1, i, j, k, l)
                        if abs(v) < 1e-16:
                            continue
                        chiJ[(J0, J1, i, j, k, l)] = v
                        nfill += 1
print(f"  nonzero={nfill}  ({time.time()-t0:.1f}s)")


def we_unpack(i, mi, j, mj, k, mk, l, ml) -> float:
    """UnitTest WE unpack. No Pauli skip (χ not AS)."""
    if abs(mi + mj - mk - ml) > 2 * lam:
        return 0.0
    M0 = (mi + mj) // 2
    M1 = (mk + ml) // 2
    mu = M0 - M1
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0:
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        for J1 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if abs(M1) > J1:
                continue
            if not tri(J0, J1, lam):
                continue
            ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            if lam == 0:
                cj = 1.0 if (J0 == J1 and M0 == M1) else 0.0
            else:
                cj = CG(J1, M1, lam, mu, J0, M0)
            if abs(cj) < 1e-15:
                continue
            v = chiJ.get((J0, J1, i, j, k, l), 0.0)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


def chi_red_from_m(J0: int, J1: int, i: int, j: int, k: int, l: int) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    sm = 0.0
    for mi in mrange(i):
        for mj in mrange(j):
            M0 = (mi + mj) // 2
            if abs(M0) > J0:
                continue
            cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
            if abs(cab) < 1e-15:
                continue
            for mk in mrange(k):
                for ml in mrange(l):
                    M1 = (mk + ml) // 2
                    if abs(M1) > J1:
                        continue
                    mu = M0 - M1
                    if abs(mu) > lam:
                        continue
                    ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J1, M1)
                    if abs(ccd) < 1e-15:
                        continue
                    if lam == 0:
                        cj = 1.0 if (J0 == J1 and M0 == M1) else 0.0
                    else:
                        cj = CG(J1, M1, lam, mu, J0, M0)
                    if abs(cj) < 1e-15:
                        continue
                    cm = chi_m(i, mi, j, mj, k, mk, l, ml)
                    if abs(cm) < 1e-16:
                        continue
                    sm += cj * cab * ccd * cm
    return sm / hat(J0)


# --- (A) m ≡ WE(AMC) ---
random.seed(seed)
cands = []
for i in orbits:
    for mi in mrange(i):
        for j in orbits:
            for mj in mrange(j):
                for k in orbits:
                    for mk in mrange(k):
                        for l in orbits:
                            for ml in mrange(l):
                                if abs(mi + mj - mk - ml) > 2 * lam:
                                    continue
                                if not any(
                                    abs(w_occ(a, b, k)) > 1e-12
                                    for a in orbits
                                    for b in orbits
                                ):
                                    continue
                                cands.append((i, mi, j, mj, k, mk, l, ml))
random.shuffle(cands)

samples = []
for tup in cands:
    if len(samples) >= nsamp:
        break
    cm = chi_m(*tup)
    if abs(cm) < 1e-8:
        continue
    samples.append((tup, cm))

print(f"\n(A) m ≡ WE(χ_AMC)  [{len(samples)} samples]")
rats_A = Counter()
max_A = 0.0
for tup, cm in samples:
    ca = we_unpack(*tup)
    max_A = max(max_A, abs(cm - ca))
    if abs(ca) > 1e-8:
        rats_A[round(cm / ca, 4)] += 1
print(f"  max|m−WE|={max_A:.3e}  m/WE={rats_A.most_common(4)}")

# --- (B) χ_red←m ≡ AMC ---
print("\n(B) χ_red←m (/Ĵ0) ≡ χ_AMC")
t0 = time.time()
rats_B = Counter()
max_B = 0.0
nB = 0
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not (
                            tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)
                        ):
                            continue
                        cr = chi_red_from_m(J0, J1, i, j, k, l)
                        ca = chiJ.get((J0, J1, i, j, k, l), 0.0)
                        if abs(cr) < 1e-10 and abs(ca) < 1e-10:
                            continue
                        nB += 1
                        max_B = max(max_B, abs(cr - ca))
                        if abs(ca) > 1e-8:
                            rats_B[round(cr / ca, 4)] += 1
                        if nB >= 40:
                            break
                    if nB >= 40:
                        break
                if nB >= 40:
                    break
            if nB >= 40:
                break
        if nB >= 40:
            break
    if nB >= 40:
        break
print(
    f"  n={nB}  max|red−A|={max_B:.3e}  "
    f"red/A={rats_B.most_common(4)}  ({time.time()-t0:.1f}s)"
)

ok = max_A < tol and max_B < tol
print()
if ok:
    print(
        "PASS — m χ^ι ≡ AMC (η-analog Ω_bjka Γ_iabl, combined-integer phases)"
    )
    sys.exit(0)
print(
    f"FAIL — max|m−WE|={max_A:.3e} max|red−A|={max_B:.3e}. "
    "If ratios ~ Ĵ^{±1}, packaging; if r=−1, check half-int phase combining."
)
sys.exit(1)
