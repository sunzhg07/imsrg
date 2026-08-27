#!/usr/bin/env python3
"""Bare χ^κ: m-scheme ↔ AMC direct (no (1−P), no suspect ±/Ĵ fudges).

Packaging — name both sides first (REDUCED_UNREDUCED.md):

  Ω  : reduced tensor (IsReduced) for **all** λ including 0 — WE path
       (λ=0 does NOT mean unreduced scalar)
  Γ  : unreduced scalar
  χ  : AMC tensor LHS (scalar=false). Formula taken **as printed**
       (input reduce=true does not change this topology’s .tex)
  m  : always physical

Correct compares (do not invent hats/signs):
  (A) m ≡ WE_unpack(χ_AMC)     # UnitTest: Σ CG_λ CG CG / Ĵ0 · χ_J
  (B) χ_red←m ≡ χ_AMC          # χ_red = Σ_m CG_λ CG CG χ(m) / Ĵ0
                                 # (same /Ĵ as scalar reduce=true; 3 CGs for WE)

Phases (locked T×S convention — χ^β / ethS / factored_fII NOTES):
  AMC may print (−1)^{J+j}(−1)^{J'+j'} with half-integer j’s.
  Evaluate as **one** integer phase: combine all half-int j via
    (j2_a + j2_b + …)//2 + (integer J’s)
  Do **not** phase each half-int with j2//2 separately — that drops
  one (−1) per pair of fermions (j2//2+j2'//2 = (j2+j2')//2 − 1).

  χ^κ print: −(−1)^{J1+j_i}(−1)^{J2+J3+J4+j_c}
  → ph = −iphase(J1+J2+J3+J4 + (j2_i+j2_c)//2)

AMC tex: learn/amc_tts/factored_GIV/output/chi_kappa_analyze.tex
  (formula as printed; no overall-sign retune)

Usage:
  PYTHONPATH=build python3 -B run/test_chi_kappa_m_vs_amc.py [emax=1] [lambda=2]
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
# Tensor path for all λ (including 0): Ω is WE-reduced. λ=0 ≠ unreduced scalar.
if not Eta.IsReduced():
    Eta.MakeReduced()

orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
print(
    f"emax={emax} λ={lam} seed={seed} nsamp={nsamp}\n"
    f"  Ω reduced={Eta.IsReduced()} Jrank={Eta.GetJRank()} ||2b||={Eta.TwoBodyNorm():.4g}"
    f"  [tensor/WE path for all λ]\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  AMC χ: as printed (no overall-sign / Ĵ retune)"
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


def w_occ(a: int, c: int, d: int) -> float:
    return nbar(c) * nbar(d) * occ(a) + occ(c) * occ(d) * nbar(a)


def chi_m(i, mi, j, mj, b, mb, d, md) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for c in orbits:
                for mc in mrange(c):
                    w = w_occ(a, c, d)
                    if abs(w) < 1e-12:
                        continue
                    o = ut.GetMschemeMatrixElement_2b(
                        Eta, a, ma, i, mi, c, mc, d, md
                    )
                    if abs(o) < 1e-16:
                        continue
                    g = ut.GetMschemeMatrixElement_2b(
                        Gamma, j, mj, c, mc, b, mb, a, ma
                    )
                    sm += w * o * g
    return sm


def chi_amc(J0: int, J1: int, i: int, j: int, b: int, d: int) -> float:
    """AMC chi_kappa_analyze.tex with χ^β combined-integer phases."""
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(b), jo(d), J1)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for c in orbits:
            w = w_occ(a, c, d)
            if abs(w) < 1e-12:
                continue
            ja, jc = jo(a), jo(c)
            for J2 in range(
                abs(j2i(a) - j2i(i)) // 2, (j2i(a) + j2i(i)) // 2 + 1
            ):
                for J3 in range(
                    abs(j2i(c) - j2i(d)) // 2, (j2i(c) + j2i(d)) // 2 + 1
                ):
                    if not tri(J2, J3, lam):
                        continue
                    om = Eta.TwoBody.GetTBME_J(J2, J3, a, i, c, d)
                    if abs(om) < 1e-16:
                        continue
                    J4min = max(abs(j2i(j) - j2i(c)), abs(j2i(b) - j2i(a))) // 2
                    J4max = min(j2i(j) + j2i(c), j2i(b) + j2i(a)) // 2
                    for J4 in range(J4min, J4max + 1):
                        gam = Gamma.TwoBody.GetTBME_J(J4, J4, j, c, b, a)
                        if abs(gam) < 1e-16:
                            continue
                        J5min = max(
                            abs(j2i(a) - j2i(c)), abs(j2i(b) - j2i(j))
                        ) // 2
                        J5max = min(j2i(a) + j2i(c), j2i(b) + j2i(j)) // 2
                        for J5 in range(J5min, J5max + 1):
                            j0max = int(max(J2, J3, J0, J1, jo(i), ja, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                six = (
                                    SixJ(J3, lam, J2, jo(i), ja, j0)
                                    * SixJ(J1, lam, J0, jo(i), jo(j), j0)
                                    * SixJ(ja, jc, J5, jo(d), j0, J3)
                                    * SixJ(jo(b), jo(j), J5, jc, ja, J4)
                                    * SixJ(jo(d), jo(b), J1, jo(j), j0, J5)
                                )
                                if abs(six) < 1e-16:
                                    continue
                                hats = (
                                    hat(J2)
                                    * hat(J3)
                                    * (2 * J4 + 1)
                                    * (2 * J5 + 1)
                                    * (2 * j0 + 1)
                                )
                                # −(−1)^{J1+j_i}(−1)^{J2+J3+J4+j_c}
                                # → single integer phase (χ^β convention)
                                ph = -iphase(
                                    J1
                                    + J2
                                    + J3
                                    + J4
                                    + (j2i(i) + j2i(c)) // 2
                                )
                                tot += ph * hats * six * w * om * gam
    return hat(J0) * hat(J1) * tot


print("Caching χ_AMC as printed ...")
t0 = time.time()
chiJ: dict = {}
nfill = 0
for i in orbits:
    for j in orbits:
        for b in orbits:
            for d in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not (
                            tri(jo(i), jo(j), J0) and tri(jo(b), jo(d), J1)
                        ):
                            continue
                        v = chi_amc(J0, J1, i, j, b, d)
                        if abs(v) < 1e-16:
                            continue
                        chiJ[(J0, J1, i, j, b, d)] = v
                        nfill += 1
print(f"  nonzero={nfill}  ({time.time()-t0:.1f}s)")


def we_unpack(i, mi, j, mj, b, mb, d, md) -> float:
    """UnitTest WE unpack. No Pauli skip (χ not AS)."""
    if abs(mi + mj - mb - md) > 2 * lam:
        return 0.0
    M0 = (mi + mj) // 2
    M1 = (mb + md) // 2
    mu = M0 - M1
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0:
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        for J1 in range(abs(j2i(b) - j2i(d)) // 2, (j2i(b) + j2i(d)) // 2 + 1):
            if abs(M1) > J1:
                continue
            if not tri(J0, J1, lam):
                continue
            ccd = CG(jo(b), mb * 0.5, jo(d), md * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            if lam == 0:
                cj = 1.0 if (J0 == J1 and M0 == M1) else 0.0
            else:
                cj = CG(J1, M1, lam, mu, J0, M0)
            if abs(cj) < 1e-15:
                continue
            v = chiJ.get((J0, J1, i, j, b, d), 0.0)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


def chi_red_from_m(J0: int, J1: int, i: int, j: int, b: int, d: int) -> float:
    """REDUCED_UNREDUCED: χ_red = Σ CG_λ CG CG χ(m) / Ĵ0."""
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(b), jo(d), J1)):
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
            for mb in mrange(b):
                for md in mrange(d):
                    M1 = (mb + md) // 2
                    if abs(M1) > J1:
                        continue
                    mu = M0 - M1
                    if abs(mu) > lam:
                        continue
                    ccd = CG(jo(b), mb * 0.5, jo(d), md * 0.5, J1, M1)
                    if abs(ccd) < 1e-15:
                        continue
                    if lam == 0:
                        cj = 1.0 if (J0 == J1 and M0 == M1) else 0.0
                    else:
                        cj = CG(J1, M1, lam, mu, J0, M0)
                    if abs(cj) < 1e-15:
                        continue
                    cm = chi_m(i, mi, j, mj, b, mb, d, md)
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
                for b in orbits:
                    for mb in mrange(b):
                        for d in orbits:
                            for md in mrange(d):
                                if abs(mi + mj - mb - md) > 2 * lam:
                                    continue
                                if not any(
                                    abs(w_occ(a, c, d)) > 1e-12
                                    for a in orbits
                                    for c in orbits
                                ):
                                    continue
                                cands.append((i, mi, j, mj, b, mb, d, md))
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
        for b in orbits:
            for d in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not (
                            tri(jo(i), jo(j), J0) and tri(jo(b), jo(d), J1)
                        ):
                            continue
                        cr = chi_red_from_m(J0, J1, i, j, b, d)
                        ca = chiJ.get((J0, J1, i, j, b, d), 0.0)
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
    print("PASS — m χ^κ ≡ AMC (combined-integer phases, /Ĵ0 packaging)")
    sys.exit(0)
print(
    f"FAIL — max|m−WE|={max_A:.3e} max|red−A|={max_B:.3e}. "
    "If ratios ~ Ĵ^{±1}, packaging; if r=−1, check half-int phase combining."
)
sys.exit(1)
