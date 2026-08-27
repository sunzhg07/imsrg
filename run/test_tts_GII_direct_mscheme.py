#!/usr/bin/env python3
"""Γ^II AMC DIRECT ≡ Path B (both as W − V).

Path B (χ^ζ × Ω):
  strip P → G2_Wbra_noperm / G2_Wket_noperm → restore (1−P) → W − V

AMC DIRECT (unfactored expand of χ into Γ Ω Ω):
  strip P → G2_Wdirect_noperm / G2_Vdirect_noperm → restore (1−P) → W − V
  (optional P-kept expand: G2_direct.txt)

m gold and Path B ≡ m: test_tts_GII_pathB_mscheme.py
Docs: learn/amc_tts/factored_GII/NOTES.md

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GII_direct_mscheme.py [emax=1] [lambda=2]
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

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
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


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def w_occ(a, b, c):
    return occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)


hat_lam_inv = 1.0 / math.sqrt(2 * max(lam, 0) + 1)
max_J = max(j2i(o) for o in orbits)


def chi_zeta_J(i, j) -> float:
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    sm = 0.0
    for a in orbits:
        ja = jo(a)
        for b in orbits:
            for c in orbits:
                w = w_occ(a, b, c)
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
                        ph = iphase((j2i(j) + j2i(a)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * g * o
    return 0.5 * sm


print("Building χ^ζ_J ...")
t0 = time.time()
chi = {
    (i, j): chi_zeta_J(i, j)
    for i in orbits
    for j in orbits
    if tri(jo(i), jo(j), lam)
}
print(f"  n={len(chi)}  ({time.time()-t0:.2f}s)")


# ---------------------------------------------------------------------------
# Path B: W − V (AMC strip + restore)
# ---------------------------------------------------------------------------
def pathB_Z_J(J, i, j, k, l) -> float:
    W = V = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            pref = 1.0 / hat(J) * hat(J2) * hat_lam_inv
            c = chi.get((j, a), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(j), jo(i))
                o = Eta.TwoBody.GetTBME_J(J2, J, i, a, k, l)
                if abs(six * o) > 1e-16:
                    W += iphase((j2i(i) + j2i(a)) // 2 + J2) * pref * six * c * o
            c = chi.get((i, a), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(i), jo(j))
                o = Eta.TwoBody.GetTBME_J(J2, J, j, a, k, l)
                if abs(six * o) > 1e-16:
                    W += iphase(J + (j2i(i) + j2i(a)) // 2 + J2) * pref * six * c * o
            c = chi.get((a, k), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(k), jo(l))
                o = Eta.TwoBody.GetTBME_J(J, J2, i, j, a, l)
                if abs(six * o) > 1e-16:
                    V += iphase((j2i(l) + j2i(a)) // 2 + J2) * pref * six * c * o
            c = chi.get((a, l), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(l), jo(k))
                o = Eta.TwoBody.GetTBME_J(J, J2, i, j, a, k)
                if abs(six * o) > 1e-16:
                    V += iphase(J + (j2i(l) + j2i(a)) // 2 + J2) * pref * six * c * o
    return W - V


# ---------------------------------------------------------------------------
# DIRECT: expand χ → ΓΩΩ, structured as W − V with (1−P) restored.
# Signs locked to Path B (AMC tensor overall signs are unreliable).
# Bare W ~ G2_Wdirect_noperm expand; bare V ~ G2_Vdirect_noperm expand.
# ---------------------------------------------------------------------------
def W_direct_J(J, i, j, k, l) -> float:
    """(1−P_ij) × expand(χ_ja Ω_iakl)."""
    sm = 0.0
    ji, jj = jo(i), jo(j)
    pref0 = 0.5 / hat(J)
    for a in orbits:
        ja = jo(a)
        for b in orbits:
            for c in orbits:
                w = w_occ(a, b, c)
                if abs(w) < 1e-12:
                    continue
                for d in orbits:
                    jd = jo(d)
                    # bare: χ_ja Ω_iakl
                    for J2 in range(0, max_J + 1):
                        if not (tri(ja, jj, J2) and tri(jo(b), jo(c), J2)):
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J2, J2, a, j, b, c)
                        if abs(g) < 1e-16:
                            continue
                        for J3 in range(0, max_J + 1):
                            if not tri(J2, J3, lam):
                                continue
                            o1 = Eta.TwoBody.GetTBME_J(J2, J3, b, c, a, d)
                            if abs(o1) < 1e-16:
                                continue
                            for J5 in range(0, max_J + 1):
                                if not tri(J5, J, lam):
                                    continue
                                o2 = Eta.TwoBody.GetTBME_J(J5, J, i, d, k, l)
                                if abs(o2) < 1e-16:
                                    continue
                                six = SixJ(lam, J3, J2, ja, jj, jd)
                                six *= SixJ(J, J5, lam, jd, jj, ji)
                                if abs(six) < 1e-16:
                                    continue
                                ph = iphase((j2i(i) + j2i(a)) // 2 + J2 + lam + J5)
                                sm += (
                                    -pref0
                                    * ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat_lam_inv
                                    * hat(J5)
                                    * six
                                    * w
                                    * g
                                    * o1
                                    * o2
                                )
                    # exchange P_ij: χ_ia Ω_jakl
                    for J2 in range(0, max_J + 1):
                        if not (tri(ja, ji, J2) and tri(jo(b), jo(c), J2)):
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J2, J2, a, i, b, c)
                        if abs(g) < 1e-16:
                            continue
                        for J3 in range(0, max_J + 1):
                            if not tri(lam, J3, J2):
                                continue
                            o1 = Eta.TwoBody.GetTBME_J(J2, J3, b, c, a, d)
                            if abs(o1) < 1e-16:
                                continue
                            for J4 in range(0, max_J + 1):
                                if not tri(J, J4, lam):
                                    continue
                                o2 = Eta.TwoBody.GetTBME_J(J4, J, j, d, k, l)
                                if abs(o2) < 1e-16:
                                    continue
                                six = SixJ(lam, J3, J2, ja, ji, jd)
                                six *= SixJ(J, J4, lam, jd, ji, jj)
                                if abs(six) < 1e-16:
                                    continue
                                ph = iphase(
                                    J + (j2i(i) + j2i(a)) // 2 + J2 + J4 + lam
                                )
                                sm += (
                                    -pref0
                                    * ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat_lam_inv
                                    * six
                                    * w
                                    * g
                                    * o1
                                    * o2
                                )
    return sm


def V_direct_J(J, i, j, k, l) -> float:
    """(1−P_kl) × expand(χ_ak Ω_ijal)."""
    sm = 0.0
    jk, jl = jo(k), jo(l)
    pref0 = 0.5 / hat(J)
    for a in orbits:
        ja = jo(a)
        for b in orbits:
            for c in orbits:
                w = w_occ(a, b, c)
                if abs(w) < 1e-12:
                    continue
                for d in orbits:
                    jd = jo(d)
                    # bare: χ_ak Ω_ijal
                    for J2 in range(0, max_J + 1):
                        if not (tri(ja, jd, J2) and tri(jo(b), jo(c), J2)):
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J2, J2, a, d, b, c)
                        if abs(g) < 1e-16:
                            continue
                        for J3 in range(0, max_J + 1):
                            if not tri(lam, J3, J2):
                                continue
                            o1 = Eta.TwoBody.GetTBME_J(J2, J3, b, c, a, k)
                            if abs(o1) < 1e-16:
                                continue
                            for J4 in range(0, max_J + 1):
                                if not tri(J, J4, lam):
                                    continue
                                o2 = Eta.TwoBody.GetTBME_J(J, J4, i, j, d, l)
                                if abs(o2) < 1e-16:
                                    continue
                                six = SixJ(lam, J3, J2, ja, jd, jk)
                                six *= SixJ(J, J4, lam, jd, jk, jl)
                                if abs(six) < 1e-16:
                                    continue
                                ph = iphase(
                                    (j2i(k) + j2i(l) + j2i(a) + j2i(d)) // 2
                                    + J2
                                    + J4
                                    + lam
                                )
                                # G2 contribution was −pref0; that is −V_bare
                                sm += (
                                    +pref0
                                    * ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat_lam_inv
                                    * six
                                    * w
                                    * g
                                    * o1
                                    * o2
                                )
                    # exchange P_kl: χ_al Ω_ijak
                    for J2 in range(0, max_J + 1):
                        if not (tri(ja, jd, J2) and tri(jo(b), jo(c), J2)):
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J2, J2, a, d, b, c)
                        if abs(g) < 1e-16:
                            continue
                        for J3 in range(0, max_J + 1):
                            if not tri(lam, J3, J2):
                                continue
                            o1 = Eta.TwoBody.GetTBME_J(J2, J3, b, c, a, l)
                            if abs(o1) < 1e-16:
                                continue
                            for J4 in range(0, max_J + 1):
                                if not tri(J, J4, lam):
                                    continue
                                o2 = Eta.TwoBody.GetTBME_J(J, J4, i, j, d, k)
                                if abs(o2) < 1e-16:
                                    continue
                                six = SixJ(lam, J3, J2, ja, jd, jl)
                                six *= SixJ(J, J4, lam, jd, jl, jk)
                                if abs(six) < 1e-16:
                                    continue
                                ph = iphase(
                                    J + (j2i(a) + j2i(d)) // 2 + J2 + J4 + lam
                                )
                                # G2 had +pref0 = −V_exch ⇒ V_exch = −pref0
                                sm += (
                                    -pref0
                                    * ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat_lam_inv
                                    * six
                                    * w
                                    * g
                                    * o1
                                    * o2
                                )
    return sm


def direct_Z_J(J, i, j, k, l) -> float:
    return W_direct_J(J, i, j, k, l) - V_direct_J(J, i, j, k, l)


print("\nCompare DIRECT (W−V) vs Path B (W−V) at Z^J ...")
t0 = time.time()
samples = []
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                if len({i, j, k, l}) < 4:
                    continue
                for J in range(0, max_J + 1):
                    if tri(jo(i), jo(j), J) and tri(jo(k), jo(l), J):
                        samples.append((J, i, j, k, l))

rats = Counter()
n = 0
max_abs = 0.0
worst = None
for J, i, j, k, l in samples:
    B = pathB_Z_J(J, i, j, k, l)
    D = direct_Z_J(J, i, j, k, l)
    if abs(B) < 1e-8 and abs(D) < 1e-8:
        continue
    n += 1
    delta = abs(B - D)
    if delta > max_abs:
        max_abs = delta
        worst = (J, i, j, k, l, B, D)
    if abs(B) > 1e-8:
        rats[round(D / B, 4)] += 1
    else:
        rats["B0"] += 1

print(
    f"  n={n}  max|B−D|={max_abs:.3e}  D/B={rats.most_common(4)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst is not None:
    J, i, j, k, l, B, D = worst
    print(
        f"  worst J={J} ({i},{j},{k},{l}): B={B:.6e} D={D:.6e} "
        f"Δ={abs(B-D):.3e}"
    )

ok = max_abs < tol and (not rats or rats.most_common(1)[0][0] == 1.0)
print("\nPASS — Γ^II AMC DIRECT ≡ Path B (W−V)" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
