#!/usr/bin/env python3
"""ethS Gamma^II ≡ Path B −W_OG − V_ζ.

Γ^II = −W − V with
  W = (1−P_ij) χ^{ΩΓ}_ja Ω_iakl
  V = (1−P_kl) χ^ζ_ak Ω_ijal
See learn/amc_tts/factored_GII/NOTES.md.

Compare on channel-ordered kets (ibra≤iket) via GetTBME_J.

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GII_eths_pathB.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time

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
print(f"emax={emax} λ={lam} seed={seed}")


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


def occ(a):
    return ms.GetOrbit(a).occ


def nbar(a):
    return 1.0 - occ(a)


hat_lam_inv = 1.0 / math.sqrt(2 * max(lam, 0) + 1)
max_J = max(j2i(o) for o in orbits)


def chi_zeta_J(i, j) -> float:
    """AMC chi_zeta_analyze: χ^ζ_{ij} from Γ_ciab Ω_abcj."""
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)
                if abs(w) < 1e-12:
                    continue
                jc = jo(c)
                for J0 in range(0, max_J + 1):
                    if not (tri(jc, ji, J0) and tri(jo(a), jo(b), J0)):
                        continue
                    g = Gamma.TwoBody.GetTBME_J(J0, J0, c, i, a, b)
                    if abs(g) < 1e-16:
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam) or not tri(jc, jj, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, a, b, c, j)
                        if abs(o) < 1e-16:
                            continue
                        six = SixJ(lam, J1, J0, jc, ji, jj)
                        if abs(six) < 1e-16:
                            continue
                        ph = iphase((j2i(j) + j2i(c)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * g * o
    return 0.5 * sm


def chi_OG_J(i, j) -> float:
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)
                if abs(w) < 1e-12:
                    continue
                jc = jo(c)
                for J0 in range(0, max_J + 1):
                    if not tri(jc, ji, J0):
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not tri(jo(a), jo(b), J1):
                            continue
                        if not tri(jc, jj, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, c, i, a, b)
                        if abs(o) < 1e-16:
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J1, J1, a, b, c, j)
                        if abs(g) < 1e-16:
                            continue
                        six = SixJ(J0, J1, lam, jj, ji, jc)
                        if abs(six) < 1e-16:
                            continue
                        ph = iphase((j2i(j) + j2i(c)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * o * g
    return 0.5 * sm


print("Building χ^ζ and χ^{ΩΓ} ...")
t0 = time.time()
chi = {}
chi_og = {}
for i in orbits:
    for j in orbits:
        if tri(jo(i), jo(j), lam):
            chi[(i, j)] = chi_zeta_J(i, j)
            chi_og[(i, j)] = chi_OG_J(i, j)
print(f"  n={len(chi)} ({time.time()-t0:.2f}s)")


def pathB_Z_J(J, i, j, k, l) -> float:
    """−W_OG − V_ζ with (1−P) restored."""
    W = V = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            pref = 1.0 / hat(J) * hat(J2) * hat_lam_inv
            c = chi_og.get((j, a), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(j), jo(i))
                o = Eta.TwoBody.GetTBME_J(J2, J, i, a, k, l)
                if abs(six * o) > 1e-16:
                    W += iphase((j2i(i) + j2i(a)) // 2 + J2) * pref * six * c * o
            c = chi_og.get((i, a), 0.0)
            if abs(c) > 1e-16:
                six = SixJ(J, J2, lam, jo(a), jo(i), jo(j))
                o = Eta.TwoBody.GetTBME_J(J2, J, j, a, k, l)
                if abs(six * o) > 1e-16:
                    W += (
                        iphase(J + (j2i(i) + j2i(a)) // 2 + J2)
                        * pref
                        * six
                        * c
                        * o
                    )
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
                    V += (
                        iphase(J + (j2i(l) + j2i(a)) // 2 + J2)
                        * pref
                        * six
                        * c
                        * o
                    )
    return -W - V


# ethS GII only
cm = Commutator.FactorizedDoubleCommutator_eths
cm.SetUse_1b_Intermediates(True)
cm.SetUse_2b_Intermediates(False)
cm.SetUse_TypeGI_2b(False)
cm.SetUse_TypeGII_2b(True)
for name in (
    "SetUse_TypeGIIIa_2b",
    "SetUse_TypeGIIIb_2b",
    "SetUse_TypeGIIIc_2b",
    "SetUse_TypeGIVa_2b",
    "SetUse_TypeGIVb_2b",
    "SetUse_TypeGIVc_2b",
    "SetUse_TypeII_2b",
    "SetUse_TypeIII_2b",
    "SetUse_GT_TypeI_2b",
    "SetUse_GT_TypeIV_2b",
):
    fn = getattr(cm, name, None)
    if fn:
        fn(False)

Z = Operator(ms, 0, 0, 0, 2)
Z.SetHermitian()
if Z.IsReduced():
    Z.MakeNotReduced()

t0 = time.time()
cm.comm223_232(Eta, Gamma, Z)
print(f"ethS GII done ({time.time()-t0:.2f}s)  ‖Z‖={Z.TwoBodyNorm():.6e}")

# Channel-ordered kets only (matches Hermitian AddToTBME upper triangle)
max_abs = 0.0
n = 0
worst = None
for ch in range(ms.GetNumberTwoBodyChannels()):
    tbc = ms.GetTwoBodyChannel(ch)
    J = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        bra = tbc.GetKet(ib)
        i, j = bra.p, bra.q
        for ik in range(ib, nk):
            ket = tbc.GetKet(ik)
            k, l = ket.p, ket.q
            zb = pathB_Z_J(J, i, j, k, l)
            ze = Z.TwoBody.GetTBME_J(J, J, i, j, k, l)
            if abs(zb) < 1e-14 and abs(ze) < 1e-14:
                continue
            n += 1
            err = abs(zb - ze)
            if err > max_abs:
                max_abs = err
                worst = (i, j, k, l, J, zb, ze, err)

print(f"n={n}  max|PathB−ethS|={max_abs:.3e}")
if worst:
    i, j, k, l, J, zb, ze, err = worst
    print(
        f"  worst J={J} ({i},{j},{k},{l}): B={zb:.6e} ethS={ze:.6e} Δ={err:.3e}"
    )

ok = max_abs < tol and Z.TwoBodyNorm() > 1e-8
print("PASS — ethS Γ^II ≡ Path B (−W−V)" if ok else "FAIL")
sys.exit(0 if ok else 1)
