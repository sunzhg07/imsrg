#!/usr/bin/env python3
"""Γ^I three-way: AMC χ^ε×Γ (J) ≡ ethS slow ≡ ethS DGEMM; χ^ε m≡AMC separate.

Locks
-----
* χ^ε m ≡ AMC: run/test_chi_epsilon_mscheme.py
* ethS ≡ tts_GI: run/test_tts_GI.py

This test: build Z from AMC χ^ε (unreduced ĵ^{-2}) × Γ via G1_from_chi
four-term fold (same as ethS write path), compare TBME to ethS.

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GI_mscheme.py [emax=1] [lambda=2]
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


def j2i(a):
    return ms.GetOrbit(a).j2


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def hat_lam_inv():
    return 1.0 / math.sqrt(2 * lam + 1)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def ket_phase(J, p, q) -> float:
    """Ket::Phase(J) = (−1)^{J + (jp+jq)/2 + 1}  (|pqJ> = Phase |qpJ>)."""
    return iphase(J + (j2i(p) + j2i(q)) // 2 + 1)


def chi_eps_J(i, j) -> float:
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
                        sm += w * iphase(J0 + J1 + lam) * hat_lam_inv() * o1 * o2
    return 0.5 * sm / (j2i(i) + 1.0)


# AMC χ^ε
chi = {}
t0 = time.time()
for i in orbits:
    for j in orbits:
        if j2i(i) == j2i(j):
            chi[(i, j)] = chi_eps_J(i, j)
print(f"χ^ε built  ({time.time()-t0:.2f}s)  n={len(chi)}")


def pathB_me(J, p, q, r, s) -> float:
    """Un-normalized TBME before √2 identical-particle factors."""
    phase_pg = ket_phase(J, p, q)
    phase_qh = ket_phase(J, r, s)
    z = 0.0
    for d in orbits:
        if j2i(d) == j2i(p):
            z += chi.get((p, d), 0.0) * Gamma.TwoBody.GetTBME_J(J, J, d, q, r, s)
        if j2i(d) == j2i(q):
            z += (
                phase_pg
                * chi.get((q, d), 0.0)
                * Gamma.TwoBody.GetTBME_J(J, J, d, p, r, s)
            )
        if j2i(d) == j2i(s):
            z += chi.get((d, s), 0.0) * Gamma.TwoBody.GetTBME_J(J, J, p, q, r, d)
        if j2i(d) == j2i(r):
            z += (
                phase_qh
                * chi.get((d, r), 0.0)
                * Gamma.TwoBody.GetTBME_J(J, J, p, q, d, s)
            )
    if p == q:
        z /= math.sqrt(2.0)
    if r == s:
        z /= math.sqrt(2.0)
    return z


cm = Commutator.FactorizedDoubleCommutator_eths
cm.SetUse_1b_Intermediates(True)
cm.SetUse_2b_Intermediates(False)
cm.SetUse_TypeGI_2b(True)
cm.SetUse_TypeGII_2b(False)
for name in (
    "SetUse_TypeGIIIa_2b",
    "SetUse_TypeGIIIb_2b",
    "SetUse_TypeGIIIc_2b",
    "SetUse_TypeGIVa_2b",
    "SetUse_TypeGIVb_2b",
    "SetUse_TypeGIVc_2b",
):
    getattr(cm, name)(False)

# Gamma^I has a single production path (no slow/fast fork), so the only
# comparison left is AMC Path B vs the ethS implementation.
Z_fast = Operator(ms, 0, 0, 0, 2)
Z_fast.SetHermitian()
t0 = time.time()
cm.comm223_232(Eta, Gamma, Z_fast)
print(f"ethS DGEMM ||2b||={Z_fast.TwoBodyNorm():.6e}  ({time.time()-t0:.2f}s)")

print("\nCompare PathB(AMC χ×Γ) vs ethS GetTBME_J ...")
t0 = time.time()
max_abs = 0.0
n = 0
n_match = 0
worst = None
nch = ms.GetNumberTwoBodyChannels()
for ch in range(nch):
    tbc = ms.GetTwoBodyChannel(ch)
    J = tbc.J
    nk = tbc.GetNumberKets()
    for ibra in range(nk):
        bra = tbc.GetKet(ibra)
        p, q = bra.p, bra.q
        for iket in range(ibra, nk):
            ket = tbc.GetKet(iket)
            r, s = ket.p, ket.q
            zb = pathB_me(J, p, q, r, s)  # normalized (after √2), like AddToTBME
            zj = Z_fast.TwoBody.GetTBME_J_norm(J, J, p, q, r, s)
            if abs(zb) < 1e-14 and abs(zj) < 1e-14:
                continue
            n += 1
            err = abs(zb - zj)
            if err < tol:
                n_match += 1
            if err > max_abs:
                max_abs = err
                worst = (p, q, r, s, J, zb, zj, err)

print(f"  n={n}  match={n_match}  max|PathB−ethS|={max_abs:.3e}  ({time.time()-t0:.2f}s)")
if worst:
    p, q, r, s, J, zb, zj, err = worst
    print(
        f"  worst pqrs=({p},{q},{r},{s}) J={J}: "
        f"PathB={zb:.6e} ethS={zj:.6e} Δ={err:.3e}"
    )

ok = max_abs < tol
print(
    "\nPASS — AMC PathB ≡ ethS DGEMM (Γ^I)"
    if ok
    else "\nFAIL"
)
sys.exit(0 if ok else 1)
