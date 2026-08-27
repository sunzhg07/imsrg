#!/usr/bin/env python3
"""Literal m-scheme (note §unfact) vs DIRECT / Path B — reduce/degeneracy only.

Rules
-----
* m-scheme is golden and **always unreduced** (physical ⟨m|…|m⟩).
* J-scheme may be reduced or not. Compare only after converting J → m via
  UnitTest.GetMschemeMatrixElement_* (WE / degeneracy). For scalar 1b unreduced
  storage that is just OneBody(i,j).
* Do **not** patch m-scheme with ad-hoc phases/partner swaps.

Tensor × tensor → scalar
------------------------
AMC declares Ω tensor and f scalar. The note writes Ω·Ω·…, but the *scalar*
LHS means only the K=0 piece of Ω^λ×Ω^λ. Numerically that is

  CG(λ μ, λ −μ; 0 0) = (−1)^{λ−μ}/λ̂

on the two Ω legs (same factor AMC inserts as λ̂^{-1} when reducing).
Without it, the product is m-dependent for j≥3/2 (K≠0 contamination) and
off by ∼λ̂ even after averaging.

Usage:
  PYTHONPATH=build python3 run/test_tts_f_mscheme.py [emax=1] [lambda=2] [I,II,IIIb]
"""

from __future__ import annotations

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
h_Omega = -1
h_Gamma = +1

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}  "
    f"||Ω||={Eta.TwoBodyNorm():.4g}  ||Γ||={Gamma.TwoBodyNorm():.4g}  "
    f"Γ.reduced={Gamma.IsReduced()}"
)
print("compare: m_literal vs GetMschemeMatrixElement_1b(Z_J)  [reduce handled by UnitTest]\n")


def occ(a: int) -> float:
    return ms.GetOrbit(a).occ


def nbar(a: int) -> float:
    return 1.0 - occ(a)


def om2(*ix) -> float:
    return ut.GetMschemeMatrixElement_2b(Eta, *ix)


def g2(*ix) -> float:
    return ut.GetMschemeMatrixElement_2b(Gamma, *ix)


def m_range(o):
    j2 = ms.GetOrbit(o).j2
    return range(-j2, j2 + 1, 2)


def same_j(i, j) -> bool:
    return ms.GetOrbit(i).j2 == ms.GetOrbit(j).j2


def cg_scalar_OmOm(mu: int) -> float:
    """K=0 projection of Ω^λ×Ω^λ (AMC: f scalar, Ω tensor)."""
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


def fI_m(i, mi, j, mj) -> float:
    # note (73–80) / input/f1.txt  + scalar projection on Ω×Ω
    if mi != mj:
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
                        for mc in m_range(c):
                            for d in orbits:
                                nd, nbd = occ(d), nbar(d)
                                w = nba * nbb * nc * nd - na * nb * nbc * nbd
                                if abs(w) < 1e-12:
                                    continue
                                for md in m_range(d):
                                    o_cdab = om2(c, mc, d, md, a, ma, b, mb)
                                    if abs(o_cdab) < 1e-16:
                                        continue
                                    mu = (mc + md - ma - mb) // 2
                                    wK = cg_scalar_OmOm(mu)
                                    if abs(wK) < 1e-16:
                                        continue
                                    for e in orbits:
                                        for me in m_range(e):
                                            o_abce = om2(
                                                a, ma, b, mb, c, mc, e, me
                                            )
                                            if abs(o_abce) < 1e-16:
                                                continue
                                            g_eidj = g2(
                                                e, me, i, mi, d, md, j, mj
                                            )
                                            g_diej = g2(
                                                d, md, i, mi, e, me, j, mj
                                            )
                                            sm += (
                                                wK
                                                * w
                                                * o_cdab
                                                * o_abce
                                                * (g_eidj + g_diej)
                                            )
    return 0.5 * sm


def fII_m(i, mi, j, mj) -> float:
    # AMC analyze f2a/f2b: Γ Ω Ω_eidj + h_Γ Γ Ω Ω_ejdi, with [Ω×Ω]_0 on the
    # two tensor Ω legs (scalar LHS). Partner is ejdi (not diej); see
    # factored_fII and factorized_code_analyze.tex §verify-chain.
    if mi != mj:
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
                        for mc in m_range(c):
                            for e in orbits:
                                ne, nbe = occ(e), nbar(e)
                                w = na * nb * nbc * nbe - nba * nbb * nc * ne
                                if abs(w) < 1e-12:
                                    continue
                                for me in m_range(e):
                                    for d in orbits:
                                        for md in m_range(d):
                                            g_cdab = g2(
                                                c, mc, d, md, a, ma, b, mb
                                            )
                                            if abs(g_cdab) < 1e-16:
                                                continue
                                            o_abce = om2(
                                                a, ma, b, mb, c, mc, e, me
                                            )
                                            if abs(o_abce) < 1e-16:
                                                continue
                                            mu = (ma + mb - mc - me) // 2
                                            wK = cg_scalar_OmOm(mu)
                                            if abs(wK) < 1e-16:
                                                continue
                                            o_eidj = om2(
                                                e, me, i, mi, d, md, j, mj
                                            )
                                            o_ejdi = om2(
                                                e, me, j, mj, d, md, i, mi
                                            )
                                            sm += wK * w * g_cdab * o_abce * (
                                                o_eidj + h_Gamma * o_ejdi
                                            )
    return 0.5 * sm


def fIIIb_m(i, mi, j, mj) -> float:
    # note (98–105) / input/f3b.txt  + scalar projection on Ω×Ω
    if mi != mj:
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
                        for mc in m_range(c):
                            for d in orbits:
                                nd, nbd = occ(d), nbar(d)
                                w = nba * nbb * nc * nd - na * nb * nbc * nbd
                                if abs(w) < 1e-12:
                                    continue
                                for md in m_range(d):
                                    o_abcd = om2(a, ma, b, mb, c, mc, d, md)
                                    if abs(o_abcd) < 1e-16:
                                        continue
                                    mu = (ma + mb - mc - md) // 2
                                    wK = cg_scalar_OmOm(mu)
                                    if abs(wK) < 1e-16:
                                        continue
                                    for e in orbits:
                                        for me in m_range(e):
                                            o_cdej = om2(
                                                c, mc, d, md, e, me, j, mj
                                            )
                                            g_eiab = g2(
                                                e, me, i, mi, a, ma, b, mb
                                            )
                                            o_eiab = om2(
                                                e, me, i, mi, a, ma, b, mb
                                            )
                                            g_cdej = g2(
                                                c, mc, d, md, e, me, j, mj
                                            )
                                            sm += wK * w * o_abcd * (
                                                o_cdej * g_eiab - o_eiab * g_cdej
                                            )
    return 0.25 * sm


def run_direct(diagram: str) -> Operator:
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    if diagram == "I":
        ReferenceImplementations.comm223_231_tts_fI(Eta, Gamma, Z)
    elif diagram == "II":
        ReferenceImplementations.comm223_231_tts_fII(Eta, Gamma, Z)
    else:
        ReferenceImplementations.comm223_231_tts_fIIIb(Eta, Gamma, Z)
    return Z


def run_pathB(diagram: str) -> Operator:
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    eth = Commutator.FactorizedDoubleCommutator_eths
    eth.SetUse_1b_Intermediates(diagram in ("I", "II"))
    eth.SetUse_2b_Intermediates(diagram == "IIIb")
    eth.SetUse_TypeI_1b(diagram == "I")
    eth.SetUse_TypeII_1b(diagram == "II")
    eth.SetUse_TypeIII_1b(diagram == "IIIb")
    eth.SetUse_TypeIIIa_1b(False)
    eth.comm223_231_st(Eta, Gamma, Z)
    return Z


def compare_m_vs_J(tag, fn, Zj) -> bool:
    """Compare literal m to J unpacked with GetMschemeMatrixElement_1b."""
    print(f"  Z_J.IsReduced={Zj.IsReduced()}  ||1b||={Zj.OneBodyNorm():.6e}")
    max_abs = 0.0
    max_spread = 0.0
    ratios = Counter()
    n = 0
    worst = None
    t0 = time.time()
    for i in orbits:
        oi = ms.GetOrbit(i)
        for j in orbits:
            if not same_j(i, j):
                continue
            j2 = oi.j2
            vals_m = []
            vals_j = []
            for mi in range(-j2, j2 + 1, 2):
                if abs(mi) > ms.GetOrbit(j).j2:
                    continue
                zm = fn(i, mi, j, mi)
                zj = ut.GetMschemeMatrixElement_1b(Zj, i, mi, j, mi)
                vals_m.append(zm)
                vals_j.append(zj)
                err = abs(zm - zj)
                n += 1
                if err > max_abs:
                    max_abs = err
                    worst = (i, j, mi, zm, zj, err)
            if vals_m:
                max_spread = max(max_spread, max(vals_m) - min(vals_m))
                am, aj = sum(vals_m) / len(vals_m), sum(vals_j) / len(vals_j)
                if abs(aj) > 1e-12:
                    ratios[round(am / aj, 6)] += 1
    ok = max_abs < tol
    print(
        f"  {tag}: n_m={n}  max|m−J_m|={max_abs:.3e}  "
        f"m_spread={max_spread:.3e}  m/J_m={ratios.most_common(4)}  "
        f"({time.time()-t0:.1f}s)  => {'PASS' if ok else 'FAIL'}"
    )
    if worst and not ok:
        i, j, mi, zm, zj, err = worst
        print(
            f"    worst ({i},{j}) mi={mi}: m={zm:.6e} J_m={zj:.6e} Δ={err:.3e}"
        )
    return ok


all_ok = True
diags = sys.argv[3].split(",") if len(sys.argv) > 3 else ["I", "II", "IIIb"]
fns = {"I": fI_m, "II": fII_m, "IIIb": fIIIb_m}

for diag in diags:
    print(f"=== f^{diag} ===")
    Zd = run_direct(diag)
    Zb = run_pathB(diag)
    ok_d = compare_m_vs_J("m vs DIRECT", fns[diag], Zd)
    ok_b = compare_m_vs_J("m vs PathB ", fns[diag], Zb)
    # DIRECT vs PathB at J level (unreduced 1b)
    max_jj = 0.0
    for i in orbits:
        for j in orbits:
            if same_j(i, j):
                max_jj = max(
                    max_jj, abs(Zd.GetOneBody(i, j) - Zb.GetOneBody(i, j))
                )
    print(f"  DIRECT vs PathB (OneBody): max|Δ|={max_jj:.3e}")
    all_ok = all_ok and ok_d and ok_b and max_jj < tol
    print()

print("ALL PASS" if all_ok else "SOME FAILED — see notes below")
print(
    """
Interpretation
  * J reduce/degeneracy: GetMschemeMatrixElement_1b(Z) ≡ Z_unred for scalar 1b.
  * Tensor×tensor→scalar: require [Ω×Ω]_0 CG (AMC scalar LHS).
  * Chain: m golden → AMC DIRECT/PathB → code, no tune
    (factorized_code_analyze.tex §verify-chain).
  * f^I / f^II / f^III_b: expect m ≡ DIRECT ≡ Path B.
"""
)
sys.exit(0 if all_ok else 1)
