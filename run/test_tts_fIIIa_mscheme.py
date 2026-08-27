#!/usr/bin/env python3
"""Step 1: m-scheme gold vs DIRECT (tts_fIIIa) for f^III_a.

Chain (factorized_code_analyze.tex §verify-chain)
-----------------------------------------------
* m-scheme unreduced equation is golden (AMC input f3a.txt / note §unfact).
* DIRECT = ReferenceImplementations.comm223_231_tts_fIIIa (AMC case2 RME).
* Compare only via UnitTest.GetMschemeMatrixElement_1b — no m-scheme patches.

m-scheme (AMC analyze):
  f_ij = Σ_abcde w(a,b,c,d) · (
      Ω_abcd Ω_idae Γ_cejb − Ω_abcd Ω_edaj Γ_cieb
  )
  w = n̄_a n̄_b n_c n_d − n_a n_b n̄_c n̄_d

Ω is tensor rank λ, f is scalar ⇒ [Ω×Ω]^(0) on the two Ω legs:
  CG(λ μ, λ −μ; 0 0) = (−1)^{λ−μ}/λ̂,  μ = (m_a+m_b−m_c−m_d)/2.

Usage:
  PYTHONPATH=build python3 run/test_tts_fIIIa_mscheme.py [emax=1] [lambda=2]
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

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}  "
    f"||Ω||={Eta.TwoBodyNorm():.4g}  ||Γ||={Gamma.TwoBodyNorm():.4g}  "
    f"Γ.reduced={Gamma.IsReduced()}"
)
print("compare: m (f3a analyze + [Ω×Ω]_0) vs GetMschemeMatrixElement_1b(DIRECT)\n")


def occ(a: int) -> float:
    return ms.GetOrbit(a).occ


def nbar(a: int) -> float:
    return 1.0 - occ(a)


def j2(a: int) -> int:
    return ms.GetOrbit(a).j2


def om2(*ix) -> float:
    return ut.GetMschemeMatrixElement_2b(Eta, *ix)


def g2(*ix) -> float:
    return ut.GetMschemeMatrixElement_2b(Gamma, *ix)


def m_range(o):
    return range(-j2(o), j2(o) + 1, 2)


def same_j(i, j) -> bool:
    return j2(i) == j2(j)


def cg0(mu: int) -> float:
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


def fIIIa_m(i, mi, j, mj) -> float:
    """Literal unfactored m-scheme f^III_a with scalar projection on Ω×Ω."""
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
                                    wK = cg0(mu)
                                    if abs(wK) < 1e-16:
                                        continue
                                    for e in orbits:
                                        for me in m_range(e):
                                            o_idae = om2(
                                                i, mi, d, md, a, ma, e, me
                                            )
                                            o_edaj = om2(
                                                e, me, d, md, a, ma, j, mj
                                            )
                                            if abs(o_idae) < 1e-16 and abs(
                                                o_edaj
                                            ) < 1e-16:
                                                continue
                                            g_cejb = g2(
                                                c, mc, e, me, j, mj, b, mb
                                            )
                                            g_cieb = g2(
                                                c, mc, i, mi, e, me, b, mb
                                            )
                                            sm += (
                                                wK
                                                * w
                                                * o_abcd
                                                * (
                                                    o_idae * g_cejb
                                                    - o_edaj * g_cieb
                                                )
                                            )
    return sm


print("=== f^III_a : m vs DIRECT (tts_fIIIa) ===")
t0 = time.time()
Z_dir = Operator(ms, 0, 0, 0, 2)
Z_dir.SetHermitian()
ReferenceImplementations.comm223_231_tts_fIIIa(Eta, Gamma, Z_dir)
print(
    f"  DIRECT: IsReduced={Z_dir.IsReduced()}  "
    f"||1b||={Z_dir.OneBodyNorm():.6e}  ({time.time()-t0:.1f}s)"
)

max_abs = 0.0
max_spread = 0.0
ratios = Counter()
n = 0
worst = None
t1 = time.time()
for i in orbits:
    for j in orbits:
        if not same_j(i, j):
            continue
        vals_m = []
        vals_j = []
        for mi in m_range(i):
            zm = fIIIa_m(i, mi, j, mi)
            zj = ut.GetMschemeMatrixElement_1b(Z_dir, i, mi, j, mi)
            vals_m.append(zm)
            vals_j.append(zj)
            err = abs(zm - zj)
            n += 1
            if err > max_abs:
                max_abs = err
                worst = (i, j, mi, zm, zj, err)
        if vals_m:
            max_spread = max(max_spread, max(vals_m) - min(vals_m))
            am = sum(vals_m) / len(vals_m)
            aj = sum(vals_j) / len(vals_j)
            if abs(aj) > 1e-12:
                ratios[round(am / aj, 6)] += 1
            # also print per-orbit diagonal snapshot
            if i == j:
                mi0 = j2(i)
                print(
                    f"  diag i={i} j={ms.GetOrbit(i).j2/2:.1f}  "
                    f"m={fIIIa_m(i, mi0, i, mi0):.6e}  "
                    f"DIRECT_m={ut.GetMschemeMatrixElement_1b(Z_dir, i, mi0, i, mi0):.6e}  "
                    f"spread_m={max(vals_m)-min(vals_m):.3e}"
                )

ok = max_abs < tol and max_spread < tol
print(
    f"\n  m vs DIRECT: n_m={n}  max|m−J_m|={max_abs:.3e}  "
    f"m_spread={max_spread:.3e}  m/J_m={ratios.most_common(4)}  "
    f"({time.time()-t1:.1f}s)  => {'PASS' if ok else 'FAIL'}"
)
if worst and not ok:
    i, j, mi, zm, zj, err = worst
    print(
        f"    worst ({i},{j}) mi={mi}: m={zm:.6e} J_m={zj:.6e} Δ={err:.3e}"
    )

print(
    """
Notes
  * m is scalar (m_spread~0) ⇒ [Ω×Ω]_0 CG is applied correctly.
  * DIRECT ≡ ethS slow (AMC W1/W2) — both J paths agree with each other.
  * m ≡ DIRECT locks the χ^γ + ladder oracle (use_TypeIIIa_slow).
  * Production CC Path B is checked separately: run/test_tts_fIIIa_pathB_cc.py.
  * Side check: unfactored m = −(χ^γ×Γ)_m with printed w_γ / eq:fIIIa
    (consistent with w_γ=−w_unfact; gold remains AMC f3a / note §unfact).
"""
)
sys.exit(0 if ok else 1)
