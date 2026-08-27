#!/usr/bin/env python3
"""χ^β (tensor 1b, rank λ): m-scheme vs AMC/ethS J-scheme.

m-scheme (golden, always unreduced) — AMC input/chi_beta.txt:
  χ_de(md,me) = 1/2 Σ_{abc} (n_a n_b n̄_c n̄_e − n̄_a n̄_b n_c n_e)
                · Γ_cdab(m) · Ω_abce(m)
  Γ scalar × Ω^λ → χ^λ (no extra [Ω×Ω]_0).

J-scheme (AMC chi_beta.tex / ethS Chi_beta), **reduced** RME:
  χ_de^λ = 1/2 Σ (−1)^{je+jc+λ+J0} Ĵ0 Ĵ1 {λ J1 J0; jc jd je}
           · w · Γ^{J0}_{cdab} · Ω^{J0 J1 λ}_{abce}

  Phase: use integer form (−1)^{(j2_e+j2_c)/2 + λ + J0} (ethS).
  Written AMC (−1)^{je+λ}(−1)^{J0+jc} is the same only when evaluated
  as a single integer phase — do not phase()-half-integers separately.

Compare: m_lit vs GetMschemeMatrixElement_1b(Chi_J)  [WE handles reduce].

Γ: unreduced (ethS MakeNotReduced if needed). Ω: reduced tensor.

Usage:
  PYTHONPATH=build python3 run/test_chi_beta_mscheme.py [emax=1] [lambda=2]
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


def tri(a, b, c):
    return abs(a - b) <= c <= a + b


def chi_beta_m(d, md, e, me) -> float:
    ne, nbe = occ(e), nbar(e)
    sm = 0.0
    for a in orbits:
        na, nba = occ(a), nbar(a)
        for ma in m_range(a):
            for b in orbits:
                nb, nbb = occ(b), nbar(b)
                for mb in m_range(b):
                    for c in orbits:
                        nc, nbc = occ(c), nbar(c)
                        w = na * nb * nbc * nbe - nba * nbb * nc * ne
                        if abs(w) < 1e-12:
                            continue
                        for mc in m_range(c):
                            g = ut.GetMschemeMatrixElement_2b(
                                Gamma, c, mc, d, md, a, ma, b, mb
                            )
                            if abs(g) < 1e-16:
                                continue
                            o = ut.GetMschemeMatrixElement_2b(
                                Eta, a, ma, b, mb, c, mc, e, me
                            )
                            sm += w * g * o
    return 0.5 * sm


def chi_beta_J_rme(d, e) -> float:
    jd, je = jo(d), jo(e)
    if not tri(jd, je, lam):
        return 0.0
    ne, nbe = occ(e), nbar(e)
    max_J = max(j2i(o) for o in orbits)
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = (
                    occ(a) * occ(b) * nbar(c) * nbe
                    - nbar(a) * nbar(b) * occ(c) * ne
                )
                if abs(w) < 1e-12:
                    continue
                jc = jo(c)
                for J0 in range(0, max_J + 1):
                    if not (tri(jc, jd, J0) and tri(jo(a), jo(b), J0)):
                        continue
                    g = Gamma.TwoBody.GetTBME_J(J0, J0, c, d, a, b)
                    if abs(g) < 1e-16:
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam) or not tri(jc, je, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, a, b, c, e)
                        if abs(o) < 1e-16:
                            continue
                        six = SixJ(lam, J1, J0, jc, jd, je)
                        if abs(six) < 1e-16:
                            continue
                        # ethS: phase((j2_e+j2_c)/2 + λ + J0)
                        ph = iphase((j2i(e) + j2i(c)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * g * o
    return 0.5 * sm


print("\nBuilding reduced χ^β_J ...")
t0 = time.time()
Chi = Operator(ms, lam, 0, 0, 1)
Chi.SetNonHermitian()
nch = 0
for d in orbits:
    for e in orbits:
        if not tri(jo(d), jo(e), lam):
            continue
        Chi.SetOneBody(d, e, chi_beta_J_rme(d, e))
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
for d in orbits:
    for e in orbits:
        if not tri(jo(d), jo(e), lam):
            continue
        for md in m_range(d):
            for me in m_range(e):
                if abs(md - me) > 2 * lam:
                    continue
                zm = chi_beta_m(d, md, e, me)
                zj = ut.GetMschemeMatrixElement_1b(Chi, d, md, e, me)
                if abs(zm) < 1e-14 and abs(zj) < 1e-14:
                    continue
                n += 1
                err = abs(zm - zj)
                if err > max_abs:
                    max_abs = err
                    worst = (d, e, md, me, zm, zj, err)
                if abs(zj) > 1e-12:
                    rats[round(zm / zj, 6)] += 1

print(
    f"  n={n}  max|m−J_m|={max_abs:.3e}  m/J_m={rats.most_common(3)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst:
    d, e, md, me, zm, zj, err = worst
    print(
        f"  worst de=({d},{e}) m=({md},{me}): "
        f"m={zm:.6e} J_m={zj:.6e} Δ={err:.3e}"
    )

ok = max_abs < tol
print("\nPASS — χ^β m ≡ J (reduce via WE unpack)" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
