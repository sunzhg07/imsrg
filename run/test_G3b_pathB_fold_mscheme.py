#!/usr/bin/env python3
"""Γ^{III_b} full Z (IIb+IId) — any-λ analyze-fold gold.

  m_4index  ≡  fold(χ^η_m)     ← this bench (any λ)
  χ^η_m_red ≡ χ^η_AMC ≡ χ^η_PathB   (test_chi_eta_mscheme.py, JT)

Why not Z[χ_PathB_WE]≡Z[χ_m] like IV_b?
  χ^η is scalar but **not antisymmetrized**. JT-reduced χ loses m-scheme
  content that the arxiv fold still uses; WE unpack ≠ χ_m, and
  Z[(1−P)^2 W[χ_AS]] ≠ Z[W[χ_m]]. IV_b’s χ^ι is a proper tensor (WE OK).

Pipeline (mirror IV_b fold packaging; χ from locked m / Path B JT):
  1) χ^η_m = Σ_ab w_η CG[λ] Ω_iabl Ω_bjka   (locked)
  2) Analyze fold (arxiv): W = −Σ_ab (χ_bkai Γ_jbla + χ_lajb Γ_aibk)
  3) Z = (1−P_ij)(1−P_kl) W
  4) lock Z_fold ≡ Z_4index (IIb+IId with CG for Ω×Ω)

Factorized Fac Pandya χ̄ remains the λ=0 code-path dual
(test_G3b_factorized_mscheme.py). Any-λ production dual is Path B pack
(test_G3b_pathB_pack_mscheme.py).

Usage:
  PYTHONPATH=build python3 -B run/test_G3b_pathB_fold_mscheme.py [emax=1] [lambda=2] [nsamp=12]
"""

from __future__ import annotations

import random
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 12
tol = 1e-5
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()
if not Omega.IsReduced():
    Omega.MakeReduced()

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed} nsamp={nsamp}\n"
    f"  WRAP-UP: fold(χ^η_m) ≡ m_4index  "
    f"(χ Path B locked at JT; Factorized RC = λ=0 dual)"
)


def j2i(a: int) -> int:
    return ms.GetOrbit(a).j2


def occ(a: int) -> float:
    return ms.GetOrbit(a).occ


def nbar(a: int) -> float:
    return 1.0 - occ(a)


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def w_eta(a: int, b: int, k: int) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def om(*x) -> float:
    return ut.GetMschemeMatrixElement_2b(Omega, *x)


def gm(*x) -> float:
    return ut.GetMschemeMatrixElement_2b(Gamma, *x)


def cg_tt_scalar(m_bra1, m_bra2, m_ket1, m_ket2) -> float:
    """CG(λ μ, λ −μ; 0 0) for T×T→S; μ from first Ω’s ΔM."""
    mu = 0.5 * (m_bra1 + m_bra2 - m_ket1 - m_ket2)
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


# ---------------------------------------------------------------------------
# χ^η_m (locked definition)
# ---------------------------------------------------------------------------
_chi: dict = {}


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    key = (i, mi, j, mj, k, mk, l, ml)
    if key in _chi:
        return _chi[key]
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    w = w_eta(a, b, k)
                    if abs(w) < 1e-12:
                        continue
                    o1 = om(i, mi, a, ma, b, mb, l, ml)
                    if abs(o1) < 1e-16:
                        continue
                    cg = cg_tt_scalar(mi, ma, mb, ml)
                    if abs(cg) < 1e-16:
                        continue
                    o2 = om(b, mb, j, mj, k, mk, a, ma)
                    sm += cg * w * o1 * o2
    _chi[key] = sm
    return sm


# ---------------------------------------------------------------------------
# Analyze fold (arxiv Γ^{III_b})
# ---------------------------------------------------------------------------
def Wm_fold(i, mi, j, mj, k, mk, l, ml) -> float:
    """W = −Σ_ab (χ_bkai Γ_jbla + χ_lajb Γ_aibk)."""
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    c1 = chi_m(b, mb, k, mk, a, ma, i, mi)
                    if abs(c1) > 1e-16:
                        g1 = gm(j, mj, b, mb, l, ml, a, ma)
                        if abs(g1) > 1e-16:
                            sm -= c1 * g1
                    c2 = chi_m(l, ml, a, ma, j, mj, b, mb)
                    if abs(c2) > 1e-16:
                        g2 = gm(a, ma, i, mi, b, mb, k, mk)
                        if abs(g2) > 1e-16:
                            sm -= c2 * g2
    return sm


def Zm_fold(i, mi, j, mj, k, mk, l, ml) -> float:
    w = Wm_fold(i, mi, j, mj, k, mk, l, ml)
    w -= Wm_fold(j, mj, i, mi, k, mk, l, ml)
    w -= Wm_fold(i, mi, j, mj, l, ml, k, mk)
    w += Wm_fold(j, mj, i, mi, l, ml, k, mk)
    return w


# ---------------------------------------------------------------------------
# m gold: IIb+IId 4-index (CG for Ω×Ω when λ≠0)
# ---------------------------------------------------------------------------
def Wm_4index(i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    for c in orbits:
                        for mc in mrange(c):
                            for d in orbits:
                                for md in mrange(d):
                                    wb = (
                                        nbar(b) * occ(c) * occ(d)
                                        + nbar(c) * nbar(d) * occ(b)
                                    )
                                    if abs(wb) > 1e-12:
                                        o1 = om(d, md, c, mc, b, mb, k, mk)
                                        if abs(o1) > 1e-16:
                                            cg = cg_tt_scalar(md, mc, mb, mk)
                                            if abs(cg) > 1e-16:
                                                o2 = om(
                                                    b, mb, i, mi, a, ma, c, mc
                                                )
                                                if abs(o2) > 1e-16:
                                                    g = gm(
                                                        j,
                                                        mj,
                                                        a,
                                                        ma,
                                                        l,
                                                        ml,
                                                        d,
                                                        md,
                                                    )
                                                    sm -= wb * cg * o1 * o2 * g
                                    wd = (
                                        nbar(c) * occ(b) * occ(d)
                                        + nbar(b) * nbar(d) * occ(c)
                                    )
                                    if abs(wd) > 1e-12:
                                        o1 = om(j, mj, c, mc, b, mb, d, md)
                                        if abs(o1) > 1e-16:
                                            cg = cg_tt_scalar(mj, mc, mb, md)
                                            if abs(cg) > 1e-16:
                                                o2 = om(
                                                    b, mb, a, ma, l, ml, c, mc
                                                )
                                                if abs(o2) > 1e-16:
                                                    g = gm(
                                                        d,
                                                        md,
                                                        i,
                                                        mi,
                                                        a,
                                                        ma,
                                                        k,
                                                        mk,
                                                    )
                                                    sm -= wd * cg * o1 * o2 * g
    return sm


def Zm_4index(i, mi, j, mj, k, mk, l, ml) -> float:
    w = Wm_4index(i, mi, j, mj, k, mk, l, ml)
    w -= Wm_4index(j, mj, i, mi, k, mk, l, ml)
    w -= Wm_4index(i, mi, j, mj, l, ml, k, mk)
    w += Wm_4index(j, mj, i, mi, l, ml, k, mk)
    return w


# ---------------------------------------------------------------------------
# Sample
# ---------------------------------------------------------------------------
random.seed(seed)
cands = []
for _ in range(4000):
    i, j, k, l = (random.choice(orbits) for _ in range(4))
    mi = random.choice(list(mrange(i)))
    mj = random.choice(list(mrange(j)))
    M = mi + mj
    mks = [mk for mk in mrange(k) if (M - mk) in mrange(l)]
    if not mks:
        continue
    mk = random.choice(mks)
    ml = M - mk
    if i == j and mi == mj:
        continue
    if k == l and mk == ml:
        continue
    cands.append((i, mi, j, mj, k, mk, l, ml))

print(f"Trying up to {len(cands)} random MEs for |Z_4|>1e-5 ...")
t0 = time.time()
hits = []
for tup in cands:
    if len(hits) >= nsamp:
        break
    z4 = Zm_4index(*tup)
    if abs(z4) < 1e-5:
        continue
    zf = Zm_fold(*tup)
    hits.append((abs(z4), z4, zf, tup))
    print(
        f"  hit {len(hits)}/{nsamp} |4|={abs(z4):.3e} "
        f"({time.time()-t0:.0f}s)",
        flush=True,
    )

if len(hits) < max(5, nsamp // 4):
    print(f"FAIL — too few hits ({len(hits)})")
    sys.exit(1)

rats: Counter = Counter()
maxd = 0.0
n1 = 0
for _, z4, zf, tup in hits:
    r = zf / z4 if abs(z4) > 1e-12 else float("nan")
    rats[round(r, 4)] += 1
    maxd = max(maxd, abs(zf - z4))
    if abs(r - 1.0) < 1e-3:
        n1 += 1
    i, mi, j, mj, k, mk, l, ml = tup
    print(
        f"  ({i},{mi})({j},{mj})({k},{mk})({l},{ml}) "
        f"4idx={z4:.5e} fold={zf:.5e} r={r:.5g}"
    )

print(
    f"n1={n1}/{len(hits)} rats={rats.most_common(4)} maxΔ={maxd:.3e} "
    f"({time.time()-t0:.1f}s)"
)
if maxd < tol and n1 == len(hits):
    print(f"PASS — m_4index ≡ fold(χ^η)  (full Z, λ={lam})")
    sys.exit(0)
print("FAIL")
sys.exit(1)
