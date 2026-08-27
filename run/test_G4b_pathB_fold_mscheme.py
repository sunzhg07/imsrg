#!/usr/bin/env python3
"""Γ^{IV_b} full Z (term1+term2) — WRAP-UP gold.

  m-scheme  ≡  AMC direct  ≡  Path B

Chain:
  χ_m ≡ χ_AMC (test_chi_iota_m_vs_amc.py)
  χ_PathB ≡ χ_AMC (test_chi_iota_pathB_vs_direct.py)
  Z[χ_PathB] ≡ Z[χ_m]   ← this bench

Pipeline (no Factorized RC — that dual is λ=0 only so far):
  1) Path B χ̄ = Γ̄·(occ⊙Ω̄) → invPlus → χ_J  (locked ≡ m ≡ AMC all λ)
  2) Unpack χ_J → m (WE)
  3) Analyze fold: W = Σ_ab (χ_aibk Ω_jbla − χ_akbi Ω_jalb) = W1 − W2
  4) Z = (1−P_ij)(1−P_kl) W
  5) lock Z(PathB) ≡ Z(m)

Usage:
  PYTHONPATH=build python3 -B run/test_G4b_pathB_fold_mscheme.py [emax=1] [lambda=2] [nsamp=20]
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
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 20
tol = 1e-5
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
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
    f"  WRAP-UP: Path B fold ≡ m  (⇒ m ≡ AMC_direct ≡ Path B for full Z)"
)


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def phase_half(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


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


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def w_AbarBC(a, b, c) -> float:
    return nbar(a) * occ(b) * occ(c) + occ(a) * nbar(b) * nbar(c)


def w_occ(a, b, k) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


# ---------------------------------------------------------------------------
# Path B χ (locked)
# ---------------------------------------------------------------------------
_barO: dict = {}
_barG: dict = {}
_barV: dict = {}
_chiJ: dict = {}


def bar_Omega(a, b, c, d, Jbra, Jket) -> float:
    key = (a, b, c, d, Jbra, Jket)
    if key in _barO:
        return _barO[key]
    if not tri(Jbra, Jket, lam):
        _barO[key] = 0.0
        return 0.0
    oa, ob, oc, od = (ms.GetOrbit(x) for x in (a, b, c, d))
    ja, jb, jc, jd = jo(a), jo(b), jo(c), jo(d)
    sm = 0.0
    for J1 in range(abs(oa.j2 - od.j2) // 2, (oa.j2 + od.j2) // 2 + 1):
        j2min = max(abs(oc.j2 - ob.j2) // 2, abs(J1 - lam))
        j2max = min((oc.j2 + ob.j2) // 2, J1 + lam)
        for J2 in range(j2min, j2max + 1):
            ninej = NineJ(ja, jd, J1, jb, jc, J2, Jbra, Jket, lam)
            if abs(ninej) < 1e-14:
                continue
            tb = Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b)
            sm -= (
                hat(J1)
                * hat(J2)
                * hat(Jbra)
                * hat(Jket)
                * phase_half((ob.j2 + od.j2) / 2 + Jket + J2)
                * ninej
                * tb
            )
    _barO[key] = sm
    return sm


def bar_Gamma(a, b, c, d, Jcc) -> float:
    key = (a, b, c, d, Jcc)
    if key in _barG:
        return _barG[key]
    oa, ob, oc, od = (ms.GetOrbit(x) for x in (a, b, c, d))
    ja, jb, jc, jd = jo(a), jo(b), jo(c), jo(d)
    jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
    jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
    dJ = 1
    if a == d or b == c:
        dJ = 2
        jmin += jmin % 2
    sm = 0.0
    for Jstd in range(jmin, jmax + 1, dJ):
        six = SixJ(ja, jb, Jcc, jc, jd, Jstd)
        if abs(six) < 1e-8:
            continue
        sm -= (2 * Jstd + 1) * six * Gamma.TwoBody.GetTBME_J(
            Jstd, Jstd, a, d, c, b
        )
    _barG[key] = sm
    return sm


def bar_chi_V(i, j, k, l, J0, J1) -> float:
    """Γ̄(il;ab)·(occ⊙Ω̄)(ab;kj) — labels as Path B locked."""
    key = (i, j, k, l, J0, J1)
    if key in _barV:
        return _barV[key]
    if not tri(J0, J1, lam):
        _barV[key] = 0.0
        return 0.0
    if not (tri(jo(i), jo(l), J0) and tri(jo(k), jo(j), J1)):
        _barV[key] = 0.0
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_AbarBC(a, b, k)
            if abs(w) < 1e-12:
                continue
            if not tri(jo(a), jo(b), J0):
                continue
            bg = bar_Gamma(i, l, a, b, J0)
            if abs(bg) < 1e-16:
                continue
            bo = bar_Omega(a, b, k, j, J0, J1)
            if abs(bo) < 1e-16:
                continue
            sm += w * bg * bo
    _barV[key] = sm
    return sm


def inv_plus(i, j, k, l, J0, J1) -> float:
    """AMC/IMSRG invPlus of bar_CHI_V (locked)."""
    key = (J0, J1, i, j, k, l)
    if key in _chiJ:
        return _chiJ[key]
    if not tri(J0, J1, lam):
        _chiJ[key] = 0.0
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        _chiJ[key] = 0.0
        return 0.0
    tot = 0.0
    for J2 in range(0, max_J + 1):
        for J3 in range(0, max_J + 1):
            if not tri(J2, J3, lam):
                continue
            if not (tri(jo(i), jo(l), J2) and tri(jo(k), jo(j), J3)):
                continue
            bc = bar_chi_V(i, j, k, l, J2, J3)
            if abs(bc) < 1e-16:
                continue
            nj = NineJ(lam, J0, J1, J3, jo(j), jo(k), J2, jo(i), jo(l))
            if abs(nj) < 1e-16:
                continue
            tot += iphase(J2) * hat(J2) * hat(J3) * nj * bc
    v = (
        iphase(J0 + (j2i(i) + j2i(k)) // 2 + lam)
        * hat(J0)
        * hat(J1)
        * tot
    )
    _chiJ[key] = v
    return v


def chi_we_pb(i, mi, j, mj, k, mk, l, ml) -> float:
    """WE unpack of Path B χ_J → m."""
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
            if abs(M1) > J1 or not tri(J0, J1, lam):
                continue
            ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            cj = (
                1.0
                if (lam == 0 and J0 == J1 and M0 == M1)
                else CG(J1, M1, lam, mu, J0, M0)
            )
            if abs(cj) < 1e-15:
                continue
            v = inv_plus(i, j, k, l, J0, J1)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


# ---------------------------------------------------------------------------
# m gold χ + fold
# ---------------------------------------------------------------------------
_chi_m: dict = {}


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    key = (i, mi, j, mj, k, mk, l, ml)
    if key in _chi_m:
        return _chi_m[key]
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
    _chi_m[key] = sm
    return sm


def Wm_from_chi(chi_fn, i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    c1 = chi_fn(a, ma, i, mi, b, mb, k, mk)
                    if abs(c1) > 1e-16:
                        o1 = ut.GetMschemeMatrixElement_2b(
                            Eta, j, mj, b, mb, l, ml, a, ma
                        )
                        if abs(o1) > 1e-16:
                            sm += c1 * o1
                    c2 = chi_fn(a, ma, k, mk, b, mb, i, mi)
                    if abs(c2) > 1e-16:
                        o2 = ut.GetMschemeMatrixElement_2b(
                            Eta, j, mj, a, ma, l, ml, b, mb
                        )
                        if abs(o2) > 1e-16:
                            sm -= c2 * o2
    return sm


def Zm_from_chi(chi_fn, i, mi, j, mj, k, mk, l, ml) -> float:
    w = Wm_from_chi(chi_fn, i, mi, j, mj, k, mk, l, ml)
    w -= Wm_from_chi(chi_fn, j, mj, i, mi, k, mk, l, ml)
    w -= Wm_from_chi(chi_fn, i, mi, j, mj, l, ml, k, mk)
    w += Wm_from_chi(chi_fn, j, mj, i, mi, l, ml, k, mk)
    return w


# ---------------------------------------------------------------------------
# Sample MEs
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

print(f"Trying up to {len(cands)} random MEs for |Zm|>1e-5 ...")
t0 = time.time()
hits = []
for tup in cands:
    if len(hits) >= nsamp:
        break
    zm = Zm_from_chi(chi_m, *tup)
    if abs(zm) < 1e-5:
        continue
    zp = Zm_from_chi(chi_we_pb, *tup)
    hits.append((abs(zm), zm, zp, tup))
    print(
        f"  hit {len(hits)}/{nsamp} |m|={abs(zm):.3e} "
        f"({time.time()-t0:.0f}s)",
        flush=True,
    )

if len(hits) < max(5, nsamp // 4):
    print(f"FAIL — too few hits ({len(hits)})")
    sys.exit(1)

rats: Counter = Counter()
maxd = 0.0
n1 = 0
for _, zm, zp, tup in hits:
    r = zp / zm if abs(zm) > 1e-12 else float("nan")
    rats[round(r, 4)] += 1
    maxd = max(maxd, abs(zp - zm))
    if abs(r - 1.0) < 1e-3:
        n1 += 1
    i, mi, j, mj, k, mk, l, ml = tup
    print(
        f"  ({i},{mi})({j},{mj})({k},{mk})({l},{ml}) "
        f"m={zm:.5e} PB={zp:.5e} r={r:.5g}"
    )

print(
    f"n1={n1}/{len(hits)} rats={rats.most_common(4)} maxΔ={maxd:.3e} "
    f"({time.time()-t0:.1f}s)"
)
if maxd < tol and n1 == len(hits):
    print(f"PASS — m ≡ AMC_direct ≡ Path B  (full Z, λ={lam})")
    sys.exit(0)
print("FAIL")
sys.exit(1)
