#!/usr/bin/env python3
"""χ^ι Path B (Pandya→DGEMM→inv) ≡ AMC direct — no RC yet.

Factorized bare bar_CHI_V (before RC / second Ω):
  barχ = Γ̄ · (occ_AbarBC ⊙ Ω̄)

Rectangular Path B (any λ; square at λ=0):
  1) Scalar Pandya(Γ) — 1×6j, equal J (IMSRG adcb)
  2) Tensor Pandya(Ω) — 1×9j, rectangular Jbra≠Jket
  3) DGEMM CHI_V (no transpose, no hΩ / (−1)^{J0+J1}):
       barχ(il;kj)^{J0,J1}
         = Σ_ab occ_AbarBC(a,b,k) · Γ̄^{J0}(il;ab) · Ω̄^{J0,J1}(ab;kj)
  4) Inv Pandya — IMSRG kernel ≡ AMC invPlus (NO printed leading minus)

Gold: AMC analyze direct Ω_bjka Γ_iabl (≡ m). Not arxiv Ω_bika.

Usage:
  PYTHONPATH=build python3 -B run/test_chi_iota_pathB_vs_direct.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
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
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  Path B: Pandya + CHI_V DGEMM + invPlus/imsrg (no RC)"
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


def w_analyze(a, b, k) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def w_AbarBC(a, b, c) -> float:
    """Factorized Pandya occ_AbarBC: n̄_a n_b n_c + n_a n̄_b n̄_c."""
    return nbar(a) * occ(b) * occ(c) + occ(a) * nbar(b) * nbar(c)


def chi_amc_direct(J0, J1, i, j, k, l) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    tot = 0.0
    for a in orbits:
        for b in orbits:
            w = w_analyze(a, b, k)
            if abs(w) < 1e-12:
                continue
            ja, jb = jo(a), jo(b)
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


_barO: dict = {}
_barG: dict = {}


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
            hats = hat(J1) * hat(J2) * hat(Jbra) * hat(Jket)
            tb = Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b)
            sm -= (
                hats
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
        sm -= (
            (2 * Jstd + 1)
            * six
            * Gamma.TwoBody.GetTBME_J(Jstd, Jstd, a, d, c, b)
        )
    _barG[key] = sm
    return sm


def bar_chi_V(i, j, k, l, J0, J1) -> float:
    """CHI_V: Γ̄ · (occ_AbarBC ⊙ Ω̄), rectangular. No hΩ / J-phase."""
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(l), J0) and tri(jo(k), jo(j), J1)):
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
    return sm


def inv_tensor_plus(i, j, k, l, J0, J1, bar_fn) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    tot = 0.0
    for J2 in range(0, max_J + 1):
        for J3 in range(0, max_J + 1):
            if not tri(J2, J3, lam):
                continue
            if not (tri(jo(i), jo(l), J2) and tri(jo(k), jo(j), J3)):
                continue
            bc = bar_fn(i, j, k, l, J2, J3)
            if abs(bc) < 1e-16:
                continue
            nj = NineJ(lam, J0, J1, J3, jo(j), jo(k), J2, jo(i), jo(l))
            if abs(nj) < 1e-16:
                continue
            tot += iphase(J2) * hat(J2) * hat(J3) * nj * bc
    return (
        iphase(J0 + (j2i(i) + j2i(k)) // 2 + lam)
        * hat(J0)
        * hat(J1)
        * tot
    )


def inv_imsrg(i, j, k, l, J0, J1, bar_fn) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    sm = 0.0
    for J3 in range(abs(j2i(i) - j2i(l)) // 2, (j2i(i) + j2i(l)) // 2 + 1):
        for J4 in range(
            max(abs(j2i(k) - j2i(j)) // 2, abs(J3 - lam)),
            min((j2i(k) + j2i(j)) // 2, J3 + lam) + 1,
        ):
            if not tri(J3, J4, lam):
                continue
            bc = bar_fn(i, j, k, l, J3, J4)
            if abs(bc) < 1e-16:
                continue
            if lam == 0:
                if J0 != J1 or J3 != J4:
                    continue
                ninej = (
                    phase_half(jj + jl + J0 + J3)
                    * SixJ(ji, jj, J0, jk, jl, J3)
                    / math.sqrt((2 * J1 + 1) * (2 * J4 + 1))
                )
            else:
                ninej = NineJ(ji, jl, J3, jj, jk, J4, J0, J1, lam)
            if abs(ninej) < 1e-14:
                continue
            hatfactor = math.sqrt(
                (2 * J0 + 1) * (2 * J1 + 1) * (2 * J3 + 1) * (2 * J4 + 1)
            )
            sm += (
                hatfactor
                * phase_half(jj + jl + J1 + J4)
                * ninej
                * bc
            )
    return sm


print("Caching AMC-direct ...")
t0 = time.time()
direct = {}
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
                        v = chi_amc_direct(J0, J1, i, j, k, l)
                        if abs(v) > 1e-12:
                            direct[(J0, J1, i, j, k, l)] = v
print(f"  nonzero={len(direct)} ({time.time()-t0:.1f}s)")

variants = [
    (
        "V_invPlus",
        lambda J0, J1, i, j, k, l: inv_tensor_plus(
            i, j, k, l, J0, J1, bar_chi_V
        ),
    ),
    (
        "V_imsrg",
        lambda J0, J1, i, j, k, l: inv_imsrg(i, j, k, l, J0, J1, bar_chi_V),
    ),
]

print(f"\nComparing ALL {len(direct)} MEs ...")
t0 = time.time()
best = None
for vname, vfn in variants:
    rats = Counter()
    n = 0
    max_abs = 0.0
    for (J0, J1, i, j, k, l), ca in direct.items():
        cb = vfn(J0, J1, i, j, k, l)
        if abs(ca) < 1e-10 and abs(cb) < 1e-10:
            continue
        n += 1
        max_abs = max(max_abs, abs(ca - cb))
        if abs(cb) > 1e-8:
            rats[round(ca / cb, 4)] += 1
    print(f"  {vname:12s} n={n:4d} maxΔ={max_abs:.3e} A/B={rats.most_common(2)}")
    if best is None or max_abs < best[1]:
        best = (vname, max_abs, n, rats)

print(f"  ({time.time()-t0:.1f}s)")
ok = best is not None and best[1] < tol and best[2] > 0
print()
if ok:
    print(f"PASS — Path B ≡ AMC direct ({best[0]}, λ={lam})")
    sys.exit(0)
print(f"FAIL — best {best[0]} maxΔ={best[1]:.3e} (λ={lam})")
sys.exit(1)
