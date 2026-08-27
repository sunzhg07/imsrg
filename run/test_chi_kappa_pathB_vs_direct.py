#!/usr/bin/env python3
"""χ^κ Path B rectangular ≡ AMC analyze direct (all λ).

Path B (analyze Chi_VI_II):
  1) Scalar Pandya(Γ) — 1×6j, equal J (IMSRG adcb)
  2) Tensor Pandya(Ω) — 1×9j, rectangular Jbra≠Jket (IMSRG adcb)
  3) DGEMM: barχ(il;kj)^{J0,J1} = hΩ (−1)^{J0+J1} Σ_ab occ_ABbarD(a,b,l)
              · Ω̄^{J1,J0}(ab;il) · Γ̄^{J1}(ab;kj)
  4) Inv Pandya tensor — IMSRG kernel or AMC via Eq4 (NO printed leading minus)

Gold: AMC analyze direct (≡ m, tensor/WE Ω).

Usage:
  PYTHONPATH=build python3 -B run/test_chi_kappa_pathB_vs_direct.py [emax=1] [lambda=2]
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

hEta = -1 if Eta.IsAntiHermitian() else (1 if Eta.IsHermitian() else -1)
orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
print(
    f"emax={emax} λ={lam} seed={seed} hΩ={hEta}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  Path B: rectangular Pandya + VI_II DGEMM + tensor invPlus"
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


def w_analyze(a, c, d) -> float:
    return nbar(c) * nbar(d) * occ(a) + occ(c) * occ(d) * nbar(a)


def w_ABbarD(a, b, d) -> float:
    return occ(a) * nbar(b) * occ(d) + nbar(a) * occ(b) * nbar(d)


# ---------------------------------------------------------------------------
# AMC direct gold
# ---------------------------------------------------------------------------
def chi_amc_direct(J0, J1, i, j, b, d) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(b), jo(d), J1)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for c in orbits:
            w = w_analyze(a, c, d)
            if abs(w) < 1e-12:
                continue
            ja, jc = jo(a), jo(c)
            for J2 in range(abs(j2i(a) - j2i(i)) // 2, (j2i(a) + j2i(i)) // 2 + 1):
                for J3 in range(
                    abs(j2i(c) - j2i(d)) // 2, (j2i(c) + j2i(d)) // 2 + 1
                ):
                    if not tri(J2, J3, lam):
                        continue
                    om = Eta.TwoBody.GetTBME_J(J2, J3, a, i, c, d)
                    if abs(om) < 1e-16:
                        continue
                    for J4 in range(
                        max(abs(j2i(j) - j2i(c)), abs(j2i(b) - j2i(a))) // 2,
                        min(j2i(j) + j2i(c), j2i(b) + j2i(a)) // 2 + 1,
                    ):
                        gam = Gamma.TwoBody.GetTBME_J(J4, J4, j, c, b, a)
                        if abs(gam) < 1e-16:
                            continue
                        for J5 in range(
                            max(abs(j2i(a) - j2i(c)), abs(j2i(b) - j2i(j))) // 2,
                            min(j2i(a) + j2i(c), j2i(b) + j2i(j)) // 2 + 1,
                        ):
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
                                ph = -iphase(
                                    J1 + J2 + J3 + J4 + (j2i(i) + j2i(c)) // 2
                                )
                                tot += ph * hats * six * w * om * gam
    return hat(J0) * hat(J1) * tot


# ---------------------------------------------------------------------------
# Pandya (IMSRG adcb) — rectangular for Ω
# ---------------------------------------------------------------------------
_barO: dict = {}
_barG: dict = {}


def bar_Omega(a, b, c, d, Jbra, Jket) -> float:
    """IMSRG DoTensorPandya: <ab Jbra|Ω̄|cd Jket> from Ω(a,d,c,b)."""
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
    """IMSRG scalar Pandya: <ab J|Γ̄|cd J>."""
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


def bar_chi_VII(i, j, k, l, J0, J1) -> float:
    """VI_II channel DGEMM (no extra 9j), rectangular.

    barχ(il; kj)^{J0,J1} = hΩ (−1)^{J0+J1} Σ_ab occ(a,b,l)
         · Ω̄^{J1,J0}(ab;il) · Γ̄^{J1}(ab;kj)

    The (−1)^{J0+J1} is required for λ≠0 (at λ=0, J0=J1 ⇒ +1).
    Scalar Γ forces mid channel = J1.
    """
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(l), J0) and tri(jo(k), jo(j), J1)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_ABbarD(a, b, l)
            if abs(w) < 1e-12:
                continue
            if not tri(jo(a), jo(b), J1):
                continue
            bg = bar_Gamma(a, b, k, j, J1)
            if abs(bg) < 1e-16:
                continue
            bo = bar_Omega(a, b, i, l, J1, J0)
            if abs(bo) < 1e-16:
                continue
            sm += hEta * w * bo * bg
    return iphase(J0 + J1) * sm


def bar_chi_VI(i, j, k, l, J0, J1) -> float:
    """VI = Γ̄·ND: Γ̄^{J0}(il;ab) Ω̄^{J0,J1}(ab;kj) · occ(a,b,j)."""
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(l), J0) and tri(jo(k), jo(j), J1)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_ABbarD(a, b, j)
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
    """AMC tensor inv WITHOUT leading minus (IMSRG / via Eq4).

    χ^{J0 J1 λ}_{ijkl} = (+1)(-1)^{J0+ji+jk+λ} Ĵ0 Ĵ1
      Σ_{J2 J3} (-1)^{J2} Ĵ2 Ĵ3 {λ J0 J1; J3 jj jk; J2 ji jl} barχ^{J2 J3}
    """
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    tot = 0.0
    for J2 in range(0, max_J + 1):
        for J3 in range(0, max_J + 1):
            if not tri(J2, J3, lam):
                continue
            # barχ_ijkl with bra (i,l) J2, ket (k,j) J3
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


def inv_tensor_minus(i, j, k, l, J0, J1, bar_fn) -> float:
    return -inv_tensor_plus(i, j, k, l, J0, J1, bar_fn)


def inv_imsrg(i, j, k, l, J0, J1, bar_fn) -> float:
    """IMSRG AddInverseTensorPandya kernel (direct term only, no AS).

    ninej(ji,jl,J3, jj,jk,J4, J0,J1,λ) with phase (jj+jl+J1+J4),
    hats √{(2J0+1)(2J1+1)(2J3+1)(2J4+1)}, bar at (J3,J4) for (il;kj).
    """
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
            # barχ_ijkl AMC/IMSRG: bra (i,l) J3, ket (k,j) J4
            bc = bar_fn(i, j, k, l, J3, J4)
            if abs(bc) < 1e-16:
                continue
            if lam == 0:
                # IMSRG λ=0 reduction of ninej
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


# ---------------------------------------------------------------------------
# Compare Path B vs AMC direct
# ---------------------------------------------------------------------------
print("Caching AMC-direct ...")
t0 = time.time()
direct = {}
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
                        v = chi_amc_direct(J0, J1, i, j, b, d)
                        if abs(v) > 1e-12:
                            direct[(J0, J1, i, j, b, d)] = v
print(f"  nonzero={len(direct)} ({time.time()-t0:.1f}s)")

variants = [
    ("VII_invPlus", lambda J0, J1, i, j, b, d: inv_tensor_plus(
        i, j, b, d, J0, J1, bar_chi_VII
    )),
    ("VII_invMinus", lambda J0, J1, i, j, b, d: inv_tensor_minus(
        i, j, b, d, J0, J1, bar_chi_VII
    )),
    ("VII_imsrg", lambda J0, J1, i, j, b, d: inv_imsrg(
        i, j, b, d, J0, J1, bar_chi_VII
    )),
    ("VII_imsrg_m", lambda J0, J1, i, j, b, d: -inv_imsrg(
        i, j, b, d, J0, J1, bar_chi_VII
    )),
    ("VI_invPlus", lambda J0, J1, i, j, b, d: inv_tensor_plus(
        i, j, b, d, J0, J1, bar_chi_VI
    )),
    ("VI_imsrg", lambda J0, J1, i, j, b, d: inv_imsrg(
        i, j, b, d, J0, J1, bar_chi_VI
    )),
]

print(f"\nComparing ALL {len(direct)} MEs ...")
t0 = time.time()
best = None
for vname, vfn in variants:
    rats = Counter()
    n = 0
    max_abs = 0.0
    worst = None
    for (J0, J1, i, j, b, d), ca in direct.items():
        cb = vfn(J0, J1, i, j, b, d)
        if abs(ca) < 1e-10 and abs(cb) < 1e-10:
            continue
        n += 1
        err = abs(ca - cb)
        if err > max_abs:
            max_abs = err
            r = ca / cb if abs(cb) > 1e-10 else float("nan")
            worst = (J0, J1, i, j, b, d, ca, cb, r)
        if abs(cb) > 1e-8:
            rats[round(ca / cb, 4)] += 1
    print(
        f"  {vname:14s} n={n:4d} maxΔ={max_abs:.3e} A/B={rats.most_common(3)}"
    )
    if worst and max_abs >= tol:
        J0, J1, i, j, b, d, ca, cb, r = worst
        print(
            f"    worst J={J0},{J1} ({i},{j}|{b},{d}) "
            f"A={ca:.4e} B={cb:.4e} r={r:.4g}"
        )
    if best is None or max_abs < best[1]:
        best = (vname, max_abs, n, rats)

print(f"  ({time.time()-t0:.1f}s)")
ok = best is not None and best[1] < tol and best[2] > 0
print(
    f"\nPASS — Path B ≡ AMC direct ({best[0]}, λ={lam})"
    if ok
    else f"\nFAIL — best {best[0]} maxΔ={best[1]:.3e} (λ={lam})"
)
sys.exit(0 if ok else 1)
