#!/usr/bin/env python3
"""Cross χ = Σ_ab Ω_ajkb Ω_ibal (no occ): m-scheme vs AMC Path A (reduce=true).

Packaging (locked — same as ring / test_z_ring_mscheme_sign.py):
  χ(m) = Σ CG(λ,μ; λ,−μ; 0,0) · Ω_ajkb(m) · Ω_ibal(m)
  χ_red(J) = Σ_m CG_ij CG_kl χ(m) / Ĵ

AMC: learn/amc_tts/factored_fIIIa/output/omega_cross_noocc_direct_plain.tex

Usage:
  PYTHONPATH=build python3 run/test_omega_cross_mscheme.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-7
seed = 7

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
if not Omega.IsReduced():
    Omega.MakeReduced()

orbits = list(ms.all_orbits)
print(
    f"Ω: λ={Omega.GetJRank()} reduced={Omega.IsReduced()} "
    f"||2b||={Omega.TwoBodyNorm():.6g} emax={emax} seed={seed}"
)


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def hat(J: float) -> float:
    return math.sqrt(2 * J + 1)


def jo(a: int) -> float:
    return ms.GetOrbit(a).j2 * 0.5


def j2i(a: int) -> int:
    return ms.GetOrbit(a).j2


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def tbme(J1, J2, a, b, c, d) -> float:
    return Omega.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


def om(*x) -> float:
    return ut.GetMschemeMatrixElement_2b(Omega, *x)


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def chi_mscheme(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    o1 = om(a, ma, j, mj, k, mk, b, mb)
                    if abs(o1) < 1e-16:
                        continue
                    mu = 0.5 * (ma + mj - mk - mb)
                    if abs(mu) > lam:
                        continue
                    w = CG(lam, mu, lam, -mu, 0, 0)
                    if abs(w) < 1e-16:
                        continue
                    sm += w * o1 * om(i, mi, b, mb, a, ma, l, ml)
    return sm


def chi_reduced_from_m(i, j, k, l, J0) -> float:
    """χ_red = 1/Ĵ Σ_m CG_ij CG_kl χ(m). Do not Pauli-skip m's — AMC χ is not AS."""
    sm = 0.0
    for mi in mrange(i):
        for mj in mrange(j):
            M = 0.5 * (mi + mj)
            if abs(M) > J0:
                continue
            c1 = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(c1) < 1e-15:
                continue
            for mk in mrange(k):
                ml = mi + mj - mk
                if ml not in mrange(l):
                    continue
                c2 = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(c2) < 1e-15:
                    continue
                sm += c1 * c2 * chi_mscheme(i, mi, j, mj, k, mk, l, ml)
    return sm / hat(J0)


def chi_path_A(i, j, k, l, J0) -> float:
    """AMC direct plain 5×6j with reduce=true (overall Ĵ0)."""
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    tot = 0.0
    maxJ = max(j2i(o) for o in orbits)
    for a in orbits:
        for b in orbits:
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(a) - j2i(j)) // 2, (j2i(a) + j2i(j)) // 2 + 1):
                for J3 in range(abs(j2i(k) - j2i(b)) // 2, (j2i(k) + j2i(b)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = tbme(J2, J3, a, j, k, b)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(i) - j2i(b)) // 2, (j2i(i) + j2i(b)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(a) - j2i(l)) // 2, (j2i(a) + j2i(l)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = tbme(J4, J5, i, b, a, l)
                            if abs(o2) < 1e-16:
                                continue
                            for J6 in range(maxJ + 1):
                                for j0_2 in range(0, 2 * (maxJ + lam + 2) + 1):
                                    j0 = 0.5 * j0_2
                                    six = (
                                        SixJ(J3, lam, J2, ja, jj, j0)
                                        * SixJ(J4, lam, J5, ja, jl, j0)
                                        * SixJ(jj, jk, J6, jb, j0, J3)
                                        * SixJ(jl, ji, J6, jb, j0, J4)
                                        * SixJ(jk, jl, J0, ji, jj, J6)
                                    )
                                    if abs(six) < 1e-16:
                                        continue
                                    tot += (
                                        phase(J3 + J5 + lam)
                                        * hat(J2)
                                        * hat(J3)
                                        * hat(J4)
                                        * hat(J5)
                                        * (2 * J6 + 1)
                                        * (2 * j0 + 1)
                                        / hat(lam)
                                        * six
                                        * o1
                                        * o2
                                    )
    return -phase(jj + jl) * hat(J0) * tot


maxd = 0.0
n = 0
rats = Counter()
worst = None

for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                oi, oj, ok, ol = (ms.GetOrbit(x) for x in (i, j, k, l))
                if (oi.l + oj.l + ok.l + ol.l) % 2 != 0:
                    continue
                if abs((oi.tz2 + oj.tz2) - (ok.tz2 + ol.tz2)) != 0:
                    continue
                Jmin = max(abs(j2i(i) - j2i(j)) // 2, abs(j2i(k) - j2i(l)) // 2)
                Jmax = min((j2i(i) + j2i(j)) // 2, (j2i(k) + j2i(l)) // 2)
                for J0 in range(Jmin, Jmax + 1):
                    Xm = chi_reduced_from_m(i, j, k, l, J0)
                    XA = chi_path_A(i, j, k, l, J0)
                    if abs(Xm) < 1e-12 and abs(XA) < 1e-12:
                        continue
                    err = abs(Xm - XA)
                    n += 1
                    maxd = max(maxd, err)
                    if abs(XA) > 1e-12:
                        rats[round(Xm / XA, 6)] += 1
                    if worst is None or err > worst[0]:
                        worst = (err, i, j, k, l, J0, Xm, XA)

ok = maxd < tol and n > 0
print("\n=== reduced χ^J  (m with [Ω×Ω]_0) vs AMC Path A ===")
print(f"n={n}  max|χ_m−χ_A|={maxd:.3e}  Xm/XA={rats.most_common(3)}")
if worst:
    err, i, j, k, l, J0, Xm, XA = worst
    print(
        f"  worst: ijkl=({i},{j},{k},{l}) J0={J0} "
        f"Xm={Xm:.6e} XA={XA:.6e} Δ={err:.3e}"
    )
print(f"\nVerdict: {'PASS' if ok else 'FAIL'}  (m ≡ AMC direct, reduced packaging)")
sys.exit(0 if ok else 1)
