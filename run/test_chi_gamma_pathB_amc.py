#!/usr/bin/env python3
"""χ^γ / cross: m ≡ AMC direct ≡ AMC Path B (literal equations).

Packaging: reduce=true on χ/barχ. Path B uses corrected inv (chi = barChi,
not the AMC-sample chi = -barChi which gives Path B = −direct).

AMC: learn/amc_tts/factored_fIIIa/output/chi_gamma_{direct,via_pandya}_plain.tex

Usage:
  PYTHONPATH=build python3 run/test_chi_gamma_pathB_amc.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-7
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
if not Omega.IsReduced():
    Omega.MakeReduced()

orbits = list(ms.all_orbits)
print(
    f"Ω: λ={Omega.GetJRank()} reduced={Omega.IsReduced()} "
    f"||2b||={Omega.TwoBodyNorm():.6g} emax={emax}"
)


def phase(x: float) -> float:
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


def tbme(*a) -> float:
    return Omega.TwoBody.GetTBME_J(*a)


def om(*x) -> float:
    return ut.GetMschemeMatrixElement_2b(Omega, *x)


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def w_gamma(a, b, j, k) -> float:
    return occ(a) * nbar(b) * occ(j) * nbar(k) - nbar(a) * occ(b) * nbar(j) * occ(k)


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    w = w_gamma(a, b, j, k)
                    if abs(w) < 1e-12:
                        continue
                    o1 = om(a, ma, j, mj, k, mk, b, mb)
                    if abs(o1) < 1e-16:
                        continue
                    mu = 0.5 * (ma + mj - mk - mb)
                    if abs(mu) > lam:
                        continue
                    wk = CG(lam, mu, lam, -mu, 0, 0)
                    if abs(wk) < 1e-16:
                        continue
                    sm += wk * w * o1 * om(i, mi, b, mb, a, ma, l, ml)
    return sm


def red_from_m(i, j, k, l, J0) -> float:
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
                sm += c1 * c2 * chi_m(i, mi, j, mj, k, mk, l, ml)
    return sm / hat(J0)


def amc_direct(i, j, k, l, J0) -> float:
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    tot = 0.0
    maxJ = max(j2i(o) for o in orbits)
    for a in orbits:
        for b in orbits:
            w = w_gamma(a, b, j, k)
            if abs(w) < 1e-12:
                continue
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(a) - j2i(j)) // 2, (j2i(a) + j2i(j)) // 2 + 1):
                for J3 in range(
                    abs(j2i(k) - j2i(b)) // 2, (j2i(k) + j2i(b)) // 2 + 1
                ):
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
                                        w
                                        * phase(J3 + J5 + lam)
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


bar_cache: dict = {}


def bar_amc(i, j, k, l, Jbra, Jket) -> float:
    """AMC Path B eq1 (tensor Pandya)."""
    key = (i, j, k, l, Jbra, Jket)
    if key in bar_cache:
        return bar_cache[key]
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not tri(Jbra, Jket, lam):
        bar_cache[key] = 0.0
        return 0.0
    sm = 0.0
    maxJ = max(j2i(o) for o in orbits)
    for J2 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        for J3 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if not tri(J2, J3, lam):
                continue
            o = tbme(J2, J3, i, j, k, l)
            if abs(o) < 1e-16:
                continue
            for j0_2 in range(0, 2 * (maxJ + lam + 2) + 1):
                j0 = 0.5 * j0_2
                six = (
                    SixJ(lam, J3, J2, ji, jj, j0)
                    * SixJ(ji, jl, Jbra, jk, j0, J3)
                    * SixJ(Jket, Jbra, lam, j0, jj, jk)
                )
                if abs(six) < 1e-16:
                    continue
                sm += phase(J2) * hat(J2) * hat(J3) * (2 * j0 + 1) * six * o
    bar_cache[key] = phase(Jbra + ji + jk + lam) * hat(Jbra) * hat(Jket) * sm
    return bar_cache[key]


def barchi_amc(i, j, k, l, Jp) -> float:
    """AMC Path B eq2 (reduce=true)."""
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_gamma(a, b, j, k)
            if abs(w) < 1e-12:
                continue
            for J2 in range(0, max(j2i(o) for o in orbits) + lam + 2):
                if not tri(J2, Jp, lam):
                    continue
                sm += (
                    phase(J2 + lam)
                    / hat(lam)
                    * w
                    * bar_amc(a, j, k, b, J2, Jp)
                    * bar_amc(i, b, a, l, Jp, J2)
                )
    return phase(Jp) / hat(Jp) * sm


def amc_pathB(i, j, k, l, J0) -> float:
    """AMC Path B eq3 corrected: no overall minus (chi = barChi)."""
    tot = 0.0
    Jpmax = max((j2i(i) + j2i(l)) // 2, (j2i(j) + j2i(k)) // 2) + lam + 1
    for Jp in range(0, Jpmax + 1):
        six = SixJ(jo(l), jo(k), J0, jo(j), jo(i), Jp)
        if abs(six) < 1e-14:
            continue
        tot += hat(Jp) * six * barchi_amc(i, j, k, l, Jp)
    return hat(J0) * tot  # corrected: + Ĵ0 Σ … (not −)


stats = {k: {"n": 0, "max": 0.0, "rats": Counter()} for k in ("mA", "mB", "AB")}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                oi, oj, ok, ol = (ms.GetOrbit(x) for x in (i, j, k, l))
                if (oi.l + oj.l + ok.l + ol.l) % 2 != 0:
                    continue
                if abs((oi.tz2 + oj.tz2) - (ok.tz2 + ol.tz2)) != 0:
                    continue
                for J0 in range(
                    max(abs(j2i(i) - j2i(j)) // 2, abs(j2i(k) - j2i(l)) // 2),
                    min((j2i(i) + j2i(j)) // 2, (j2i(k) + j2i(l)) // 2) + 1,
                ):
                    Xm = red_from_m(i, j, k, l, J0)
                    XA = amc_direct(i, j, k, l, J0)
                    XB = amc_pathB(i, j, k, l, J0)
                    for tag, va, vb in (("mA", Xm, XA), ("mB", Xm, XB), ("AB", XA, XB)):
                        if abs(va) < 1e-12 and abs(vb) < 1e-12:
                            continue
                        err = abs(va - vb)
                        stats[tag]["n"] += 1
                        stats[tag]["max"] = max(stats[tag]["max"], err)
                        if abs(vb) > 1e-12:
                            stats[tag]["rats"][round(va / vb, 6)] += 1

print("\n=== χ^γ (occ): reduced χ^J ===")
ok = True
for tag, lab in (
    ("mA", "m vs AMC direct"),
    ("mB", "m vs AMC Path B (corrected inv)"),
    ("AB", "AMC direct vs Path B"),
):
    st = stats[tag]
    passed = st["max"] < tol
    ok = ok and passed
    print(
        f"{lab}: n={st['n']} max|Δ|={st['max']:.3e} "
        f"ratios={st['rats'].most_common(2)} => {'PASS' if passed else 'FAIL'}"
    )

print(f"\nVerdict: {'PASS' if ok else 'FAIL'}  (m ≡ AMC direct ≡ AMC Path B)")
sys.exit(0 if ok else 1)
