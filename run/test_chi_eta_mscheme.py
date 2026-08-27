#!/usr/bin/env python3
"""χ^η (Γ^{III_a/b} intermediate): m ≡ AMC direct ≡ Pandya Path B.

Physics (analyze):
  χ^η_ijkl = Σ_ab (n̄_a n_b n̄_k + n_a n̄_b n_k) Ω_iabl Ω_bjka
  T×T→S ring (same class as χ^γ / z_ring; different legs + occ).

AMC: learn/amc_tts/factored_GIIIa/input/chi_eta_{direct,via_pandya}.txt
  Path A: chi_eta_direct_ninej.tex
  Path B: fwd Pandya → mid RME (DGEMM) → inv Pandya; chi = +barChi
          (drop AMC-sample overall minus — tts_ring / chi_gamma rule)

Packaging: Ω reduced tensor; χ_red = S/Ĵ (reduce=true).

Usage:
  PYTHONPATH=build python3 -B run/test_chi_eta_mscheme.py [emax=1] [lambda=2]
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


def occ(a: int) -> float:
    return ms.GetOrbit(a).occ


def nbar(a: int) -> float:
    return 1.0 - occ(a)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def tbme(J1, J2, a, b, c, d) -> float:
    return Omega.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


def om(*x) -> float:
    return ut.GetMschemeMatrixElement_2b(Omega, *x)


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def w_eta(a: int, b: int, k: int) -> float:
    """(n̄_a n_b n̄_k + n_a n̄_b n_k)."""
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


# ---------------------------------------------------------------------------
# m-scheme oracle
# ---------------------------------------------------------------------------
def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    w = w_eta(a, b, k)
                    if abs(w) < 1e-12:
                        continue
                    o1 = om(i, mi, a, ma, b, mb, l, ml)  # Ω_iabl
                    if abs(o1) < 1e-16:
                        continue
                    mu = 0.5 * (mi + ma - mb - ml)
                    if abs(mu) > lam:
                        continue
                    cg = CG(lam, mu, lam, -mu, 0, 0)
                    if abs(cg) < 1e-16:
                        continue
                    o2 = om(b, mb, j, mj, k, mk, a, ma)  # Ω_bjka
                    sm += cg * w * o1 * o2
    return sm


def chi_red_from_m(i, j, k, l, J0) -> float:
    """χ_red = Ĵ^{-1} Σ CG CG χ(m). Do not Pauli-skip — χ^η is not AS."""
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


# ---------------------------------------------------------------------------
# Path A — AMC direct (chi_eta_direct_ninej.tex)
# ---------------------------------------------------------------------------
def chi_path_A(i, j, k, l, J0) -> float:
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    tot = 0.0
    maxJ = max(j2i(o) for o in orbits)
    for a in orbits:
        for b in orbits:
            w = w_eta(a, b, k)
            if abs(w) < 1e-12:
                continue
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(i) - j2i(a)) // 2, (j2i(i) + j2i(a)) // 2 + 1):
                for J3 in range(abs(j2i(b) - j2i(l)) // 2, (j2i(b) + j2i(l)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = tbme(J2, J3, i, a, b, l)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(b) - j2i(j)) // 2, (j2i(b) + j2i(j)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(k) - j2i(a)) // 2, (j2i(k) + j2i(a)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = tbme(J4, J5, b, j, k, a)
                            if abs(o2) < 1e-16:
                                continue
                            j0max = int(max(J2, J3, J4, J5, ji, jk, ja, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                s1 = SixJ(J3, lam, J2, ja, ji, j0)
                                s2 = SixJ(J4, lam, J5, ja, jk, j0)
                                n9 = NineJ(jl, jb, J3, jk, J4, j0, J0, jj, ji)
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                tot += (
                                    w
                                    * phase(J2 + J4 + lam)
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (2 * j0 + 1)
                                    / hat(lam)
                                    * s1
                                    * s2
                                    * n9
                                    * o1
                                    * o2
                                )
    return -phase(ji + jk) * hat(J0) * tot


# ---------------------------------------------------------------------------
# Path B — AMC Pandya → mid RME → inv (corrected chi = +barChi)
# ---------------------------------------------------------------------------
def pandya_bar_amc(i, j, k, l, Jbra, Jket) -> float:
    """AMC same-label fwd Pandya (eq1 of via_pandya ninej)."""
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not tri(Jbra, Jket, lam):
        return 0.0
    sm = 0.0
    for J2 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        for J3 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if not tri(J2, J3, lam):
                continue
            n9 = NineJ(lam, Jbra, Jket, J3, jl, jk, J2, ji, jj)
            if abs(n9) < 1e-16:
                continue
            sm += phase(J2) * hat(J2) * hat(J3) * n9 * tbme(J2, J3, i, j, k, l)
    return -phase(Jbra + ji + jk + lam) * hat(Jbra) * hat(Jket) * sm


bar_cache: dict = {}


def get_bar(i, j, k, l, Jbra, Jket):
    key = (i, j, k, l, Jbra, Jket)
    if key not in bar_cache:
        bar_cache[key] = pandya_bar_amc(i, j, k, l, Jbra, Jket)
    return bar_cache[key]


def bar_chi_mid(i, j, k, l, J0) -> float:
    """barχ^η mid RME (eq2): Ĵ^{-1}(-1)^J Σ w (−1)^{J2+λ}/λ̂ barΩ_iabl barΩ_bjka."""
    sm = 0.0
    maxJ = max(j2i(o) for o in orbits) + lam + 2
    for J2 in range(0, maxJ):
        if not tri(J0, J2, lam):
            continue
        for a in orbits:
            for b in orbits:
                w = w_eta(a, b, k)
                if abs(w) < 1e-12:
                    continue
                sm += (
                    w
                    * phase(J2 + lam)
                    / hat(lam)
                    * get_bar(i, a, b, l, J0, J2)
                    * get_bar(b, j, k, a, J2, J0)
                )
    return phase(J0) / hat(J0) * sm


def chi_path_B(i, j, k, l, J0) -> float:
    """Corrected inv (eq3): χ = +Ĵ0 Σ Ĵ2 {jl jk J0; jj ji J2} barχ^{J2}."""
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    sm = 0.0
    Jpmax = max((j2i(i) + j2i(j)) // 2, (j2i(k) + j2i(l)) // 2) + 1
    for Jp in range(0, Jpmax + 1):
        six = SixJ(jl, jk, J0, jj, ji, Jp)
        if abs(six) < 1e-16:
            continue
        sm += hat(Jp) * six * bar_chi_mid(i, j, k, l, Jp)
    return hat(J0) * sm  # corrected (no overall minus)


# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
stats = {
    k: {"max": 0.0, "n": 0, "rats": Counter()} for k in ("A", "B", "A_vs_B")
}
worst = {"A": None, "B": None}

for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                Jmin = max(abs(j2i(i) - j2i(j)) // 2, abs(j2i(k) - j2i(l)) // 2)
                Jmax = min((j2i(i) + j2i(j)) // 2, (j2i(k) + j2i(l)) // 2)
                for J0 in range(Jmin, Jmax + 1):
                    if not (
                        tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                    ):
                        continue
                    Zm = chi_red_from_m(i, j, k, l, J0)
                    ZA = chi_path_A(i, j, k, l, J0)
                    ZB = chi_path_B(i, j, k, l, J0)

                    for tag, Zj in (("A", ZA), ("B", ZB)):
                        if abs(Zm) < 1e-12 and abs(Zj) < 1e-12:
                            continue
                        err = abs(Zm - Zj)
                        stats[tag]["n"] += 1
                        stats[tag]["max"] = max(stats[tag]["max"], err)
                        if abs(Zj) > 1e-12:
                            stats[tag]["rats"][round(Zm / Zj, 6)] += 1
                        if worst[tag] is None or err > worst[tag][0]:
                            worst[tag] = (err, i, j, k, l, J0, Zm, Zj)

                    if abs(ZA) > 1e-12 or abs(ZB) > 1e-12:
                        stats["A_vs_B"]["n"] += 1
                        stats["A_vs_B"]["max"] = max(
                            stats["A_vs_B"]["max"], abs(ZA - ZB)
                        )
                        if abs(ZB) > 1e-12:
                            stats["A_vs_B"]["rats"][round(ZA / ZB, 6)] += 1

print("\n=== χ^η_red  m vs AMC Path A / Path B ===")
for tag in ("A", "B"):
    st = stats[tag]
    print(
        f"Path {tag}: n={st['n']}  max|χ_m−χ_{tag}|={st['max']:.3e}  "
        f"Zm/Zj={st['rats'].most_common(3)}"
    )
    if worst[tag]:
        err, i, j, k, l, J0, Zm, Zj = worst[tag]
        print(
            f"  worst: ijkl=({i},{j},{k},{l}) J0={J0} "
            f"Zm={Zm:.6e} Z{tag}={Zj:.6e} Δ={err:.3e}"
        )

st = stats["A_vs_B"]
print(
    f"\nPath A vs B: n={st['n']}  max|A−B|={st['max']:.3e}  "
    f"A/B={st['rats'].most_common(3)}"
)

pass_A = stats["A"]["max"] < tol
pass_B = stats["B"]["max"] < tol
pass_AB = st["max"] < tol

print("\nVerdict:")
print(f"  m ≡ AMC direct (Path A):  {'PASS' if pass_A else 'FAIL'}")
print(f"  m ≡ Pandya+mid+inv (B):   {'PASS' if pass_B else 'FAIL'}")
print(f"  Path A ≡ Path B:          {'PASS' if pass_AB else 'FAIL'}")
ok = pass_A and pass_B and pass_AB
print("\nPASS — χ^η m ≡ AMC direct ≡ Pandya Path B" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
