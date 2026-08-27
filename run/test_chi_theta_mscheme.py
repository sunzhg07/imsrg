#!/usr/bin/env python3
"""χ^θ (T×T→S scalar 2b): m-scheme ≡ AMC.

AMC: learn/amc_tts/factored_GIIIc/input/chi_theta{,_reduced}.txt
  Unreduced: χ = (−1)^{J0} Ĵ^{-2} Σ (−1)^{J2+λ} λ̂^{-1} Ω_ijab Ω_abkl
  Reduced:   χ = (−1)^{J0} Ĵ^{-1} Σ (−1)^{J2+λ} λ̂^{-1} Ω_ijab Ω_abkl
  w = f(a,b,k)+f(a,b,j), f(a,b,x)=n_a n_b nbar_x + nbar_a nbar_b n_x

m (physical):
  χ(m) = Σ_ab w CG(λ,μ; λ,−μ; 0,0) Ω_ijab(m) Ω_abkl(m)

Compare:
  χ_red ← Σ_m CG_bra CG_ket χ(m) / Ĵ   ≡ AMC reduce=true  (**lock this**)
  χ_unred = χ_red / Ĵ                  ≡ AMC default
  Note: bare Σ CG CG χ(m) = S = Ĵ χ_red = Ĵ² χ_unred (code ChiTab),
  not χ_unred — do not compare that sum to AMC unreduced.
  Do **not** Pauli-skip equal-m (χ is not AS in AMC).

Code ChiTab stores S = Ĵ χ_red = Ĵ² χ_unred (bare product).

Usage:
  PYTHONPATH=build python3 -B run/test_chi_theta_mscheme.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
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
if not Eta.IsReduced():
    Eta.MakeReduced()

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}"
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


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def cg0(mu: float) -> float:
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1)
max_J = max(j2i(o) for o in orbits)


def w_theta(a, b, j, k):
    return (
        occ(a) * occ(b) * nbar(k)
        + nbar(a) * nbar(b) * occ(k)
        + occ(a) * occ(b) * nbar(j)
        + nbar(a) * nbar(b) * occ(j)
    )


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """Physical m: [Ω×Ω]^(0) with χ^θ occupation."""
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            for b in orbits:
                for mb in m_range(b):
                    w = w_theta(a, b, j, k)
                    if abs(w) < 1e-12:
                        continue
                    o1 = ut.GetMschemeMatrixElement_2b(
                        Eta, i, mi, j, mj, a, ma, b, mb
                    )
                    if abs(o1) < 1e-16:
                        continue
                    # μ from first Ω (tensor): Δm / 2
                    mu = 0.5 * (mi + mj - ma - mb)
                    cg = cg0(mu)
                    if abs(cg) < 1e-16:
                        continue
                    o2 = ut.GetMschemeMatrixElement_2b(
                        Eta, a, ma, b, mb, k, mk, l, ml
                    )
                    sm += w * cg * o1 * o2
    return sm


def project_from_m(i, j, k, l, J0, reduced: bool) -> float:
    """Σ CG CG χ(m) [/ Ĵ if reduced]. No Pauli skip."""
    sm = 0.0
    for mi in m_range(i):
        for mj in m_range(j):
            M = 0.5 * (mi + mj)
            if abs(M) > J0:
                continue
            c1 = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(c1) < 1e-15:
                continue
            for mk in m_range(k):
                ml = mi + mj - mk
                if ml not in m_range(l):
                    continue
                c2 = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(c2) < 1e-15:
                    continue
                sm += c1 * c2 * chi_m(i, mi, j, mj, k, mk, l, ml)
    if reduced:
        return sm / hat(J0)
    return sm


def chi_amc(i, j, k, l, J0, reduced: bool) -> float:
    """AMC χ^θ.tex / chi_theta_reduced.tex."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_theta(a, b, j, k)
            if abs(w) < 1e-12:
                continue
            for J2 in range(0, max_J + 1):
                if not tri(J0, J2, lam):
                    continue
                o1 = Eta.TwoBody.GetTBME_J(J0, J2, i, j, a, b)
                o2 = Eta.TwoBody.GetTBME_J(J2, J0, a, b, k, l)
                if abs(o1 * o2) < 1e-16:
                    continue
                ph = iphase(J0 + J2 + lam)
                sm += ph * hat_lam_inv * w * o1 * o2
    # bare S; then Ĵ^{-1} (red) or Ĵ^{-2} (unred)
    if reduced:
        return iphase(J0) * sm / hat(J0)  # wait: ph already has J0
    # Fix: AMC has (−1)^{J0} * (−1)^{J2+λ} = (−1)^{J0+J2+λ} already in sm via ph
    # so sm already includes (−1)^{J0}. Prefactor is only Ĵ power.
    # Recalculate carefully:
    return None  # replaced below


def chi_amc_S(i, j, k, l, J0) -> float:
    """Bare S = Σ (−1)^{J0+J2+λ} λ̂^{-1} w Ω Ω  (code ChiTab / AMC numerator)."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = w_theta(a, b, j, k)
            if abs(w) < 1e-12:
                continue
            for J2 in range(0, max_J + 1):
                if not tri(J0, J2, lam):
                    continue
                o1 = Eta.TwoBody.GetTBME_J(J0, J2, i, j, a, b)
                o2 = Eta.TwoBody.GetTBME_J(J2, J0, a, b, k, l)
                if abs(o1 * o2) < 1e-16:
                    continue
                sm += iphase(J0 + J2 + lam) * hat_lam_inv * w * o1 * o2
    return sm


print("Compare m→χ_red vs AMC reduced (unred = red/Ĵ) ...")
t0 = time.time()
rats_r = Counter()
rats_u = Counter()
n = 0
max_r = 0.0
max_u = 0.0
worst = None

for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                if len({i, j, k, l}) < 3:
                    continue
                for J0 in range(0, max_J + 1):
                    if not (
                        tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                    ):
                        continue
                    S = chi_amc_S(i, j, k, l, J0)
                    if abs(S) < 1e-10:
                        continue
                    amc_r = S / hat(J0)
                    amc_u = S / (hat(J0) ** 2)
                    m_r = project_from_m(i, j, k, l, J0, reduced=True)
                    m_u = m_r / hat(J0)  # packaging: unred = red / Ĵ
                    n += 1
                    dr = abs(m_r - amc_r)
                    du = abs(m_u - amc_u)
                    if dr > max_r:
                        max_r = dr
                        worst = ("red", i, j, k, l, J0, m_r, amc_r)
                    if du > max_u:
                        max_u = du
                    if abs(amc_r) > 1e-10:
                        rats_r[round(m_r / amc_r, 4)] += 1
                    if abs(amc_u) > 1e-10:
                        rats_u[round(m_u / amc_u, 4)] += 1

print(
    f"  n={n}  max|m−AMC|_red={max_r:.3e}  red ratios={rats_r.most_common(4)}\n"
    f"         max|m−AMC|_unred={max_u:.3e}  unred ratios={rats_u.most_common(4)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst is not None:
    tag, i, j, k, l, J0, m, a = worst
    print(
        f"  worst {tag} ({i},{j},{k},{l}) J={J0}: m={m:.6e} AMC={a:.6e} "
        f"r={m/a if abs(a)>1e-14 else float('nan'):g}"
    )

ok = (
    max_r < tol
    and max_u < tol
    and (not rats_r or rats_r.most_common(1)[0][0] == 1.0)
    and (not rats_u or rats_u.most_common(1)[0][0] == 1.0)
)
print("\nPASS — χ^θ m ≡ AMC (red & unred)" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
