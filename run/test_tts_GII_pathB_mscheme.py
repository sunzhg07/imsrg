#!/usr/bin/env python3
"""Γ^II gold chain: m-scheme ≡ Path B (−W_OG − V_ζ).

M-scheme unfact:
  χ^{ΩΓ}_ij = 1/2 Σ_abc w Ω_ciab Γ_abcj
  χ^ζ_ij    = 1/2 Σ_abc w Γ_aibc Ω_bcaj
  Γ^II = −(1−P_ij)[χ^{ΩΓ}_ja × Ω_iakl]^(0) − (1−P_kl)[χ^ζ_ak × Ω_ijal]^(0)

Path B: AMC χ^{ΩΓ} / χ^ζ, same Wbra/Wket 6j, assemble −W − V.

Docs: learn/amc_tts/factored_GII/NOTES.md

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GII_pathB_mscheme.py [emax=1] [lambda=2]
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
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}"
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


def cg0(mu: int) -> float:
    if abs(mu) > lam:
        return 0.0
    return CG(lam, mu, lam, -mu, 0, 0)


hat_lam_inv = 1.0 / math.sqrt(2 * max(lam, 0) + 1)
max_J = max(j2i(o) for o in orbits)


# ---------------------------------------------------------------------------
# χ^ζ (AMC chi_zeta_analyze.tex) — reduced tensor 1b
# ---------------------------------------------------------------------------
def chi_zeta_J(i, j) -> float:
    """AMC chi_zeta_analyze: χ^ζ_{ij} from Γ_ciab Ω_abcj."""
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)
                if abs(w) < 1e-12:
                    continue
                jc = jo(c)
                for J0 in range(0, max_J + 1):
                    if not (tri(jc, ji, J0) and tri(jo(a), jo(b), J0)):
                        continue
                    g = Gamma.TwoBody.GetTBME_J(J0, J0, c, i, a, b)
                    if abs(g) < 1e-16:
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam) or not tri(jc, jj, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, a, b, c, j)
                        if abs(o) < 1e-16:
                            continue
                        six = SixJ(lam, J1, J0, jc, ji, jj)
                        if abs(six) < 1e-16:
                            continue
                        ph = iphase((j2i(j) + j2i(c)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * g * o
    return 0.5 * sm


def chi_OG_J(i, j) -> float:
    """AMC chi_omega_gamma_analyze: χ^{ΩΓ}_{ij} from Ω_ciab Γ_abcj."""
    ji, jj = jo(i), jo(j)
    if not tri(ji, jj, lam):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            for c in orbits:
                w = occ(a) * occ(b) * nbar(c) + nbar(a) * nbar(b) * occ(c)
                if abs(w) < 1e-12:
                    continue
                jc = jo(c)
                for J0 in range(0, max_J + 1):
                    if not tri(jc, ji, J0):
                        continue
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not tri(jo(a), jo(b), J1):
                            continue
                        if not tri(jc, jj, J1):
                            continue
                        o = Eta.TwoBody.GetTBME_J(J0, J1, c, i, a, b)
                        if abs(o) < 1e-16:
                            continue
                        g = Gamma.TwoBody.GetTBME_J(J1, J1, a, b, c, j)
                        if abs(g) < 1e-16:
                            continue
                        six = SixJ(J0, J1, lam, jj, ji, jc)
                        if abs(six) < 1e-16:
                            continue
                        ph = iphase((j2i(j) + j2i(c)) // 2 + lam + J0)
                        sm += ph * hat(J0) * hat(J1) * six * w * o * g
    return 0.5 * sm


print("Building χ^ζ_J and χ^{ΩΓ}_J ...")
t0 = time.time()
chi = {
    (i, j): chi_zeta_J(i, j)
    for i in orbits
    for j in orbits
    if tri(jo(i), jo(j), lam)
}
chi_og = {
    (i, j): chi_OG_J(i, j)
    for i in orbits
    for j in orbits
    if tri(jo(i), jo(j), lam)
}
print(f"  n={len(chi)}  ({time.time()-t0:.2f}s)")

ChiOp = Operator(ms, lam, 0, 0, 1)
ChiOp.SetNonHermitian()
ChiOG = Operator(ms, lam, 0, 0, 1)
ChiOG.SetNonHermitian()
for (i, j), v in chi.items():
    ChiOp.SetOneBody(i, j, v)
for (i, j), v in chi_og.items():
    ChiOG.SetOneBody(i, j, v)


# ---------------------------------------------------------------------------
# Path B J-scheme: AMC bare + restore (1−P); Γ = W − V
# ---------------------------------------------------------------------------
def W_bare_J(J, i, j, k, l) -> float:
    """AMC G2_Wbra with χ^{ΩΓ}_ja Ω_iakl."""
    sm = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            c = chi_og.get((j, a), 0.0)
            if abs(c) < 1e-16:
                continue
            six = SixJ(J, J2, lam, jo(a), jo(j), jo(i))
            o = Eta.TwoBody.GetTBME_J(J2, J, i, a, k, l)
            if abs(six * o) < 1e-16:
                continue
            ph = iphase((j2i(i) + j2i(a)) // 2 + J2)
            sm += ph * (1.0 / hat(J)) * hat(J2) * hat_lam_inv * six * c * o
    return sm


def V_bare_J(J, i, j, k, l) -> float:
    """AMC G2_Wket_noperm: χ_ak Ω_ijal."""
    sm = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            c = chi.get((a, k), 0.0)
            if abs(c) < 1e-16:
                continue
            six = SixJ(J, J2, lam, jo(a), jo(k), jo(l))
            o = Eta.TwoBody.GetTBME_J(J, J2, i, j, a, l)
            if abs(six * o) < 1e-16:
                continue
            ph = iphase((j2i(l) + j2i(a)) // 2 + J2)
            sm += ph * (1.0 / hat(J)) * hat(J2) * hat_lam_inv * six * c * o
    return sm


def W_J(J, i, j, k, l) -> float:
    """(1−P_ij) W_bare — bra exchange phase (−1)^{J+j_i}."""
    # bare + exchange (i↔j on bra / χ / Ω)
    exch = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            c = chi_og.get((i, a), 0.0)
            if abs(c) < 1e-16:
                continue
            six = SixJ(J, J2, lam, jo(a), jo(i), jo(j))
            o = Eta.TwoBody.GetTBME_J(J2, J, j, a, k, l)
            if abs(six * o) < 1e-16:
                continue
            ph = iphase(J + (j2i(i) + j2i(a)) // 2 + J2)
            exch += ph * (1.0 / hat(J)) * hat(J2) * hat_lam_inv * six * c * o
    return W_bare_J(J, i, j, k, l) + exch


def V_J(J, i, j, k, l) -> float:
    """(1−P_kl) V_bare — ket exchange phase (−1)^{J+j_l}."""
    exch = 0.0
    for a in orbits:
        for J2 in range(0, max_J + 1):
            if not tri(J, J2, lam):
                continue
            c = chi.get((a, l), 0.0)
            if abs(c) < 1e-16:
                continue
            six = SixJ(J, J2, lam, jo(a), jo(l), jo(k))
            o = Eta.TwoBody.GetTBME_J(J, J2, i, j, a, k)
            if abs(six * o) < 1e-16:
                continue
            ph = iphase(J + (j2i(l) + j2i(a)) // 2 + J2)
            exch += ph * (1.0 / hat(J)) * hat(J2) * hat_lam_inv * six * c * o
    return V_bare_J(J, i, j, k, l) + exch


def pathB_Z_J(J, i, j, k, l) -> float:
    return -W_J(J, i, j, k, l) - V_J(J, i, j, k, l)


def pathB_to_m(i, mi, j, mj, k, mk, l, ml) -> float:
    if mi + mj != mk + ml:
        return 0.0
    if i == j and mi == mj:
        return 0.0
    if k == l and mk == ml:
        return 0.0
    M = (mi + mj) // 2
    sm = 0.0
    Jmin = max(abs(j2i(i) - j2i(j)), abs(j2i(k) - j2i(l))) // 2
    Jmax = min(j2i(i) + j2i(j), j2i(k) + j2i(l)) // 2
    for J in range(Jmin, Jmax + 1):
        if abs(M) > J:
            continue
        cab = CG(0.5 * j2i(i), 0.5 * mi, 0.5 * j2i(j), 0.5 * mj, J, M)
        ccd = CG(0.5 * j2i(k), 0.5 * mk, 0.5 * j2i(l), 0.5 * ml, J, M)
        if abs(cab * ccd) < 1e-16:
            continue
        sm += cab * ccd * pathB_Z_J(J, i, j, k, l)
    return sm


# ---------------------------------------------------------------------------
# m-scheme gold: W − V with (1−P) and [χ × Ω]^(0)
# ---------------------------------------------------------------------------
def chi_x_ome(chi_op, chi_i, chi_mi, chi_j, chi_mj, ome) -> float:
    """Scalar product [χ^λ × Ω^λ]^(0)."""
    if abs(ome) < 1e-16:
        return 0.0
    chi = ut.GetMschemeMatrixElement_1b(chi_op, chi_i, chi_mi, chi_j, chi_mj)
    if abs(chi) < 1e-16:
        return 0.0
    cg = cg0((chi_mi - chi_mj) // 2)
    if abs(cg) < 1e-16:
        return 0.0
    return cg * chi * ome


def W_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """(1−P_ij) [χ^{ΩΓ}_ja × Ω_iakl]^(0)."""
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            ome = ut.GetMschemeMatrixElement_2b(Eta, i, mi, a, ma, k, mk, l, ml)
            sm += chi_x_ome(ChiOG, j, mj, a, ma, ome)
            ome = ut.GetMschemeMatrixElement_2b(Eta, j, mj, a, ma, k, mk, l, ml)
            sm -= chi_x_ome(ChiOG, i, mi, a, ma, ome)
    return sm


def V_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """(1−P_kl) [χ^ζ_ak × Ω_ijal]^(0)."""
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            ome = ut.GetMschemeMatrixElement_2b(Eta, i, mi, j, mj, a, ma, l, ml)
            sm += chi_x_ome(ChiOp, a, ma, k, mk, ome)
            ome = ut.GetMschemeMatrixElement_2b(Eta, i, mi, j, mj, a, ma, k, mk)
            sm -= chi_x_ome(ChiOp, a, ma, l, ml, ome)
    return sm


def gII_m(i, mi, j, mj, k, mk, l, ml) -> float:
    return -W_m(i, mi, j, mj, k, mk, l, ml) - V_m(i, mi, j, mj, k, mk, l, ml)


print("\nCompare m_gold (−W_OG−V) vs PathB→m vs Mscheme_fact_GII ...")
t0 = time.time()
max_abs = 0.0
n = 0
rats = Counter()
worst = None

for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for mi in m_range(i):
                    for mj in m_range(j):
                        for mk in m_range(k):
                            for ml in m_range(l):
                                if mi + mj != mk + ml:
                                    continue
                                if i == j and mi == mj:
                                    continue
                                if k == l and mk == ml:
                                    continue
                                zm = gII_m(i, mi, j, mj, k, mk, l, ml)
                                zj = pathB_to_m(i, mi, j, mj, k, mk, l, ml)
                                zf = ut.Mscheme_fact_GII(
                                    Eta, Gamma, i, mi, j, mj, k, mk, l, ml
                                )
                                if abs(zm) < 1e-14 and abs(zj) < 1e-14 and abs(zf) < 1e-14:
                                    continue
                                n += 1
                                err = max(abs(zm - zj), abs(zm - zf), abs(zj - zf))
                                if err > max_abs:
                                    max_abs = err
                                    ratio = (
                                        zm / zj if abs(zj) > 1e-12 else float("nan")
                                    )
                                    worst = (
                                        i, j, k, l, mi, mj, mk, ml, zm, zj, err, ratio
                                    )
                                if abs(zj) > 1e-12:
                                    rats[round(zm / zj, 6)] += 1

print(
    f"  n={n}  max|m−PathB|={max_abs:.3e}  m/B={rats.most_common(5)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst:
    i, j, k, l, mi, mj, mk, ml, zm, zj, err, ratio = worst
    print(
        f"  worst ({i},{j},{k},{l}) m=({mi},{mj},{mk},{ml}): "
        f"m={zm:.6e} B={zj:.6e} Δ={err:.3e} r={ratio:.6g}"
    )

ok = max_abs < tol
print("\nPASS — Γ^II Path B ≡ m-scheme ≡ Mscheme_fact_GII" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
