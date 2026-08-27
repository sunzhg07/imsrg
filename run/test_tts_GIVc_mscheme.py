#!/usr/bin/env python3
"""Γ^{IV_c} / χ^λ: m ≡ ring fold ≡ ethS Path B (Pandya→DGEMM).

Path A TTS strips removed from ethS. Gold = m / AMC direct / Pandya+DGEMM.

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIVc_mscheme.py [emax=1] [lambda=2]
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


def chi_amc_J(J0, J1, i, j, k, l) -> float:
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J1)):
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w_l = nbar(a) * nbar(b) * occ(l) + occ(a) * occ(b) * nbar(l)
            w_j = nbar(a) * nbar(b) * occ(j) + occ(a) * occ(b) * nbar(j)
            if abs(w_l) > 1e-12 and tri(jo(a), jo(b), J0):
                sm += w_l * (
                    Gamma.TwoBody.GetTBME_J(J0, J0, i, j, a, b)
                    * Eta.TwoBody.GetTBME_J(J0, J1, a, b, k, l)
                )
            if abs(w_j) > 1e-12 and tri(jo(a), jo(b), J1):
                sm += w_j * (
                    Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b)
                    * Gamma.TwoBody.GetTBME_J(J1, J1, a, b, k, l)
                )
    return sm


print("Building χ^λ cache ...")
t0 = time.time()
chiJ: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        v = chi_amc_J(J0, J1, i, j, k, l)
                        if abs(v) > 1e-16:
                            chiJ[(J0, J1, i, j, k, l)] = v
print(f"  nonzero={len(chiJ)}  ({time.time()-t0:.2f}s)")


def chi_m_WE(i, mi, j, mj, k, mk, l, ml) -> float:
    if i == j and mi == mj:
        return 0.0
    if k == l and mk == ml:
        return 0.0
    if abs(mi + mj - mk - ml) > 2 * lam:
        return 0.0
    M0 = (mi + mj) // 2
    M1 = (mk + ml) // 2
    mu = M0 - M1
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0 or (i == j and J0 % 2 > 0):
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        for J1 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if abs(M1) > J1 or (k == l and J1 % 2 > 0):
                continue
            if not tri(J0, J1, lam):
                continue
            ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            cj = CG(J1, M1, lam, mu, J0, M0)
            if abs(cj) < 1e-15:
                continue
            v = chiJ.get((J0, J1, i, j, k, l), 0.0)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


def phase(x: float) -> float:
    """(-1)^x for half-integer/integer x (same as test_z_ring_mscheme_sign)."""
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def ring_X(p, q, s, r, J0) -> float:
    """Bare T×T→S ring X_pqsr = Σ χ_pbar Ω_aqsb (tts_ring Path A, reduced)."""
    jp, jq, js, jr = jo(p), jo(q), jo(s), jo(r)
    if not (tri(jp, jq, J0) and tri(js, jr, J0)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for b in orbits:
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(p) - j2i(b)) // 2, (j2i(p) + j2i(b)) // 2 + 1):
                for J3 in range(abs(j2i(a) - j2i(r)) // 2, (j2i(a) + j2i(r)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = chiJ.get((J2, J3, p, b, a, r), 0.0)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(a) - j2i(q)) // 2, (j2i(a) + j2i(q)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(s) - j2i(b)) // 2, (j2i(s) + j2i(b)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = Eta.TwoBody.GetTBME_J(J4, J5, a, q, s, b)
                            if abs(o2) < 1e-16:
                                continue
                            j0max = int(max(J2, J3, J4, J5, jp, js, jb, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                s1 = SixJ(J3, lam, J2, jb, jp, j0)
                                s2 = SixJ(J4, lam, J5, jb, js, j0)
                                n9 = NineJ(jr, ja, J3, js, J4, j0, J0, jq, jp)
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                tot += (
                                    phase(J2 + J4 + lam)
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (2 * j0 + 1)
                                    * hat_lam_inv
                                    * s1
                                    * s2
                                    * n9
                                    * o1
                                    * o2
                                )
    return -phase(jp + js) * hat(J0) * tot


def amc_fold(i, j, k, l, J0, reduced: bool) -> float:
    """Ring X + fermionic AS → Z_red; Z_unred = Z_red/Ĵ."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    x = ring_X(i, j, k, l, J0)
    xkl = ring_X(i, j, l, k, J0)
    xij = ring_X(j, i, k, l, J0)
    xijkl = ring_X(j, i, l, k, J0)
    pkl = phase(jo(k) + jo(l) - J0)
    pij = phase(jo(i) + jo(j) - J0)
    z_red = 0.5 * (x - pkl * xkl - pij * xij + pij * pkl * xijkl)
    return z_red if reduced else z_red / hat(J0)

def Z_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            for b in orbits:
                for mb in m_range(b):
                    for sign, cx, ox in [
                        (
                            +1,
                            (i, mi, a, ma, l, ml, b, mb),
                            (b, mb, j, mj, a, ma, k, mk),
                        ),
                        (
                            -1,
                            (i, mi, a, ma, k, mk, b, mb),
                            (b, mb, j, mj, a, ma, l, ml),
                        ),
                        (
                            -1,
                            (j, mj, a, ma, l, ml, b, mb),
                            (b, mb, i, mi, a, ma, k, mk),
                        ),
                        (
                            +1,
                            (j, mj, a, ma, k, mk, b, mb),
                            (b, mb, i, mi, a, ma, l, ml),
                        ),
                    ]:
                        chi = chi_m_WE(*cx)
                        if abs(chi) < 1e-16:
                            continue
                        mu = 0.5 * (cx[1] + cx[3] - cx[5] - cx[7])
                        cg = cg0(mu)
                        if abs(cg) < 1e-16:
                            continue
                        om = ut.GetMschemeMatrixElement_2b(Eta, *ox)
                        sm += sign * cg * chi * om
    return 0.5 * sm


def project_S(i, j, k, l, J0) -> float:
    """Bare S = Σ CG CG Z(m)."""
    sm = 0.0
    for mi in m_range(i):
        for mj in m_range(j):
            if i == j and mi == mj:
                continue
            M = (mi + mj) // 2
            if abs(M) > J0:
                continue
            c1 = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(c1) < 1e-15:
                continue
            for mk in m_range(k):
                ml = mi + mj - mk
                if ml not in m_range(l):
                    continue
                if k == l and mk == ml:
                    continue
                c2 = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(c2) < 1e-15:
                    continue
                sm += c1 * c2 * Z_m(i, mi, j, mj, k, mk, l, ml)
    return sm


def fresh_Z():
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    return Z


print("ethS Path B (Pandya→DGEMM) ...")
t0 = time.time()
cm = Commutator.FactorizedDoubleCommutator_eths
Z_B = fresh_Z()
cm.comm223_232_GIVc(Eta, Gamma, Z_B)
print(f"  ‖Z_B‖={Z_B.TwoBodyNorm():.6g}  ({time.time()-t0:.2f}s)")

print("\nm ≡ ring fold (S/Ĵ = red, S/Ĵ² = unred), exhaustive ...")
t0 = time.time()
rats_r = Counter()
rats_u = Counter()
n = 0
max_r = 0.0
max_u = 0.0
for i in orbits:
    for j in orbits:
        if i >= j:
            continue
        for k in orbits:
            for l in orbits:
                if k >= l:
                    continue
                for J0 in range(0, max_J + 1):
                    if not (
                        tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                    ):
                        continue
                    S = project_S(i, j, k, l, J0)
                    aR = amc_fold(i, j, k, l, J0, reduced=True)
                    aU = amc_fold(i, j, k, l, J0, reduced=False)
                    mR = S / hat(J0)
                    mU = S / (hat(J0) ** 2)
                    if abs(S) < 1e-10 and abs(aR) < 1e-10:
                        continue
                    n += 1
                    max_r = max(max_r, abs(mR - aR))
                    max_u = max(max_u, abs(mU - aU))
                    if abs(aR) > 1e-8:
                        rats_r[round(mR / aR, 4)] += 1
                    if abs(aU) > 1e-8:
                        rats_u[round(mU / aU, 4)] += 1

print(
    f"  n={n}  max|S/Ĵ−red|={max_r:.3e}  mR/aR={rats_r.most_common(3)}\n"
    f"       max|S/Ĵ²−unred|={max_u:.3e}  mU/aU={rats_u.most_common(3)}  "
    f"({time.time()-t0:.2f}s)"
)

# ethS Path B vs m (via ring red → unred packaging)
max_be = 0.0
n_be = 0
for i in orbits:
    for j in orbits:
        if i >= j:
            continue
        for k in orbits:
            for l in orbits:
                if k >= l:
                    continue
                for J0 in range(0, max_J + 1):
                    if not (
                        tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                    ):
                        continue
                    aU = amc_fold(i, j, k, l, J0, reduced=False)
                    bU = Z_B.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
                    if abs(aU) < 1e-10 and abs(bU) < 1e-10:
                        continue
                    n_be += 1
                    max_be = max(max_be, abs(aU - bU))

ok_m = (
    n > 0
    and max_r < tol
    and max_u < tol
    and (not rats_r or rats_r.most_common(1)[0][0] == 1.0)
    and (not rats_u or rats_u.most_common(1)[0][0] == 1.0)
)
ok_be = n_be > 0 and max_be < tol

print(
    f"\nm ≡ ring fold (χ^λ):     {'PASS' if ok_m else 'FAIL'}\n"
    f"ethS Path B ≡ ring:      {'PASS' if ok_be else 'FAIL'}  "
    f"(n={n_be} max|Δ|={max_be:.3e})"
)
ok = ok_m and ok_be
print("\nPASS — Γ^{IV_c} m ≡ ring ≡ ethS Path B" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
