#!/usr/bin/env python3
"""Γ^{IV_a} from Path B χ^κ + ladder RME; (1−P) on W only.

Packaging (REDUCED_UNREDUCED.md — name both sides first):
  Ω  : reduced tensor for **all** λ (WE path; λ=0 ≠ unreduced scalar)
  Γ  : unreduced scalar
  χ^κ: Path B tensor ≡ AMC analyze (WE-reduced); **no P / no herm fill**
  W,Z: reduced scalar — X_red = Σ CG CG X(m) / Ĵ
       AMC ladder: W_red = −(−1)^{J0} Ĵ0^{-1} Σ (−1)^{jb+jd+λ} λ̂^{-1} χ Ω

Pipeline:
  1. Path B χ^κ = invPlus(VI_II)          ← bare χ
  2. W_ijkl  = −Σ_bd χ_ijbd Ω_dbkl        ← DGEMM intermediate
     W_klij  = −Σ_bd χ_klbd Ω_dbij        ← same kernel, legs swapped
  3. Z = (1−P_ij) W + (1−P_kl) W_klij     ← P on W only (not χ, not late Z)

Hermitian assemble: reuse χ^κ via W_klij — do **not** multiply h_Ω again
(that would give A−A^T). Gold: Z = A + A^T in equal-J blocks.

Usage:
  PYTHONPATH=build python3 -B run/test_G4a_pathB_mscheme.py [emax=1] [lambda=2]
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
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 30

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
hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1) if lam else 1.0
print(
    f"emax={emax} λ={lam} seed={seed} hΩ={hEta} nsamp={nsamp}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  Path B χ → W=−χΩ → (1−P) on W → Z (Hermitian W+W_klij)"
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


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def w_analyze(a, c, d) -> float:
    return nbar(c) * nbar(d) * occ(a) + occ(c) * occ(d) * nbar(a)


def w_ABbarD(a, b, d) -> float:
    return occ(a) * nbar(b) * occ(d) + nbar(a) * occ(b) * nbar(d)


# ---------------------------------------------------------------------------
# Path B: Pandya + VI_II DGEMM + invPlus
# ---------------------------------------------------------------------------
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


def bar_chi_VII(i, j, k, l, J0, J1) -> float:
    """VI_II: barχ(il;kj)^{J0,J1} = hΩ (−1)^{J0+J1} Σ occ Ω̄^{J1 J0} Γ̄^{J1}."""
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


def inv_tensor_plus(i, j, k, l, J0, J1) -> float:
    """AMC tensor inv WITHOUT leading minus (Path B gold)."""
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
            bc = bar_chi_VII(i, j, k, l, J2, J3)
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
            for J2 in range(
                abs(j2i(a) - j2i(i)) // 2, (j2i(a) + j2i(i)) // 2 + 1
            ):
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
                            max(abs(j2i(a) - j2i(c)), abs(j2i(b) - j2i(j)))
                            // 2,
                            min(j2i(a) + j2i(c), j2i(b) + j2i(j)) // 2 + 1,
                        ):
                            j0max = int(
                                max(J2, J3, J0, J1, jo(i), ja, lam) + 2
                            )
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
                                    J1
                                    + J2
                                    + J3
                                    + J4
                                    + (j2i(i) + j2i(c)) // 2
                                )
                                tot += ph * hats * six * w * om * gam
    return hat(J0) * hat(J1) * tot


print("Caching Path B χ^κ ...")
t0 = time.time()
chiJ: dict = {}
n_chi = 0
max_chi_vs_amc = 0.0
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
                        v = inv_tensor_plus(i, j, b, d, J0, J1)
                        if abs(v) > 1e-16:
                            chiJ[(J0, J1, i, j, b, d)] = v
                            n_chi += 1
                            # spot-check vs AMC direct
                            if n_chi <= 40:
                                va = chi_amc_direct(J0, J1, i, j, b, d)
                                max_chi_vs_amc = max(
                                    max_chi_vs_amc, abs(v - va)
                                )
print(
    f"  nonzero χ={len(chiJ)}  PathB↔AMC maxΔ={max_chi_vs_amc:.3e}"
    f"  ({time.time()-t0:.1f}s)"
)
if max_chi_vs_amc > tol:
    print(f"FAIL — Path B χ ≢ AMC direct (maxΔ={max_chi_vs_amc:.3e})")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Ladder RME: W = −χΩ  (DGEMM intermediate; no P yet)
# ---------------------------------------------------------------------------
def w_rme(J0: int, i: int, j: int, k: int, l: int) -> float:
    """W_red^{J0}_{ijkl} from Path B χ (G4a_Wbra_noperm)."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    pref = -iphase(J0) / hat(J0)
    tot = 0.0
    for b in orbits:
        for d in orbits:
            ph_bd = iphase((j2i(b) + j2i(d)) // 2 + lam)
            for J2 in range(0, max_J + 1):
                if not tri(J0, J2, lam):
                    continue
                if not tri(jo(b), jo(d), J2):
                    continue
                ch = chiJ.get((J0, J2, i, j, b, d), 0.0)
                if abs(ch) < 1e-16:
                    continue
                om = Eta.TwoBody.GetTBME_J(J2, J0, d, b, k, l)
                if abs(om) < 1e-16:
                    continue
                tot += ph_bd * hat_lam_inv * ch * om
    return pref * tot


def w_klij_rme(J0: int, i: int, j: int, k: int, l: int) -> float:
    """W_red^{J0}_{klij} — same kernel, externals swapped."""
    return w_rme(J0, k, l, i, j)


def z_from_W(J0: int, i: int, j: int, k: int, l: int) -> float:
    """(1−P) on W intermediates only → Z_red."""
    # bra: (1−P_ij) W_ijkl
    w1 = w_rme(J0, i, j, k, l)
    w2 = w_rme(J0, j, i, k, l)
    bra = w1 - iphase((j2i(i) + j2i(j)) // 2 - J0) * w2
    # ket: (1−P_kl) W_klij  (Hermitian twin — no extra h_Ω)
    wk1 = w_klij_rme(J0, i, j, k, l)
    wk2 = w_klij_rme(J0, i, j, l, k)
    ket = wk1 - iphase((j2i(k) + j2i(l)) // 2 - J0) * wk2
    return bra + ket


# ---------------------------------------------------------------------------
# m-scheme gold
# ---------------------------------------------------------------------------
def chi_we(i, mi, j, mj, b, mb, d, md) -> float:
    """WE unpack of Path B χ_J (≡ χ(m)). No Pauli skip — χ not AS."""
    if abs(mi + mj - mb - md) > 2 * lam:
        return 0.0
    M0 = (mi + mj) // 2
    M1 = (mb + md) // 2
    mu = M0 - M1
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0:
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        for J1 in range(abs(j2i(b) - j2i(d)) // 2, (j2i(b) + j2i(d)) // 2 + 1):
            if abs(M1) > J1 or not tri(J0, J1, lam):
                continue
            ccd = CG(jo(b), mb * 0.5, jo(d), md * 0.5, J1, M1)
            if abs(ccd) < 1e-15:
                continue
            cj = (
                1.0
                if (lam == 0 and J0 == J1 and M0 == M1)
                else CG(J1, M1, lam, mu, J0, M0)
            )
            if abs(cj) < 1e-15:
                continue
            v = chiJ.get((J0, J1, i, j, b, d), 0.0)
            if abs(v) < 1e-16:
                continue
            sm += cj * cab * ccd / hat(J0) * v
    return sm


def w_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """Bare W(m)=−Σ [χ×Ω]^(0); no P."""
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for b in orbits:
        for mb in mrange(b):
            for d in orbits:
                for md in mrange(d):
                    cm = chi_we(i, mi, j, mj, b, mb, d, md)
                    if abs(cm) < 1e-16:
                        continue
                    om = ut.GetMschemeMatrixElement_2b(
                        Eta, d, md, b, mb, k, mk, l, ml
                    )
                    if abs(om) < 1e-16:
                        continue
                    M_ij = (mi + mj) // 2
                    M_bd = (mb + md) // 2
                    mu = M_ij - M_bd
                    if abs(mu) > lam:
                        continue
                    cg0 = CG(lam, mu, lam, -mu, 0, 0)
                    if abs(cg0) < 1e-16:
                        continue
                    sm += cg0 * cm * om
    return -sm


def z_m(i, mi, j, mj, k, mk, l, ml) -> float:
    """(1−P_ij)W + (1−P_kl)W_klij in m — P on W only."""
    wm = w_m(i, mi, j, mj, k, mk, l, ml)
    wp = w_m(j, mj, i, mi, k, mk, l, ml)
    # W_klij(m): same contraction with legs (kl|ij)
    wk = w_m(k, mk, l, ml, i, mi, j, mj)
    wkp = w_m(l, ml, k, mk, i, mi, j, mj)
    return (wm - wp) + (wk - wkp)


def z_red_from_m(J0: int, i: int, j: int, k: int, l: int) -> float:
    """Z_red = Σ CG CG Z(m) / Ĵ  (AMC reduce=true packaging)."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    # Pauli: final Z is AS — skip i==j or k==l same-m only when projecting
    # AS operators; for compare to Z_red from (1−P)W we still sum all m
    # (W path already applied P). Use full sum without Pauli skip.
    sm = 0.0
    for mi in mrange(i):
        for mj in mrange(j):
            M = (mi + mj) // 2
            if abs(M) > J0:
                continue
            cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(cab) < 1e-15:
                continue
            for mk in mrange(k):
                ml = mi + mj - mk
                if abs(ml) > j2i(l) or ((ml + j2i(l)) % 2) != 0:
                    continue
                ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(ccd) < 1e-15:
                    continue
                zm = z_m(i, mi, j, mj, k, mk, l, ml)
                if abs(zm) < 1e-16:
                    continue
                sm += cab * ccd * zm
    return sm / hat(J0)


def w_red_from_m(J0: int, i: int, j: int, k: int, l: int) -> float:
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    sm = 0.0
    for mi in mrange(i):
        for mj in mrange(j):
            M = (mi + mj) // 2
            if abs(M) > J0:
                continue
            cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(cab) < 1e-15:
                continue
            for mk in mrange(k):
                ml = mi + mj - mk
                if abs(ml) > j2i(l) or ((ml + j2i(l)) % 2) != 0:
                    continue
                ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(ccd) < 1e-15:
                    continue
                wm = w_m(i, mi, j, mj, k, mk, l, ml)
                if abs(wm) < 1e-16:
                    continue
                sm += cab * ccd * wm
    return sm / hat(J0)


# ---------------------------------------------------------------------------
# Compares
# ---------------------------------------------------------------------------
def sample_compare(label, j_fn, m_fn):
    rats = Counter()
    n = 0
    max_abs = 0.0
    worst = None
    t0 = time.time()
    for i in orbits:
        for j in orbits:
            for k in orbits:
                for l in orbits:
                    for J0 in range(0, max_J + 1):
                        if not (
                            tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                        ):
                            continue
                        ja = j_fn(J0, i, j, k, l)
                        if abs(ja) < 1e-10:
                            continue
                        jm = m_fn(J0, i, j, k, l)
                        n += 1
                        err = abs(jm - ja)
                        if err > max_abs:
                            max_abs = err
                            r = jm / ja if abs(ja) > 1e-10 else float("nan")
                            worst = (J0, i, j, k, l, jm, ja, r)
                        rats[round(jm / ja, 4)] += 1
                        if n >= nsamp:
                            break
                    if n >= nsamp:
                        break
                if n >= nsamp:
                    break
            if n >= nsamp:
                break
        if n >= nsamp:
            break
    print(
        f"\n{label}\n"
        f"  n={n}  max|Δ|={max_abs:.3e}  m/J={rats.most_common(5)}"
        f"  ({time.time()-t0:.1f}s)"
    )
    if worst:
        J0, i, j, k, l, jm, ja, r = worst
        print(
            f"  worst J={J0} ({i},{j}|{k},{l}): m={jm:.6e} J={ja:.6e} r={r:.6g}"
        )
    return n > 0 and max_abs < tol, max_abs, rats


ok_W, dW, _ = sample_compare(
    "(A) bare W=−χΩ (no P): m→W_red ≡ Path B RME",
    w_rme,
    w_red_from_m,
)

ok_Z, dZ, ratsZ = sample_compare(
    "(B) full Z=(1−P_ij)W+(1−P_kl)W_klij: m→Z_red ≡ Path B",
    z_from_W,
    z_red_from_m,
)

print()
if ok_W and ok_Z:
    print(
        "PASS — Path B χ → ladder W → (1−P) on W ≡ m  "
        f"(Z_red packaging, λ={lam})"
    )
    sys.exit(0)
if ok_W:
    print(f"PARTIAL — bare W PASS; full Z maxΔ={dZ:.3e} m/J={ratsZ.most_common(3)}")
    sys.exit(1)
print(f"FAIL — bare W maxΔ={dW:.3e}")
sys.exit(1)
