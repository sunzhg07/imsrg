#!/usr/bin/env python3
"""Γ^{III_b} full Z — Path B pack (Pandya/RC + DGEMM) ≡ fold ≡ 4-index, any λ.

Path B χ̄ (same as term1/term2):
  λ=0:  Fac Pandya Ω̄ @ (occ⊙Ω̄)
  λ≠0:  AMC mid-J DGEMM → CHI = (−1)^J/Ĵ² · P[(a,b),(d,c)]

Then production pack:
  RC = code RC[χ̄_bc,ad + χ̄_ad,bc]
  W̄ = Γ̄ @ RC
  Z = InvPandya_noperm → (1−P_ij)(1−P_kl)

Gold: WE(Z) ≡ Zm_fold ≡ Z_4index
  (fold ≡ 4-index locked in test_G3b_pathB_fold_mscheme.py)

Usage:
  PYTHONPATH=build python3 -B run/test_G3b_pathB_pack_mscheme.py [emax=1] [lambda=2] [nsamp=6]
"""

from __future__ import annotations

import math
import random
import sys
import time
from collections import Counter

import numpy as np
from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 6
tol = 1e-5
seed = 11

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()
# λ=0 Fac: unreduced Ω. λ≠0 AMC Path B: reduced Ω.
if lam == 0:
    if Omega.IsReduced():
        Omega.MakeNotReduced()
else:
    if not Omega.IsReduced():
        Omega.MakeReduced()

hEta = -1 if Omega.IsAntiHermitian() else (1 if Omega.IsHermitian() else -1)
hGamma = 1 if Gamma.IsHermitian() else (-1 if Gamma.IsAntiHermitian() else 1)
orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)

print(
    f"emax={emax} λ_Ω={lam} hΩ={hEta} hΓ={hGamma} seed={seed}\n"
    f"  PathB pack: χ̄ DGEMM → RC(bcad+adbc) → Γ̄·RC → Inv → (1−P)² ≡ Zm_fold"
)


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


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


def w_eta(a, b, k) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def tbme5(Op, J, a, d, c, b) -> float:
    try:
        return Op.TwoBody.GetTBME_J(J, a, d, c, b)
    except TypeError:
        return Op.TwoBody.GetTBME_J(J, J, a, d, c, b)


def tbme_O(J1, J2, a, b, c, d) -> float:
    return Omega.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


# ---------------------------------------------------------------------------
# Pandya helpers
# ---------------------------------------------------------------------------
_fac_pandya: dict = {}
_amc_bar: dict = {}


def fac_pandya_Omega(a, b, c, d, Jcc) -> float:
    """Factorized scalar Pandya (legs adcb), λ=0 equal-J."""
    key = (a, b, c, d, Jcc)
    if key in _fac_pandya:
        return _fac_pandya[key]
    oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
    oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
    ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
    jc, jd = oc.j2 * 0.5, od.j2 * 0.5
    jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
    jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
    dJ = 1
    if a == d or b == c:
        dJ = 2
        jmin += jmin % 2
    X = 0.0
    for Jstd in range(jmin, jmax + 1, dJ):
        six = SixJ(ja, jb, Jcc, jc, jd, Jstd)
        if abs(six) > 1e-8:
            X -= (2 * Jstd + 1) * six * tbme5(Omega, Jstd, a, d, c, b)
    _fac_pandya[key] = X
    return X


def amc_bar_Omega(i, j, k, l, Jbra, Jket) -> float:
    """AMC / Path B Pandya of reduced Ω (any λ)."""
    key = (i, j, k, l, Jbra, Jket)
    if key in _amc_bar:
        return _amc_bar[key]
    if not tri(Jbra, Jket, lam):
        _amc_bar[key] = 0.0
        return 0.0
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    sm = 0.0
    for J2 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        for J3 in range(abs(j2i(k) - j2i(l)) // 2, (j2i(k) + j2i(l)) // 2 + 1):
            if not tri(J2, J3, lam):
                continue
            n9 = NineJ(lam, Jbra, Jket, J3, jl, jk, J2, ji, jj)
            if abs(n9) < 1e-16:
                continue
            sm += phase(J2) * hat(J2) * hat(J3) * n9 * tbme_O(J2, J3, i, j, k, l)
    v = -phase(Jbra + ji + jk + lam) * hat(Jbra) * hat(Jket) * sm
    _amc_bar[key] = v
    return v


def orient_pq(tbc, idx: int):
    """Oriented (p,q) for 2n index in a CC channel."""
    nK = tbc.GetNumberKets()
    if idx < nK:
        ket = tbc.GetKet(idx)
        return ket.p, ket.q
    ket = tbc.GetKet(idx - nK)
    return ket.q, ket.p


def orient_idx(tbc, p: int, q: int) -> int:
    nK = tbc.GetNumberKets()
    loc = tbc.GetLocalIndex(min(p, q), max(p, q))
    if loc < 0:
        return -1
    return loc + (nK if p > q else 0)


# ---------------------------------------------------------------------------
# Build Γ̄ + χ̄ (DGEMM)
# ---------------------------------------------------------------------------
t0 = time.time()
n_cc = ms.GetNumberTwoBodyChannels_CC()
bar_Gamma: list = []
barCHI: list = []


def fill_Gamma():
    for ch in range(n_cc):
        tbc = ms.GetTwoBodyChannel_CC(ch)
        nK = tbc.GetNumberKets()
        Jcc = tbc.J
        if nK < 1:
            bar_Gamma.append(None)
            continue
        n2 = 2 * nK
        bG = np.zeros((n2, n2))
        for ibra in range(nK):
            bra = tbc.GetKet(ibra)
            a, b = bra.p, bra.q
            oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
            ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
            for iket in range(n2):
                if (iket % nK) < ibra:
                    continue
                if iket < nK:
                    ket = tbc.GetKet(iket)
                    c, d = ket.p, ket.q
                else:
                    ket = tbc.GetKet(iket - nK)
                    d, c = ket.p, ket.q
                oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
                jc, jd = oc.j2 * 0.5, od.j2 * 0.5
                jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
                jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
                dJ = 1
                if a == d or b == c:
                    dJ = 2
                    jmin += jmin % 2
                Gammabar = 0.0
                for Jstd in range(jmin, jmax + 1, dJ):
                    six = SixJ(ja, jb, Jcc, jc, jd, Jstd)
                    if abs(six) > 1e-8:
                        Gammabar -= (2 * Jstd + 1) * six * tbme5(
                            Gamma, Jstd, a, d, c, b
                        )
                flip = iphase((oa.j2 + ob.j2 + oc.j2 + od.j2) // 2)
                if iket < nK or (iket >= nK and c != d):
                    bG[ibra, iket] = Gammabar
                    if iket != ibra:
                        bG[iket, ibra] = hGamma * Gammabar
                if a != b:
                    bG[ibra + nK, (iket + nK) % n2] = Gammabar * flip * hGamma
                if iket >= nK or (iket < nK and c != d):
                    bG[(iket + nK) % n2, ibra + nK] = Gammabar * flip
        bar_Gamma.append(bG)


def fill_chi_equalJ_fac():
    """λ=0: χ̄ = Ω̄ @ (occ⊙Ω̄) Factorized DGEMM."""
    for ch in range(n_cc):
        tbc = ms.GetTwoBodyChannel_CC(ch)
        nK = tbc.GetNumberKets()
        Jcc = tbc.J
        if nK < 1:
            barCHI.append(None)
            continue
        n2 = 2 * nK
        bE = np.zeros((n2, n2))
        nnn = np.zeros((n2, n2))
        for ibra in range(nK):
            bra = tbc.GetKet(ibra)
            a, b = bra.p, bra.q
            oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
            na, nb = oa.occ, ob.occ
            nba, nbb = 1.0 - na, 1.0 - nb
            for iket in range(n2):
                if (iket % nK) < ibra:
                    continue
                if iket < nK:
                    ket = tbc.GetKet(iket)
                    c, d = ket.p, ket.q
                else:
                    ket = tbc.GetKet(iket - nK)
                    d, c = ket.p, ket.q
                oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
                nc, nd = oc.occ, od.occ
                nbc, nbd = 1.0 - nc, 1.0 - nd
                occ_AbarBC = nba * nb * nc + na * nbb * nbc
                occ_ABbarD = na * nbb * nd + nba * nb * nbd
                occ_BCDbar = nb * nc * nbd + nbb * nbc * nd
                occ_ACbarD = na * nbc * nd + nba * nc * nbd
                Etabar = fac_pandya_Omega(a, b, c, d, Jcc)
                flip = iphase((oa.j2 + ob.j2 + oc.j2 + od.j2) // 2)
                if iket < nK or (iket >= nK and c != d):
                    bE[ibra, iket] = Etabar
                    nnn[ibra, iket] = Etabar * occ_AbarBC
                    if iket != ibra:
                        bE[iket, ibra] = hEta * Etabar
                        nnn[iket, ibra] = hEta * Etabar * occ_ACbarD
                if a != b:
                    bE[ibra + nK, (iket + nK) % n2] = Etabar * flip * hEta
                    nnn[ibra + nK, (iket + nK) % n2] = (
                        Etabar * flip * hEta * occ_ABbarD
                    )
                if iket >= nK or (iket < nK and c != d):
                    bE[(iket + nK) % n2, ibra + nK] = Etabar * flip
                    nnn[(iket + nK) % n2, ibra + nK] = (
                        Etabar * flip * occ_BCDbar
                    )
        barCHI.append(bE @ nnn)


def fill_chi_midJ_amc():
    """λ≠0: AMC-bar mid-J DGEMM → CHI ≡ mid(adcb)/Ĵ."""
    hatlam = hat(lam)
    shared = [(a, b) for a in orbits for b in orbits]
    ns = len(shared)

    for ch in range(n_cc):
        tbc = ms.GetTwoBodyChannel_CC(ch)
        nK = tbc.GetNumberKets()
        if nK < 1:
            barCHI.append(None)
            continue
        J0 = tbc.J
        n2 = 2 * nK
        P = np.zeros((n2, n2))
        for ch2 in range(n_cc):
            t2 = ms.GetTwoBodyChannel_CC(ch2)
            if t2.GetNumberKets() < 1:
                continue
            J2 = t2.J
            if not tri(J0, J2, lam):
                continue
            if tbc.parity != t2.parity or tbc.Tz != t2.Tz:
                continue
            L = np.zeros((n2, ns))
            R = np.zeros((ns, n2))
            for ib in range(n2):
                i, l = orient_pq(tbc, ib)
                for s, (a, b) in enumerate(shared):
                    L[ib, s] = amc_bar_Omega(i, a, b, l, J0, J2)
            for s, (a, b) in enumerate(shared):
                for jk in range(n2):
                    j, k = orient_pq(tbc, jk)
                    R[s, jk] = w_eta(a, b, k) * amc_bar_Omega(
                        b, j, k, a, J2, J0
                    )
            P += (phase(J2 + lam) / hatlam) * (L @ R)

        # CHI[a,b;c,d] = phase(J0)/Ĵ² · P[(a,b),(d,c)] ≡ mid(a,d,c,b)/Ĵ
        CHI = np.zeros((n2, n2))
        sc = phase(J0) / (hat(J0) ** 2)
        for ib in range(n2):
            for ik in range(n2):
                c, d = orient_pq(tbc, ik)
                jk = orient_idx(tbc, d, c)
                if jk < 0:
                    continue
                CHI[ib, ik] = sc * P[ib, jk]
        barCHI.append(CHI)


fill_Gamma()
if lam == 0:
    fill_chi_equalJ_fac()
else:
    fill_chi_midJ_amc()

print(f"  Pandya Γ + χ̄ DGEMM done ({time.time()-t0:.1f}s)")


# ---------------------------------------------------------------------------
# RC₂ (term2: χ̄_ad,bc) + GEMM; also RC₁ for pack sanity
# ---------------------------------------------------------------------------
def build_RC(which: str) -> list:
    """which='bcad' → term1 RC₁; which='adbc' → term2 RC₂."""
    out = []
    for ch in range(n_cc):
        tbc = ms.GetTwoBodyChannel_CC(ch)
        nK = tbc.GetNumberKets()
        Jcc = tbc.J
        if nK < 1 or barCHI[ch] is None:
            out.append(None)
            continue
        n2 = 2 * nK
        RC = np.zeros((n2, n2))
        for ibra in range(n2):
            a, b = orient_pq(tbc, ibra)
            if ibra >= nK and a == b:
                continue
            oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
            ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
            for iket in range(n2):
                c, d = orient_pq(tbc, iket)
                if iket >= nK and c == d:
                    continue
                oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
                jc, jd = oc.j2 * 0.5, od.j2 * 0.5
                jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
                jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
                X = 0.0
                for Jp in range(jmin, jmax + 1):
                    six = SixJ(ja, jb, Jcc, jc, jd, Jp)
                    if abs(six) < 1e-8:
                        continue
                    parity = (oa.l + od.l) % 2
                    Tz = abs(oa.tz2 - od.tz2) // 2
                    ch_old = ms.GetTwoBodyChannelIndex(Jp, parity, Tz)
                    if ch_old < 0 or barCHI[ch_old] is None:
                        continue
                    t_old = ms.GetTwoBodyChannel_CC(ch_old)
                    if t_old.J != Jp or t_old.parity != parity or t_old.Tz != Tz:
                        continue
                    nk = t_old.GetNumberKets()
                    if nk < 1:
                        continue
                    indx_ad = t_old.GetLocalIndex(min(a, d), max(a, d))
                    indx_bc = t_old.GetLocalIndex(min(b, c), max(b, c))
                    if indx_ad < 0 or indx_bc < 0:
                        continue
                    if a > d:
                        indx_ad += nk
                    if b > c:
                        indx_bc += nk
                    Mchi = barCHI[ch_old]
                    if which == "bcad":
                        me = Mchi[indx_bc, indx_ad]
                    elif which == "adbc":
                        me = Mchi[indx_ad, indx_bc]
                    else:
                        raise ValueError(which)
                    X -= (
                        iphase((ob.j2 + oc.j2) // 2 + Jp)
                        * (2 * Jp + 1)
                        * six
                        * me
                    )
                RC[ibra, iket] = X
        out.append(RC)
    return out


RC2 = build_RC("adbc")
RC1 = build_RC("bcad")
Wmats = [
    None if RC2[ch] is None else bar_Gamma[ch] @ RC2[ch] for ch in range(n_cc)
]
Wmats1 = [
    None if RC1[ch] is None else bar_Gamma[ch] @ RC1[ch] for ch in range(n_cc)
]
print(f"  RC₂ + GEMM done ({time.time()-t0:.1f}s)")


def readW_from(mats, a, b, c, d, J, parity, Tz) -> float:
    ch = ms.GetTwoBodyChannelIndex(J, parity, Tz)
    if ch < 0 or mats[ch] is None:
        return 0.0
    tbc = ms.GetTwoBodyChannel_CC(ch)
    if tbc.J != J or tbc.parity != parity or tbc.Tz != Tz:
        return 0.0
    nK = tbc.GetNumberKets()
    if nK < 1:
        return 0.0
    iab = tbc.GetLocalIndex(min(a, b), max(a, b))
    icd = tbc.GetLocalIndex(min(c, d), max(c, d))
    if iab < 0 or icd < 0:
        return 0.0
    iab += nK if a > b else 0
    icd += nK if c > d else 0
    M = mats[ch]
    if iab >= M.shape[0] or icd >= M.shape[1]:
        return 0.0
    return float(M[iab, icd])


def readW(a, b, c, d, J, parity, Tz) -> float:
    return readW_from(Wmats, a, b, c, d, J, parity, Tz)


def z_noperm_from(mats, J0, i, j, k, l) -> float:
    oi, oj = ms.GetOrbit(i), ms.GetOrbit(j)
    ok, ol = ms.GetOrbit(k), ms.GetOrbit(l)
    ji, jj, jk, jl = oi.j2, oj.j2, ok.j2, ol.j2
    if not (tri(ji * 0.5, jj * 0.5, J0) and tri(jk * 0.5, jl * 0.5, J0)):
        return 0.0
    sm = 0.0
    parity = (oi.l + ok.l) % 2
    Tz = abs(oi.tz2 - ok.tz2) // 2
    for Jp in range(
        max(abs(jj - jl), abs(ji - jk)) // 2,
        min(jj + jl, ji + jk) // 2 + 1,
    ):
        six = SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jp)
        if abs(six) < 1e-8:
            continue
        me = readW_from(mats, j, l, i, k, Jp, parity, Tz)
        sm -= iphase(Jp + (ji + jk) // 2) * (2 * Jp + 1) * six * me
    return iphase(J0 + (ji + jj) // 2) * sm


def z_noperm(J0, i, j, k, l) -> float:
    return z_noperm_from(Wmats, J0, i, j, k, l)


def z_we_from(mats, i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    M = (mi + mj) // 2
    sm = 0.0
    for J0 in range(0, max_J + 1):
        if abs(M) > J0:
            continue
        if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
        ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
        if abs(cab * ccd) < 1e-20:
            continue
        sm += cab * ccd * z_noperm_from(mats, J0, i, j, k, l)
    return sm


def z_we(i, mi, j, mj, k, mk, l, ml) -> float:
    return z_we_from(Wmats, i, mi, j, mj, k, mk, l, ml)


# ---------------------------------------------------------------------------
# Gold W2_m  (−Σ χ_lajb Γ_aibk)
# ---------------------------------------------------------------------------
Omega_m = Omega
if lam != 0 and not Omega.IsReduced():
    Omega_m = Operator(Omega)
    Omega_m.MakeReduced()

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
                    w = w_eta(a, b, k)
                    if abs(w) < 1e-12:
                        continue
                    o1 = ut.GetMschemeMatrixElement_2b(
                        Omega_m, i, mi, a, ma, b, mb, l, ml
                    )
                    if abs(o1) < 1e-16:
                        continue
                    mu = 0.5 * (mi + ma - mb - ml)
                    if abs(mu) > lam:
                        continue
                    cg = CG(lam, mu, lam, -mu, 0, 0)
                    if abs(cg) < 1e-16:
                        continue
                    o2 = ut.GetMschemeMatrixElement_2b(
                        Omega_m, b, mb, j, mj, k, mk, a, ma
                    )
                    sm += cg * w * o1 * o2
    _chi_m[key] = sm
    return sm


def W2_m(i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    c = chi_m(l, ml, a, ma, j, mj, b, mb)
                    if abs(c) < 1e-16:
                        continue
                    g = ut.GetMschemeMatrixElement_2b(
                        Gamma, a, ma, i, mi, b, mb, k, mk
                    )
                    sm -= c * g
    return sm


def W1_m(i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    c = chi_m(b, mb, k, mk, a, ma, i, mi)
                    if abs(c) < 1e-16:
                        continue
                    g = ut.GetMschemeMatrixElement_2b(
                        Gamma, j, mj, b, mb, l, ml, a, ma
                    )
                    sm -= c * g
    return sm


# Pack sanity: RC₁ + RC₂ ≡ RC(χ̄_bc,ad + χ̄_ad,bc)
max_pack_d = 0.0
n_pack = 0
RC_pack = []
for ch in range(n_cc):
    if RC1[ch] is None:
        RC_pack.append(None)
        continue
    RC_pack.append(RC1[ch] + RC2[ch])


def build_RC_pack() -> list:
    out = []
    for ch in range(n_cc):
        tbc = ms.GetTwoBodyChannel_CC(ch)
        nK = tbc.GetNumberKets()
        Jcc = tbc.J
        if nK < 1 or barCHI[ch] is None:
            out.append(None)
            continue
        n2 = 2 * nK
        RC = np.zeros((n2, n2))
        for ibra in range(n2):
            a, b = orient_pq(tbc, ibra)
            if ibra >= nK and a == b:
                continue
            oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
            ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
            for iket in range(n2):
                c, d = orient_pq(tbc, iket)
                if iket >= nK and c == d:
                    continue
                oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
                jc, jd = oc.j2 * 0.5, od.j2 * 0.5
                jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
                jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
                X = 0.0
                for Jp in range(jmin, jmax + 1):
                    six = SixJ(ja, jb, Jcc, jc, jd, Jp)
                    if abs(six) < 1e-8:
                        continue
                    parity = (oa.l + od.l) % 2
                    Tz = abs(oa.tz2 - od.tz2) // 2
                    ch_old = ms.GetTwoBodyChannelIndex(Jp, parity, Tz)
                    if ch_old < 0 or barCHI[ch_old] is None:
                        continue
                    t_old = ms.GetTwoBodyChannel_CC(ch_old)
                    if t_old.J != Jp or t_old.parity != parity or t_old.Tz != Tz:
                        continue
                    nk = t_old.GetNumberKets()
                    if nk < 1:
                        continue
                    indx_ad = t_old.GetLocalIndex(min(a, d), max(a, d))
                    indx_bc = t_old.GetLocalIndex(min(b, c), max(b, c))
                    if indx_ad < 0 or indx_bc < 0:
                        continue
                    if a > d:
                        indx_ad += nk
                    if b > c:
                        indx_bc += nk
                    Mchi = barCHI[ch_old]
                    me = Mchi[indx_bc, indx_ad] + Mchi[indx_ad, indx_bc]
                    X -= (
                        iphase((ob.j2 + oc.j2) // 2 + Jp)
                        * (2 * Jp + 1)
                        * six
                        * me
                    )
                RC[ibra, iket] = X
        out.append(RC)
    return out


print("  building pack RC for sanity…", flush=True)
RCpack_exp = build_RC_pack()
for ch in range(n_cc):
    if RC_pack[ch] is None or RCpack_exp[ch] is None:
        continue
    d = float(np.max(np.abs(RC_pack[ch] - RCpack_exp[ch])))
    max_pack_d = max(max_pack_d, d)
    n_pack += 1
print(f"  pack sanity RC₁+RC₂ vs RC(χ̄_bc+χ̄_ad): maxΔ={max_pack_d:.3e} ({n_pack} ch)")

Wmats_pack = [
    None
    if RC_pack[ch] is None or bar_Gamma[ch] is None
    else bar_Gamma[ch] @ RC_pack[ch]
    for ch in range(n_cc)
]


def antisym_we(fn, i, mi, j, mj, k, mk, l, ml):
    w = fn(i, mi, j, mj, k, mk, l, ml)
    w -= fn(j, mj, i, mi, k, mk, l, ml)
    w -= fn(i, mi, j, mj, l, ml, k, mk)
    w += fn(j, mj, i, mi, l, ml, k, mk)
    return w


def Wm_fold(i, mi, j, mj, k, mk, l, ml) -> float:
    return W1_m(i, mi, j, mj, k, mk, l, ml) + W2_m(i, mi, j, mj, k, mk, l, ml)


def Zm_fold(i, mi, j, mj, k, mk, l, ml) -> float:
    return antisym_we(Wm_fold, i, mi, j, mj, k, mk, l, ml)


def Z_pack(i, mi, j, mj, k, mk, l, ml) -> float:
    return antisym_we(lambda *x: z_we_from(Wmats_pack, *x), i, mi, j, mj, k, mk, l, ml)


random.seed(seed)
cands = []
for _ in range(8000):
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

hits = []
t1 = time.time()
for tup in cands:
    if len(hits) >= nsamp:
        break
    zf = Zm_fold(*tup)
    if abs(zf) < 1e-5:
        continue
    zp = Z_pack(*tup)
    hits.append((zf, zp, tup))
    print(
        f"  hit {len(hits)}/{nsamp} |Zfold|={abs(zf):.3e} ({time.time()-t1:.0f}s)",
        flush=True,
    )

if len(hits) < max(3, nsamp // 4):
    print(f"FAIL — too few hits ({len(hits)})")
    sys.exit(1)

rats: Counter = Counter()
n1 = 0
maxd = 0.0
for zf, zp, tup in hits:
    r = zp / zf if abs(zf) > 1e-12 else float("nan")
    rats[round(r, 4)] += 1
    maxd = max(maxd, abs(zp - zf))
    if abs(r - 1.0) < 1e-3:
        n1 += 1
    i, mi, j, mj, k, mk, l, ml = tup
    print(
        f"  ({i},{mi})({j},{mj})({k},{mk})({l},{ml}) "
        f"fold={zf:.5e} pack={zp:.5e} r={r:.5g}"
    )

print(f"n1={n1}/{len(hits)} rats={rats.most_common(5)} maxΔ={maxd:.3e}")
print(f"RC pack linearity maxΔ={max_pack_d:.3e}")
if max_pack_d >= tol:
    print("FAIL — RC pack linearity broken")
    sys.exit(1)
if maxd < tol and n1 == len(hits):
    print(
        f"PASS — PathB pack(χ̄ DGEMM→RC→Γ̄→Inv→(1−P)²) ≡ Zm_fold  (λ_Ω={lam})"
    )
    sys.exit(0)
print("FAIL — Path B pack not locked to fold")
sys.exit(1)
