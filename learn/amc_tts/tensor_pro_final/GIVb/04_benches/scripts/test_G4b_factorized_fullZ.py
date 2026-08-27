#!/usr/bin/env python3
"""Γ^{IV_b} Factorized Path B full Z ≡ m (λ=0).

Faithful Python port of FactorizedDoubleCommutator.cc χ^ι branch:
  1) Scalar Pandya → 2n×2n CC mats (bar_Eta, bar_Gamma, nnnbar_Eta)
  2) bar_CHI_V = bar_Gamma @ nnnbar_Eta
  3) Factorized RC (pack χ̄_adbc − h_Z χ̄_bcad) → bar_CHI_V_RC
  4) CHI_V_final = bar_Eta @ bar_CHI_V_RC
  5) Inv Pandya with (1−P)^2 only (L1884–2013)
  6) Compare magnetic MEs to m-scheme fold

Gold: m-scheme analyze fold (no TTS).
  χ^ι_ijkl = Σ_ab (n̄_a n_b n̄_k + n_a n̄_b n_k) Ω_bjka Γ_iabl
  Γ^{IV_b} = (1−P_ij)(1−P_kl) Σ_ab (χ_aibk Ω_jbla − χ_akbi Ω_jalb)

Usage:
  PYTHONPATH=build python3 -B run/test_G4b_factorized_fullZ.py [emax=1] [nsamp=20]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

import numpy as np
from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 20
tol = 1e-5
seed = 11
lam = 0  # scalar Factorized dual

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

hEta = -1 if Eta.IsAntiHermitian() else (1 if Eta.IsHermitian() else -1)
hGamma = 1 if Gamma.IsHermitian() else (-1 if Gamma.IsAntiHermitian() else 1)
hZ = hGamma
orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)

print(
    f"emax={emax} λ={lam} seed={seed} hEta={hEta} hGamma={hGamma} hZ={hZ}\n"
    f"  Factorized: Pandya→RC→Ω̄·RC→InvPandya(1−P)^2 ≡ m\n"
    f"  gold: m-scheme IV_b fold (no TTS)"
)


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


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


def tbme5(Op, J, a, d, c, b) -> float:
    try:
        return Op.TwoBody.GetTBME_J(J, a, d, c, b)
    except TypeError:
        return Op.TwoBody.GetTBME_J(J, J, a, d, c, b)


# ---------------------------------------------------------------------------
# Factorized CC mats (L1210–1350, L1756–1869)
# ---------------------------------------------------------------------------
t0 = time.time()
n_cc = ms.GetNumberTwoBodyChannels_CC()
bar_Eta: list = []
bar_CHI_V: list = []
CHI_V_final: list = []

for ch in range(n_cc):
    tbc = ms.GetTwoBodyChannel_CC(ch)
    nK = tbc.GetNumberKets()
    Jcc = tbc.J
    if nK < 1:
        bar_Eta.append(None)
        bar_CHI_V.append(None)
        continue
    n2 = 2 * nK
    bE = np.zeros((n2, n2))
    bG = np.zeros((n2, n2))
    nnn = np.zeros((n2, n2))
    for ibra in range(nK):
        bra = tbc.GetKet(ibra)
        a, b = bra.p, bra.q
        oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
        ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
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
            jc, jd = oc.j2 * 0.5, od.j2 * 0.5
            nc, nd = oc.occ, od.occ
            nbc, nbd = 1.0 - nc, 1.0 - nd
            occ_AbarBC = nba * nb * nc + na * nbb * nbc
            occ_ABbarD = na * nbb * nd + nba * nb * nbd
            occ_BCDbar = nb * nc * nbd + nbb * nbc * nd
            occ_ACbarD = na * nbc * nd + nba * nc * nbd
            jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
            jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
            dJ = 1
            if a == d or b == c:
                dJ = 2
                jmin += jmin % 2
            Etabar = Gammabar = 0.0
            for Jstd in range(jmin, jmax + 1, dJ):
                six = SixJ(ja, jb, Jcc, jc, jd, Jstd)
                if abs(six) > 1e-8:
                    Etabar -= (2 * Jstd + 1) * six * tbme5(Eta, Jstd, a, d, c, b)
                    Gammabar -= (2 * Jstd + 1) * six * tbme5(
                        Gamma, Jstd, a, d, c, b
                    )
            flip = iphase((oa.j2 + ob.j2 + oc.j2 + od.j2) // 2)
            if iket < nK or (iket >= nK and c != d):
                bG[ibra, iket] = Gammabar
                bE[ibra, iket] = Etabar
                nnn[ibra, iket] = Etabar * occ_AbarBC
                if iket != ibra:
                    bG[iket, ibra] = hGamma * Gammabar
                    bE[iket, ibra] = hEta * Etabar
                    nnn[iket, ibra] = hEta * Etabar * occ_ACbarD
            if a != b:
                bG[ibra + nK, (iket + nK) % n2] = Gammabar * flip * hGamma
                bE[ibra + nK, (iket + nK) % n2] = Etabar * flip * hEta
                nnn[ibra + nK, (iket + nK) % n2] = (
                    Etabar * flip * hEta * occ_ABbarD
                )
            if iket >= nK or (iket < nK and c != d):
                bG[(iket + nK) % n2, ibra + nK] = Gammabar * flip
                bE[(iket + nK) % n2, ibra + nK] = Etabar * flip
                nnn[(iket + nK) % n2, ibra + nK] = Etabar * flip * occ_BCDbar
    bar_Eta.append(bE)
    bar_CHI_V.append(bG @ nnn)

for ch in range(n_cc):
    tbc = ms.GetTwoBodyChannel_CC(ch)
    nK = tbc.GetNumberKets()
    Jcc = tbc.J
    if nK < 1 or bar_CHI_V[ch] is None:
        CHI_V_final.append(None)
        continue
    n2 = 2 * nK
    RC = np.zeros((n2, n2))
    for ibra in range(n2):
        if ibra < nK:
            bra = tbc.GetKet(ibra)
            a, b = bra.p, bra.q
        else:
            bra = tbc.GetKet(ibra - nK)
            b, a = bra.p, bra.q
        if ibra >= nK and a == b:
            continue
        oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
        ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
        for iket in range(n2):
            if iket < nK:
                ket = tbc.GetKet(iket)
                c, d = ket.p, ket.q
            else:
                ket = tbc.GetKet(iket - nK)
                d, c = ket.p, ket.q
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
                if ch_old < 0 or bar_CHI_V[ch_old] is None:
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
                mat = bar_CHI_V[ch_old]
                me = mat[indx_ad, indx_bc] - hZ * mat[indx_bc, indx_ad]
                X += (
                    iphase((ob.j2 + oc.j2) // 2 + Jp)
                    * (2 * Jp + 1)
                    * six
                    * me
                )
            RC[ibra, iket] = X
    CHI_V_final.append(bar_Eta[ch] @ RC)

print(f"  built CC mats in {time.time() - t0:.1f}s")


def chi_final_read(a, b, c, d, J, parity, Tz) -> float:
    ch = ms.GetTwoBodyChannelIndex(J, parity, Tz)
    if ch < 0 or CHI_V_final[ch] is None:
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
    M = CHI_V_final[ch]
    if iab >= M.shape[0] or icd >= M.shape[1]:
        return 0.0
    return float(M[iab, icd])


def z_pb(J0, i, j, k, l) -> float:
    """Factorized InvPandya of CHI_V_final (L1884–2013)."""
    oi, oj = ms.GetOrbit(i), ms.GetOrbit(j)
    ok, ol = ms.GetOrbit(k), ms.GetOrbit(l)
    ji, jj, jk, jl = oi.j2, oj.j2, ok.j2, ol.j2
    if not (tri(ji * 0.5, jj * 0.5, J0) and tri(jk * 0.5, jl * 0.5, J0)):
        return 0.0
    cijkl = cjikl = cijlk = cjilk = 0.0
    parity = (oi.l + ok.l) % 2
    Tz = abs(oi.tz2 - ok.tz2) // 2
    Jpmin = max(abs(jj - jl), abs(ji - jk)) // 2
    Jpmax = min(jj + jl, ji + jk) // 2
    for Jp in range(Jpmin, Jpmax + 1):
        six = SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jp)
        if abs(six) < 1e-8:
            continue
        me1 = chi_final_read(j, l, i, k, Jp, parity, Tz)
        me2 = chi_final_read(i, k, j, l, Jp, parity, Tz)
        cjikl -= iphase(Jp + (ji + jk) // 2) * (2 * Jp + 1) * six * me1
        cijlk -= iphase(Jp + (jj + jl) // 2) * (2 * Jp + 1) * six * me2
    parity = (oi.l + ol.l) % 2
    Tz = abs(oi.tz2 - ol.tz2) // 2
    Jpmin = max(abs(ji - jl), abs(jj - jk)) // 2
    Jpmax = min(ji + jl, jj + jk) // 2
    for Jp in range(Jpmin, Jpmax + 1):
        six = SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jp)
        if abs(six) < 1e-8:
            continue
        me1 = chi_final_read(j, k, i, l, Jp, parity, Tz)
        me2 = chi_final_read(i, l, j, k, Jp, parity, Tz)
        cjilk -= iphase(Jp + (ji + jl) // 2) * (2 * Jp + 1) * six * me1
        cijkl -= iphase(Jp + (jj + jk) // 2) * (2 * Jp + 1) * six * me2
    z = cjikl - iphase((ji + jj) // 2 - J0) * cijkl
    z += -iphase((jl + jk) // 2 - J0) * cjilk + iphase(
        (jk + jl + ji + jj) // 2
    ) * cijlk
    if i == j:
        z /= math.sqrt(2.0)
    if k == l:
        z /= math.sqrt(2.0)
    return iphase(J0 + (ji + jj) // 2) * z


def z_we(i, mi, j, mj, k, mk, l, ml) -> float:
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
        if abs(cab) < 1e-15:
            continue
        ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
        if abs(ccd) < 1e-15:
            continue
        zj = z_pb(J0, i, j, k, l)
        if abs(zj) < 1e-16:
            continue
        sm += cab * ccd * zj
    return sm


# ---------------------------------------------------------------------------
# m-scheme gold
# ---------------------------------------------------------------------------
_chi: dict = {}


def chi_m(i, mi, j, mj, k, mk, l, ml) -> float:
    key = (i, mi, j, mj, k, mk, l, ml)
    if key in _chi:
        return _chi[key]
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    w = nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)
                    if abs(w) < 1e-12:
                        continue
                    g = ut.GetMschemeMatrixElement_2b(
                        Gamma, i, mi, a, ma, b, mb, l, ml
                    )
                    if abs(g) < 1e-16:
                        continue
                    o = ut.GetMschemeMatrixElement_2b(
                        Eta, b, mb, j, mj, k, mk, a, ma
                    )
                    sm += w * g * o
    _chi[key] = sm
    return sm


def Wm(i, mi, j, mj, k, mk, l, ml) -> float:
    if (mi + mj) != (mk + ml):
        return 0.0
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    c1 = chi_m(a, ma, i, mi, b, mb, k, mk)
                    if abs(c1) > 1e-16:
                        o1 = ut.GetMschemeMatrixElement_2b(
                            Eta, j, mj, b, mb, l, ml, a, ma
                        )
                        if abs(o1) > 1e-16:
                            sm += c1 * o1
                    c2 = chi_m(a, ma, k, mk, b, mb, i, mi)
                    if abs(c2) > 1e-16:
                        o2 = ut.GetMschemeMatrixElement_2b(
                            Eta, j, mj, a, ma, l, ml, b, mb
                        )
                        if abs(o2) > 1e-16:
                            sm -= c2 * o2
    return sm


def Zm(i, mi, j, mj, k, mk, l, ml) -> float:
    w = Wm(i, mi, j, mj, k, mk, l, ml)
    w -= Wm(j, mj, i, mi, k, mk, l, ml)
    w -= Wm(i, mi, j, mj, l, ml, k, mk)
    w += Wm(j, mj, i, mi, l, ml, k, mk)
    return w


t1 = time.time()
cands = []
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for mi in mrange(i):
                    for mj in mrange(j):
                        Mtot = mi + mj
                        for mk in mrange(k):
                            ml = Mtot - mk
                            if ml not in mrange(l):
                                continue
                            if i == j and mi == mj:
                                continue
                            if k == l and mk == ml:
                                continue
                            zm = Zm(i, mi, j, mj, k, mk, l, ml)
                            if abs(zm) > 1e-6:
                                zp = z_we(i, mi, j, mj, k, mk, l, ml)
                                cands.append(
                                    (abs(zm), zm, zp, i, mi, j, mj, k, mk, l, ml)
                                )
cands.sort(reverse=True)
print(f"  compared {len(cands)} MEs in {time.time() - t1:.1f}s")

rats: Counter = Counter()
maxd = 0.0
n1 = 0
for _, zm, zp, i, mi, j, mj, k, mk, l, ml in cands[:nsamp]:
    r = zp / zm if abs(zm) > 1e-12 else float("nan")
    rats[round(r, 4)] += 1
    maxd = max(maxd, abs(zp - zm))
    if abs(r - 1.0) < 1e-3:
        n1 += 1
    print(
        f"  ({i},{mi})({j},{mj})({k},{mk})({l},{ml}) "
        f"m={zm:.5e} Fac={zp:.5e} r={r:.5g}"
    )

print(f"n1={n1}/{nsamp} rats={rats.most_common(6)} maxΔ={maxd:.3e}")
if maxd < tol and n1 == nsamp:
    print("PASS")
    sys.exit(0)
print("FAIL")
sys.exit(1)
