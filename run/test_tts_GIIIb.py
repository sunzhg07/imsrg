#!/usr/bin/env python3
"""Γ^{III_b} ethS Path B vs Python Path B (z_pb).

ethS: locked χ^η (normal, not AS) → Fac Pandya (no hermiticity) → code RC
      → Γ̄ DGEMM → Inv + (1−P_ij)(1−P_kl)
Gold: test_G3b_normal_to_RC.py (Python Path B ≡ fold ≡ m). Do not use tts_GIIIb.

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIIIb.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

import numpy as np
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

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if not Omega.IsReduced():
    Omega.MakeReduced()
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

hGamma = 1 if Gamma.IsHermitian() else (-1 if Gamma.IsAntiHermitian() else 1)
orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)

print(
    f"emax={emax} λ_Ω={lam} seed={seed}\n"
    f"  ethS GIIIb Path B vs Python Path B / fold"
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


def tbme_O(J1, J2, a, b, c, d) -> float:
    return Omega.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


def tbme5(Op, J, a, d, c, b) -> float:
    try:
        return Op.TwoBody.GetTBME_J(J, a, d, c, b)
    except TypeError:
        return Op.TwoBody.GetTBME_J(J, J, a, d, c, b)


_barO: dict = {}
_barChi: dict = {}
_chiB: dict = {}


def pandya_bar_amc(i, j, k, l, Jbra, Jket) -> float:
    key = (i, j, k, l, Jbra, Jket)
    if key in _barO:
        return _barO[key]
    if not tri(Jbra, Jket, lam):
        _barO[key] = 0.0
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
    _barO[key] = v
    return v


def bar_chi_mid(i, j, k, l, J0) -> float:
    key = (i, j, k, l, J0)
    if key in _barChi:
        return _barChi[key]
    sm = 0.0
    maxJ = max_J + lam + 2
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
                    * pandya_bar_amc(i, a, b, l, J0, J2)
                    * pandya_bar_amc(b, j, k, a, J2, J0)
                )
    v = phase(J0) / hat(J0) * sm
    _barChi[key] = v
    return v


def chi_path_B_red(i, j, k, l, J0) -> float:
    key = (i, j, k, l, J0)
    if key in _chiB:
        return _chiB[key]
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not (tri(ji, jj, J0) and tri(jk, jl, J0)):
        _chiB[key] = 0.0
        return 0.0
    sm = 0.0
    Jpmax = max((j2i(i) + j2i(j)) // 2, (j2i(k) + j2i(l)) // 2) + 1
    for Jp in range(0, Jpmax + 1):
        six = SixJ(jl, jk, J0, jj, ji, Jp)
        if abs(six) < 1e-16:
            continue
        sm += hat(Jp) * six * bar_chi_mid(i, j, k, l, Jp)
    v = hat(J0) * sm
    _chiB[key] = v
    return v


def chi_unnorm(J, i, j, k, l) -> float:
    return hat(J) * chi_path_B_red(i, j, k, l, J)


t0 = time.time()
n_cc = ms.GetNumberTwoBodyChannels_CC()
bar_Gamma: list = []
barCHI: list = []
CHI_final: list = []

for ch in range(n_cc):
    tbc = ms.GetTwoBodyChannel_CC(ch)
    nK = tbc.GetNumberKets()
    Jcc = tbc.J
    if nK < 1:
        bar_Gamma.append(None)
        barCHI.append(None)
        continue
    n2 = 2 * nK
    bG = np.zeros((n2, n2))
    bC = np.zeros((n2, n2))

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

    for ibra in range(n2):
        if ibra < nK:
            bra = tbc.GetKet(ibra)
            a, b = bra.p, bra.q
        else:
            bra = tbc.GetKet(ibra - nK)
            a, b = bra.q, bra.p
        oa, ob = ms.GetOrbit(a), ms.GetOrbit(b)
        ja, jb = oa.j2 * 0.5, ob.j2 * 0.5
        for iket in range(n2):
            if iket < nK:
                ket = tbc.GetKet(iket)
                c, d = ket.p, ket.q
            else:
                ket = tbc.GetKet(iket - nK)
                c, d = ket.q, ket.p
            oc, od = ms.GetOrbit(c), ms.GetOrbit(d)
            jc, jd = oc.j2 * 0.5, od.j2 * 0.5
            jmin = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
            jmax = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
            X = 0.0
            for Jp in range(jmin, jmax + 1):
                six = SixJ(ja, jb, Jcc, jc, jd, Jp)
                if abs(six) < 1e-8:
                    continue
                X -= (2 * Jp + 1) * six * chi_unnorm(Jp, a, d, c, b)
            bC[ibra, iket] = X

    bar_Gamma.append(bG)
    barCHI.append(bC)

for ch in range(n_cc):
    tbc = ms.GetTwoBodyChannel_CC(ch)
    nK = tbc.GetNumberKets()
    Jcc = tbc.J
    if nK < 1 or barCHI[ch] is None or bar_Gamma[ch] is None:
        CHI_final.append(None)
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
                mat = barCHI[ch_old]
                me = mat[indx_bc, indx_ad] + mat[indx_ad, indx_bc]
                X -= (
                    iphase((ob.j2 + oc.j2) // 2 + Jp)
                    * (2 * Jp + 1)
                    * six
                    * me
                )
            RC[ibra, iket] = X
    CHI_final.append(bar_Gamma[ch] @ RC)

print(f"  Python Path B pack {time.time() - t0:.1f}s")


def chi_final_read(a, b, c, d, J, parity, Tz) -> float:
    ch = ms.GetTwoBodyChannelIndex(J, parity, Tz)
    if ch < 0 or CHI_final[ch] is None:
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
    M = CHI_final[ch]
    if iab >= M.shape[0] or icd >= M.shape[1]:
        return 0.0
    return float(M[iab, icd])


def z_pb(J0, i, j, k, l) -> float:
    oi, oj = ms.GetOrbit(i), ms.GetOrbit(j)
    ok, ol = ms.GetOrbit(k), ms.GetOrbit(l)
    ji, jj, jk, jl = oi.j2, oj.j2, ok.j2, ol.j2
    if not (tri(ji * 0.5, jj * 0.5, J0) and tri(jk * 0.5, jl * 0.5, J0)):
        return 0.0
    cijkl = cjikl = cijlk = cjilk = 0.0
    parity = (oi.l + ok.l) % 2
    Tz = abs(oi.tz2 - ok.tz2) // 2
    for Jp in range(
        max(abs(jj - jl), abs(ji - jk)) // 2, min(jj + jl, ji + jk) // 2 + 1
    ):
        six = SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jp)
        if abs(six) < 1e-8:
            continue
        me1 = chi_final_read(j, l, i, k, Jp, parity, Tz)
        me2 = chi_final_read(i, k, j, l, Jp, parity, Tz)
        cjikl -= iphase(Jp + (ji + jk) // 2) * (2 * Jp + 1) * six * me1
        cijlk -= iphase(Jp + (jj + jl) // 2) * (2 * Jp + 1) * six * me2
    parity = (oi.l + ol.l) % 2
    Tz = abs(oi.tz2 - ol.tz2) // 2
    for Jp in range(
        max(abs(ji - jl), abs(jj - jk)) // 2, min(ji + jl, jj + jk) // 2 + 1
    ):
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
    return iphase(J0 + (ji + jj) // 2) * z


print("ethS GIIIb ...")
t1 = time.time()
Z_fac = Operator(ms, 0, 0, 0, 2)
Z_fac.SetHermitian()
Commutator.FactorizedDoubleCommutator_eths.comm223_232_GIIIb(Omega, Gamma, Z_fac)
print(f"  ethS wall {time.time() - t1:.2f}s  ‖Z‖={Z_fac.TwoBodyNorm():.6g}")

rats = Counter()
n = 0
max_ad = 0.0
fails = []
for ch in range(ms.GetNumberTwoBodyChannels()):
    tbc = ms.GetTwoBodyChannel(ch)
    J0 = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        i, j = tbc.GetKet(ib).p, tbc.GetKet(ib).q
        for ik in range(nk):
            k, l = tbc.GetKet(ik).p, tbc.GetKet(ik).q
            zc = Z_fac.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
            zg = z_pb(J0, i, j, k, l)
            if abs(zc) < 1e-12 and abs(zg) < 1e-12:
                continue
            n += 1
            err = abs(zc - zg)
            max_ad = max(max_ad, err)
            r = zc / zg if abs(zg) > 1e-8 else float("nan")
            if abs(zg) > 1e-8:
                rats[round(r, 6)] += 1
            if err > tol:
                fails.append((err, J0, i, j, k, l, zc, zg, r, i == j, k == l))

print(
    f"ethS ≡ Python Path B (JT): n={n}  max|Δ|={max_ad:.3e}  "
    f"eth/gold={rats.most_common(8)}"
)
fails.sort(reverse=True)
print(f"  nfail={len(fails)}")
for f in fails[:12]:
    err, J0, i, j, k, l, zc, zg, r, iij, kkl = f
    print(
        f"  J={J0} ({i},{j})({k},{l}) i=j={iij} k=l={kkl} "
        f"eth={zc:.5e} py={zg:.5e} r={r:.5g} Δ={err:.3e}"
    )
ok = n > 5 and max_ad < tol
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
