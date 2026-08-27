#!/usr/bin/env python3
"""ethS Γ^{IV_a} ≡ Path B gold (χ→W→(1−P) on W), packaging locked.

Packaging:
  Path B Python: Z_red
  ethS store:    Z_unred = Z_red / Ĵ  (+ √2 for i=j / k=l)
  Compare: ethS GetTBME ≡ Z_red / Ĵ / √2…

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIVa_eths_vs_pathB.py [emax=1] [lambda=2]
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

if lam == 0:
    print(
        "SKIP — λ=0 uses Factorized scalar CHI_VI (not rectangular Path B).\n"
        "       Compare λ≠0, or use test_chi_kappa_pandya_dgemm.py for λ=0 twin."
    )
    sys.exit(0)

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

hEta = -1 if Eta.IsAntiHermitian() else 1
orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1) if lam else 1.0


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


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def w_ABbarD(a, b, d) -> float:
    return occ(a) * nbar(b) * occ(d) + nbar(a) * occ(b) * nbar(d)


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


print(f"emax={emax} λ={lam} seed={seed}")
print("Caching Path B χ ...")
t0 = time.time()
chiJ: dict = {}
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
print(f"  nonzero={len(chiJ)} ({time.time()-t0:.1f}s)")


def w_rme(J0, i, j, k, l) -> float:
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    pref = -iphase(J0) / hat(J0)
    tot = 0.0
    for b in orbits:
        for d in orbits:
            ph_bd = iphase((j2i(b) + j2i(d)) // 2 + lam)
            for J2 in range(0, max_J + 1):
                if not tri(J0, J2, lam) or not tri(jo(b), jo(d), J2):
                    continue
                ch = chiJ.get((J0, J2, i, j, b, d), 0.0)
                if abs(ch) < 1e-16:
                    continue
                om = Eta.TwoBody.GetTBME_J(J2, J0, d, b, k, l)
                if abs(om) < 1e-16:
                    continue
                tot += ph_bd * hat_lam_inv * ch * om
    return pref * tot


def z_red(J0, i, j, k, l) -> float:
    w1 = w_rme(J0, i, j, k, l)
    w2 = w_rme(J0, j, i, k, l)
    bra = w1 - iphase((j2i(i) + j2i(j)) // 2 - J0) * w2
    wk1 = w_rme(J0, k, l, i, j)
    wk2 = w_rme(J0, l, k, i, j)
    ket = wk1 - iphase((j2i(k) + j2i(l)) // 2 - J0) * wk2
    return bra + ket


def z_unred_GetTBME(J0, i, j, k, l) -> float:
    """Match ethS GetTBME_J: Z_red/Ĵ (GetTBME undoes √2 store norms)."""
    return z_red(J0, i, j, k, l) / hat(J0)


print("ethS GIVa ...")
cm = Commutator.FactorizedDoubleCommutator_eths
cm.SetUse_1b_Intermediates(False)
cm.SetUse_2b_Intermediates(False)
cm.SetUse_TypeGIVa_2b(True)
for name in (
    "SetUse_TypeGI_2b",
    "SetUse_TypeGII_2b",
    "SetUse_TypeGIIIa_2b",
    "SetUse_TypeGIIIb_2b",
    "SetUse_TypeGIIIc_2b",
    "SetUse_TypeGIVb_2b",
    "SetUse_TypeGIVc_2b",
):
    fn = getattr(cm, name, None)
    if fn:
        fn(False)

Z = Operator(ms, 0, 0, 0, 2)
Z.SetHermitian()
t0 = time.time()
cm.comm223_232_GIVa(Eta, Gamma, Z)
print(f"  ‖Z‖={Z.TwoBodyNorm():.6e} ({time.time()-t0:.2f}s)")

print("Compare ethS ≡ Path B store ...")
rats = Counter()
n = 0
max_abs = 0.0
worst = None
for ch in range(ms.GetNumberTwoBodyChannels()):
    tbc = ms.GetTwoBodyChannel(ch)
    J = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        i, j = tbc.GetKet(ib).p, tbc.GetKet(ib).q
        for ik in range(ib, nk):
            k, l = tbc.GetKet(ik).p, tbc.GetKet(ik).q
            zp = z_unred_GetTBME(J, i, j, k, l)
            ze = Z.TwoBody.GetTBME_J(J, J, i, j, k, l)
            if abs(zp) < 1e-10 and abs(ze) < 1e-10:
                continue
            n += 1
            err = abs(zp - ze)
            if err > max_abs:
                max_abs = err
                r = zp / ze if abs(ze) > 1e-10 else float("nan")
                worst = (J, i, j, k, l, zp, ze, r)
            if abs(ze) > 1e-10:
                rats[round(zp / ze, 4)] += 1

print(f"  n={n} max|Δ|={max_abs:.3e} PathB/ethS={rats.most_common(5)}")
if worst:
    J, i, j, k, l, zp, ze, r = worst
    print(f"  worst J={J} ({i},{j}|{k},{l}): PB={zp:.6e} ethS={ze:.6e} r={r:.6g}")

ok = n > 0 and max_abs < tol
print(
    "\nPASS — ethS ≡ Path B Γ^{IV_a}"
    if ok
    else f"\nFAIL maxΔ={max_abs:.3e}"
)
sys.exit(0 if ok else 1)
