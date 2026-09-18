#!/usr/bin/env python3
"""λ=3: lock AMC Path B Term1 / Path A vs CG-gold Z (kind=0 ethS).

Gold: ethS givc_leftover_kind=0 (CG) → Z_unred GetTBME_J.
AMC Path A: G4c_from_chi_reduced Term1 (printed ninej) → Z_red = A; compare A/Ĵ.
AMC Path B: same-label Pandya → Eq3 Term1 → inv Eq4 (no sample minus).

Packaging:
  χ, Ω WE-reduced; AMC reduce=true → Z_red; gold store Z_unred = Z_red/Ĵ.

Usage:
  PYTHONPATH=build python3 -B run/test_givc_pathB_term1_vs_cg.py [emax=1] [lambda=3]
"""
from __future__ import annotations

import math
import os
import sys
from collections import Counter

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 3
seed = 17
tol = 1e-6

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
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
max_j = max(ms.GetOrbit(o).j2 * 0.5 for o in orbits)
hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1)


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


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()}  Γ reduced={Gamma.IsReduced()}"
)


def chi_amc(J0, J1, i, j, k, l) -> float:
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


print("χ cache ...")
chiJ: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        v = chi_amc(J0, J1, i, j, k, l)
                        if abs(v) > 1e-16:
                            chiJ[(J0, J1, i, j, k, l)] = v
print(f"  nonzero={len(chiJ)}")


# ---- CG gold (ethS kind=0), T1 only ----
cm = Commutator.FactorizedDoubleCommutator_eths
cm.SetGIVcChiWhich(1)
cm.SetGIVcFoldAS(True)
cm.SetGIVcLeftoverKind(0)
Zg = Operator(ms, 0, 0, 0, 2)
Zg.SetHermitian()
Zg *= 0.0
cm.comm223_232_GIVc(Eta, Gamma, Zg)
print(f"  CG gold T1 ||Z||={Zg.TwoBodyNorm():.6g} reduced={Zg.IsReduced()}")


# ---- Path A Term1 (G4c_from_chi_reduced print) ----
def pathA_term1(i, j, k, l, J0) -> float:
    """Printed Term1 → Z_red."""
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    p = iphase((j2i(i) + j2i(l)) // 2)
    sm = 0.0
    for a in orbits:
        ja = jo(a)
        for b in orbits:
            jb = jo(b)
            for J2 in range(0, max_J + 1):
                for J3 in range(0, max_J + 1):
                    if not tri(J2, J3, lam):
                        continue
                    chi = chiJ.get((J2, J3, i, a, l, b), 0.0)
                    if abs(chi) < 1e-16:
                        continue
                    for J4 in range(0, max_J + 1):
                        for J5 in range(0, max_J + 1):
                            if not tri(J4, J5, lam):
                                continue
                            om = Eta.TwoBody.GetTBME_J(J4, J5, b, j, a, k)
                            if abs(om) < 1e-16:
                                continue
                            for two_j0 in range(0, int(2 * max_j) + 3):
                                j0 = 0.5 * two_j0
                                # Printed: 6j(J3 λ J2; ja ji j0), 6j(J4 λ J5; ja jk j0)
                                s1 = SixJ(J3, lam, J2, ja, jo(i), j0)
                                s2 = SixJ(J4, lam, J5, ja, jo(k), j0)
                                n9 = NineJ(
                                    jo(l), jb, J3, jo(k), J4, j0, J0, jo(j), jo(i)
                                )
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                ph = iphase(
                                    J2
                                    + J3
                                    + J4
                                    + J5
                                    + (j2i(a) + j2i(b)) // 2
                                    + lam
                                )
                                sm += (
                                    ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (hat(j0) ** 2)
                                    * hat_lam_inv
                                    * s1
                                    * s2
                                    * n9
                                    * chi
                                    * om
                                )
    return 0.5 * p * hat(J0) * sm


# ---- Path B same-label Pandya + Term1 product + inv ----
def amc_fwd(get_me, Jbra, Jket, i, j, k, l) -> float:
    if not tri(Jbra, Jket, lam):
        return 0.0
    if not (tri(jo(i), jo(l), Jbra) and tri(jo(k), jo(j), Jket)):
        return 0.0
    sm = 0.0
    for J2 in range(0, max_J + 1):
        for J3 in range(0, max_J + 1):
            if not tri(J2, J3, lam):
                continue
            if not (tri(jo(i), jo(j), J2) and tri(jo(k), jo(l), J3)):
                continue
            me = get_me(J2, J3, i, j, k, l)
            if abs(me) < 1e-16:
                continue
            n9 = NineJ(lam, Jbra, Jket, J3, jo(l), jo(k), J2, jo(i), jo(j))
            if abs(n9) < 1e-16:
                continue
            sm += iphase(J2) * hat(J2) * hat(J3) * n9 * me
    pref = -iphase(Jbra + (j2i(i) + j2i(k)) // 2 + lam) * hat(Jbra) * hat(Jket)
    return pref * sm


def get_chi(J2, J3, a, b, c, d):
    return chiJ.get((J2, J3, a, b, c, d), 0.0)


def get_om(J2, J3, a, b, c, d):
    return Eta.TwoBody.GetTBME_J(J2, J3, a, b, c, d)


print("Path B fwd Pandya cache ...")
barChi: dict = {}
barO: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        bc = amc_fwd(get_chi, J0, J1, i, j, k, l)
                        bo = amc_fwd(get_om, J0, J1, i, j, k, l)
                        if abs(bc) > 1e-16:
                            barChi[(J0, J1, i, j, k, l)] = bc
                        if abs(bo) > 1e-16:
                            barO[(J0, J1, i, j, k, l)] = bo
print(f"  barχ={len(barChi)} barΩ={len(barO)}")


def pathB_bar_term1(i, j, k, l, J0) -> float:
    """AMC Eq3 Term1 → barG reduced (before inv)."""
    sm = 0.0
    p = iphase((j2i(i) + j2i(l)) // 2)
    for a in orbits:
        for b in orbits:
            for J2 in range(0, max_J + 1):
                for J3 in range(0, max_J + 1):
                    if not tri(J2, J3, lam):
                        continue
                    bc = barChi.get((J2, J3, i, a, l, b), 0.0)
                    if abs(bc) < 1e-16:
                        continue
                    for J4 in range(0, max_J + 1):
                        for J5 in range(0, max_J + 1):
                            if not tri(J4, J5, lam):
                                continue
                            bo = barO.get((J4, J5, b, j, a, k), 0.0)
                            if abs(bo) < 1e-16:
                                continue
                            for two_j0 in range(0, int(2 * max_j) + 3):
                                j0 = 0.5 * two_j0
                                s1 = SixJ(lam, J3, J2, jo(i), jo(b), j0)
                                s2 = SixJ(jo(i), jo(l), J0, jo(a), j0, J3)
                                s3 = SixJ(jo(a), jo(j), J5, jo(k), j0, J0)
                                s4 = SixJ(J5, lam, J4, jo(b), jo(k), j0)
                                if abs(s1 * s2 * s3 * s4) < 1e-16:
                                    continue
                                ph = iphase(J2 + J3 + (j2i(a) + j2i(b)) // 2)
                                sm += (
                                    ph
                                    * hat(J2)
                                    * hat(J3)
                                    * hat(J4)
                                    * hat(J5)
                                    * (hat(j0) ** 2)
                                    * hat_lam_inv
                                    * s1
                                    * s2
                                    * s3
                                    * s4
                                    * bc
                                    * bo
                                )
    return 0.5 * p * hat(J0) * sm


def scalar_inv(get_bar, i, j, k, l, J0) -> float:
    """AMC Eq4 corrected (no overall minus) → Z_red."""
    sm = 0.0
    for J2 in range(0, max_J + 1):
        bz = get_bar(i, j, k, l, J2)
        if abs(bz) < 1e-16:
            continue
        s = SixJ(jo(l), jo(k), J0, jo(j), jo(i), J2)
        if abs(s) < 1e-16:
            continue
        sm += hat(J2) * s * bz
    return hat(J0) * sm


barZ1: dict = {}


def ensure_barZ1(i, j, k, l, J0):
    key = (i, j, k, l, J0)
    if key not in barZ1:
        barZ1[key] = pathB_bar_term1(i, j, k, l, J0)
    return barZ1[key]


def cmp(tag, pred_red):
    """pred_red → Z_red; gold is Z_unred GetTBME."""
    rats = Counter()
    n = 0
    max_abs = 0.0
    worst = None
    nch = ms.GetNumberTwoBodyChannels()
    for ch in range(nch):
        tbc = ms.GetTwoBodyChannel(ch)
        J0 = tbc.J
        hatJ = hat(J0)
        nk = tbc.GetNumberKets()
        for ib in range(nk):
            bra = tbc.GetKet(ib)
            i, j = int(bra.p), int(bra.q)
            for ik in range(ib, nk):
                ket = tbc.GetKet(ik)
                k, l = int(ket.p), int(ket.q)
                g = Zg.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
                # identical-orbit √2 already in GetTBME store
                pr = pred_red(i, j, k, l, J0)
                # convert reduced → unreduced; match identical √2
                pu = pr / hatJ
                if i == j:
                    pu /= math.sqrt(2.0)
                if k == l:
                    pu /= math.sqrt(2.0)
                if abs(g) < 1e-10 and abs(pu) < 1e-10:
                    continue
                n += 1
                d = abs(g - pu)
                if d > max_abs:
                    max_abs = d
                    worst = (i, j, k, l, J0, g, pu, pr)
                if abs(g) > 1e-10:
                    rats[round(pu / g, 4)] += 1
    top = rats.most_common(5)
    ok = n > 0 and max_abs < tol and (top and top[0][0] == 1.0)
    print(
        f"  {tag}: n={n} max|Δ|={max_abs:.3e} pu/gold={top} => "
        f"{'PASS' if ok else 'FAIL'}"
    )
    if worst and not ok:
        i, j, k, l, J0, g, pu, pr = worst
        print(
            f"    worst ({i},{j},{k},{l})J={J0}: gold={g:.6e} "
            f"pred_unred={pu:.6e} pred_red={pr:.6e} r={pu/g if abs(g)>1e-14 else float('nan'):g}"
        )
    return ok


print("\n--- vs CG gold T1 (Z_unred) ---")
ok_a = cmp("PathA Term1 (print)", pathA_term1)
ok_b = cmp(
    "PathB Term1→inv",
    lambda i, j, k, l, J0: scalar_inv(ensure_barZ1, i, j, k, l, J0),
)

# Path A vs Path B (both reduced, no √2)
print("\n--- PathA Term1 ≡ PathB Term1→inv (Z_red, no √2) ---")
rats = Counter()
n = 0
max_abs = 0.0
worst = None
for i in orbits:
    for j in orbits:
        if i > j:
            continue
        for k in orbits:
            for l in orbits:
                if k > l:
                    continue
                for J0 in range(0, max_J + 1):
                    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
                        continue
                    a = pathA_term1(i, j, k, l, J0)
                    b = scalar_inv(ensure_barZ1, i, j, k, l, J0)
                    if abs(a) < 1e-8 and abs(b) < 1e-8:
                        continue
                    n += 1
                    d = abs(a - b)
                    if d > max_abs:
                        max_abs = d
                        worst = (i, j, k, l, J0, a, b)
                    if abs(a) > 1e-8:
                        rats[round(b / a, 4)] += 1
print(f"  n={n} max|Δ|={max_abs:.3e} B/A={rats.most_common(4)}")
if worst:
    i, j, k, l, J0, a, b = worst
    print(
        f"  worst ({i},{j},{k},{l})J={J0}: A={a:.6e} B={b:.6e} "
        f"r={b/a if abs(a)>1e-14 else float('nan'):g}"
    )

cm.SetGIVcChiWhich(0)
print(
    "\nPASS — Path B Term1 ≡ CG"
    if ok_b
    else "\nFAIL — Path B Term1 not locked to CG gold at this λ"
)
print(
    "PASS — Path A Term1 ≡ CG"
    if ok_a
    else "FAIL — Path A Term1 (known often wrong vs m; diagnostic only)"
)
sys.exit(0 if ok_b else 1)
