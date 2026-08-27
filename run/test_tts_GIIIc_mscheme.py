#!/usr/bin/env python3
"""Γ^{III_c}: m ≡ AMC fold ≡ DIRECT ≡ Path B (after χ packaging fix).

Packaging (locked)
------------------
χ^θ T×T→S: S = Σ (−1)^{J0+J2+λ} λ̂^{-1} w ΩΩ
  χ_red = S/Ĵ, χ_unred = S/Ĵ²
  m: χ(m) = Σ w CG(λμ;λ−μ;00) ΩΩ ; χ_red = Σ CG CG χ(m)/Ĵ

Fold: Z = −½(1−Pij)(1−Pkl) Σ χ_iabl Γ_bjka
  Unreduced AMC: Ĵ2² Ĵ3² χ_unred Γ_unred  (G3c_from_chi_ninej.tex)
  Case-2: Ĵ0 Ĵ2 Ĵ3 χ_red Γ_red → Z_unred = Z_red/Ĵ0
  Code DIRECT uses Case-2; Path B uses χ_unred via Pandya.

m gold: −½(R1−R2−R3+R4) with χ from unpack(χ_unred), Γ(m)
Compare via 2-CG unpack of GetTBME_J on ordered bras (i<j, k<l).

Usage:
  PYTHONPATH=build python3 -B run/test_tts_GIIIc_mscheme.py [emax=1] [lambda=2]
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


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def m_range(o):
    return range(-j2i(o), j2i(o) + 1, 2)


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1)
max_J = max(j2i(o) for o in orbits)

chiS: dict = {}


def chi_S(i, j, k, l, J0) -> float:
    key = (i, j, k, l, J0)
    if key in chiS:
        return chiS[key]
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        chiS[key] = 0.0
        return 0.0
    sm = 0.0
    for a in orbits:
        for b in orbits:
            w = (
                occ(a) * occ(b) * nbar(k)
                + nbar(a) * nbar(b) * occ(k)
                + occ(a) * occ(b) * nbar(j)
                + nbar(a) * nbar(b) * occ(j)
            )
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
    chiS[key] = sm
    return sm


def chi_to_m(i, mi, a, ma, b, mb, l, ml) -> float:
    if mi + ma != mb + ml:
        return 0.0
    M = (mi + ma) // 2
    sm = 0.0
    Jmin = max(abs(j2i(i) - j2i(a)), abs(j2i(b) - j2i(l))) // 2
    Jmax = min(j2i(i) + j2i(a), j2i(b) + j2i(l)) // 2
    for J2 in range(Jmin, Jmax + 1):
        if abs(M) > J2:
            continue
        c1 = CG(0.5 * j2i(i), 0.5 * mi, 0.5 * j2i(a), 0.5 * ma, J2, M)
        c2 = CG(0.5 * j2i(b), 0.5 * mb, 0.5 * j2i(l), 0.5 * ml, J2, M)
        if abs(c1 * c2) < 1e-16:
            continue
        sm += c1 * c2 * chi_S(i, a, b, l, J2) / (hat(J2) ** 2)
    return sm


def gIIIc_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in m_range(a):
            for b in orbits:
                for mb in m_range(b):
                    chi = chi_to_m(i, mi, a, ma, b, mb, l, ml)
                    if abs(chi) > 1e-16:
                        sm += chi * ut.GetMschemeMatrixElement_2b(
                            Gamma, b, mb, j, mj, k, mk, a, ma
                        )
                    chi = chi_to_m(i, mi, a, ma, b, mb, k, mk)
                    if abs(chi) > 1e-16:
                        sm -= chi * ut.GetMschemeMatrixElement_2b(
                            Gamma, b, mb, j, mj, l, ml, a, ma
                        )
                    chi = chi_to_m(j, mj, a, ma, b, mb, l, ml)
                    if abs(chi) > 1e-16:
                        sm -= chi * ut.GetMschemeMatrixElement_2b(
                            Gamma, b, mb, i, mi, k, mk, a, ma
                        )
                    chi = chi_to_m(j, mj, a, ma, b, mb, k, mk)
                    if abs(chi) > 1e-16:
                        sm += chi * ut.GetMschemeMatrixElement_2b(
                            Gamma, b, mb, i, mi, l, ml, a, ma
                        )
    return -0.5 * sm


def unpack_Z(Z, i, mi, j, mj, k, mk, l, ml) -> float:
    if mi + mj != mk + ml:
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
        sm += cab * ccd * Z.TwoBody.GetTBME_J(J, J, i, j, k, l)
    return sm


def fresh_Z():
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    return Z


print("DIRECT ...")
t0 = time.time()
Z_D = fresh_Z()
ReferenceImplementations.comm223_232_tts_GIIIc(Eta, Gamma, Z_D, 0)
print(f"  ‖Z_D‖={Z_D.TwoBodyNorm():.6g}  ({time.time()-t0:.2f}s)")

print("Path B ...")
t0 = time.time()
cm = Commutator.FactorizedDoubleCommutator_eths
# comm223_232_GIIIc is Path B only -- the old factorized/slow/which_term
# toggles no longer exist, so call it directly.
Z_B = fresh_Z()
cm.comm223_232_GIIIc(Eta, Gamma, Z_B)
diff_DB = (Z_D - Z_B).TwoBodyNorm()
print(f"  ‖Z_B‖={Z_B.TwoBodyNorm():.6g}  ‖D−B‖={diff_DB:.3e}")

print("\nCompare m vs unpack(DIRECT/PathB) [ordered i<j, k<l] ...")
t0 = time.time()
rats_D = Counter()
rats_B = Counter()
n = 0
max_D = 0.0
max_B = 0.0
worst = None

for i in orbits:
    for j in orbits:
        if i >= j:
            continue
        for k in orbits:
            for l in orbits:
                if k >= l:
                    continue
                if len({i, j, k, l}) < 4:
                    continue
                for mi in m_range(i):
                    for mj in m_range(j):
                        for mk in m_range(k):
                            ml = mi + mj - mk
                            if ml not in m_range(l):
                                continue
                            m = gIIIc_m(i, mi, j, mj, k, mk, l, ml)
                            d = unpack_Z(Z_D, i, mi, j, mj, k, mk, l, ml)
                            b = unpack_Z(Z_B, i, mi, j, mj, k, mk, l, ml)
                            if abs(m) < 1e-8 and abs(d) < 1e-8:
                                continue
                            n += 1
                            dd = abs(m - d)
                            db = abs(m - b)
                            if dd > max_D:
                                max_D = dd
                                worst = (i, j, k, l, (mi, mj, mk, ml), m, d, b)
                            if db > max_B:
                                max_B = db
                            if abs(d) > 1e-8:
                                rats_D[round(m / d, 4)] += 1
                            if abs(b) > 1e-8:
                                rats_B[round(m / b, 4)] += 1

print(
    f"  n={n}  max|m−D|={max_D:.3e}  m/D={rats_D.most_common(4)}\n"
    f"         max|m−B|={max_B:.3e}  m/B={rats_B.most_common(4)}  "
    f"({time.time()-t0:.2f}s)"
)
if worst is not None:
    i, j, k, l, ms_, m, d, b = worst
    print(
        f"  worst ({i},{j},{k},{l}) m={ms_}: m={m:.6e} D={d:.6e} B={b:.6e} "
        f"rD={m/d if abs(d)>1e-14 else float('nan'):g}"
    )

ok = (
    diff_DB < tol
    and max_D < tol
    and max_B < tol
    and (not rats_D or rats_D.most_common(1)[0][0] == 1.0)
    and (not rats_B or rats_B.most_common(1)[0][0] == 1.0)
)
print(
    "\nPASS — Γ^{III_c} m ≡ DIRECT ≡ Path B"
    if ok
    else "\nFAIL"
)
sys.exit(0 if ok else 1)
