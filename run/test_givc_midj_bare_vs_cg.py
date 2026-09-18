#!/usr/bin/env python3
"""λ=3: lock bare mid-J leftover vs ethS Path B (Pandya→DGEMM→inv).

Gold: ethS production leftover, FoldAS=False, T1 only → Z_unred = K recoupled
      K ≡ tts_ring X_pqsr = Σ χ_pbar Ω_aqsb

Mid-J candidates:
  Path A: ring_X (tts_ring)
  Path B: adcb Pandya mid+inv (same ring; drop AMC-sample minus)
  Diagnostic: printed AMC G4c_from_chi_noperm (known wrong vs m)

See learn/amc_tts/tensor_pro_final/GIVc/03_pathB/LESSON.md

Usage:
  PYTHONPATH=build python3 -B run/test_givc_midj_bare_vs_cg.py [emax=1] [lambda=3]
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
tol = 1e-7

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
hat_lam_inv = 1.0 / math.sqrt(2 * lam + 1)
max_j = max(ms.GetOrbit(o).j2 * 0.5 for o in orbits)


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


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def iphase(n: int) -> float:
    return 1.0 if int(n) % 2 == 0 else -1.0


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


print(
    f"emax={emax} λ={lam} seed={seed}\n"
    f"  Ω reduced={Eta.IsReduced()}  Γ reduced={Gamma.IsReduced()}\n"
    f"  lock: bare K (FoldAS=False), T1 only"
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
            if abs(w_l) > 1e-12 and tri(jo(a), jo(b), J0):
                sm += w_l * (
                    Gamma.TwoBody.GetTBME_J(J0, J0, i, j, a, b)
                    * Eta.TwoBody.GetTBME_J(J0, J1, a, b, k, l)
                )
    return sm


print("χ T1 cache ...")
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

# ---- CG gold: FoldAS=False, kind=0, T1 ----
cm = Commutator.FactorizedDoubleCommutator_eths
cm.SetGIVcChiWhich(1)
cm.SetGIVcFoldAS(False)
# Production leftover = Path B only; gold here is ethS bare K.
Zg = Operator(ms, 0, 0, 0, 2)
Zg.SetHermitian()
Zg *= 0.0
cm.comm223_232_GIVc(Eta, Gamma, Zg)
print(f"  ethS bare PathB ||Z||={Zg.TwoBodyNorm():.6g}")


def ring_X_red(p, q, s, r, J0) -> float:
    """tts_ring Path A: X_pqsr = Σ χ_pbar Ω_aqsb → reduced."""
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


def amc_noperm_unred(i, j, k, l, J0) -> float:
    """Printed G4c_from_chi_noperm (unreduced LHS) — diagnostic only."""
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
    return p * sm


def cmp(tag, pred_unred_from_ij):
    """pred returns unreduced GetTBME-like (incl √2 for identical)."""
    rats = Counter()
    n = 0
    max_abs = 0.0
    worst = None
    samples = []
    nch = ms.GetNumberTwoBodyChannels()
    for ch in range(nch):
        tbc = ms.GetTwoBodyChannel(ch)
        J0 = tbc.J
        nk = tbc.GetNumberKets()
        for ib in range(nk):
            bra = tbc.GetKet(ib)
            i, j = int(bra.p), int(bra.q)
            for ik in range(ib, nk):
                ket = tbc.GetKet(ik)
                k, l = int(ket.p), int(ket.q)
                g = Zg.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
                pu = pred_unred_from_ij(i, j, k, l, J0)
                if abs(g) < 1e-10 and abs(pu) < 1e-10:
                    continue
                n += 1
                d = abs(g - pu)
                if d > max_abs:
                    max_abs = d
                    worst = (i, j, k, l, J0, g, pu)
                if abs(g) > 1e-10:
                    rats[round(pu / g, 4)] += 1
                if len(samples) < 4 and abs(g) > 1e-8:
                    samples.append((i, j, k, l, J0, g, pu))
    top = rats.most_common(6)
    ok = n > 0 and max_abs < tol and (top and abs(top[0][0] - 1.0) < 1e-3)
    print(
        f"  {tag}: n={n} max|Δ|={max_abs:.3e} pu/gold={top} => "
        f"{'PASS' if ok else 'FAIL'}"
    )
    if worst and not ok:
        i, j, k, l, J0, g, pu = worst
        print(
            f"    worst ({i},{j}|{k},{l})J={J0}: gold={g:.6e} "
            f"pred={pu:.6e} r={pu/g if abs(g)>1e-14 else float('nan'):g}"
        )
    for i, j, k, l, J0, g, pu in samples:
        print(
            f"    sample ({i},{j}|{k},{l})J={J0}: gold={g:.6e} "
            f"pred={pu:.6e} r={pu/g if abs(g)>1e-14 else float('nan'):g}"
        )
    return ok


def unred_from_red(J0, z_red):
    """GetTBME_J compares to unreduced tilde (√2 already restored on identical)."""
    return z_red / hat(J0)


# ---- Path B ring: IMSRG adcb Pandya + AMC↔adcb map (locked) ----
# AMC bar(i,j,k,l) ≡ IMSRG adcb(i,l,k,j).
# Mid: AMC barχ(p,b,a,r)·barΩ(a,q,s,b)
#    ≡ adcbχ(p,r,a,b)·adcbΩ(a,b,s,q)
def pandya_adcb(get_me, a, b, c, d, Jbra, Jket) -> float:
    """IMSRG DoTensorPandya adcb: bar(a,b,c,d) ← O(a,d,c,b)."""
    if not tri(Jbra, Jket, lam):
        return 0.0
    if not (tri(jo(a), jo(b), Jbra) and tri(jo(c), jo(d), Jket)):
        return 0.0
    sm = 0.0
    for J1 in range(abs(j2i(a) - j2i(d)) // 2, (j2i(a) + j2i(d)) // 2 + 1):
        for J2 in range(
            max(abs(j2i(c) - j2i(b)) // 2, abs(J1 - lam)),
            min((j2i(c) + j2i(b)) // 2, J1 + lam) + 1,
        ):
            if not tri(J1, J2, lam):
                continue
            n9 = NineJ(
                jo(a), jo(d), J1, jo(b), jo(c), J2, Jbra, Jket, lam
            )
            if abs(n9) < 1e-16:
                continue
            tbme = get_me(J1, J2, a, d, c, b)
            sm -= (
                hat(J1)
                * hat(J2)
                * hat(Jbra)
                * hat(Jket)
                * phase(jo(b) + jo(d) + Jket + J2)
                * n9
                * tbme
            )
    return sm


def get_chi(J2, J3, w, x, y, z):
    return chiJ.get((J2, J3, w, x, y, z), 0.0)


def get_om(J2, J3, w, x, y, z):
    return Eta.TwoBody.GetTBME_J(J2, J3, w, x, y, z)


bar_chi: dict = {}
bar_om: dict = {}


def get_bar_chi_adcb(p, r, a, b, Jpr, Jab):
    """AMC barχ(p,b,a,r) via adcb(p,r,a,b)."""
    key = (p, r, a, b, Jpr, Jab)
    if key not in bar_chi:
        bar_chi[key] = pandya_adcb(get_chi, p, r, a, b, Jpr, Jab)
    return bar_chi[key]


def get_bar_om_adcb(a, b, s, q, Jab, Jsq):
    """AMC barΩ(a,q,s,b) via adcb(a,b,s,q)."""
    key = (a, b, s, q, Jab, Jsq)
    if key not in bar_om:
        bar_om[key] = pandya_adcb(get_om, a, b, s, q, Jab, Jsq)
    return bar_om[key]


def ring_barz_mid(p, q, s, r, Jp) -> float:
    """barX mid at CC J=Jp on (p,r)×(s,q)."""
    sm = 0.0
    for J2 in range(0, max_J + lam + 2):
        if not tri(Jp, J2, lam):
            continue
        for a in orbits:
            for b in orbits:
                bc = get_bar_chi_adcb(p, r, a, b, Jp, J2)
                if abs(bc) < 1e-16:
                    continue
                bo = get_bar_om_adcb(a, b, s, q, J2, Jp)
                if abs(bo) < 1e-16:
                    continue
                sm += phase(J2 + lam) * hat_lam_inv * bc * bo
    return phase(Jp) / hat(Jp) * sm


def ring_pathB_red(p, q, s, r, J0) -> float:
    """Corrected inv (drop AMC-sample overall minus) → Z_red ≡ Path A."""
    jp, jq, js, jr = jo(p), jo(q), jo(s), jo(r)
    if not (tri(jp, jq, J0) and tri(js, jr, J0)):
        return 0.0
    sm = 0.0
    Jpmax = max((j2i(p) + j2i(r)) // 2, (j2i(q) + j2i(s)) // 2) + 1
    for Jp in range(0, Jpmax + 1):
        six = SixJ(jr, js, J0, jq, jp, Jp)
        if abs(six) < 1e-16:
            continue
        sm += hat(Jp) * six * ring_barz_mid(p, q, s, r, Jp)
    return hat(J0) * sm


print("\n--- bare mid-J vs CG gold (FoldAS=False) ---")
# GetTBME_J = unreduced tilde (identical √2 restored). Do NOT ÷√2 again.
ok_ring = cmp(
    "tts_ring PathA X_pqsr",
    lambda i, j, k, l, J0: unred_from_red(J0, ring_X_red(i, j, k, l, J0)),
)
ok_amc = cmp(
    "AMC noperm print (χ_ialb)",
    lambda i, j, k, l, J0: amc_noperm_unred(i, j, k, l, J0),
)
print("Path B ring (adcb map, minus dropped) ...")
ok_b = cmp(
    "tts_ring PathB adcb-map",
    lambda i, j, k, l, J0: unred_from_red(J0, ring_pathB_red(i, j, k, l, J0)),
)

# Single stretched ME dump for the first nonzero channel
print("\n--- single ME dump (first |gold|>1e-8, i≠j,k≠l) ---")
found = False
nch = ms.GetNumberTwoBodyChannels()
for ch in range(nch):
    tbc = ms.GetTwoBodyChannel(ch)
    J0 = tbc.J
    nk = tbc.GetNumberKets()
    for ib in range(nk):
        bra = tbc.GetKet(ib)
        i, j = int(bra.p), int(bra.q)
        if i == j:
            continue
        for ik in range(ib, nk):
            ket = tbc.GetKet(ik)
            k, l = int(ket.p), int(ket.q)
            if k == l:
                continue
            g = Zg.TwoBody.GetTBME_J(J0, J0, i, j, k, l)
            if abs(g) < 1e-8:
                continue
            xr = ring_X_red(i, j, k, l, J0)
            xu = xr / hat(J0)
            an = amc_noperm_unred(i, j, k, l, J0)
            xb = ring_pathB_red(i, j, k, l, J0) / hat(J0)
            print(
                f"  ({i},{j}|{k},{l}) J={J0}\n"
                f"    CG gold unred = {g:.10e}\n"
                f"    ring X_red     = {xr:.10e}\n"
                f"    ring X_unred   = {xu:.10e}  r={xu/g if abs(g)>1e-14 else float('nan'):.6g}\n"
                f"    AMC noperm     = {an:.10e}  r={an/g if abs(g)>1e-14 else float('nan'):.6g}\n"
                f"    PathB red/Ĵ    = {xb:.10e}  r={xb/g if abs(g)>1e-14 else float('nan'):.6g}"
            )
            found = True
            break
        if found:
            break
    if found:
        break

cm.SetGIVcChiWhich(0)
cm.SetGIVcFoldAS(True)

ok = ok_ring and ok_b
print(
    "\nPASS — Path A + Path B bare ≡ ethS Path B"
    if ok
    else "\nFAIL — bare mid-J not locked"
)
if ok_ring and not ok_b:
    print("  (Path A OK; Path B still open)")
elif not ok_ring:
    print("  (Path A FAIL — unexpected)")
sys.exit(0 if ok else 1)
