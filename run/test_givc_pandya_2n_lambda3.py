#!/usr/bin/env python3
"""GIVc Path B step-1 lock (λ=3): fwd Pandya of χ^λ and Ω into 2n×2n CC.

Packaging (named first):
  χ^λ  : tensor → WE-reduced (ChiTab = GetTBME_J units; always reduced)
  Ω    : tensor → WE-reduced (Eta.IsReduced() True)
  barχ, barΩ : tensors in Pandya/CC → also WE-reduced (AMC Eq1/2)

AMC Eq1/2 (G4c_chi_omega_pandya_ninej.tex), same-label ME:
  barO_ijkl^{J0 J1 λ} = −(−1)^{J0+ji+jk+λ} Ĵ0 Ĵ1
      Σ_{J2 J3} (−1)^{J2} Ĵ2 Ĵ3  NineJ(λ J0 J1; J3 jl jk; J2 ji jj) O_ijkl^{J2 J3}

IMSRG DoTensorPandya is adcb: bar(a,b,c,d) ← O(a,d,c,b). Do not mix without a map.

Usage:
  PYTHONPATH=build python3 -B run/test_givc_pandya_2n_lambda3.py [emax=1]
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
lam = 3
seed = 17
tol = 1e-8

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
    f"  Ω IsReduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ IsReduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}\n"
    f"  χ^λ packaging: WE-reduced (tensor; ChiTab ≡ GetTBME_J)"
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


print("χ^λ cache (ordinary WE-reduced) ...")
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


def amc_same_label(get_me, Jbra, Jket, i, j, k, l) -> float:
    """AMC Eq1/2: barO_ijkl from O_ijkl (same labels).

    scheme=((1,-4),(3,-2)) ⇒ LHS couples (iℓ)(kj):
      Triangle(ji,jl,Jbra) and Triangle(jk,jj,Jket).
    RHS ordinary ME still uses (ij)(kl) with J2,J3.
    """
    if not tri(Jbra, Jket, lam):
        return 0.0
    # LHS Pandya channel triangles — NOT ordinary (ij)(kl)
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


def imsrg_adcb(get_me, Jbra, Jket, a, b, c, d) -> float:
    """IMSRG DoTensorPandya(adcb): bar(a,b,c,d) ← O(a,d,c,b)."""
    if not tri(Jbra, Jket, lam):
        return 0.0
    if not (tri(jo(a), jo(b), Jbra) and tri(jo(c), jo(d), Jket)):
        return 0.0
    sm = 0.0
    j1min = abs(j2i(a) - j2i(d)) // 2
    j1max = (j2i(a) + j2i(d)) // 2
    for J1 in range(j1min, j1max + 1):
        j2min = max(abs(j2i(c) - j2i(b)) // 2, abs(J1 - lam))
        j2max = min((j2i(c) + j2i(b)) // 2, J1 + lam)
        for J2 in range(j2min, j2max + 1):
            n9 = NineJ(
                jo(a), jo(d), J1, jo(b), jo(c), J2, Jbra, Jket, lam
            )
            if abs(n9) < 1e-14:
                continue
            me = get_me(J1, J2, a, d, c, b)
            if abs(me) < 1e-16:
                continue
            hats = hat(J1) * hat(J2) * hat(Jbra) * hat(Jket)
            sm -= hats * iphase((j2i(b) + j2i(d)) // 2 + Jket + J2) * n9 * me
    return sm


def get_chi(J2, J3, a, b, c, d):
    return chiJ.get((J2, J3, a, b, c, d), 0.0)


def get_om(J2, J3, a, b, c, d):
    return Eta.TwoBody.GetTBME_J(J2, J3, a, b, c, d)


def cmp_maps(tag, amc_fn, maps):
    """Compare AMC same-label bar(i,j,k,l) to each candidate map."""
    print(f"\n=== {tag}: AMC same-label vs candidates (λ={lam}) ===")
    best = None
    for name, fn in maps:
        rats = Counter()
        n = 0
        max_abs = 0.0
        worst = None
        for i in orbits:
            for j in orbits:
                for k in orbits:
                    for l in orbits:
                        for J0 in range(0, max_J + 1):
                            for J1 in range(0, max_J + 1):
                                if not tri(J0, J1, lam):
                                    continue
                                a = amc_fn(J0, J1, i, j, k, l)
                                b = fn(J0, J1, i, j, k, l)
                                if abs(a) < 1e-10 and abs(b) < 1e-10:
                                    continue
                                n += 1
                                d = abs(a - b)
                                if d > max_abs:
                                    max_abs = d
                                    worst = (i, j, k, l, J0, J1, a, b)
                                if abs(a) > 1e-10:
                                    rats[round(b / a, 4)] += 1
                                if n >= 200:
                                    break
                            if n >= 200:
                                break
                        if n >= 200:
                            break
                    if n >= 200:
                        break
                if n >= 200:
                    break
            if n >= 200:
                break
        top = rats.most_common(3)
        ok = n > 0 and max_abs < tol and (top and top[0][0] == 1.0)
        print(
            f"  {name:28s} n={n:4d} max|Δ|={max_abs:.3e} B/A={top} "
            f"=> {'PASS' if ok else 'FAIL'}"
        )
        if worst and not ok:
            i, j, k, l, J0, J1, a, b = worst
            print(
                f"    worst ({i},{j},{k},{l})J={J0},{J1}: "
                f"AMC={a:.6e} cand={b:.6e} r={b/a if abs(a)>1e-14 else float('nan'):g}"
            )
        if ok:
            best = name
    return best


# Candidates: how to get AMC bar(i,j,k,l) from IMSRG adcb
om_maps = [
    ("adcb(i,j,k,l)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_om, J0, J1, i, j, k, l)),
    ("adcb(i,l,k,j)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_om, J0, J1, i, l, k, j)),
    ("adcb(i,l,j,k)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_om, J0, J1, i, l, j, k)),
    ("adcb(i,j,l,k)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_om, J0, J1, i, j, l, k)),
    (
        "NOTES AMC(p,b,a,r)=adcb(p,r,a,b)",
        lambda J0, J1, i, j, k, l: imsrg_adcb(get_om, J0, J1, i, l, k, j),
    ),
]
chi_maps = [
    ("adcb(i,j,k,l)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_chi, J0, J1, i, j, k, l)),
    ("adcb(i,l,k,j)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_chi, J0, J1, i, l, k, j)),
    ("adcb(i,l,j,k)", lambda J0, J1, i, j, k, l: imsrg_adcb(get_chi, J0, J1, i, l, j, k)),
]

# Also: ethS helpers aimed at χ_ialb / Ω_bjak product topology
def pandya_chi_ialb_as_bar_il_ab(J0, J1, i, l, a, b):
    """Current ethS PandyaChiIalb intent: χ̄(il;ab) from χ(i,a,l,b)."""
    # Rebuild with same formula as C++ PandyaChiIalb
    if not tri(J0, J1, lam):
        return 0.0
    if not (tri(jo(i), jo(l), J0) and tri(jo(a), jo(b), J1)):
        return 0.0
    sm = 0.0
    for Jia in range(abs(j2i(i) - j2i(a)) // 2, (j2i(i) + j2i(a)) // 2 + 1):
        for Jlb in range(
            max(abs(j2i(l) - j2i(b)) // 2, abs(Jia - lam)),
            min((j2i(l) + j2i(b)) // 2, Jia + lam) + 1,
        ):
            n9 = NineJ(jo(i), jo(a), Jia, jo(l), jo(b), Jlb, J0, J1, lam)
            if abs(n9) < 1e-14:
                continue
            me = get_chi(Jia, Jlb, i, a, l, b)
            if abs(me) < 1e-16:
                continue
            hats = hat(Jia) * hat(Jlb) * hat(J0) * hat(J1)
            sm -= hats * iphase((j2i(l) + j2i(b)) // 2 + J1 + Jlb) * n9 * me
    return sm


print("\n--- Map search: AMC barΩ_ijkl vs IMSRG adcb ---")
best_om = cmp_maps(
    "Ω",
    lambda J0, J1, i, j, k, l: amc_same_label(get_om, J0, J1, i, j, k, l),
    om_maps,
)

print("\n--- Map search: AMC barχ_ijkl vs IMSRG adcb ---")
best_chi = cmp_maps(
    "χ",
    lambda J0, J1, i, j, k, l: amc_same_label(get_chi, J0, J1, i, j, k, l),
    chi_maps,
)

# Product topology: AMC uses barChi_ialb (= χ̄ scheme-label i,a,l,b ⇒ CC (ib)(la)?
# scheme (1,-4),(3,-2) on (i,a,l,b): bra=(i,-b), ket=(l,-a) → (ib)(la)
# Wanted for DGEMM: χ̄(il;ab). Map: AMC(i,j,k,l) = IMSRG(i,l,k,j)
#   so χ̄_IMSRG(il;ab) = AMC(i,b,a,l)   [j=b,k=a,l=l]
print("\n--- Product topology: χ̄(il;ab) from AMC map vs PandyaChiIalb ---")
rats = Counter()
n = 0
max_abs = 0.0
worst = None
for i in orbits:
    for l in orbits:
        for a in orbits:
            for b in orbits:
                for J0 in range(0, max_J + 1):
                    for J1 in range(0, max_J + 1):
                        if not tri(J0, J1, lam):
                            continue
                        if not (tri(jo(i), jo(l), J0) and tri(jo(a), jo(b), J1)):
                            continue
                        # AMC(i,b,a,l) ↔ IMSRG(i,l,a,b)
                        amc = amc_same_label(get_chi, J0, J1, i, b, a, l)
                        eth = pandya_chi_ialb_as_bar_il_ab(J0, J1, i, l, a, b)
                        ad = imsrg_adcb(get_chi, J0, J1, i, l, a, b)
                        if abs(amc) < 1e-10 and abs(eth) < 1e-10 and abs(ad) < 1e-10:
                            continue
                        n += 1
                        d = abs(amc - ad)
                        if d > max_abs:
                            max_abs = d
                            worst = (i, l, a, b, J0, J1, amc, eth, ad)
                        if abs(amc) > 1e-10:
                            rats[("ad/AMC", round(ad / amc, 4))] += 1
                            rats[("eth/AMC", round(eth / amc, 4))] += 1
                        if n >= 200:
                            break
                    if n >= 200:
                        break
                if n >= 200:
                    break
            if n >= 200:
                break
        if n >= 200:
            break
    if n >= 200:
        break
print(f"  n={n} max|AMC−adcb|={max_abs:.3e}")
print(f"  ratios: {rats.most_common(8)}")
if worst:
    i, l, a, b, J0, J1, amc, eth, ad = worst
    print(
        f"  worst (il;ab)=({i},{l};{a},{b}) J={J0},{J1}: "
        f"AMC={amc:.6e} adcb={ad:.6e} ChiIalb={eth:.6e}"
    )

# 2n×2n pack in IMSRG adcb layout: bar(p,q,r,s) = AMC(p,s,r,q)
print("\n--- 2n×2n CC pack: AMC map ≡ adcb (Ω, λ=3) ---")
n_cmp = 0
n_bad = 0
maxd = 0.0
n_cc = ms.GetNumberTwoBodyChannels_CC()
for ch_b in range(n_cc):
    tb = ms.GetTwoBodyChannel_CC(ch_b)
    nb = tb.GetNumberKets()
    if nb < 1:
        continue
    for ch_k in range(n_cc):
        tk = ms.GetTwoBodyChannel_CC(ch_k)
        nk = tk.GetNumberKets()
        if nk < 1 or not tri(tb.J, tk.J, lam):
            continue
        if (tb.parity + tk.parity) % 2 != Eta.GetParity():
            continue
        for ib in range(2 * nb):
            if ib < nb:
                ket = tb.GetKet(ib)
                p, q = int(ket.p), int(ket.q)
            else:
                ket = tb.GetKet(ib - nb)
                if ket.p == ket.q:
                    continue
                p, q = int(ket.q), int(ket.p)
            for ik in range(2 * nk):
                if ik < nk:
                    ket = tk.GetKet(ik)
                    r, s = int(ket.p), int(ket.q)
                else:
                    ket = tk.GetKet(ik - nk)
                    if ket.p == ket.q:
                        continue
                    r, s = int(ket.q), int(ket.p)
                # IMSRG bar(p,q,r,s) ↔ AMC bar(p,s,r,q)
                amc = amc_same_label(get_om, tb.J, tk.J, p, s, r, q)
                ad = imsrg_adcb(get_om, tb.J, tk.J, p, q, r, s)
                if abs(amc) < 1e-12 and abs(ad) < 1e-12:
                    continue
                n_cmp += 1
                d = abs(amc - ad)
                if d > maxd:
                    maxd = d
                if d > tol:
                    n_bad += 1
print(
    f"  n={n_cmp} nbad={n_bad} max|Δ|={maxd:.3e}  => "
    f"{'PASS' if n_bad == 0 and n_cmp > 0 else 'FAIL'}"
)

print("\nSummary:")
print(f"  best Ω map : {best_om}")
print(f"  best χ map : {best_chi}")
print(
    "  Locked packaging: χ/Ω WE-reduced; AMC bar(i,j,k,l)=IMSRG adcb(i,l,k,j).\n"
    "  2n×2n store IMSRG layout; DGEMM uses barχ(il;ab)·barΩ(ab;kj)."
)
ok = best_om is not None and best_chi is not None and n_bad == 0 and n_cmp > 0
print("PASS — Pandya map locked" if ok else "FAIL — need map / packaging fix")
sys.exit(0 if ok else 1)
