#!/usr/bin/env python3
"""Ring Z = Σ_ab Ω_pbar Ω_aqsb: m-scheme oracle vs AMC Path A / Path B.

m-scheme (no occ). Two reduced tensors → scalar needs the rank coupling:

  Z(m) = Σ_{a ma b mb} CG(λ,μ; λ,−μ; 0,0) · Ω_{pb,ar}(m) · Ω_{aq,sb}(m)
       = Σ (−1)^{λ−μ}/λ̂ · Ω · Ω
  with μ = (mp+mb)−(ma+mr).

Then project to reduced scalar Z^J (AMC reduce=true):
  Z_red(J) = Σ_m CG_pq CG_sr Z(m) / Ĵ

J-scheme:
  Path A: z_pbar_aqsb_direct_ninej.tex
  Path B: z_pbar_aqsb_pandya_ninej.tex (composed)

Usage:
  PYTHONPATH=build python3 run/test_z_ring_mscheme_sign.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-7
seed = 7

ms = ModelSpace(emax, "He4", "He4")
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
if not Omega.IsReduced():
    Omega.MakeReduced()

orbits = list(ms.all_orbits)
print(
    f"Ω: λ={Omega.GetJRank()} reduced={Omega.IsReduced()} "
    f"||2b||={Omega.TwoBodyNorm():.6g} emax={emax} seed={seed}"
)


def phase(x: float) -> float:
    return 1.0 if int(round(2 * x)) % 4 == 0 else -1.0


def hat(J: float) -> float:
    return math.sqrt(2 * J + 1)


def jo(a: int) -> float:
    return ms.GetOrbit(a).j2 * 0.5


def j2(a: int) -> int:
    return ms.GetOrbit(a).j2


def tri(a, b, c) -> bool:
    return abs(a - b) <= c <= a + b


def om_tbme(J1, J2, a, b, c, d) -> float:
    return Omega.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


# ---------------------------------------------------------------------------
# m-scheme oracle
# ---------------------------------------------------------------------------
def z_mscheme(p, mp, q, mq, s, ms_, r, mr) -> float:
    if p == q and mp == mq:
        return 0.0
    if s == r and ms_ == mr:
        return 0.0
    sm = 0.0
    for a in orbits:
        oa = ms.GetOrbit(a)
        for ma in range(-oa.j2, oa.j2 + 1, 2):
            for b in orbits:
                ob = ms.GetOrbit(b)
                for mb in range(-ob.j2, ob.j2 + 1, 2):
                    o1 = ut.GetMschemeMatrixElement_2b(
                        Omega, p, mp, b, mb, a, ma, r, mr
                    )
                    if abs(o1) < 1e-16:
                        continue
                    o2 = ut.GetMschemeMatrixElement_2b(
                        Omega, a, ma, q, mq, s, ms_, b, mb
                    )
                    if abs(o2) < 1e-16:
                        continue
                    mu = (mp + mb - ma - mr) // 2
                    if abs(mu) > lam:
                        continue
                    # [Ω^λ × Ω^λ]_0 : CG(λ μ, λ −μ; 0 0) = (−1)^{λ−μ}/λ̂
                    w = CG(lam, mu, lam, -mu, 0, 0)
                    sm += w * o1 * o2
    return sm


def z_reduced_from_m(p, q, s, r, J0) -> float:
    """Reduced scalar TBME (AMC reduce=true / IsReduced)."""
    op, oq, os, or_ = (ms.GetOrbit(x) for x in (p, q, s, r))
    jp, jq, js, jr = (0.5 * o.j2 for o in (op, oq, os, or_))
    sm = 0.0
    for mp in range(-op.j2, op.j2 + 1, 2):
        for mq in range(-oq.j2, oq.j2 + 1, 2):
            if p == q and mp == mq:
                continue
            M = (mp + mq) // 2
            if abs(M) > J0:
                continue
            cg_pq = CG(jp, 0.5 * mp, jq, 0.5 * mq, J0, M)
            if abs(cg_pq) < 1e-15:
                continue
            for ms_ in range(-os.j2, os.j2 + 1, 2):
                for mr in range(-or_.j2, or_.j2 + 1, 2):
                    if s == r and ms_ == mr:
                        continue
                    if (ms_ + mr) // 2 != M:
                        continue
                    cg_sr = CG(js, 0.5 * ms_, jr, 0.5 * mr, J0, M)
                    if abs(cg_sr) < 1e-15:
                        continue
                    sm += cg_pq * cg_sr * z_mscheme(p, mp, q, mq, s, ms_, r, mr)
    return sm / hat(J0)


def z_unreduced_from_m(p, q, s, r, J0) -> float:
    return z_reduced_from_m(p, q, s, r, J0) / hat(J0)


# ---------------------------------------------------------------------------
# Path A — AMC direct (reduced Z)
# ---------------------------------------------------------------------------
def z_path_A(p, q, s, r, J0) -> float:
    jp, jq, js, jr = jo(p), jo(q), jo(s), jo(r)
    tot = 0.0
    for a in orbits:
        for b in orbits:
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2(p) - j2(b)) // 2, (j2(p) + j2(b)) // 2 + 1):
                for J3 in range(abs(j2(a) - j2(r)) // 2, (j2(a) + j2(r)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = om_tbme(J2, J3, p, b, a, r)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(abs(j2(a) - j2(q)) // 2, (j2(a) + j2(q)) // 2 + 1):
                        for J5 in range(abs(j2(s) - j2(b)) // 2, (j2(s) + j2(b)) // 2 + 1):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = om_tbme(J4, J5, a, q, s, b)
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
                                    / hat(lam)
                                    * s1
                                    * s2
                                    * n9
                                    * o1
                                    * o2
                                )
    return -phase(jp + js) * hat(J0) * tot


# ---------------------------------------------------------------------------
# Path B — AMC Pandya composed (reduced Z)
# ---------------------------------------------------------------------------
def pandya_bar_amc(i, j, k, l, Jbra, Jket) -> float:
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not tri(Jbra, Jket, lam):
        return 0.0
    sm = 0.0
    for J2 in range(abs(j2(i) - j2(j)) // 2, (j2(i) + j2(j)) // 2 + 1):
        if i == j and J2 % 2 != 0:
            continue
        for J3 in range(abs(j2(k) - j2(l)) // 2, (j2(k) + j2(l)) // 2 + 1):
            if k == l and J3 % 2 != 0:
                continue
            if not tri(J2, J3, lam):
                continue
            n9 = NineJ(lam, Jbra, Jket, J3, jl, jk, J2, ji, jj)
            if abs(n9) < 1e-16:
                continue
            sm += phase(J2) * hat(J2) * hat(J3) * n9 * om_tbme(J2, J3, i, j, k, l)
    return -phase(Jbra + ji + jk + lam) * hat(Jbra) * hat(Jket) * sm


bar_cache: dict = {}


def get_bar(i, j, k, l, Jbra, Jket):
    key = (i, j, k, l, Jbra, Jket)
    if key not in bar_cache:
        bar_cache[key] = pandya_bar_amc(i, j, k, l, Jbra, Jket)
    return bar_cache[key]


def barz_mid(p, q, s, r, J0) -> float:
    sm = 0.0
    for J2 in range(0, max(j2(x) for x in orbits) + lam + 2):
        if not tri(J0, J2, lam):
            continue
        for a in orbits:
            for b in orbits:
                sm += (
                    phase(J2 + lam)
                    / hat(lam)
                    * get_bar(p, b, a, r, J0, J2)
                    * get_bar(a, q, s, b, J2, J0)
                )
    return phase(J0) / hat(J0) * sm


def z_path_B(p, q, s, r, J0) -> float:
    jp, jq, js, jr = jo(p), jo(q), jo(s), jo(r)
    sm = 0.0
    Jpmax = max((j2(p) + j2(r)) // 2, (j2(q) + j2(s)) // 2) + 1
    for Jp in range(0, Jpmax + 1):
        six = SixJ(jr, js, J0, jq, jp, Jp)
        if abs(six) < 1e-16:
            continue
        sm += hat(Jp) * six * barz_mid(p, q, s, r, Jp)
    return -hat(J0) * sm


# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
def channel_ok(p, q, s, r, J0) -> bool:
    if p == q and J0 % 2 != 0:
        return False
    if s == r and J0 % 2 != 0:
        return False
    op, oq, os, or_ = (ms.GetOrbit(x) for x in (p, q, s, r))
    if (op.l + oq.l + os.l + or_.l) % 2 != 0:
        return False
    if abs((op.tz2 + oq.tz2) - (os.tz2 + or_.tz2)) != 0:
        return False
    if not (
        abs(j2(p) - j2(q)) // 2 <= J0 <= (j2(p) + j2(q)) // 2
        and abs(j2(s) - j2(r)) // 2 <= J0 <= (j2(s) + j2(r)) // 2
    ):
        return False
    return True


stats = {k: {"max": 0.0, "n": 0, "ratios": Counter()} for k in ("A", "B", "A_vs_B")}
worst = {"A": None, "B": None}

for p in orbits:
    for q in orbits:
        if q < p:
            continue
        for s in orbits:
            for r in orbits:
                if r < s:
                    continue
                Jmin = max(abs(j2(p) - j2(q)) // 2, abs(j2(s) - j2(r)) // 2)
                Jmax = min((j2(p) + j2(q)) // 2, (j2(s) + j2(r)) // 2)
                for J0 in range(Jmin, Jmax + 1):
                    if not channel_ok(p, q, s, r, J0):
                        continue
                    Zm = z_reduced_from_m(p, q, s, r, J0)
                    ZA = z_path_A(p, q, s, r, J0)
                    ZB = z_path_B(p, q, s, r, J0)

                    for tag, Zj in (("A", ZA), ("B", ZB)):
                        if abs(Zm) < 1e-12 and abs(Zj) < 1e-12:
                            continue
                        err = abs(Zm - Zj)
                        stats[tag]["n"] += 1
                        stats[tag]["max"] = max(stats[tag]["max"], err)
                        if abs(Zj) > 1e-12:
                            stats[tag]["ratios"][round(Zm / Zj, 6)] += 1
                        if worst[tag] is None or err > worst[tag][0]:
                            worst[tag] = (err, p, q, s, r, J0, Zm, Zj)

                    if abs(ZA) > 1e-12 or abs(ZB) > 1e-12:
                        stats["A_vs_B"]["n"] += 1
                        stats["A_vs_B"]["max"] = max(
                            stats["A_vs_B"]["max"], abs(ZA + ZB)
                        )
                        if abs(ZB) > 1e-12:
                            stats["A_vs_B"]["ratios"][round(ZA / ZB, 6)] += 1

print("\n=== reduced Z^J  (m-scheme with [Ω×Ω]_0) vs J-scheme ===")
for tag in ("A", "B"):
    st = stats[tag]
    print(
        f"Path {tag}: n={st['n']}  max|Z_m−Z_{tag}|={st['max']:.3e}  "
        f"Zm/Zj={st['ratios'].most_common(3)}"
    )
    if worst[tag]:
        err, p, q, s, r, J0, Zm, Zj = worst[tag]
        print(
            f"  worst: pqsr=({p},{q},{s},{r}) J0={J0} "
            f"Zm={Zm:.6e} Z{tag}={Zj:.6e} Δ={err:.3e}"
        )

st = stats["A_vs_B"]
print(
    f"\nPath A vs B: n={st['n']}  max|A+B|={st['max']:.3e}  "
    f"A/B={st['ratios'].most_common(3)}"
)

# reduce vs unreduce spot check
for p, q, s, r, J0 in [(0, 1, 0, 1, 1), (0, 0, 0, 0, 0)]:
    if not channel_ok(p, q, s, r, J0):
        continue
    Zr = z_reduced_from_m(p, q, s, r, J0)
    Zu = z_unreduced_from_m(p, q, s, r, J0)
    print(
        f"\nreduce: pqsr=({p},{q},{s},{r}) J0={J0}  "
        f"Z_red={Zr:.6e}  Z_unred={Zu:.6e}  Z_red/Ĵ={Zr/hat(J0):.6e}"
    )
    break

pass_A = stats["A"]["max"] < tol
pass_B = stats["B"]["max"] < tol

print("\nVerdict:")
if pass_A and not pass_B:
    print("  CORRECT: Path A (AMC direct).  Path B = −Path A (wrong overall sign).")
elif pass_B and not pass_A:
    print("  CORRECT: Path B (AMC Pandya).  Path A = −Path B (wrong overall sign).")
elif pass_A and pass_B:
    print("  Both match (unexpected if A=−B).")
else:
    print("  Neither matches by exact equality; see ratios above.")

sys.exit(0 if (pass_A or pass_B) else 1)
