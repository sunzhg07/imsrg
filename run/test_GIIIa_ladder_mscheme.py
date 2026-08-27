#!/usr/bin/env python3
"""Γ^{III_a} ladder: m ≡ AMC from_chi ≡ χ^η × Γ (channel DGEMM).

Physics (analyze eq:GIIIa):
  Γ^{III_a} = −Σ_ab { (1−Pij) χ^η_ijab Γ_abkl + (1−Pkl) Γ_ijab χ^η_klab }

χ^η locked: run/test_chi_eta_mscheme.py (m ≡ AMC direct ≡ Pandya).
AMC ladder: learn/amc_tts/factored_GIIIa/input/G3a_from_chi.txt
  reduce=true → Z_red = −χ_red Γ_unred (no Ĵ between factors).

DGEMM (per J): Z = −Chi_AS @ Gam − Gam @ Chi_AS.T
  Chi_AS[ij,ab] = χ_ijab − (−1)^{J+ji+jj} χ_jiab

Gold is m / from_chi / χ×Γ — not ethS TTS Path A strips (different χ routing).

Usage:
  PYTHONPATH=build python3 -B run/test_GIIIa_ladder_mscheme.py [emax=1] [lambda=2]
"""

from __future__ import annotations

import math
import sys
import time
from collections import Counter

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lam = int(sys.argv[2]) if len(sys.argv) > 2 else 2
tol = 1e-6
seed = 11
# m projection is cheap at emax=1; set >0 to subsample largest |Z| only
n_m_sample = 0

ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = UnitTest(ms)
ut.SetRandomSeed(seed)

Eta = ut.RandomOp(ms, lam, 0, 0, 2, -1)
if not Eta.IsReduced():
    Eta.MakeReduced()
Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
if Gamma.IsReduced():
    Gamma.MakeNotReduced()

orbits = list(ms.all_orbits)
max_J = max(ms.GetOrbit(o).j2 for o in orbits)
print(
    f"emax={emax} λ={lam} seed={seed} norb={len(orbits)}\n"
    f"  Ω reduced={Eta.IsReduced()} ||2b||={Eta.TwoBodyNorm():.4g}\n"
    f"  Γ reduced={Gamma.IsReduced()} ||2b||={Gamma.TwoBodyNorm():.4g}"
)


def phase(x: float) -> float:
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


def mrange(o: int):
    return range(-j2i(o), j2i(o) + 1, 2)


def w_eta(a, b, k) -> float:
    return nbar(a) * occ(b) * nbar(k) + occ(a) * nbar(b) * occ(k)


def tbme(J1, J2, a, b, c, d) -> float:
    return Eta.TwoBody.GetTBME_J(J1, J2, a, b, c, d)


def gam_J(J0, a, b, c, d) -> float:
    return Gamma.TwoBody.GetTBME_J(J0, J0, a, b, c, d)


# ---------------------------------------------------------------------------
# χ^η reduced — AMC Path A (locked ≡ m)
# ---------------------------------------------------------------------------
def chi_eta_red(i, j, k, l, J0) -> float:
    ji, jj, jk, jl = jo(i), jo(j), jo(k), jo(l)
    if not (tri(ji, jj, J0) and tri(jk, jl, J0)):
        return 0.0
    tot = 0.0
    for a in orbits:
        for b in orbits:
            w = w_eta(a, b, k)
            if abs(w) < 1e-12:
                continue
            ja, jb = jo(a), jo(b)
            for J2 in range(abs(j2i(i) - j2i(a)) // 2, (j2i(i) + j2i(a)) // 2 + 1):
                for J3 in range(abs(j2i(b) - j2i(l)) // 2, (j2i(b) + j2i(l)) // 2 + 1):
                    if not tri(J2, J3, lam):
                        continue
                    o1 = tbme(J2, J3, i, a, b, l)
                    if abs(o1) < 1e-16:
                        continue
                    for J4 in range(
                        abs(j2i(b) - j2i(j)) // 2, (j2i(b) + j2i(j)) // 2 + 1
                    ):
                        for J5 in range(
                            abs(j2i(k) - j2i(a)) // 2, (j2i(k) + j2i(a)) // 2 + 1
                        ):
                            if not tri(J4, J5, lam):
                                continue
                            o2 = tbme(J4, J5, b, j, k, a)
                            if abs(o2) < 1e-16:
                                continue
                            j0max = int(max(J2, J3, J4, J5, ji, jk, ja, lam) + 2)
                            for j0_2 in range(0, 2 * j0max + 1):
                                j0 = 0.5 * j0_2
                                s1 = SixJ(J3, lam, J2, ja, ji, j0)
                                s2 = SixJ(J4, lam, J5, ja, jk, j0)
                                n9 = NineJ(jl, jb, J3, jk, J4, j0, J0, jj, ji)
                                if abs(s1 * s2 * n9) < 1e-16:
                                    continue
                                tot += (
                                    w
                                    * phase(J2 + J4 + lam)
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
    return -phase(ji + jk) * hat(J0) * tot


print("Building χ^η_red (Path A) ...")
t0 = time.time()
chi_red: dict = {}
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
                        continue
                    v = chi_eta_red(i, j, k, l, J0)
                    if abs(v) > 1e-16:
                        chi_red[(i, j, k, l, J0)] = v
print(f"  nonzero={len(chi_red)}  ({time.time()-t0:.2f}s)")


def chi_m_from_red(i, mi, j, mj, k, mk, l, ml) -> float:
    """Reduced scalar unpack: Σ_J CG CG χ_red / Ĵ."""
    M0 = 0.5 * (mi + mj)
    sm = 0.0
    for J0 in range(abs(j2i(i) - j2i(j)) // 2, (j2i(i) + j2i(j)) // 2 + 1):
        if abs(M0) > J0:
            continue
        if not tri(jo(k), jo(l), J0):
            continue
        cab = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M0)
        if abs(cab) < 1e-15:
            continue
        ccd = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M0)
        if abs(ccd) < 1e-15:
            continue
        v = chi_red.get((i, j, k, l, J0), 0.0)
        sm += cab * ccd / hat(J0) * v
    return sm


def gam_m(i, mi, j, mj, k, mk, l, ml) -> float:
    return ut.GetMschemeMatrixElement_2b(Gamma, i, mi, j, mj, k, mk, l, ml)


# ---------------------------------------------------------------------------
# AMC ladder (G3a_from_chi.tex, reduce=true)
# ---------------------------------------------------------------------------
def amc_ladder_red(i, j, k, l, J0) -> float:
    if not (tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)):
        return 0.0
    pij = phase(J0 + jo(i) + jo(j))
    pkl = phase(J0 + jo(k) + jo(l))
    sm = 0.0
    for a in orbits:
        for b in orbits:
            gab = gam_J(J0, a, b, k, l)
            if abs(gab) > 1e-16:
                sm += -chi_red.get((i, j, a, b, J0), 0.0) * gab
                sm += pij * chi_red.get((j, i, a, b, J0), 0.0) * gab
            gij = gam_J(J0, i, j, a, b)
            if abs(gij) > 1e-16:
                sm += -gij * chi_red.get((k, l, a, b, J0), 0.0)
                sm += pkl * gij * chi_red.get((l, k, a, b, J0), 0.0)
    return sm


# ---------------------------------------------------------------------------
# Channel DGEMM: Z = −Chi_AS @ Gam − Gam @ Chi_AS.T
# ---------------------------------------------------------------------------
def matmul(A, B):
    """Pure-Python C = A @ B (lists of lists)."""
    n, m, p = len(A), len(A[0]), len(B[0])
    C = [[0.0] * p for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for k in range(m):
            aik = Ai[k]
            if abs(aik) < 1e-16:
                continue
            Bk = B[k]
            Ci = C[i]
            for j in range(p):
                Ci[j] += aik * Bk[j]
    return C


def transpose(A):
    return [list(row) for row in zip(*A)]


print("Channel DGEMM χ^η × Γ ...")
t0 = time.time()
z_dgemm: dict = {}
for J0 in range(0, max_J + 1):
    pairs = [
        (p, q)
        for p in orbits
        for q in orbits
        if tri(jo(p), jo(q), J0)
    ]
    if not pairs:
        continue
    n = len(pairs)
    idx = {pq: i for i, pq in enumerate(pairs)}
    Chi = [[0.0] * n for _ in range(n)]
    Gam = [[0.0] * n for _ in range(n)]
    for ia, (i, j) in enumerate(pairs):
        for ib, (a, b) in enumerate(pairs):
            Chi[ia][ib] = chi_red.get((i, j, a, b, J0), 0.0)
            Gam[ia][ib] = gam_J(J0, i, j, a, b)
    # bra AS: Chi_AS[ij,ab] = χ_ijab − (−1)^{J+ji+jj} χ_jiab
    Chi_AS = [[0.0] * n for _ in range(n)]
    for ia, (i, j) in enumerate(pairs):
        pij = phase(J0 + jo(i) + jo(j))
        ji_idx = idx.get((j, i))
        for ib in range(n):
            if ji_idx is None:
                Chi_AS[ia][ib] = Chi[ia][ib]
            else:
                Chi_AS[ia][ib] = Chi[ia][ib] - pij * Chi[ji_idx][ib]
    # Z = −Chi_AS @ Gam − Gam @ Chi_AS.T
    T1 = matmul(Chi_AS, Gam)
    T2 = matmul(Gam, transpose(Chi_AS))
    for ia, (i, j) in enumerate(pairs):
        for ib, (k, l) in enumerate(pairs):
            v = -(T1[ia][ib] + T2[ia][ib])
            if abs(v) > 1e-16:
                z_dgemm[(i, j, k, l, J0)] = v
print(f"  nonzero={len(z_dgemm)}  ({time.time()-t0:.2f}s)")


# ---------------------------------------------------------------------------
# AMC ≡ DGEMM (exhaustive)
# ---------------------------------------------------------------------------
print("\nAMC ≡ DGEMM (exhaustive) ...")
t0 = time.time()
rats_ad = Counter()
n_ad = 0
max_ad = 0.0
channels = []
for i in orbits:
    for j in orbits:
        for k in orbits:
            for l in orbits:
                for J0 in range(0, max_J + 1):
                    if not (
                        tri(jo(i), jo(j), J0) and tri(jo(k), jo(l), J0)
                    ):
                        continue
                    aR = amc_ladder_red(i, j, k, l, J0)
                    dR = z_dgemm.get((i, j, k, l, J0), 0.0)
                    if abs(aR) < 1e-12 and abs(dR) < 1e-12:
                        continue
                    n_ad += 1
                    max_ad = max(max_ad, abs(aR - dR))
                    if abs(dR) > 1e-8:
                        rats_ad[round(aR / dR, 6)] += 1
                    channels.append((abs(aR), i, j, k, l, J0, aR))

print(
    f"  n={n_ad}  max|AMC−DGEMM|={max_ad:.3e}  "
    f"AMC/DGEMM={rats_ad.most_common(3)}  ({time.time()-t0:.2f}s)"
)


# ---------------------------------------------------------------------------
# m-scheme ladder → Z_red = S/Ĵ  (sampled largest |AMC|)
# ---------------------------------------------------------------------------
def Z_m(i, mi, j, mj, k, mk, l, ml) -> float:
    sm = 0.0
    for a in orbits:
        for ma in mrange(a):
            for b in orbits:
                for mb in mrange(b):
                    g = gam_m(a, ma, b, mb, k, mk, l, ml)
                    if abs(g) > 1e-16:
                        sm += chi_m_from_red(i, mi, j, mj, a, ma, b, mb) * g
                        sm -= chi_m_from_red(j, mj, i, mi, a, ma, b, mb) * g
                    g2 = gam_m(i, mi, j, mj, a, ma, b, mb)
                    if abs(g2) > 1e-16:
                        sm += g2 * chi_m_from_red(k, mk, l, ml, a, ma, b, mb)
                        sm -= g2 * chi_m_from_red(l, ml, k, mk, a, ma, b, mb)
    return -sm


def project_red(i, j, k, l, J0) -> float:
    sm = 0.0
    for mi in mrange(i):
        for mj in mrange(j):
            M = 0.5 * (mi + mj)
            if abs(M) > J0:
                continue
            c1 = CG(jo(i), mi * 0.5, jo(j), mj * 0.5, J0, M)
            if abs(c1) < 1e-15:
                continue
            for mk in mrange(k):
                ml = mi + mj - mk
                if ml not in mrange(l):
                    continue
                c2 = CG(jo(k), mk * 0.5, jo(l), ml * 0.5, J0, M)
                if abs(c2) < 1e-15:
                    continue
                sm += c1 * c2 * Z_m(i, mi, j, mj, k, mk, l, ml)
    return sm / hat(J0)


channels.sort(reverse=True)
sample = channels if n_m_sample <= 0 else channels[:n_m_sample]
print(
    f"\nm ≡ AMC/DGEMM "
    f"({'exhaustive' if n_m_sample <= 0 else f'top-{len(sample)} |Z|'}) ..."
)
t0 = time.time()
rats_ma = Counter()
rats_md = Counter()
n_m = 0
max_ma = 0.0
max_md = 0.0
fails = []
for _, i, j, k, l, J0, aR in sample:
    mR = project_red(i, j, k, l, J0)
    dR = z_dgemm.get((i, j, k, l, J0), 0.0)
    n_m += 1
    max_ma = max(max_ma, abs(mR - aR))
    max_md = max(max_md, abs(mR - dR))
    if abs(aR) > 1e-8:
        rats_ma[round(mR / aR, 4)] += 1
    if abs(dR) > 1e-8:
        rats_md[round(mR / dR, 4)] += 1
    if abs(mR - aR) > 1e-5:
        fails.append(((i, j, k, l, J0), mR, aR, dR))

print(
    f"  n={n_m}  max|m−AMC|={max_ma:.3e}  m/AMC={rats_ma.most_common(3)}\n"
    f"       max|m−DGEMM|={max_md:.3e}  m/DGEMM={rats_md.most_common(3)}\n"
    f"       ({time.time()-t0:.2f}s)"
)
if fails[:3]:
    print("  sample fails:", fails[:3])

ok_ad = n_ad > 0 and max_ad < tol and (
    not rats_ad or rats_ad.most_common(1)[0][0] == 1.0
)
ok_ma = n_m > 0 and max_ma < tol and (
    not rats_ma or rats_ma.most_common(1)[0][0] == 1.0
)
ok_md = n_m > 0 and max_md < tol and (
    not rats_md or rats_md.most_common(1)[0][0] == 1.0
)

print(
    f"\nm ≡ AMC from_chi:  {'PASS' if ok_ma else 'FAIL'}\n"
    f"m ≡ χ^η×Γ DGEMM:  {'PASS' if ok_md else 'FAIL'}\n"
    f"AMC ≡ DGEMM:      {'PASS' if ok_ad else 'FAIL'}"
)
ok = ok_ad and ok_ma and ok_md
print("\nPASS — Γ^{III_a} m ≡ AMC ≡ χ^η×Γ DGEMM" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
