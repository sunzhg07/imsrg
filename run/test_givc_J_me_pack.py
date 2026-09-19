#!/usr/bin/env python3
"""Quick J-scheme ME pack compare for GIVc.

Unfactored gold (ethS leftover packaging):
  S = Σ_m CG(ji mi jj mj; J M) CG(jk mk jl ml; J M) * Mscheme_fact_GIVc(...)
  tilde = S / (2J+1)   # = Z_red/Ĵ with Z_red = S/Ĵ; GetTBME_J is unreduced tilde

Factorized J: Zf.TwoBody.GetTBME_J(J, J, i, j, k, l)

Grid (convincing subset; full JT-parity m-gate is 30/30 elsewhere):
  λ∈{0,1,2,3} × T∈{0,1} × π∈{0,1}  +  λ=4 T=0 both π

PASS if max|Δ| < TOL * scale (default 1e-6).

Usage:
  PYTHONPATH=build OMP_NUM_THREADS=1 python3 -B run/test_givc_J_me_pack.py [emax=1]
  N2B=8 TOL=1e-6 python3 -B run/test_givc_J_me_pack.py 1
"""
from __future__ import annotations

import math
import os
import sys
import time
import traceback
from collections import Counter

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
tol = float(os.environ.get("TOL", "1e-6"))
n2b = int(os.environ.get("N2B", "8"))

# Explicit case list (do not expand to full T=2 / λ=4 grid).
CASES = []
for lam in (0, 1, 2, 3):
    for T in (0, 1):
        for parity in (0, 1):
            CASES.append((lam, T, parity))
for parity in (0, 1):
    CASES.append((4, 0, parity))

eth = Commutator.FactorizedDoubleCommutator_eths

_FLAGS = (
    "SetUse_TypeI_1b",
    "SetUse_TypeII_1b",
    "SetUse_TypeIII_1b",
    "SetUse_TypeIIIa_1b",
    "SetUse_TypeGI_2b",
    "SetUse_TypeGII_2b",
    "SetUse_TypeGIIIa_2b",
    "SetUse_TypeGIIIb_2b",
    "SetUse_TypeGIIIc_2b",
    "SetUse_TypeGIVa_2b",
    "SetUse_TypeGIVb_2b",
    "SetUse_TypeGIVc_2b",
)


def reset():
    eth.SetUse_1b_Intermediates(False)
    eth.SetUse_2b_Intermediates(False)
    for name in _FLAGS:
        fn = getattr(eth, name, None)
        if fn:
            fn(False)


def j2i(ms, a):
    return ms.GetOrbit(a).j2


def mrange(ms, o):
    return range(-j2i(ms, o), j2i(ms, o) + 1, 2)


def make_Z(ms):
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    if Z.IsReduced():
        Z.MakeNotReduced()
    Z *= 0.0
    return Z


def collect_tbmes(ms, Z, nkeep):
    """Largest |GetTBME_J| samples from Z channels."""
    pts = []
    nch = ms.GetNumberTwoBodyChannels()
    for ch in range(nch):
        tbc = ms.GetTwoBodyChannel(ch)
        nk = tbc.GetNumberKets()
        J = tbc.J
        for ib in range(nk):
            bra = tbc.GetKet(ib)
            i, j = int(bra.p), int(bra.q)
            for ik in range(nk):
                ket = tbc.GetKet(ik)
                k, l = int(ket.p), int(ket.q)
                zj = Z.TwoBody.GetTBME_J(J, J, i, j, k, l)
                pts.append((abs(zj), i, j, k, l, J, zj))
    pts.sort(key=lambda t: t[0], reverse=True)
    out, seen = [], set()
    for p in pts:
        key = p[1:6]
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
        if len(out) >= nkeep:
            break
    return out


def project_tilde(ut, Eta, Gamma, ms, i, j, k, l, J):
    """Unreduced tilde: S/(2J+1) with Km = Mscheme_fact_GIVc."""
    oi, oj = ms.GetOrbit(i), ms.GetOrbit(j)
    ok, ol = ms.GetOrbit(k), ms.GetOrbit(l)
    ji, jj = 0.5 * oi.j2, 0.5 * oj.j2
    jk, jl = 0.5 * ok.j2, 0.5 * ol.j2
    cache = {}
    sm = 0.0
    for mi in mrange(ms, i):
        for mj in mrange(ms, j):
            if i == j and mi == mj:
                continue
            twoM = mi + mj
            if abs(twoM) > 2 * J:
                continue
            cgij = CG(ji, 0.5 * mi, jj, 0.5 * mj, J, 0.5 * twoM)
            if abs(cgij) < 1e-12:
                continue
            for mk in mrange(ms, k):
                ml = twoM - mk
                if ml not in mrange(ms, l):
                    continue
                if k == l and mk == ml:
                    continue
                cgkl = CG(jk, 0.5 * mk, jl, 0.5 * ml, J, 0.5 * twoM)
                if abs(cgkl) < 1e-12:
                    continue
                key = (i, mi, j, mj, k, mk, l, ml)
                if key not in cache:
                    cache[key] = ut.Mscheme_fact_GIVc(Eta, Gamma, *key)
                sm += cgij * cgkl * cache[key]
    return sm / (2.0 * J + 1.0)


def ratio_bucket(r, J):
    if r is None or not math.isfinite(r):
        return "nan"
    hatJ = math.sqrt(2.0 * J + 1.0)
    targets = [
        ("+1", 1.0),
        ("−1", -1.0),
        ("+Ĵ", hatJ),
        ("−Ĵ", -hatJ),
        ("+1/Ĵ", 1.0 / hatJ),
        ("−1/Ĵ", -1.0 / hatJ),
        ("+Ĵ²", hatJ * hatJ),
        ("−Ĵ²", -hatJ * hatJ),
        ("+1/Ĵ²", 1.0 / (hatJ * hatJ)),
        ("−1/Ĵ²", -1.0 / (hatJ * hatJ)),
        ("+√2", math.sqrt(2.0)),
        ("−√2", -math.sqrt(2.0)),
        ("+1/√2", 1.0 / math.sqrt(2.0)),
        ("−1/√2", -1.0 / math.sqrt(2.0)),
    ]
    for name, target in targets:
        if abs(r - target) < 0.08 * max(1.0, abs(target)):
            return name
    return "other"


def close(a, b):
    d = abs(a - b)
    scale = max(abs(a), abs(b), 1.0)
    return d <= tol * scale or d <= 1e-8


def run_case(ms, ut, lam, T, parity, seed):
    ut.SetRandomSeed(seed)
    Eta = ut.RandomOp(ms, lam, T, parity, 2, -1)
    Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
    if Gamma.IsReduced():
        Gamma.MakeNotReduced()
    if not Eta.IsReduced():
        Eta.MakeReduced()

    Zf = make_Z(ms)
    t0 = time.time()
    reset()
    eth.comm223_232_GIVc(Eta, Gamma, Zf)
    dt_e = time.time() - t0

    pts = collect_tbmes(ms, Zf, n2b)
    n, nfail, maxd = 0, 0, 0.0
    buckets = Counter()
    worst = None
    t0 = time.time()
    for _a, i, j, k, l, J, zj in pts:
        # Skip near-zero samples when Zf is tiny overall
        if abs(zj) < 1e-14 and Zf.TwoBodyNorm() < 1e-10:
            continue
        zm = project_tilde(ut, Eta, Gamma, ms, i, j, k, l, J)
        n += 1
        d = abs(zj - zm)
        if d > maxd:
            maxd = d
            worst = (i, j, k, l, J, zj, zm)
        if not close(zj, zm):
            nfail += 1
            r = zj / zm if abs(zm) > 1e-14 else None
            buckets[ratio_bucket(r, J)] += 1
            print(
                f"    fail λ={lam} T={T} π={parity}  "
                f"ijkl=({i},{j},{k},{l}) J={J}  "
                f"GetTBME={zj:.6e}  tilde={zm:.6e}  "
                f"ratio={zj / zm if abs(zm) > 1e-14 else float('nan'):.6g}",
                flush=True,
            )
    dt_m = time.time() - t0
    ok = nfail == 0 and n > 0
    return {
        "ok": ok,
        "n": n,
        "nfail": nfail,
        "maxd": maxd,
        "dt_e": dt_e,
        "dt_m": dt_m,
        "nrm_f": Zf.TwoBodyNorm(),
        "buckets": buckets,
        "worst": worst,
    }


def main():
    ms = ModelSpace(emax, "He4", "He4")
    ms.SetHbarOmega(20.0)
    ms.PreCalculateSixJ()
    ms.PreCalculateNineJ()
    ut = UnitTest(ms)

    print(
        f"GIVc J-ME pack  emax={emax} He4  N2B={n2b}  tol={tol}\n"
        f"  gold: tilde = Σ CG CG Mscheme_fact_GIVc / (2J+1)\n"
        f"  code: GetTBME_J(ethS GIVc)  [unreduced tilde; √2 restored]\n"
        f"  cases: {len(CASES)}  "
        f"(λ=0..3 × T=0,1 × π=0,1 + λ=4 T=0 both π)",
        flush=True,
    )
    print(
        f"{'λ':>2} {'T':>2} {'π':>2} {'seed':>5}  {'gate':6s}  "
        f"{'nOK/n':>8s}  {'maxΔ':>10s}  {'||Zf||':>10s}  dt",
        flush=True,
    )

    n_pass = n_tot = n_err = 0
    failures = []
    pack_hint = Counter()
    t_all = time.time()

    for lam, T, parity in CASES:
        n_tot += 1
        seed = 11 + 100 * lam + 10 * T + parity
        try:
            r = run_case(ms, ut, lam, T, parity, seed)
            if r["ok"]:
                n_pass += 1
                gate = "PASS"
            else:
                gate = "FAIL"
                failures.append((lam, T, parity, seed, r))
                pack_hint.update(r["buckets"])

            print(
                f"{lam:2d} {T:2d} {parity:2d} {seed:5d}  {gate:6s}  "
                f"{r['n']-r['nfail']:3d}/{r['n']:<3d}  "
                f"{r['maxd']:10.3e}  {r['nrm_f']:10.4g}  "
                f"ethS {r['dt_e']:.2f}s fold {r['dt_m']:.2f}s",
                flush=True,
            )
        except Exception as ex:
            n_err += 1
            failures.append((lam, T, parity, seed, None))
            print(
                f"{lam:2d} {T:2d} {parity:2d} {seed:5d}  ERROR  "
                f"{type(ex).__name__}: {ex}",
                flush=True,
            )
            traceback.print_exc()

    dt = time.time() - t_all
    print(flush=True)
    print(
        f"SUMMARY  J-ME pack: {n_pass}/{n_tot} PASS  errors={n_err}  ({dt:.1f}s)",
        flush=True,
    )
    if failures:
        print("FAILURES:", flush=True)
        for lam, T, parity, seed, r in failures:
            if r is None:
                print(f"  λ={lam} T={T} π={parity} seed={seed}  ERROR", flush=True)
            else:
                print(
                    f"  λ={lam} T={T} π={parity} seed={seed}  "
                    f"fail={r['nfail']}/{r['n']} maxΔ={r['maxd']:.3e}  "
                    f"ratio buckets={dict(r['buckets'])}",
                    flush=True,
                )
        if pack_hint:
            print(
                f"  packaging hint (fail ratios): {dict(pack_hint)}",
                flush=True,
            )
            print(
                "  If ratios cluster near Ĵ^{±1}/√2: fix gold fold once "
                "(do not retune ethS).",
                flush=True,
            )
    else:
        print("FAILURES: (none)", flush=True)

    verdict = n_pass == n_tot and n_err == 0
    print(
        "\nVERDICT: J-scheme (m-gold packed) "
        + ("≡" if verdict else "≢")
        + " ethS GetTBME_J across sampled grid  →  "
        + ("PASS" if verdict else "FAIL"),
        flush=True,
    )
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
