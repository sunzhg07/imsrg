#!/usr/bin/env python3
"""Final GIVc grid: ethS factorized J vs m-scheme gold, over (λ, T, π).

Primary PASS gate (same as bench_eths_pathB_vs_mscheme):
  GetMscheme(Zf) ≡ UnitTest.Mscheme_fact_GIVc
  where Zf = ethS.comm223_232_GIVc (unreduced scalar leftover).

DIAG (not a gate; tts_GIVc ≠ m historically):
  Zu = Ref.comm223_232_tts_GIVc; report (Zu−Zf).TwoBodyNorm().
  For even π only, also note whether J TTS ≡ ethS (norm ≈ 0).

Usage:
  PYTHONPATH=build OMP_NUM_THREADS=1 python3 -B run/test_givc_JT_parity_grid.py [emax=1]
  N2B=10 TOL=1e-6 python3 -B run/test_givc_JT_parity_grid.py 1
"""
from __future__ import annotations

import os
import sys
import time
import traceback

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
tol = float(os.environ.get("TOL", "1e-6"))
n2b = int(os.environ.get("N2B", "10"))
jlist = [int(x) for x in os.environ.get("JLIST", "0,1,2,3,4").split(",") if x.strip()]
tlist = [int(x) for x in os.environ.get("TLIST", "0,1,2").split(",") if x.strip()]
plist = [int(x) for x in os.environ.get("PLIST", "0,1").split(",") if x.strip()]

eth = Commutator.FactorizedDoubleCommutator_eths
Ref = ReferenceImplementations
tts_GIVc = getattr(Ref, "comm223_232_tts_GIVc", None)

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


def collect_2b_cands(ms, nkeep):
    pts = []
    nch = ms.GetNumberTwoBodyChannels()
    for ch in range(nch):
        tbc = ms.GetTwoBodyChannel(ch)
        nk = tbc.GetNumberKets()
        J = tbc.J
        M2 = 2 * J
        for ib in range(nk):
            bra = tbc.GetKet(ib)
            i, j = int(bra.p), int(bra.q)
            jj = j2i(ms, j)
            for mi in mrange(ms, i):
                mj = M2 - mi
                if mj < -jj or mj > jj or ((mj + jj) % 2):
                    continue
                if i == j and mi == mj:
                    continue
                for ik in range(nk):
                    ket = tbc.GetKet(ik)
                    k, l = int(ket.p), int(ket.q)
                    mk, ml = mi, mj
                    if abs(mk) > j2i(ms, k) or abs(ml) > j2i(ms, l):
                        mk = j2i(ms, k)
                        ml = M2 - mk
                        if abs(ml) > j2i(ms, l):
                            continue
                    if k == l and mk == ml:
                        continue
                    pts.append((i, mi, j, mj, k, mk, l, ml))
                    break
                break
            if len(pts) > 8 * nkeep:
                break
        if len(pts) > 8 * nkeep:
            break
    return pts


def close(a, b):
    d = abs(a - b)
    scale = max(abs(a), abs(b), 1.0)
    return d <= tol * scale or d <= 1e-8


def run_case(ms, ut, lam, T, parity, seed, cands):
    ut.SetRandomSeed(seed)
    # RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)
    Eta = ut.RandomOp(ms, lam, T, parity, 2, -1)
    Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)
    if Gamma.IsReduced():
        Gamma.MakeNotReduced()
    # λ=0 RandomOp quirk: may come unreduced; ethS wants reduced Ω
    if not Eta.IsReduced():
        Eta.MakeReduced()

    Zf = make_Z(ms)
    t0 = time.time()
    reset()
    eth.comm223_232_GIVc(Eta, Gamma, Zf)
    dt_e = time.time() - t0

    Zu = None
    dt_u = 0.0
    jdiff = float("nan")
    tts_eq = None
    if tts_GIVc is not None:
        Zu = make_Z(ms)
        t0 = time.time()
        try:
            tts_GIVc(Eta, Gamma, Zu)
            dt_u = time.time() - t0
            Diff = Zu - Zf
            jdiff = Diff.TwoBodyNorm()
            if parity == 0:
                tts_eq = jdiff <= max(1e-8, tol * max(Zf.TwoBodyNorm(), 1.0))
        except Exception as ex:
            dt_u = time.time() - t0
            jdiff = float("nan")
            tts_eq = None
            print(
                f"    TTS diag skip ({type(ex).__name__}: {ex})",
                flush=True,
            )
            Zu = None

    # Sample largest |Zf| m-scheme MEs vs Mscheme_fact_GIVc
    ranked = []
    for tup in cands:
        zj = ut.GetMschemeMatrixElement_2b(Zf, *tup)
        ranked.append((abs(zj), tup, zj))
    ranked.sort(key=lambda t: t[0], reverse=True)
    pts = ranked[:n2b]

    n, nfail, maxd = 0, 0, 0.0
    for _a, tup, zj in pts:
        zm = ut.Mscheme_fact_GIVc(Eta, Gamma, *tup)
        n += 1
        d = abs(zj - zm)
        if d > maxd:
            maxd = d
        if not close(zj, zm):
            nfail += 1

    ok = nfail == 0 and n > 0
    return {
        "ok": ok,
        "n": n,
        "nfail": nfail,
        "maxd": maxd,
        "jdiff": jdiff,
        "tts_eq": tts_eq,
        "dt_e": dt_e,
        "dt_u": dt_u,
        "nrm_f": Zf.TwoBodyNorm(),
        "nrm_u": Zu.TwoBodyNorm() if Zu is not None else float("nan"),
        "eta_red": Eta.IsReduced(),
        "gam_red": Gamma.IsReduced(),
    }


def main():
    ms = ModelSpace(emax, "He4", "He4")
    ms.SetHbarOmega(20.0)
    ms.PreCalculateSixJ()
    ms.PreCalculateNineJ()
    ut = UnitTest(ms)
    cands = collect_2b_cands(ms, n2b)

    print(
        f"GIVc JT-parity grid  emax={emax} He4  λ={jlist}  T={tlist}  π={plist}  "
        f"N2B={n2b}  tol={tol}  Ω herm=-1  Γ herm=+1",
        flush=True,
    )
    print(
        "  GATE: GetMscheme(ethS GIVc) ≡ Mscheme_fact_GIVc\n"
        "  DIAG: (tts_GIVc − ethS).TwoBodyNorm()  [may FAIL; not a gate]\n"
        "  even-π note: TTS≡ethS when DIAG norm ≈ 0",
        flush=True,
    )
    print(
        f"{'λ':>2} {'T':>2} {'π':>2} {'seed':>5}  {'gate':6s}  "
        f"{'m nOK/n':>10s}  {'maxΔ':>10s}  {'||Zf||':>10s}  "
        f"{'DIAG||Zu-Zf||':>14s}  {'even TTS≡?':>10s}  dt",
        flush=True,
    )

    n_pass = n_tot = n_err = 0
    failures = []
    t_all = time.time()

    for lam in jlist:
        for T in tlist:
            for parity in plist:
                n_tot += 1
                seed = 11 + 100 * lam + 10 * T + parity
                try:
                    r = run_case(ms, ut, lam, T, parity, seed, cands)
                    if r["ok"]:
                        n_pass += 1
                        gate = "PASS"
                    else:
                        gate = "FAIL"
                        failures.append((lam, T, parity, seed, r))

                    tts_note = "—"
                    if parity == 0:
                        if r["tts_eq"] is True:
                            tts_note = "yes"
                        elif r["tts_eq"] is False:
                            tts_note = "no"
                        else:
                            tts_note = "?"

                    jdiff_s = (
                        f"{r['jdiff']:.3e}"
                        if r["jdiff"] == r["jdiff"]
                        else "nan"
                    )
                    print(
                        f"{lam:2d} {T:2d} {parity:2d} {seed:5d}  {gate:6s}  "
                        f"{r['n']-r['nfail']:3d}/{r['n']:<3d}     "
                        f"{r['maxd']:10.3e}  {r['nrm_f']:10.4g}  "
                        f"{jdiff_s:>14s}  {tts_note:>10s}  "
                        f"ethS {r['dt_e']:.2f}s"
                        + (f" TTS {r['dt_u']:.2f}s" if r["dt_u"] else ""),
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
        f"SUMMARY  m-gold gate: {n_pass}/{n_tot} PASS  errors={n_err}  ({dt:.1f}s)",
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
                    f"m fail={r['nfail']}/{r['n']} maxΔ={r['maxd']:.3e}  "
                    f"DIAG||Zu-Zf||={r['jdiff']:.3e}",
                    flush=True,
                )
    else:
        print("FAILURES: (none)", flush=True)

    return 0 if n_pass == n_tot and n_err == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
