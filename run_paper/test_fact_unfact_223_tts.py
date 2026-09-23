#!/usr/bin/env python3
"""J-scheme leftover: nested 3b vs ethS factorized.

Unfactorized (builds W3):
  W = −comm223st(H, Ω)
  Z = comm231tts(Ω,W) + comm232tts(Ω,W) + comm132tts(Ω,W)

Factorized (no 3b store):
  Z = comm223_231_st + comm223_232 + comm223_132_tts

One script, grid over λ, T, P (env: LAMBDAS, TS, PARITIES).

He8: MS = He4 core + p-shell. Default emax=2, λ=0,4 × T=0,1,2 × even/odd.

From run_paper/:
  PYTHONPATH=../build python3 test_fact_unfact_223_tts.py
  LAMBDAS=4 TS=1 PARITIES=0,1 python3 test_fact_unfact_223_tts.py
"""
from __future__ import annotations

import os
import sys
import time
from contextlib import contextmanager

from pyIMSRG import *

EMAX = int(os.environ.get("EMAX", "2"))
LAMBDAS = [int(x) for x in os.environ.get("LAMBDAS", "0,1,2,3,4").split(",")]
TS = [int(x) for x in os.environ.get("TS", "0,1,2").split(",")]
PARITIES = [int(x) for x in os.environ.get("PARITIES", "0,1").split(",")]
TOL = float(os.environ.get("FACT_TOL", "1e-6"))
SEED = int(os.environ.get("SEED", "17"))
HUSH = os.environ.get("HUSH", "1") != "0"

eth = Commutator.FactorizedDoubleCommutator_eths
cm = Commutator


@contextmanager
def hush():
    if not HUSH:
        yield
        return
    devnull = os.open(os.devnull, os.O_WRONLY)
    old_out, old_err = os.dup(1), os.dup(2)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(old_out, 1)
        os.dup2(old_err, 2)
        os.close(devnull)
        os.close(old_out)
        os.close(old_err)


def make_Z(ms):
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    Z *= 0.0
    return Z


def make_W3(ms, lam, tz, parity):
    W = Operator(ms, lam, tz, parity, 3)
    W.SetHermitian()
    W.ThreeBody.SetMode("pn")
    W *= 0.0
    return W


def copy_omega_rank3(ms, omega, lam, tz, parity):
    eta3 = Operator(ms, lam, tz, parity, 3)
    if omega.IsHermitian():
        eta3.SetHermitian()
    else:
        eta3.SetAntiHermitian()
    eta3.ThreeBody.SetMode("pn")
    eta3 *= 0.0
    eta3.OneBody = omega.OneBody
    eta3.TwoBody = omega.TwoBody
    if lam == 0 and not eta3.IsReduced():
        eta3.MakeReduced()
    return eta3


def leftover_eths(ms, omega, H):
    Z = make_Z(ms)
    eth.SetUse_1b_Intermediates(True)
    eth.SetUse_2b_Intermediates(True)
    eth.comm223_231_st(omega, H, Z)
    Z232 = make_Z(ms)
    eth.comm223_232(omega, H, Z232)
    Z132 = make_Z(ms)
    eth.comm223_132_tts(omega, H, Z132)
    Z += Z232
    Z += Z132
    return Z


def leftover_nested(ms, omega3, H, lam, tz, parity):
    W = make_W3(ms, lam, tz, parity)
    cm.comm223st(H, omega3, W)
    W *= -1.0
    Z231 = make_Z(ms)
    cm.comm231tts(omega3, W, Z231)
    Z232 = make_Z(ms)
    cm.comm232tts(omega3, W, Z232)
    Z132 = make_Z(ms)
    cm.comm132tts(omega3, W, Z132)
    Z = make_Z(ms)
    Z += Z231
    Z += Z232
    Z += Z132
    return Z, W, {"231": Z231, "232": Z232, "132": Z132}


def omega_for_fact(omega, omega3, lam):
    return omega3 if lam == 0 else omega


def run_cell(ms, H, lam, tz, parity):
    unt = UnitTest(ms)
    unt.SetRandomSeed(SEED + 17 * lam + 5 * tz + 3 * parity)
    omega = unt.RandomOp(ms, lam, tz, parity, 2, -1)
    omega3 = copy_omega_rank3(ms, omega, lam, tz, parity)
    of = omega_for_fact(omega, omega3, lam)
    t0 = time.perf_counter()
    Z_un, W, pieces_un = leftover_nested(ms, omega3, H, lam, tz, parity)
    t_un = time.perf_counter() - t0
    t0 = time.perf_counter()
    Z_fa = leftover_eths(ms, of, H)
    t_fa = time.perf_counter() - t0
    dZ = (Z_un - Z_fa).Norm()
    n_un = Z_un.Norm()
    n_fa = Z_fa.Norm()
    scale = max(1.0, n_un, n_fa)
    rel = dZ / scale
    ok = dZ <= TOL * scale
    split = {}
    for name in ("231", "232", "132"):
        Zu = pieces_un[name]
        Zf = make_Z(ms)
        if name == "231":
            eth.comm223_231_st(of, H, Zf)
        elif name == "232":
            eth.comm223_232(of, H, Zf)
        else:
            eth.comm223_132_tts(of, H, Zf)
        dn = (Zu - Zf).Norm()
        sn = max(1.0, Zu.Norm(), Zf.Norm())
        split[name] = (Zu.Norm(), Zf.Norm(), dn / sn, dn <= TOL * sn)
    return dict(
        ok=ok,
        n_un=n_un,
        n_fa=n_fa,
        dZ=dZ,
        rel=rel,
        w3=W.ThreeBodyNorm(),
        t_un=t_un,
        t_fa=t_fa,
        split=split,
    )


def main():
    print(
        f"He8  MS=He4+p-shell  emax={EMAX}  "
        f"λ={LAMBDAS}  T={TS}  P={PARITIES}  tol={TOL}  seed={SEED}"
    )
    print("unfact:  W=−comm223st(H,Ω) ; Z=231tts+232tts+132tts")
    print("fact:    comm223_231_st + comm223_232 + comm223_132_tts")
    print("H, Ω random (H scalar Hermitian 2b; Ω antiherm 2b, same λTP as W3)")

    ms = ModelSpace(EMAX, "He4", "p-shell")
    ms.SetHbarOmega(16)
    ms.PreCalculateSixJ()
    ms.PreCalculateNineJ()
    unt = UnitTest(ms)
    unt.SetRandomSeed(SEED)
    H = unt.RandomOp(ms, 0, 0, 0, 2, 1)
    print(f"  ||H1||={H.OneBodyNorm():.4g}  ||H2||={H.TwoBodyNorm():.4g}")

    print(
        f"\n{'λ':>3} {'T':>3} {'P':>3}  {'||Z_un||':>12} {'||Z_fa||':>12} "
        f"{'||ΔZ||':>10} {'rel':>10}  "
        f"{'231':>5} {'232':>5} {'132':>5}  "
        f"{'t_un':>6} {'t_fa':>6}  status"
    )
    rows = []
    for lam in LAMBDAS:
        for tz in TS:
            for parity in PARITIES:
                try:
                    with hush():
                        rec = run_cell(ms, H, lam, tz, parity)
                except Exception as exc:
                    print(
                        f"{lam:3d} {tz:3d} {parity:3d}  "
                        f"{'':>12} {'':>12} {'':>10} {'':>10}  "
                        f"{'':>5} {'':>5} {'':>5}  "
                        f"{'':>6} {'':>6}  ERR  {exc}",
                        flush=True,
                    )
                    rows.append((lam, tz, parity, dict(ok=False, status="ERR")))
                    continue
                rec["status"] = "PASS" if rec["ok"] else "FAIL"
                rows.append((lam, tz, parity, rec))
                sp = rec["split"]
                print(
                    f"{lam:3d} {tz:3d} {parity:3d}  "
                    f"{rec['n_un']:12.5e} {rec['n_fa']:12.5e} "
                    f"{rec['dZ']:10.3e} {rec['rel']:10.3e}  "
                    f"{'OK' if sp['231'][3] else 'FAIL':>5} "
                    f"{'OK' if sp['232'][3] else 'FAIL':>5} "
                    f"{'OK' if sp['132'][3] else 'FAIL':>5}  "
                    f"{rec['t_un']:6.2f} {rec['t_fa']:6.2f}  {rec['status']}"
                    f"  ||W3||={rec['w3']:.3g}",
                    flush=True,
                )
                if not rec["ok"]:
                    for name in ("231", "232", "132"):
                        nu, nf, rel, okp = rec["split"][name]
                        print(
                            f"      {name}: ||Z_un||={nu:.5e}  ||Z_fa||={nf:.5e}  "
                            f"rel={rel:.3e}  {'OK' if okp else 'FAIL'}",
                            flush=True,
                        )

    n_pass = sum(1 for *_, r in rows if r.get("status") == "PASS")
    n_fail = sum(1 for *_, r in rows if r.get("status") == "FAIL")
    n_err = sum(1 for *_, r in rows if r.get("status") == "ERR")
    print(f"\n{n_pass} PASS  {n_fail} FAIL  {n_err} ERR  / {len(rows)} cells")
    ok = n_fail == 0 and n_err == 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
