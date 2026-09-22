#!/usr/bin/env python3
"""Tensor ComputeNorm benchmark: v^T N v vs leftover ½⟨[Qa, Qb⁻]⟩_ρ.

One script, grid over λ, T, P (env: LAMBDAS, TS, PARITIES).

He8: MS = He4 core + p-shell. Default emax=2, λ=1,2,3,4 × T=0,1,2 × even/odd.
Gold is leftover TTS (NormMultiref), not production Commutator() on T×T.

From run_paper/:
  PYTHONPATH=../build python3 test_computenorm_tensor.py
  LAMBDAS=3 TS=1 PARITIES=0,1 python3 test_computenorm_tensor.py
"""
from __future__ import annotations

import os
import sys
from contextlib import contextmanager

import numpy as np
from pyIMSRG import *

from he8_io import RUN2

EMAX = int(os.environ.get("EMAX", "2"))
N_PAIR = int(os.environ.get("N_PAIR", "2"))
TOL = float(os.environ.get("NORM_TOL", "1e-6"))
LAMBDAS = [int(x) for x in os.environ.get("LAMBDAS", "1,2,3,4").split(",")]
TS = [int(x) for x in os.environ.get("TS", "0,1,2").split(",")]
PARITIES = [int(x) for x in os.environ.get("PARITIES", "0,1").split(",")]
HS_FILE = os.environ.get("HS_FILE", os.path.join(RUN2, "Hs_back"))
REF_FILE = os.environ.get("REF_FILE", os.path.join(RUN2, "he8.ref"))


@contextmanager
def hush():
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


def build_ms():
    ms = ModelSpace(EMAX, "He4", "p-shell")
    ms.SetHbarOmega(16)
    Hs = Operator(ms, 0, 0, 0, 3)
    if os.path.isfile(HS_FILE):
        try:
            ReadWrite().ReadOperator(Hs, HS_FILE)
        except Exception:
            Hs *= 0.0
    return ms, Hs


def run_one(ms, Hs, lam, tz, parity):
    tensor = not (lam == 0 and tz == 0 and parity == 0)
    eom = EOM(Hs, REF_FILE, lam, parity, tz)
    eom.SetUseRdm3(False)
    eom.SetIncludeConfigs(True, True, True, True, True)
    if tensor:
        eom.ConstructConfigs_tensor()
        eom.ConstructNormMatrix_tensor()
    else:
        eom.ConstructConfigs()
        eom.ConstructNormMatrix()
    dim = int(eom.eom_dims)
    if dim == 0:
        return dict(status="SKIP", dim=0, reason="empty")
    N = np.array(eom.Nkernel, dtype=float)
    chis, vecs = [], []
    for k in range(N_PAIR):
        unt = UnitTest(ms)
        unt.SetRandomSeed(1000 + 17 * lam + 5 * tz + 3 * parity + k)
        h = unt.RandomOp(ms, lam, tz, parity, 2, 1)
        chi = eom.GetVSEOM_ladder_multiref(h, 1)
        chis.append(chi)
        vecs.append(np.array(eom.FlattenOperator(chi), dtype=float))
    worst = 0.0
    n_ok = 0
    n_tot = 0
    n_zero = 0
    for i in range(N_PAIR):
        for j in range(N_PAIR):
            matmul = float(vecs[i] @ (N @ vecs[j]))
            cn = float(eom.ComputeNorm(chis[i], chis[j]))
            gold = float(eom.NormMultiref(chis[i], chis[j]))
            d_cn = abs(matmul - cn)
            d = abs(cn - gold)
            scale = max(1.0, abs(cn), abs(gold))
            rel = d / scale
            n_tot += 1
            if abs(cn) < 1e-14 and abs(gold) < 1e-14:
                n_zero += 1
            ok = (d <= TOL or rel <= TOL) and d_cn < 1e-8
            n_ok += int(ok)
            worst = max(worst, rel)
    status = "PASS" if n_ok == n_tot else "FAIL"
    return dict(
        status=status,
        dim=dim,
        tensor=tensor,
        qv=eom.qv_dim,
        ph=eom.ph_dim,
        ppvv=eom.ppvv_dim,
        pphv=eom.pphv_dim,
        pphh=eom.pphh_dim,
        n_ok=n_ok,
        n_tot=n_tot,
        n_zero=n_zero,
        worst=worst,
    )


def main():
    print(
        f"He8 nmax={EMAX}  MS=He4+p-shell  RDM={REF_FILE}  "
        f"n_pair={N_PAIR}  tol={TOL}"
    )
    print(f"grid  λ={LAMBDAS}  T={TS}  P={PARITIES}")
    print("ComputeNorm = v^T N v   gold = 1/2 ⟨[Qa, Qb^-]⟩_ρ  (TTS leftovers)")
    if not os.path.isfile(REF_FILE):
        print(f"missing RDM: {REF_FILE}", file=sys.stderr)
        return 2
    ms, Hs = build_ms()
    rows = []
    print(
        f"\n{'λ':>3} {'T':>3} {'P':>3}  {'kern':>6} {'dim':>6}  "
        f"{'qv':>4} {'ph':>4} {'ppvv':>5} {'pphv':>5} {'pphh':>5}  "
        f"{'pairs':>7}  {'worst rel':>10}  status"
    )
    for lam in LAMBDAS:
        for tz in TS:
            for parity in PARITIES:
                try:
                    with hush():
                        rec = run_one(ms, Hs, lam, tz, parity)
                except Exception as exc:
                    rec = dict(status="ERR", dim=0, reason=str(exc)[:80])
                    print(
                        f"{lam:3d} {tz:3d} {parity:3d}  {'':>6} {0:6d}  "
                        f"{'':>4} {'':>4} {'':>5} {'':>5} {'':>5}  "
                        f"{'':>7}  {'':>10}  ERR  {rec['reason']}",
                        flush=True,
                    )
                    rows.append((lam, tz, parity, rec))
                    continue
                rec["lam"] = lam
                rec["tz"] = tz
                rec["parity"] = parity
                rows.append((lam, tz, parity, rec))
                if rec["status"] == "SKIP":
                    print(
                        f"{lam:3d} {tz:3d} {parity:3d}  {'':>6} {0:6d}  "
                        f"{'':>4} {'':>4} {'':>5} {'':>5} {'':>5}  "
                        f"{'':>7}  {'':>10}  SKIP empty",
                        flush=True,
                    )
                    continue
                kern = "tensor" if rec["tensor"] else "scalar"
                print(
                    f"{lam:3d} {tz:3d} {parity:3d}  {kern:>6} {rec['dim']:6d}  "
                    f"{rec['qv']:4d} {rec['ph']:4d} {rec['ppvv']:5d} "
                    f"{rec['pphv']:5d} {rec['pphh']:5d}  "
                    f"{rec['n_ok']:2d}/{rec['n_tot']:<4d}  "
                    f"{rec['worst']:10.3e}  {rec['status']}",
                    flush=True,
                )

    n_pass = sum(1 for *_, r in rows if r["status"] == "PASS")
    n_fail = sum(1 for *_, r in rows if r["status"] == "FAIL")
    n_skip = sum(1 for *_, r in rows if r["status"] == "SKIP")
    n_err = sum(1 for *_, r in rows if r["status"] == "ERR")
    print(f"\n{n_pass} PASS  {n_fail} FAIL  {n_skip} SKIP  {n_err} ERR  "
          f"/ {len(rows)} cells")
    ok = n_fail == 0 and n_err == 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
