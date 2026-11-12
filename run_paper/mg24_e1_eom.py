#!/usr/bin/env python3
"""Mg24: 0+ -> 1- MR-EOM and E1 strength of each 1- state.

Reference: Mg24 0+_1 density, run2/mg24/mg24_gs_e{E}.ref (OSM format,
written by mg24_kshell_rdm.py from the KSHELL transit densities).
Hs and the evolved E1 (recoil-corrected, e_p=N/A, e_n=-Z/A) from
mg24_vs_imsrg.py, both NO wrt the inert O16 core.

1- states: tensor MR-EOM (J=1, parity -) with the Arnoldi started from the
excitation part of E1 acting on the reference (Lanczos strength method), so
the states that carry E1 strength converge first.

Transition amplitude (scalar-coupled, ethS leftover-CG convention):
  T_nu = NormMultiref_tensor(E1, chi_nu) = <[E1, chi_nu - chi_nu^dag]>_rho / 2
  N_nu = ComputeNorm(chi_nu, chi_nu)
  B(E1; 0+ -> 1-_nu) = -lambda_hat^2 * T_nu^2 / N_nu     [e^2 fm^2]
  (code N, T carry -1/lambda_hat^2; TEST=1 pins this.)
TEST=1 checks this formula on a pure particle-hole chi against the exact
|<a||E1||i>|^2 with the bare E1.

From run_paper/:
  TEST=1 python3 mg24_e1_eom.py
  MAX_ITER=60 STATE_WANT=40 python3 mg24_e1_eom.py
"""
from __future__ import annotations

import json
import math
import os
import sys

from he8_io import BUILD, RUN2, eom_result_path, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import (  # noqa: E402
    EOM,
    ElectricMultipoleOp,
    ModelSpace,
    Operator,
    ReadWrite,
    UnitTest,
)

emax = int(os.environ.get("EMAX", "4"))
hw = 16
nucleus = "Mg24"
val = "sd-shell"
inert_core = "O16"
J2, parity, itz = 1, 1, 0          # 1- states
LAM = 1
HAT = math.sqrt(2 * LAM + 1)

Hs_file = hs_path(nucleus, hw, emax, val=val)
E1_TAG = os.environ.get("E1_OP", "E1c")
E1_file = os.path.join(os.path.dirname(Hs_file), E1_TAG + os.path.basename(Hs_file)[len("Hs"):])
ref_file = os.environ.get("REF_FILE", os.path.join(RUN2, "mg24", f"mg24_gs_e{emax}.ref"))


def setup():
    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    rw = ReadWrite()
    Hs = Operator(ms, 0, 0, 0, 3)
    rw.ReadOperator(Hs, Hs_file)
    E1 = Operator(ms, LAM, 0, 1, 2)
    rw.ReadOperator(E1, E1_file)
    eom = EOM(Hs, ref_file, J2, parity, itz)
    eom.SetUseRdm3(False)
    eom.PrepareHamiltonianForArnoldi()
    eom.ConstructConfigs_tensor()
    eom.ConstructNormMatrix_tensor()
    eom.ConstructProjectMatrix()
    eref = eom.GetVSEOM_Overlap_multiref(Hs)
    return ms, Hs, E1, eom, eref


def strength(eom, E1, chi):
    N = eom.ComputeNorm(chi, chi)
    T = eom.NormMultiref_tensor(E1, chi)
    # Code N, T = -(1/lambda_hat^2) x physical reduced quantities (checked with
    # TEST=1: unit core->q particle-hole gives N=-1/3), so B = -lambda_hat^2 T^2/N.
    return N, T, (-HAT * HAT * T * T / N if N != 0 else float("nan"))


def run_test(ms, eom):
    """Single particle-hole chi = (a^dag i)^1 with the bare proton E1.

    Prints raw N = ComputeNorm(chi,chi), T = NormMultiref_tensor(E1,chi),
    lambda_hat*T^2/N and |<a||E1||i>|^2 for E1-allowed pairs (Delta N = 1),
    to pin the tensor normalization convention.
    """
    E1bare = ElectricMultipoleOp(ms, 1)
    pairs = (("0p3/2 core -> 0d5/2 val", (0, 1, 3), (0, 2, 5)),
             ("0p1/2 core -> 1s1/2 val", (0, 1, 1), (1, 0, 1)),
             ("0d5/2 val  -> 0f7/2 q  ", (0, 2, 5), (0, 3, 7)),
             ("1s1/2 val  -> 1p3/2 q  ", (1, 0, 1), (1, 1, 3)))
    for lab, (ni, li, ji), (na, la, ja) in pairs:
        i = ms.GetOrbitIndex(ni, li, ji, -1)
        a = ms.GetOrbitIndex(na, la, ja, -1)
        x = Operator(ms, LAM, 0, 1, 2)
        x.SetOneBodyME(a, i, 1.0)
        chi = eom.GetVSEOM_ladder_multiref(x, 1)
        N = eom.ComputeNorm(chi, chi)
        T = eom.NormMultiref_tensor(E1bare, chi)
        e1 = E1bare.GetOneBody(a, i)
        B = HAT * T * T / N if N != 0 else float("nan")
        print(f"TEST {lab}: occ_i={ms.GetOrbit(i).occ:.4f} occ_a={ms.GetOrbit(a).occ:.4f}  "
              f"N={N:+.6e}  T={T:+.6e}  <a||E1||i>={e1:+.6e}  "
              f"lam_hat*T^2/N={B:+.6e}  |<a||E1||i>|^2={e1 * e1:.6e}  "
              f"T/<a||E1||i>={T / e1 if e1 else float('nan'):+.6f}", flush=True)
    return True


def main():
    ms, Hs, E1, eom, eref = setup()
    print(f"Mg24 1- MR-EOM  emax={emax}  ref={ref_file}\n  E_ref=<H>_rho={eref:.6f} MeV  "
          f"E1 op={E1_file}", flush=True)
    if os.environ.get("TEST", "0") != "0":
        sys.exit(0 if run_test(ms, eom) else 1)

    max_iter = int(os.environ.get("MAX_ITER", "60"))
    state_want = int(os.environ.get("STATE_WANT", "40"))
    e1seed = eom.GetVSEOM_ladder_multiref(E1, 1)
    eom.ProjectOprator(e1seed)
    N0, T0, B0 = strength(eom, E1, e1seed)
    print(f"E1 excitation part (projected):  N={N0:.6e}  -lam^2 T^2/N={B0:.6e}", flush=True)
    if os.environ.get("SEED", "random") == "e1":
        seed = e1seed
    else:
        # Same start as EOM::RunMR: random rank-J op -> ladder(+1) -> projector.
        ut = UnitTest(ms)
        ut.SetRandomSeed(int(os.environ.get("RSEED", "11")))
        seed = eom.GetVSEOM_ladder_multiref(ut.RandomOp(ms, J2, itz, parity, 2, 1), 1)
        eom.ProjectOprator(seed)
    n0 = eom.ComputeNorm(seed, seed)
    print(f"Arnoldi start ({os.environ.get('SEED', 'random')}): N={n0:.6e}", flush=True)
    if n0 <= 0:
        raise RuntimeError("start vector has non-positive N-norm")
    seed = seed * (1.0 / math.sqrt(n0))

    res = eom.ArnoldiSolveH2(seed, max_iter, state_want)
    energies = [float(e) for e in res.energies]
    rows = []
    print(f"\n{'k':>3} {'E_x (MeV)':>10} {'N':>10} {'T':>12} {'B(E1)up':>12} {'<E1|nu>':>10}",
          flush=True)
    for k, (e, chi) in enumerate(zip(energies, res.ritz)):
        N, T, B = strength(eom, E1, chi)
        ov = eom.ComputeNorm(e1seed, chi) / math.sqrt(abs(N)) if N != 0 else float("nan")
        rows.append({"k": k, "E_x": e, "E_abs": e + eref, "N": N, "T": T,
                     "B_E1_up": B, "seed_overlap": ov})
        print(f"{k:3d} {e:10.4f} {N:10.4e} {T:+12.5e} {B:12.5e} {ov:+10.4f}", flush=True)
    s0 = sum(r["B_E1_up"] for r in rows if r["B_E1_up"] == r["B_E1_up"])
    s1 = sum(r["B_E1_up"] * r["E_x"] for r in rows if r["B_E1_up"] == r["B_E1_up"])
    print(f"\nsum B(E1) = {s0:.5e} e^2fm^2   (seed total {B0:.5e})", flush=True)
    print(f"sum E_x B(E1) = {s1:.5e} MeV e^2fm^2", flush=True)
    path = eom_result_path(nucleus, hw, emax, "e1_1m")
    with open(path, "w") as f:
        json.dump({"emax": emax, "ref_file": ref_file, "hs_file": Hs_file,
                   "e1_file": E1_file, "E_ref": eref, "max_iter": max_iter,
                   "stop_reason": str(res.stop_reason), "converged": bool(res.converged),
                   "seed_B_total": B0, "sum_B": s0, "sum_EB": s1, "states": rows},
                  f, indent=2)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
