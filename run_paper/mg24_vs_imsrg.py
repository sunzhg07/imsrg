#!/usr/bin/env python3
"""Mg24 VS-IMSRG: reference Mg24, valence sd-shell, inert O16 core.

Same interaction and solver settings as he8_mr_eom.py (EM1.8/2.0, hw=16,
NO2B 3N, Magnus + atan then shell-model-atan).

Outputs (hs_cache/):
  Hs_Mg24_sdshell_{LECs}_hw16_emax{E}.dat   Hs, NO wrt inert O16 core
  Hs_Mg24_sdshell_{LECs}_hw16_emax{E}.snt   sd-shell interaction for KSHELL

From run_paper/:
  python3 mg24_vs_imsrg.py              # EMAX=4
  EMAX=6 python3 mg24_vs_imsrg.py
"""
from __future__ import annotations

import os
import sys
import time

from he8_io import BUILD, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import (  # noqa: E402
    ElectricMultipoleOp,
    HartreeFock,
    IMSRGSolver,
    ModelSpace,
    Operator,
    OperatorFromString,
    ReadWrite,
)

emax = int(os.environ.get("EMAX", "4"))
hw = 16
nucleus = "Mg24"
ref = nucleus          # IMSRG reference (ensemble NO for open-shell Mg24)
val = "sd-shell"       # valence space decoupled by shell-model-atan
inert_core = "O16"     # core for the packaged Hs / KSHELL interaction

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"

Hs_file = hs_path(nucleus, hw, emax, val=val)
snt_file = Hs_file[: -len(".dat")] + ".snt"
Z_nuc, N_nuc = 12, 12
A_nuc = Z_nuc + N_nuc


def op_file(tag):
    """hs_cache/<tag>_Mg24_sdshell_..._emax{E}.dat next to the Hs file."""
    d, f = os.path.split(Hs_file)
    return os.path.join(d, tag + f[len("Hs"):])


def e1_operators(ms):
    """Bare E1 (rank 1, odd parity, reduced, e fm).

    E1p : protons only (imsrg_util::ElectricMultipoleOp(ms, 1)).
    E1c : recoil-corrected, e_p = N/A, e_n = -Z/A, which removes the
          spurious c.m. dipole to first order. The neutron part copies the
          proton one-body MEs onto the mirror neutron orbits (same radial
          and angular factors).
    """
    E1p = ElectricMultipoleOp(ms, 1)
    E1n = 0.0 * E1p
    for i in range(ms.GetNumberOrbits()):
        oi = ms.GetOrbit(i)
        if oi.tz2 != 1:
            continue
        ip = ms.GetOrbitIndex(oi.n, oi.l, oi.j2, -1)
        for j in range(ms.GetNumberOrbits()):
            oj = ms.GetOrbit(j)
            if oj.tz2 != 1:
                continue
            jp = ms.GetOrbitIndex(oj.n, oj.l, oj.j2, -1)
            v = E1p.GetOneBody(ip, jp)
            if v != 0.0:
                E1n.SetOneBodyME(i, j, v)
    E1c = (N_nuc / A_nuc) * E1p - (Z_nuc / A_nuc) * E1n
    return {"E1p": E1p, "E1c": E1c}


def main():
    t0 = time.time()
    print(f"Mg24 VS-IMSRG  emax={emax}  hw={hw}  ref={ref}  val={val}  "
          f"core={inert_core}", flush=True)
    rw = ReadWrite()
    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)

    H = Operator(ms, 0, 0, 0, 3)
    rw.ReadBareTBME_Darmstadt(f2b, H, 14, 28, 14)
    H.ThreeBody.SetMode("no2b")
    H.ThreeBody.ReadFile([f3b], [18, 36, 18])
    H += OperatorFromString(ms, "Trel")
    print(f"read interaction  ({time.time() - t0:.0f} s)", flush=True)

    hf = HartreeFock(H)
    hf.Solve()
    hf.PrintSPEandWF()
    HNO = hf.GetNormalOrderedH(2)
    print(f"E_HF (ensemble ref {ref}) = {HNO.ZeroBody:.6f} MeV", flush=True)

    solver = IMSRGSolver(HNO)
    solver.SetMethod("magnus")
    solver.SetHunterGatherer(True)
    solver.SetOmegaNormMax(0.1)
    solver.SetGenerator("atan")
    solver.SetSmax(50)
    solver.Solve()
    print(f"after core decoupling: E0 = {solver.GetH_s().ZeroBody:.6f} MeV  "
          f"({time.time() - t0:.0f} s)", flush=True)
    solver.SetGenerator("shell-model-atan")
    solver.SetSmax(100)
    solver.Solve()
    Hs = solver.GetH_s()
    print(f"VS Hs (NO wrt {ref}): E0 = {Hs.ZeroBody:.6f} MeV  "
          f"({time.time() - t0:.0f} s)", flush=True)

    # Repackage NO wrt the inert O16 core (as he8: UndoNO -> DoNormalOrderingCore).
    Hcore = Hs.UndoNormalOrdering()
    Hcore = Hcore.DoNormalOrderingCore()
    ms_core = ModelSpace(emax, inert_core, val)
    ms_core.SetHbarOmega(hw)
    Hcore.SetModelSpace(ms_core)
    print(f"core-NO ({inert_core} inert): E_core = {Hcore.ZeroBody:.6f} MeV",
          flush=True)

    rw.WriteOperator(Hcore, Hs_file)
    rw.WriteTokyo(Hcore, snt_file, "")
    print(f"wrote {Hs_file}\nwrote {snt_file}", flush=True)

    # E1 evolved with the same Omega, then repackaged NO wrt the O16 core.
    for tag, op in e1_operators(ms).items():
        bare_norm = op.OneBodyNorm()
        op_s = solver.Transform(op)
        op_c = op_s.UndoNormalOrdering().DoNormalOrderingCore()
        op_c.SetModelSpace(ms_core)
        path = op_file(tag)
        rw.WriteOperator(op_c, path)
        back = 0.0 * op_c
        rw.ReadOperator(back, path)
        diff = (back - op_c).Norm()
        print(f"{tag}: rank_J={op_c.GetJRank()} parity={op_c.GetParity()}  "
              f"||1b|| bare={bare_norm:.5f} evolved={op_c.OneBodyNorm():.5f}  "
              f"||2b||={op_c.TwoBodyNorm():.5f}  wrote {path} "
              f"(read-back |diff|={diff:.1e})", flush=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
