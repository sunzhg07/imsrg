#!/usr/bin/env python3
"""He8 closed-shell SR-IMSRG + SR-EOM for 0+ and 2+.

Hs cache (parameters → filename):
  hs_cache/Hs_He8_EM1.8_2.0_hw16_emax3.dat

FCI .snt (JOB=fci; UndoNO after decouple, no EOM):
  hs_cache/Hs_He8_SR_FCI_EM1.8_2.0_hw16_emax3.snt

EOM JSON:
  eom_result/He8_EM1.8_2.0_hw16_emax3_sr0.json

From run_paper/:
  PYTHONPATH=../build python3 he8_sr_eom.py              # 0+
  PYTHONPATH=../build JRANK=2 python3 he8_sr_eom.py      # 2+
  PHASE=1 python3 he8_sr_eom.py                          # rebuild Hs
  python3 he8_sr_eom.py JOB=fci                          # decouple + UndoNO snt
"""
from __future__ import annotations

import json
import os
import sys

from he8_io import BUILD, eom_result_path, fci_snt_path, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import (  # noqa: E402
    EOM,
    HartreeFock,
    IMSRGSolver,
    ModelSpace,
    Operator,
    OperatorFromString,
    ReadWrite,
)

emax = int(os.environ.get("EMAX", "3"))
hw = 16
nucleus = "He8"
ref = nucleus
val = nucleus  # closed shell
Hs_file = hs_path(nucleus, hw, emax, val=val)

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"


def make_ms():
    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)
    return ms


def write_fci_hamiltonian(Hs):
    """Undo He8 normal ordering of SR Hs and write a vacuum FCI .snt."""
    Hvac = Hs.UndoNormalOrdering()
    msfci = ModelSpace(emax, "vacuum", "FCI")
    msfci.SetHbarOmega(hw)
    Hvac.SetModelSpace(msfci)
    path = fci_snt_path(nucleus, hw, emax, tag="SR")
    ReadWrite().WriteTokyo(Hvac, path, "")
    print(f"FCI UndoNO ZeroBody={Hvac.ZeroBody:.8f} MeV", flush=True)
    print(f"wrote {path}", flush=True)
    return path


def decouple_hs():
    """Closed-shell SR-IMSRG. Returns (ms, Hs, solver, HNO); Hs NO wrt He8."""
    rw = ReadWrite()
    ms = make_ms()
    H = Operator(ms, 0, 0, 0, 3)
    rw.ReadBareTBME_Darmstadt(f2b, H, 14, 28, 14)
    H.ThreeBody.SetMode("no2b")
    H.ThreeBody.ReadFile([f3b], [18, 36, 18])
    H += OperatorFromString(ms, "Trel")
    hf = HartreeFock(H)
    hf.Solve()
    hf.PrintSPEandWF()
    # HNO must outlive the solver: IMSRGSolver stores Operator* H_0.
    HNO = hf.GetNormalOrderedH(2)
    solver = IMSRGSolver(HNO)
    solver.SetMethod("magnus")
    solver.SetHunterGatherer(True)
    solver.SetOmegaNormMax(0.1)
    solver.SetGenerator("atan")
    solver.SetSmax(50)
    solver.Solve()
    Hs = solver.GetH_s()
    print(f"SR Hs (NO wrt {ref}): ZeroBody={Hs.ZeroBody:.8f} MeV", flush=True)
    return ms, Hs, solver, HNO


def build_hs():
    rw = ReadWrite()
    phase = int(os.environ.get("PHASE", "2" if os.path.isfile(Hs_file) else "1"))
    print(f"PHASE={phase}  Hs={Hs_file}  exists={os.path.isfile(Hs_file)}",
          flush=True)

    if phase != 1:
        if not os.path.isfile(Hs_file):
            raise FileNotFoundError(Hs_file)
        ms = make_ms()
        Hs = Operator(ms, 0, 0, 0, 3)
        rw.ReadOperator(Hs, Hs_file)
        print(f"read {Hs_file}  ZeroBody={Hs.ZeroBody:.8f} MeV", flush=True)
        return ms, Hs

    ms, Hs, solver, HNO = decouple_hs()
    rw.WriteOperator(Hs, Hs_file)
    print(f"wrote {Hs_file}", flush=True)
    del solver, HNO
    Hs = Operator(ms, 0, 0, 0, 3)
    rw.ReadOperator(Hs, Hs_file)
    return ms, Hs


def dump_result(name, lab, J2, res):
    energies = [float(x) for x in res.arnoldi.energies]
    eref = float(res.eref)
    rec = {
        "name": name,
        "lab": lab,
        "kind": name,
        "mode": "sr",
        "J2": J2,
        "ref_file": None,
        "nucleus": nucleus,
        "lecs": os.environ.get("LECS", "EM1.8_2.0"),
        "hw": hw,
        "emax": emax,
        "hs_file": Hs_file,
        "eref": eref,
        "E_ex": energies,
        "E_abs": [e + eref for e in energies],
        "converged": bool(getattr(res.arnoldi, "converged", False)),
        "stop_reason": str(getattr(res.arnoldi, "stop_reason", "")),
        "steps": int(getattr(res.arnoldi, "steps", 0)),
        "residuals": [float(x) for x in getattr(res.arnoldi, "residuals", [])],
    }
    path = eom_result_path(nucleus, hw, emax, name)
    with open(path, "w") as f:
        json.dump(rec, f, indent=2)
    print(f"wrote {path}", flush=True)
    return rec


def requested_job():
    job = os.environ.get("JOB", "")
    if not job and len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        job = sys.argv[1]
    if job.startswith("JOB="):
        job = job.split("=", 1)[1]
    return job


def main():
    if requested_job() == "fci":
        print("JOB=fci  (decouple + write snt, no EOM)", flush=True)
        ms, Hs, solver, HNO = decouple_hs()
        write_fci_hamiltonian(Hs)
        del solver, HNO, ms
        return

    J2 = int(os.environ.get("JRANK", "0"))
    parity = int(os.environ.get("PARITY", "0"))
    itz = int(os.environ.get("ITZ", "0"))
    max_iter = int(os.environ.get("MAX_ITER", "80"))
    state_want = int(os.environ.get("STATE_WANT", "5"))
    name = "sr0" if J2 == 0 else f"sr{J2}"
    lab = rf"SR ${J2}^+$"

    ms, Hs = build_hs()
    print(f"ZeroBody (E_ref) = {Hs.ZeroBody:.6f} MeV", flush=True)
    print(f"SR-EOM  JπTz = {J2} {parity} {itz}  max_iter={max_iter}  "
          f"state_want={state_want}", flush=True)

    eom = EOM(Hs, J2, parity, itz)
    res = eom.Run(max_iter, state_want)

    print(f"\n=== He8 closed-shell SR-EOM {J2}+ ===", flush=True)
    print(f"eref (ZeroBody) = {res.eref:.6f} MeV", flush=True)
    print(f"stop={res.arnoldi.stop_reason}  steps={res.arnoldi.steps}  "
          f"converged={res.arnoldi.converged}", flush=True)
    for k, e in enumerate(res.arnoldi.energies):
        print(f"  E({k}): excitation={e:.6f}  absolute={e + res.eref:.6f} MeV",
              flush=True)
    dump_result(name, lab, J2, res)


if __name__ == "__main__":
    main()
