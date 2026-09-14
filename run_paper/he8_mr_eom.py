#!/usr/bin/env python3
"""He8 MR-EOM on the VS Hamiltonian (He4-core NO) + OSM RDM.

VS-IMSRG model space (decoupling):
  ModelSpace(emax, ref=He8, val=p-shell)
    ref  = core_ref for IMSRG: defines particle / hole occupations
           (He8 → 0s1/2 occupied; neutron 0p3/2 occupied / hole)
    val  = valence orbits (p-shell → 0p3/2, 0p1/2)

After core + valence decoupling, the same Hs (still NO wrt He8) is used for:
  1. write_pshell_hamiltonian — UndoNO + DoNormalOrderingCore(He4)
     → KSHELL p-shell effective interaction (.snt)
  2. write_fci_hamiltonian — UndoNO to true vacuum + FCI model space
     → full-space FCI interaction (.snt)

MR-EOM itself uses the He4-core NO Hs (.dat in hs_cache/).

Hs cache:
  hs_cache/Hs_He8_pshell_EM1.8_2.0_hw16_emax3.dat

Jobs (JOB=...):
  mr0gs    scalar  0+_gs → 0+     (he8.ref)
  mr0x     scalar  0+_1  → 0+     (he8/he8_exc0.ref)
  mr2gs    scalar  2+_gs → 2+     (he8/he8_2p_gs.ref)
  mr2x     scalar  2+_1  → 2+     (he8/he8_2p_exc.ref)
  t02      tensor  0+_gs → 2+     (he8.ref, J=2)
  t20      tensor  2+_gs → 0+     (he8/he8_2p_gs.ref, J=2)
  pshell   decouple → write_pshell_hamiltonian (no EOM)
  fci      decouple → write_fci_hamiltonian (no EOM)
  snt      decouple → both writers from the same Hs (no EOM)

From run_paper/:
  PYTHONPATH=../build python3 he8_mr_eom.py JOB=mr0gs
  PYTHONPATH=../build MAX_ITER=100 python3 he8_mr_eom.py JOB=t02
  PHASE=1 python3 he8_mr_eom.py JOB=mr0gs     # rebuild Hs .dat
  python3 he8_mr_eom.py JOB=pshell            # He4-core p-shell .snt
  python3 he8_mr_eom.py JOB=fci               # vacuum FCI .snt
  python3 he8_mr_eom.py JOB=snt               # both .snt from one Hs
"""
from __future__ import annotations

import json
import os
import sys

from he8_io import (
    BUILD,
    RUN2,
    eom_result_path,
    fci_snt_path,
    hs_path,
    lecs,
    pshell_snt_path,
)

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
# IMSRG ref: particle/hole occupations (He8 holes, not the inert He4 core).
ref = nucleus
# Valence orbits for VS decoupling / off-diagonal generator.
val = "p-shell"
# Inert core used only when packaging the p-shell interaction / MR-EOM Hs.
inert_core = "He4"
Hs_file = hs_path(nucleus, hw, emax, val=val)

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"

JOBS = {
    "mr0gs": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}0^+$",
        "kind": "mr0",
        "J2": 0,
        "ref": "he8.ref",
    },
    "mr0x": {
        "lab": r"MR $0^+_1{\to}0^+$",
        "kind": "mr0x",
        "J2": 0,
        "ref": os.path.join("he8", "he8_exc0.ref"),
    },
    "mr2gs": {
        "lab": r"MR $2^+_{\mathrm{gs}}{\to}2^+$",
        "kind": "mr2",
        "J2": 0,
        "ref": os.path.join("he8", "he8_2p_gs.ref"),
    },
    "mr2x": {
        "lab": r"MR $2^+_1{\to}2^+$",
        "kind": "mr2x",
        "J2": 0,
        "ref": os.path.join("he8", "he8_2p_exc.ref"),
    },
    "t02": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}2^+$",
        "kind": "t02",
        "J2": 2,
        "ref": "he8.ref",
    },
    "t20": {
        "lab": r"MR $2^+_{\mathrm{gs}}{\to}0^+$",
        "kind": "t20",
        "J2": 2,
        "ref": os.path.join("he8", "he8_2p_gs.ref"),
    },
}


def make_vs_ms():
    """Model space for VS-IMSRG: ref=He8 (holes), val=p-shell."""
    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)
    return ms


def make_core_ms():
    """Model space for He4-core NO Hs / MR-EOM / KSHELL p-shell."""
    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    return ms


def to_he4_core_hs(Hs):
    """Repackage VS Hs (NO wrt He8) as NO wrt inert He4 core."""
    Hs = Hs.UndoNormalOrdering()
    Hs = Hs.DoNormalOrderingCore()
    ms_core = make_core_ms()
    Hs.SetModelSpace(ms_core)
    print(f"core-NO ({inert_core} inert): ZeroBody = {Hs.ZeroBody:.8f} MeV",
          flush=True)
    return ms_core, Hs


def write_pshell_hamiltonian(Hs):
    """Write KSHELL p-shell effective interaction from VS-decoupled Hs.

    Input Hs is still normal-ordered wrt the IMSRG ref (He8). We UndoNO,
    then DoNormalOrderingCore so the interaction is NO wrt the inert He4
    core (holes = 0s1/2 only; valence = p-shell).
    """
    Hcore = Hs.UndoNormalOrdering()
    Hcore = Hcore.DoNormalOrderingCore()
    ms_core = make_core_ms()
    Hcore.SetModelSpace(ms_core)
    path = pshell_snt_path(nucleus, hw, emax)
    ReadWrite().WriteTokyo(Hcore, path, "")
    print(f"p-shell core-NO ({inert_core}) ZeroBody={Hcore.ZeroBody:.8f} MeV",
          flush=True)
    print(f"wrote {path}", flush=True)
    return path


def write_fci_hamiltonian(Hs):
    """Write vacuum FCI interaction from the same VS-decoupled Hs.

    UndoNormalOrdering removes the He8 reference so the Hamiltonian is
    normal-ordered wrt the true vacuum; ModelSpace(vacuum, FCI) includes
    every orbit for WriteTokyo.
    """
    Hvac = Hs.UndoNormalOrdering()
    msfci = ModelSpace(emax, "vacuum", "FCI")
    msfci.SetHbarOmega(hw)
    Hvac.SetModelSpace(msfci)
    path = fci_snt_path(nucleus, hw, emax)
    ReadWrite().WriteTokyo(Hvac, path, "")
    print(f"FCI vacuum ZeroBody={Hvac.ZeroBody:.8f} MeV", flush=True)
    print(f"wrote {path}", flush=True)
    return path


def decouple_hs():
    """Core + valence VS-IMSRG. Returns Hs still NO wrt ref=He8."""
    rw = ReadWrite()
    ms_vs = make_vs_ms()
    H = Operator(ms_vs, 0, 0, 0, 3)
    rw.ReadBareTBME_Darmstadt(f2b, H, 14, 28, 14)
    H.ThreeBody.SetMode("no2b")
    H.ThreeBody.ReadFile([f3b], [18, 36, 18])
    H += OperatorFromString(ms_vs, "Trel")
    hf = HartreeFock(H)
    hf.Solve()
    hf.PrintSPEandWF()
    HNO = hf.GetNormalOrderedH(2)
    solver = IMSRGSolver(HNO)
    solver.SetMethod("magnus")
    solver.SetHunterGatherer(True)
    solver.SetOmegaNormMax(0.1)
    solver.SetGenerator("atan")
    solver.SetSmax(50)
    solver.Solve()
    solver.SetGenerator("shell-model-atan")
    solver.SetSmax(100)
    solver.Solve()
    Hs = solver.GetH_s()
    print(f"VS Hs (NO wrt {ref}, val={val}): ZeroBody={Hs.ZeroBody:.8f} MeV",
          flush=True)
    return ms_vs, Hs, solver, HNO


def build_hs():
    rw = ReadWrite()
    phase = int(os.environ.get("PHASE", "2" if os.path.isfile(Hs_file) else "1"))
    print(f"PHASE={phase}  Hs={Hs_file}  exists={os.path.isfile(Hs_file)}",
          flush=True)

    if phase != 1:
        if not os.path.isfile(Hs_file):
            raise FileNotFoundError(Hs_file)
        ms_core = make_core_ms()
        Hs = Operator(ms_core, 0, 0, 0, 3)
        rw.ReadOperator(Hs, Hs_file)
        print(f"read {Hs_file}  ZeroBody={Hs.ZeroBody:.8f} MeV", flush=True)
        return ms_core, Hs

    ms_vs, Hs, solver, HNO = decouple_hs()
    ms_core, Hs = to_he4_core_hs(Hs)
    rw.WriteOperator(Hs, Hs_file)
    print(f"wrote {Hs_file}", flush=True)
    del solver, HNO, ms_vs
    Hs = Operator(ms_core, 0, 0, 0, 3)
    rw.ReadOperator(Hs, Hs_file)
    return ms_core, Hs


def dump_result(name, spec, eref, res):
    energies = [float(x) for x in res.arnoldi.energies]
    rec = {
        "name": name,
        "lab": spec["lab"],
        "kind": spec["kind"],
        "mode": "mr",
        "J2": spec["J2"],
        "ref_file": spec["ref"],
        "nucleus": nucleus,
        "val": val,
        "lecs": lecs(),
        "hw": hw,
        "emax": emax,
        "hs_file": Hs_file,
        "eref": float(eref),
        "E_ex": energies,
        "E_abs": [e + float(eref) for e in energies],
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


def resolve_job():
    job = requested_job()
    if job and job in JOBS:
        spec = dict(JOBS[job])
        spec["ref"] = os.path.join(RUN2, spec["ref"])
        return job, spec
    J2 = int(os.environ.get("JRANK", "0"))
    ref_file = os.environ.get("REF_FILE", os.path.join(RUN2, "he8.ref"))
    if not os.path.isabs(ref_file):
        ref_file = os.path.join(RUN2, ref_file)
    name = os.environ.get("NAME", f"mr_J{J2}")
    spec = {
        "lab": rf"MR $J={J2}$",
        "kind": name,
        "J2": J2,
        "ref": ref_file,
    }
    return name, spec


def run_snt_jobs(which: str):
    """Decouple once; write p-shell and/or FCI .snt from that Hs (no EOM)."""
    print(f"JOB={which}  (decouple + write interaction .snt, no EOM)", flush=True)
    ms_vs, Hs, solver, HNO = decouple_hs()
    if which in ("pshell", "snt"):
        write_pshell_hamiltonian(Hs)
    if which in ("fci", "snt"):
        write_fci_hamiltonian(Hs)
    del solver, HNO, ms_vs


def main():
    job = requested_job()
    if job in ("pshell", "fci", "snt"):
        run_snt_jobs(job)
        return

    name, spec = resolve_job()

    max_iter = int(os.environ.get("MAX_ITER", "50"))
    state_want = int(os.environ.get("STATE_WANT", "5"))
    use_h3 = os.environ.get("USE_H3", "1") != "0"
    use_projection = os.environ.get("USE_PROJECTION", "1") != "0"
    print_timing = os.environ.get("PRINT_TIMING", "0") != "0"
    J2 = spec["J2"]
    ref_file = spec["ref"]
    if not os.path.isfile(ref_file):
        raise FileNotFoundError(ref_file)

    ms, Hs = build_hs()
    print(f"JOB={name}  RDM={ref_file}  JπTz={J2} 0 0  "
          f"use_h3={use_h3}  max_iter={max_iter}  state_want={state_want}",
          flush=True)

    eom = EOM(Hs, ref_file, J2, 0, 0)
    eom.SetUseRdm3(False)
    eom.SetArnoldiUseProjection(use_projection)
    eom.SetArnoldiUseH3(use_h3)
    eom.SetArnoldiPrintTiming(print_timing)

    run_result = eom.Run(max_iter, state_want)
    result = run_result.arnoldi
    eref = run_result.eref
    print(f"\n=== He8 MR-EOM {name} ===", flush=True)
    print(f"E_ref (<H>_RDM) = {eref:.8f} MeV", flush=True)
    print(f"stop={result.stop_reason}  steps={result.steps}  "
          f"converged={result.converged}", flush=True)
    for k, e in enumerate(result.energies):
        print(f"  E({k}): excitation={e:.6f}  absolute={e + eref:.6f} MeV",
              flush=True)
    dump_result(name, spec, eref, run_result)


if __name__ == "__main__":
    main()
