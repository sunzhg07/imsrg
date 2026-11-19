#!/usr/bin/env python3
"""Ne20 0+ -> 0+ MR-EOM with and without the rho3 metric terms.

PHASE=imsrg  VS-IMSRG (ref Ne20, sd-shell, inert O16 core), EM1.8/2.0 hw16,
             same flow as mg24_vs_imsrg.py -> hs_cache/Hs_Ne20_sdshell_*.dat/.snt
PHASE=kshell KSHELL (kshell_threebody build) Ne20 M=0 states, transit
             OBTD/TBTD (-> WriteTdm .ref) and TDTD of 0+_1 (-> TRBTD, appended)
             -> run2/ne20/ne20_gs_e{E}.ref (1b+2b+3b densities, OSM format)
PHASE=eom    scalar MR-EOM J=0+ Tz=0 (EOM.Run = RunMR, H3 on), use_rdm3 off and
             on; prints E_x, the metric spectrum and the KSHELL 0+ E_x.
PHASE=all    (default) the three in order, skipping finished steps.

From run_paper/:  EMAX=4 python3 ne20_rho3_eom.py
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

from he8_io import BUILD, PAPER, RUN2, eom_result_path, hs_path, read_em18_hw16

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

emax = int(os.environ.get("EMAX", "4"))
hw = 16
nucleus = "Ne20"
val = "sd-shell"
inert_core = "O16"
mass = 20
Hs_file = hs_path(nucleus, hw, emax, val=val)
snt_src = Hs_file[: -len(".dat")] + ".snt"

KBIN = os.environ.get("KSHELL3_BIN", os.path.join(PAPER, "rho3_test", "kshell", "bin"))
WORK = os.path.join(PAPER, "kshell_ne20", f"e{emax}")
PTN_SRC = os.path.join(PAPER, "rho3_test", "ne20", "ne20.ptn")
SNT, PTN = "ne20.snt", "ne20.ptn"
N_EIGEN = int(os.environ.get("N_EIGEN", "12"))
REFDIR = os.path.join(RUN2, "ne20")
REF = os.path.join(REFDIR, f"ne20_gs_e{emax}.ref")
KJSON = os.path.join(REFDIR, f"ne20_kshell_e{emax}.json")

H_LINE = re.compile(r"^\s*(\d+)\s+<H>:\s+([-+0-9.eE]+)\s+<JJ>:\s+[-+0-9.eE]+\s+J:\s+(\d+)/2\s+prty\s+(-?\d+)")
ZB_LINE = re.compile(r"Zero body term:\s*([-+0-9.eE]+)")


# ----------------------------------------------------------------- IMSRG
def phase_imsrg():
    if os.path.isfile(Hs_file) and os.path.isfile(snt_src):
        print(f"[imsrg] have {Hs_file}", flush=True)
        return
    t0 = time.time()
    rw = ReadWrite()
    ms = ModelSpace(emax, nucleus, val)
    ms.SetHbarOmega(hw)
    H = Operator(ms, 0, 0, 0, 3)
    read_em18_hw16(H, rw)
    H += OperatorFromString(ms, "Trel")
    hf = HartreeFock(H)
    hf.Solve()
    HNO = hf.GetNormalOrderedH(2)
    print(f"[imsrg] E_HF (ensemble ref {nucleus}) = {HNO.ZeroBody:.6f}", flush=True)
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
    print(f"[imsrg] VS E0 = {Hs.ZeroBody:.6f}  ({time.time() - t0:.0f} s)", flush=True)
    Hcore = Hs.UndoNormalOrdering().DoNormalOrderingCore()
    ms_core = ModelSpace(emax, inert_core, val)
    ms_core.SetHbarOmega(hw)
    Hcore.SetModelSpace(ms_core)
    print(f"[imsrg] E_core(O16) = {Hcore.ZeroBody:.6f}", flush=True)
    rw.WriteOperator(Hcore, Hs_file)
    rw.WriteTokyo(Hcore, snt_src, "")
    print(f"[imsrg] wrote {Hs_file} and .snt", flush=True)


# ----------------------------------------------------------------- KSHELL
def run_cmd(cmd, log):
    print("[kshell] running", " ".join(cmd), flush=True)
    with open(log, "w") as f:
        rc = subprocess.call(cmd, cwd=WORK, stdout=f, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"{cmd} failed ({rc}), see {log}")


def nml(path, body):
    with open(path, "w") as f:
        f.write("&input\n" + "".join(f"  {k} = {v}\n" for k, v in body) + "/\n")


def read_hs(ms):
    Hs = Operator(ms, 0, 0, 0, 3)
    ReadWrite().ReadOperator(Hs, Hs_file)
    return Hs


def phase_kshell():
    if os.path.isfile(REF) and os.path.isfile(KJSON):
        print(f"[kshell] have {REF}", flush=True)
        return
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(REFDIR, exist_ok=True)
    shutil.copy2(snt_src, os.path.join(WORK, SNT))
    shutil.copy2(PTN_SRC, os.path.join(WORK, PTN))
    with open(os.path.join(WORK, SNT)) as f:
        e_core = float(ZB_LINE.search("".join(next(f) for _ in range(8))).group(1))
    kx, tx = os.path.join(KBIN, "kshell.exe"), os.path.join(KBIN, "transit.exe")
    # (1) all low M=0 states (J=0 ones are the 0+ spectrum), (2) gs only wave
    nml(os.path.join(WORK, "ne20.input"), [
        ("fn_int", f'"{SNT}"'), ("fn_ptn", f'"{PTN}"'), ("fn_save_wave", '"ne20.wav"'),
        ("mtot", 0), ("n_eigen", N_EIGEN), ("n_restart_vec", 20), ("max_lanc_vec", 200),
        ("maxiter", 300), ("tol", "1.d-10"), ("is_double_j", ".false."), ("mass", mass)])
    run_cmd([kx, "ne20.input"], os.path.join(WORK, "log_ne20.txt"))
    states = []
    for line in open(os.path.join(WORK, "log_ne20.txt")):
        m = H_LINE.match(line)
        if m:
            states.append({"index": int(m.group(1)), "E_sm": float(m.group(2)),
                           "twoJ": int(m.group(3)), "prty": int(m.group(4))})
    if not states or states[0]["twoJ"] != 0:
        raise RuntimeError(f"KSHELL ground state not 0+: {states[:2]}")
    gs = states[0]
    print(f"[kshell] E_core={e_core:.6f}", flush=True)
    for st in states:
        print(f"  {st['index']:3d} J={st['twoJ'] / 2:g}  E_SM={st['E_sm']:.5f}  "
              f"E_x={st['E_sm'] - gs['E_sm']:.5f}", flush=True)
    # transit 1b/2b and 3b densities of state 1
    nml(os.path.join(WORK, "tr12.input"), [
        ("fn_int", f'"{SNT}"'), ("fn_ptn_l", f'"{PTN}"'), ("fn_ptn_r", f'"{PTN}"'),
        ("fn_load_wave_l", '"ne20.wav"'), ("fn_load_wave_r", '"ne20.wav"'),
        ("is_tbtd", ".true."), ("is_obtd", ".true."), ("n_eig_lr_pair", "1, 1"),
        ("mass", mass)])
    tr12 = os.path.join(WORK, "transit_ne20.log")
    run_cmd([tx, "tr12.input"], tr12)
    nml(os.path.join(WORK, "tr3.input"), [
        ("fn_int", f'"{SNT}"'), ("fn_ptn_l", f'"{PTN}"'), ("fn_ptn_r", f'"{PTN}"'),
        ("fn_load_wave_l", '"ne20.wav"'), ("fn_load_wave_r", '"ne20.wav"'),
        ("is_tdtd", ".true."), ("n_eig_lr_pair", "1, 1"), ("mass", mass)])
    tr3 = os.path.join(WORK, "transit_ne20_tdtd.log")
    run_cmd([tx, "tr3.input"], tr3)

    # 1b+2b .ref via LoadKshellRdm -> WriteTdm, check <H>_rho
    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    Hs = read_hs(ms)
    eom = EOM(Hs, 0, 0, 0)
    eom.SetUseRdm3(False)
    eom.LoadKshellRdm(tr12, os.path.join(WORK, SNT), 1, 1)
    e_rho = eom.GetVSEOM_Overlap_multiref(Hs)
    ref12 = os.path.join(WORK, "ne20_gs_12.ref")
    eom.WriteTdm(eom.rdm, ref12)
    del eom
    e_abs = gs["E_sm"] + e_core
    print(f"[kshell] <H>_rho={e_rho:.6f}  KSHELL+E_core={e_abs:.6f}  delta={e_rho - e_abs:.2e}",
          flush=True)
    # 3b: TDTD -> TRBTD (orbit table of the .snt), appended to the 1b+2b .ref
    ref3 = os.path.join(WORK, "ne20_gs_3.ref")
    subprocess.check_call([sys.executable, os.path.join(PAPER, "rho3_test", "tdtd2ref.py"),
                           os.path.join(WORK, SNT), tr3, ref3])
    l12 = open(ref12).read().split("\n")
    l3 = [l for l in open(ref3).read().split("\n") if l.strip()]
    norb = int(l12[1])
    orb12 = [l.split()[1:5] for l in l12[2:2 + norb]]
    orb3 = [l.split()[1:5] for l in l3[2:2 + norb]]
    if orb12 != orb3:
        raise RuntimeError(f"orbit tables differ: {orb12} vs {orb3}")
    while l12 and not l12[-1].strip():
        l12.pop()
    if l12[-1].strip() != "0":
        raise RuntimeError("expected an empty TRBTD block at the end of " + ref12)
    trb = [l for l in l3 if l.startswith("TRBTD")]
    with open(REF, "w") as f:
        f.write("\n".join(l12[:-1] + [str(len(trb))] + trb) + "\n")
    print(f"[kshell] wrote {REF}  ({len(trb)} TRBTD)", flush=True)
    with open(KJSON, "w") as f:
        json.dump({"E_core": e_core, "E_rho": e_rho, "states": states}, f, indent=2)


# ----------------------------------------------------------------- EOM
def phase_eom():
    ks = json.load(open(KJSON))
    st = ks["states"]
    e0 = st[0]["E_sm"]
    ks0 = [s["E_sm"] - e0 for s in st if s["twoJ"] == 0][1:]
    max_iter = int(os.environ.get("MAX_ITER", "60"))
    state_want = int(os.environ.get("STATE_WANT", "10"))
    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    Hs = read_hs(ms)
    out = {"emax": emax, "ref": REF, "kshell_0p_Ex": ks0, "runs": {}}
    flags = [f == "1" for f in os.environ.get("RDM3", "0 1").split()]
    for flag in flags:
        t0 = time.time()
        print(f"\n===== Ne20 0+ MR-EOM  emax={emax}  use_rdm3={flag} =====", flush=True)
        eom = EOM(Hs, REF, 0, 0, 0)
        eom.SetUseRdm3(flag)
        eom.SetArnoldiUseProjection(True)
        eom.SetArnoldiUseH3(True)
        rr = eom.Run(max_iter, state_want)
        ar = rr.arnoldi
        en = [float(x) for x in ar.energies]
        N = np.array(eom.Nkernel)
        w = np.linalg.eigvalsh(0.5 * (N + N.T))
        tol = 1e-8 * max(1.0, abs(w).max())
        print(f"E_ref=<H>_rho={rr.eref:.6f}  steps={ar.steps}  converged={ar.converged}  "
              f"({time.time() - t0:.0f} s)", flush=True)
        print(f"metric: dim={len(w)}  min={w.min():+.4e}  max={w.max():+.4e}  "
              f"#neg(<-{tol:.0e})={(w < -tol).sum()}  #|w|<tol={(abs(w) <= tol).sum()}  "
              f"lowest 8: {np.array2string(w[:8], precision=4)}", flush=True)
        for k, e in enumerate(en):
            print(f"  E({k}) = {e:10.4f}   absolute {e + rr.eref:10.4f}", flush=True)
        out["runs"]["rdm3" if flag else "rdm2"] = {
            "E_ref": rr.eref, "E_x": en, "steps": int(ar.steps),
            "converged": bool(ar.converged), "metric_min": float(w.min()),
            "metric_nneg": int((w < -tol).sum()), "metric_lowest": w[:20].tolist()}
        del eom
    print("\nKSHELL 0+ E_x (in-valence):", " ".join(f"{e:.4f}" for e in ks0), flush=True)
    path = eom_result_path(nucleus, hw, emax, "0p_rho3")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    ph = os.environ.get("PHASE", "all")
    if ph in ("imsrg", "all"):
        phase_imsrg()
    if ph in ("kshell", "all"):
        phase_kshell()
    if ph in ("eom", "all"):
        phase_eom()
