#!/usr/bin/env python3
"""O14 (0+) -> N14 (1+) Gamow-Teller in the p-shell.

One VS-IMSRG decoupling (reference N14, valence p-shell, inert He4 core)
gives the interaction for both nuclei and the consistently evolved GT
operator <f||sigma tau_pm||i> (no g_A). Then KSHELL:
  O14 0+  (6 valence protons, 4 neutrons, 2M=0)
  N14 1+  (5+5, 2M=2: lowest states with J>=1)
  transit O14 -> N14: KSHELL bare B(GT) and the charge-changing OBTDs.

Outputs:
  hs_cache/Hs_N14_pshell_{LECs}_hw16_emax{E}.{dat,snt}
  hs_cache/GT_N14_pshell_{LECs}_hw16_emax{E}.dat   (evolved GT, He4-core NO)
  run_paper/kshell_o14n14/e{E}/   KSHELL inputs, logs, transit log

From run_paper/:
  python3 o14_n14_gt.py            # EMAX=3
  SKIP_IMSRG=1 python3 o14_n14_gt.py   # reuse Hs/.snt, rerun KSHELL only
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time

from he8_io import BUILD, PAPER, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import (  # noqa: E402
    HartreeFock,
    IMSRGSolver,
    ModelSpace,
    Operator,
    OperatorFromString,
    ReadWrite,
)

emax = int(os.environ.get("EMAX", "3"))
hw = 16
ref = "N14"
val = "p-shell"
inert_core = "He4"
f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"

Hs_file = hs_path(ref, hw, emax, val=val)
snt_src = Hs_file[: -len(".dat")] + ".snt"
GT_file = os.path.join(os.path.dirname(Hs_file), "GT" + os.path.basename(Hs_file)[len("Hs"):])

WORK = os.path.join(PAPER, "kshell_o14n14", f"e{emax}")
KSHELL_BIN = "/Users/wolf/work/imsrg/kshell-20230714-ver4/bin"
KSHELL_EXE = os.path.join(PAPER, "kshell_he8", "kshell.exe")
TRANSIT_EXE = os.path.join(PAPER, "kshell_he8", "transit.exe")
SNT = "pshell.snt"
NUCLEI = {  # name: (valence Z, valence N, 2M, n_eigen)
    "O14": (6, 4, 0, 3),   # Z=8,N=6 minus He4 core (2p,2n)
    "N14": (5, 5, 2, 8),   # Z=7,N=7 minus He4 core
}
H_LINE = re.compile(r"^\s*(\d+)\s+<H>:\s+([-+0-9.eE]+)\s+<JJ>:\s+[-+0-9.eE]+\s+J:\s+(\d+)/2")


def run_imsrg():
    t0 = time.time()
    rw = ReadWrite()
    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)
    H = Operator(ms, 0, 0, 0, 3)
    rw.ReadBareTBME_Darmstadt(f2b, H, 14, 28, 14)
    H.ThreeBody.SetMode("no2b")
    H.ThreeBody.ReadFile([f3b], [18, 36, 18])
    H += OperatorFromString(ms, "Trel")
    hf = HartreeFock(H)
    hf.Solve()
    HNO = hf.GetNormalOrderedH(2)
    print(f"E_HF (ensemble ref {ref}) = {HNO.ZeroBody:.6f} MeV", flush=True)
    GT = OperatorFromString(ms, "GamowTeller")
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
    print(f"VS Hs (NO wrt {ref}): E0 = {Hs.ZeroBody:.6f} MeV  ({time.time() - t0:.0f} s)", flush=True)

    ms_core = ModelSpace(emax, inert_core, val)
    ms_core.SetHbarOmega(hw)
    Hc = Hs.UndoNormalOrdering().DoNormalOrderingCore()
    Hc.SetModelSpace(ms_core)
    GTs = solver.Transform(GT)
    GTc = GTs.UndoNormalOrdering().DoNormalOrderingCore()
    GTc.SetModelSpace(ms_core)
    rw.WriteOperator(Hc, Hs_file)
    rw.WriteTokyo(Hc, snt_src, "")
    rw.WriteOperator(GTc, GT_file)
    print(f"core-NO ({inert_core}): E_core = {Hc.ZeroBody:.6f} MeV", flush=True)
    print(f"GT: ||1b|| bare={GT.OneBodyNorm():.5f}  evolved={GTc.OneBodyNorm():.5f}  "
          f"||2b|| evolved={GTc.TwoBodyNorm():.5f}", flush=True)
    print(f"wrote {Hs_file}\nwrote {snt_src}\nwrote {GT_file}", flush=True)


def run(cmd, log, stdin=None):
    with open(log, "w") as f:
        rc = subprocess.call(cmd, cwd=WORK, stdout=f, stderr=subprocess.STDOUT,
                             stdin=stdin)
    if rc != 0:
        raise RuntimeError(f"{cmd[0]} failed ({rc}); see {log}")


def kshell_input(path, ptn, wav, mtot, n_eigen):
    with open(path, "w") as f:
        f.write("&input\n"
                f'  fn_int = "{SNT}"\n  fn_ptn = "{ptn}"\n  fn_save_wave = "{wav}"\n'
                f"  mtot = {mtot}\n  n_eigen = {n_eigen}\n  n_restart_vec = 10\n"
                "  max_lanc_vec = 100\n  maxiter = 300\n  is_double_j = .true.\n"
                "  mass = 14\n/\n")


def run_kshell():
    os.makedirs(WORK, exist_ok=True)
    shutil.copy2(snt_src, os.path.join(WORK, SNT))
    states = {}
    for nuc, (zv, nv, mtot, n_eigen) in NUCLEI.items():
        ptn = f"{nuc.lower()}.ptn"
        if not os.path.isfile(os.path.join(WORK, ptn)):
            with open(os.path.join(WORK, "gen.in"), "w") as f:
                f.write("\n")
            with open(os.path.join(WORK, "gen.in")) as f:
                run([sys.executable, os.path.join(KSHELL_BIN, "gen_partition.py"),
                     SNT, ptn, str(zv), str(nv), "+"], os.path.join(WORK, f"gen_{nuc}.log"), stdin=f)
        wav = f"{nuc.lower()}.wav"
        inp = f"{nuc.lower()}.input"
        kshell_input(os.path.join(WORK, inp), ptn, wav, mtot, n_eigen)
        log = os.path.join(WORK, f"log_{nuc}.txt")
        run([KSHELL_EXE, inp], log)
        st = []
        with open(log) as f:
            for line in f:
                m = H_LINE.match(line)
                if m:
                    st.append((int(m.group(1)), float(m.group(2)), int(m.group(3))))
        states[nuc] = st
        print(f"KSHELL {nuc} (2M={mtot}):", flush=True)
        for i, e, tj in st:
            print(f"  {i}: J={tj / 2:g}+  E_SM={e:.6f}  Ex={e - st[0][1]:.4f}", flush=True)

    with open(os.path.join(WORK, "gt_transit.input"), "w") as f:
        f.write("&input\n"
                f'  fn_int = "{SNT}"\n'
                '  fn_ptn_l = "o14.ptn"\n  fn_load_wave_l = "o14.wav"\n'
                '  fn_ptn_r = "n14.ptn"\n  fn_load_wave_r = "n14.wav"\n'
                "  is_obtd = .true.\n  mass = 14\n/\n")
    tlog = os.path.join(WORK, "transit_O14_N14.log")
    run([TRANSIT_EXE, "gt_transit.input"], tlog)
    print(f"transit O14 -> N14: {tlog}", flush=True)
    return states, tlog


def main():
    print(f"O14 -> N14 GT  emax={emax}  ref={ref}  val={val}  core={inert_core}", flush=True)
    if os.environ.get("SKIP_IMSRG", "0") == "0" or not os.path.isfile(snt_src):
        run_imsrg()
    run_kshell()


if __name__ == "__main__":
    main()
