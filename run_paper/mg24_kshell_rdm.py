#!/usr/bin/env python3
"""Mg24 sd-shell KSHELL states -> reference densities for EOM.

Needs the VS Hs from mg24_vs_imsrg.py (same EMAX). Mirrors run_rdm_jobs in
he8_mr_eom.py:
  KSHELL J=0 (3 states) and J=2 (2 states) with the O16-core sd .snt
  transit OBTD/TBTD
  EOM.LoadKshellRdm  ->  <H>_rho check against E_SM + E_core
  WriteTdm           ->  run2/mg24/mg24_{key}_e{EMAX}.ref (OSM format)

Refs written: gs = 0+_1, exc0 = 0+_2, 2p_gs = 2+_1.
Work dir: run_paper/kshell_mg24/e{EMAX}/ (needs mg24_sd.ptn there).

From run_paper/:
  python3 mg24_kshell_rdm.py               # EMAX=4
  SKIP_KSHELL=1 python3 mg24_kshell_rdm.py # reuse KSHELL/transit logs
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys

from he8_io import BUILD, PAPER, RUN2, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import EOM, ModelSpace, Operator, ReadWrite  # noqa: E402

emax = int(os.environ.get("EMAX", "4"))
hw = 16
nucleus = "Mg24"
val = "sd-shell"
inert_core = "O16"
mass = 24

Hs_file = hs_path(nucleus, hw, emax, val=val)
snt_src = Hs_file[: -len(".dat")] + ".snt"
WORK = os.path.join(PAPER, "kshell_mg24", f"e{emax}")
SNT_NAME = "mg24.snt"
PTN_NAME = "mg24_sd.ptn"
KSHELL_EXE = os.environ.get("KSHELL_EXE", os.path.join(PAPER, "kshell_he8", "kshell.exe"))
TRANSIT_EXE = os.environ.get("TRANSIT_EXE", os.path.join(PAPER, "kshell_he8", "transit.exe"))

H_LINE = re.compile(
    r"^\s*(\d+)\s+<H>:\s+([-+0-9.eE]+)\s+<JJ>:\s+[-+0-9.eE]+\s+J:\s+(\d+)/2"
)
ZB_LINE = re.compile(r"Zero body term:\s*([-+0-9.eE]+)")

# (channel tag, 2*M, n_eigen, [(ref key, state index)])
CHANNELS = (
    ("0p", 0, 3, (("gs", 1), ("exc0", 2))),
    ("2p", 4, 2, (("2p_gs", 1),)),
)


def write_kshell_input(path, mtot, wav, n_eigen):
    with open(path, "w") as f:
        f.write("&input\n")
        f.write(f'  fn_int = "{SNT_NAME}"\n')
        f.write(f'  fn_ptn = "{PTN_NAME}"\n')
        f.write(f'  fn_save_wave = "{wav}"\n')
        f.write(f"  mtot = {mtot}\n")
        f.write(f"  n_eigen = {n_eigen}\n")
        f.write("  n_restart_vec = 10\n")
        f.write("  max_lanc_vec = 100\n")
        f.write("  maxiter = 300\n")
        f.write("  is_double_j = .true.\n")
        f.write("  is_calc_tbdm = .false.\n")
        f.write(f"  mass = {mass}\n")
        f.write("/\n")


def write_transit_input(path, wav):
    with open(path, "w") as f:
        f.write("&input\n")
        f.write(f'  fn_int = "{SNT_NAME}"\n')
        f.write(f'  fn_ptn_l = "{PTN_NAME}"\n')
        f.write(f'  fn_ptn_r = "{PTN_NAME}"\n')
        f.write(f'  fn_load_wave_l = "{wav}"\n')
        f.write(f'  fn_load_wave_r = "{wav}"\n')
        f.write("  is_tbtd = .true.\n")
        f.write("  is_obtd = .true.\n")
        f.write(f"  mass = {mass}\n")
        f.write("/\n")


def run_cmd(cmd, log_path):
    print("running:", " ".join(cmd), flush=True)
    with open(log_path, "w") as log:
        rc = subprocess.call(cmd, cwd=WORK, stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"{cmd[0]} failed ({rc}); see {log_path}")


def parse_states(log_path):
    out = []
    with open(log_path) as f:
        for line in f:
            m = H_LINE.match(line)
            if m:
                out.append({"index": int(m.group(1)), "E_sm": float(m.group(2)),
                            "twoJ": int(m.group(3))})
    return out


def main():
    for p in (Hs_file, snt_src, KSHELL_EXE, TRANSIT_EXE, os.path.join(WORK, PTN_NAME)):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
    shutil.copy2(snt_src, os.path.join(WORK, SNT_NAME))
    with open(os.path.join(WORK, SNT_NAME)) as f:
        e_core = float(ZB_LINE.search("".join(next(f) for _ in range(8))).group(1))
    print(f"Mg24 sd-shell  emax={emax}  E_core({inert_core})={e_core:.6f} MeV", flush=True)

    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    Hs = Operator(ms, 0, 0, 0, 3)
    ReadWrite().ReadOperator(Hs, Hs_file)
    print(f"read {Hs_file}  ZeroBody={Hs.ZeroBody:.6f}", flush=True)

    refdir = os.path.join(RUN2, "mg24")
    os.makedirs(refdir, exist_ok=True)
    skip = os.environ.get("SKIP_KSHELL", "0") != "0"
    recs, spectrum = [], {}
    for tag, mtot, n_eigen, dests in CHANNELS:
        wav = f"mg24_m{mtot}.wav"
        ks_in, tr_in = f"mg24_m{mtot}.input", f"mg24_m{mtot}_transit.input"
        ks_log = os.path.join(WORK, f"log_Mg24_m{mtot}.txt")
        tr_log = os.path.join(WORK, f"transit_Mg24_m{mtot}.log")
        write_kshell_input(os.path.join(WORK, ks_in), mtot, wav, n_eigen)
        write_transit_input(os.path.join(WORK, tr_in), wav)
        if not (skip and os.path.isfile(tr_log)):
            run_cmd([KSHELL_EXE, ks_in], ks_log)
            run_cmd([TRANSIT_EXE, tr_in], tr_log)
        states = parse_states(ks_log)
        if not states:
            raise RuntimeError(f"no KSHELL states in {ks_log}")
        spectrum[tag] = states
        print(f"KSHELL 2M={mtot}:", flush=True)
        for st in states:
            print(f"  state {st['index']}: J={st['twoJ'] / 2:g}  E_SM={st['E_sm']:.6f}  "
                  f"E_SM+E_core={st['E_sm'] + e_core:.6f}", flush=True)
        by_idx = {st["index"]: st for st in states}
        for key, idx in dests:
            st = by_idx[idx]
            dest = os.path.join(refdir, f"mg24_{key}_e{emax}.ref")
            eom = EOM(Hs, 0, 0, 0)
            eom.SetUseRdm3(False)
            eom.LoadKshellRdm(tr_log, os.path.join(WORK, SNT_NAME), idx, idx)
            e_rdm = eom.GetVSEOM_Overlap_multiref(Hs)
            eom.WriteTdm(eom.rdm, dest)
            del eom
            e_abs = st["E_sm"] + e_core
            print(f"RDM {key} (J={st['twoJ'] / 2:g}, state {idx}) -> {dest}\n"
                  f"  <H>_rho={e_rdm:.6f}  KSHELL+E_core={e_abs:.6f}  "
                  f"delta={e_rdm - e_abs:.2e}", flush=True)
            recs.append({"key": key, "J": st["twoJ"] / 2, "state": idx, "file": dest,
                         "transit": tr_log, "E_sm": st["E_sm"], "E_abs": e_abs,
                         "E_rdm": e_rdm, "delta": e_rdm - e_abs})

    jpath = os.path.join(refdir, f"mg24_kshell_e{emax}.json")
    with open(jpath, "w") as f:
        json.dump({"E_core": e_core, "hs_file": Hs_file, "emax": emax,
                   "spectrum": spectrum, "refs": recs}, f, indent=2)
    print(f"wrote {jpath}", flush=True)


if __name__ == "__main__":
    main()
