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
  t02      tensor  0+_gs → 2+     (he8.ref, J=2)
  t02x     tensor  0+_1  → 2+     (he8/he8_exc0.ref, J=2)
  mr2gs    scalar  2+_gs → 2+     (he8/he8_2p_gs.ref)
  mr2x     scalar  2+_1  → 2+     (he8/he8_2p_exc.ref)
  t20      tensor  2+_gs → 0+     (he8/he8_2p_gs.ref, J=2)
  t20x     tensor  2+_1  → 0+     (he8/he8_2p_exc.ref, J=2)
  p00gs    tensor  0+_gs → 0-     (he8.ref, J=0, π=−)
  p00x     tensor  0+_1  → 0-     (he8/he8_exc0.ref, J=0, π=−)
  p20gs    tensor  2+_gs → 0-     (he8/he8_2p_gs.ref, J=2, π=−)
  p20x     tensor  2+_1  → 0-     (he8/he8_2p_exc.ref, J=2, π=−)
  p02gs    tensor  0+_gs → 2-     (he8.ref, J=2, π=−)
  p02x     tensor  0+_1  → 2-     (he8/he8_exc0.ref, J=2, π=−)
  0m       the four 0- jobs
  2m       0+_gs/x → 2-
  all      KSHELL RDMs (if needed) + the eight even-parity EOM jobs
  rdm      KSHELL 0+ / 2+ + transit OBTD/TBTD → OSM .ref (no EOM)
  hs       decouple (if needed) + write Hs .dat and p-shell .snt
  pshell   decouple → write_pshell_hamiltonian (no EOM)
  fci      decouple → write_fci_hamiltonian (no EOM)
  snt      decouple → both writers from the same Hs (no EOM)

From run_paper/:
  PYTHONPATH=../build python3 he8_mr_eom.py JOB=mr0gs
  PYTHONPATH=../build MAX_ITER=100 python3 he8_mr_eom.py JOB=t02
  PHASE=1 python3 he8_mr_eom.py JOB=mr0gs     # rebuild Hs .dat
  EMAX=4 python3 he8_mr_eom.py JOB=hs         # VS-IMSRG + .snt
  EMAX=4 python3 he8_mr_eom.py JOB=rdm        # KSHELL densities
  EMAX=4 python3 he8_mr_eom.py JOB=all        # RDM + eight EOM jobs
  python3 he8_mr_eom.py JOB=pshell            # He4-core p-shell .snt
  python3 he8_mr_eom.py JOB=fci               # vacuum FCI .snt
  python3 he8_mr_eom.py JOB=snt               # both .snt from one Hs
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys

from he8_io import (
    BUILD,
    PAPER,
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
Hs_file = os.environ.get("HS_FILE", hs_path(nucleus, hw, emax, val=val))
if Hs_file and not os.path.isabs(Hs_file):
    Hs_file = os.path.join(RUN2, Hs_file)

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"

KSHELL_HE8 = os.path.join(PAPER, "kshell_he8")
SNT_NAME = "he8.snt"
PTN_NAME = "he8_00.ptn"
H_LINE = re.compile(
    r"^\s*(\d+)\s+<H>:\s+([-+0-9.eE]+)\s+<JJ>:\s+[-+0-9.eE]+\s+J:\s+(\d+)/2"
)
ZB_LINE = re.compile(r"Zero body term:\s*([-+0-9.eE]+)")

RDM_KEYS = ("gs", "exc0", "2p_gs", "2p_exc")

JOBS = {
    "mr0gs": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}0^+$",
        "kind": "mr0",
        "J2": 0,
        "parity": 0,
        "rdm_key": "gs",
    },
    "mr0x": {
        "lab": r"MR $0^+_1{\to}0^+$",
        "kind": "mr0x",
        "J2": 0,
        "parity": 0,
        "rdm_key": "exc0",
    },
    "mr2gs": {
        "lab": r"MR $2^+_{\mathrm{gs}}{\to}2^+$",
        "kind": "mr2",
        "J2": 0,
        "parity": 0,
        "rdm_key": "2p_gs",
    },
    "mr2x": {
        "lab": r"MR $2^+_1{\to}2^+$",
        "kind": "mr2x",
        "J2": 0,
        "parity": 0,
        "rdm_key": "2p_exc",
    },
    "t02": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}2^+$",
        "kind": "t02",
        "J2": 2,
        "parity": 0,
        "rdm_key": "gs",
    },
    "t02x": {
        "lab": r"MR $0^+_1{\to}2^+$",
        "kind": "t02x",
        "J2": 2,
        "parity": 0,
        "rdm_key": "exc0",
    },
    "t20": {
        "lab": r"MR $2^+_{\mathrm{gs}}{\to}0^+$",
        "kind": "t20",
        "J2": 2,
        "parity": 0,
        "rdm_key": "2p_gs",
    },
    "t20x": {
        "lab": r"MR $2^+_1{\to}0^+$",
        "kind": "t20x",
        "J2": 2,
        "parity": 0,
        "rdm_key": "2p_exc",
    },
    "p00gs": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}0^-$",
        "kind": "p00",
        "J2": 0,
        "parity": 1,
        "rdm_key": "gs",
    },
    "p00x": {
        "lab": r"MR $0^+_1{\to}0^-$",
        "kind": "p00x",
        "J2": 0,
        "parity": 1,
        "rdm_key": "exc0",
    },
    "p20gs": {
        "lab": r"MR $2^+_{\mathrm{gs}}{\to}0^-$",
        "kind": "p20",
        "J2": 2,
        "parity": 1,
        "rdm_key": "2p_gs",
    },
    "p20x": {
        "lab": r"MR $2^+_1{\to}0^-$",
        "kind": "p20x",
        "J2": 2,
        "parity": 1,
        "rdm_key": "2p_exc",
    },
    "p02gs": {
        "lab": r"MR $0^+_{\mathrm{gs}}{\to}2^-$",
        "kind": "p02",
        "J2": 2,
        "parity": 1,
        "rdm_key": "gs",
    },
    "p02x": {
        "lab": r"MR $0^+_1{\to}2^-$",
        "kind": "p02x",
        "J2": 2,
        "parity": 1,
        "rdm_key": "exc0",
    },
}

EOM_JOB_ORDER = (
    "mr0gs", "mr0x", "t02", "t02x", "mr2gs", "mr2x", "t20", "t20x",
)
ODD_JOB_ORDER = (
    "p00gs", "p00x", "p20gs", "p20x",
)
P02_JOB_ORDER = (
    "p02gs", "p02x",
)


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


def rdm_path(key: str) -> str:
    """Per-emax OSM .ref written from KSHELL transit (gs, exc0, 2p_gs, 2p_exc)."""
    os.makedirs(os.path.join(RUN2, "he8"), exist_ok=True)
    return os.path.join(RUN2, "he8", f"he8_{key}_e{emax}.ref")


def rdm_paths() -> dict[str, str]:
    return {key: rdm_path(key) for key in RDM_KEYS}


def transit_log(key: str) -> str:
    work = kshell_work()
    if key in ("gs", "exc0"):
        return os.path.join(work, "transit_He8_j0p.log")
    return os.path.join(work, "transit_He8_j4p.log")


def transit_state(key: str) -> int:
    return 1 if key in ("gs", "2p_gs") else 2


def refs_ready() -> bool:
    work = kshell_work()
    return (
        os.path.isfile(os.path.join(work, "transit_He8_j0p.log"))
        and os.path.isfile(os.path.join(work, "transit_He8_j4p.log"))
        and os.path.isfile(os.path.join(work, SNT_NAME))
    )


def kshell_work() -> str:
    d = os.path.join(KSHELL_HE8, f"e{emax}")
    os.makedirs(d, exist_ok=True)
    return d


def write_core_snt(Hs):
    """WriteTokyo of an already He4-core-NO Hs (do not UndoNO again)."""
    path = pshell_snt_path(nucleus, hw, emax)
    ReadWrite().WriteTokyo(Hs, path, "")
    print(f"p-shell core-NO ZeroBody={Hs.ZeroBody:.8f} MeV", flush=True)
    print(f"wrote {path}", flush=True)
    return path


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


def parse_kshell_states(log_path):
    states = []
    with open(log_path) as f:
        for line in f:
            m = H_LINE.match(line)
            if m:
                states.append({
                    "index": int(m.group(1)),
                    "E_sm": float(m.group(2)),
                    "twoJ": int(m.group(3)),
                })
    return states


def write_kshell_input(path, mtot, wav, n_eigen):
    with open(path, "w") as f:
        f.write("&input\n")
        f.write(f'  fn_int = "{SNT_NAME}"\n')
        f.write(f'  fn_ptn = "{PTN_NAME}"\n')
        f.write(f'  fn_save_wave = "{wav}"\n')
        f.write(f"  mtot = {mtot}\n")
        f.write(f"  n_eigen = {n_eigen}\n")
        f.write("  n_restart_vec = 10\n")
        f.write("  max_lanc_vec = 40\n")
        f.write("  maxiter = 200\n")
        f.write("  is_double_j = .true.\n")
        f.write("  is_calc_tbdm = .false.\n")
        f.write("  mass = 8\n")
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
        f.write("  mass = 8\n")
        f.write("/\n")


def run_cmd(cmd, cwd, log_path):
    print("running:", " ".join(cmd), flush=True)
    with open(log_path, "w") as log:
        rc = subprocess.call(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
    print(f"  exit={rc}  log={log_path}", flush=True)
    if rc != 0:
        raise RuntimeError(f"{cmd[0]} failed ({rc}); see {log_path}")
    return rc


def load_kshell_rdm_energy(Hs, transit, snt, dest, state=1):
    """⟨H⟩_ρ from KSHELL transit; write an OSM-style .ref for EOM(Hs, file)."""
    print(f"  load RDM {transit}  state={state}", flush=True)
    eom = EOM(Hs, 0, 0, 0)
    eom.SetUseRdm3(False)
    eom.LoadKshellRdm(transit, snt, state, state)
    print("  compute ⟨H⟩_RDM", flush=True)
    e_rdm = eom.GetVSEOM_Overlap_multiref(Hs)
    print(f"  ⟨H⟩_RDM={e_rdm:.8f}  write {dest}", flush=True)
    eom.WriteTdm(eom.rdm, dest)
    del eom
    return e_rdm


def run_rdm_jobs(Hs):
    """KSHELL 0+ / 2+ + transit OBTD/TBTD → per-emax OSM refs."""
    work = kshell_work()
    snt_src = pshell_snt_path(nucleus, hw, emax)
    if not os.path.isfile(snt_src):
        snt_src = write_core_snt(Hs)
    snt_dst = os.path.join(work, SNT_NAME)
    shutil.copy2(snt_src, snt_dst)
    ptn_src = os.path.join(KSHELL_HE8, PTN_NAME)
    if not os.path.isfile(ptn_src):
        raise FileNotFoundError(ptn_src)
    shutil.copy2(ptn_src, os.path.join(work, PTN_NAME))

    with open(snt_dst) as f:
        head = "".join(next(f) for _ in range(8))
    e0_snt = float(ZB_LINE.search(head).group(1))
    print(f"snt {snt_src}  ZeroBody={e0_snt:.8f} MeV", flush=True)

    kshell_exe = os.environ.get("KSHELL_EXE", os.path.join(KSHELL_HE8, "kshell.exe"))
    transit_exe = os.environ.get("TRANSIT_EXE", os.path.join(KSHELL_HE8, "transit.exe"))
    for exe in (kshell_exe, transit_exe):
        if not os.path.isfile(exe):
            raise FileNotFoundError(exe)

    n_eigen = int(os.environ.get("N_EIGEN", "3"))
    recs = []
    channels = (
        ("0p", 0, (("gs", 1), ("exc0", 2))),
        ("2p", 4, (("2p_gs", 1), ("2p_exc", 2))),
    )
    skip = os.environ.get("SKIP_KSHELL", "0") != "0"
    for tag, mtot, dests in channels:
        wav = f"he8_j{mtot}p.wav"
        ks_in = f"he8_j{mtot}.input"
        tr_in = f"he8_j{mtot}_transit.input"
        ks_log = os.path.join(work, f"log_He8_j{mtot}p.txt")
        tr_log = os.path.join(work, f"transit_He8_j{mtot}p.log")
        write_kshell_input(os.path.join(work, ks_in), mtot, wav, n_eigen)
        write_transit_input(os.path.join(work, tr_in), wav)
        if skip and os.path.isfile(tr_log):
            print(f"SKIP_KSHELL=1  reuse {ks_log} / {tr_log}", flush=True)
        else:
            run_cmd([kshell_exe, ks_in], work, ks_log)
            run_cmd([transit_exe, tr_in], work, tr_log)

        states = parse_kshell_states(ks_log)
        if not states:
            raise RuntimeError(f"no KSHELL states in {ks_log}")
        print(f"KSHELL He8 J={mtot // 2}+ valence eigenvalues:", flush=True)
        for st in states:
            print(f"  state {st['index']}: E_SM={st['E_sm']:.8f}  "
                  f"E_SM+E0={st['E_sm'] + e0_snt:.8f}  "
                  f"Ex={st['E_sm'] - states[0]['E_sm']:.8f}", flush=True)

        by_index = {st["index"]: st for st in states}
        for key, state in dests:
            if state not in by_index:
                raise RuntimeError(f"KSHELL {tag} missing state {state}")
            dest = rdm_path(key)
            e_rdm = load_kshell_rdm_energy(Hs, tr_log, snt_dst, dest, state=state)
            e_abs = by_index[state]["E_sm"] + e0_snt
            dE = e_rdm - e_abs
            print(f"RDM {key} → {dest}", flush=True)
            print(f"  ⟨H⟩_RDM={e_rdm:.8f}  KSHELL+E0={e_abs:.8f}  Δ={dE:.3e}",
                  flush=True)
            recs.append({
                "key": key,
                "tag": tag,
                "file": dest,
                "J": mtot // 2,
                "state": state,
                "E_sm": by_index[state]["E_sm"],
                "E_abs": e_abs,
                "E_rdm": e_rdm,
                "delta": dE,
            })

    jpath = os.path.join(RUN2, "he8", f"he8_kshell_e{emax}.json")
    with open(jpath, "w") as f:
        json.dump({
            "E0_snt": e0_snt,
            "hs_file": Hs_file,
            "emax": emax,
            "states": recs,
        }, f, indent=2)
    print(f"wrote {jpath}", flush=True)
    return recs


def dump_result(name, spec, eref, res):
    energies = [float(x) for x in res.arnoldi.energies]
    rec = {
        "name": name,
        "lab": spec["lab"],
        "kind": spec["kind"],
        "mode": "mr",
        "J2": spec["J2"],
        "parity": int(spec.get("parity", 0)),
        "itz": int(spec.get("itz", 0)),
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


def resolve_ref(spec):
    spec = dict(spec)
    key = spec.get("rdm_key")
    if key:
        spec["ref"] = transit_log(key)
        spec["snt"] = os.path.join(kshell_work(), SNT_NAME)
        spec["rdm_state"] = transit_state(key)
        spec["rdm_format"] = "kshell"
        return spec
    ref_file = spec.get("ref", os.path.join(RUN2, "he8.ref"))
    if not os.path.isabs(ref_file):
        ref_file = os.path.join(RUN2, ref_file)
    spec["ref"] = ref_file
    spec.setdefault("rdm_format", "osm")
    return spec


def resolve_job():
    job = requested_job()
    if job and job in JOBS:
        return job, resolve_ref(JOBS[job])
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


def run_eom_job(name, spec, Hs):
    max_iter = int(os.environ.get("MAX_ITER", "50"))
    state_want = int(os.environ.get("STATE_WANT", "5"))
    use_h3 = os.environ.get("USE_H3", "1") != "0"
    use_projection = os.environ.get("USE_PROJECTION", "1") != "0"
    print_timing = os.environ.get("PRINT_TIMING", "0") != "0"
    J2 = spec["J2"]
    parity = int(spec.get("parity", 0))
    itz = int(spec.get("itz", 0))
    ref_file = spec["ref"]
    if not os.path.isfile(ref_file):
        raise FileNotFoundError(ref_file)

    print(f"JOB={name}  RDM={ref_file}  JπTz={J2} {parity} {itz}  "
          f"use_h3={use_h3}  max_iter={max_iter}  state_want={state_want}",
          flush=True)

    fmt = spec.get("rdm_format", "osm")
    if fmt == "kshell":
        state = int(spec.get("rdm_state", 1))
        eom = EOM(Hs, J2, parity, itz)
        eom.LoadKshellRdm(ref_file, spec.get("snt", ""), state, state)
    else:
        eom = EOM(Hs, ref_file, J2, parity, itz)
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
    rec = dump_result(name, spec, eref, run_result)
    del eom
    return rec


def run_snt_jobs(which: str):
    """Decouple once; write p-shell and/or FCI .snt from that Hs (no EOM)."""
    print(f"JOB={which}  (decouple + write interaction .snt, no EOM)", flush=True)
    ms_vs, Hs, solver, HNO = decouple_hs()
    if which in ("pshell", "snt"):
        write_pshell_hamiltonian(Hs)
    if which in ("fci", "snt"):
        write_fci_hamiltonian(Hs)
    ms_core, Hs_core = to_he4_core_hs(Hs)
    ReadWrite().WriteOperator(Hs_core, Hs_file)
    print(f"wrote {Hs_file}", flush=True)
    del solver, HNO, ms_vs, ms_core, Hs_core


def want_rdm() -> bool:
    flag = os.environ.get("RUN_RDM")
    if flag is None:
        return not refs_ready()
    return flag != "0"


def main():
    job = requested_job()
    if job in ("pshell", "fci", "snt"):
        run_snt_jobs(job)
        return

    ms, Hs = build_hs()

    if job == "hs":
        write_core_snt(Hs)
        return

    do_rdm = job in ("rdm",) or want_rdm()
    if do_rdm:
        write_core_snt(Hs)
        print(f"JOB=rdm  KSHELL He8 0+ / 2+ RDMs  emax={emax}", flush=True)
        run_rdm_jobs(Hs)
        print(flush=True)
        if job == "rdm":
            return

    names = None
    if job in ("all", ""):
        names = list(EOM_JOB_ORDER)
    elif job in ("0m", "odd"):
        names = list(ODD_JOB_ORDER)
    elif job == "2m":
        names = list(P02_JOB_ORDER)
    if names is None:
        name, spec = resolve_job()
        run_eom_job(name, spec, Hs)
        return

    print(f"JOB=all  {len(names)} He8 MR-EOM channels", flush=True)
    recs = []
    for name in names:
        spec = resolve_ref(JOBS[name])
        recs.append(run_eom_job(name, spec, Hs))
        print(flush=True)

    print("=== He8 MR-EOM summary ===", flush=True)
    for rec in recs:
        e0 = rec["E_abs"][0] if rec["E_abs"] else float("nan")
        print(f"  {rec['name']:<6}  E_ref={rec['eref']:10.4f}  "
              f"E0={e0:10.4f}  steps={rec['steps']}  "
              f"{rec['stop_reason']}", flush=True)


if __name__ == "__main__":
    main()
