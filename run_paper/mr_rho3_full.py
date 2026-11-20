#!/usr/bin/env python3
"""MR-EOM full pipeline (NUCLEUS=He8 or Be8), generalized from he8_rho3_full.py.

NUCLEUS=Be8: VS-IMSRG with reference Be8 (p-shell, He4 core), KSHELL 2p2n.
Original He8 notes follow.

He8 MR-EOM, full pipeline, with and without the rho3 metric terms.

PHASE=imsrg  VS-IMSRG (ref He8, p-shell, inert He4 core), EM1.8/2.0 hw16,
             same flow as he8_mr_eom.py -> hs_cache/Hs_He8_pshell_*_emax{E}.dat/.snt
PHASE=kshell KSHELL (kshell_threebody build): He8 M=0 states -> 1b/2b/3b
             densities of 0+_1, 0+_2, 2+_1, 2+_2 -> run2/he8_rho3/he8_{key}_e{E}.ref
             (3b trace checked against C(4,3)); Li8 (1p3n) spectrum for the
             charge-changing comparison.
PHASE=eom    the eight he8_mr_eom.py jobs (mr0gs mr0x t02 t02x mr2gs mr2x t20 t20x)
             and He8 0+ -> Li8 (itz=1, J=1+, 2+, 3+), each with use_rdm3 off and on.
             Prints E_x and the metric spectrum (negative modes).
PHASE=all    (default) all three, skipping finished steps.

  JOBS="mr0gs t02" RDM3="0 1" MAX_ITER=60 python3 he8_rho3_full.py
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

emax = int(os.environ.get("EMAX", "3"))
hw = 16
nucleus = os.environ.get("NUCLEUS", "He8")
low = nucleus.lower()
val = "p-shell"
inert_core = "He4"
mass = {"He8": 8, "Be8": 8, "He6": 6}[nucleus]
NVAL = {"He8": 4, "Be8": 4, "He6": 2}[nucleus]  # valence nucleons above He4
# charge-changing partner (n -> p) and its p-shell partition
CX = {"He8": "li8", "He6": "li6"}.get(nucleus)
PTN_TEXT = {
    "he6": "# partition file  Z=0  N=2  parity=+1\n 0 2 1\n# num. of  proton partition, neutron partition\n 1 3\n# proton partition\n     1     0  0\n# neutron partition\n     1     0  2\n     2     1  1\n     3     2  0\n# partition of proton and neutron\n3\n    1     1\n    1     2\n    1     3\n",
    "li6": "# partition file  Z=1  N=1  parity=+1\n 1 1 1\n# num. of  proton partition, neutron partition\n 2 2\n# proton partition\n     1     0  1\n     2     1  0\n# neutron partition\n     1     0  1\n     2     1  0\n# partition of proton and neutron\n4\n    1     1\n    1     2\n    2     1\n    2     2\n",
}
Hs_file = hs_path(nucleus, hw, emax, val=val)
snt_src = Hs_file[: -len(".dat")] + ".snt"

KBIN = os.environ.get("KSHELL3_BIN", os.path.join(PAPER, "rho3_test", "kshell", "bin"))
WORK = os.path.join(PAPER, f"kshell_{low}_rho3", f"e{emax}")
PTN_DIR = os.path.join(PAPER, "kshell_he8_rho3", "e3")  # he8.ptn, li8.ptn (p-shell)
BE8_PTN = os.path.join(PAPER, "kshell_be8", "be8_p.ptn")
SNT = f"{low}.snt"
REFDIR = os.path.join(RUN2, f"{low}_rho3")
KJSON = os.path.join(REFDIR, f"{low}_kshell_e{emax}.json")
TDTD2REF = os.path.join(PAPER, "rho3_test", "tdtd2ref.py")
# reference states: key -> (2J, n-th state of that J)
REFS = {"gs": (0, 1), "exc0": (0, 2), "2p_gs": (4, 1), "2p_exc": (4, 2)}
# he8_mr_eom.py jobs: name -> (ref key, J2, parity, itz)
JOBS = {
    "mr0gs": ("gs", 0, 0, 0), "mr0x": ("exc0", 0, 0, 0),
    "t02": ("gs", 2, 0, 0), "t02x": ("exc0", 2, 0, 0),
    "mr2gs": ("2p_gs", 0, 0, 0), "mr2x": ("2p_exc", 0, 0, 0),
    "t20": ("2p_gs", 2, 0, 0), "t20x": ("2p_exc", 2, 0, 0),
    # charge-changing He8 0+ -> Li8 (n -> p); X^dag kills He8 (no valence p)
    "cx1": ("gs", 1, 0, 1), "cx2": ("gs", 2, 0, 1), "cx3": ("gs", 3, 0, 1),
}
if nucleus == "He6":
    JOBS["cx0"] = ("gs", 0, 0, 1)  # isobaric analog in Li6
if CX is None:
    JOBS = {k: v for k, v in JOBS.items() if not k.startswith("cx")}

H_LINE = re.compile(r"^\s*(\d+)\s+<H>:\s+([-+0-9.eE]+)\s+<JJ>:\s+[-+0-9.eE]+\s+J:\s+(\d+)/2\s+prty\s+(-?\d+)")
ZB_LINE = re.compile(r"Zero body term:\s*([-+0-9.eE]+)")


def ref_path(key):
    return os.path.join(REFDIR, f"{low}_{key}_e{emax}.ref")


def core_ms():
    ms = ModelSpace(emax, inert_core, val)
    ms.SetHbarOmega(hw)
    return ms


def read_hs(ms):
    Hs = Operator(ms, 0, 0, 0, 3)
    ReadWrite().ReadOperator(Hs, Hs_file)
    return Hs


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
    print(f"[imsrg] E_HF (ref {nucleus}) = {HNO.ZeroBody:.6f}", flush=True)
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
    ms_core = core_ms()  # keep alive: the Operator only stores a pointer
    Hcore.SetModelSpace(ms_core)
    print(f"[imsrg] E_core(He4) = {Hcore.ZeroBody:.6f}", flush=True)
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


def kshell_states(tag, ptn, n_eigen):
    nml(os.path.join(WORK, f"{tag}.input"), [
        ("fn_int", f'"{SNT}"'), ("fn_ptn", f'"{ptn}"'), ("fn_save_wave", f'"{tag}.wav"'),
        ("mtot", 0), ("n_eigen", n_eigen), ("n_restart_vec", 10), ("max_lanc_vec", 100),
        ("maxiter", 300), ("tol", "1.d-10"), ("is_double_j", ".false."), ("mass", mass)])
    log = os.path.join(WORK, f"log_{tag}.txt")
    run_cmd([os.path.join(KBIN, "kshell.exe"), f"{tag}.input"], log)
    states = []
    for line in open(log):
        m = H_LINE.match(line)
        if m:
            states.append({"index": int(m.group(1)), "E_sm": float(m.group(2)),
                           "twoJ": int(m.group(3)), "prty": int(m.group(4))})
    if not states:
        raise RuntimeError(f"no KSHELL states in {log}")
    return states


def trbtd_trace(lines):
    """sum over diagonal TRBTD of value/g; = C(N,3) for the scalar density."""
    t = 0.0
    for x in lines:
        a, b, c, d, e, f, jab, jde, j3 = map(int, x.split()[1:10])
        if (a, b, c, jab) == (d, e, f, jde):
            t += float(x.split()[10]) / {3: 1, 2: 2, 1: 6}[len({a, b, c})]
    return t


def phase_kshell():
    if os.path.isfile(KJSON) and all(os.path.isfile(ref_path(k)) for k in REFS):
        print(f"[kshell] have {KJSON}", flush=True)
        return
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(REFDIR, exist_ok=True)
    shutil.copy2(snt_src, os.path.join(WORK, SNT))
    if nucleus == "He8":
        for p in ("he8.ptn", "li8.ptn"):
            if os.path.abspath(PTN_DIR) != os.path.abspath(WORK):
                shutil.copy2(os.path.join(PTN_DIR, p), os.path.join(WORK, p))
    elif nucleus == "Be8":
        shutil.copy2(BE8_PTN, os.path.join(WORK, f"{low}.ptn"))
    else:
        for tag in (low, CX):
            with open(os.path.join(WORK, f"{tag}.ptn"), "w") as f:
                f.write(PTN_TEXT[tag])
    with open(os.path.join(WORK, SNT)) as f:
        e_core = float(ZB_LINE.search("".join(next(f) for _ in range(8))).group(1))
    he8 = kshell_states(low, f"{low}.ptn", {"He8": 5, "Be8": 15, "He6": 5}[nucleus])
    li8 = kshell_states(CX, f"{CX}.ptn", {"li8": 12, "li6": 6}[CX]) if CX else []
    print(f"[kshell] E_core={e_core:.6f}", flush=True)
    for name, sts in ((nucleus, he8), (str(CX).capitalize(), li8)):
        for st in sts:
            print(f"  {name} {st['index']:3d} J={st['twoJ'] / 2:g}  E_SM={st['E_sm']:.5f}  "
                  f"E-E(He8 gs)={st['E_sm'] - he8[0]['E_sm']:.5f}", flush=True)
    ms = core_ms()
    Hs = read_hs(ms)
    recs = {}
    for key, (twoJ, nth) in REFS.items():
        cand = [s for s in he8 if s["twoJ"] == twoJ]
        if len(cand) < nth:
            raise RuntimeError(f"no {nth}-th state with 2J={twoJ}")
        st = cand[nth - 1]
        i = st["index"]
        base = [("fn_int", f'"{SNT}"'), ("fn_ptn_l", f'"{low}.ptn"'), ("fn_ptn_r", f'"{low}.ptn"'),
                ("fn_load_wave_l", f'"{low}.wav"'), ("fn_load_wave_r", f'"{low}.wav"'),
                ("n_eig_lr_pair", f"{i}, {i}"), ("mass", mass)]
        nml(os.path.join(WORK, f"tr12_{key}.input"),
            base + [("is_tbtd", ".true."), ("is_obtd", ".true.")])
        tr12 = os.path.join(WORK, f"transit_{key}.log")
        run_cmd([os.path.join(KBIN, "transit.exe"), f"tr12_{key}.input"], tr12)
        nml(os.path.join(WORK, f"tr3_{key}.input"), base + [("is_tdtd", ".true.")])
        tr3 = os.path.join(WORK, f"transit_{key}_tdtd.log")
        run_cmd([os.path.join(KBIN, "transit.exe"), f"tr3_{key}.input"], tr3)

        eom = EOM(Hs, 0, 0, 0)
        eom.SetUseRdm3(False)
        eom.LoadKshellRdm(tr12, os.path.join(WORK, SNT), i, i)
        e_rho = eom.GetVSEOM_Overlap_multiref(Hs)
        ref12 = os.path.join(WORK, f"{low}_{key}_12.ref")
        eom.WriteTdm(eom.rdm, ref12)
        del eom
        ref3 = os.path.join(WORK, f"{low}_{key}_3.ref")
        subprocess.check_call([sys.executable, TDTD2REF, os.path.join(WORK, SNT), tr3, ref3,
                               str(i), str(i), str(twoJ)])
        l12 = open(ref12).read().split("\n")
        l3 = [l for l in open(ref3).read().split("\n") if l.strip()]
        norb = int(l12[1])
        if [l.split()[1:5] for l in l12[2:2 + norb]] != [l.split()[1:5] for l in l3[2:2 + norb]]:
            raise RuntimeError("orbit tables differ between WriteTdm and the .snt")
        while l12 and not l12[-1].strip():
            l12.pop()
        if l12[-1].strip() != "0":
            raise RuntimeError("expected an empty TRBTD block at the end of " + ref12)
        trb = [l for l in l3 if l.startswith("TRBTD")]
        tr = trbtd_trace(trb)
        with open(ref_path(key), "w") as f:
            f.write("\n".join(l12[:-1] + [str(len(trb))] + trb) + "\n")
        e_abs = st["E_sm"] + e_core
        print(f"[kshell] ref {key}: state {i} J={twoJ / 2:g}  <H>_rho={e_rho:.6f}  "
              f"KSHELL+E_core={e_abs:.6f}  delta={e_rho - e_abs:.1e}  "
              f"rho3 trace={tr:.6f} (C({NVAL},3)={math.comb(NVAL, 3)})  -> {ref_path(key)}",
              flush=True)
        if abs(e_rho - e_abs) > 1e-4 or abs(tr - math.comb(NVAL, 3)) > 1e-4:
            raise RuntimeError(f"density check failed for {key}")
        recs[key] = {"state": i, "twoJ": twoJ, "E_sm": st["E_sm"], "E_abs": e_abs,
                     "E_rho": e_rho, "rho3_trace": tr}
    with open(KJSON, "w") as f:
        json.dump({"E_core": e_core, "he8": he8, "li8": li8, "refs": recs}, f, indent=2)


# ----------------------------------------------------------------- EOM
def phase_eom():
    ks = json.load(open(KJSON))
    e_gs = ks["he8"][0]["E_sm"]
    max_iter = int(os.environ.get("MAX_ITER", "60"))
    state_want = int(os.environ.get("STATE_WANT", "8"))
    flags = [f == "1" for f in os.environ.get("RDM3", "0 1").split()]
    names = os.environ.get("JOBS", " ".join(JOBS)).split()
    ms = core_ms()
    Hs = read_hs(ms)
    allres = {}
    for name in names:
        key, J2, par, itz = JOBS[name]
        rec = ks["refs"][key]
        for flag in flags:
            t0 = time.time()
            print(f"\n===== {nucleus} {name}: ref {key} (J={rec['twoJ'] / 2:g})  "
                  f"EOM J2={J2} P={par} itz={itz}  emax={emax}  use_rdm3={flag} =====",
                  flush=True)
            eom = EOM(Hs, ref_path(key), J2, par, itz)
            if itz:
                eom.SetDeltaTz(-itz)  # n -> p (He8 -> Li8)
            eom.SetUseRdm3(flag)
            eom.SetArnoldiUseProjection(True)
            eom.SetArnoldiUseH3(True)
            try:
                rr = eom.Run(max_iter, state_want)
            except Exception as exc:  # keep going with the other jobs
                print(f"FAILED: {exc}", flush=True)
                allres[f"{name}_{'rdm3' if flag else 'rdm2'}"] = {"error": str(exc)}
                continue
            ar = rr.arnoldi
            en = [float(x) for x in ar.energies]
            N = np.array(eom.Nkernel)
            w = np.linalg.eigvalsh(0.5 * (N + N.T))
            tol = 1e-8 * max(1.0, abs(w).max())
            print(f"E_ref=<H>_rho={rr.eref:.6f}  steps={ar.steps}  converged={ar.converged}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
            print(f"metric: dim={len(w)}  min={w.min():+.4e}  #neg={(w < -tol).sum()}  "
                  f"#null={(abs(w) <= tol).sum()}", flush=True)
            for k, e in enumerate(en):
                extra = ""
                if itz:
                    extra = f"   E - E(He8 gs) = {e + rr.eref - (e_gs + ks['E_core']):9.4f}"
                print(f"  E({k}) = {e:10.4f}   absolute {e + rr.eref:10.4f}{extra}", flush=True)
            chk = []
            for k in range(len(en)):
                rk = ar.ritz[k]
                chk.append((float(eom.ExpectationValue(rk)), float(eom.ComputeNorm(rk, rk))))
            for k, (ev, nk) in enumerate(chk):
                print(f"  check E({k}): Ritz={en[k]:.6f}  <H>/N={ev:.6f}  diff={ev - en[k]:+.1e}  N={nk:.8f}", flush=True)
            allres[f"{name}_{'rdm3' if flag else 'rdm2'}"] = {
                "ref": key, "J2": J2, "itz": itz, "E_ref": rr.eref, "E_x": en,
                "E_abs": [e + rr.eref for e in en],
                "steps": int(ar.steps), "converged": bool(ar.converged),
                "stop_reason": str(ar.stop_reason),
                "residuals": [float(x) for x in ar.residuals],
                "expect": [c[0] for c in chk], "ritz_norm": [c[1] for c in chk],
                "metric_dim": int(len(w)), "metric_min": float(w.min()),
                "metric_nneg": int((w < -tol).sum()), "metric_nnull": int((abs(w) <= tol).sum())}
            del eom
    print("\nKSHELL p-shell (E - E(He8 gs), MeV):", flush=True)
    for name in ("he8", "li8") if CX else ("he8",):
        print(f"  {name}: " + "  ".join(f"{s['twoJ'] / 2:g}+ {s['E_sm'] - e_gs:.3f}"
                                         for s in ks[name]), flush=True)
    path = eom_result_path(nucleus, hw, emax, "rho3_full")
    with open(path, "w") as f:
        json.dump({"emax": emax, "kshell": ks, "runs": allres}, f, indent=2)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    ph = os.environ.get("PHASE", "all")
    if ph in ("imsrg", "all"):
        phase_imsrg()
    if ph in ("kshell", "all"):
        phase_kshell()
    if ph in ("eom", "all"):
        phase_eom()
