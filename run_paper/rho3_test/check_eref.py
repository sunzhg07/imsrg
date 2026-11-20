#!/usr/bin/env python3
"""<Hbar>_rho vs KSHELL eigenvalue + E_core, for every reference, rho3 off/on."""
import os
import sys, json, os
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "build")); sys.path.insert(0, os.path.join(ROOT, "run_paper"))
from pyIMSRG import EOM, ModelSpace, Operator, ReadWrite
from he8_io import hs_path, RUN2
for nuc in ("He8", "Be8"):
    low = nuc.lower()
    ms = ModelSpace(3, "He4", "p-shell"); ms.SetHbarOmega(16)
    Hs = Operator(ms, 0, 0, 0, 3)
    ReadWrite().ReadOperator(Hs, hs_path(nuc, 16, 3, val="p-shell"))
    n3 = Hs.ThreeBodyNorm() if Hs.ThreeBody.IsAllocated() else 0.0
    print(f"\n{nuc}: Hs file {os.path.basename(hs_path(nuc,16,3,val='p-shell'))}  particle rank {Hs.GetParticleRank()}  "
          f"3B allocated={Hs.ThreeBody.IsAllocated()}  |H_3b|={n3:.3e}  ZeroBody(He4 core)={Hs.ZeroBody:.6f}")
    ks = json.load(open(os.path.join(RUN2, f"{low}_rho3", f"{low}_kshell_e3.json")))
    print(f"  KSHELL .snt zero-body (E_core) = {ks['E_core']:.6f}")
    print(f"  {'ref':7s} {'J':>2s} {'KSHELL+E_core':>14s} {'<H>_rho (rho3 off)':>19s} {'<H>_rho (rho3 on)':>18s} {'diff':>9s} {'rho3 effect':>12s}")
    for key, rec in ks["refs"].items():
        eom = EOM(Hs, os.path.join(RUN2, f"{low}_rho3", f"{low}_{key}_e3.ref"), 0, 0, 0)
        eom.SetUseRdm3(False); e2 = eom.GetVSEOM_Overlap_multiref(Hs)
        eom.SetUseRdm3(True);  e3 = eom.GetVSEOM_Overlap_multiref(Hs)
        print(f"  {key:7s} {rec['twoJ']//2:2d} {rec['E_abs']:14.6f} {e2:19.6f} {e3:18.6f} {e2-rec['E_abs']:9.1e} {e3-e2:12.1e}")
        del eom
