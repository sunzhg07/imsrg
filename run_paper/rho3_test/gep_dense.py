#!/usr/bin/env python3
"""Dense generalized eigenproblem (canonical orthogonalization) for scalar MR-EOM, rho3 on.
  gep.py NUCLEUS key [eps nsub]"""
import os
import sys, os, json, time
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "build")); sys.path.insert(0, os.path.join(ROOT, "run_paper"))
from pyIMSRG import EOM, ModelSpace, Operator, ReadWrite
from he8_io import hs_path, RUN2
nuc, key = sys.argv[1], sys.argv[2]
eps = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
nsub = int(sys.argv[4]) if len(sys.argv) > 4 else 80
low = nuc.lower()
ms = ModelSpace(3, "He4", "p-shell"); ms.SetHbarOmega(16)
Hs = Operator(ms, 0, 0, 0, 3); ReadWrite().ReadOperator(Hs, hs_path(nuc, 16, 3, val="p-shell"))
eom = EOM(Hs, os.path.join(RUN2, f"{low}_rho3", f"{low}_{key}_e3.ref"), 0, 0, 0)
eom.SetUseRdm3(True); eom.SetArnoldiUseH3(True)
eom.PrepareHamiltonianForArnoldi()
eref = eom.GetVSEOM_Overlap_multiref(eom.GetHsFull())
eom.ConstructConfigs(); eom.ConstructNormMatrix()
t0 = time.time()
out = {}
for h3 in (False, True):
    eom.SetArnoldiUseH3(h3)
    r = eom.SolveGEPCanonical(12, eps, nsub)
    en = [float(x) for x in r.energies]
    tag = "H2+H3" if h3 else "H2"
    print(f"GEP {nuc} {key} {tag} eps={eps} nsub={nsub} E_ref={eref:.6f} ({time.time()-t0:.0f}s)")
    print("  dE :", " ".join(f"{e:9.4f}" for e in en))
    print("  abs:", " ".join(f"{e+eref:9.4f}" for e in en), flush=True)
    out[tag] = {"E_x": en, "E_abs": [e + eref for e in en]}
out["E_ref"] = eref
json.dump(out, open(os.path.join(ROOT, "eom_result", f"{nuc}_EM1.8_2.0_hw16_emax3_gep_{key}.json"), "w"), indent=1)
