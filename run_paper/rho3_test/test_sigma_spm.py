#!/usr/bin/env python3
"""Does the matvec depend on using S+ or S- ?  sigma(S) = ladder_{+1}([H, S]).
S+ - S- = 2 chi^dagger, so the difference is the pool part of [H, chi^dagger]_{<=2b}.
Exact theory: [H,chi^dag]|Psi> = 0, but only for the full (all-rank) commutator."""
import os
import sys, os
import numpy as np
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "build")); sys.path.insert(0, os.path.join(ROOT, "run_paper"))
from pyIMSRG import EOM, ModelSpace, Operator, ReadWrite, Commutator, UnitTest
from he8_io import hs_path, RUN2
nuc = sys.argv[1]; low = nuc.lower()
J2, par, itz = (int(x) for x in sys.argv[2:5])
ms = ModelSpace(3, "He4", "p-shell"); ms.SetHbarOmega(16)
Hs = Operator(ms, 0, 0, 0, 3); ReadWrite().ReadOperator(Hs, hs_path(nuc, 16, 3, val="p-shell"))
eom = EOM(Hs, os.path.join(RUN2, f"{low}_rho3", f"{low}_gs_e3.ref"), J2, par, itz)
if itz: eom.SetDeltaTz(-itz)
eom.SetUseRdm3(True); eom.SetArnoldiUseProjection(False)
eom.PrepareHamiltonianForArnoldi(); H = eom.GetHsFull()
tensor = (J2, par, itz) != (0, 0, 0)
(eom.ConstructConfigs_tensor if tensor else eom.ConstructConfigs)()
(eom.ConstructNormMatrix_tensor if tensor else eom.ConstructNormMatrix)()
rng = np.random.default_rng(5)
def rnd():
    op = Operator(ms, J2, itz, par, 2); eom.UnflattenOperator(op, rng.normal(size=eom.eom_dims))
    return eom.GetVSEOM_ladder_multiref(op, 1)          # S+ storage
A, B = rnd(), rnd()
Bm = eom.GetVSEOM_ladder_multiref(B, -1)                 # S- of the same chi
sig_p = eom.GetVSEOM_ladder_multiref(Commutator.Commutator(H, B), 1)
sig_m = eom.GetVSEOM_ladder_multiref(Commutator.Commutator(H, Bm), 1)
vp, vm = eom.FlattenOperator(sig_p), eom.FlattenOperator(sig_m)
print(f"{nuc} J2={J2} P={par} itz={itz} dim={eom.eom_dims}")
print(f"  |sigma(S+)| = {np.linalg.norm(vp):.6e}   |sigma(S-)| = {np.linalg.norm(vm):.6e}   |difference| = {np.linalg.norm(vp-vm):.3e}")
print(f"  H2_AB with S+ : {eom.ComputeNorm(A, sig_p):+.10f}")
print(f"  H2_AB with S- : {eom.ComputeNorm(A, sig_m):+.10f}")
print(f"  H2_AA with S+ : {eom.ComputeNorm(A, eom.GetVSEOM_ladder_multiref(Commutator.Commutator(H, A), 1)):+.10f}")
print(f"  H2_AA with S- : {eom.ComputeNorm(A, eom.GetVSEOM_ladder_multiref(Commutator.Commutator(H, eom.GetVSEOM_ladder_multiref(A,-1)), 1)):+.10f}")
