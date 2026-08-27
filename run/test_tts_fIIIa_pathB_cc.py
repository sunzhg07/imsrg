#!/usr/bin/env python3
"""f^III_a: CC/Pandya production path vs the AMC-direct W1/W2 + ladder gold.

Production (default) keeps chi^gamma in the Pandya representation and contracts
it with bar{Gamma} there -- the scalar Factorized II_a/II_c structure. The gold
oracle (use_TypeIIIa_slow) is the AMC-direct W1/W2 chi + ordinary-channel
ladder, which is locked against m-scheme.

Usage:
  PYTHONPATH=build python3 run/test_tts_fIIIa_pathB_cc.py [emax=1] [lambdas=0,1,2]
"""

from __future__ import annotations

import sys

from pyIMSRG import *

emax = int(sys.argv[1]) if len(sys.argv) > 1 else 1
lams = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0, 1, 2]
tol = 1e-6
seed = 11

FD = Commutator.FactorizedDoubleCommutator_eths


def run_one(lam: int):
    ms = ModelSpace(emax, "He4", "He4")
    ms.SetHbarOmega(20.0)
    ms.PreCalculateSixJ()
    ms.PreCalculateNineJ()
    ut = UnitTest(ms)
    ut.SetRandomSeed(seed)

    Omega = ut.RandomOp(ms, lam, 0, 0, 2, -1)
    if not Omega.IsReduced():
        Omega.MakeReduced()
    Gamma = ut.RandomOp(ms, 0, 0, 0, 2, +1)

    # Only f^III_a
    FD.SetUse_1b_Intermediates(True)
    FD.SetUse_2b_Intermediates(True)
    FD.SetUse_TypeI_1b(False)
    FD.SetUse_TypeII_1b(False)
    FD.SetUse_TypeIII_1b(False)
    FD.SetUse_TypeIIIa_1b(True)

    out = {}
    for tag, slow in (("gold", True), ("cc", False)):
        FD.SetUse_TypeIIIa_slow(slow)
        Z = Operator(ms, 0, 0, 0, 2)
        Z.SetHermitian()
        FD.comm223_231_st(Omega, Gamma, Z)
        out[tag] = Z

    Zg, Zc = out["gold"], out["cc"]
    ng = Zg.OneBodyNorm()
    D = Operator(Zc)
    D -= Zg
    nd = D.OneBodyNorm()
    rel = nd / ng if ng > 1e-14 else nd
    ok = rel < tol
    print(
        f"lambda={lam}: ||f_gold||={ng:.6g}  ||f_cc||={Zc.OneBodyNorm():.6g}  "
        f"||diff||={nd:.3e}  rel={rel:.3e}  => {'PASS' if ok else 'FAIL'}"
    )
    return ok


print(f"=== f^III_a CC/Pandya vs W1/W2 gold (emax={emax}) ===")
allok = all(run_one(l) for l in lams)
print(f"\nVerdict: {'PASS' if allok else 'FAIL'}")
sys.exit(0 if allok else 1)
