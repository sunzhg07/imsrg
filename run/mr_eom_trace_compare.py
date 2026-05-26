#!/usr/bin/env python3
from pyIMSRG import *


def build_problem():
    emax = 3
    ref = "He8"
    val = "p-shell"
    core_generator = "atan"
    valence_generator = "shell-model-atan"
    smax_core = 50
    smax_valence = 100

    f2b = "input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
    f2e1, f2e2, f2l = 14, 28, 14
    f3b = "input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
    f3e1, f3e2, f3e3 = 18, 36, 18
    hw = 16
    mode3n = "no2b"
    lecs = "EM1820"

    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)
    rw = ReadWrite()

    particle_rank = 3 if f3b != "none" else 2
    H = Operator(ms, 0, 0, 0, particle_rank)

    if lecs == "Minnesota":
        H += OperatorFromString(ms, "VMinnesota")
    else:
        rw.ReadBareTBME_Darmstadt(f2b, H, f2e1, f2e2, f2l)
        if f3b != "none":
            if mode3n == "no2b":
                H.ThreeBody.SetMode("no2b")
                H.ThreeBody.ReadFile([f3b], [f3e1, f3e2, f3e3])
            else:
                rw.Read_Darmstadt_3body(f3b, H, f3e1, f3e2, f3e3)

    H += OperatorFromString(ms, "Trel")

    hf = HartreeFock(H)
    hf.Solve()
    hf.PrintSPEandWF()

    HNO = hf.GetNormalOrderedH(2)
    imsrgsolver = IMSRGSolver(HNO)
    imsrgsolver.SetMethod("magnus")
    imsrgsolver.SetGenerator(core_generator)
    imsrgsolver.SetSmax(smax_core)
    imsrgsolver.Solve()
    imsrgsolver.SetGenerator(valence_generator)
    imsrgsolver.SetSmax(smax_valence)
    imsrgsolver.Solve()

    Hs = imsrgsolver.GetH_s()
    Hs.UndoNormalOrdering()
    ms2 = ModelSpace(emax, "He4", val)
    Hs.SetModelSpace(ms2)
    Hs.DoNormalOrdering()
    return Hs, ms2


def make_eom(Hs):
    eom = EOM(Hs, "he8.ref", 0, 0, 0)
    eom.ConstructConfigs()
    eom.ConstructNormMatrix()
    eom.ConstructProjectMatrix()
    return eom


def main():
    Hs, ms2 = build_problem()
    unt = UnitTest(ms2)
    h1 = unt.RandomOp(ms2, 0, 0, 0, 2, 1)

    seed_eom = make_eom(Hs)
    chi_seed = seed_eom.GetVSEOM_ladder_multiref(h1, 1)
    diff = seed_eom.CompareArnoldiHallBuild(chi_seed, 20, 1e-10)

    if diff.found:
        print("first divergence:")
        print("  step =", diff.step, "i =", diff.i, "j =", diff.j)
        print("  new hall   =", diff.hall_new)
        print("  old hall   =", diff.hall_old)
        print("  delta hall =", diff.delta_hall)
        print("  new h1_sym =", diff.h1_sym_new)
        print("  old h1_sym =", diff.h1_sym_old)
        print("  new h2_cross =", diff.h2_cross_new)
        print("  old h2_cross =", diff.h2_cross_old)
        print("  new h2_sym =", diff.h2_sym_new)
        print("  old h2_sym =", diff.h2_sym_old)
        print("  new h2_diag_i/j =", diff.h2_diag_i_new, diff.h2_diag_j_new)
        print("  old h2_diag_i/j =", diff.h2_diag_i_old, diff.h2_diag_j_old)
    else:
        print("no per-entry divergence found in the traced hall build")

    print("max abs diff =", diff.max_abs_diff)
    print("max rel diff =", diff.max_rel_diff)


if __name__ == "__main__":
    main()