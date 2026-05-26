#!/usr/bin/env python3
import os
import numpy as np
from pyIMSRG import *

emax = 3
ref = "He4"
val = "p-shell"

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
f3e1, f3e2, f3e3 = 18, 36, 18
hw = 16
mode3n = "no2b"
LECs = "EM1820"


def build_hs():
    ms = ModelSpace(emax, ref, val)
    ms.SetHbarOmega(hw)
    rw = ReadWrite()

    rank_j, parity, rank_Tz, particle_rank = 0, 0, 0, 2
    if f3b != "none":
        particle_rank = 3

    H = Operator(ms, rank_j, parity, rank_Tz, particle_rank)

    if LECs == "Minnesota":
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

    phase = int(os.environ.get("PHASE", "2"))
    if phase == 1:
        imsrgsolver.SetGenerator("atan")
        imsrgsolver.SetSmax(50)
        imsrgsolver.Solve()
        imsrgsolver.SetGenerator("shell-model-atan")
        imsrgsolver.SetSmax(100)
        imsrgsolver.Solve()
        Hs = imsrgsolver.GetH_s()
        Hs.ZeroBody = 0.0
        rw.WriteOperator(Hs, "Hs_back")
    else:
        Hs = Operator(ms, rank_j, parity, rank_Tz, particle_rank)
        rw.ReadOperator(Hs, "Hs_back")

    return ms, Hs


def build_eom(Hs, use_projection, check_expectation, print_timing):
    eom = EOM(Hs, "Be10.ref", 0, 0, 0)
    eom.SetArnoldiUseProjection(use_projection)
    eom.SetArnoldiCheckExpectation(check_expectation)
    eom.SetArnoldiPrintTiming(print_timing)
    return eom


def main():
    use_projection = os.environ.get("USE_PROJECTION", "1") != "0"
    check_expectation = os.environ.get("CHECK_EXPECTATION", "0") != "0"
    print_timing = os.environ.get("PRINT_TIMING", "0") != "0"
    max_iter = int(os.environ.get("MAX_ITER", "80"))
    state_want = int(os.environ.get("STATE_WANT", "3"))

    ms, Hs = build_hs()
    eom = build_eom(Hs, use_projection, check_expectation, print_timing)
    print("use_projection =", use_projection)
    print("check_expectation =", check_expectation)
    print("print_timing =", print_timing)
    print("max_iter =", max_iter, "state_want =", state_want)

    run_result = eom.Run(max_iter, state_want)
    result = run_result.arnoldi

    hall_cpp = np.array(result.hall, dtype=float)
    eigvecs = np.array(result.eigvecs, dtype=float)
    energies = np.array(result.energies, dtype=float)
    ritz = result.ritz
    dim = hall_cpp.shape[0]
    reference_energy = run_result.eref

    print("reference_energy =", reference_energy)
    print("reference_energy_shift =", eom.GetReferenceEnergyShift())

    print("\n=== C++ Arnoldi Result ===")
    print("hall shape =", hall_cpp.shape)
    print("cpp energies =", energies)
    if hall_cpp.shape[0] and eigvecs.size:
        hall_cpp_sym = 0.5 * (hall_cpp + hall_cpp.T)
        hall_is_finite = np.all(np.isfinite(hall_cpp_sym))
        if not hall_is_finite:
            print("skipping cT hall c check: hall contains non-finite entries")
        else:
            for k in range(min(state_want, eigvecs.shape[1], len(energies))):
                coeff = eigvecs[:dim, k]
                if coeff.shape[0] != dim or not np.all(np.isfinite(coeff)):
                    print("state", k)
                    print("  eigval    =", energies[k])
                    print("  cT hall c = skipped: eigenvector contains non-finite entries")
                    continue
                hall_num = float(coeff.T @ hall_cpp_sym @ coeff)
                if not np.isfinite(hall_num):
                    print("state", k)
                    print("  eigval    =", energies[k])
                    print("  cT hall c = skipped: result is non-finite")
                    continue
                print("state", k)
                print("  eigval    =", energies[k])
                print("  cT hall c =", hall_num)
                print("  delta     =", hall_num - energies[k])

    if check_expectation:
        print("\n=== C++ Ritz Vector Expectation Check ===")
        max_delta_e = 0.0
        for k in range(min(state_want, len(ritz), len(energies))):
            norm = eom.ComputeNorm(ritz[k], ritz[k])
            expectation_value = eom.ExpectationValue(ritz[k])
            delta_e = expectation_value - energies[k]
            max_delta_e = max(max_delta_e, abs(delta_e))
            print("state", k)
            print("  eigval             =", energies[k])
            print("  norm               =", norm)
            print("  ExpectationValue   =", expectation_value)
            print("  expect-eigval      =", delta_e)

        print("\nmax |ExpectationValue-eigval| =", max_delta_e)


if __name__ == "__main__":
    main()
