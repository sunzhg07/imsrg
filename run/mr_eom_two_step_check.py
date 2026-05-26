#!/usr/bin/env python3
import math
import numpy as np
from pyIMSRG import *


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
LECs = "EM1820"

def build_hamiltonian():
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


def build_eom(Hs):
    eom = EOM(Hs, "he8.ref", 0, 0, 0)
    eom.ConstructConfigs()
    eom.ConstructNormMatrix()
    eom.ConstructProjectMatrix()
    return eom


def generate_start_vector(eom, ms2):
    rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1
    unt = UnitTest(ms2)
    h1 = unt.RandomOp(ms2, rank_j, rank_Tz, parity, particle_rank, herm)
    chi_b = eom.GetVSEOM_ladder_multiref(h1, 1)
    return chi_b


def normalize(eom, vec):
    nn = eom.ComputeNorm(vec, vec)
    return vec / math.sqrt(nn)


def build_next_basis_vector(eom, Hs, basis_vec, h1v_basis):
    w = h1v_basis * 1.0
    eom.ProjectOprator(w)
    for _ in range(2):
        cij = eom.ComputeNorm(basis_vec, w)
        w = w - cij * basis_vec
        eom.ProjectOprator(w)
    bj = eom.ComputeNorm(w, w)
    return bj, w / math.sqrt(bj)


def print_start_vector_checks(eom, Hs, v0):
    h1v0 = eom.HtcMultiref(Hs, v0)
    h1_diag = eom.ComputeNorm(v0, h1v0)
    h2_diag_direct = eom.DcomMultiref(Hs, v0)
    h2_diag_cross = eom.DcomMultiref(Hs, v0 + v0)
    h2_diag_polar = 0.5 * (h2_diag_cross - 2.0 * h2_diag_direct)

    print("\n=== Check 1: Start Vector Expectations ===")
    print("norm(v0)                =", eom.ComputeNorm(v0, v0))
    print("<v0|H1|v0>              =", h1_diag)
    print("<v0|H2|v0> direct       =", h2_diag_direct)
    print("<v0|H2|v0> polarization =", h2_diag_polar)
    print("H2 diagonal delta       =", h2_diag_polar - h2_diag_direct)
    print("<v0|H|v0> direct        =", h1_diag + h2_diag_direct)
    print("<v0|H|v0> polarization  =", h1_diag + h2_diag_polar)
    return h1v0, h1_diag, h2_diag_direct, h2_diag_polar


def print_manual_two_step_hall(eom, Hs, v0, h1v0, h1_00, h2_00_direct, h2_00_polar):
    bj, v1 = build_next_basis_vector(eom, Hs, v0, h1v0)
    h1v1 = eom.HtcMultiref(Hs, v1)

    h1_11 = eom.ComputeNorm(v1, h1v1)
    h2_11_direct = eom.DcomMultiref(Hs, v1)
    h2_11_cross = eom.DcomMultiref(Hs, v1 + v1)
    h2_11_polar = 0.5 * (h2_11_cross - 2.0 * h2_11_direct)

    h1_01 = eom.ComputeNorm(v0, h1v1)
    h1_10 = eom.ComputeNorm(v1, h1v0)
    h1_sym = 0.5 * (h1_01 + h1_10)

    h2_01_cross = eom.DcomMultiref(Hs, v0 + v1)
    h2_01_sym = 0.5 * (h2_01_cross - h2_00_direct - h2_11_direct)

    hall_manual = np.array(
        [
            [h1_00 + h2_00_polar, h1_sym + h2_01_sym],
            [h1_sym + h2_01_sym, h1_11 + h2_11_polar],
        ],
        dtype=float,
    )

    print("\n=== Check 2: Manual Two-Step Hall Build ===")
    print("bj(step 0)              =", bj)
    print("<v1|H2|v1> direct       =", h2_11_direct)
    print("<v1|H2|v1> polarization =", h2_11_polar)
    print("v1 H2 diagonal delta    =", h2_11_polar - h2_11_direct)
    print("<v0|H1|v1>              =", h1_01)
    print("<v1|H1|v0>              =", h1_10)
    print("H1 sym(0,1)             =", h1_sym)
    print("<v0+v1|H2|v0+v1>        =", h2_01_cross)
    print("H2 sym(0,1)             =", h2_01_sym)
    print("manual hall(0,1)        =", hall_manual[0, 1])
    print("manual 2x2 hall =")
    print(hall_manual)
    return hall_manual


def run_two_step_solver_checks(Hs, ms2):
    eom_new = build_eom(Hs)
    chi_new = generate_start_vector(eom_new, ms2)
    ar_new = eom_new.ArnoldiSolve(chi_new, 3, 2)

    eom_old = build_eom(Hs)
    chi_old = generate_start_vector(eom_old, ms2)
    ar_old = eom_old.ArnoldiSolve_old(chi_old, 3, 2)

    hall_new = np.array(ar_new.hall, dtype=float)
    hall_old = np.array(ar_old.hall, dtype=float)

    print("\n=== Check 3: Solver Hall Comparison For Two Steps ===")
    print("new hall =")
    print(hall_new)
    print("old hall =")
    print(hall_old)
    print("difference =")
    print(hall_new - hall_old)

    if hall_new.shape == hall_old.shape and hall_new.size:
        print("max |new-old| =", np.max(np.abs(hall_new - hall_old)))
        if hall_new.shape[0] > 1:
            print("offdiag new hall(0,1) =", hall_new[0, 1])
            print("offdiag old hall(0,1) =", hall_old[0, 1])


def main():
    Hs, ms2 = build_hamiltonian()

    eom = build_eom(Hs)
    chi_b = generate_start_vector(eom, ms2)
    v0 = normalize(eom, chi_b)

    h1v0, h1_00, h2_00_direct, h2_00_polar = print_start_vector_checks(eom, Hs, v0)
    print_manual_two_step_hall(eom, Hs, v0, h1v0, h1_00, h2_00_direct, h2_00_polar)
    run_two_step_solver_checks(Hs, ms2)


if __name__ == "__main__":
    main()