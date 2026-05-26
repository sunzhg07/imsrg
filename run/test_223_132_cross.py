#!/usr/bin/env python3
import os
import time

from pyIMSRG import *


emax = int(os.environ.get("EMAX", "3"))
ref = "He4"
val = "p-shell"
hw = 16

rank_j, parity, rank_Tz, particle_rank = 0, 0, 0, 2
phase = int(os.environ.get("PHASE", "2"))
nrepeat = int(os.environ.get("NREPEAT", "1"))
herm = int(os.environ.get("HERM", "1"))

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
f3e1, f3e2, f3e3 = 18, 36, 18
mode3n = "no2b"

core_generator = "atan"
valence_generator = "shell-model-atan"
smax_core = 50
smax_valence = 100


def build_or_read_hamiltonian(ms):
    rw = ReadWrite()
    if phase == 1:
        H = Operator(ms, rank_j, parity, rank_Tz, 3)
        rw.ReadBareTBME_Darmstadt(f2b, H, f2e1, f2e2, f2l)
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
        Hs.ZeroBody = 0.0
        rw.WriteOperator(Hs, "Hs_back")
        return Hs

    Hs = Operator(ms, rank_j, parity, rank_Tz, particle_rank)
    rw.ReadOperator(Hs, "Hs_back")
    return Hs


def zero_like(op):
    z = op * 0.0
    z.SetHermitian()
    return z


def operator_norms(op):
    return {
        "Norm": op.Norm(),
        "OneBodyNorm": op.OneBodyNorm(),
        "TwoBodyNorm": op.TwoBodyNorm(),
        "ThreeBodyNorm": op.ThreeBodyNorm(),
    }


def print_norms(label, op):
    norms = operator_norms(op)
    print(
        f"{label:>14}: "
        f"Norm={norms['Norm']:.12e} "
        f"1b={norms['OneBodyNorm']:.12e} "
        f"2b={norms['TwoBodyNorm']:.12e} "
        f"3b={norms['ThreeBodyNorm']:.12e}"
    )


ms = ModelSpace(emax, ref, val)
ms.SetHbarOmega(hw)
Hs = build_or_read_hamiltonian(ms)

eom = EOM(Hs, "he8.ref", rank_j, parity, rank_Tz)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()
eom.force_decouple(Hs)
eom.EraseValence(Hs)

unt = UnitTest(ms)
h1 = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)
chi_in = eom.GetVSEOM_ladder_multiref(h1, 1)
eta = eom.GetVSEOM_ladder_multiref(chi_in, -1)

print("benchmarking comm223_132_cross only")
print(f"EMAX={emax} PHASE={phase} NREPEAT={nrepeat} HERM={herm}")
print_norms("eta", eta)
print_norms("Hs", Hs)

for irepeat in range(nrepeat):
    z_reference = zero_like(chi_in)
    t0 = time.perf_counter()
    ReferenceImplementations.comm223_132_cross(eta, Hs, z_reference)
    t_reference = time.perf_counter() - t0

    z_factorized = zero_like(chi_in)
    t0 = time.perf_counter()
    Commutator.FactorizedDoubleCommutator.comm223_132_cross(eta, Hs, z_factorized)
    t_factorized = time.perf_counter() - t0

    zdiff = z_factorized - z_reference
    rel = zdiff.Norm() / max(z_reference.Norm(), 1.0e-30)

    print(f"\nrepeat {irepeat + 1}/{nrepeat}")
    print(f"reference time:  {t_reference:.6f} s")
    print(f"factorized time: {t_factorized:.6f} s")
    print_norms("reference", z_reference)
    print_norms("factorized", z_factorized)
    print_norms("difference", zdiff)
    print(f"relative difference Norm: {rel:.12e}")