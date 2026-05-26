#!/usr/bin/env python3
import os

import numpy as np
from pyIMSRG import *

from lanczos import *

emax = 3  # maximum number of oscillator quanta in the model space
ref = "He8"  # reference used for normal ordering
val = "p-shell"  # valence space

core_generator = "atan"  # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = "shell-model-atan"  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 50  # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 100  # limit of s for second stage of decoupling

f2b = "input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b = "input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
f3e1, f3e2, f3e3 = 18, 36, 18
hw = 16
mode3n = "no2b"
LECs = "EM1820"

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

eom = EOM(Hs, "he8.ref", 0, 0, 0)
use_projection = os.environ.get("USE_PROJECTION", "1") != "0"
eom.SetArnoldiUseProjection(use_projection)
print("use arnoldi projection =", use_projection)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()

rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1
unt = UnitTest(ms2)

fixed_start_op = os.environ.get("FIXED_START_OP", "output/arnoldi_start.op")
fixed_start_vector_op = os.environ.get("FIXED_START_VECTOR_OP", "output/arnoldi_start_vector.op")
start_vector_op = os.environ.get("START_VECTOR_OP")

if start_vector_op:
    chi_b = Operator(ms2, rank_j, parity, rank_Tz, particle_rank)
    rw.ReadOperator(chi_b, start_vector_op)
    print("loaded Arnoldi start vector", start_vector_op)
else:
    h1 = unt.RandomOp(ms2, rank_j, rank_Tz, parity, particle_rank, herm)
    rw.WriteOperator(h1, fixed_start_op)
    print("saved raw fixed start operator", fixed_start_op)
    chi_b = eom.GetVSEOM_ladder_multiref(h1, 1)
    rw.WriteOperator(chi_b, fixed_start_vector_op)
    print("saved laddered fixed start vector", fixed_start_vector_op)

max_iter = int(os.environ.get("MAX_ITER", "20"))
state_want = int(os.environ.get("STATE_WANT", "5"))

result = eom.ArnoldiSolve(chi_b, max_iter, state_want)

hall = np.array(result.hall, dtype=float)
energies = np.array(result.energies, dtype=float)

print("solver = new")
print("hall shape =", hall.shape)
print("energies =", energies)

np.save("output/arnoldi_hall_new.npy", hall)
np.save("output/arnoldi_energies_new.npy", energies)
print("saved output/arnoldi_hall_new.npy")
print("saved output/arnoldi_energies_new.npy")