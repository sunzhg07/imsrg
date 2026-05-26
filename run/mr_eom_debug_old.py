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
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()

rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1
unt = UnitTest(ms2)

fixed_start_op = os.environ.get("FIXED_START_OP", "output/arnoldi_start.op")

h1 = Operator(ms2, rank_j, parity, rank_Tz, particle_rank)
rw.ReadOperator(h1, fixed_start_op)
print("loaded fixed start operator", fixed_start_op)

chi_b = eom.GetVSEOM_ladder_multiref(h1, 1)

max_iter = 20
state_want = 5

result = eom.ArnoldiSolve_old(chi_b, max_iter, state_want)

hall = np.array(result.hall, dtype=float)
energies = np.array(result.energies, dtype=float)

print("solver = old")
print("hall shape =", hall.shape)
print("energies =", energies)

np.save("output/arnoldi_hall_old.npy", hall)
np.save("output/arnoldi_energies_old.npy", energies)
print("saved output/arnoldi_hall_old.npy")
print("saved output/arnoldi_energies_old.npy")