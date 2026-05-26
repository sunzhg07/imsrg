#!/usr/bin/env python3
import numpy as np
from pyIMSRG import *

from lanczos import *

emax = 3  # maximum number of oscillator quanta in the model space
ref = "He4"  # reference used for normal ordering
val = "p-shell"  # valence space

core_generator = "atan"  # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = "shell-model-atan"  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 50  # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 100  # limit of s for second stage of decoupling

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
# f2b = 'input/chi2b_srg0800_eMax16_EMax16_hwHO020.me2j.gz'
# f2e1,f2e2,f2l = 16,16,16
# f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
# f3e1,f3e2,f3e3 = 12,24,12
# LECs = 'srg0800'

f2b = "/Users/wolf/work/srg_io/input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b = "/Users/wolf/work/srg_io/input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
f3e1, f3e2, f3e3 = 18, 36, 18
hw = 16
mode3n = "no2b"
LECs = "EM1820"

#### Otherwise, we use the Minnesota NN potential
# LECs = 'Minnesota'
# f3b = 'none'
# hw = 20    # harmonic oscillator basis frequency

### name of file to write resulting shell model effective interaction.
### *.snt is the exension used with KSHELL
valence_fname = "/Users/wolf/work/srg_io/output/{}_{}_{}_e{}_hw{}.snt".format(
    val, ref, LECs, emax, hw
)


##########################################################################
###  END PARAMETER SETTING. BEGIN ACTUALLY DOING STUFF ##################
##########################################################################


### Create an instance of the ModelSpace class
ms = ModelSpace(emax, ref, val)
ms.SetHbarOmega(hw)

### the ReadWrite object handles reading and writing of files
rw = ReadWrite()

rank_j, parity, rank_Tz, particle_rank = 0, 0, 0, 2
if f3b != "none":
    particle_rank = 3

### Create an instance of the Operator class, representing the Hamiltonian
H = Operator(ms, rank_j, parity, rank_Tz, particle_rank)


### Either generate the matrix elements of the Minnesota potential, or read in matrix elements from file
if LECs == "Minnesota":
    H += OperatorFromString(ms, "VMinnesota")


else:
    ### Read Two-body matrix elements
    rw.ReadBareTBME_Darmstadt(f2b, H, f2e1, f2e2, f2l)
    ### Read Three-body matrix elements
    if f3b != "none":
        if mode3n == "no2b":
            H.ThreeBody.SetMode("no2b")
            H.ThreeBody.ReadFile([f3b], [f3e1, f3e2, f3e3])
        else:
            rw.Read_Darmstadt_3body(f3b, H, f3e1, f3e2, f3e3)


### Add the relative kinetic energy, so H = Trel + V
H += OperatorFromString(ms, "Trel")

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis
HNO = hf.GetNormalOrderedH(2)
### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod(
    "magnus"
)  # Solve using the Magnus formulation. Could also be 'flow_RK4'

rw = ReadWrite()

phase = 2
if phase == 1:
    imsrgsolver.SetGenerator(core_generator)
    imsrgsolver.SetSmax(smax_core)

    ### Do the first stage of integration to decouple the core
    imsrgsolver.Solve()

    ### Now set the generator for the second stage to decouple the valence space
    imsrgsolver.SetGenerator(valence_generator)
    imsrgsolver.SetSmax(smax_valence)

    imsrgsolver.Solve()

    ### Hs is the IMSRG-evolved Hamiltonian
    Hs = imsrgsolver.GetH_s()
    Hs.ZeroBody = 0.0
    rw.WriteOperator(Hs, "Hs_back")
    rw.WriteTokyo(Hs, valence_fname, "")
else:
    Hs = Operator(ms, rank_j, parity, rank_Tz, particle_rank)
    rw.ReadOperator(Hs, "Hs_back")


# tdm_op = read_tdm("he8.ref", ms)


cm = Commutator

eom = EOM(Hs, "he8.ref", rank_j, parity, rank_Tz)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()
# eom.PrintConfigs()
nm = eom.GetVSEOM_Overlap_multiref(Hs)
print(
    "reference energy: ",
    nm,
    "vs [ -4.19234292  9.75174484 ] in nmax2 for he6_2",
)

##
##
##
##
eom.force_decouple(Hs)
# print("norm here a: ", Hs.Norm())
# eom.EraseQspace(Hs)
# print("norm here b: ", Hs.Norm())
eom.EraseValence(Hs)

unt = UnitTest(ms)
### set anti hermitian
rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1

h1 = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)
h2 = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)


chi_a = eom.GetVSEOM_ladder_multiref(h1, 1)
chi_b = eom.GetVSEOM_ladder_multiref(h2, 1)


print("test herm of norm")
nm = norm_multiref(eom, chi_a, chi_b)
print("original py ab: ", nm)
nm = norm_multiref(eom, chi_b, chi_a)
print("original py ba: ", nm)
nm = eom.NormMultiref(chi_a, chi_b)
print("nab c++: ", nm)
nm = eom.NormMultiref(chi_b, chi_a)
print("nab c++: ", nm)


print("test herm of norm c++")
nm = eom.ComputeNorm(chi_a, chi_b)
print("nab c++: ", nm)
nm = eom.ComputeNorm(chi_b, chi_a)
print("nab c++ direct: ", nm)


print("test norm3 implement")

h3ab = norm3_multiref_new(eom, chi_b, chi_b, Hs, ms)
print("norm3 ab py: ", h3ab)

nm2 = dcom222312(eom, Hs, chi_b)
print("factorized: ", nm2)

print("test final:")

hb = htc_multiref(eom, Hs, chi_b)
hab_2b_py = norm_multiref(eom, chi_a, hb)
print("nab py: ", hab_2b_py)


hb = htc_multiref(eom, Hs, chi_a)
hba_2b_py = norm_multiref(eom, chi_b, hb)
print("nba py: ", hba_2b_py)


chi_c = chi_a + chi_b

nm1 = eom.DcomMultiref(Hs, chi_a)

nm2 = eom.DcomMultiref(Hs, chi_b)

nm3 = eom.DcomMultiref(Hs, chi_c)


hab_plus_hba_3b_decomposed = nm3 - nm1 - nm2
print("dcomposed 3b ab+ba: ", hab_plus_hba_3b_decomposed)


hb = eom.HtcMultiref(Hs, chi_b)
hab_2b = eom.NormMultiref(chi_a, hb)


ha = eom.HtcMultiref(Hs, chi_a)
hba_2b = eom.NormMultiref(chi_b, ha)

hab_minus_hba_2b = hab_2b - hba_2b
hab_minus_hba_3b_decomposed = -hab_minus_hba_2b

hab_3b_decomposed = (
    hab_plus_hba_3b_decomposed + hab_minus_hba_3b_decomposed
) / 2
hba_3b_decomposed = (
    hab_plus_hba_3b_decomposed - hab_minus_hba_3b_decomposed
) / 2


print("2b hab, hba: ", hab_2b, hba_2b)
print("2bdiff hab-hba: ", hab_minus_hba_2b)
print("decomposed 3bdiff hab-hba: ", hab_minus_hba_3b_decomposed)


hab_3b_direct = eom.Norm3Multiref(chi_a, chi_b, Hs)
print("norm3 ab c++ direct/decomposed: ", hab_3b_direct, hab_3b_decomposed)
hba_3b_direct = eom.Norm3Multiref(chi_b, chi_a, Hs)
print("norm3 ba c++ direct/decomposed: ", hba_3b_direct, hba_3b_decomposed)

print("norm3 ab direct-decomposed: ", hab_3b_direct - hab_3b_decomposed)
print("norm3 ba direct-decomposed: ", hba_3b_direct - hba_3b_decomposed)
print(
    "3bdiff direct/decomposed: ",
    hab_3b_direct - hba_3b_direct,
    hab_minus_hba_3b_decomposed,
)
print("bruteforce 3b ab+ba: ", hab_3b_direct + hba_3b_direct)

np.testing.assert_allclose(hab_3b_direct, hab_3b_decomposed, rtol=1e-8, atol=1e-8)
np.testing.assert_allclose(hba_3b_direct, hba_3b_decomposed, rtol=1e-8, atol=1e-8)
np.testing.assert_allclose(
    hab_3b_direct + hba_3b_direct,
    hab_plus_hba_3b_decomposed,
    rtol=1e-8,
    atol=1e-8,
)
np.testing.assert_allclose(
    hab_3b_direct - hba_3b_direct,
    hab_minus_hba_3b_decomposed,
    rtol=1e-8,
    atol=1e-8,
)


nm3 = eom.Norm3Multiref(chi_a, chi_a, Hs)
print("norm3 aa c++: ", nm3)

nm4 = eom.DcomMultiref(Hs, chi_a)
print("norm3 decom: ", nm4)

hb = htc_multiref(eom, Hs, chi_b, eom.ProjectOprator)
hab = norm_multiref(eom, chi_a, hb)
print("hab py: ", hab)

hb = eom.HtcMultiref(Hs, chi_b)
hab = eom.NormMultiref(chi_a, hb)
print("hab c++: ", hab)
