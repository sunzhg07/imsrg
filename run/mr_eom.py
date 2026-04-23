#!/usr/bin/env python3
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

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
# f2b = 'input/chi2b_srg0800_eMax16_EMax16_hwHO020.me2j.gz'
# f2e1,f2e2,f2l = 16,16,16
# f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
# f3e1,f3e2,f3e3 = 12,24,12
# LECs = 'srg0800'

f2b = "input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b = "input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"
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
valence_fname = "output/{}_{}_{}_e{}_hw{}.snt".format(val, ref, LECs, emax, hw)
hbar_fname = "output/hbar_{}_{}_{}_e{}_hw{}.dat".format(val, ref, LECs, emax, hw)
fci_fname = "output/fci_{}_{}_{}_e{}_hw{}.snt".format(val, ref, LECs, emax, hw)


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

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

stage = 2
if stage == 1:
    ### Do the first stage of integration to decouple the core
    imsrgsolver.Solve()

    ### Now set the generator for the second stage to decouple the valence space
    imsrgsolver.SetGenerator(valence_generator)
    imsrgsolver.SetSmax(smax_valence)

    imsrgsolver.Solve()

    # Hs is the IMSRG-evolved Hamiltonian
    # rw.WriteTokyo( HNO, valence_fname,'')
    Hs = imsrgsolver.GetH_s()
    Hs.UndoNormalOrdering()
    ms2 = ModelSpace(emax, "He4", val)
    Hs.SetModelSpace(ms2)
    Hs.DoNormalOrdering()

    rw.WriteOperator(Hs, hbar_fname)
else:
    ms2 = ModelSpace(emax, "He4", val)
    Hs = Operator(ms2, rank_j, parity, rank_Tz, particle_rank)
    Hs = 0.0 * Hs
    rw.ReadOperator(Hs, hbar_fname)


## ---------------------------------------------------------------
eom = EOM(Hs, "he8.ref", 0, 0, 0)
eom.ConstructConfigs()
eom.PrintConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()

##########################################################################
###  Loop over configs from cfs file, run comm223ss for each
##########################################################################
configs = []
with open("cfs") as f:
    for line in f:
        parts = line.split()
        if len(parts) < 6:
            continue
        typ, idx, c0, c1, c2, c3 = (
            parts[0],
            int(parts[1]),
            int(parts[2]),
            int(parts[3]),
            int(parts[4]),
            int(parts[5]),
        )
        configs.append((typ, idx, c0, c1, c2, c3))

# for typ, idx, c0, c1, c2, c3 in configs[0]:

# typ, idx, c0, c1, c2, c3 = configs[0]
typ, idx, c0, c1, c2, c3 = configs[1]
print(c0, c1, c2)
Z_test = Operator(ms, rank_j, rank_Tz, parity, 2)
Z_test.SetHermitian()
Z_test.ThreeBody.SetMode("pn")
X_test = 0.0 * Hs
X_test.SetHermitian()
typ, idx, c0, c1, c2, c3 = configs[2]
X_test.TwoBody.SetTBME_chij(c2, c2, c0, c1, 1.0)
print(X_test.Norm())
Y_test = 0.0 * Hs
Y_test.SetAntiHermitian()
typ, idx, c0, c1, c2, c3 = configs[3]
Y_test.TwoBody.SetTBME_chij(c2, c2, c0, c1, 1.0)

Commutator.comm223ss(Y_test, X_test, Z_test)
norm_z = Z_test.Norm()
print(norm_z)
# val_3b = Z_test.ThreeBody.GetME_pn(0, 2, 3, 2, 2, 2, 2, 2, 2)
# val_3b = Z_test.ThreeBody.GetME_pn(0, 2, 3, 2, 2, 2, 2, 2, 2)
# print(f"Z_test 3b ME (a=2,b=2,Jab=0,c=2,d=2,e=2,Jde=2,f=2,2J=3) = {val_3b:.8f}")


# Print all non-zero 3b matrix elements in Z_test
for ch_bra, ch_ket in Z_test.ThreeBody.Get_ch_start_keys():
    tbc_bra = ms.GetThreeBodyChannel(ch_bra)
    tbc_ket = ms.GetThreeBodyChannel(ch_ket)
    nk_bra = tbc_bra.GetNumber3bKets()
    nk_ket = tbc_ket.GetNumber3bKets()
    for ibra in range(nk_bra):
        kb = tbc_bra.GetKet(ibra)
        iket_start = ibra if ch_bra == ch_ket else 0
        for iket in range(iket_start, nk_ket):
            v = Z_test.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket)
            if abs(v) > 1e-10:
                kk = tbc_ket.GetKet(iket)
                print(
                    f"  3b nonzero: ch_bra={ch_bra}(2J={tbc_bra.twoJ}) ch_ket={ch_ket}(2J={tbc_ket.twoJ})"
                    f" ibra={ibra}({kb.p},{kb.q},{kb.r},Jpq={kb.Jpq}) iket={iket}({kk.p},{kk.q},{kk.r},Jpq={kk.Jpq}) val={v:.8f}"
                )

# Check diagonal (Jpq=0,Jpq=0) and cross (Jpq=0,Jpq=2) elements explicitly
for ch_bra, ch_ket in Z_test.ThreeBody.Get_ch_start_keys():
    tbc_bra = ms.GetThreeBodyChannel(ch_bra)
    tbc_ket = ms.GetThreeBodyChannel(ch_ket)
    if ch_bra != ch_ket:
        continue
    for ibra in range(tbc_bra.GetNumber3bKets()):
        kb = tbc_bra.GetKet(ibra)
        for iket in range(ibra, tbc_ket.GetNumber3bKets()):
            kk = tbc_ket.GetKet(iket)
            if not (kb.p == kk.p == 2 and kb.q == kk.q == 2 and kb.r == kk.r == 2):
                continue
            v = Z_test.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket)
            print(
                f"  diag check: ibra={ibra}(Jpq={kb.Jpq}) iket={iket}(Jpq={kk.Jpq}) val={v:.8f}"
            )

# old signature: (a, b, c, d, e, f, g, j0, j2)
# val = eom.ThreeBody_Diagram(a, b, c, d, e, f, g, j0, j2)

val = eom.ThreeBody_Diagram(2, 2, 4, 4, 4, 4, 18, 0, 0)
# ThreeBME_type GetME_pn(int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f)
val_pn = Z_test.ThreeBody.GetME_pn(0, 1, 3, 4, 4, 2, 2, 4, 4)
print(f"Z_test(4,4,2,2,4,4, jpq_bra=1, jpq_ket=1, jtot=3) = {val_pn:.8f}")
val_pn = Z_test.ThreeBody.GetME_pn(0, 2, 3, 4, 4, 2, 2, 4, 4)
print(f"Z_test(4,4,2,2,4,4, jpq_bra=1, jpq_ket=1, jtot=3) = {val_pn:.8f}")
# val = eom.ThreeBody_Diagram(2,2,2,2,2,2,16, 0, 0)
# Compute ThreeBody_Diagram for all (i,j) pairs of ppvv configs
# For each config: c[0]=ibra -> (a=ket.p, b=ket.q), c[1]=iket -> (c=ket.p, d=ket.q), c[2]=ch
# ms2 = Hs.GetModelSpace()
# for typ_j, idx_j, c0_j, c1_j, c2_j, c3_j in configs:
#     if typ_j != "ppvv" or typ != "ppvv":
#         continue
#     tbc_bra = ms2.GetTwoBodyChannel(c2)
#     tbc_ket = ms2.GetTwoBodyChannel(c2_j)
#     ket_bra_pp = tbc_bra.GetKet(c0)  # ibra of bra config
#     ket_bra_vv = tbc_bra.GetKet(c1)  # iket of bra config
#     ket_ket_pp = tbc_ket.GetKet(c0_j)  # ibra of ket config
#     ket_ket_vv = tbc_ket.GetKet(c1_j)  # iket of ket config
#     a1, b1 = ket_bra_pp.p, ket_bra_pp.q
#     c_orb1, d1 = ket_bra_vv.p, ket_bra_vv.q
#     a2, b2 = ket_ket_pp.p, ket_ket_pp.q
#     c_orb2, d2 = ket_ket_vv.p, ket_ket_vv.q
#     # only when b1==b2 and both valence (matches C++ condition)
#     ob1 = ms2.GetOrbit(b1)
#     ob2 = ms2.GetOrbit(b2)
#     if b1 == b2 and ob1.cvq == 1 and ob2.cvq == 1:
#         J_bra = tbc_bra.J
#         J_ket = tbc_ket.J
#         val = eom.ThreeBody_Diagram(c_orb1, d1, a2, a1, d2, c_orb2, b1, J_bra, J_ket)
#         print(f"ThreeBody_Diagram ({typ} idx={idx}, {typ_j} idx={idx_j}): {val:.8f}")

# print(f"{typ} idx={idx} c=({c0},{c1},{c2},{c3}) -> Z.Norm={norm_z:.8f}")
