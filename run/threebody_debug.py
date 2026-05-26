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
imsrgsolver.SetHunterGatherer(True)
imsrgsolver.SetOmegaNormMax(0.1)
imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

stage = 1
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
    Omega = imsrgsolver.GetOmega(0)
    Htilde = Hs + 0.5 * Commutator.Commutator(Omega, Hs)
    Hs.SetParticleRank(3)
    Hs.ThreeBody.SetMode("pn")
    Commutator.comm223ss(Omega, Htilde, Hs)
    valence3n_filename = "threebd_forbench.snt"
    valence2n_filename = "twobd_forbench.snt"
    rw.WriteValence3body(Hs.ThreeBody, valence3n_filename)

    Hss = Hs.ReNormalOrderCore()
    print("computing the norm of Hs 3body: ", Hs.ThreeBodyNorm())
    Hs.OneBody = Hss.OneBody
    Hs.TwoBody = Hss.TwoBody
    Hs.ThreeBody = Hs.ThreeBody 

    rw.WriteTokyo(Hs, valence2n_filename, "")
    # rw.WriteOperator(Hs, hbar_fname)
else:
    ms2 = ModelSpace(emax, "He4", val)
    Hs = Operator(ms2, rank_j, parity, rank_Tz, particle_rank)
    Hs = 0.0 * Hs
    rw.ReadOperator(Hs, hbar_fname)


## ---------------------------------------------------------------
eom = EOM(Hs, "Be10.ref", 0, 0, 0)
eom.ConstructConfigs()
# eom.PrintConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()
#print(eom.GetVSEOM_Overlap_multiref(Hs))

unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
h4= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)

chia= eom.GetVSEOM_ladder_single(h3,0)
chib= eom.GetVSEOM_ladder_single(h4,0)
nm=norm_multiref(eom,chia, chib)
print('norm of ladder single: ', nm)


# now we benchmark the new nrom function


###########################################################################
####  Loop over configs from cfs file, run comm223ss for each
###########################################################################
# configs = []
# with open("cfs") as f:
#    for line in f:
#        parts = line.split()
#        if len(parts) < 6:
#            continue
#        typ, idx, c0, c1, c2, c3 = (
#            parts[0],
#            int(parts[1]),
#            int(parts[2]),
#            int(parts[3]),
#            int(parts[4]),
#            int(parts[5]),
#        )
#        configs.append((typ, idx, c0, c1, c2, c3))
#
# configs2 = []
# with open("cfs2") as f:
#    for line in f:
#        parts = line.split()
#        if len(parts) < 7:
#            continue
#        typ, idx, c0, c1, c2, c3, c4 = (
#            parts[0],
#            int(parts[1]),
#            int(parts[2]),
#            int(parts[3]),
#            int(parts[4]),
#            int(parts[5]),
#            int(parts[6]),
#        )
#        configs2.append((typ, idx, c0, c1, c2, c3, c4))
#
#
# Z_test = Operator(ms, rank_j, rank_Tz, parity, 3)
# Z_test.SetHermitian()
# Z_test.ThreeBody.SetMode("pn")
#
# X_test = 0.0 * Hs
# X_test.SetHermitian()
#
# Y_test = 0.0 * Hs
# Y_test.SetAntiHermitian()
#
## for i, cfsi in enumerate(configs):
##    for j, cfsj in enumerate(configs):
##        if i>j:
##            continue
##        stx, idx, d,g1,a,b, j0=configs2[i]
##        stx, idx, c,g2,f,e, j2=configs2[j]
##        if(g1!=g2):
##            continue
##        print(i,j,'exist')
##        typ_i, idx_i, c0_i, c1_i, c2_i, c3_i = cfsi
##        typ_j, idx_j, c0_j, c1_j, c2_j, c3_j = cfsj
##
##        X_test=0.0*X_test
##        Y_test=0.0*Y_test
##        Z_test=0.0*Z_test
##
##        X_test.TwoBody.SetTBME_chij(c2_i, c2_i, c0_i, c1_i, 1.0)
##        Y_test.TwoBody.SetTBME_chij(c2_j, c2_j, c0_j, c1_j, 1.0)
##
##        Commutator.comm223ss(X_test, Y_test, Z_test)
##
##        result = eom.ThreeBody_Diagram_Entries(a,b,c,d,e,f,g1,j0,j2)
##
##        for jab, jde, jtot, diag_val in result:
##            z_val = Z_test.ThreeBody.GetME_pn(jab, jde, jtot, a, b, c, d, e, f)
##            ratio = diag_val / z_val if abs(z_val) > 1e-14 else float('nan')
##            print(f"  {a} {b} {c} {d} {e} {f} {g1}  j0={jab} j1={jde} jtot={jtot}  diagram={diag_val:.8f}  Z_test={z_val:.8f}  ratio={ratio:.6f}")
#
#
# if 1 == 1:
#    i = 0
#    j = 6
#    cfsi = configs[i]
#    cfsj = configs[j]
#    stx, idx, d, g1, a, b, j0 = configs2[i]
#    stx, idx, c, g2, f, e, j2 = configs2[j]
#    print(i, j, "exist")
#    typ_i, idx_i, c0_i, c1_i, c2_i, c3_i = cfsi
#    typ_j, idx_j, c0_j, c1_j, c2_j, c3_j = cfsj
#    X_test = 0.0 * X_test
#    Y_test = 0.0 * Y_test
#    Z_test = 0.0 * Z_test
#    X_test.TwoBody.SetTBME_chij(c2_i, c2_i, c0_i, c1_i, 1.0)
#    Y_test.TwoBody.SetTBME_chij(c2_j, c2_j, c0_j, c1_j, 1.0)
#    Commutator.comm223ss(X_test, Y_test, Z_test)
#    print(f"  Z_test.ThreeBodyNorm() = {Z_test.ThreeBodyNorm():.12f}")
#
#    result = eom.ThreeBody_Diagram_Entries(a, b, c, d, e, f, g1, j0, j2)
#
#    for jab, jde, jtot, diag_val in result:
#        z_val = Z_test.ThreeBody.GetME_pn(jab, jde, jtot, a, b, c, d, e, f)
#        ratio = diag_val / z_val if abs(z_val) > 1e-14 else float("nan")
#        print(
#            f"  {a} {b} {c} {d} {e} {f} {g1}  j0={jab} j1={jde} jtot={jtot}  diagram={diag_val:.8f}  Z_test={z_val:.8f}  ratio={ratio:.6f}"
#        )
