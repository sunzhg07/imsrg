#!/usr/bin/env python3
import numpy as np
from pyIMSRG import *
import sys
#emax = int(sys.argv[1])    # maximum number of oscillator quanta in the model space
#ref = sys.argv[2]  # reference used for normal ordering
emax = 4
ref = "O16"
val = ref  # valence space

core_generator = "white"  # definition of generator eta for decoupling the core (could also use 'white')
smax_core = 0  # limit of integration in flow parameter s for first stage of decoupling

f2b = "../../input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
#f2b = "../../input/TwBME-HO_NN-only_N2LO_sat_bare_hw16_emax16_e2max32.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
f3b='../../input/NO2B_ThBME_EM7.5_1.8_2.0_IS_hw16from16_ms14_28_18.me3j.gz'
#3f3b = '../../input/NO2B_ThBME_N2LOsat_3NFJmax15_IS_hw16_ms18_36_24.stream.bin'
#f3b = "none"
mode3n = "no2b"
f3e1, f3e2, f3e3 = 14, 28, 18
LECs = "EM7.5_1820"
hw = 16


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
# Choose a beta 1-5
#H += 1.5*OperatorFromString(ms, "HCM")
#H -= 1.5*1.5*hw

#Add beta(Hcm - 3/2homega)
print("after reading files, 3-body norm is", H.ThreeBodyNorm())

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis, and retain only up to 2-body operators
HNO = hf.GetNormalOrderedH(2)

### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod(
    "magnus"
)  # Solve using the Magnus formulation. Could also be 'flow_RK4'

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetDenominatorPartitioning("MP")
imsrgsolver.SetSmax(smax_core)
imsrgsolver.Solve()


### Hs is the IMSRG-evolved Hamiltonian
Hs = imsrgsolver.GetH_s()
#Hs = imsrgsolver.Transform(HNO)
#Hs.ZeroBody += triples


IVD = OperatorFromString(ms, "IVD")
IVD = hf.TransformToHFBasis(IVD)
IVDT = imsrgsolver.Transform(IVD)
com1 = Commutator.Commutator(Hs,IVDT)
print(com1.IsAntiHermitian())
Z1 = IVDT*0.0
Z2 = IVDT*0.0 
Z1.SetHermitian()
Z2.SetHermitian()

print("Commutator after first part norm:",com1.Norm())
# IVD is only a one body operator so we can ignore the 220 part

ReferenceImplementations.comm110tt(IVDT, com1, Z1)
ReferenceImplementations.comm220tt(IVDT, com1, Z2)

# passed = True
# ut = UnitTest(ms)
# passed &= ut.Mscheme_Test_comm220tt(IVDT,com1)
# print('passed?', passed)
# exit()
denom = 8
full_com = -(Z1 *np.sqrt(3)) - Z2
print("Comm110:",Z1.ZeroBody * np.sqrt(3)/denom)
print("Comm220:",Z2.ZeroBody * np.sqrt(3)/denom)
print(full_com.ZeroBody/denom)
exit()
eom=EOM(Hs,1,1,0)



## My work for Dipole Polarizability is done below
## Obtaining E1 operator and performing similarity transofrmation
E1 = OperatorFromString(ms, "E1")
#E1.PrintOneBody()
E1 = hf.TransformToHFBasis(E1)
E1T = imsrgsolver.Transform(E1)


print("Hs:")
H_d = imsrgsolver.GetH_sDiagonal(Hs)
H_od = Hs - H_d
print("H_od is:\n\n")
Gen = imsrgsolver.generator
omega_k = E1T*0.0
omega_k.SetAntiHermitian()



Niter = 1000
diis_conv_tol = 1e-10


for it in range(Niter):
    omega_k = eom.GetVSEOM_ladder_single(omega_k,-1)
    omega_old = Operator(omega_k) * 1.0
    A = imsrgsolver.GetA(E1T, H_od, omega_old)
    omega_next = Operator(omega_k) * 0.0
    Gen.UpdateGeneral(A, H_d, omega_next)
    residual = eom.GetVSEOM_ladder_single(omega_next,-1) - eom.GetVSEOM_ladder_single(omega_old,-1)
    res_norm = residual.OneBodyNorm()
    
    sp = 0
    omega_k = omega_next * (1-sp) + sp*omega_old
    if it % 10 == 0:
        print(f"Iter {it+1}, res_norm: {res_norm:.4e}")
    if res_norm < diis_conv_tol:
        print(f"Converged at iteration {it+1}.")
        break

# Final check and comparison
final_eta = Operator(omega_next)
O = imsrgsolver.CheckWork(omega_next, Hs)

O2 = eom.GetVSEOM_ladder_single(O, 1)
O1 = eom.GetVSEOM_ladder_single(E1T, 1)
O3 = O1 + O2
print("Norm of the components", O3.Norm())

print("Eta has converged, calculating dipole polarizability")
Z1 = E1T*0.0
Z2 = E1T*0.0 
ReferenceImplementations.comm110tt(final_eta,E1T, Z1)
ReferenceImplementations.comm220tt(final_eta,E1T, Z2)


dp = 0.5 * (Z1.ZeroBody + Z2.ZeroBody)
print("Dipole Polarizability One Body Part: ",Z1.ZeroBody)
print("Dipole Polarizability Two Body Part: ",Z2.ZeroBody)
print("Dipole Polarizability:", -Z1.ZeroBody - Z2.ZeroBody)




# def etafunc(h, denom, denominator_cutoff=1e-6):
#     if abs(denom) < denominator_cutoff:
#         denom = denominator_cutoff if denom >= 0 else -denominator_cutoff
#     return h / denom


# def UpdateOmega(ms, A, H_denom, omega,
#                 denominator_cutoff=1e-6,
#                 denominator_delta=0.0,
#                 use_denominator_delta=False):

#     # One-body
#     for i in ms.all_orbits:
#         for j in ms.all_orbits:
#             a_val = A.OneBody(i, j)
#             if abs(a_val) < 1e-16:
#                 continue
#             denom = H_denom.OneBody(i, i) - H_denom.OneBody(j, j)
#             if use_denominator_delta:
#                 denom += denominator_delta
#             omega.SetOneBody(i, j, etafunc(a_val, denom, denominator_cutoff))

#     # Two-body
#     num_channels = ms.GetNumberTwoBodyChannels()
#     for ch_bra in range(num_channels):
#         for ch_ket in range(num_channels):
#             tbc_bra = ms.GetTwoBodyChannel(ch_bra)
#             tbc_ket = ms.GetTwoBodyChannel(ch_ket)
#             num_bra = tbc_bra.GetNumberKets()
#             num_ket = tbc_ket.GetNumberKets()

#             if num_bra == 0 or num_ket == 0:
#                 continue

#             for ibra in range(num_bra):
#                 bra = tbc_bra.GetKet(ibra)
#                 i = bra.p
#                 j = bra.q

#                 for iket in range(num_ket):
#                     try:
#                         a_val = A.TwoBody.GetTBMEnorm_chij(ch_bra, ch_ket, ibra, iket)
#                     except:
#                         continue
#                     if abs(a_val) < 1e-16:
#                         continue

#                     ket = tbc_ket.GetKet(iket)
#                     k = ket.p
#                     l = ket.q

#                     denom = (H_denom.OneBody(i, i) + H_denom.OneBody(j, j)
#                            - H_denom.OneBody(k, k) - H_denom.OneBody(l, l))
#                     print(denom)
#                     if use_denominator_delta:
#                         denom += denominator_delta

#                     omega.TwoBody.SetTBME_chij(ch_bra, ch_ket, ibra, iket,
#                                                etafunc(a_val, denom, denominator_cutoff))

#     print(f"omega 1b norm: {omega.OneBodyNorm():.6e}")
#     print(f"omega 2b norm: {omega.TwoBodyNorm():.6e}")
#     return omega

