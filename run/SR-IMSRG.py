#!/usr/bin/env python3
import numpy as np
from pyIMSRG import *



emax = 3       # maximum number of oscillator quanta in the model space
ref = 'O22'     # reference used for normal ordering
val = ref # valence space

core_generator = 'atan'   # definition of generator eta for decoupling the core (could also use 'white')
smax_core = 50       # limit of integration in flow parameter s for first stage of decoupling

f2b='../../input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz'
f2e1,f2e2,f2l = 14,28,14
f3b='../../input/NO2B_ThBME_EM7.5_1.8_2.0_IS_hw16from16_ms14_28_18.me3j.gz'
f3e1,f3e2,f3e3 = 14,28,18
mode3n='no2b'
LECs = 'EM7.5_1820'
hw=16




### Create an instance of the ModelSpace class
ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

### the ReadWrite object handles reading and writing of files
rw = ReadWrite()

rank_j, parity, rank_Tz, particle_rank = 0,0,0,2
if f3b != 'none':
   particle_rank = 3

### Create an instance of the Operator class, representing the Hamiltonian
H = Operator(ms,rank_j, parity, rank_Tz, particle_rank)

### Either generate the matrix elements of the Minnesota potential, or read in matrix elements from file
if LECs == 'Minnesota':
    H += OperatorFromString(ms,'VMinnesota')

else:
  ### Read Two-body matrix elements
  rw.ReadBareTBME_Darmstadt(f2b,H,f2e1,f2e2,f2l)
  ### Read Three-body matrix elements
  if f3b != 'none':
     if mode3n == 'no2b':
        H.ThreeBody.SetMode('no2b')
        H.ThreeBody.ReadFile([f3b],[f3e1,f3e2,f3e3])
     else:
        rw.Read_Darmstadt_3body(f3b,H,f3e1,f3e2,f3e3)


### Add the relative kinetic energy, so H = Trel + V
H += OperatorFromString(ms,'Trel')
print('after reading files, 3-body norm is',H.ThreeBodyNorm())

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis, and retain only up to 2-body operators
HNO = hf.GetNormalOrderedH(2)

### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod('magnus')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
imsrgsolver.Solve()


### Hs is the IMSRG-evolved Hamiltonian
Hs = imsrgsolver.GetH_s()


## My work for Dipole Polarizability is done below
## Obtaining E1 operator and performing similarity transofrmation
E1 = OperatorFromString(ms,'E1')
E1T = imsrgsolver.Transform(E1)
E1T.PrintOneBody()
print("Before confusion")
H_d = imsrgsolver.GetH_sDiagonal(Hs)

H_od = Hs - H_d
#Hs.PrintOneBody()
#H_od.PrintOneBody()
#H_d.PrintOneBody()
#print(H_d.PrintTwoBody())
Gen = imsrgsolver.generator 
A = imsrgsolver.GetA(E1T, H_od, H_d) 
#A.PrintOneBody()

E1Tin = E1T
print("\n")
Eta_in = Gen.GetEta()
Gen.UpdateGeneral(A, H_d, E1Tin)
Eta_1up = Gen.GetEta()
Eta_1up.PrintOneBody()
print("\n")
A2 = imsrgsolver.GetA(E1T, H_od, Eta_1up)
Gen.UpdateGeneral(A2,H_d,E1T)
Eta_2up = Gen.GetEta()
Eta_2up.PrintOneBody()

print("\n")
print("\n")
E1T.PrintOneBody()
Hs.PrintOneBody()
