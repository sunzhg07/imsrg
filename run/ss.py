#!/usr/bin/env python3
from pyIMSRG import *
import numpy as np
from lanczos import *



emax =2         # maximum number of oscillator quanta in the model space
ref = 'He4'     # reference used for normal ordering
val = 'p-shell' # valence space

core_generator = 'atan'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model-atan'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 50       # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 100   # limit of s for second stage of decoupling

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
#f2b = 'input/chi2b_srg0800_eMax16_EMax16_hwHO020.me2j.gz'
#f2e1,f2e2,f2l = 16,16,16
#f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
#f3e1,f3e2,f3e3 = 12,24,12
#LECs = 'srg0800'

f2b='input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz'
f2e1,f2e2,f2l = 14,28,14
f3b='input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin'
f3e1,f3e2,f3e3 = 18,36,18
hw=16
mode3n = 'no2b'
LECs = 'EM1820'

#### Otherwise, we use the Minnesota NN potential
#LECs = 'Minnesota'
#f3b = 'none'
#hw = 20    # harmonic oscillator basis frequency

### name of file to write resulting shell model effective interaction.
### *.snt is the exension used with KSHELL
valence_fname = 'output/{}_{}_{}_e{}_hw{}.snt'.format(val,ref,LECs,emax,hw)


##########################################################################
###  END PARAMETER SETTING. BEGIN ACTUALLY DOING STUFF ##################
##########################################################################


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

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis
HNO = hf.GetNormalOrderedH(2)
ips1 = ms.GetOrbitIndex(0,0,1,-1)
ins1 = ms.GetOrbitIndex(0,0,1,+1)
HNO.SetOneBody(ips1,ips1, HNO.GetOneBody(ips1,ips1)-10)
HNO.SetOneBody(ins1,ins1, HNO.GetOneBody(ins1,ins1)-10)
### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod('magnus')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
imsrgsolver.Solve()

### Now set the generator for the second stage to decouple the valence space
imsrgsolver.SetGenerator(valence_generator)
imsrgsolver.SetSmax(smax_valence)

imsrgsolver.Solve()

### Hs is the IMSRG-evolved Hamiltonian
#rw.WriteTokyo( HNO, valence_fname,'')
Hs = imsrgsolver.GetH_s()
Hs.ZeroBody=0.



tdm_op=read_tdm("he6.ref",ms)


#cm.SetIMSRG3Noqqq(True)
#cm.SetIMSRG3Onlyvvv(True)
gm=Generator()
cm=Commutator


#rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
#
#
eom=EOM(Hs, tdm_op,rank_j, parity, rank_Tz)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()
#nm=eom.GetVSEOM_Overlap_multiref(Hs)
#print('reference energy: ', nm, 'vs [ -4.19234292  9.75174484 ] in nmax2 for he6_2',)
##
##
##
##
unt = UnitTest(ms)
### set anti hermitian
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
h1= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
h2= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
t3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, 3,herm)
t3.ThreeBody.SetMode("pn")

chi_a=eom.GetVSEOM_ladder_multiref(h1,1)
print(chi_a.IsAntiHermitian())
opo=chi_a*0.
opo.SetHermitian()
chi_a.ThreeBody.Erase()

cm.comm223ss(chi_a, Hs, t3)
cm.comm231ss(chi_a, t3, opo )
print('direct norm: ', opo.Norm())
cm.comm232ss(chi_a, t3, opo )
print('direct norm: ', opo.Norm())

ops=chi_a*0.
ops.SetHermitian()

cm.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True);
cm.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True);
cm.FactorizedDoubleCommutator.comm223_231(chi_a, Hs, ops);
print('fact norm: ', ops.Norm())
cm.FactorizedDoubleCommutator.comm223_232(chi_a, Hs, ops);
print('fact norm: ', ops.Norm())

#chi_b=eom.GetVSEOM_ladder_multiref(h2,0)

#print(chi_b.Norm())

#hb=htc_multiref(eom,Hs, chi_b)
#hab = norm_multiref(eom,chi_a, hb)
#print('<a|H|b>= ',hab)

#ha=htc_multiref(eom,Hs, chi_a)
#hba = norm_multiref(eom,chi_b, ha)
#print('<b|H|a>= ',hba)
#
#
#
#nma, chi3aa = dcom222312(eom, Hs, chi_a)
#print('<a|h3|a> facto= ',nma)

#nmb, chi3bb = dcom222312(eom, Hs, chi_b)
#print('<b|h3|b>= ',nmb, chi3bb.Norm())
#
#chiab=chi_a+chi_b
#nmab, chi3apb = dcom222312(eom, Hs, chiab)
#print('<a+b|h3|a+b>= ',nmab, chi3apb.Norm())
#
#nmab=nmab-nma-nmb
#
#h3ab=(hba-hab+nmab)/2.
#h3ba=-(hba-hab-nmab)/2.
#
#print('<a|h3|b> in direct= ', h3ab/2)
#print('<b|h3|a> in direct= ', h3ba/2)
#unt.TestFactorizedDoubleCommutators()

#nm3=norm3_multiref_fact(eom, chi_a,chi_a, Hs, ms)
#print('<a|h3|a> direct= ',nm3)
#nmaa, chi3apb = dcom222312(eom, Hs, chi_a)
#print('<a|h3|a> factor= ',nmaa)

#
#nm3=norm3_multiref(eom, chi_b,chi_b, Hs, ms)
#print('<b|h3|b> direct= ',nm3)

#nm3=norm3_multiref(eom, chi_a,chi_b, Hs, ms)
#print('<a|h3|b> direct= ',nm3)
#
#nm3=norm3_multiref(eom, chi_b,chi_a, Hs, ms)
#print('<b|h3|a> direct= ',nm3)





#e,v1,v2=lanczos_proc(htc_multiref, norm_multiref, Hs, chi_a, 140, 8,ms,eom,  prjop=eom.ProjectOprator)

#e,v1,v2=arnoldi_proc_new(htc_multiref, norm_multiref, norm3_multiref, Hs, chi, 140, 4,ms,eom, prjop=eom.ProjectOprator)


jrank = 0;
tz = 0;
parity = 0;
particle_rank = 2;

eta = chi_a
H = Hs
t3=Operator(ms, jrank, tz, parity, 3);
opo=chi_a*0.
ops=chi_a*0.
opo.SetHermitian()
ops.SetHermitian()

#opo=Operator(ms, jrank, tz, parity, 2);
#ops=Operator(ms, jrank, tz, parity, 2);
t3.ThreeBody.SetMode("pn");


print('eta norm: ', eta.Norm())
print('H norm: ', H.Norm())
cm.comm223ss(eta, H, t3);
print('t3 norm: ', t3.Norm(),t3.IsHermitian())
cm.comm231ss(eta, t3, opo);
cm.comm232ss(eta, t3, opo);


cm.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True);
cm.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True);
cm.FactorizedDoubleCommutator.comm223_231(eta, H, ops);
cm.FactorizedDoubleCommutator.comm223_232(eta, H, ops);


print(opo.Norm(),ops.Norm())
