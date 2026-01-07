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
#
#### Do the first stage of integration to decouple the core
imsrgsolver.Solve()
#
#### Now set the generator for the second stage to decouple the valence space
imsrgsolver.SetGenerator(valence_generator)
imsrgsolver.SetSmax(smax_valence)

imsrgsolver.Solve()

### Hs is the IMSRG-evolved Hamiltonian
#rw.WriteTokyo( HNO, valence_fname,'')
Hs = imsrgsolver.GetH_s()
Hs.ZeroBody=0.
#gm.force_decouple(Hs)

tdm_op=rdm_el();

tdm_op.read_tdm("he6.ref",ms)

nm=tdm_op.GetVSEOM_Overlap_rd(Hs)

print('reference energy: ', nm, 'vs [ -4.19234292  9.75174484 ] in nmax2 for he6_2',)
eom = EOM(ms,tdm_op)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()

max_iter=3
unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,3,1
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
h4= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)

chi= gm.GetEOM_ladder(h3,0)
chi2= gm.GetEOM_ladder(h4,0)


eom.ProjectOprator(chi)
eom.ProjectOprator(chi2)

cnorm=Norm_vs_new(chi,chi,tdm_op)
chi=chi/np.sqrt(cnorm)

cnorm=Norm_vs_new(chi2,chi2,tdm_op)
chi2=chi2/np.sqrt(cnorm)

cnorm=Norm_vs_new(chi,chi2,tdm_op)
chi2=chi2 - cnorm*chi

cnorm=Norm_vs_new(chi2,chi2,tdm_op)
chi2=chi2/np.sqrt(cnorm)


vec_a = htc_vs(Hs, chi2,eom.ProjectOprator)
cnorm1=Norm_vs_new(chi,vec_a,tdm_op)
print("<1h2> 2b : ",cnorm1)

vec_a = htc_vs(Hs, chi,eom.ProjectOprator)
cnorm2=Norm_vs_new(chi2,vec_a,tdm_op)
print("<2h1> 2b : ",cnorm2)

cnorm1+=Norm3_new(chi,chi2,Hs,ms,tdm_op)
cnorm2+=Norm3_new(chi2,chi,Hs,ms,tdm_op)

print("<1h2/2h1> diff 2b+3b : ",cnorm1-cnorm2,cnorm1, cnorm2)

#print(chi.Norm())
#print(chi2.Norm())
#
chi_back =chi*1.
chi2_back =chi2*1.

Hs.EraseThreeBody()

h3=h3*0.
opa=chi*0.
h3.SetAntiHermitian()
cm.comm223ss(Hs,chi,h3)
cm.comm231ss(chi2,h3,opa)
rst=gm.GetVSEOM_Overlap_rd(opa,tdm_op)
a231a=rst
print('comm231: ',rst)

h3=h3*0.
opa=chi*0.
h3.SetAntiHermitian()
cm.comm223ss(Hs,chi2,h3)
cm.comm231ss(chi,h3,opa)
rst=gm.GetVSEOM_Overlap_rd(opa,tdm_op)
a231b=rst
print('comm231: ',rst)
print('Diff 231: ', a231a-a231b)

Hs.EraseOneBody()

nm1=dcom22232(chi2,chi,Hs,tdm_op)
print(nm1)
nm2=dcom22232(chi,chi2,Hs,tdm_op)
print(nm2)

print('Diff 231: ', nm1-nm2)


#e,v1,v2=lanczos_proc(htc_vs, Norm_vs_new, Norm3_new,Hs, chi, 4, 1,ms,tdm_op,eom.ProjectOprator)
#e,vs,v2=arnoldi_proc_new(htc_vs, Norm_vs_new, Norm3_new, Hs, chi, 17,4,ms,tdm_op,eom.ProjectOprator)
#e,vs,v2=arnoldi_proc(htc_vs, Norm_vs_new, Norm3_new, Hs, chi, 24,2,ms,tdm_op,eom.ProjectOprator)
