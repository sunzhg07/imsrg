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
nm=eom.GetVSEOM_Overlap_multiref(Hs)
print('reference energy: ', nm, 'vs [ -4.19234292  9.75174484 ] in nmax2 for he6_2',)
##
##
##
##

unt = UnitTest(ms)
### set anti hermitian
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1

h1= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
h2= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
##t3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, 3,herm)
##t3.ThreeBody.SetMode("pn")
#
chi_a=eom.GetVSEOM_ladder_multiref(h1,1)
chi_b=eom.GetVSEOM_ladder_multiref(h2,1)
eom.ProjectOprator(chi_a)
eom.ProjectOprator(chi_b)

#nm=norm_multiref(eom,chi_a,chi_b)
#print('original py ab: ',nm)
#
#nm=norm_multiref(eom,chi_b,chi_a)
#print('original py ba: ',nm)
#nm=eom.ComputeNorm(chi_a,chi_b)
#print('original c++ ab: ', nm)
#nm=eom.ComputeNorm(chi_b,chi_a)
#print('original c++ ba: ', nm)
#
#
#nm=eom.ComputeNorm(chi_a,chi_a)
#print('original c++ aa: ', nm)
#nm=eom.ComputeNorm(chi_a,chi_a)
#print('original c++ aa: ', nm)
#
#chi_b=1.0*chi_a
#
#eom.ProjectOprator(chi_a)
#chi_b=chi_b-chi_a
#
#nm=eom.ComputeNorm(chi_a,chi_a)
#print('c++ after proj aa: ', nm)
#
#nm=eom.ComputeNorm(chi_b,chi_b)
#print('c++ after proj 1-aa: ', nm)

#print("vecryfy p*p=p")
#eom.ProjectOprator(chi_a)
#eom.ProjectOprator(chi_b)
#
#nm=eom.ComputeNorm(chi_a,chi_b)
#print('c++ after proj: ', nm)
#nm=norm_multiref(eom,chi_a,chi_b)
#print('python after proj: ',nm)


#eom.ProjectOprator(chi_b)
#
#
#hb=htc_multiref(eom, Hs, chi_b)
#hab = norm_multiref(eom,chi_a,hb)
#h3ab = norm3_multiref(eom,chi_a,chi_b,Hs,ms)
#print("hab: ", hab, "h3ab: ", h3ab)
#print("hab+h3ab: ", hab+h3ab)
#
#
#ha = htc_multiref(eom, Hs, chi_a)
#hba = norm_multiref(eom, chi_b,ha)
#h3ba = norm3_multiref(eom,chi_b,chi_a,Hs,ms)
#print("hba: ", hba, "h3ba: ", h3ba)
#print("hba+h3ba: ", hba+h3ba)
#print("diff:", hab-hba)
#
## ---------------------------------------------------------------
## Test: does projecting the Arnoldi vector break hermiticity?
## Take chi_a (projected), compute w = H1*chi_a (raw, no prjop).
## Then compare H[chi_a, w] vs H[chi_a, w_proj] where w_proj = P*w.
## If projection breaks hermiticity, H[chi_a,w_proj] != H[w_proj,chi_a].
## ---------------------------------------------------------------
##print("\n--- hermiticity test: unprojected vs projected Arnoldi vector ---")
##w_raw  = htc_multiref(eom, Hs, chi_a)          # H1*chi_a, no projection
##w_proj = w_raw * 1.0
##eom.ProjectOprator(w_proj)                       # P * (H1*chi_a)
##nm=norm_multiref(eom, chi_a, w_raw)
##print(nm)
##nm=norm_multiref(eom, chi_a, w_proj)
##print(nm)
#
## matrix elements with raw w
##H_ab_raw = norm_multiref(eom, chi_a, w_raw) # + norm3_multiref(eom, chi_a, w_raw,  Hs, ms)
##H_ba_raw = norm_multiref(eom, w_raw,  chi_a)# + norm3_multiref(eom, w_raw,  chi_a, Hs, ms)
##print(f"  raw  w:  H[chi_a,w]={H_ab_raw:.8f}  H[w,chi_a]={H_ba_raw:.8f}  diff={H_ab_raw-H_ba_raw:.3e}")
##
### matrix elements with projected w
##H_ab_prj = norm_multiref(eom, chi_a, w_proj) #+ norm3_multiref(eom, chi_a, w_proj, Hs, ms)
##H_ba_prj = norm_multiref(eom, w_proj, chi_a) #+ norm3_multiref(eom, w_proj, chi_a, Hs, ms)
##print(f"  proj w:  H[chi_a,w]={H_ab_prj:.8f}  H[w,chi_a]={H_ba_prj:.8f}  diff={H_ab_prj-H_ba_prj:.3e}")
## ---------------------------------------------------------------
#
##n3a =norm3_multiref(eom, chi_a,chi_a,Hs,ms)
##n3b =norm3_multiref(eom, chi_b,chi_b,Hs,ms)
##n3ab =norm3_multiref(eom, chi_a+chi_b,chi_a+chi_b,Hs,ms)
#
##n3a , oprs = dcom222312(eom, Hs,chi_a)
##n3b , oprs = dcom222312(eom, Hs,chi_b)
##n3ab , oprs = dcom222312(eom, Hs,chi_a+chi_b)
##
##print("n3a: ",n3a )
##print("n3b: ",n3b )
##print("n3ab: ",n3ab )
##
##fab=(n3ab-n3a-n3b-hab+hba)/2
##print("fab: ", fab)
##
##fba=(n3ab-n3a-n3b+hab-hba)/2
##print("fba: ", fba)
#
#
#
e,vs,v2=arnoldi_proc(
   htc_multiref,
   norm_multiref,
   Hs,
   chi_b,
   max_iter=200,
   state_want=6,
   ms=ms,
   eom=eom,
   norm_three=lambda eom,a,b,h,ms: dcom222312(eom,h,a)[0],  # diagonal H2 via polarization identity
   rdmat=tdm_op,
   prjop=eom.ProjectOprator)


