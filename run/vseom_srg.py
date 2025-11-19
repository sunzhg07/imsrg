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


tdm_op=read_tdm("he6.ref",ms)

nm=gm.GetVSEOM_Overlap_rd(Hs,tdm_op)

print('reference energy: ', nm, 'vs [ -4.19234292  9.75174484 ] in nmax2 for he6_2',)
eom = EOM(ms,tdm_op)
eom.ConstructConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()

max_iter=3
unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
chi= gm.GetEOM_ladder(h3,0)


cnorm=Norm_vs_new(chi,chi,tdm_op)
print('Original norm: ', cnorm)
eom.ProjectOprator(chi)
cnorm=Norm_vs_new(chi,chi,tdm_op)
print('Norm without null vectors: ', cnorm)
chi=chi/np.sqrt(cnorm)

e,v1,v2=lanczos_proc(htc_vs, Norm_vs_new, Norm3,Hs, chi, 24, 1,ms,tdm_op)
#e,vs,v2=arnoldi_proc(htc_vs, Norm_vs_new5,Norm3, Hs, chi, 14,2,ms,tdm_op)
### Generate all configurations
#
##cfs=[]
##
#### 
##for itz in [1]:
##    for pp in [0]:
##        for jj in [0,2]:
##            ch=ms.GetTwoBodyChannelIndex(jj,pp,itz)
##            print(jj,pp,itz,ch,'channel')
##            kcf=ms.GetTwoBodyChannel(ch)
##            bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()
##            kets=kcf.GetKetIndex_vv()
##            for ibra in bras:
##                dbra=kcf.GetKet(ibra)
##                for iket in kets:
##                    dket=kcf.GetKet(iket)
##                    cfs.append([dbra.p,dbra.q,dket.p,dket.q,jj,pp,itz])
##print(len(cfs))
######## ppvh
#####nch=ms.GetNumberTwoBodyChannels()
#####
#####for ich in range(nch):
#####    dch = ms.GetTwoBodyChannel(ich)
#####    jj=dch.J
#####    pp=dch.parity
#####    tt=dch.Tz
#####    ch=ich
#####
#####    kcf=dch
#####    if(kcf.GetNumberKets() == 0):
#####        continue
#####    bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()+kcf.GetKetIndex_vv()
#####    kets=kcf.GetKetIndex_vc()
#####    for ibra in bras:
#####        dbra=kcf.GetKet(ibra)
#####        for iket in kets:
#####            dket=kcf.GetKet(iket)
#####            cfs.append([dbra.p,dbra.q,dket.p,dket.q,jj,pp,tt])
##dims=len(cfs)
##
##print(dims)
##
##h3mat=np.zeros([dims,dims])
##nmat=np.zeros([dims,dims])
##print(cfs)
##
##
##for i,cfl in enumerate(cfs):
##    ql=Hs*0.
##    ql.SetTwoBody(cfl[4],cfl[5],cfl[6],cfl[4],cfl[5],cfl[6], cfl[0],cfl[1],cfl[2],cfl[3],1.)
##    for j,cfr in enumerate(cfs):
##        qr=Hs*0.
##        qr.SetTwoBody(cfr[4],cfr[5],cfr[6],cfr[4],cfr[5],cfr[6], cfr[0],cfr[1],cfr[2],cfr[3],1.)
##        nmat[i,j] = Norm_vs_new4(ql,qr,rdms)
##        qv=Hs*0.
##        qv=htc_vs(Hs, qr)
##        h3mat[i,j]=Norm4(ql,qr,Hs,ms,tdm_op)
##np.save("nmat_he6_3b.npy", nmat)
##np.save("h3mat_he6_3b.npy",h3mat)
###
##hmatbare=[]
##for i,cfl in enumerate(cfs):
##    ql=Hs*0.
##    ql.SetTwoBody(cfl[4],cfl[5],cfl[6],cfl[4],cfl[5],cfl[6], cfl[0],cfl[1],cfl[2],cfl[3],1.)
##    qv=htc_vs(Hs, ql)
##    qout=print_op(qv,ms)
##    hmatbare.append(qout)
##hmatbare=np.array(hmatbare)
##np.save("hmat_he6_3b.npy", hmatbare)
#
#
