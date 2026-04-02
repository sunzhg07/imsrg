#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'build'))
from pyIMSRG import *
import numpy as np
from lanczos import *



emax =3         # maximum number of oscillator quanta in the model space
ref = 'He8'     # reference used for normal ordering
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
Hs.UndoNormalOrdering()
ms2 = ModelSpace(emax,'He4',val)
Hs.SetModelSpace(ms2)
Hs.DoNormalOrdering()


#tdm_op = read_tdm('he8_pshell_e3_em1290_gs.dat', ms2)
#
#rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1
#unt = UnitTest(ms2)
#h1    = unt.RandomOp(ms2, rank_j, rank_Tz, parity, particle_rank, herm)
#h2    = unt.RandomOp(ms2, rank_j, rank_Tz, parity, particle_rank, herm)
#
#
#eom = EOM(Hs, tdm_op, rank_j, parity, rank_Tz)
#eom.ConstructConfigs()
#eom.ConstructNormMatrix()
#eom.ConstructProjectMatrix()
#
#chi_a = eom.GetVSEOM_ladder_multiref(h1, 1)
#chi_b = eom.GetVSEOM_ladder_multiref(h2, 1)
#
#nm=norm_multiref(eom,chi_a,chi_b)
#print('norm comm: ', nm)
#nm=eom.ComputeNorm(chi_a,chi_b)
#print('norm matmul: ', nm)

#eom.ProjectOprator(chi_a)
#eom.ProjectOprator(chi_b)


## ---------------------------------------------------------------
eom = EOM(Hs, 'he8.ref', 0, 0, 0)
#print("3b allocated:", eom.rdm.ThreeBody.IsAllocated())
eom.WriteTdm(eom.rdm,'ss')
#res = eom.Run(26, 4)
## ---------------------------------------------------------------
#print(f'\n  [he6.ref]  E_ref = {res.eref:.6f} MeV')
#for k, e in enumerate(res.arnoldi.energies):
#    print(f'    E({k}): excitation={e:.4f}  absolute={e+res.eref:.4f} MeV')

#cm.SetIMSRG3Noqqq(True)
###cm.SetIMSRG3Onlyvvv(True)
##
##cm  = Commutator
##unt = UnitTest(ms)
##
##ref_files = [f'he6_t{round(i*0.1,1):.1f}.ref' for i in range(-13, 14)]
##results   = {}   # {ref_file: (eref, e_array, hall_nb, v2, eom)}
##
##for ref_file in ref_files:
##    print(f'\n{"="*60}')
##    print(f'  Reference state: {ref_file}')
##    print(f'{"="*60}')
##
##    tdm_op = read_tdm(ref_file, ms)
##
##    rank_j, parity, rank_Tz, particle_rank, herm = 0, 0, 0, 2, 1
##    eom = EOM(Hs, tdm_op, rank_j, parity, rank_Tz)
##    eom.ConstructConfigs()
##    eom.ConstructNormMatrix()
##    eom.ConstructProjectMatrix()
##    nm   = eom.GetVSEOM_Overlap_multiref(Hs)
##    eref = nm 
##    print(f'  E_ref (valence) = {nm:.6f}   ZeroBody = {Hs.ZeroBody:.6f}   E_ref total = {eref:.6f} MeV')
##
##    h1    = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)
##    h2    = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, herm)
##    chi_a = eom.GetVSEOM_ladder_multiref(h1, 1)
##    chi_b = eom.GetVSEOM_ladder_multiref(h2, 1)
##    eom.ProjectOprator(chi_a)
##    eom.ProjectOprator(chi_b)
##    e, vs, v2, hall_nb = arnoldi_proc(
##        htc_multiref, norm_multiref, Hs, chi_b,
##        max_iter=200, state_want=6, ms=ms, eom=eom,
##        norm_three=lambda eom, a, b, h, ms: dcom222312(eom, h, a)[0],
##        rdmat=tdm_op, prjop=eom.ProjectOprator,
##        restart_gen=lambda: eom.GetVSEOM_ladder_multiref(unt.RandomOp(ms, 0, 0, 0, 2, 1), 1))
##
##    print(f'\n  [{ref_file}] Arnoldi eigenvalues:')
##    print(f'  E_ref = {eref:.6f} MeV')
##    for k in range(len(e)):
##        print(f'    E({k}): excitation={e[k]:.4f}  absolute={e[k]+eref:.4f} MeV')
##
##    results[ref_file] = (eref, e.copy(), hall_nb, v2, eom)
##
### ========================================================
###  Cross-reference comparison (EOM reference-independence)
###  For an exact EOM, E_ref1 + omega_k1 == E_ref2 + omega_k2
###  Deviations show the truncation error of the 2-body EOM.
### ========================================================
##print(f'\n{"="*60}')
##print('  COMPARISON: absolute energies (E_ref + excitation)')
##print('  EOM is exact iff both columns are identical')
##print(f'{"="*60}')
##refs = list(results.keys())
##n_states = min(len(results[r][1]) for r in refs)
##header = f'  {"State":<8}' + ''.join(f'  {r:<22}' for r in refs) + '  diff (MeV)'
##print(header)
##print('  ' + '-' * (len(header) - 2))
##for k in range(n_states):
##    vals = [results[r][0] + results[r][1][k] for r in refs]
##    diff = vals[-1] - vals[0]
##    row  = f'  E({k})    ' + ''.join(f'  {v:<22.4f}' for v in vals) + f'  {diff:+.4f}'
##    print(row)
##
##print(f'\n  E_ref values:')
##for r in refs:
##    print(f'    {r}: {results[r][0]:.6f} MeV')



