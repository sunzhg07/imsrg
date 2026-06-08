#!/usr/bin/env python3
import numpy as np
from pyIMSRG import *

emax = 2  # maximum number of oscillator quanta in the model space
ref = "H2"  # reference used for normal ordering
val = ref  # valence space

core_generator = "white"  # definition of generator eta for decoupling the core (could also use 'white')
smax_core = 50  # limit of integration in flow parameter s for first stage of decoupling

f2b = "../../input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
f2e1, f2e2, f2l = 14, 28, 14
# f3b='../../input/NO2B_ThBME_EM7.5_1.8_2.0_IS_hw16from16_ms14_28_18.me3j.gz'
f3b = "none"
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
imsrgsolver.SetDenominatorPartitioning("M")
imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
imsrgsolver.Solve()


### Hs is the IMSRG-evolved Hamiltonian
Hs = imsrgsolver.GetH_s()


## My work for Dipole Polarizability is done below
## Obtaining E1 operator and performing similarity transofrmation
E1 = OperatorFromString(ms, "E1")
E1.PrintOneBody()
E1T = imsrgsolver.Transform(E1)
E1T.PrintOneBody()
print(
    "E1 Info: Parity",
    E1T.GetParity(),
    "\tJRank",
    E1T.GetJRank(),
    "\tTRank",
    E1T.GetTRank(),
    "\tParticleRank",
    E1T.GetParticleRank(),
)
# E1T.PrintOneBody()

print("Hs:")
# Hs.PrintOneBody()


H_d = imsrgsolver.GetH_sDiagonal(Hs)
# H_d = Hs
# EOM.force_decouple(H_d)

H_od = Hs - H_d
print("H_od is:\n\n")
# H_od.PrintOneBody()
Gen = imsrgsolver.generator

omega_k = E1T*0.0
omega_k.SetAntiHermitian()



Niter = 100

# Brody-DIIS controls
diis_start = 3
diis_dim = 6
diis_eps = 1e-12
diis_conv_tol = 1e-12
omega_hist = []
res_hist = []


def op_inner(op_a, op_b):
    # Use polarization identity to get an inner product from operator norms.
    apb = op_a + op_b
    return 0.5 * (apb.Norm() ** 2 - op_a.Norm() ** 2 - op_b.Norm() ** 2)


def diis_extrapolate(ops, residuals):
    m = len(residuals)
    if m < 2:
        return None

    B = np.empty((m + 1, m + 1), dtype=float)
    B.fill(-1.0)
    B[m, m] = 0.0
    rhs = np.zeros(m + 1, dtype=float)
    rhs[m] = -1.0

    for i in range(m):
        for j in range(m):
            B[i, j] = op_inner(residuals[i], residuals[j])

    # Regularize nearly singular DIIS systems.
    for i in range(m):
        B[i, i] += diis_eps

    try:
        coeff = np.linalg.solve(B, rhs)[:m]
    except np.linalg.LinAlgError:
        return None

    mixed = ops[0] * 0.0
    for c, op in zip(coeff, ops):
        mixed += op * float(c)
    mixed.SetAntiHermitian()
    return mixed


for it in range(Niter):
    omega_old = Operator(omega_k) * 1.0
    A = imsrgsolver.GetA(E1T, H_od, omega_old)
    omega_next = Operator(omega_k) * 0.0
    Gen.UpdateGeneral(A, H_d, omega_next)
    
    residual = omega_next - omega_old
    res_norm = residual.Norm()
    
    omega_k = (-1)*omega_next
    #print("Norm of eta at iteration ", it, "is:  ", res_norm)
    if it % 10 == 0:
        print(f"Iter {it+1}, res_norm: {res_norm:.4e}")
    if res_norm < diis_conv_tol:
        print(f"Converged at iteration {it+1}.")
        break
    
    
# This will now contain the cleanly evolved generator matrix elements
final_eta = Operator(omega_k)
O = imsrgsolver.CheckWork(omega_k, Hs)
print("Final Norm of commutator: ", O.Norm())
O.PrintOneBody()
print("Norm we started with: ", E1T.Norm())
E1T.PrintOneBody()

