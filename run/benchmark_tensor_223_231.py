#!/usr/bin/env python3
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))

from pyIMSRG import *


emax = int(os.environ.get("EMAX", "1"))
ref = os.environ.get("REF", "He4")
val = os.environ.get("VAL", "p-shell")
hw = int(os.environ.get("HW", "16"))
rank_j = int(os.environ.get("JRANK", "1"))
rank_t = int(os.environ.get("TRANK", "0"))
parity = int(os.environ.get("PARITY", "0"))
particle_rank = int(os.environ.get("PARTICLE_RANK", "2"))
herm = int(os.environ.get("HERM", "1"))
nrepeat = int(os.environ.get("NREPEAT", "1"))


def zero_scalar_like(op):
    z = op * 0.0
    z.SetHermitian()
    return z


def print_norms(label, op):
    print(
        f"{label:>18}: "
        f"Norm={op.Norm():.12e} "
        f"1b={op.OneBodyNorm():.12e} "
        f"2b={op.TwoBodyNorm():.12e} "
        f"3b={op.ThreeBodyNorm():.12e}"
    )


ms = ModelSpace(emax, ref, val)
ms.SetHbarOmega(hw)
unt = UnitTest(ms)

eta = unt.RandomOp(ms, rank_j, rank_t, parity, particle_rank, herm)
gamma = unt.RandomOp(ms, 0, 0, 0, particle_rank, herm)

print("benchmarking tensor eta comm223_231st against comm223st -> comm231st")
print(
    f"EMAX={emax} REF={ref} VAL={val} JRANK={rank_j} "
    f"TRANK={rank_t} PARITY={parity} NREPEAT={nrepeat}"
)
print_norms("eta", eta)
print_norms("gamma", gamma)


def comm223_231st_chain_reference(eta, gamma, out):
    t3 = Operator(ms, rank_j, rank_t, parity, 3)
    t3.ThreeBody.SetMode("pn")
    t3 *= 0.0
    t3.SetHermitian()

    Commutator.comm223st(gamma, eta, t3)
    Commutator.comm231st(eta, t3, out)
    return t3


for irepeat in range(nrepeat):
    z_chain = zero_scalar_like(gamma)
    t0 = time.perf_counter()
    t3 = comm223_231st_chain_reference(eta, gamma, z_chain)
    t_chain = time.perf_counter() - t0

    z_cpp_chain = zero_scalar_like(gamma)
    t0 = time.perf_counter()
    ReferenceImplementations.comm223_231st(eta, gamma, z_cpp_chain)
    t_cpp_chain = time.perf_counter() - t0

    z_public = zero_scalar_like(gamma)
    t_public = None
    public_error = None
    t0 = time.perf_counter()
    try:
        Commutator.FactorizedDoubleCommutator.comm223_231st(eta, gamma, z_public)
        t_public = time.perf_counter() - t0
    except RuntimeError as err:
        public_error = str(err)

    zdiff_cpp_chain = z_cpp_chain - z_chain
    zdiff_public = z_public - z_chain
    rel_cpp_chain = zdiff_cpp_chain.Norm() / max(z_chain.Norm(), 1.0e-30)
    rel_public = zdiff_public.Norm() / max(z_chain.Norm(), 1.0e-30)

    print(f"\nrepeat {irepeat + 1}/{nrepeat}")
    print(f"python chain time:       {t_chain:.6f} s")
    print(f"c++ chain time:          {t_cpp_chain:.6f} s")
    if t_public is not None:
        print(f"public implementation time: {t_public:.6f} s")
    else:
        print("public implementation: skipped (disabled)")
    print_norms("threebody int", t3)
    print_norms("python chain", z_chain)
    print_norms("c++ chain", z_cpp_chain)
    print_norms("public", z_public)
    print_norms("cpp-python", zdiff_cpp_chain)
    print_norms("public-chain", zdiff_public)
    print(f"cpp-python relative Norm: {rel_cpp_chain:.12e}")
    print(f"public-chain relative Norm: {rel_public:.12e}")
    if public_error is not None:
        print(f"public disable reason: {public_error}")
