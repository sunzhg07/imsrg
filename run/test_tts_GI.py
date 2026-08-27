#!/usr/bin/env python3
"""Benchmark ethS factorized Gamma^I (chi^epsilon) vs unfactorized TTS reference."""

from pyIMSRG import *

emax = 2
ms = ModelSpace(emax, "He4", "He4")
ms.SetHbarOmega(20.0)

unt = UnitTest(ms)

# Tensor anti-Hermitian Omega/Eta, scalar Hermitian Gamma, scalar Z
rank_j, parity, rank_Tz, particle_rank = 2, 0, 0, 2
Eta = unt.RandomOp(ms, rank_j, rank_Tz, parity, particle_rank, -1)  # antiherm
Gamma = unt.RandomOp(ms, 0, 0, 0, particle_rank, +1)  # scalar herm

Z_ref = Operator(ms, 0, 0, 0, 2)
Z_fac = Operator(ms, 0, 0, 0, 2)
Z_ref.SetHermitian()
Z_fac.SetHermitian()

# Unfactorized TTS Gamma^I only
ReferenceImplementations.comm223_232_tts_GI(Eta, Gamma, Z_ref)

# Factorized ethS: chi1b Gamma^I only
cm = Commutator
cm.FactorizedDoubleCommutator_eths.SetUse_1b_Intermediates(True)
cm.FactorizedDoubleCommutator_eths.SetUse_2b_Intermediates(False)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeGI_2b(True)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeGII_2b(False)
cm.FactorizedDoubleCommutator_eths.comm223_232(Eta, Gamma, Z_fac)

diff = Z_ref - Z_fac
print(f"Z_ref 2b norm = {Z_ref.TwoBodyNorm():.8e}")
print(f"Z_fac 2b norm = {Z_fac.TwoBodyNorm():.8e}")
print(f"diff  2b norm = {diff.TwoBodyNorm():.8e}")
print(f"diff  total   = {diff.Norm():.8e}")
tol = 1e-6
ok = diff.TwoBodyNorm() < tol
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
