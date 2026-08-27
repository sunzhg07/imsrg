#!/usr/bin/env python3
"""Benchmark ethS factorized f^III_b (chi^delta RME+DGEMM) vs Case-2 TTS reference.

Reference: reduced Γ (MakeReduced), reduced Ω, Z stored unreduced (÷ĵ²).
Factorized: χ^δ = [Ω⊗Ω]^(0) DGEMM; fold χ·Γ with unreduced Γ (≡ G_red/Ĵ).
"""

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

# Unfactorized TTS f^III_b only
ReferenceImplementations.comm223_231_tts_fIIIb(Eta, Gamma, Z_ref)

# Factorized ethS: chi2b f^III_b only (scalar chi^delta); disable 1b and IIIa
cm = Commutator
cm.FactorizedDoubleCommutator_eths.SetUse_1b_Intermediates(False)
cm.FactorizedDoubleCommutator_eths.SetUse_2b_Intermediates(True)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeI_1b(False)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeII_1b(False)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeIII_1b(True)
cm.FactorizedDoubleCommutator_eths.SetUse_TypeIIIa_1b(False)
cm.FactorizedDoubleCommutator_eths.comm223_231_st(Eta, Gamma, Z_fac)

diff = Z_ref - Z_fac
print(f"Z_ref 1b norm = {Z_ref.OneBodyNorm():.8e}")
print(f"Z_fac 1b norm = {Z_fac.OneBodyNorm():.8e}")
print(f"diff  1b norm = {diff.OneBodyNorm():.8e}")
print(f"diff  total   = {diff.Norm():.8e}")
tol = 1e-6
ok = diff.OneBodyNorm() < tol
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
