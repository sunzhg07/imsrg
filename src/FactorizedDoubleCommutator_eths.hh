///////////////////////////////////////////////////////////////////////////////////
//    FactorizedDoubleCommutator.hh, part of  imsrg++
//    Copyright (C) 2023 Bingcheng He and Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef FactorizedDoubleCommutator_eths_hh
#define FactorizedDoubleCommutator_eths_hh 1

#include "Operator.hh"
#include <vector>

namespace Commutator {

namespace FactorizedDoubleCommutator_eths {
extern bool use_1b_intermediates;
extern bool use_2b_intermediates;

// Per-diagram switches (default on, every λ). Equations and contraction
// paths: learn/factorized_code_analyze.tex, Sec. jscheme and Sec. opt.
extern bool use_TypeI_1b;     // f^I     / χ^α  (scalar 1b)
extern bool use_TypeII_1b;    // f^II    / χ^β  (tensor 1b)
extern bool use_TypeIII_1b;   // f^III_b / χ^δ  (scalar 2b)
extern bool use_TypeIIIa_1b;  // f^III_a / χ^γ  (scalar 2b, Pandya)
extern bool use_TypeGI_2b;    // Γ^I     / χ^ε  (scalar 1b)
extern bool use_TypeGII_2b;   // Γ^II    / χ^ΩΓ (bra), χ^ζ (ket), tensor 1b
extern bool use_TypeGIIIa_2b; // Γ^III_a / χ^η  ladder
extern bool use_TypeGIIIb_2b; // Γ^III_b / χ^η  Pandya → RC → Γ̄·RC → inverse
extern bool use_TypeGIIIc_2b; // Γ^III_c / χ^θ  Pandya → χ̄·Γ̄ → inverse
extern bool use_TypeGIVa_2b;  // Γ^IV_a  / χ^κ  (tensor 2b) ladder
extern bool use_TypeGIVb_2b;  // Γ^IV_b  / χ^ι  (tensor 2b) leftover Pandya GEMM
extern bool use_TypeGIVc_2b;  // Γ^IV_c  / χ^λ  (tensor 2b) Pandya → DGEMM
/// Isolated-test helper: 0=both χ^λ terms, 1=T1 (ΓΩ, occ on l), 2=T2 (ΩΓ, occ on j).
extern int givc_chi_which;
extern bool givc_fold_as;
void SetGIVcChiWhich(int which);
void SetGIVcFoldAS(bool tf);
void SetUse_1b_Intermediates(bool tf);
void SetUse_2b_Intermediates(bool tf);

void SetUse_TypeI_1b(bool tf);
void SetUse_TypeII_1b(bool tf);
void SetUse_TypeIII_1b(bool tf);
void SetUse_TypeIIIa_1b(bool tf);
void SetUse_TypeGI_2b(bool tf);
void SetUse_TypeGII_2b(bool tf);
void SetUse_TypeGIIIa_2b(bool tf);
void SetUse_TypeGIIIb_2b(bool tf);
void SetUse_TypeGIIIc_2b(bool tf);
/// Isolated-test helper: treat MEs as reduced-convention → MakeNotReduced.
void ForceScalarMakeNotReduced(Operator &Z);
void SetUse_TypeGIVa_2b(bool tf);
void SetUse_TypeGIVb_2b(bool tf);
void SetUse_TypeGIVc_2b(bool tf);
// Leftover [Ω,[Ω,Γ]_3b] for tensor Ω (reduced; λ=0 converted on entry),
// scalar unreduced Γ and Z. Together with comm223_132_tts this equals
// comm231tts + comm232tts + comm132tts applied to W = −comm223st(Γ, Ω).
void comm223_231_st(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z);

void comm223_231_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
							  Operator &Z);
void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
							  Operator &Z);
void comm223_232_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
							  Operator &Z);
void comm223_232_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z);
// Single-diagram entry points (used by tests); each applies ReducedRMEView.
void comm223_232_GIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_GIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_GIIIc(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_GIVa(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_GIVb(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_GIVc(const Operator &Eta, const Operator &Gamma, Operator &Z);
/// Tensor ω, scalar H: leftover 2b of [ω1,[ω2,H2]_3]_2 without storing 3b.
/// Nested gold: W = −comm223st(H, ω), then comm132tts(ω, W).
/// λ=0 unreduced Ω is MakeReduced on entry (ReducedRMEView), then this kernel.
void comm223_132_tts(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_132_tts_ladder(const Operator &Eta, const Operator &Gamma,
                            Operator &Z);
void comm223_132_tts_onebody(const Operator &Eta, const Operator &Gamma,
                             Operator &Z);
void comm223_132_tts_cross(const Operator &Eta, const Operator &Gamma,
                           Operator &Z);
} // namespace FactorizedDoubleCommutator_eths
} // namespace Commutator

#endif
