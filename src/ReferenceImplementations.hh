
///////////////////////////////////////////////////////////////////////////////////
//    ReferenceImplementations.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
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

#ifndef ReferenceImplementations_hh
#define ReferenceImplementations_hh 1

#include "Operator.hh"

namespace ReferenceImplementations
{

  void comm110ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm220ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm111ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm121ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm221ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm122ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_phss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hh_221ss(const Operator &X, const Operator &Y, Operator &Z);

  void comm330ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm331ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm231ss(const Operator &X, const Operator &Y, Operator &Z);

  void comm132ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm232ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm332_ppph_hhhpss(const Operator &X, const Operator &Y, Operator &Z);
  void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z);

  void comm133ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm223ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm233_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm233_phss(const Operator &X, const Operator &Y, Operator &Z);

  void comm333_ppp_hhhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm333_pph_hhpss(const Operator &X, const Operator &Y, Operator &Z);

  // scalar-tensor commutators
  void comm111st(const Operator &X, const Operator &Y, Operator &Z);
  void comm121st(const Operator &X, const Operator &Y, Operator &Z);
  void comm122st(const Operator &X, const Operator &Y, Operator &Z);
  void comm221st(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_phst(const Operator &X, const Operator &Y, Operator &Z);


  // scalar-tensor with a 3b operator
  void comm331st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm223st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm231st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm232st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  // AMC print of comm232st seeds + restore ×1/2 and (1-P). See learn/amc_tts/comm232st/
  void comm232st_amc(const Operator &X, const Operator &Y, Operator &Z);
  // Diagnostic: AMC tex Eq1 only (X_icab Y_abjklc), no 1/2 and no (1-P).
  void comm232st_amc_eq1(const Operator &X, const Operator &Y, Operator &Z);
  void comm133st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm132st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)

  // tensor × tensor → scalar leftover. Naive J-scheme (and production wraps
  // for diagrams with no ss-level BLAS: 121, 132, 231, 223, 232).
  // Production GEMM: 111, 122, 221, 220, 222. Gold: UnitTest::Mscheme_Test_comm*tts.
  // Not-production tested: comm232tts_bare, comm223_231_tts*, comm223_232_tts*.
  // AMC: learn/amc_tts/comm_tts/input/comm*tts_unred.txt
  void comm110tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm220tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm111tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm121tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm122tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm221tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hhtts(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_phtts(const Operator &X, const Operator &Y, Operator &Z);
  void comm231tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm132tts(const Operator &X, const Operator &Y, Operator &Z);
  void comm232tts(const Operator &X, const Operator &Y, Operator &Z);
  // AMC comm232tts Eq. eq (1..4) as printed: no 1/2, no (1-P).
  void comm232tts_bare(const Operator &X, const Operator &Y, Operator &Z, int eq);
  // Tensor 2b × tensor 2b → scalar 3b. AMC primitive + comm223ss 9 perms.
  void comm223tts(const Operator &X, const Operator &Y, Operator &Z);
  
  void comm332_ppph_hhhpst(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test (J and T)
  void comm332_pphhst(const Operator &X, const Operator &Y, Operator &Z);       // PASS the unit test (J and T)
  void comm233_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);      // PASS the unit test (J and T)
  void comm233_phst(const Operator &X, const Operator &Y, Operator &Z);         // PASS the unit test (J and T)
  void comm333_ppp_hhhst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test (J and T)
  void comm333_pph_hhpst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test (J and T)


  /// Two-nested-commutator expressions Z = [X,[X,Y]_3]  where X and Y are 2-body.
  void diagram_CIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIb(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIb(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIc(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIId(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIIb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_DIb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIVa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_DIVb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIVb_intermediate(const Operator &X, const Operator &Y, Operator &Z);

  // The commutators for [Omega, [Omega, Gamma]]
  void comm223_231_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_tts_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_tts(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts_GI(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts_GII(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts_GIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z);
  /// which_term: 0=both, 1=G3c Term1 (χ_k), 2=G3c Term2 (χ_j)
  void comm223_232_tts_GIIIc(const Operator &Eta, const Operator &Gamma, Operator &Z,
                             int which_term = 0);
  void comm223_232_tts_GIIIc_term1(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_tts_GIIIc_term2(const Operator &Eta, const Operator &Gamma, Operator &Z);
  /// Tensor-convention DIRECT: MakeReduced(Γ), reduced χ^θ, AMC reduced 9j on
  /// diagonal (λ_Z=0) channels only, then MakeNotReduced → scalar Z.
  /// AMC: learn/amc_tts/output/G3c_chi_theta_reduced_gamma_ninej.tex
  void comm223_232_tts_GIIIc_tensor_red(const Operator &Eta, const Operator &Gamma,
                                        Operator &Z, int which_term = 0);
  void comm223_232_tts_GIVa(const Operator &Eta, const Operator &Gamma, Operator &Z);
  /// which_term: 0 both, 1 AMC term1 (W1 / Ω_dibc Γ_acdk Ω_jbla), 2 AMC term2 (Ω_dkbc Γ_acdi Ω_lbja).
  void comm223_232_tts_GIVb(const Operator &Eta, const Operator &Gamma, Operator &Z,
                           int which_term = 0);
  void comm223_232_tts_GIVc(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_tts_fI(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_tts_fII(const Operator &Eta, const Operator &Gamma, Operator &Z);
  /// Unreduced χ^γ tab (AMC direct / W1·W2). Gold for f^III_a ladder.
  void build_chi_gamma_tab(const Operator &Eta, std::vector<double> &chi_tab);
  /// Scalar ladder: f^III_a from unreduced χ^γ × Γ.
  void fold_fIIIa_ladder(const Operator &Gamma,
                         const std::vector<double> &chi_tab, Operator &Z);
  void comm223_231_tts_fIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_tts_fIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231sst(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232sst(const Operator &Eta, const Operator &Gamma, Operator &Z);

  void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_132(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_132_ladder(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_132_cross(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_132_onebody(const Operator &Eta, const Operator &Gamma, Operator &Z);

  double TriplesGuess(const Operator &Omega, const Operator &H);
  void comm223_231_fI(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_fIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z);

  // EVC Z1 m-scheme gold (unreduced). t,z excitation-stored. dz.OneBody(a,i) = d z_i^a / dλ.
  void evc_z1_mscheme(const Operator &t, const Operator &z, Operator &dz);
  void evc_z2_mscheme(const Operator &t, const Operator &z, Operator &dz);
  double evc_z0_mscheme(const Operator &t, const Operator &z);

} // namespace ReferenceImplementations

#endif
