//////////////////////////////////////////////////////////////////////////////////
//    EOM.hh, part of  imsrg++
//    Copyright (C) 2021  Ragnar Stroberg
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

#ifndef EOM_hh
#define EOM_hh 1

#include "Generator.hh"
#include "ModelSpace.hh"
#include "Operator.hh"
#include <string>
#include <tuple>
#include <vector>

class EOM {
public:
  // essential fields
  ModelSpace *modelspace; ///< Pointer to the associated modelspace
  Operator Hs;
  // JPT for the EOM, angular momentum, parity and isospin change
  int J2 = 0;
  int parity = 0;
  int itz = 0;

  /// true when an rdm / tdm_file was supplied at construction (MR mode).
  /// false for plain single-reference EOM.
  bool is_multiref = false;

  // The following are used for the multi-reference EOM only
  ModelSpace rdm_modelspace; ///< Owned modelspace for rdm when built from a tdm_file
  ModelSpace *rdm_ms = nullptr; ///< Points to rdm_modelspace (file) or modelspace (Hs)
  std::vector<int> hs_to_rdm_orb; ///< Map: Hs orbit index -> rdm orbit index (-1 if absent)
  Operator rdm;
  arma::sp_mat Nkernel;
  arma::sp_mat Prj_kernel;
  // arma::mat Nkernel;

  arma::vec Energies;
  // configurations
  std::vector<std::array<index_t, 4>> eom_confs;
  index_t qv_start, qv_end, qv_dim;
  index_t ph_start, ph_end, ph_dim;
  index_t ppvv_start, ppvv_end, ppvv_dim;
  index_t pphv_start, pphv_end, pphv_dim;
  index_t pphh_start, pphh_end, pphh_dim;
  size_t channel;
  index_t eom_dims = 0;

  // Methods
  EOM(Operator &Hs, Operator &rdm, int J2, int parity, int itz);
  EOM(Operator &Hs, const std::string &tdm_file, int J2, int parity, int itz);
  EOM(Operator &Hs, int J2, int parity, int itz);

  void force_decouple(Operator &H);

  double GetVSEOM_Overlap_single(Operator &H1, Operator &H2);
  double GetVSEOM_Overlap_multiref(Operator &H);
  Operator GetVSEOM_ladder_single(Operator &H, int herm);
  Operator GetVSEOM_ladder_multiref(Operator &H, int herm);

  void ConstructConfigs();
  void PrintConfigs();
  void ConstructNormMatrix();
  void ConstructProjectMatrix();
  double Core_Diagram(size_t a, size_t b, size_t c, size_t d, size_t e,
                      size_t f, double j1, double j2);
  // Debug/inspection helper: return the uncoupled 3-body diagram entries before RDM contraction.
  std::vector<std::tuple<size_t,size_t,size_t,double>> ThreeBody_Diagram_Entries(size_t a, size_t b, size_t c,
                           size_t d, size_t e, size_t f, size_t g, double j0,
                           double j2);
  // Main API: compute the full diagram and immediately contract it with the 3-body RDM.
  double ThreeBody_Diagram(size_t a, size_t b, size_t c,
                           size_t d, size_t e, size_t f,
                           size_t g, double j0, double j2);

  // Compute the comm223ss diagram for a fixed 3B bra (i,j,k) and ket (l,m,n) with
  // total angular momentum twoJ, using the bra/ket ordering exactly as provided
  // by the caller rather than requiring any canonical a<=b<=c or d<=e<=f sorting.
  // Comparison against the full comm223ss result should be done through GetME_pn,
  // which applies the built-in recoupling for arbitrary external ordering.
  // Sparse seed operators are:
  //   X^{jab}_{sa sb, sf sg} = 1  (all other X elements zero)
  //   Y^{jde}_{sg sc, sd se} = 1  (all other Y elements zero)
  // where sg is the shared intermediate orbit.
  // Returns an arma::mat of size (nJ1 x nJ2) where row r <-> J1 = J1min_out + r
  // and col c <-> J2 = J2min_out + c.
  // herm_X, herm_Y = +1 (hermitian) or -1 (antihermitian).
  arma::mat threebody_diagram_comm(
      size_t i, size_t j, size_t k,
      size_t l, size_t m, size_t n,
      int twoJ,
      size_t sa, size_t sb, size_t sf, size_t sg, int jab,
      size_t sc, size_t sd, size_t se, int jde,
      int herm_X, int herm_Y,
      int &J1min_out, int &J2min_out);
  arma::vec GetEnergies();
  void SqrtMat(arma::mat &Amat, size_t n);
  void ProjectOprator(Operator &Qin);
  void block_svd(std::vector<int> &coupled_vector);
  double ComputeNorm(Operator &Op1, Operator &Op2);

  // -----------------------------------------------
  // Lanczos / Arnoldi eigensolver infrastructure
  // (C++ ports of run/lanczos.py)
  // -----------------------------------------------

  // --- norm helpers ---
  double NormSingle(Operator &T1, Operator &T2);
  double NormMultiref(Operator &T1, Operator &T2);
  double Norm3Multiref(Operator &t1, Operator &t2, Operator &haml);

  // --- H * v helpers (direct member calls, no callbacks) ---
  Operator HtcSingle(Operator &haml, Operator &chi);   ///< [Haml,chi] → ladder_single
  Operator HtcMultiref(Operator &haml, Operator &chi); ///< [Haml,chi] → ladder_multiref + ProjectOprator

  // --- double commutator diagonal (223_231 + 223_232 + 223_132) ---
  std::pair<double, Operator> DcomMultiref(Operator &haml, Operator &chi);

  // --- eigensolvers ---
  // SR Lanczos (mirrors sr_eom.py lanczos_proc with htc_single/norm_single)
  std::pair<arma::vec, std::vector<Operator>>
  LanczosSolve(Operator &vi, int max_iter, int state_want);

  struct ArnoldiResult {
    arma::vec   energies;
    arma::mat   eigvecs;
    std::vector<Operator> ritz;
    arma::mat   hall;
  };

  // MR Arnoldi (mirrors mr_eom.py arnoldi_proc with htc_multiref/norm_multiref/dcom)
  ArnoldiResult ArnoldiSolve(Operator &vi, int max_iter, int state_want);

  /// Read a transition density matrix file and populate a scalar 2-body Operator.
  /// File format (mirrors the Python read_tdm in run/lanczos.py):
  ///   line 0 : J_total (float)
  ///   line 1 : norb (int)
  ///   lines 2..norb+1 : orbit table  "idx n l j2 tz2"
  ///   next line : n_obtd
  ///   n_obtd lines : "_ a b ... rd"  (1-body density matrix elements)
  ///   next line : n_tbtd
  ///   n_tbtd lines : "_ a b c d J ... rd"  (2-body density matrix elements)
  ///   next line : n_3btd
  ///   n_3btd lines : "_ a b c d e f jab jef jtot ... rd"  (3-body density matrix elements)
  Operator ReadTdm(const std::string &tdm_file);

  /// Write the rdm operator to a file in the exact format ReadTdm reads,
  /// iterating over 3-body channels and kets in native memory order.
  void WriteTdm(const Operator &op, const std::string &filename) const;

  /// Build hs_to_rdm_orb map from modelspace -> rdm_ms (must be called after rdm_ms is set).
  void BuildOrbMap();
  /// 1-body rdm element: translates Hs orbit indices to rdm indices, returns 0 if absent.
  double RdmOB(size_t i_hs, size_t j_hs) const;
  /// 2-body rdm element via J + Hs orbit indices, returns 0 if any orbit absent in rdm.
  double RdmTB_J(double J, size_t a_hs, size_t b_hs, size_t c_hs, size_t d_hs) const;
  /// 3-body rdm element <(ab)^Jab c ; J | rho3 | (de)^Jed f ; J> via Hs orbit indices.
  /// Jab, Jed are integer (not 2*J), twoJ = 2*J_total.
  /// Returns 0 if any orbit is absent in the rdm model space.
  double RdmThreeBody_J(int Jab, size_t a_hs, size_t b_hs, size_t c_hs,
                        int Jed, size_t d_hs, size_t e_hs, size_t f_hs, int twoJ) const;

  /// Convenience result bundle returned by Run().
  struct RunResult {
    double      eref;    ///< reference-state energy = <Hs> in the valence space
    ArnoldiResult arnoldi; ///< full Arnoldi output (energies, Ritz vecs, ...)
  };

  /// Unified EOM solve — dispatches on is_multiref:
  ///
  ///  SR (is_multiref == false, mirrors run/sr_eom.py):
  ///    Calls force_decouple, builds a random single-ref initial vector,
  ///    then runs LanczosSolve with HtcSingle / NormSingle.
  ///    RunResult::eref = Hs.ZeroBody.
  ///
  ///  MR (is_multiref == true, mirrors run/mr_eom.py lines 133-173):
  ///    Calls ConstructConfigs / ConstructNormMatrix / ConstructProjectMatrix,
  ///    builds a random projected multiref initial vector, then runs
  ///    ArnoldiSolve with HtcMultiref / NormMultiref / DcomMultiref wired in.
  ///    RunResult::eref = <Hs> in the valence space.
  ///
  /// @param max_iter   maximum Lanczos/Arnoldi steps  (default 200)
  /// @param state_want number of lowest eigenvalues to target (default 6)
  RunResult Run(int max_iter = 200, int state_want = 6);

private:
  std::vector<std::tuple<size_t,size_t,size_t,double>> ThreeBody_Diagram_Entries_Internal(size_t a, size_t b, size_t c,
                           size_t d, size_t e, size_t f, size_t g, double j0,
                           double j2);
  double ThreeBody_Diagram_Internal(size_t a, size_t b, size_t c,
                                    size_t d, size_t e, size_t f,
                                    size_t g, double j0, double j2);
  RunResult RunSR(int max_iter, int state_want);  ///< called by Run() when !is_multiref
  RunResult RunMR(int max_iter, int state_want);  ///< called by Run() when  is_multiref
};

#endif
