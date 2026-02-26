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
  void Setup_rdm();
  void ConstructNormMatrix();
  void ConstructProjectMatrix();
  void SolveEOM();
  double Core_Diagram(size_t a, size_t b, size_t c, size_t d, size_t e,
                      size_t f, double j1, double j2);
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
  RunResult RunSR(int max_iter, int state_want);  ///< called by Run() when !is_multiref
  RunResult RunMR(int max_iter, int state_want);  ///< called by Run() when  is_multiref
};

#endif
