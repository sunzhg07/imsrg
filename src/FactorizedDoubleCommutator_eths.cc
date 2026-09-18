
#include "FactorizedDoubleCommutator_eths.hh"
#include "Commutator.hh"
#include "ReferenceImplementations.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <omp.h>
#include <map>
#include <cstdint>
#include <array>
#include <deque>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
namespace Commutator {

namespace FactorizedDoubleCommutator_eths {

void comm223_231_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z);
void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z);

bool use_1b_intermediates = true;
bool use_2b_intermediates = true;
bool use_TypeI_1b = true;
bool use_TypeII_1b = true;
bool use_TypeIII_1b = true;
bool use_TypeIIIa_1b = true; // Pandya+DGEMM f^III_a (paper B4c/B5c)
bool use_TypeGI_2b = true;  // Gamma^I / chi^epsilon
bool use_TypeGII_2b = true; // Gamma^II / chi^zeta (universal AMC Path B)
bool use_TypeGIIIa_2b = true; // Gamma^III_a Path B χ^η (2n×2n) → Chi_AS×Γ
bool use_TypeGIIIb_2b = true; // Gamma^III_b Path B Fac Pandya χ → RC → DGEMM
bool use_TypeGIIIc_2b = true; // Gamma^III_c Factorized IIe/IIf (λ=0)
bool use_TypeGIVa_2b = true; // Gamma^IV_a / chi^kappa Factorized DGEMM (λ=0)
bool use_TypeGIVb_2b = true;  // Gamma^IV_b / chi^iota Factorized DGEMM (λ=0)
bool use_TypeGIVc_2b = true;  // Gamma^IV_c / chi^lambda (enable diagram)
int givc_chi_which = 0;       // 0 both, 1 T1 ΓΩ, 2 T2 ΩΓ
bool givc_fold_as = true;     // false: recouple K only (no ½(1−P)²)
void SetUse_1b_Intermediates(bool tf) { use_1b_intermediates = tf; }
void SetUse_2b_Intermediates(bool tf) { use_2b_intermediates = tf; }

void SetUse_TypeI_1b(bool tf) { use_TypeI_1b = tf; }

void SetUse_TypeII_1b(bool tf) { use_TypeII_1b = tf; }

void SetUse_TypeIII_1b(bool tf) { use_TypeIII_1b = tf; }

void SetUse_TypeIIIa_1b(bool tf) { use_TypeIIIa_1b = tf; }
void SetUse_TypeGI_2b(bool tf) { use_TypeGI_2b = tf; }

void SetUse_TypeGII_2b(bool tf) { use_TypeGII_2b = tf; }

void SetUse_TypeGIIIa_2b(bool tf) { use_TypeGIIIa_2b = tf; }
void SetUse_TypeGIIIb_2b(bool tf) { use_TypeGIIIb_2b = tf; }
void SetUse_TypeGIIIc_2b(bool tf) { use_TypeGIIIc_2b = tf; }
void ForceScalarMakeNotReduced(Operator &Z) {
  if (Z.GetJRank() != 0)
    return;
  if (not Z.IsReduced())
    Z.is_reduced = true;
  Z.MakeNotReduced();
}
void SetUse_TypeGIVa_2b(bool tf) { use_TypeGIVa_2b = tf; }
void SetUse_TypeGIVb_2b(bool tf) { use_TypeGIVb_2b = tf; }
void SetUse_TypeGIVc_2b(bool tf) { use_TypeGIVc_2b = tf; }
void SetGIVcChiWhich(int which) { givc_chi_which = which; }
void SetGIVcFoldAS(bool tf) { givc_fold_as = tf; }

// Pandya CC Tz (TensorCommutators / ModelSpace::CalculatePandyaLookup).
// CC Tz = |tz1−tz2|/2 ∈ {0,1}. Rank T:
//   0 → (0,0),(1,1);  1 → (0,1),(1,0);  2 → (1,1).
static bool CcTzCouples(int tz_bra, int tz_ket, int rank_T) {
  return (tz_bra + tz_ket == rank_T) or
         (std::abs(tz_bra - tz_ket) == rank_T);
}

// FDC.cc comm223_232_chi2b ~1227–1393: one Pandya of Ω and Γ, then occupancy
// nnn (not a second transform). Tensor leftover is the same except Ω is
// rank-λ (rectangular CC, 9j) and some χ need two copies (non-Hermitian).
namespace {
using CcPair = std::array<int, 2>;
using CcMatMap = std::map<CcPair, arma::mat>;
struct Tensor232Bars {
  std::vector<CcPair> cc_pairs;
  CcMatMap bar_Omega;   // IMSRG DoTensorPandya(adcb)
  CcMatMap nnn_AbarBC;  // FDC nnnbar_Eta  → χ^ι / GIVb
  CcMatMap nnn_ABbarD;  // FDC nnnbar_Eta_d → χ^κ / GIVa
  const std::deque<arma::mat> *bar_Gamma = nullptr;
};
void FillTensor232Bars(const Operator &Eta, const Operator &Gamma, Operator &Z,
                       Tensor232Bars &B);
const std::deque<arma::mat> &CachedBarGammaScalarCC(const Operator &Gamma,
                                                    Operator &Z, int hGamma);
void comm223_232_GIVa_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const Tensor232Bars &B);
void comm223_232_GIVb_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const Tensor232Bars &B);
} // namespace
// factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
// Eta is tensor (reduced). Gamma and Z are scalar unreduced.
void comm223_231_st(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  if (use_1b_intermediates)
    comm223_231_chi1b_tensor(Eta, Gamma, Z);
  if (use_2b_intermediates)
    comm223_231_chi2b_tensor(Eta, Gamma, Z);
} // comm223_231_st

////////////////////////////////////////////////////////////////////////////
/// ethS one-body intermediates:
///   chi^alpha (scalar) -> f^I   [use_TypeI_1b]
///   chi^beta  (tensor) -> f^II  [use_TypeII_1b]
/// M-scheme: factorized_code_analyze.tex
/// AMC factorized: learn/amc_tts/factored_fII/output/
////////////////////////////////////////////////////////////////////////////
void comm223_231_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {

  double t_internal = omp_get_wtime();
  double t_start = omp_get_wtime();

  Z.modelspace->PreCalculateSixJ();

  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);

  int max_j2 = 0;
  for (auto x : Z.modelspace->all_orbits)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
  int max_J = max_j2;

  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };

  int norbits = Z.modelspace->all_orbits.size();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),
                                  Z.modelspace->all_orbits.end());

  // ==================================================================
  // f^I via scalar chi^alpha
  // ==================================================================
  if (use_TypeI_1b) {
    auto Chi_alpha = Z.OneBody;
    Chi_alpha.zeros();

    // χ^α: [Ω⊗Ω]^(0) via RME+DGEMM → scalar 2b on equal-J channels,
    // then trace spectator c:  (c,d)×(c,e) → 1b (d,e), ÷ĵ_d² (unreduced).
    // occ: n̄_a n̄_b n_c n_d − n_a n_b n̄_c n̄_d
    //   = n_c n_d T_pp − n̄_c n̄_d T_hh,  T = Σ_{J'} (−1)^{J+J'+λ} λ̂^{-1} Ω W Ω
    for (int J0 = 0; J0 <= max_J; ++J0) {
      std::vector<std::array<index_t, 2>> pairs0;
      pairs0.reserve((size_t)norbits * norbits);
      for (auto i : allorb_vec) {
        Orbit &oi = Z.modelspace->GetOrbit(i);
        const double ji = oi.j2 * 0.5;
        for (auto j : allorb_vec) {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          if (not AngMom::Triangle(ji, oj.j2 * 0.5, (double)J0))
            continue;
          pairs0.push_back({i, j});
        }
      }
      const int n0 = (int)pairs0.size();
      if (n0 < 1)
        continue;

      arma::mat T_pp(n0, n0, arma::fill::zeros);
      arma::mat T_hh(n0, n0, arma::fill::zeros);

      for (int J1 = 0; J1 <= max_J; ++J1) {
        if (not AngMom::Triangle(J0, J1, lambda))
          continue;

        std::vector<std::array<index_t, 2>> pairs1;
        pairs1.reserve((size_t)norbits * norbits);
        for (auto a : allorb_vec) {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          const double ja = oa.j2 * 0.5;
          for (auto b : allorb_vec) {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            if (not AngMom::Triangle(ja, ob.j2 * 0.5, (double)J1))
              continue;
            pairs1.push_back({a, b});
          }
        }
        const int n1 = (int)pairs1.size();
        if (n1 < 1)
          continue;

        arma::vec w_pp(n1), w_hh(n1);
        for (int m = 0; m < n1; ++m) {
          Orbit &oa = Z.modelspace->GetOrbit(pairs1[m][0]);
          Orbit &ob = Z.modelspace->GetOrbit(pairs1[m][1]);
          const double na = oa.occ, nb = ob.occ;
          w_pp(m) = (1.0 - na) * (1.0 - nb);
          w_hh(m) = na * nb;
        }

        arma::mat Om(n0, n1, arma::fill::zeros);
        arma::mat Om_r(n1, n0, arma::fill::zeros);
        for (int n = 0; n < n0; ++n) {
          const index_t i = pairs0[n][0], j = pairs0[n][1];
          for (int m = 0; m < n1; ++m) {
            const index_t a = pairs1[m][0], b = pairs1[m][1];
            Om(n, m) = Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b);
          }
        }
        for (int m = 0; m < n1; ++m) {
          const index_t a = pairs1[m][0], b = pairs1[m][1];
          for (int n = 0; n < n0; ++n) {
            const index_t k = pairs0[n][0], l = pairs0[n][1];
            Om_r(m, n) = Eta.TwoBody.GetTBME_J(J1, J0, a, b, k, l);
          }
        }

        const double ang =
            Z.modelspace->phase(J0 + J1 + lambda) * hat_lambda_inv;
        T_pp += ang * (Om * arma::diagmat(w_pp) * Om_r);
        T_hh += ang * (Om * arma::diagmat(w_hh) * Om_r);
      }

      // Map (c,d) → row index for fast trace of spectator c
      std::map<std::pair<index_t, index_t>, int> pair_index;
      for (int n = 0; n < n0; ++n)
        pair_index[{pairs0[n][0], pairs0[n][1]}] = n;

      for (int n = 0; n < n0; ++n) {
        const index_t c = pairs0[n][0];
        const index_t d = pairs0[n][1];
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        const double n_c = oc.occ, nbar_c = 1.0 - n_c;
        const double n_d = od.occ, nbar_d = 1.0 - n_d;
        const double hatj2_inv = 1.0 / (od.j2 + 1.0);

        for (auto e : allorb_vec) {
          Orbit &oe = Z.modelspace->GetOrbit(e);
          if (oe.j2 != od.j2)
            continue;
          auto it = pair_index.find({c, e});
          if (it == pair_index.end())
            continue;
          const int m = it->second;
          const double chi_c =
              n_c * n_d * T_pp(n, m) - nbar_c * nbar_d * T_hh(n, m);
          Chi_alpha(d, e) += chi_c * hatj2_inv;
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chialpha"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // ------------------------------------------------------------------
    // f^I (scalar final: j_p == j_q)
    //   Z_pq += 1/2/(2jp+1) * sum_deJ (2J+1) Chi_de (Gamma_epdq + Gamma_depq)
    // Loop J over 0..max_J like the unfactored reference (do not share one
    // triangular bound for both Gamma strings — their J windows can differ).
    // ------------------------------------------------------------------
#pragma omp parallel for schedule(dynamic, 1)
    for (int indexp = 0; indexp < norbits; ++indexp) {
      auto p = allorb_vec[indexp];
      Orbit &op = Z.modelspace->GetOrbit(p);
      for (auto q : Z.modelspace->all_orbits) {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if (oq.j2 != op.j2)
          continue;

        double zij = 0.0;
        for (auto d : Z.modelspace->all_orbits) {
          Orbit &od = Z.modelspace->GetOrbit(d);
          for (auto e : Z.modelspace->all_orbits) {
            Orbit &oe = Z.modelspace->GetOrbit(e);
            if (oe.j2 != od.j2)
              continue;
            double chi = Chi_alpha(d, e);
            if (std::abs(chi) < 1e-14)
              continue;

            for (int J2 = 0; J2 <= max_J; ++J2) {
              double g = Gamma.TwoBody.GetTBME_J(J2, J2, e, p, d, q)
                       + Gamma.TwoBody.GetTBME_J(J2, J2, d, p, e, q);
              if (std::abs(g) < 1e-14)
                continue;
              zij += (2 * J2 + 1.0) * chi * g;
            }
          }
        }

        Z.OneBody(p, q) += 0.5 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.5 * hZ * zij / (op.j2 + 1.0);
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fI"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
  }

  // ==================================================================
  // f^II via tensor chi^beta  (Omega*Gamma -> chi^beta; chi^beta*Omega -> f)
  // AMC: learn/amc_tts/factored_fII/output/chi_beta.tex, f2a_from_chi.tex,
  //      f2b_from_chi.tex. Assemble f = f2a + hGamma * f2b.
  // Overall 1/2 lives in chi (do not repeat in the contraction).
  //
  // Convention: Ω tensor reduced; Γ scalar unreduced (st driver); χ^β
  // reduced tensor 1b; Z scalar unreduced (÷ĵ_p²). Equivalent to TTS Case2
  // (MakeReduced Γ, drop Ĵ0) because G_red = Ĵ0 G_unred.
  // ==================================================================
  if (use_TypeII_1b) {
    auto Chi_beta = Z.OneBody;
    Chi_beta.zeros();

    // χ^β: rectangular ladder DGEMM  Γ^{J0} W Ω^{J0 J1 λ}  then 6j + trace c
    //   χ_de^λ = 1/2 Σ_{c J0 J1} (−1)^{je+λ+J0+jc} Ĵ0 Ĵ1
    //            {λ J1 J0; jc jd je} (n̄_c n̄_e T_hh − n_c n_e T_pp)
    //   T_hh = Γ diag(n_a n_b) Ω,  T_pp = Γ diag(n̄_a n̄_b) Ω
    // No ÷ĵ_d²: reduced tensor 1b.
    for (int J0 = 0; J0 <= max_J; ++J0) {
      std::vector<std::array<index_t, 2>> pairs0;
      pairs0.reserve((size_t)norbits * norbits);
      for (auto i : allorb_vec) {
        Orbit &oi = Z.modelspace->GetOrbit(i);
        const double ji = oi.j2 * 0.5;
        for (auto j : allorb_vec) {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          if (not AngMom::Triangle(ji, oj.j2 * 0.5, (double)J0))
            continue;
          pairs0.push_back({i, j});
        }
      }
      const int n0 = (int)pairs0.size();
      if (n0 < 1)
        continue;

      // Γ^{J0} on pairs0 × pairs_ab (ab also at J0)
      arma::mat Gmat(n0, n0, arma::fill::zeros);
      for (int n = 0; n < n0; ++n) {
        const index_t c = pairs0[n][0], d = pairs0[n][1];
        for (int m = 0; m < n0; ++m) {
          const index_t a = pairs0[m][0], b = pairs0[m][1];
          Gmat(n, m) = Gamma.TwoBody.GetTBME_J(J0, J0, c, d, a, b);
        }
      }

      for (int J1 = 0; J1 <= max_J; ++J1) {
        if (not AngMom::Triangle(J0, J1, lambda))
          continue;

        std::vector<std::array<index_t, 2>> pairs1;
        pairs1.reserve((size_t)norbits * norbits);
        for (auto i : allorb_vec) {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          const double ji = oi.j2 * 0.5;
          for (auto j : allorb_vec) {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            if (not AngMom::Triangle(ji, oj.j2 * 0.5, (double)J1))
              continue;
            pairs1.push_back({i, j});
          }
        }
        const int n1 = (int)pairs1.size();
        if (n1 < 1)
          continue;

        // Mid (a,b) at J0 — same pair list as pairs0
        arma::vec w_hh(n0), w_pp(n0);
        for (int m = 0; m < n0; ++m) {
          Orbit &oa = Z.modelspace->GetOrbit(pairs0[m][0]);
          Orbit &ob = Z.modelspace->GetOrbit(pairs0[m][1]);
          const double na = oa.occ, nb = ob.occ;
          w_hh(m) = na * nb;
          w_pp(m) = (1.0 - na) * (1.0 - nb);
        }

        arma::mat Om(n0, n1, arma::fill::zeros); // (a,b)_{J0} × (c,e)_{J1}
        for (int m = 0; m < n0; ++m) {
          const index_t a = pairs0[m][0], b = pairs0[m][1];
          for (int k = 0; k < n1; ++k) {
            const index_t c = pairs1[k][0], e = pairs1[k][1];
            Om(m, k) = Eta.TwoBody.GetTBME_J(J0, J1, a, b, c, e);
          }
        }

        // T = Γ W Ω : rows (c,d)@J0, cols (c',e)@J1
        const double hats01 = hat(J0) * hat(J1);
        arma::mat T_hh = hats01 * (Gmat * arma::diagmat(w_hh) * Om);
        arma::mat T_pp = hats01 * (Gmat * arma::diagmat(w_pp) * Om);

        std::map<std::pair<index_t, index_t>, int> idx1;
        for (int k = 0; k < n1; ++k)
          idx1[{pairs1[k][0], pairs1[k][1]}] = k;

        for (int n = 0; n < n0; ++n) {
          const index_t c = pairs0[n][0];
          const index_t d = pairs0[n][1];
          Orbit &oc = Z.modelspace->GetOrbit(c);
          Orbit &od = Z.modelspace->GetOrbit(d);
          const double jc = oc.j2 / 2.0, jd = od.j2 / 2.0;
          const double n_c = oc.occ, nbar_c = 1.0 - n_c;

          for (auto e : allorb_vec) {
            Orbit &oe = Z.modelspace->GetOrbit(e);
            const double je = oe.j2 / 2.0;
            if (not AngMom::Triangle(jd, je, (double)lambda))
              continue;
            auto it = idx1.find({c, e});
            if (it == idx1.end())
              continue;
            const int k = it->second;
            const double n_e = oe.occ, nbar_e = 1.0 - n_e;
            // occ: n_a n_b n̄_c n̄_e − n̄_a n̄_b n_c n_e
            const double block =
                nbar_c * nbar_e * T_hh(n, k) - n_c * n_e * T_pp(n, k);
            if (std::abs(block) < 1e-14)
              continue;
            const double sixj = AngMom::SixJ(lambda, J1, J0, jc, jd, je);
            if (std::abs(sixj) < 1e-14)
              continue;
            const double phase =
                Z.modelspace->phase((oe.j2 + oc.j2) / 2 + lambda + J0);
            Chi_beta(d, e) += 0.5 * phase * sixj * block;
          }
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chibeta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // f^II (scalar final): tensor χ^β × tensor Ω → scalar via RME [χ⊗Ω]^(0)
    //   Z_pq += 1/(2jp+1) * sum_deJ3J4 (−1)^{jp+J4+jd} * Ĵ3 Ĵ4 / λ̂
    //           * SixJ(J4,λ,J3; je,jp,jd) * Chi_de
    //           * (Ω_eidj^{J3 J4 λ} + hΓ Ω_ejdi^{J3 J4 λ})
#pragma omp parallel for schedule(dynamic, 1)
    for (int indexp = 0; indexp < norbits; ++indexp) {
      auto p = allorb_vec[indexp];
      Orbit &op = Z.modelspace->GetOrbit(p);
      double jp = op.j2 / 2.0;
      for (auto q : Z.modelspace->all_orbits) {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if (oq.j2 != op.j2)
          continue;

        double zij = 0.0;
        for (auto d : Z.modelspace->all_orbits) {
          Orbit &od = Z.modelspace->GetOrbit(d);
          double jd = od.j2 / 2.0;
          for (auto e : Z.modelspace->all_orbits) {
            Orbit &oe = Z.modelspace->GetOrbit(e);
            double je = oe.j2 / 2.0;
            if (not AngMom::Triangle(jd, je, (double)lambda))
              continue;
            double chi = Chi_beta(d, e);
            if (std::abs(chi) < 1e-14)
              continue;

            for (int J3 = 0; J3 <= max_J; ++J3)
            for (int J4 = 0; J4 <= max_J; ++J4) {
              if (not AngMom::Triangle(J3, J4, lambda))
                continue;
              double sixj = AngMom::SixJ(J4, lambda, J3, je, jp, jd);
              if (std::abs(sixj) < 1e-14)
                continue;
              double eta2a = Eta.TwoBody.GetTBME_J(J3, J4, e, p, d, q);
              double eta2b = Eta.TwoBody.GetTBME_J(J3, J4, e, q, d, p);
              double omega = eta2a + hGamma * eta2b;
              if (std::abs(omega) < 1e-14)
                continue;
              // AMC f2a: combined (−1)^{j_i+J_4+j_d} ≡ phase((jp2+jd2)/2+J4).
              // Do NOT multiply an extra (−1)^{j_i}: that double-counts and was
              // only a J↔J tune vs tts_fII (see factorized_code_analyze.tex).
              double phase = Z.modelspace->phase((op.j2 + od.j2) / 2 + J4);
              zij += phase * hat(J3) * hat(J4) * hat_lambda_inv * sixj * chi *
                     omega;
            }
          }
        }

        Z.OneBody(p, q) += zij / (op.j2 + 1.0); // unreduced scalar 1b
        if (p != q)
          Z.OneBody(q, p) += hZ * zij / (op.j2 + 1.0);
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fII"] += omp_get_wtime() - t_internal;
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// ethS two-body intermediates for one-body final:
///   chi^delta (scalar 2b from tensor×tensor Omega*Omega) -> f^III_b
///   chi^gamma (scalar 2b, Pandya topology; AMC — not Neithan) -> f^III_a
///
/// χ^γ gold: m ≡ AMC direct ≡ AMC Path B (corrected inv). Do NOT use
/// neithan.tex TT→0 here — it disagrees with m/AMC by ~λ̂.
/// Docs: learn/amc_tts/factored_fIIIa/OMEGA_TT_TO_SCALAR.md
///
/// Both chi^delta and chi^gamma are NON-HERMITIAN (occ-weighted). Store full
/// chi(ch,ch) squares only — no hermiticity fill / no ch_bra>ch_ket conjugate.
///   f^III_b: M = chi*Gamma; same-ch M+=M.t(); unequal M-=Gamma*chi_ket
///   f^III_a: M = hat(J)^2 (Gamma*chi - chi*Gamma)  [no .t() shortcut]
////////////////////////////////////////////////////////////////////////////
/// Leftover χ^γ × Γ ladder (any λ). Same W1·W2 as tts_fIIIa, (a,b) by DGEMM.
static void comm223_231_fIIIa_leftover_dgemm(const Operator &Eta,
                                            const Operator &Gamma,
                                            Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  const int lambda = Eta.GetJRank();
  const double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  const int hZ = Gamma.IsHermitian() ? 1 : -1;
  std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                              Z.modelspace->all_orbits.end());
  const int n_orb = (int)allorb.size();
  int max_j2 = 0;
  int n_orb_alloc = 0;
  for (auto o : allorb) {
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(o).j2);
    n_orb_alloc = std::max(n_orb_alloc, (int)o + 1);
  }
  const int max_J = max_j2;
  const int max_j0_2 = 2 * (max_J + lambda + 2);
  const int n_J = max_J + 1;
  const int n_j0 = max_j0_2 + 1;

  std::vector<int> compact(n_orb_alloc, -1);
  std::vector<int> oj2(n_orb);
  std::vector<double> ojh(n_orb), nocc(n_orb), nnocc(n_orb);
  for (int t = 0; t < n_orb; ++t) {
    compact[(int)allorb[t]] = t;
    Orbit &o = Z.modelspace->GetOrbit(allorb[t]);
    oj2[t] = o.j2;
    ojh[t] = o.j2 * 0.5;
    nocc[t] = o.occ;
    nnocc[t] = 1.0 - o.occ;
  }

  const size_t w_s_j0 = 1;
  const size_t w_s_J6 = (size_t)n_j0 * w_s_j0;
  const size_t w_s_b = (size_t)n_J * w_s_J6;
  const size_t w_s_k = (size_t)n_orb * w_s_b;
  const size_t w_s_j = (size_t)n_orb * w_s_k;
  const size_t w_s_a = (size_t)n_orb * w_s_j;
  std::vector<double> W1((size_t)n_orb * w_s_a, 0.0),
      W2((size_t)n_orb * w_s_a, 0.0);
  auto w_index = [&](int a, int j, int k, int b, int J6, int j02) -> size_t {
    return (size_t)a * w_s_a + (size_t)j * w_s_j + (size_t)k * w_s_k +
           (size_t)b * w_s_b + (size_t)J6 * w_s_J6 + (size_t)j02 * w_s_j0;
  };

  double t_w1 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ia = 0; ia < n_orb; ++ia) {
    const int a = (int)allorb[ia];
    const double ja = ojh[ia];
    const int ja2 = oj2[ia];
    std::vector<double> s1((size_t)n_j0, 0.0);
    for (int ij = 0; ij < n_orb; ++ij) {
      const int j = (int)allorb[ij];
      const double jj = ojh[ij];
      const int jj2 = oj2[ij];
      const int J2min = std::abs(ja2 - jj2) / 2;
      const int J2max = (ja2 + jj2) / 2;
      for (int J2 = J2min; J2 <= J2max; ++J2) {
        for (int J3 = 0; J3 <= max_J; ++J3) {
          if (not AngMom::Triangle(J2, J3, lambda))
            continue;
          int s1lo = std::max(std::abs(2 * J3 - jj2), std::abs(ja2 - 2 * lambda));
          int s1hi = std::min({max_j0_2, 2 * J3 + jj2, ja2 + 2 * lambda});
          if ((s1lo & 1) == 0)
            ++s1lo;
          if (s1lo > s1hi)
            continue;
          std::fill(s1.begin(), s1.end(), 0.0);
          bool any_s1 = false;
          for (int j02 = s1lo; j02 <= s1hi; j02 += 2) {
            const double sj =
                AngMom::SixJ(J3, lambda, J2, ja, jj, 0.5 * j02);
            if (std::abs(sj) < 1e-14)
              continue;
            s1[(size_t)j02] = sj;
            any_s1 = true;
          }
          if (not any_s1)
            continue;
          const double hats23 =
              std::sqrt((2.0 * J2 + 1.0) * (2.0 * J3 + 1.0));
          const double phase3 = Z.modelspace->phase(J3);
          for (int ik = 0; ik < n_orb; ++ik) {
            const int k = (int)allorb[ik];
            const double jk = ojh[ik];
            const int jk2 = oj2[ik];
            const int J6min = std::abs(jj2 - jk2) / 2;
            const int J6max = (jj2 + jk2) / 2;
            for (int ib = 0; ib < n_orb; ++ib) {
              const int jb2 = oj2[ib];
              const int J3min = std::abs(jk2 - jb2) / 2;
              const int J3max = (jk2 + jb2) / 2;
              if (J3 < J3min or J3 > J3max)
                continue;
              const int b = (int)allorb[ib];
              const double o1 = Eta.TwoBody.GetTBME_J(J2, J3, a, j, k, b);
              if (std::abs(o1) < 1e-14)
                continue;
              const double jb = ojh[ib];
              const double pref = phase3 * hats23 * o1;
              for (int J6 = J6min; J6 <= J6max; ++J6) {
                int lo = std::max(s1lo, std::abs(jb2 - 2 * J6));
                int hi = std::min(s1hi, jb2 + 2 * J6);
                if ((lo & 1) == 0)
                  ++lo;
                for (int j02 = lo; j02 <= hi; j02 += 2) {
                  const double s1v = s1[(size_t)j02];
                  if (std::abs(s1v) < 1e-14)
                    continue;
                  const double s2 =
                      AngMom::SixJ(jj, jk, J6, jb, 0.5 * j02, J3);
                  if (std::abs(s2) < 1e-14)
                    continue;
                  W1[w_index(ia, ij, ik, ib, J6, j02)] += pref * s1v * s2;
                }
              }
            }
          }
        }
      }
    }
  }
  Z.profiler.timer["_fIIIa_W1"] += omp_get_wtime() - t_w1;

  double t_w2 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ia = 0; ia < n_orb; ++ia) {
    const int a = (int)allorb[ia];
    const double ja = ojh[ia];
    const int ja2 = oj2[ia];
    std::vector<double> s1((size_t)n_j0, 0.0);
    for (int il = 0; il < n_orb; ++il) {
      const int l = (int)allorb[il];
      const double jl = ojh[il];
      const int jl2 = oj2[il];
      const int J5min_al = std::abs(ja2 - jl2) / 2;
      const int J5max_al = (ja2 + jl2) / 2;
      for (int J4 = 0; J4 <= max_J; ++J4) {
        for (int J5 = J5min_al; J5 <= J5max_al; ++J5) {
          if (not AngMom::Triangle(J4, J5, lambda))
            continue;
          int s1lo = std::max(std::abs(2 * J4 - jl2), std::abs(ja2 - 2 * lambda));
          int s1hi = std::min({max_j0_2, 2 * J4 + jl2, ja2 + 2 * lambda});
          if ((s1lo & 1) == 0)
            ++s1lo;
          if (s1lo > s1hi)
            continue;
          std::fill(s1.begin(), s1.end(), 0.0);
          bool any_s1 = false;
          for (int j02 = s1lo; j02 <= s1hi; j02 += 2) {
            const double sj =
                AngMom::SixJ(J4, lambda, J5, ja, jl, 0.5 * j02);
            if (std::abs(sj) < 1e-14)
              continue;
            s1[(size_t)j02] = sj;
            any_s1 = true;
          }
          if (not any_s1)
            continue;
          const double hats45 =
              std::sqrt((2.0 * J4 + 1.0) * (2.0 * J5 + 1.0));
          const double phase5 = Z.modelspace->phase(J5);
          for (int ii = 0; ii < n_orb; ++ii) {
            const int i = (int)allorb[ii];
            const double ji = ojh[ii];
            const int ji2 = oj2[ii];
            const int J6min = std::abs(jl2 - ji2) / 2;
            const int J6max = (jl2 + ji2) / 2;
            for (int ib = 0; ib < n_orb; ++ib) {
              const int jb2 = oj2[ib];
              const int J4min = std::abs(ji2 - jb2) / 2;
              const int J4max = (ji2 + jb2) / 2;
              if (J4 < J4min or J4 > J4max)
                continue;
              const int b = (int)allorb[ib];
              const double o2 = Eta.TwoBody.GetTBME_J(J4, J5, i, b, a, l);
              if (std::abs(o2) < 1e-14)
                continue;
              const double jb = ojh[ib];
              const double pref = phase5 * hats45 * o2;
              for (int J6 = J6min; J6 <= J6max; ++J6) {
                int lo = std::max(s1lo, std::abs(jb2 - 2 * J6));
                int hi = std::min(s1hi, jb2 + 2 * J6);
                if ((lo & 1) == 0)
                  ++lo;
                for (int j02 = lo; j02 <= hi; j02 += 2) {
                  const double s1v = s1[(size_t)j02];
                  if (std::abs(s1v) < 1e-14)
                    continue;
                  const double s2 =
                      AngMom::SixJ(jl, ji, J6, jb, 0.5 * j02, J4);
                  if (std::abs(s2) < 1e-14)
                    continue;
                  W2[w_index(ii, ib, ia, il, J6, j02)] += pref * s1v * s2;
                }
              }
            }
          }
        }
      }
    }
  }
  Z.profiler.timer["_fIIIa_W2"] += omp_get_wtime() - t_w2;

  const size_t chi_stride_l = (size_t)n_J;
  const size_t chi_stride_k = (size_t)n_orb * chi_stride_l;
  const size_t chi_stride_j = (size_t)n_orb * chi_stride_k;
  const size_t chi_stride_i = (size_t)n_orb * chi_stride_j;
  std::vector<double> chi_tab((size_t)n_orb * chi_stride_i, 0.0);
  auto chi_index = [&](int i, int j, int k, int l, int J0) -> size_t {
    return (size_t)i * chi_stride_i + (size_t)j * chi_stride_j +
           (size_t)k * chi_stride_k + (size_t)l * chi_stride_l + (size_t)J0;
  };

  const double phase_lambda = Z.modelspace->phase(lambda);
  const int n2 = n_orb * n_orb;
  double t_dgemm = omp_get_wtime();
  if (n2 > 0) {
    arma::mat X(n2, n2), Y(n2, n2);
    for (int J6 = 0; J6 <= max_J; ++J6) {
      for (int j02 = 1; j02 <= max_j0_2; j02 += 2) {
        X.zeros();
        Y.zeros();
#pragma omp parallel for schedule(dynamic, 1)
        for (int col = 0; col < n2; ++col) {
          const int ia = col / n_orb;
          const int ib = col - ia * n_orb;
          const double occ_ab_1 = nocc[ia] * nnocc[ib];
          const double occ_ab_2 = nnocc[ia] * nocc[ib];
          for (int ij = 0; ij < n_orb; ++ij) {
            for (int ik = 0; ik < n_orb; ++ik) {
              const double w1 = W1[w_index(ia, ij, ik, ib, J6, j02)];
              if (std::abs(w1) < 1e-16)
                continue;
              const double occ = nocc[ij] * nnocc[ik] * occ_ab_1 -
                                 nnocc[ij] * nocc[ik] * occ_ab_2;
              X(ij * n_orb + ik, col) = occ * w1;
            }
          }
          for (int ii = 0; ii < n_orb; ++ii) {
            for (int il = 0; il < n_orb; ++il) {
              Y(ii * n_orb + il, col) = W2[w_index(ii, ib, ia, il, J6, j02)];
            }
          }
        }
        const arma::mat M = X * Y.t();
        const double sc = -phase_lambda * hat_lambda_inv * (2.0 * J6 + 1.0) *
                          (j02 + 1.0);
#pragma omp parallel for schedule(dynamic, 1)
        for (int ii = 0; ii < n_orb; ++ii) {
          const double ji = ojh[ii];
          for (int ij = 0; ij < n_orb; ++ij) {
            const double jj = ojh[ij];
            const int jj2 = oj2[ij];
            for (int ik = 0; ik < n_orb; ++ik) {
              const double jk = ojh[ik];
              for (int il = 0; il < n_orb; ++il) {
                const double jl = ojh[il];
                const double ph_jl =
                    Z.modelspace->phase((jj2 + oj2[il]) / 2);
                const double m = M(ij * n_orb + ik, ii * n_orb + il);
                if (std::abs(m) < 1e-16)
                  continue;
                for (int J0 = 0; J0 <= max_J; ++J0) {
                  if (not AngMom::Triangle(ji, jj, (double)J0) or
                      not AngMom::Triangle(jk, jl, (double)J0))
                    continue;
                  const double sixj =
                      AngMom::SixJ(jk, jl, J0, ji, jj, J6);
                  if (std::abs(sixj) < 1e-14)
                    continue;
                  chi_tab[chi_index(ii, ij, ik, il, J0)] +=
                      sc * ph_jl * sixj * m;
                }
              }
            }
          }
        }
      }
    }
  }
  Z.profiler.timer["_fIIIa_dgemm"] += omp_get_wtime() - t_dgemm;
  W1.clear();
  W2.clear();

  double t_fold = omp_get_wtime();
  std::vector<double> gam((size_t)n_orb * chi_stride_i, 0.0);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const int i = (int)allorb[ii];
    const double ji = ojh[ii];
    for (int ij = 0; ij < n_orb; ++ij) {
      const int j = (int)allorb[ij];
      const double jj = ojh[ij];
      for (int ik = 0; ik < n_orb; ++ik) {
        const int k = (int)allorb[ik];
        const double jk = ojh[ik];
        for (int il = 0; il < n_orb; ++il) {
          const int l = (int)allorb[il];
          const double jl = ojh[il];
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(ji, jj, (double)J0) or
                not AngMom::Triangle(jk, jl, (double)J0))
              continue;
            gam[chi_index(ii, ij, ik, il, J0)] =
                Gamma.TwoBody.GetTBME_J(J0, J0, i, j, k, l);
          }
        }
      }
    }
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < n_orb; ++ip) {
    const int p = (int)allorb[ip];
    const int jp2 = oj2[ip];
    for (int iq = 0; iq <= ip; ++iq) {
      if (oj2[iq] != jp2)
        continue;
      const int q = (int)allorb[iq];
      double zpq = 0.0;
      for (int ic = 0; ic < n_orb; ++ic) {
        for (int ia = 0; ia < n_orb; ++ia) {
          for (int ib = 0; ib < n_orb; ++ib) {
            for (int J0 = 0; J0 <= max_J; ++J0) {
              const double g1 = gam[chi_index(ic, ip, ia, ib, J0)];
              const double g2 = gam[chi_index(ia, ib, iq, ic, J0)];
              if (std::abs(g1) < 1e-14 && std::abs(g2) < 1e-14)
                continue;
              zpq += (2.0 * J0 + 1.0) *
                     (g1 * chi_tab[chi_index(ia, ib, ic, iq, J0)] -
                      chi_tab[chi_index(ip, ic, ia, ib, J0)] * g2);
            }
          }
        }
      }
      Z.OneBody(p, q) += zpq / (jp2 + 1.0);
      if (p != q)
        Z.OneBody(q, p) += hZ * zpq / (jp2 + 1.0);
    }
  }
  Z.profiler.timer["_fIIIa_fold"] += omp_get_wtime() - t_fold;
  Z.profiler.timer["comm223_231_tts_fIIIa"] += omp_get_wtime() - t_start;
}

void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {

  double t_internal = omp_get_wtime();
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();

  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);

  int norbits = Z.modelspace->all_orbits.size();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),
                                  Z.modelspace->all_orbits.end());

  int max_j2 = 0;
  for (auto x : Z.modelspace->all_orbits)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
  int max_J = max_j2;

  // ==================================================================
  // f^III_b via scalar chi^delta  (use_TypeIII_1b)
  // Case 2: Ω tensor reduced; Γ scalar unreduced here ≡ G_red/Ĵ from TTS;
  // Z stored unreduced (÷ĵ²).  χ^δ: RME [Ω⊗Ω]^(0) via DGEMM on each J channel.
  // ==================================================================
  if (use_TypeIII_1b) {
    // ------------------------------------------------------------------
    // χ^δ: non-Hermitian scalar, FULL chi(ch,ch) only.
    //   T_pp/hh = Σ_{J2} (-1)^{J+J2+λ} λ̂^{-1}  Ω^{J J2} W_{pp/hh} Ω^{J2 J}
    //   χ_phys = n̄_k n̄_l T_hh − n_k n_l T_pp
    // No 1/4 in χ (4×0.25 with Γ). Store NORMALIZED GetMatrix entries.
    // ------------------------------------------------------------------
    TwoBodyME Chi_delta = Z.TwoBody;
    Chi_delta.Erase();
    Chi_delta.SetNonHermitian();

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch; ++ch) {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      if (nKets < 1)
        continue;
      arma::mat &Chi_mat = Chi_delta.GetMatrix(ch, ch);

      arma::mat T_pp(nKets, nKets, arma::fill::zeros);
      arma::mat T_hh(nKets, nKets, arma::fill::zeros);

      for (int J2 = 0; J2 <= max_J; ++J2) {
        if (not AngMom::Triangle(J0, J2, lambda))
          continue;

        std::vector<std::array<index_t, 2>> pairs2;
        pairs2.reserve((size_t)norbits * norbits);
        for (auto a : allorb_vec) {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          const double ja = oa.j2 * 0.5;
          for (auto b : allorb_vec) {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            if (not AngMom::Triangle(ja, ob.j2 * 0.5, (double)J2))
              continue;
            pairs2.push_back({a, b});
          }
        }
        const int n2 = (int)pairs2.size();
        if (n2 < 1)
          continue;

        arma::vec w_pp(n2), w_hh(n2);
        for (int m = 0; m < n2; ++m) {
          Orbit &oa = Z.modelspace->GetOrbit(pairs2[m][0]);
          Orbit &ob = Z.modelspace->GetOrbit(pairs2[m][1]);
          const double na = oa.occ, nb = ob.occ;
          w_pp(m) = (1.0 - na) * (1.0 - nb);
          w_hh(m) = na * nb;
        }

        arma::mat Om(nKets, n2, arma::fill::zeros);
        arma::mat Om_r(n2, nKets, arma::fill::zeros);
        for (int ibra = 0; ibra < nKets; ++ibra) {
          Ket &bra = tbc.GetKet(ibra);
          const index_t i = bra.p, j = bra.q;
          for (int m = 0; m < n2; ++m) {
            const index_t a = pairs2[m][0], b = pairs2[m][1];
            Om(ibra, m) = Eta.TwoBody.GetTBME_J(J0, J2, i, j, a, b);
          }
        }
        for (int m = 0; m < n2; ++m) {
          const index_t a = pairs2[m][0], b = pairs2[m][1];
          for (int iket = 0; iket < nKets; ++iket) {
            Ket &ket = tbc.GetKet(iket);
            Om_r(m, iket) = Eta.TwoBody.GetTBME_J(J2, J0, a, b, ket.p, ket.q);
          }
        }

        const double ang =
            Z.modelspace->phase(J0 + J2 + lambda) * hat_lambda_inv;
        T_pp += ang * (Om * arma::diagmat(w_pp) * Om_r);
        T_hh += ang * (Om * arma::diagmat(w_hh) * Om_r);
      }

      for (int ibra = 0; ibra < nKets; ++ibra) {
        Ket &bra = tbc.GetKet(ibra);
        const index_t i = bra.p, j = bra.q;
        for (int iket = 0; iket < nKets; ++iket) {
          Ket &ket = tbc.GetKet(iket);
          const index_t k = ket.p, l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          const double n_k = ok.occ, nbar_k = 1.0 - n_k;
          const double n_l = ol.occ, nbar_l = 1.0 - n_l;
          // occ = n_a n_b n̄_k n̄_l − n̄_a n̄_b n_k n_l
          const double chi_phys =
              nbar_k * nbar_l * T_hh(ibra, iket) - n_k * n_l * T_pp(ibra, iket);
          double nrm = 1.0;
          if (i == j)
            nrm *= PhysConst::SQRT2;
          if (k == l)
            nrm *= PhysConst::SQRT2;
          Chi_mat(ibra, iket) = chi_phys / nrm;
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chidelta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // ------------------------------------------------------------------
    // Fold χ^δ × Γ (normalized DGEMM). Case-2 TTS uses G_red/Ĵ ≡ G_unred.
    // No extra hat(J)^2 (not scalar's 4*(2J+1)). M*=2 then 0.25 on 1b.
    // ------------------------------------------------------------------
    TwoBodyME intermediateTB = Z.TwoBody;
    intermediateTB.Erase();

    std::vector<int> bra_channels;
    std::vector<int> ket_channels;
    for (auto &itmat : Z.TwoBody.MatEl) {
      bra_channels.push_back(itmat.first[0]);
      ket_channels.push_back(itmat.first[1]);
    }
    int nbra_ket_ch = bra_channels.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nbra_ket_ch; ich++) {
      size_t ch_bra = bra_channels[ich];
      size_t ch_ket = ket_channels[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      if (tbc_bra.J != tbc_ket.J)
        continue;

      const arma::mat &Chi_bra = Chi_delta.GetMatrix(ch_bra, ch_bra);
      const arma::mat &Chi_ket = Chi_delta.GetMatrix(ch_ket, ch_ket);
      const arma::mat &Gamma_mat = Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);

      // Normalized-space product (like scalar Eta * Eta_nnnn * Gamma)
      arma::mat M = Chi_bra * Gamma_mat;
      if (ch_bra == ch_ket) {
        // Square block: A + A^T  (h_Omega=-1; analyze §code-fIIIb)
        M += M.t();
      } else {
        // Rectangular: explicit partner — do NOT .t() into this MatEl key
        M -= Gamma_mat * Chi_ket;
      }
      // 1/4 from chi^delta: use *=2 here because A+A^T already supplies
      // the second topology; scalar uses *=4*(2J+1) with a different ME
      // convention (same-J Omega only). TTS unfactored has no net hat(J)^2.
      M *= 2.0;
      intermediateTB.GetMatrix(ch_bra, ch_ket) = M;
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fIIIb_chiG"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // One-body reduction: sum_c M_{cp,cq}  (same as scalar)
#pragma omp parallel for schedule(dynamic, 1)
    for (int indexp = 0; indexp < norbits; ++indexp) {
      auto p = allorb_vec[indexp];
      Orbit &op = Z.modelspace->GetOrbit(p);
      for (auto q : Z.modelspace->all_orbits) {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if (oq.j2 != op.j2)
          continue;

        double zpq = 0.0;
        for (auto &c : Z.modelspace->all_orbits) {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          int J0min = std::abs(oc.j2 - op.j2) / 2;
          int J0max = (oc.j2 + op.j2) / 2;
          for (int J0 = J0min; J0 <= J0max; ++J0) {
            zpq += intermediateTB.GetTBME_J(J0, J0, c, p, c, q);
          }
        }
        Z.OneBody(p, q) += 0.25 * zpq / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.25 * hZ * zpq / (op.j2 + 1.0);
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fIIIb"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
  }


  // ==================================================================
  // f^III_a via χ^γ (scalar, NON-HERMITIAN)  (use_TypeIIIa_1b)
  // Leftover χ^γ × Γ is AMC Path B / leftover_dgemm for any λ
  // (reduced Ω; hats and Triangle(J,J′,λ) already cover λ=0).
  // ==================================================================
  if (use_TypeIIIa_1b) {
    comm223_231_fIIIa_leftover_dgemm(Eta, Gamma, Z);
    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fIIIa"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
  }


  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}


void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  // Ω is reduced RME (including λ=0). Gamma and Z are scalar unreduced.
  if (use_1b_intermediates) {
    comm223_232_chi1b_tensor(Eta, Gamma,
                             Z); // topology with 1-body intermediate (fast)
  }
  if (use_2b_intermediates) {
    comm223_232_chi2b(Eta, Gamma,
                      Z); // topology with 2-body intermediate (slow)
  }
}

////////////////////////////////////////////////////////////////////////////
/// factorized 223_232 double commutator with 1b intermediate
///
/// Gamma^I  : chi^epsilon (Omega x Omega -> scalar 1b) x Gamma  -- any lambda
/// Gamma^II : bra χ^{ΩΓ} (Ω×Γ → tensor 1b) × Ω, ket χ^ζ (Γ×Ω) × Ω
///            Unfact bra is ΩΓΩ; χ^ζ on the bra matches only at λ=0.
///            Docs: learn/amc_tts/factored_GII/NOTES.md
////////////////////////////////////////////////////////////////////////////
void comm223_232_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {
  double t_start = omp_get_wtime();
  double t_internal = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  size_t nch = ch_bra_list.size();
  auto &Z2 = Z.TwoBody;

  bool Z_is_scalar = (Z.TwoBody.rank_T == 0);
  int hEta = Eta.IsHermitian() ? 1 : -1;
  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  (void)hEta;
  (void)hZ;
  bool tensor_case = (Gamma.GetJRank() == 0) && Z_is_scalar;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv =
      1.0 / std::sqrt(2.0 * std::max(lambda, 0) + 1.0);

  bool do_GI = use_TypeGI_2b;
  bool do_GII = use_TypeGII_2b;
  if (not do_GI and not do_GII)
    return;

  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };

  // ######################################################################
  // Gamma^I -- chi^epsilon (scalar 1b) x Gamma  (ordinary-channel DGEMM)
  // ######################################################################
  if (do_GI) {
    arma::mat CHI_I = Gamma.OneBody * 0;

    int nch_tb = Z.modelspace->GetNumberTwoBodyChannels();

#pragma omp parallel
    {
      arma::mat CHI_I_loc(CHI_I.n_rows, CHI_I.n_cols, arma::fill::zeros);

#pragma omp for schedule(dynamic, 1)
      for (int ch0 = 0; ch0 < nch_tb; ++ch0) {
        TwoBodyChannel &tbc0 = Z.modelspace->GetTwoBodyChannel(ch0);
        int J0 = tbc0.J;
        int n0 = tbc0.GetNumberKets();
        if (n0 == 0)
          continue;

        int ch1_lo = tensor_case ? 0 : ch0;
        int ch1_hi = tensor_case ? nch_tb : ch0 + 1;
        for (int ch1 = ch1_lo; ch1 < ch1_hi; ++ch1) {
          TwoBodyChannel &tbc1 = Z.modelspace->GetTwoBodyChannel(ch1);
          int J1 = tbc1.J;
          int n1 = tbc1.GetNumberKets();
          if (n1 == 0)
            continue;
          if (tensor_case and not AngMom::Triangle(J0, J1, lambda))
            continue;

          arma::vec w_pp(n1, arma::fill::zeros);
          arma::vec w_hh(n1, arma::fill::zeros);
          for (int iab = 0; iab < n1; ++iab) {
            Ket &kab = tbc1.GetKet(iab);
            double na = kab.op->occ;
            double nb = kab.oq->occ;
            w_pp(iab) = (1.0 - na) * (1.0 - nb);
            w_hh(iab) = na * nb;
          }

          double ang =
              tensor_case
                  ? (Z.modelspace->phase(J0 + J1 + lambda) * hat_lambda_inv)
                  : (2.0 * J0 + 1.0);
          double pref = 2.0 * ang;

          auto fill_norm = [&](const Operator &Op, int ch_bra, int ch_ket,
                               TwoBodyChannel &tb_bra, TwoBodyChannel &tb_ket,
                               int nbra, int nket) {
            arma::mat O(nbra, nket, arma::fill::zeros);
            auto &MatEl = Op.TwoBody.MatEl;
            size_t clo = (size_t)std::min(ch_bra, ch_ket);
            size_t chi = (size_t)std::max(ch_bra, ch_ket);
            auto it = MatEl.find({clo, chi});
            if (it == MatEl.end())
              return O;
            if (ch_bra <= ch_ket) {
              O = it->second;
            } else {
              double jph = Z.modelspace->phase(tb_bra.J - tb_ket.J);
              if (Op.IsHermitian())
                O = jph * it->second.t();
              else if (Op.IsAntiHermitian())
                O = -jph * it->second.t();
              else
                O = jph * it->second.t();
            }
            return O;
          };

          auto has_matel = [](const Operator &Op, int cha, int chb) {
            size_t clo = (size_t)std::min(cha, chb);
            size_t chi = (size_t)std::max(cha, chb);
            return Op.TwoBody.MatEl.find({clo, chi}) != Op.TwoBody.MatEl.end();
          };
          if (not has_matel(Eta, ch0, ch1))
            continue;

          arma::mat O10 = fill_norm(Eta, ch1, ch0, tbc1, tbc0, n1, n0);
          arma::mat O01 = fill_norm(Eta, ch0, ch1, tbc0, tbc1, n0, n1);

          arma::mat T_pp = pref * O01 * arma::diagmat(w_pp) * O10;
          arma::mat T_hh = pref * O01 * arma::diagmat(w_hh) * O10;
          for (auto c : Z.modelspace->all_orbits) {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;
            for (auto p : Z.modelspace->all_orbits) {
              Orbit &op = Z.modelspace->GetOrbit(p);
              int icp = static_cast<int>(
                  tbc0.GetLocalIndex(std::min(c, p), std::max(c, p)));
              if (icp < 0 or icp >= n0)
                continue;
              Ket &kcp = tbc0.GetKet(icp);
              if (kcp.p != std::min(c, p) or kcp.q != std::max(c, p))
                continue;
              double ph_cp = (c > p) ? kcp.Phase(J0) : 1.0;
              double N_cp = (c == p) ? PhysConst::SQRT2 : 1.0;
              for (auto q : Z.modelspace->all_orbits) {
                Orbit &oq = Z.modelspace->GetOrbit(q);
                if (oq.j2 != op.j2)
                  continue;
                int icq = static_cast<int>(
                    tbc0.GetLocalIndex(std::min(c, q), std::max(c, q)));
                if (icq < 0 or icq >= n0)
                  continue;
                Ket &kcq = tbc0.GetKet(icq);
                if (kcq.p != std::min(c, q) or kcq.q != std::max(c, q))
                  continue;
                double ph_cq = (c > q) ? kcq.Phase(J0) : 1.0;
                double N_cq = (c == q) ? PhysConst::SQRT2 : 1.0;
                double me = n_c * T_pp(icp, icq) + nbar_c * T_hh(icp, icq);
                if (std::abs(me) < 1e-16)
                  continue;
                double hatj2 =
                    tensor_case ? (op.j2 + 1.0) : (oq.j2 + 1.0);
                CHI_I_loc(p, q) +=
                    0.5 / hatj2 * ph_cp * ph_cq * N_cp * N_cq * me;
              }
            }
          }
        }
      }
#pragma omp critical
      { CHI_I += CHI_I_loc; }
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nch; ich++) {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras; ibra++) {
        Ket &bra = tbc_bra.GetKet(ibra);
        index_t p = bra.p;
        index_t q = bra.q;
        Orbit &op = Z.modelspace->GetOrbit(p);
        Orbit &oq = Z.modelspace->GetOrbit(q);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (size_t iket = ketmin; iket < nkets; iket++) {
          Ket &ket = tbc_ket.GetKet(iket);
          index_t r = ket.p;
          index_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double zpqrs = 0;

          // Leftover Z is scalar: χ^ε (T×T→S) × Γ is the scalar (1−P) fold.
          // Loop Z/Gamma 1b channels (same l,j,tz), not Eta's tensor channels:
          // Triangle(j,λ,j) fails for j=1/2, λ=2, so Eta.OneBodyChannels
          // drops the s1/2 partners leftover χ actually needs.
          if (tensor_case and tbc_ket.J != J)
            continue;
          {
            for (auto b : Z.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
              auto ibra_bq =
                  tbc_bra.GetLocalIndex(std::min(b, q), std::max(b, q));
              if (ibra_bq < 0 or ibra_bq > (int)nbras)
                continue;
              double norm = (b == q ? PhysConst::SQRT2 : 1) *
                            (p == q ? 1 / PhysConst::SQRT2 : 1);
              if (b > q)
                norm *= bra.Phase(tbc_bra.J);
              zpqrs += norm * CHI_I(p, b) *
                       Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq,
                                                  iket);
            }
            for (auto b : Z.OneBodyChannels.at({oq.l, oq.j2, oq.tz2})) {
              auto ibra_pb =
                  tbc_bra.GetLocalIndex(std::min(p, b), std::max(p, b));
              if (ibra_pb < 0 or ibra_pb > (int)nbras)
                continue;
              double norm = (b == p ? PhysConst::SQRT2 : 1) *
                            (p == q ? 1 / PhysConst::SQRT2 : 1);
              if (p > b)
                norm *= bra.Phase(tbc_bra.J);
              zpqrs += norm * CHI_I(q, b) *
                       Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb,
                                                  iket);
            }
            for (auto b : Z.OneBodyChannels.at({oR.l, oR.j2, oR.tz2})) {
              auto iket_bs =
                  tbc_ket.GetLocalIndex(std::min(b, s), std::max(b, s));
              if (iket_bs < 0 or iket_bs > (int)nkets)
                continue;
              double norm = (b == s ? PhysConst::SQRT2 : 1) *
                            (r == s ? 1 / PhysConst::SQRT2 : 1);
              if (b > s)
                norm *= ket.Phase(tbc_ket.J);
              zpqrs +=
                  norm *
                  Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) *
                  CHI_I(b, r);
            }
            for (auto b : Z.OneBodyChannels.at({os.l, os.j2, os.tz2})) {
              auto iket_rb =
                  tbc_ket.GetLocalIndex(std::min(r, b), std::max(r, b));
              if (iket_rb < 0 or iket_rb > (int)nkets)
                continue;
              double norm = (b == r ? PhysConst::SQRT2 : 1) *
                            (r == s ? 1 / PhysConst::SQRT2 : 1);
              if (r > b)
                norm *= ket.Phase(tbc_ket.J);
              zpqrs +=
                  norm *
                  Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) *
                  CHI_I(b, s);
            }
          }

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
        }
      }
    }

    CHI_I.clear();
    if (Commutator::verbose) {
      Z.profiler.timer["_232_eths_GI"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
  }

  // ######################################################################
  // Gamma^II -- bra χ^{ΩΓ} (Ω×Γ), ket χ^ζ (Γ×Ω)
  //
  // Unfact: Γ^II = −1/2 Σ w {(1−P_ij) Ω_cjab Γ_abcd Ω_idkl
  //                           + (1−P_kl) Γ_cdab Ω_abcl Ω_ijkd}
  // Bra: χ^{ΩΓ}_ij = 1/2 Σ w Ω_ciab Γ_abcj  (AMC chi_omega_gamma_analyze)
  //      → −(1−P_ij) χ^{ΩΓ}_ja Ω_iakl
  // Ket: χ^ζ_ij = 1/2 Σ w Γ_ciab Ω_abcj
  //      → −(1−P_kl) χ^ζ_ak Ω_ijal
  // Fold (AMC G2_Wbra / G2_Wket) is the same 6j for any tensor 1b χ.
  // ######################################################################
  if (do_GII) {
    int max_j2 = 0;
    for (auto x : Z.modelspace->all_orbits)
      max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
    const int max_J = max_j2;
    std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                                Z.modelspace->all_orbits.end());
    const int norb = (int)allorb.size();

    arma::mat Chi_zeta(Gamma.OneBody.n_rows, Gamma.OneBody.n_cols,
                       arma::fill::zeros);
    arma::mat Chi_OG(Gamma.OneBody.n_rows, Gamma.OneBody.n_cols,
                     arma::fill::zeros);

    // χ^ζ / χ^{ΩΓ}: same dummy (a,b) DGEMM as leftover f^II (Γ W Ω), then 6j on c.
    auto pairs_at_J = [&](int J) {
      std::vector<std::array<index_t, 2>> pairs;
      pairs.reserve((size_t)norb * norb);
      for (auto x : allorb) {
        Orbit &ox = Z.modelspace->GetOrbit(x);
        const double jx = ox.j2 * 0.5;
        for (auto y : allorb) {
          Orbit &oy = Z.modelspace->GetOrbit(y);
          if (not AngMom::Triangle(jx, oy.j2 * 0.5, (double)J))
            continue;
          pairs.push_back({x, y});
        }
      }
      return pairs;
    };

    for (int J0 = 0; J0 <= max_J; ++J0) {
      auto pairs0 = pairs_at_J(J0);
      const int n0 = (int)pairs0.size();
      if (n0 < 1)
        continue;
      arma::mat G0(n0, n0, arma::fill::zeros);
      arma::vec w_hh0(n0), w_pp0(n0);
      for (int n = 0; n < n0; ++n) {
        Orbit &oa = Z.modelspace->GetOrbit(pairs0[n][0]);
        Orbit &ob = Z.modelspace->GetOrbit(pairs0[n][1]);
        const double na = oa.occ, nb = ob.occ;
        w_hh0(n) = na * nb;
        w_pp0(n) = (1.0 - na) * (1.0 - nb);
        for (int m = 0; m < n0; ++m)
          G0(n, m) = Gamma.TwoBody.GetTBME_J(J0, J0, pairs0[n][0],
                                             pairs0[n][1], pairs0[m][0],
                                             pairs0[m][1]);
      }

      for (int J1 = 0; J1 <= max_J; ++J1) {
        if (not AngMom::Triangle(J0, J1, lambda))
          continue;
        auto pairs1 = pairs_at_J(J1);
        const int n1 = (int)pairs1.size();
        if (n1 < 1)
          continue;
        const double hats01 = hat(J0) * hat(J1);

        // χ^ζ: Γ^{J0}_{ci,ab} W_ab Ω^{J0 J1}_{ab,cj}
        arma::mat OmZ(n0, n1, arma::fill::zeros);
        for (int m = 0; m < n0; ++m)
          for (int k = 0; k < n1; ++k)
            OmZ(m, k) = Eta.TwoBody.GetTBME_J(
                J0, J1, pairs0[m][0], pairs0[m][1], pairs1[k][0],
                pairs1[k][1]);
        arma::mat TZ_hh = hats01 * (G0 * arma::diagmat(w_hh0) * OmZ);
        arma::mat TZ_pp = hats01 * (G0 * arma::diagmat(w_pp0) * OmZ);

        // χ^{ΩΓ}: Ω^{J0 J1}_{ci,ab} W_ab Γ^{J1}_{ab,cj}
        arma::mat OmOG(n0, n1, arma::fill::zeros);
        arma::mat G1(n1, n1, arma::fill::zeros);
        arma::vec w_hh1(n1), w_pp1(n1);
        for (int m = 0; m < n1; ++m) {
          Orbit &oa = Z.modelspace->GetOrbit(pairs1[m][0]);
          Orbit &ob = Z.modelspace->GetOrbit(pairs1[m][1]);
          const double na = oa.occ, nb = ob.occ;
          w_hh1(m) = na * nb;
          w_pp1(m) = (1.0 - na) * (1.0 - nb);
          for (int k = 0; k < n1; ++k)
            G1(m, k) = Gamma.TwoBody.GetTBME_J(J1, J1, pairs1[m][0],
                                               pairs1[m][1], pairs1[k][0],
                                               pairs1[k][1]);
        }
        for (int n = 0; n < n0; ++n)
          for (int m = 0; m < n1; ++m)
            OmOG(n, m) = Eta.TwoBody.GetTBME_J(
                J0, J1, pairs0[n][0], pairs0[n][1], pairs1[m][0],
                pairs1[m][1]);
        arma::mat TOG_hh = hats01 * (OmOG * arma::diagmat(w_hh1) * G1);
        arma::mat TOG_pp = hats01 * (OmOG * arma::diagmat(w_pp1) * G1);

        for (int n = 0; n < n0; ++n) {
          const index_t c = pairs0[n][0], i = pairs0[n][1];
          Orbit &oc = Z.modelspace->GetOrbit(c);
          Orbit &oi = Z.modelspace->GetOrbit(i);
          const double jc = oc.j2 * 0.5, ji = oi.j2 * 0.5;
          const double n_c = oc.occ, nbar_c = 1.0 - n_c;
          for (int k = 0; k < n1; ++k) {
            if (pairs1[k][0] != c)
              continue;
            const index_t j = pairs1[k][1];
            Orbit &oj = Z.modelspace->GetOrbit(j);
            const double jj = oj.j2 * 0.5;
            const double ph = Z.modelspace->phase(
                (oj.j2 + oc.j2) / 2 + lambda + J0);

            const double blkZ = nbar_c * TZ_hh(n, k) + n_c * TZ_pp(n, k);
            if (std::abs(blkZ) > 1e-16) {
              const double sixj =
                  AngMom::SixJ(lambda, J1, J0, jc, ji, jj);
              if (std::abs(sixj) > 1e-16)
                Chi_zeta(i, j) += 0.5 * ph * sixj * blkZ;
            }
            const double blkOG = nbar_c * TOG_hh(n, k) + n_c * TOG_pp(n, k);
            if (std::abs(blkOG) > 1e-16) {
              const double sixj =
                  AngMom::SixJ(J0, J1, lambda, jj, ji, jc);
              if (std::abs(sixj) > 1e-16)
                Chi_OG(i, j) += 0.5 * ph * sixj * blkOG;
            }
          }
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_232_eths_chizeta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // Fold: Z = −W_OG − V_ζ. W: AMC Wbra with χ^{ΩΓ}; V: Wket with χ^ζ.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nch; ++ich) {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      const int J = tbc_bra.J;
      if (tbc_ket.J != J)
        continue;
      const double hatJ_inv = 1.0 / hat(J);
      const size_t nbras = tbc_bra.GetNumberKets();
      const size_t nkets = tbc_ket.GetNumberKets();

      for (size_t ibra = 0; ibra < nbras; ++ibra) {
        Ket &bra = tbc_bra.GetKet(ibra);
        const index_t i = bra.p, j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = (int)ibra;
        for (size_t iket = ketmin; iket < nkets; ++iket) {
          Ket &ket = tbc_ket.GetKet(iket);
          const index_t k = ket.p, l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;

          double W = 0.0; // bra (1-P_ij) χ^{ΩΓ}_ja Ω_iakl
          double V = 0.0; // ket (1-P_kl) χ^ζ_ak Ω_ijal
          for (auto a : allorb) {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            const double ja = oa.j2 * 0.5;
            for (int J2 = 0; J2 <= max_J; ++J2) {
              if (not AngMom::Triangle(J, J2, lambda))
                continue;
              const double pref = hatJ_inv * hat(J2) * hat_lambda_inv;

              {
                const double c = Chi_OG(j, a);
                if (std::abs(c) > 1e-16) {
                  const double sixj =
                      AngMom::SixJ(J, J2, lambda, ja, jj, ji);
                  const double om =
                      Eta.TwoBody.GetTBME_J(J2, J, i, a, k, l);
                  if (std::abs(sixj * om) > 1e-16) {
                    const double ph =
                        Z.modelspace->phase((oi.j2 + oa.j2) / 2 + J2);
                    W += ph * pref * sixj * c * om;
                  }
                }
              }
              {
                const double c = Chi_OG(i, a);
                if (std::abs(c) > 1e-16) {
                  const double sixj =
                      AngMom::SixJ(J, J2, lambda, ja, ji, jj);
                  const double om =
                      Eta.TwoBody.GetTBME_J(J2, J, j, a, k, l);
                  if (std::abs(sixj * om) > 1e-16) {
                    const double ph = Z.modelspace->phase(
                        J + (oi.j2 + oa.j2) / 2 + J2);
                    W += ph * pref * sixj * c * om;
                  }
                }
              }
              {
                const double c = Chi_zeta(a, k);
                if (std::abs(c) > 1e-16) {
                  const double sixj =
                      AngMom::SixJ(J, J2, lambda, ja, jk, jl);
                  const double om =
                      Eta.TwoBody.GetTBME_J(J, J2, i, j, a, l);
                  if (std::abs(sixj * om) > 1e-16) {
                    const double ph =
                        Z.modelspace->phase((ol.j2 + oa.j2) / 2 + J2);
                    V += ph * pref * sixj * c * om;
                  }
                }
              }
              {
                const double c = Chi_zeta(a, l);
                if (std::abs(c) > 1e-16) {
                  const double sixj =
                      AngMom::SixJ(J, J2, lambda, ja, jl, jk);
                  const double om =
                      Eta.TwoBody.GetTBME_J(J, J2, i, j, a, k);
                  if (std::abs(sixj * om) > 1e-16) {
                    const double ph = Z.modelspace->phase(
                        J + (ol.j2 + oa.j2) / 2 + J2);
                    V += ph * pref * sixj * c * om;
                  }
                }
              }
            }
          }

          double z = -W - V; // Γ^II = −W_OG − V_ζ
          if (i == j)
            z /= PhysConst::SQRT2;
          if (k == l)
            z /= PhysConst::SQRT2;
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, z);
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_232_eths_GII"] += omp_get_wtime() - t_internal;
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
} // comm223_232_chi1b_tensor


////////////////////////////////////////////////////////////////////////////
/// factorized 223_232 double commutator with 2b intermediate
////////////////////////////////////////////////////////////////////////////
void comm223_232_chi2b(const Operator &Eta, const Operator &Gamma,
                       Operator &Z) {
  const double t_start = omp_get_wtime();

  // FDC.cc comm223_232_chi2b block order (tensor Ω / scalar Γ,Z):
  //   ~1227 allocate bar_Omega, bar_Gamma
  //   ~1243 combined Pandya + nnn occupancy
  //   ~1385 DGEMM χ  (ι,κ here; η is AMC Path B; θ,λ stay diagram-local)
  //   ~1413 GIIIa  χ^η ladder
  //   ~1528 GIVa   χ^κ
  //   ~1906 GIVb   χ^ι
  //   ~2393 GIIIc  χ^θ
  //   ~2739 GIIIb  (ethS gold is leftover IIb/IId, not FDC recouple of χ^η)
  //   ~2916 GIVc   χ^λ  (Pandya of χ^λ and Ω together, same 9j)
  Tensor232Bars bars;

  if (use_TypeGIIIa_2b)
    comm223_232_GIIIa(Eta, Gamma, Z);
  if (use_TypeGIVa_2b)
    comm223_232_GIVa_from_bars(Eta, Gamma, Z, bars);
  if (use_TypeGIVb_2b)
    comm223_232_GIVb_from_bars(Eta, Gamma, Z, bars);
  if (use_TypeGIIIc_2b)
    comm223_232_GIIIc(Eta, Gamma, Z);
  if (use_TypeGIIIb_2b)
    comm223_232_GIIIb(Eta, Gamma, Z);
  if (use_TypeGIVc_2b)
    comm223_232_GIVc(Eta, Gamma, Z);

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// χ^η → Γ^{III_a}  (ladder → ordinary-channel DGEMM once χ is scalar)
///
/// Rule: [Ω×Ω]^0 → scalar χ (λ_χ=0). Ladder = Chi*Gamma DGEMM.
/// Fill is still G3a Term strips until Pandya χ̄ matches; contraction stays
/// explicit (1-P) until χ is in normalized ket basis.
/// GIIIb: non-ladder → Pandya/CC → DGEMM → inverse (TTS stub for now).
////////////////////////////////////////////////////////////////////////////
namespace {

struct ChiTab {
  int n_orb = 0, max_J = 0;
  size_t sL = 0, sK = 0, sJ = 0, sI = 0;
  std::vector<double> data;

  void allocate(int n, int Jmax) {
    n_orb = n;
    max_J = Jmax;
    sL = (size_t)(max_J + 1);
    sK = (size_t)n_orb * sL;
    sJ = (size_t)n_orb * sK;
    sI = (size_t)n_orb * sJ;
    data.assign((size_t)n_orb * sI, 0.0);
  }
  double &at(index_t p, index_t q, index_t r, index_t s, int J) {
    return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                (size_t)s * sL + (size_t)J];
  }
  double operator()(index_t p, index_t q, index_t r, index_t s, int J) const {
    return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                (size_t)s * sL + (size_t)J];
  }
};




/// Factorized 2n layout from ChiTab (unreduced χ = S/Ĵ²).
/// Do not go through a ket-ordered Operator + fermionic AS: χ^θ has occupancy
/// on j and k, so χ_{jikl} ≠ −phase χ_{ijkl} when n_i ≠ n_j. That happens for
/// even-parity ph kets (0s+1s/0d) at emax≥2; emax=1 only has even kets with
/// n_i = n_j (hh 0s0s, pp 0p0p).
void PackChiThetaFactLayout(const ChiTab &chi_k, const ChiTab &chi_j,
                            Operator &Z, int ch, int J0, arma::mat &Out) {
  TwoBodyChannel &tb = Z.modelspace->GetTwoBodyChannel(ch);
  const int nb = tb.GetNumberKets();
  Out = arma::mat(2 * nb, 2 * nb, arma::fill::zeros);
  if (nb < 1)
    return;
  const double inv_hat2 = 1.0 / (2.0 * J0 + 1.0);
  auto pq = [&](int idx) -> std::array<int, 2> {
    Ket &ket = tb.GetKet(idx % nb);
    if (idx < nb)
      return {(int)ket.p, (int)ket.q};
    if (ket.p == ket.q)
      return {-1, -1};
    return {(int)ket.q, (int)ket.p};
  };
  for (int ibra = 0; ibra < 2 * nb; ++ibra) {
    auto ij = pq(ibra);
    if (ij[0] < 0)
      continue;
    for (int iket = 0; iket < 2 * nb; ++iket) {
      auto kl = pq(iket);
      if (kl[0] < 0)
        continue;
      if (ij[0] >= chi_k.n_orb or ij[1] >= chi_k.n_orb or
          kl[0] >= chi_k.n_orb or kl[1] >= chi_k.n_orb or J0 > chi_k.max_J)
        continue;
      const double S = chi_k(ij[0], ij[1], kl[0], kl[1], J0) +
                       chi_j(ij[0], ij[1], kl[0], kl[1], J0);
      Out(ibra, iket) = S * inv_hat2;
    }
  }
}

// chi^theta (G3c): reduced [Ω⊗Ω]^(0) via ordinary-channel DGEMM (production).
// RME: ⟨J‖[Ω⊗Ω]^(0)‖J⟩ = Σ_{J'} (-1)^{λ+J-J'}/(Ĵ λ̂) Ω^{JJ'} (occ⊙Ω^{J'J})
// AMC/TTS packaging (matches prior orbit loops / G3c): ang = (-1)^{J+J'+λ} λ̂^{-1}
//   (Ĵ absorbed in downstream sixj hats). Two occ-sided tables:
//   χ_k: weight n_k W_pp + n̄_k W_hh on ket; χ_j: same T, weight on bra j.
void FillChiThetaG3c_DGEMM(const Operator &Eta, Operator &Z, ChiTab &chi_k,
                           ChiTab &chi_j, double hat_lambda_inv, int lambda,
                           const std::vector<index_t> &allorb) {
  const int n_orb = (int)allorb.size();
  const int max_J = chi_k.max_J;
  // zero
  std::fill(chi_k.data.begin(), chi_k.data.end(), 0.0);
  std::fill(chi_j.data.begin(), chi_j.data.end(), 0.0);

  for (int J0 = 0; J0 <= max_J; ++J0) {
    std::vector<std::array<index_t, 2>> pairs0;
    pairs0.reserve((size_t)n_orb * n_orb);
    for (auto i : allorb) {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      const double ji = oi.j2 * 0.5;
      for (auto j : allorb) {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        if (not AngMom::Triangle(ji, oj.j2 * 0.5, (double)J0))
          continue;
        pairs0.push_back({i, j});
      }
    }
    const int n0 = (int)pairs0.size();
    if (n0 < 1)
      continue;

    arma::mat T_pp(n0, n0, arma::fill::zeros);
    arma::mat T_hh(n0, n0, arma::fill::zeros);

    for (int J1 = 0; J1 <= max_J; ++J1) {
      if (not AngMom::Triangle(J0, J1, lambda))
        continue;

      std::vector<std::array<index_t, 2>> pairs1;
      pairs1.reserve((size_t)n_orb * n_orb);
      for (auto a : allorb) {
        Orbit &oa = Z.modelspace->GetOrbit(a);
        const double ja = oa.j2 * 0.5;
        for (auto b : allorb) {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          if (not AngMom::Triangle(ja, ob.j2 * 0.5, (double)J1))
            continue;
          pairs1.push_back({a, b});
        }
      }
      const int n1 = (int)pairs1.size();
      if (n1 < 1)
        continue;

      arma::vec w_pp(n1), w_hh(n1);
      for (int m = 0; m < n1; ++m) {
        Orbit &oa = Z.modelspace->GetOrbit(pairs1[m][0]);
        Orbit &ob = Z.modelspace->GetOrbit(pairs1[m][1]);
        const double na = oa.occ, nb = ob.occ;
        w_pp(m) = (1.0 - na) * (1.0 - nb);
        w_hh(m) = na * nb;
      }

      arma::mat Om(n0, n1, arma::fill::zeros);
      arma::mat Om_r(n1, n0, arma::fill::zeros);
      for (int n = 0; n < n0; ++n) {
        const index_t i = pairs0[n][0], j = pairs0[n][1];
        for (int m = 0; m < n1; ++m) {
          const index_t a = pairs1[m][0], b = pairs1[m][1];
          Om(n, m) = Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b);
        }
      }
      for (int m = 0; m < n1; ++m) {
        const index_t a = pairs1[m][0], b = pairs1[m][1];
        for (int n = 0; n < n0; ++n) {
          const index_t k = pairs0[n][0], l = pairs0[n][1];
          Om_r(m, n) = Eta.TwoBody.GetTBME_J(J1, J0, a, b, k, l);
        }
      }

      const double ang =
          Z.modelspace->phase(J0 + J1 + lambda) * hat_lambda_inv;
      T_pp += ang * (Om * arma::diagmat(w_pp) * Om_r);
      T_hh += ang * (Om * arma::diagmat(w_hh) * Om_r);
    }

    for (int n = 0; n < n0; ++n) {
      const index_t i = pairs0[n][0], j = pairs0[n][1];
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const double n_j = oj.occ, nbar_j = 1.0 - n_j;
      for (int m = 0; m < n0; ++m) {
        const index_t k = pairs0[m][0], l = pairs0[m][1];
        Orbit &ok = Z.modelspace->GetOrbit(k);
        const double n_k = ok.occ, nbar_k = 1.0 - n_k;
        // occ = W_hh n̄ + W_pp n  →  χ = n̄ T_hh + n T_pp
        const double ck = n_k * T_pp(n, m) + nbar_k * T_hh(n, m);
        const double cj = n_j * T_pp(n, m) + nbar_j * T_hh(n, m);
        chi_k.at(i, j, k, l, J0) = ck;
        chi_j.at(i, j, k, l, J0) = cj;
      }
    }
  }

}



std::array<int, 2> PqCC(TwoBodyChannel_CC &tbc, int idx, int nK) {
  Ket &ket = tbc.GetKet(idx % nK);
  if (idx < nK)
    return {(int)ket.p, (int)ket.q};
  if (ket.p == ket.q)
    return {-1, -1};
  return {(int)ket.q, (int)ket.p};
}

std::array<int, 2> PqTB(TwoBodyChannel &tbc, int idx, int nK) {
  Ket &ket = tbc.GetKet(idx % nK);
  if (idx < nK)
    return {(int)ket.p, (int)ket.q};
  if (ket.p == ket.q)
    return {-1, -1};
  return {(int)ket.q, (int)ket.p};
}

// All-orbit χ table. Occupancy-asymmetric χ is not a fermionic TBME: do not
// pack it in TwoBodyChannel 2n (Pauli / p≤q). Same layout as GIVc ChiTabJJ.
struct ChiTabJJ {
  int max_J = 0;
  size_t sJ1 = 0, sL = 0, sK = 0, sJ = 0, sI = 0;
  std::vector<double> data;
  void allocate(int n, int Jmax) {
    max_J = Jmax;
    sJ1 = (size_t)(max_J + 1);
    sL = sJ1 * (size_t)(max_J + 1);
    sK = (size_t)n * sL;
    sJ = (size_t)n * sK;
    sI = (size_t)n * sJ;
    data.assign((size_t)n * sI, 0.0);
  }
  double &at(index_t p, index_t q, index_t r, index_t s, int J0, int J1) {
    return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                (size_t)s * sL + (size_t)J0 * sJ1 + (size_t)J1];
  }
  double operator()(index_t p, index_t q, index_t r, index_t s, int J0,
                    int J1) const {
    return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                (size_t)s * sL + (size_t)J0 * sJ1 + (size_t)J1];
  }
};

struct OrbPack {
  std::vector<index_t> allorb;
  int n_orb = 0, n_alloc = 0, max_J = 0;
  std::vector<int> compact, oj2;
  std::vector<double> ojh, nocc;
};

OrbPack MakeOrbPack(ModelSpace *ms) {
  OrbPack P;
  P.allorb.assign(ms->all_orbits.begin(), ms->all_orbits.end());
  P.n_orb = (int)P.allorb.size();
  int max_j2 = 0;
  for (auto o : P.allorb) {
    max_j2 = std::max(max_j2, ms->GetOrbit(o).j2);
    P.n_alloc = std::max(P.n_alloc, (int)o + 1);
  }
  P.max_J = max_j2;
  P.compact.assign(P.n_alloc, -1);
  P.oj2.resize(P.n_orb);
  P.ojh.resize(P.n_orb);
  P.nocc.resize(P.n_orb);
  for (int t = 0; t < P.n_orb; ++t) {
    P.compact[(int)P.allorb[t]] = t;
    Orbit &o = ms->GetOrbit(P.allorb[t]);
    P.oj2[t] = o.j2;
    P.ojh[t] = o.j2 * 0.5;
    P.nocc[t] = o.occ;
  }
  return P;
}

double HatJ(int J) { return std::sqrt(2.0 * J + 1.0); }

// IMSRG tensor Pandya adcb (FillTensor232Bars / Python Path B).
double PandyaOmegaAdcb(const Operator &Eta, int a, int b, int c, int d, int Jbra,
                       int Jket) {
  const int lambda = Eta.GetJRank();
  if (not AngMom::Triangle(Jbra, Jket, lambda))
    return 0.0;
  ModelSpace *ms = Eta.modelspace;
  Orbit &oa = ms->GetOrbit(a);
  Orbit &ob = ms->GetOrbit(b);
  Orbit &oc = ms->GetOrbit(c);
  Orbit &od = ms->GetOrbit(d);
  const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5, jc = oc.j2 * 0.5,
               jd = od.j2 * 0.5;
  if (not AngMom::Triangle(ja, jb, (double)Jbra) or
      not AngMom::Triangle(jc, jd, (double)Jket))
    return 0.0;
  double sm = 0.0;
  const int j1min = std::abs(oa.j2 - od.j2) / 2;
  const int j1max = (oa.j2 + od.j2) / 2;
  for (int J1 = j1min; J1 <= j1max; ++J1) {
    const int j2min =
        std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
    const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
    for (int J2 = j2min; J2 <= j2max; ++J2) {
      const double ninej =
          ms->GetNineJ(ja, jd, J1, jb, jc, J2, Jbra, Jket, lambda);
      if (std::abs(ninej) < 1e-14)
        continue;
      const double hats = HatJ(J1) * HatJ(J2) * HatJ(Jbra) * HatJ(Jket);
      const double tbme = Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
      sm -= hats * ms->phase((ob.j2 + od.j2) / 2 + Jket + J2) * ninej * tbme;
    }
  }
  return sm;
}

// Tensor Pandya adcb of occupancy-AS χ^λ (ChiTab, not a fermionic TBME).
double PandyaChiAdcb(const ChiTabJJ &Chi, ModelSpace *ms, int lambda, int a,
                     int b, int c, int d, int Jbra, int Jket) {
  if (not AngMom::Triangle(Jbra, Jket, lambda))
    return 0.0;
  Orbit &oa = ms->GetOrbit(a);
  Orbit &ob = ms->GetOrbit(b);
  Orbit &oc = ms->GetOrbit(c);
  Orbit &od = ms->GetOrbit(d);
  const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5, jc = oc.j2 * 0.5,
               jd = od.j2 * 0.5;
  if (not AngMom::Triangle(ja, jb, (double)Jbra) or
      not AngMom::Triangle(jc, jd, (double)Jket))
    return 0.0;
  double sm = 0.0;
  const int j1min = std::abs(oa.j2 - od.j2) / 2;
  const int j1max = (oa.j2 + od.j2) / 2;
  for (int J1 = j1min; J1 <= j1max; ++J1) {
    const int j2min =
        std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
    const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
    for (int J2 = j2min; J2 <= j2max; ++J2) {
      const double ninej =
          ms->GetNineJ(ja, jd, J1, jb, jc, J2, Jbra, Jket, lambda);
      if (std::abs(ninej) < 1e-14)
        continue;
      const double hats = HatJ(J1) * HatJ(J2) * HatJ(Jbra) * HatJ(Jket);
      const double tbme = Chi(a, d, c, b, J1, J2);
      sm -= hats * ms->phase((ob.j2 + od.j2) / 2 + Jket + J2) * ninej * tbme;
    }
  }
  return sm;
}

// Scalar Pandya adcb (Python bar_Gamma). Identical legs: even Jstd only.
double PandyaGammaAdcb(const Operator &Gamma, int a, int b, int c, int d,
                       int Jcc) {
  ModelSpace *ms = Gamma.modelspace;
  Orbit &oa = ms->GetOrbit(a);
  Orbit &ob = ms->GetOrbit(b);
  Orbit &oc = ms->GetOrbit(c);
  Orbit &od = ms->GetOrbit(d);
  const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5, jc = oc.j2 * 0.5,
               jd = od.j2 * 0.5;
  if (not AngMom::Triangle(ja, jb, (double)Jcc) or
      not AngMom::Triangle(jc, jd, (double)Jcc))
    return 0.0;
  int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
  int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
  int dJ = 1;
  if (a == d or b == c) {
    dJ = 2;
    jmin += jmin % 2;
  }
  double sm = 0.0;
  for (int Jstd = jmin; Jstd <= jmax; Jstd += dJ) {
    const double six = ms->GetSixJ(ja, jb, Jcc, jc, jd, Jstd);
    if (std::abs(six) < 1e-8)
      continue;
    sm -= (2.0 * Jstd + 1.0) * six *
          Gamma.TwoBody.GetTBME_J(Jstd, Jstd, a, d, c, b);
  }
  return sm;
}

// AMC tensor inv WITHOUT printed leading minus (Path B / IMSRG via Eq4).
// barChi(i,l,k,j,J2,J3) = χ̄(il;kj)^{J2 J3}.
double InvTensorPlus(const ChiTabJJ &barChi, ModelSpace *ms, int lambda,
                     int max_J, int i, int j, int k, int l, int J0, int J1) {
  Orbit &oi = ms->GetOrbit(i);
  Orbit &oj = ms->GetOrbit(j);
  Orbit &ok = ms->GetOrbit(k);
  Orbit &ol = ms->GetOrbit(l);
  const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5, jk = ok.j2 * 0.5,
               jl = ol.j2 * 0.5;
  if (not AngMom::Triangle(J0, J1, lambda))
    return 0.0;
  if (not AngMom::Triangle(ji, jj, (double)J0) or
      not AngMom::Triangle(jk, jl, (double)J1))
    return 0.0;
  double tot = 0.0;
  for (int J2 = 0; J2 <= max_J; ++J2) {
    for (int J3 = 0; J3 <= max_J; ++J3) {
      if (not AngMom::Triangle(J2, J3, lambda))
        continue;
      if (not AngMom::Triangle(ji, jl, (double)J2) or
          not AngMom::Triangle(jk, jj, (double)J3))
        continue;
      const double bc = barChi(i, l, k, j, J2, J3);
      if (std::abs(bc) < 1e-16)
        continue;
      const double nj =
          ms->GetNineJ(lambda, J0, J1, J3, jj, jk, J2, ji, jl);
      if (std::abs(nj) < 1e-16)
        continue;
      tot += ms->phase(J2) * HatJ(J2) * HatJ(J3) * nj * bc;
    }
  }
  return ms->phase(J0 + (oi.j2 + ok.j2) / 2 + lambda) * HatJ(J0) * HatJ(J1) *
         tot;
}


// Reduced tensor CC flip (LESSONS.md Γ^{IV_a} / comm222_phst):
//   Ω̄^{J1 J0}(cd;ab) = hΩ (−1)^{J1−J0} Ω̄^{J0 J1}(ab;cd)
// Occupancy-weighted χ does not obey this — only Ω̄.
arma::mat TensorBarOmegaPartner(const arma::mat &Om_J0J1, int J0, int J1,
                                int hOmega, ModelSpace *ms) {
  return (hOmega * ms->phase(J0 - J1)) * Om_J0J1.t();
}

/// AMC inv: χ_red = Ĵ_0 Σ_J' Ĵ' 6j χ̄_CC. No (1−P). Gold: test_chi_eta_mscheme.py.
double InvChiEtaRed(Operator &Z, const std::deque<arma::mat> &barCHI, int i,
                    int j, int k, int l, int J0) {
  const int n_cc = (int)barCHI.size();
  Orbit &oi = Z.modelspace->GetOrbit(i);
  Orbit &oj = Z.modelspace->GetOrbit(j);
  Orbit &ok = Z.modelspace->GetOrbit(k);
  Orbit &ol = Z.modelspace->GetOrbit(l);
  const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
  const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
  double sm = 0.0;
  const int parity_cc = (oi.l + ol.l) % 2;
  const int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
  const int Jpmin =
      std::max(std::abs(oi.j2 - ol.j2), std::abs(ok.j2 - oj.j2)) / 2;
  const int Jpmax = std::min(oi.j2 + ol.j2, ok.j2 + oj.j2) / 2;
  for (int Jp = Jpmin; Jp <= Jpmax; ++Jp) {
    const double sixj = Z.modelspace->GetSixJ(jl, jk, (double)J0, jj, ji, Jp);
    if (std::abs(sixj) < 1e-16)
      continue;
    const size_t ch_cc =
        Z.modelspace->GetTwoBodyChannelIndex(Jp, parity_cc, Tz_cc);
    if ((int)ch_cc >= n_cc)
      continue;
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    const int nkets_cc = tbc_cc.GetNumberKets();
    if (nkets_cc < 1)
      continue;
    int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
    int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
    if (indx_il < 0 or indx_kj < 0)
      continue;
    indx_il += (i > l ? nkets_cc : 0);
    indx_kj += (k > j ? nkets_cc : 0);
    if (indx_il >= (int)barCHI[ch_cc].n_rows or
        indx_kj >= (int)barCHI[ch_cc].n_cols)
      continue;
    sm += std::sqrt(2.0 * Jp + 1.0) * sixj * barCHI[ch_cc](indx_il, indx_kj);
  }
  return std::sqrt(2.0 * J0 + 1.0) * sm;
}

/// χ^η Path B (GIIIa / GIIIb): AMC same-label Pandya of Ω → occ DGEMM.
/// χ is scalar and not AS. Store 2n ordinary (not n×n GetMatrix).
/// GIIIb Fac Pandya of χ inverts on the fly — CC legs ≠ ordinary kets.
void BuildChiEtaPathB(const Operator &Eta, Operator &Z,
                      std::deque<arma::mat> &Chi2,
                      std::deque<arma::mat> *barCHI_out, bool fill_chi2) {
  const int lambda = Eta.GetJRank();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();
  const int nch_eta = Z.modelspace->GetNumberTwoBodyChannels();
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  const double hat_lambda = std::sqrt(2.0 * lambda + 1.0);
  double t_start = omp_get_wtime();

  auto amc_bar_omega = [&](int i, int j, int k, int l, int Jbra,
                           int Jket) -> double {
    if (not AngMom::Triangle(Jbra, Jket, lambda))
      return 0.0;
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
    const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
    double sm = 0.0;
    const int J2min = std::abs(oi.j2 - oj.j2) / 2;
    const int J2max = (oi.j2 + oj.j2) / 2;
    const int J3min = std::abs(ok.j2 - ol.j2) / 2;
    const int J3max = (ok.j2 + ol.j2) / 2;
    for (int J2 = J2min; J2 <= J2max; ++J2) {
      for (int J3 = J3min; J3 <= J3max; ++J3) {
        if (not AngMom::Triangle(J2, J3, lambda))
          continue;
        const double n9 = Z.modelspace->GetNineJ(
            (double)lambda, (double)Jbra, (double)Jket, (double)J3, jl, jk,
            (double)J2, ji, jj);
        if (std::abs(n9) < 1e-16)
          continue;
        const double hats = std::sqrt((2.0 * J2 + 1.0) * (2.0 * J3 + 1.0));
        sm += Z.modelspace->phase(J2) * hats * n9 *
              Eta.TwoBody.GetTBME_J(J2, J3, i, j, k, l);
      }
    }
    const double pref =
        -Z.modelspace->phase(Jbra + (oi.j2 + ok.j2) / 2 + lambda) *
        std::sqrt((2.0 * Jbra + 1.0) * (2.0 * Jket + 1.0));
    return pref * sm;
  };

  // 9j Pandya of Ω̄ only for ch_b ≤ ch_k. Partner block is reduced transpose
  // hΩ (−1)^{J0−J1} Ω̄^T (LESSONS.md). Occ χ̄ is applied after, on that layout.
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  std::vector<std::array<int, 2>> cc_all, cc_canon;
  for (int ch_b = 0; ch_b < n_cc; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    const int nb = tb.GetNumberKets();
    if (nb < 1)
      continue;
    for (int ch_k = 0; ch_k < n_cc; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      if (tk.GetNumberKets() < 1)
        continue;
      if (not AngMom::Triangle(tb.J, tk.J, lambda))
        continue;
      if ((tb.parity + tk.parity + Eta.GetParity()) % 2 != 0)
        continue;
      if (not CcTzCouples(tb.Tz, tk.Tz, Eta.GetTRank()))
        continue;
      cc_all.push_back({ch_b, ch_k});
      if (ch_b <= ch_k)
        cc_canon.push_back({ch_b, ch_k});
    }
  }
  const int nc = (int)cc_canon.size();
  std::vector<arma::mat> tmpOm(nc);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ic = 0; ic < nc; ++ic) {
    const int ch_b = cc_canon[ic][0], ch_k = cc_canon[ic][1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    const int Jb = tb.J, Jk = tk.J;
    arma::mat Om(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto il = PqCC(tb, ibra, nb);
      if (il[0] < 0)
        continue;
      const int i = il[0], lorb = il[1];
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto kj = PqCC(tk, iket, nk);
        if (kj[0] < 0)
          continue;
        Om(ibra, iket) = amc_bar_omega(i, kj[1], kj[0], lorb, Jb, Jk);
      }
    }
    tmpOm[ic] = std::move(Om);
  }
  std::map<std::array<int, 2>, arma::mat> bar_Omega;
  for (int ic = 0; ic < nc; ++ic)
    bar_Omega[cc_canon[ic]] = std::move(tmpOm[ic]);
  for (const auto &pr : cc_all) {
    const int ch_b = pr[0], ch_k = pr[1];
    if (ch_b <= ch_k)
      continue;
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    auto it = bar_Omega.find({ch_k, ch_b});
    if (it == bar_Omega.end())
      continue;
    bar_Omega[pr] =
        TensorBarOmegaPartner(it->second, tk.J, tb.J, hEta, Z.modelspace);
  }

  if (Commutator::verbose) {
    Z.profiler.timer["BuildChiEtaPathB_pandya"] += omp_get_wtime() - t_start;
  }
  double t_mid = omp_get_wtime();

  std::deque<arma::mat> barCHI(n_cc);
  for (int ch_b = 0; ch_b < n_cc; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    const int nb = tb.GetNumberKets();
    barCHI[ch_b] = arma::mat(2 * nb, 2 * nb, arma::fill::zeros);
    if (nb < 1)
      continue;
    const int Jb = tb.J;
    const double wL = Z.modelspace->phase(Jb) / std::sqrt(2.0 * Jb + 1.0);
    for (int ch_k = 0; ch_k < n_cc; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      const int nk = tk.GetNumberKets();
      if (nk < 1)
        continue;
      if (not AngMom::Triangle(Jb, tk.J, lambda))
        continue;
      auto itL = bar_Omega.find({ch_b, ch_k});
      auto itR = bar_Omega.find({ch_k, ch_b});
      if (itL == bar_Omega.end() or itR == bar_Omega.end())
        continue;
      const int Jk = tk.J;
      arma::mat Rocc = itR->second;
      for (int ibra = 0; ibra < 2 * nk; ++ibra) {
        auto ba = PqCC(tk, ibra, nk);
        if (ba[0] < 0)
          continue;
        const int b = ba[0], a = ba[1];
        Orbit &oa = Z.modelspace->GetOrbit(a);
        Orbit &ob = Z.modelspace->GetOrbit(b);
        const double na = oa.occ, nbar_a = 1.0 - na;
        const double nb_occ = ob.occ, nbar_b = 1.0 - nb_occ;
        const double w_ab_hh = nbar_a * nb_occ;
        const double w_ab_pp = na * nbar_b;
        for (int iket = 0; iket < 2 * nb; ++iket) {
          auto kj = PqCC(tb, iket, nb);
          if (kj[0] < 0)
            continue;
          Orbit &ok = Z.modelspace->GetOrbit(kj[0]);
          const double nk_occ = ok.occ;
          const double nbar_k = 1.0 - nk_occ;
          Rocc(ibra, iket) *= (w_ab_hh * nbar_k + w_ab_pp * nk_occ);
        }
      }
      const double wR = Z.modelspace->phase(Jk + lambda) / hat_lambda;
      barCHI[ch_b] += (wL * wR) * (itL->second * Rocc);
    }
  }
  bar_Omega.clear();

  if (Commutator::verbose) {
    Z.profiler.timer["BuildChiEtaPathB_dgemm_chi"] += omp_get_wtime() - t_mid;
  }
  t_mid = omp_get_wtime();

  if (fill_chi2) {
    Chi2.assign(nch_eta, arma::mat());
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch_eta; ++ch) {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      const int nK = tbc.GetNumberKets();
      const int J0 = tbc.J;
      Chi2[ch] = arma::mat(2 * nK, 2 * nK, arma::fill::zeros);
      if (nK < 1)
        continue;
      for (int ibra = 0; ibra < 2 * nK; ++ibra) {
        auto ij = PqTB(tbc, ibra, nK);
        if (ij[0] < 0)
          continue;
        for (int iket = 0; iket < 2 * nK; ++iket) {
          auto kl = PqTB(tbc, iket, nK);
          if (kl[0] < 0)
            continue;
          Chi2[ch](ibra, iket) =
              InvChiEtaRed(Z, barCHI, ij[0], ij[1], kl[0], kl[1], J0);
        }
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler.timer["BuildChiEtaPathB_inv"] += omp_get_wtime() - t_mid;
  }

  if (barCHI_out)
    *barCHI_out = std::move(barCHI);
  else {
    for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc)
      barCHI[ch_cc].clear();
    barCHI.clear();
  }
}

} // namespace

void comm223_232_GIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  // Γ^{III_a}: AMC Path B χ^η (same-label Pandya → occ DGEMM → inv)
  // then ordinary-channel Chi_AS×Γ ladder. χ is scalar and not AS.
  // Gold: test_chi_eta_mscheme.py + test_GIIIa_ladder_mscheme.py
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  auto &Z2 = Z.TwoBody;

  std::deque<arma::mat> Chi2;
  BuildChiEtaPathB(Eta, Z, Chi2, nullptr, true);

  double t_mid = omp_get_wtime();

  auto make_chi_as = [&](const arma::mat &Chi, TwoBodyChannel &tbc, int nK,
                         int J0) -> arma::mat {
    arma::mat ChiAS = Chi;
    for (int ibra = 0; ibra < 2 * nK; ++ibra) {
      auto ij = PqTB(tbc, ibra, nK);
      if (ij[0] < 0)
        continue;
      Orbit &oi = Z.modelspace->GetOrbit(ij[0]);
      Orbit &oj = Z.modelspace->GetOrbit(ij[1]);
      const double pij =
          Z.modelspace->phase(J0 + (oi.j2 + oj.j2) / 2);
      int ib_ex = ibra;
      if (ij[0] != ij[1])
        ib_ex = (ibra < nK) ? ibra + nK : ibra - nK;
      for (int ik = 0; ik < (int)Chi.n_cols; ++ik)
        ChiAS(ibra, ik) = Chi(ibra, ik) - pij * Chi(ib_ex, ik);
    }
    return ChiAS;
  };

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    const size_t ch_bra = ch_bra_list[ch];
    const size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;
    const int J0 = tbc_bra.J;
    if (tbc_ket.J != J0)
      continue;
    const double hatJ = std::sqrt(2.0 * J0 + 1.0);

    arma::mat Gam2(2 * nbras, 2 * nkets, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
      auto ij = PqTB(tbc_bra, ibra, nbras);
      if (ij[0] < 0)
        continue;
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = PqTB(tbc_ket, iket, nkets);
        if (kl[0] < 0)
          continue;
        Gam2(ibra, iket) =
            Gamma.TwoBody.GetTBME_J(J0, J0, ij[0], ij[1], kl[0], kl[1]);
      }
    }

    const arma::mat ChiAS_b = make_chi_as(Chi2[ch_bra], tbc_bra, nbras, J0);
    const arma::mat ChiAS_k = make_chi_as(Chi2[ch_ket], tbc_ket, nkets, J0);
    const arma::mat Z2n = -(ChiAS_b * Gam2 + Gam2 * ChiAS_k.t());

    arma::mat &Zmat = Z2.GetMatrix(ch_bra, ch_ket);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        double z = Z2n(ibra, iket) / hatJ;
        if (i == j)
          z /= PhysConst::SQRT2;
        if (k == l)
          z /= PhysConst::SQRT2;
        Zmat(ibra, iket) += z;
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler.timer[std::string(__func__) + "_ladder"] +=
        omp_get_wtime() - t_mid;
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// Gamma^III_b — leftover AMC as two terms (right first).
/// IIb and IId separately (P stripped). AMC: IIb_leftover.tex / IId_leftover.tex.
/// Z = (1−Pij)(1−Pkl)(IIb+IId). Unreduced leftover. Do not dagger χ^η.
/// RC of a combined χ is an optimization, not this gold path.
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  const double hat_lam_inv =
      1.0 / std::sqrt(2.0 * lambda + 1.0);
  auto hatJ = [](int J) { return std::sqrt(2.0 * J + 1.0); };
  auto &Z2 = Z.TwoBody;

  std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                              Z.modelspace->all_orbits.end());
  int max_j2 = 0;
  for (auto o : allorb)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(o).j2);
  const int max_J = max_j2;

  int n_orb_alloc = 0;
  for (auto o : allorb)
    n_orb_alloc = std::max(n_orb_alloc, (int)o + 1);
  const int n_orb = (int)allorb.size();

  struct ChiGIIIbTab {
    int n = 0;
    int max_J = 0;
    size_t sJ = 0, sL = 0, sK = 0, sJp = 0, sI = 0;
    std::vector<double> data;
    void allocate(int n_in, int Jmax) {
      n = n_in;
      max_J = Jmax;
      sJ = (size_t)(max_J + 1);
      sL = sJ * (size_t)n;
      sK = sL * (size_t)n;
      sJp = sK * (size_t)n;
      sI = sJp * (size_t)n;
      data.assign((size_t)n * sJp, 0.0);
    }
    double &at(index_t p, index_t q, index_t r, index_t s, int J) {
      return data[(size_t)p * sJp + (size_t)q * sK + (size_t)r * sL +
                  (size_t)s * sJ + (size_t)J];
    }
    double operator()(index_t p, index_t q, index_t r, index_t s, int J) const {
      return data[(size_t)p * sJp + (size_t)q * sK + (size_t)r * sL +
                  (size_t)s * sJ + (size_t)J];
    }
  };

  ChiGIIIbTab chi_iib, chi_iid;
  chi_iib.allocate(n_orb_alloc, max_J);
  chi_iid.allocate(n_orb_alloc, max_J);

  const double t_chi = omp_get_wtime();

  struct OmTabJJ {
    int max_J = 0;
    size_t sJ1 = 0, sL = 0, sK = 0, sJ = 0, sI = 0;
    std::vector<double> data;
    void allocate(int n, int Jmax) {
      max_J = Jmax;
      sJ1 = (size_t)(max_J + 1);
      sL = sJ1 * (size_t)(max_J + 1);
      sK = (size_t)n * sL;
      sJ = (size_t)n * sK;
      sI = (size_t)n * sJ;
      data.assign((size_t)n * sI, 0.0);
    }
    double &at(index_t p, index_t q, index_t r, index_t s, int J0, int J1) {
      return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                  (size_t)s * sL + (size_t)J0 * sJ1 + (size_t)J1];
    }
    double operator()(index_t p, index_t q, index_t r, index_t s, int J0,
                      int J1) const {
      return data[(size_t)p * sI + (size_t)q * sJ + (size_t)r * sK +
                  (size_t)s * sL + (size_t)J0 * sJ1 + (size_t)J1];
    }
  };
  OmTabJJ Om;
  Om.allocate(n_orb_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const index_t i = allorb[ii];
    Orbit &oi = Z.modelspace->GetOrbit(i);
    const double ji = oi.j2 * 0.5;
    for (auto j : allorb) {
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const double jj = oj.j2 * 0.5;
      for (auto k : allorb) {
        Orbit &ok = Z.modelspace->GetOrbit(k);
        const double jk = ok.j2 * 0.5;
        for (auto l : allorb) {
          Orbit &ol = Z.modelspace->GetOrbit(l);
          const double jl = ol.j2 * 0.5;
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(ji, jj, (double)J0))
              continue;
            for (int J1 = 0; J1 <= max_J; ++J1) {
              if (not AngMom::Triangle(jk, jl, (double)J1))
                continue;
              if (not AngMom::Triangle(J0, J1, lambda))
                continue;
              Om.at(i, j, k, l, J0, J1) =
                  Eta.TwoBody.GetTBME_J(J0, J1, (int)i, (int)j, (int)k, (int)l);
            }
          }
        }
      }
    }
  }

  std::vector<int> oj2(n_orb);
  std::vector<double> ojh(n_orb), nocc(n_orb), nnocc(n_orb);
  for (int t = 0; t < n_orb; ++t) {
    Orbit &o = Z.modelspace->GetOrbit(allorb[t]);
    oj2[t] = o.j2;
    ojh[t] = o.j2 * 0.5;
    nocc[t] = o.occ;
    nnocc[t] = 1.0 - o.occ;
  }

  // χ^IIb / χ^IId: same 4-6j leftover as orbit loops, contract (b,c) by DGEMM.
  const int n2 = n_orb * n_orb;
  const int j0_2max = max_j2 + 2 * lambda;
  const double ph_lam = Z.modelspace->phase(lambda);
  if (n2 > 0) {
    arma::mat Xb(n2, n2), Yb(n2, n2), Xd(n2, n2), Yd(n2, n2);
    for (int J7 = 0; J7 <= max_J; ++J7) {
      arma::mat accb(n2, n2, arma::fill::zeros);
      arma::mat accd(n2, n2, arma::fill::zeros);
      for (int j0_2 = 0; j0_2 <= j0_2max; ++j0_2) {
        const double j0 = 0.5 * j0_2;
        Xb.zeros();
        Yb.zeros();
        Xd.zeros();
        Yd.zeros();
#pragma omp parallel for schedule(dynamic, 1)
        for (int col = 0; col < n2; ++col) {
          const int ib = col / n_orb;
          const int ic = col - ib * n_orb;
          const int b = (int)allorb[ib];
          const int c = (int)allorb[ic];
          const int jb2 = oj2[ib], jc2 = oj2[ic];
          const double jb = ojh[ib], jc = ojh[ic];
          const double nb = nocc[ib], nnb = nnocc[ib];
          const double nc = nocc[ic], nnc = nnocc[ic];
          if ((j0_2 + jc2) % 2 != 0)
            continue;
          if (not AngMom::Triangle(jc, (double)lambda, j0))
            continue;

          for (int id = 0; id < n_orb; ++id) {
            const int d = (int)allorb[id];
            const double jd = ojh[id];
            const int jd2 = oj2[id];
            const double nd = nocc[id], nnd = nnocc[id];
            const double wb = nnb * nc * nd + nb * nnc * nnd;
            const double wd = nnc * nb * nd + nc * nnb * nnd;
            for (int ik = 0; ik < n_orb; ++ik) {
              const int k = (int)allorb[ik];
              const double jk = ojh[ik];
              const int jk2 = oj2[ik];
              const int row = id * n_orb + ik;
              if (AngMom::Triangle(jd, jk, (double)J7) and std::abs(wb) > 1e-12) {
                double s = 0.0;
                const int j2min = std::abs(jd2 - jc2) / 2;
                const int j2max = (jd2 + jc2) / 2;
                const int j3min = std::abs(jb2 - jk2) / 2;
                const int j3max = (jb2 + jk2) / 2;
                for (int J2 = j2min; J2 <= j2max; ++J2) {
                  for (int J3 = j3min; J3 <= j3max; ++J3) {
                    if (not AngMom::Triangle(J2, J3, lambda))
                      continue;
                    const double o1 = Om(d, c, b, k, J2, J3);
                    if (std::abs(o1) < 1e-16)
                      continue;
                    s += hatJ(J2) * hatJ(J3) * Z.modelspace->phase(J2) *
                         Z.modelspace->GetSixJ(J3, lambda, J2, jc, jd, j0) *
                         Z.modelspace->GetSixJ(jd, jk, J7, jb, j0, J3) * o1;
                  }
                }
                Xb(row, col) = wb * s;
              }
              if (AngMom::Triangle(ojh[ik], jd, (double)J7) and
                  std::abs(wd) > 1e-12) {
                // Xd row (j,d) with j←ik: Ω_{j c b d}
                const int j = k;
                const double jj = jk;
                const int jj2 = jk2;
                double s = 0.0;
                const int j2min = std::abs(jj2 - jc2) / 2;
                const int j2max = (jj2 + jc2) / 2;
                const int j3min = std::abs(jb2 - jd2) / 2;
                const int j3max = (jb2 + jd2) / 2;
                for (int J2 = j2min; J2 <= j2max; ++J2) {
                  for (int J3 = j3min; J3 <= j3max; ++J3) {
                    if (not AngMom::Triangle(J2, J3, lambda))
                      continue;
                    const double o1 = Om(j, c, b, d, J2, J3);
                    if (std::abs(o1) < 1e-16)
                      continue;
                    s += hatJ(J2) * hatJ(J3) * Z.modelspace->phase(J2) *
                         Z.modelspace->GetSixJ(J3, lambda, J2, jc, jj, j0) *
                         Z.modelspace->GetSixJ(jj, jd, J7, jb, j0, J3) * o1;
                  }
                }
                Xd(ik * n_orb + id, col) = wd * s;
              }
            }
          }

          for (int ia = 0; ia < n_orb; ++ia) {
            const int a = (int)allorb[ia];
            const double ja = ojh[ia];
            const int ja2 = oj2[ia];
            for (int ii = 0; ii < n_orb; ++ii) {
              const int i = (int)allorb[ii];
              const double ji = ojh[ii];
              const int ji2 = oj2[ii];
              const int row_ai = ia * n_orb + ii;
              if (AngMom::Triangle(ja, ji, (double)J7)) {
                double s = 0.0;
                const int j4min = std::abs(jb2 - ji2) / 2;
                const int j4max = (jb2 + ji2) / 2;
                const int j5min = std::abs(ja2 - jc2) / 2;
                const int j5max = (ja2 + jc2) / 2;
                for (int J4 = j4min; J4 <= j4max; ++J4) {
                  for (int J5 = j5min; J5 <= j5max; ++J5) {
                    if (not AngMom::Triangle(J4, J5, lambda))
                      continue;
                    const double o2 = Om(b, i, a, c, J4, J5);
                    if (std::abs(o2) < 1e-16)
                      continue;
                    s += hatJ(J4) * hatJ(J5) * Z.modelspace->phase(J4) *
                         Z.modelspace->GetSixJ(J4, lambda, J5, jc, ja, j0) *
                         Z.modelspace->GetSixJ(ja, ji, J7, jb, j0, J4) * o2;
                  }
                }
                Yb(row_ai, col) = s;
              }
              // Yd row (l,a) with l←ii, a←ia: Ω_{b a l c}
              if (AngMom::Triangle(ojh[ii], ja, (double)J7)) {
                const int l = i;
                const double jl = ji;
                const int jl2 = ji2;
                double s = 0.0;
                const int j4min = std::abs(jb2 - ja2) / 2;
                const int j4max = (jb2 + ja2) / 2;
                const int j5min = std::abs(jl2 - jc2) / 2;
                const int j5max = (jl2 + jc2) / 2;
                for (int J4 = j4min; J4 <= j4max; ++J4) {
                  for (int J5 = j5min; J5 <= j5max; ++J5) {
                    if (not AngMom::Triangle(J4, J5, lambda))
                      continue;
                    const double o2 = Om(b, a, l, c, J4, J5);
                    if (std::abs(o2) < 1e-16)
                      continue;
                    s += hatJ(J4) * hatJ(J5) * Z.modelspace->phase(J4) *
                         Z.modelspace->GetSixJ(J4, lambda, J5, jc, jl, j0) *
                         Z.modelspace->GetSixJ(jl, ja, J7, jb, j0, J4) * o2;
                  }
                }
                Yd(ii * n_orb + ia, col) = s;
              }
            }
          }
        }
        const double sc = (j0_2 + 1.0) * hat_lam_inv * ph_lam;
        accb += sc * (Xb * Yb.t());
        accd += sc * (Xd * Yd.t());
      }
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
      for (int id = 0; id < n_orb; ++id) {
        for (int ik = 0; ik < n_orb; ++ik) {
          const int d = (int)allorb[id];
          const int k = (int)allorb[ik];
          const int jd2 = oj2[id];
          for (int ia = 0; ia < n_orb; ++ia) {
            const int a = (int)allorb[ia];
            const int ja2 = oj2[ia];
            const double ph_ad = Z.modelspace->phase((ja2 + jd2) / 2);
            for (int ii = 0; ii < n_orb; ++ii) {
              const int i = (int)allorb[ii];
              chi_iib.at(d, k, a, i, J7) =
                  ph_ad * accb(id * n_orb + ik, ia * n_orb + ii);
              const int j = k;
              const int l = i;
              chi_iid.at(j, d, l, a, J7) =
                  accd(ik * n_orb + id, ii * n_orb + ia);
            }
          }
        }
      }
    }
  }
  Z.profiler.timer["_GIIIb_chi"] += omp_get_wtime() - t_chi;
  const double t_fold = omp_get_wtime();

  std::vector<int> compact(n_orb_alloc, -1);
  for (int t = 0; t < n_orb; ++t)
    compact[(int)allorb[t]] = t;

  ChiGIIIbTab Gam;
  Gam.allocate(n_orb_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const int i = (int)allorb[ii];
    const double ji = ojh[ii];
    for (int jj = 0; jj < n_orb; ++jj) {
      const int j = (int)allorb[jj];
      const double jjv = ojh[jj];
      for (int kk = 0; kk < n_orb; ++kk) {
        const int k = (int)allorb[kk];
        const double jk = ojh[kk];
        for (int ll = 0; ll < n_orb; ++ll) {
          const int l = (int)allorb[ll];
          const double jl = ojh[ll];
          for (int J6 = 0; J6 <= max_J; ++J6) {
            if (not AngMom::Triangle(ji, jjv, (double)J6) or
                not AngMom::Triangle(jk, jl, (double)J6))
              continue;
            Gam.at(i, j, k, l, J6) =
                Gamma.TwoBody.GetTBME_J(J6, J6, i, j, k, l);
          }
        }
      }
    }
  }

  auto iib_unred = [&](int i, int j, int k, int l, int J0) -> double {
    const int ci = compact[i], cj = compact[j], ck = compact[k], cl = compact[l];
    if (ci < 0 or cj < 0 or ck < 0 or cl < 0)
      return 0.0;
    const double ji = ojh[ci], jjv = ojh[cj];
    const double jk = ojh[ck], jl = ojh[cl];
    const int ji2 = oj2[ci], jj2 = oj2[cj], jk2 = oj2[ck], jl2 = oj2[cl];
    if (not AngMom::Triangle(ji, jjv, (double)J0) or
        not AngMom::Triangle(jk, jl, (double)J0))
      return 0.0;
    double tot = 0.0;
    for (int id = 0; id < n_orb; ++id) {
      const int d = (int)allorb[id];
      const double jd = ojh[id];
      const int jd2 = oj2[id];
      for (int ia = 0; ia < n_orb; ++ia) {
        const int a = (int)allorb[ia];
        const double ja = ojh[ia];
        const int ja2 = oj2[ia];
        const int J7min = std::max(std::abs(jd2 - jk2), std::abs(ja2 - ji2)) / 2;
        const int J7max = std::min(jd2 + jk2, ja2 + ji2) / 2;
        const int J6min = std::max(std::abs(jj2 - ja2), std::abs(jl2 - jd2)) / 2;
        const int J6max = std::min(jj2 + ja2, jl2 + jd2) / 2;
        for (int J7 = J7min; J7 <= J7max; ++J7) {
          const double chi = chi_iib(d, k, a, i, J7);
          if (std::abs(chi) < 1e-16)
            continue;
          const double hat7 = (2.0 * J7 + 1.0) * chi;
          for (int J6 = J6min; J6 <= J6max; ++J6) {
            const double g = Gam(j, a, l, d, J6);
            if (std::abs(g) < 1e-16)
              continue;
            const double six =
                Z.modelspace->GetSixJ(J0, J7, J6, jd, jl, jk) *
                Z.modelspace->GetSixJ(J0, J6, J7, ja, ji, jjv);
            if (std::abs(six) < 1e-16)
              continue;
            tot += (2.0 * J6 + 1.0) * hat7 * six * g;
          }
        }
      }
    }
    return tot;
  };

  auto iid_unred = [&](int i, int j, int k, int l, int J0) -> double {
    const int ci = compact[i], cj = compact[j], ck = compact[k], cl = compact[l];
    if (ci < 0 or cj < 0 or ck < 0 or cl < 0)
      return 0.0;
    const double ji = ojh[ci], jjv = ojh[cj];
    const double jk = ojh[ck], jl = ojh[cl];
    const int ji2 = oj2[ci], jj2 = oj2[cj], jk2 = oj2[ck], jl2 = oj2[cl];
    if (not AngMom::Triangle(ji, jjv, (double)J0) or
        not AngMom::Triangle(jk, jl, (double)J0))
      return 0.0;
    double tot = 0.0;
    for (int id = 0; id < n_orb; ++id) {
      const int d = (int)allorb[id];
      const double jd = ojh[id];
      const int jd2 = oj2[id];
      for (int ia = 0; ia < n_orb; ++ia) {
        const int a = (int)allorb[ia];
        const double ja = ojh[ia];
        const int ja2 = oj2[ia];
        const int J7min = std::max(std::abs(jj2 - jd2), std::abs(jl2 - ja2)) / 2;
        const int J7max = std::min(jj2 + jd2, jl2 + ja2) / 2;
        const int J6min = std::max(std::abs(jd2 - ji2), std::abs(ja2 - jk2)) / 2;
        const int J6max = std::min(jd2 + ji2, ja2 + jk2) / 2;
        for (int J7 = J7min; J7 <= J7max; ++J7) {
          const double chi = chi_iid(j, d, l, a, J7);
          if (std::abs(chi) < 1e-16)
            continue;
          const double hat7 = (2.0 * J7 + 1.0) * chi;
          for (int J6 = J6min; J6 <= J6max; ++J6) {
            const double g = Gam(d, i, a, k, J6);
            if (std::abs(g) < 1e-16)
              continue;
            const double six =
                Z.modelspace->GetSixJ(J0, J7, J6, jd, ji, jjv) *
                Z.modelspace->GetSixJ(J7, J6, J0, jk, jl, ja);
            if (std::abs(six) < 1e-16)
              continue;
            tot += (2.0 * J6 + 1.0) * hat7 * six * g;
          }
        }
      }
    }
    return Z.modelspace->phase((jj2 + jl2) / 2) * tot;
  };

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

  for (int ch = 0; ch < nch; ++ch) {
    const int ch_bra = (int)ch_bra_list[ch];
    const int ch_ket = (int)ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int J0 = tbc_bra.J;
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;

    arma::mat W2n(2 * nbras, 2 * nkets, arma::fill::zeros);
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
    for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto ij = PqTB(tbc_bra, ibra, nbras);
        if (ij[0] < 0)
          continue;
        auto kl = PqTB(tbc_ket, iket, nkets);
        if (kl[0] < 0)
          continue;
        W2n(ibra, iket) = iib_unred(ij[0], ij[1], kl[0], kl[1], J0) +
                          iid_unred(ij[0], ij[1], kl[0], kl[1], J0);
      }
    }

    arma::mat &Zmat = Z2.GetMatrix(ch_bra, ch_ket);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = bra.p, j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ibra_s = (i == j) ? ibra : ibra + nbras;
      const double Pij = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = ket.p, l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int iket_s = (k == l) ? iket : iket + nkets;
        const double Pkl = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J0);
        double z = W2n(ibra, iket) - Pij * W2n(ibra_s, iket) -
                   Pkl * W2n(ibra, iket_s) + Pij * Pkl * W2n(ibra_s, iket_s);
        if (i == j)
          z /= PhysConst::SQRT2;
        if (k == l)
          z /= PhysConst::SQRT2;
        Zmat(ibra, iket) += z;
      }
    }
  }

  Z.profiler.timer["_GIIIb_fold"] += omp_get_wtime() - t_fold;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_232_GIIIc(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  // χ^θ = T×T→S. Path B only: MakeNotReduced(χ_k)+MakeNotReduced(χ_j) → Pandya×Γ̄.
  // Factorized IIe/IIf χ^θ → Pandya × Γ̄.
  const int lambda = Eta.GetJRank();

  // ---- Path B: MakeNotReduced(χ_k)+MakeNotReduced(χ_j) → Pandya×Γ̄ ----
  // χ^θ = χ_k + χ_j (same slots, no transpose). Z ~ IIe-style inv Pandya.
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const double hat_lambda_inv =
      1.0 / std::sqrt(2.0 * std::max(lambda, 0) + 1.0);

  int max_j2 = 0;
  for (auto x : Z.modelspace->all_orbits)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
  std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                              Z.modelspace->all_orbits.end());
  const int n_orb = (int)allorb.size();

  // (2) RME + DGEMM → ChiTabs
  ChiTab chi_k, chi_j;
  chi_k.allocate(n_orb, max_j2);
  chi_j.allocate(n_orb, max_j2);
  FillChiThetaG3c_DGEMM(Eta, Z, chi_k, chi_j, hat_lambda_inv, lambda, allorb);

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  const int nch_eta = Z.modelspace->GetNumberTwoBodyChannels();
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  auto &Z2 = Z.TwoBody;

  // CHI_IV 2n from ChiTab (not ket-ordered Operator + AS)
  std::deque<arma::mat> CHI_IV(nch_eta);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch_eta; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    PackChiThetaFactLayout(chi_k, chi_j, Z, ch, tbc.J, CHI_IV[ch]);
  }

  // FDC ~1235 / ~2411: reuse scalar Γ̄; Pandya χ^θ in the same CC layout.
  const std::deque<arma::mat> &bar_Gamma =
      CachedBarGammaScalarCC(Gamma, Z, hGamma);
  std::deque<arma::mat> bar_CHI_IV(n_nonzero);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    bar_CHI_IV[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    if (nKets_cc < 1)
      continue;

    // Pandya of CHI_IV → bar_CHI_IV
    for (int ibra_cc = 0; ibra_cc < 2 * nKets_cc; ++ibra_cc) {
      int a, b;
      if (ibra_cc < nKets_cc) {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
        a = bra_cc.p;
        b = bra_cc.q;
      } else {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
        b = bra_cc.p;
        a = bra_cc.q;
      }
      if (ibra_cc >= nKets_cc and a == b)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      for (int iket_cc = 0; iket_cc < 2 * nKets_cc; ++iket_cc) {
        int c, d;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc);
          c = ket_cc.p;
          d = ket_cc.q;
        } else {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc.p;
          c = ket_cc.q;
        }
        if (iket_cc >= nKets_cc and c == d)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
        int Tz_J2_bc = (ob.tz2 + oc.tz2) / 2;
        int parity_J2 = (ob.l + oc.l) % 2;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        double Xbar = 0.0;
        for (int J_std = jmin; J_std <= jmax; ++J_std) {
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
          // Scalar χ^θ: χ_{ad,cb}=0 unless Tz_ad=Tz_bc (and parity). Do not
          // read CHI_IV[ch_bc](indx_ad from ch_ad, …) across channels.
          if (Tz_J2_ad != Tz_J2_bc)
            continue;
          if ((oa.l + od.l) % 2 != parity_J2)
            continue;
          int ch_J2_bc = Z.modelspace->GetTwoBodyChannelIndex(
              J_std, parity_J2, Tz_J2_bc);
          TwoBodyChannel &tbc_J2_bc = Z.modelspace->GetTwoBodyChannel(ch_J2_bc);
          int nkets_bc = tbc_J2_bc.GetNumberKets();
          if (nkets_bc < 1)
            continue;
          int indx_ad =
              tbc_J2_bc.GetLocalIndex(std::min(a, d), std::max(a, d));
          int indx_bc =
              tbc_J2_bc.GetLocalIndex(std::min(b, c), std::max(b, c));
          if (indx_ad < 0 or indx_bc < 0)
            continue;
          if (a > d)
            indx_ad += nkets_bc;
          int indx_cb = indx_bc;
          if (b > c)
            indx_bc += nkets_bc;
          if (c > b)
            indx_cb += nkets_bc;
          if (indx_ad >= (int)CHI_IV[ch_J2_bc].n_rows or
              indx_cb >= (int)CHI_IV[ch_J2_bc].n_cols)
            continue;
          Xbar -= (2 * J_std + 1) * sixj1 * CHI_IV[ch_J2_bc](indx_ad, indx_cb);
        }
        bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = Xbar;
      }
    }
  }
  for (int ch = 0; ch < nch_eta; ++ch)
    CHI_IV[ch].clear();
  CHI_IV.clear();

  // DGEMM χ̄^θ · Γ̄
  std::deque<arma::mat> bar_CHI_gamma(n_nonzero);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    if (bar_CHI_IV[ch_cc].n_rows < 1) {
      bar_CHI_gamma[ch_cc] = arma::mat(0, 0);
      continue;
    }
    bar_CHI_gamma[ch_cc] = bar_CHI_IV[ch_cc] * bar_Gamma[ch_cc];
  }
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    bar_CHI_IV[ch_cc].clear();
  bar_CHI_IV.clear();

  // Inverse Pandya (+½)
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    size_t ch_bra = ch_bra_list[ch];
    size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nKets = tbc_ket.GetNumberKets();
    if (nbras == 0 or nKets == 0)
      continue;
    int J0 = tbc_bra.J;
    for (int ibra = 0; ibra < (int)nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p, j = bra.q;
      Orbit &oi = *(bra.op), &oj = *(bra.oq);
      int ji = oi.j2, jj = oj.j2;
      int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < (int)nKets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p, l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2, jl = ol.j2;
        double commijkl = 0, commjikl = 0, commijlk = 0, commjilk = 0;

        int parity_cc = (oi.l + ol.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        int Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5,
                                              jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
          if (indx_il < 0 or indx_kj < 0)
            continue;
          int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
          int indx_li = indx_il + (l > i ? nkets_cc : 0);
          commjilk -=
              (2 * Jprime + 1) * sixj * bar_CHI_gamma[ch_cc](indx_jk, indx_li);
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          commijkl -=
              (2 * Jprime + 1) * sixj * bar_CHI_gamma[ch_cc](indx_il, indx_kj);
        }

        parity_cc = (oi.l + ok.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        Jpmin = std::max(std::abs(int(jj - jl)), std::abs(int(jk - ji))) / 2;
        Jpmax = std::min(int(jj + jl), int(jk + ji)) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5,
                                              jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_ki = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
          int indx_jl = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
          if (indx_ki < 0 or indx_jl < 0)
            continue;
          int indx_ik = indx_ki + (i > k ? nkets_cc : 0);
          int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
          commijlk -=
              (2 * Jprime + 1) * sixj * bar_CHI_gamma[ch_cc](indx_ik, indx_lj);
          indx_ki += (k > i ? nkets_cc : 0);
          indx_jl += (j > l ? nkets_cc : 0);
          commjikl -=
              (2 * Jprime + 1) * sixj * bar_CHI_gamma[ch_cc](indx_jl, indx_ki);
        }

        double zijkl =
            (commijkl - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commijlk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commjilk);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zijkl);
      }
    }
  }
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    bar_CHI_gamma[ch_cc].clear();
  bar_CHI_gamma.clear();
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}


////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_a / chi^kappa — Path B any λ (Pandya / DGEMM / inv)
///
/// T×S → tensor χ^κ. Gold: m ≡ AMC analyze ≡ this extract
///   (test_chi_kappa_pathB_vs_direct.py, test_G4a_pathB_mscheme.py).
///
/// 1. Scalar Pandya Γ (6j) + tensor Pandya Ω (9j, IMSRG adcb, no extra scale)
/// 2. VI_II DGEMM: χ̄^{J0 J1} = hΩ (−1)^{J0+J1} (occ⊙Ω̄^{J1 J0})^T Γ̄^{J1}
/// 3. InvPlus (AMC Eq4 without printed leading minus)
/// 4. W = −χΩ pair-channel DGEMM; (1−P) on W only; Hermitian W + W_klij
///      (no extra hΩ on W_klij). Store Z_unred = Z_red / Ĵ (+√2 i=j / k=l).
/// Ω is WE-reduced at every λ (λ=0 is the equal-J limit, not a Factorized fork).
////////////////////////////////////////////////////////////////////////////
namespace {
// Scalar Pandya of Γ (H). Identical 2×nKets CC layout in GIVa/GIVb.
// Arnoldi keeps H fixed, so cache by pointer across leftover diagrams / Dcoms.
struct BarGammaCCCache {
  const Operator *gamma = nullptr;
  int nch = -1;
  std::deque<arma::mat> bar;
};
BarGammaCCCache g_bar_gamma_cc;

void FillBarGammaScalarCC(const Operator &Gamma, Operator &Z, int hGamma,
                          std::deque<arma::mat> &bar_Gamma) {
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  bar_Gamma.assign(n_nonzero, arma::mat());
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    const int nKets_cc = tbc_cc.GetNumberKets();
    const int J_cc = tbc_cc.J;
    bar_Gamma[ch_cc] = arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    if (nKets_cc < 1)
      continue;
    for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc) {
      Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
      const int a = bra_cc.p, b = bra_cc.q;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      for (int iket_cc = 0; iket_cc < 2 * nKets_cc; ++iket_cc) {
        if ((iket_cc % nKets_cc) < ibra_cc)
          continue;
        int c, d;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc);
          c = ket_cc.p;
          d = ket_cc.q;
        } else {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc.p;
          c = ket_cc.q;
        }
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
        double Gammabar = 0.0;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        int dJ_std = 1;
        if ((a == d or b == c)) {
          dJ_std = 2;
          jmin += jmin % 2;
        }
        for (int J_std = jmin; J_std <= jmax; J_std += dJ_std) {
          const double sixj1 =
              Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8)
            Gammabar -= (2 * J_std + 1) * sixj1 *
                        Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
        }
        const double flip_phase =
            Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);
        if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d)) {
          bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
          if (iket_cc != ibra_cc)
            bar_Gamma[ch_cc](iket_cc, ibra_cc) = hGamma * Gammabar;
        }
        if (a != b)
          bar_Gamma[ch_cc](ibra_cc + nKets_cc,
                           (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Gammabar * flip_phase * hGamma;
        if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d))
          bar_Gamma[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc),
                           ibra_cc + nKets_cc) = Gammabar * flip_phase;
      }
    }
  }
}

const std::deque<arma::mat> &CachedBarGammaScalarCC(const Operator &Gamma,
                                                    Operator &Z, int hGamma) {
  const int nch = Z.modelspace->GetNumberTwoBodyChannels_CC();
  if (g_bar_gamma_cc.gamma == &Gamma && g_bar_gamma_cc.nch == nch &&
      (int)g_bar_gamma_cc.bar.size() == nch)
    return g_bar_gamma_cc.bar;
  FillBarGammaScalarCC(Gamma, Z, hGamma, g_bar_gamma_cc.bar);
  g_bar_gamma_cc.gamma = &Gamma;
  g_bar_gamma_cc.nch = nch;
  return g_bar_gamma_cc.bar;
}

// FDC.cc ~1243: one ket loop writes Ω̄ and occupancy nnn. Tensor Ω uses 9j
// DoTensorPandya(adcb); scalar Γ is 6j phss (cached). nnn is not a 2nd Pandya.
void FillTensor232Bars(const Operator &Eta, const Operator &Gamma, Operator &Z,
                       Tensor232Bars &B) {
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();
  const int lambda = Eta.GetJRank();
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  auto hatJ = [](int J) { return std::sqrt(2.0 * J + 1.0); };

  B.bar_Gamma = &CachedBarGammaScalarCC(Gamma, Z, hGamma);
  B.cc_pairs.clear();
  for (int ch_b = 0; ch_b < n_cc; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    if (tb.GetNumberKets() < 1)
      continue;
    for (int ch_k = 0; ch_k < n_cc; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      if (tk.GetNumberKets() < 1)
        continue;
      if ((tb.parity + tk.parity) % 2 != 0)
        continue;
      if (not CcTzCouples(tb.Tz, tk.Tz, Eta.GetTRank()))
        continue;
      if (not AngMom::Triangle(tb.J, tk.J, lambda))
        continue;
      B.cc_pairs.push_back({ch_b, ch_k});
    }
  }

  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const int np = (int)B.cc_pairs.size();
  std::vector<int> canon;
  canon.reserve((size_t)np);
  for (int ip = 0; ip < np; ++ip)
    if (B.cc_pairs[ip][0] <= B.cc_pairs[ip][1])
      canon.push_back(ip);

  std::vector<arma::mat> tmpO(np), tmpBC(np), tmpBD(np);
  const int nc = (int)canon.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ic = 0; ic < nc; ++ic) {
    const int ip = canon[ic];
    const int ch_b = B.cc_pairs[ip][0], ch_k = B.cc_pairs[ip][1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    const int Jb = tb.J, Jk = tk.J;
    arma::mat Om(2 * nb, 2 * nk, arma::fill::zeros);
    arma::mat nBC(2 * nb, 2 * nk, arma::fill::zeros);
    arma::mat nBD(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto ab = PqCC(tb, ibra, nb);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      const double n_a = oa.occ, nbar_a = 1.0 - n_a;
      const double n_b = ob.occ, nbar_b = 1.0 - n_b;
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto cd = PqCC(tk, iket, nk);
        if (cd[0] < 0)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(cd[0]);
        Orbit &od = Z.modelspace->GetOrbit(cd[1]);
        const double n_c = oc.occ, nbar_c = 1.0 - n_c;
        const double n_d = od.occ, nbar_d = 1.0 - n_d;
        const double occ_AbarBC = nbar_a * n_b * n_c + n_a * nbar_b * nbar_c;
        const double occ_ABbarD = n_a * nbar_b * n_d + nbar_a * n_b * nbar_d;
        if (std::abs(occ_AbarBC) < 1e-12 and std::abs(occ_ABbarD) < 1e-12)
          continue;
        const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
        double Xbar = 0.0;
        const int j1min = std::abs(oa.j2 - od.j2) / 2;
        const int j1max = (oa.j2 + od.j2) / 2;
        for (int J1 = j1min; J1 <= j1max; ++J1) {
          const int j2min =
              std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
          const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
          for (int J2 = j2min; J2 <= j2max; ++J2) {
            const double ninej = Z.modelspace->GetNineJ(
                ja, jd, J1, jb, jc, J2, Jb, Jk, lambda);
            if (std::abs(ninej) < 1e-14)
              continue;
            const double hats = hatJ(J1) * hatJ(J2) * hatJ(Jb) * hatJ(Jk);
            const double tbme =
                Eta.TwoBody.GetTBME_J(J1, J2, ab[0], cd[1], cd[0], ab[1]);
            Xbar -= hats *
                    Z.modelspace->phase((ob.j2 + od.j2) / 2 + Jk + J2) *
                    ninej * tbme;
          }
        }
        Om(ibra, iket) = Xbar;
        if (std::abs(occ_AbarBC) >= 1e-12)
          nBC(ibra, iket) = Xbar * occ_AbarBC;
        if (std::abs(occ_ABbarD) >= 1e-12)
          nBD(ibra, iket) = Xbar * occ_ABbarD;
      }
    }
    tmpO[ip] = std::move(Om);
    tmpBC[ip] = std::move(nBC);
    tmpBD[ip] = std::move(nBD);
  }

  std::map<CcPair, int> ip_canon;
  for (int ic = 0; ic < nc; ++ic)
    ip_canon[B.cc_pairs[canon[ic]]] = canon[ic];
  for (int ip = 0; ip < np; ++ip) {
    const int ch_b = B.cc_pairs[ip][0], ch_k = B.cc_pairs[ip][1];
    if (ch_b <= ch_k)
      continue;
    auto it = ip_canon.find({ch_k, ch_b});
    if (it == ip_canon.end())
      continue;
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    tmpO[ip] = TensorBarOmegaPartner(tmpO[it->second], tk.J, tb.J, hEta,
                                     Z.modelspace);
    tmpBC[ip] = arma::mat(2 * nb, 2 * nk, arma::fill::zeros);
    tmpBD[ip] = arma::mat(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto ab = PqCC(tb, ibra, nb);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double n_a = oa.occ, nbar_a = 1.0 - n_a;
      const double n_b = ob.occ, nbar_b = 1.0 - n_b;
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto cd = PqCC(tk, iket, nk);
        if (cd[0] < 0)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(cd[0]);
        Orbit &od = Z.modelspace->GetOrbit(cd[1]);
        const double n_c = oc.occ, nbar_c = 1.0 - n_c;
        const double n_d = od.occ, nbar_d = 1.0 - n_d;
        const double occ_AbarBC = nbar_a * n_b * n_c + n_a * nbar_b * nbar_c;
        const double occ_ABbarD = n_a * nbar_b * n_d + nbar_a * n_b * nbar_d;
        const double Xbar = tmpO[ip](ibra, iket);
        if (std::abs(occ_AbarBC) >= 1e-12)
          tmpBC[ip](ibra, iket) = Xbar * occ_AbarBC;
        if (std::abs(occ_ABbarD) >= 1e-12)
          tmpBD[ip](ibra, iket) = Xbar * occ_ABbarD;
      }
    }
  }
  B.bar_Omega.clear();
  B.nnn_AbarBC.clear();
  B.nnn_ABbarD.clear();
  for (int ip = 0; ip < np; ++ip) {
    B.bar_Omega[B.cc_pairs[ip]] = std::move(tmpO[ip]);
    B.nnn_AbarBC[B.cc_pairs[ip]] = std::move(tmpBC[ip]);
    B.nnn_ABbarD[B.cc_pairs[ip]] = std::move(tmpBD[ip]);
  }
}
} // namespace

void comm223_232_GIVa(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  Tensor232Bars B;
  comm223_232_GIVa_from_bars(Eta, Gamma, Z, B);
}

namespace {
void comm223_232_GIVa_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const Tensor232Bars &B) {
  (void)B;
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const double hat_lam_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  auto &Z2 = Z.TwoBody;
  const OrbPack P = MakeOrbPack(Z.modelspace);
  const int n_orb = P.n_orb, max_J = P.max_J;

  // χ̄^κ(il;kj) = hΩ (−1)^{J0+J1} Σ_ab occ_ABbarD(a,b,l) Ω̄^{J1 J0}(ab;il) Γ̄^{J1}(ab;kj)
  ChiTabJJ barChi;
  barChi.allocate(P.n_alloc, max_J);
  {
    ChiTabJJ barO, barG;
    barO.allocate(P.n_alloc, max_J);
    barG.allocate(P.n_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ia = 0; ia < n_orb; ++ia) {
      const int a = (int)P.allorb[ia];
      for (int ib = 0; ib < n_orb; ++ib) {
        const int b = (int)P.allorb[ib];
        for (int ic = 0; ic < n_orb; ++ic) {
          const int c = (int)P.allorb[ic];
          for (int id = 0; id < n_orb; ++id) {
            const int d = (int)P.allorb[id];
            for (int Jb = 0; Jb <= max_J; ++Jb) {
              if (not AngMom::Triangle(P.ojh[ia], P.ojh[ib], (double)Jb))
                continue;
              barG.at(a, b, c, d, Jb, Jb) =
                  PandyaGammaAdcb(Gamma, a, b, c, d, Jb);
              for (int Jk = 0; Jk <= max_J; ++Jk) {
                if (not AngMom::Triangle(Jb, Jk, lambda))
                  continue;
                if (not AngMom::Triangle(P.ojh[ic], P.ojh[id], (double)Jk))
                  continue;
                barO.at(a, b, c, d, Jb, Jk) =
                    PandyaOmegaAdcb(Eta, a, b, c, d, Jb, Jk);
              }
            }
          }
        }
      }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < n_orb; ++ii) {
      const int i = (int)P.allorb[ii];
      for (int il = 0; il < n_orb; ++il) {
        const int l = (int)P.allorb[il];
        const double nl = P.nocc[il], nbar_l = 1.0 - nl;
        for (int ik = 0; ik < n_orb; ++ik) {
          const int k = (int)P.allorb[ik];
          for (int ij = 0; ij < n_orb; ++ij) {
            const int j = (int)P.allorb[ij];
            for (int J0 = 0; J0 <= max_J; ++J0) {
              if (not AngMom::Triangle(P.ojh[ii], P.ojh[il], (double)J0))
                continue;
              for (int J1 = 0; J1 <= max_J; ++J1) {
                if (not AngMom::Triangle(J0, J1, lambda))
                  continue;
                if (not AngMom::Triangle(P.ojh[ik], P.ojh[ij], (double)J1))
                  continue;
                double sm = 0.0;
                for (int ia = 0; ia < n_orb; ++ia) {
                  const int a = (int)P.allorb[ia];
                  const double na = P.nocc[ia], nbar_a = 1.0 - na;
                  for (int ib = 0; ib < n_orb; ++ib) {
                    const int b = (int)P.allorb[ib];
                    if (not AngMom::Triangle(P.ojh[ia], P.ojh[ib],
                                             (double)J1))
                      continue;
                    const double nb = P.nocc[ib], nbar_b = 1.0 - nb;
                    const double w =
                        na * nbar_b * nl + nbar_a * nb * nbar_l;
                    if (std::abs(w) < 1e-12)
                      continue;
                    const double bo = barO(a, b, i, l, J1, J0);
                    const double bg = barG(a, b, k, j, J1, J1);
                    sm += w * bo * bg;
                  }
                }
                barChi.at(i, l, k, j, J0, J1) =
                    hEta * Z.modelspace->phase(J0 + J1) * sm;
              }
            }
          }
        }
      }
    }
  }

  ChiTabJJ Chi;
  Chi.allocate(P.n_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const int i = (int)P.allorb[ii];
    for (int ij = 0; ij < n_orb; ++ij) {
      const int j = (int)P.allorb[ij];
      for (int ib = 0; ib < n_orb; ++ib) {
        const int b = (int)P.allorb[ib];
        for (int id = 0; id < n_orb; ++id) {
          const int d = (int)P.allorb[id];
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(P.ojh[ii], P.ojh[ij], (double)J0))
              continue;
            for (int J2 = 0; J2 <= max_J; ++J2) {
              if (not AngMom::Triangle(J0, J2, lambda))
                continue;
              if (not AngMom::Triangle(P.ojh[ib], P.ojh[id], (double)J2))
                continue;
              Chi.at(i, j, b, d, J0, J2) = InvTensorPlus(
                  barChi, Z.modelspace, lambda, max_J, i, j, b, d, J0, J2);
            }
          }
        }
      }
    }
  }
  barChi.data.clear();
  barChi.data.shrink_to_fit();

  // AMC G4a_Wbra: W_red = −(−1)^{J0}/Ĵ0 λ̂^{-1} Σ (−1)^{jb+jd+λ} χ_{ijbd} Ω_{dbkl}
  auto W_red = [&](int i, int j, int k, int l, int J0) -> double {
    if (not AngMom::Triangle(Z.modelspace->GetOrbit(i).j2 * 0.5,
                             Z.modelspace->GetOrbit(j).j2 * 0.5, (double)J0) or
        not AngMom::Triangle(Z.modelspace->GetOrbit(k).j2 * 0.5,
                             Z.modelspace->GetOrbit(l).j2 * 0.5, (double)J0))
      return 0.0;
    double tot = 0.0;
    for (int ib = 0; ib < n_orb; ++ib) {
      const int b = (int)P.allorb[ib];
      for (int id = 0; id < n_orb; ++id) {
        const int d = (int)P.allorb[id];
        const double ph_bd =
            Z.modelspace->phase((P.oj2[ib] + P.oj2[id]) / 2 + lambda);
        for (int J2 = 0; J2 <= max_J; ++J2) {
          if (not AngMom::Triangle(J0, J2, lambda))
            continue;
          if (not AngMom::Triangle(P.ojh[ib], P.ojh[id], (double)J2))
            continue;
          const double ch = Chi(i, j, b, d, J0, J2);
          if (std::abs(ch) < 1e-16)
            continue;
          tot += ph_bd * hat_lam_inv * ch *
                 Eta.TwoBody.GetTBME_J(J2, J0, d, b, k, l);
        }
      }
    }
    return -Z.modelspace->phase(J0) / HatJ(J0) * tot;
  };

  auto fold_red = [&](int i, int j, int k, int l, int J0) -> double {
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    const double Pij = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J0);
    const double Pkl = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J0);
    const double w = W_red(i, j, k, l, J0);
    const double wji = W_red(j, i, k, l, J0);
    const double wkl = W_red(k, l, i, j, J0);
    const double wlk = W_red(l, k, i, j, J0);
    return (w - Pij * wji) + (wkl - Pkl * wlk);
  };

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    const size_t ch_bra = ch_bra_list[ch];
    const size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int J0 = tbc_bra.J;
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    arma::mat &Zmat = Z2.GetMatrix(ch_bra, ch_ket);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        double z = fold_red(i, j, k, l, J0) / HatJ(J0);
        if (i == j)
          z /= PhysConst::SQRT2;
        if (k == l)
          z /= PhysConst::SQRT2;
        Zmat(ibra, iket) += z;
      }
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}
} // namespace

////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_b / chi^iota — leftover AMC as two terms (right first).
/// χ^ι: Path B Pandya→DGEMM→invPlus (≡ m ≡ AMC). Do not dagger χ.
/// Folds (P stripped): W1 = Σ CG χ_aibk Ω_jbla, W2 = −hΩ Σ CG χ_akbi Ω_lbja.
/// AMC: G4b_W{1,2}_leftover.tex (unreduced leftover). Z = (1−Pij)(1−Pkl) W.
/// Do not pack χ_adbc−hZ χ_bcad with one Ω̄ (λ=0 identity Ω_jalb ≡ hΩ Ω_lbja).
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVb(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  Tensor232Bars B;
  comm223_232_GIVb_from_bars(Eta, Gamma, Z, B);
}

namespace {
void comm223_232_GIVb_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const Tensor232Bars &B) {
  (void)B;
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const double hat_lam_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  auto hatJ = [](int J) { return std::sqrt(2.0 * J + 1.0); };
  auto &Z2 = Z.TwoBody;
  const OrbPack P = MakeOrbPack(Z.modelspace);
  const int n_orb = P.n_orb, max_J = P.max_J;
  const int n_orb_alloc = P.n_alloc;
  const int max_j2 = max_J;
  const std::vector<index_t> &allorb = P.allorb;

  // χ̄^ι(il;kj) = Σ_ab occ_AbarBC(a,b,k) Γ̄^{J0}(il;ab) Ω̄^{J0 J1}(ab;kj)
  // No hΩ / (−1)^{J0+J1}. Occupancy-AS: all-orbit, not CC 2n.
  ChiTabJJ barChi;
  barChi.allocate(P.n_alloc, max_J);
  {
    ChiTabJJ barO, barG;
    barO.allocate(P.n_alloc, max_J);
    barG.allocate(P.n_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ia = 0; ia < n_orb; ++ia) {
      const int a = (int)P.allorb[ia];
      for (int ib = 0; ib < n_orb; ++ib) {
        const int b = (int)P.allorb[ib];
        for (int ic = 0; ic < n_orb; ++ic) {
          const int c = (int)P.allorb[ic];
          for (int id = 0; id < n_orb; ++id) {
            const int d = (int)P.allorb[id];
            for (int Jb = 0; Jb <= max_J; ++Jb) {
              if (not AngMom::Triangle(P.ojh[ia], P.ojh[ib], (double)Jb))
                continue;
              barG.at(a, b, c, d, Jb, Jb) =
                  PandyaGammaAdcb(Gamma, a, b, c, d, Jb);
              for (int Jk = 0; Jk <= max_J; ++Jk) {
                if (not AngMom::Triangle(Jb, Jk, lambda))
                  continue;
                if (not AngMom::Triangle(P.ojh[ic], P.ojh[id], (double)Jk))
                  continue;
                barO.at(a, b, c, d, Jb, Jk) =
                    PandyaOmegaAdcb(Eta, a, b, c, d, Jb, Jk);
              }
            }
          }
        }
      }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < n_orb; ++ii) {
      const int i = (int)P.allorb[ii];
      for (int il = 0; il < n_orb; ++il) {
        const int l = (int)P.allorb[il];
        for (int ik = 0; ik < n_orb; ++ik) {
          const int k = (int)P.allorb[ik];
          const double nk = P.nocc[ik], nbar_k = 1.0 - nk;
          for (int ij = 0; ij < n_orb; ++ij) {
            const int j = (int)P.allorb[ij];
            for (int J0 = 0; J0 <= max_J; ++J0) {
              if (not AngMom::Triangle(P.ojh[ii], P.ojh[il], (double)J0))
                continue;
              for (int J1 = 0; J1 <= max_J; ++J1) {
                if (not AngMom::Triangle(J0, J1, lambda))
                  continue;
                if (not AngMom::Triangle(P.ojh[ik], P.ojh[ij], (double)J1))
                  continue;
                double sm = 0.0;
                for (int ia = 0; ia < n_orb; ++ia) {
                  const int a = (int)P.allorb[ia];
                  const double na = P.nocc[ia], nbar_a = 1.0 - na;
                  for (int ib = 0; ib < n_orb; ++ib) {
                    const int b = (int)P.allorb[ib];
                    if (not AngMom::Triangle(P.ojh[ia], P.ojh[ib],
                                             (double)J0))
                      continue;
                    const double nb = P.nocc[ib], nbar_b = 1.0 - nb;
                    const double w =
                        nbar_a * nb * nk + na * nbar_b * nbar_k;
                    if (std::abs(w) < 1e-12)
                      continue;
                    sm += w * barG(i, l, a, b, J0, J0) *
                          barO(a, b, k, j, J0, J1);
                  }
                }
                barChi.at(i, l, k, j, J0, J1) = sm;
              }
            }
          }
        }
      }
    }
  }

  ChiTabJJ Chi;
  Chi.allocate(P.n_alloc, max_J);
  const double t_chi = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const int i = (int)P.allorb[ii];
    for (int ij = 0; ij < n_orb; ++ij) {
      const int j = (int)P.allorb[ij];
      for (int ik = 0; ik < n_orb; ++ik) {
        const int k = (int)P.allorb[ik];
        for (int il = 0; il < n_orb; ++il) {
          const int l = (int)P.allorb[il];
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(P.ojh[ii], P.ojh[ij], (double)J0))
              continue;
            for (int J1 = 0; J1 <= max_J; ++J1) {
              if (not AngMom::Triangle(J0, J1, lambda))
                continue;
              if (not AngMom::Triangle(P.ojh[ik], P.ojh[il], (double)J1))
                continue;
              Chi.at(i, j, k, l, J0, J1) = InvTensorPlus(
                  barChi, Z.modelspace, lambda, max_J, i, j, k, l, J0, J1);
            }
          }
        }
      }
    }
  }
  barChi.data.clear();
  barChi.data.shrink_to_fit();
  Z.profiler.timer["_GIVb_chi_inv"] += omp_get_wtime() - t_chi;
  const double t_left = omp_get_wtime();

  // Leftover W1/W2: same 5-6j AMC as G4b_W{1,2}_leftover.tex, but contract
  // (a,b) by DGEMM for each (J6,j0). W2n = w1 − hΩ w2, then (1−P)² + NAS.
  ChiTabJJ Om;
  Om.allocate(n_orb_alloc, max_J);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    const index_t i = allorb[ii];
    Orbit &oi = Z.modelspace->GetOrbit(i);
    const double ji = oi.j2 * 0.5;
    for (auto j : allorb) {
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const double jj = oj.j2 * 0.5;
      for (auto k : allorb) {
        Orbit &ok = Z.modelspace->GetOrbit(k);
        const double jk = ok.j2 * 0.5;
        for (auto l : allorb) {
          Orbit &ol = Z.modelspace->GetOrbit(l);
          const double jl = ol.j2 * 0.5;
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(ji, jj, (double)J0))
              continue;
            for (int J1 = 0; J1 <= max_J; ++J1) {
              if (not AngMom::Triangle(jk, jl, (double)J1))
                continue;
              if (not AngMom::Triangle(J0, J1, lambda))
                continue;
              Om.at(i, j, k, l, J0, J1) =
                  Eta.TwoBody.GetTBME_J(J0, J1, (int)i, (int)j, (int)k, (int)l);
            }
          }
        }
      }
    }
  }

  std::vector<int> compact(n_orb_alloc, -1);
  std::vector<int> oj2(n_orb);
  std::vector<double> ojh(n_orb);
  for (int t = 0; t < n_orb; ++t) {
    compact[(int)allorb[t]] = t;
    Orbit &o = Z.modelspace->GetOrbit(allorb[t]);
    oj2[t] = o.j2;
    ojh[t] = o.j2 * 0.5;
  }

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  std::vector<arma::mat> W2n_list(nch);
  for (int ch = 0; ch < nch; ++ch) {
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel((int)ch_bra_list[ch]);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel((int)ch_ket_list[ch]);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;
    W2n_list[ch] = arma::mat(2 * nbras, 2 * nkets, arma::fill::zeros);
  }

  const int n2 = n_orb * n_orb;
  const int j0_2max = max_j2 + 2 * lambda;
  if (n2 > 0) {
    arma::mat X1(n2, n2), Y1(n2, n2), X2(n2, n2), Y2(n2, n2);
    for (int J6 = 0; J6 <= max_J; ++J6) {
      arma::mat acc1(n2, n2, arma::fill::zeros);
      arma::mat acc2(n2, n2, arma::fill::zeros);
      for (int j0_2 = 0; j0_2 <= j0_2max; ++j0_2) {
        const double j0 = 0.5 * j0_2;
        X1.zeros();
        Y1.zeros();
        X2.zeros();
        Y2.zeros();
#pragma omp parallel for schedule(dynamic, 1)
        for (int col = 0; col < n2; ++col) {
          const int ia = col / n_orb;
          const int ib = col - ia * n_orb;
          const int a = (int)allorb[ia];
          const int b = (int)allorb[ib];
          const int ja2 = oj2[ia], jb2 = oj2[ib];
          const double ja = ojh[ia], jb = ojh[ib];
          if ((j0_2 + ja2) % 2 != 0)
            continue;
          if (not AngMom::Triangle(ja, (double)lambda, j0))
            continue;
          const double ph_ab =
              Z.modelspace->phase((ja2 + jb2) / 2 + lambda);

          for (int ii = 0; ii < n_orb; ++ii) {
            const int i = (int)allorb[ii];
            const double ji = ojh[ii];
            const int ji2 = oj2[ii];
            for (int kk = 0; kk < n_orb; ++kk) {
              const int k = (int)allorb[kk];
              const double jk = ojh[kk];
              const int jk2 = oj2[kk];
              if (not AngMom::Triangle(ji, jk, (double)J6))
                continue;
              const int row = ii * n_orb + kk;
              double s1 = 0.0, s2 = 0.0;
              const int j2min1 = std::abs(ja2 - ji2) / 2;
              const int j2max1 = (ja2 + ji2) / 2;
              const int j3min1 = std::abs(jb2 - jk2) / 2;
              const int j3max1 = (jb2 + jk2) / 2;
              for (int J2 = j2min1; J2 <= j2max1; ++J2) {
                for (int J3 = j3min1; J3 <= j3max1; ++J3) {
                  if (not AngMom::Triangle(J2, J3, lambda))
                    continue;
                  const double chi = Chi(a, i, b, k, J2, J3);
                  if (std::abs(chi) < 1e-16)
                    continue;
                  s1 += hatJ(J2) * hatJ(J3) *
                        Z.modelspace->GetSixJ(J3, lambda, J2, ja, ji, j0) *
                        Z.modelspace->GetSixJ(ji, jk, J6, jb, j0, J3) * chi;
                }
              }
              const int j2min2 = std::abs(ja2 - jk2) / 2;
              const int j2max2 = (ja2 + jk2) / 2;
              const int j3min2 = std::abs(jb2 - ji2) / 2;
              const int j3max2 = (jb2 + ji2) / 2;
              for (int J2 = j2min2; J2 <= j2max2; ++J2) {
                for (int J3 = j3min2; J3 <= j3max2; ++J3) {
                  if (not AngMom::Triangle(J2, J3, lambda))
                    continue;
                  const double chi = Chi(a, k, b, i, J2, J3);
                  if (std::abs(chi) < 1e-16)
                    continue;
                  s2 += hatJ(J2) * hatJ(J3) *
                        Z.modelspace->GetSixJ(J3, lambda, J2, ja, jk, j0) *
                        Z.modelspace->GetSixJ(jk, ji, J6, jb, j0, J3) * chi;
                }
              }
              X1(row, col) = ph_ab * s1;
              X2(row, col) = ph_ab * s2;
            }
          }

          for (int jj = 0; jj < n_orb; ++jj) {
            const int j = (int)allorb[jj];
            const double jjv = ojh[jj];
            const int jj2 = oj2[jj];
            for (int ll = 0; ll < n_orb; ++ll) {
              const int l = (int)allorb[ll];
              const double jl = ojh[ll];
              const int jl2 = oj2[ll];
              if (not AngMom::Triangle(jjv, jl, (double)J6))
                continue;
              const int row = jj * n_orb + ll;
              double t1 = 0.0, t2 = 0.0;
              const int j4min1 = std::abs(jj2 - jb2) / 2;
              const int j4max1 = (jj2 + jb2) / 2;
              const int j5min1 = std::abs(jl2 - ja2) / 2;
              const int j5max1 = (jl2 + ja2) / 2;
              for (int J4 = j4min1; J4 <= j4max1; ++J4) {
                for (int J5 = j5min1; J5 <= j5max1; ++J5) {
                  if (not AngMom::Triangle(J4, J5, lambda))
                    continue;
                  const double om = Om(j, b, l, a, J4, J5);
                  if (std::abs(om) < 1e-16)
                    continue;
                  t1 += hatJ(J4) * hatJ(J5) *
                        Z.modelspace->GetSixJ(J4, lambda, J5, ja, jl, j0) *
                        Z.modelspace->GetSixJ(jl, jjv, J6, jb, j0, J4) * om;
                }
              }
              const int j4min2 = std::abs(jl2 - jb2) / 2;
              const int j4max2 = (jl2 + jb2) / 2;
              const int j5min2 = std::abs(jj2 - ja2) / 2;
              const int j5max2 = (jj2 + ja2) / 2;
              for (int J4 = j4min2; J4 <= j4max2; ++J4) {
                for (int J5 = j5min2; J5 <= j5max2; ++J5) {
                  if (not AngMom::Triangle(J4, J5, lambda))
                    continue;
                  const double om = Om(l, b, j, a, J4, J5);
                  if (std::abs(om) < 1e-16)
                    continue;
                  t2 += hatJ(J4) * hatJ(J5) *
                        Z.modelspace->GetSixJ(J4, lambda, J5, ja, jjv, j0) *
                        Z.modelspace->GetSixJ(jjv, jl, J6, jb, j0, J4) * om;
                }
              }
              Y1(row, col) = t1;
              Y2(row, col) = t2;
            }
          }
        }
        const double sc = (j0_2 + 1.0) * hat_lam_inv;
        acc1 += sc * (X1 * Y1.t());
        acc2 += sc * (X2 * Y2.t());
      }

      const double hatJ6sq = 2.0 * J6 + 1.0;
#pragma omp parallel for schedule(dynamic, 1)
      for (int ch = 0; ch < nch; ++ch) {
        if (W2n_list[ch].n_rows < 1)
          continue;
        const int ch_bra = (int)ch_bra_list[ch];
        const int ch_ket = (int)ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        const int J0 = tbc_bra.J;
        const int nbras = tbc_bra.GetNumberKets();
        const int nkets = tbc_ket.GetNumberKets();
        for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
          auto ij = PqTB(tbc_bra, ibra, nbras);
          if (ij[0] < 0)
            continue;
          const int i = ij[0], j = ij[1];
          const int ci = compact[i], cj = compact[j];
          if (ci < 0 or cj < 0)
            continue;
          const double ji = ojh[ci], jjv = ojh[cj];
          const int ji2 = oj2[ci], jj2 = oj2[cj];
          for (int iket = 0; iket < 2 * nkets; ++iket) {
            auto kl = PqTB(tbc_ket, iket, nkets);
            if (kl[0] < 0)
              continue;
            const int k = kl[0], l = kl[1];
            const int ck = compact[k], cl = compact[l];
            if (ck < 0 or cl < 0)
              continue;
            const double jk = ojh[ck], jl = ojh[cl];
            const int row = ci * n_orb + ck;
            const int col = cj * n_orb + cl;
            const double six5 =
                Z.modelspace->GetSixJ(jk, jl, J0, jjv, ji, J6);
            if (std::abs(six5) < 1e-16)
              continue;
            const double pref = hatJ6sq * six5;
            const double w1 =
                -Z.modelspace->phase(J0 + (jj2 + oj2[ck]) / 2) * pref *
                acc1(row, col);
            const double w2 =
                -Z.modelspace->phase(J0 + (ji2 + oj2[cl]) / 2) * pref *
                acc2(row, col);
            W2n_list[ch](ibra, iket) += w1 - hEta * w2;
          }
        }
      }
    }
  }

  for (int ch = 0; ch < nch; ++ch) {
    if (W2n_list[ch].n_rows < 1)
      continue;
    const int ch_bra = (int)ch_bra_list[ch];
    const int ch_ket = (int)ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const int J0 = tbc_bra.J;
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    arma::mat &W2n = W2n_list[ch];
    arma::mat &Zmat = Z2.GetMatrix(ch_bra, ch_ket);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = bra.p, j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ibra_s = (i == j) ? ibra : ibra + nbras;
      const double Pij = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = ket.p, l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int iket_s = (k == l) ? iket : iket + nkets;
        const double Pkl = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J0);
        double z = W2n(ibra, iket) - Pij * W2n(ibra_s, iket) -
                   Pkl * W2n(ibra, iket_s) + Pij * Pkl * W2n(ibra_s, iket_s);
        if (i == j)
          z /= PhysConst::SQRT2;
        if (k == l)
          z /= PhysConst::SQRT2;
        Zmat(ibra, iket) += z;
      }
    }
  }

  Z.profiler.timer["_GIVb_leftover"] += omp_get_wtime() - t_left;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}
} // namespace

////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_c / chi^lambda — χ DGEMM then leftover of χ_{ialb} Ω_{bjak}.
/// Rank: T×S / S×T → tensor χ^λ. Fold T×T→S.
///
/// χ^λ is occupancy-AS (no h_χ). Default leftover (kind=0) is the CG gold:
/// WE-unpack of ChiTab × Ω, CG(λμ,λ−μ;00), ½(1−P)² in m, two leftover CGs.
/// kind=1 is the speed form: all-orbit Pandya(χ,Ω) → mid-J DGEMM → inv 6j.
////////////////////////////////////////////////////////////////////////////
static void comm223_232_GIVc_pathB(const Operator &Eta, const Operator &Gamma,
                                   Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };

  int max_j2 = 0;
  for (auto x : Z.modelspace->all_orbits)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
  const int max_J = max_j2;

  std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                              Z.modelspace->all_orbits.end());
  int n_orb_alloc = 0;
  for (auto o : allorb)
    n_orb_alloc = std::max(n_orb_alloc, (int)o + 1);

  // χ^λ table: χ(i,j,k,l; J0,J1) in GetTBME_J units (AMC chi_lambda.tex)
  ChiTabJJ Chi;
  Chi.allocate(n_orb_alloc, max_J);

  auto pq_tb = [](TwoBodyChannel &tbc, int idx, int nK) -> std::array<int, 2> {
    Ket &ket = tbc.GetKet(idx % nK);
    if (idx < nK)
      return {(int)ket.p, (int)ket.q};
    if (ket.p == ket.q)
      return {-1, -1};
    return {(int)ket.q, (int)ket.p};
  };

  const int nch_tb = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (int c0 = 0; c0 < nch_tb; ++c0) {
    TwoBodyChannel &t0 = Z.modelspace->GetTwoBodyChannel(c0);
    const int J0 = t0.J;
    const int n0 = t0.GetNumberKets();
    if (n0 < 1)
      continue;
    arma::mat G0(2 * n0, 2 * n0, arma::fill::zeros);
    for (int ib = 0; ib < 2 * n0; ++ib) {
      auto ij = pq_tb(t0, ib, n0);
      if (ij[0] < 0)
        continue;
      for (int ik = 0; ik < 2 * n0; ++ik) {
        auto ab = pq_tb(t0, ik, n0);
        if (ab[0] < 0)
          continue;
        G0(ib, ik) =
            Gamma.TwoBody.GetTBME_J(J0, J0, ij[0], ij[1], ab[0], ab[1]);
      }
    }
    for (int c1 = 0; c1 < nch_tb; ++c1) {
      TwoBodyChannel &t1 = Z.modelspace->GetTwoBodyChannel(c1);
      const int J1 = t1.J;
      const int n1 = t1.GetNumberKets();
      if (n1 < 1 or not AngMom::Triangle(J0, J1, lambda))
        continue;
      arma::mat Om(2 * n0, 2 * n1, arma::fill::zeros);
      arma::mat G1(2 * n1, 2 * n1, arma::fill::zeros);
      for (int ib = 0; ib < 2 * n0; ++ib) {
        auto r0 = pq_tb(t0, ib, n0);
        if (r0[0] < 0)
          continue;
        for (int ik = 0; ik < 2 * n1; ++ik) {
          auto r1 = pq_tb(t1, ik, n1);
          if (r1[0] < 0)
            continue;
          Om(ib, ik) =
              Eta.TwoBody.GetTBME_J(J0, J1, r0[0], r0[1], r1[0], r1[1]);
        }
      }
      for (int ib = 0; ib < 2 * n1; ++ib) {
        auto ab = pq_tb(t1, ib, n1);
        if (ab[0] < 0)
          continue;
        for (int ik = 0; ik < 2 * n1; ++ik) {
          auto kl = pq_tb(t1, ik, n1);
          if (kl[0] < 0)
            continue;
          G1(ib, ik) =
              Gamma.TwoBody.GetTBME_J(J1, J1, ab[0], ab[1], kl[0], kl[1]);
        }
      }

      // T1 = Γ^{J0} · (w(a,b,l) ⊙ Ω^{J0 J1})
      arma::mat Om_w1 = Om;
      for (int ik = 0; ik < 2 * n1; ++ik) {
        auto kl = pq_tb(t1, ik, n1);
        if (kl[0] < 0)
          continue;
        const double n_l = Z.modelspace->GetOrbit(kl[1]).occ;
        const double nbar_l = 1.0 - n_l;
        for (int ib = 0; ib < 2 * n0; ++ib) {
          auto ab = pq_tb(t0, ib, n0);
          if (ab[0] < 0) {
            Om_w1(ib, ik) = 0.0;
            continue;
          }
          const double na = Z.modelspace->GetOrbit(ab[0]).occ;
          const double nb = Z.modelspace->GetOrbit(ab[1]).occ;
          Om_w1(ib, ik) *= (1.0 - na) * (1.0 - nb) * n_l + na * nb * nbar_l;
        }
      }
      arma::mat T(2 * n0, 2 * n1, arma::fill::zeros);
      if (givc_chi_which != 2)
        T += G0 * Om_w1;

      // T2 = (w(a,b,j) ⊙ Ω^{J0 J1}) · Γ^{J1}
      if (givc_chi_which != 1) {
        arma::mat Om_w2 = Om;
        for (int ib = 0; ib < 2 * n0; ++ib) {
          auto ij = pq_tb(t0, ib, n0);
          if (ij[0] < 0)
            continue;
          const double n_j = Z.modelspace->GetOrbit(ij[1]).occ;
          const double nbar_j = 1.0 - n_j;
          for (int ik = 0; ik < 2 * n1; ++ik) {
            auto ab = pq_tb(t1, ik, n1);
            if (ab[0] < 0) {
              Om_w2(ib, ik) = 0.0;
              continue;
            }
            const double na = Z.modelspace->GetOrbit(ab[0]).occ;
            const double nb = Z.modelspace->GetOrbit(ab[1]).occ;
            Om_w2(ib, ik) *= (1.0 - na) * (1.0 - nb) * n_j + na * nb * nbar_j;
          }
        }
        T += Om_w2 * G1;
      }

      for (int ib = 0; ib < 2 * n0; ++ib) {
        auto ij = pq_tb(t0, ib, n0);
        if (ij[0] < 0)
          continue;
        for (int ik = 0; ik < 2 * n1; ++ik) {
          auto kl = pq_tb(t1, ik, n1);
          if (kl[0] < 0)
            continue;
          Chi.at(ij[0], ij[1], kl[0], kl[1], J0, J1) = T(ib, ik);
        }
      }
    }
  }

  Z.profiler.timer["_GIVc_chi"] += omp_get_wtime() - t_start;


  // Speed form = tts_ring Path B (locked run/test_givc_midj_bare_vs_cg.py):
  //   X_pqsr = Σ_ab χ_pbar Ω_aqsb  (Path A) ≡ CG bare K
  //   Path B: adcb Pandya → mid → inv (drop AMC-sample overall minus)
  //
  // Index map (AMC bar(i,j,k,l) ≡ IMSRG adcb(i,l,k,j)):
  //   mid Σ_ab λ̂^{-1}(−1)^{J_ab+λ} χ̄_adcb(p,r;a,b) Ω̄_adcb(a,b;s,q)
  //   then × (−1)^{Jp}/Ĵp on CC channel Jp=(p,r)=(s,q)
  //   inv: Ĵ0 Σ_Jp Ĵp {jr js J0; jq jp Jp} barG  (no overall minus)
  //
  // Packaging: χ/Ω WE-reduced; barG/Z reduce=true; store Z_unred=Z_red/Ĵ.
  // 2n unnormalized (GetTBME_J); ÷√2 only on AddToTBME.
  const double hat_lam_inv =
      1.0 / std::sqrt(2.0 * std::max(lambda, 0) + 1.0);
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  const int parity_eta = Eta.GetParity();
  const int rank_T = Eta.GetTRank();
  std::vector<CcPair> cc_pairs;
  for (int ch_b = 0; ch_b < n_cc; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    if (tb.GetNumberKets() < 1)
      continue;
    for (int ch_k = 0; ch_k < n_cc; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      if (tk.GetNumberKets() < 1)
        continue;
      if ((tb.parity + tk.parity) % 2 != parity_eta)
        continue;
      if (not CcTzCouples(tb.Tz, tk.Tz, rank_T))
        continue;
      if (not AngMom::Triangle(tb.J, tk.J, lambda))
        continue;
      cc_pairs.push_back({ch_b, ch_k});
    }
  }
  const int np = (int)cc_pairs.size();
  using CcKey = std::array<int, 2>;
  std::map<CcKey, arma::mat> barChi, barO;
  std::vector<arma::mat *> chi_ptr(np, nullptr), o_ptr(np, nullptr);
  for (int ip = 0; ip < np; ++ip) {
    TwoBodyChannel_CC &tb =
        Z.modelspace->GetTwoBodyChannel_CC(cc_pairs[ip][0]);
    TwoBodyChannel_CC &tk =
        Z.modelspace->GetTwoBodyChannel_CC(cc_pairs[ip][1]);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    barChi[cc_pairs[ip]] = arma::mat(2 * nb, 2 * nk, arma::fill::zeros);
    barO[cc_pairs[ip]] = arma::mat(2 * nb, 2 * nk, arma::fill::zeros);
    chi_ptr[ip] = &barChi[cc_pairs[ip]];
    o_ptr[ip] = &barO[cc_pairs[ip]];
  }
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < np; ++ip) {
    const int ch_b = cc_pairs[ip][0], ch_k = cc_pairs[ip][1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    const int Jb = tb.J, Jk = tk.J;
    arma::mat &bChi = *chi_ptr[ip];
    arma::mat &bO = *o_ptr[ip];
    // adcb: bar(a,b,c,d) ← O(a,d,c,b). (ch_il,ch_ab) → χ̄(i,l;a,b);
    // (ch_ab,ch_il) → Ω̄(a,b;k,j) for all (k,j) in ch_il.
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto pq_b = PqCC(tb, ibra, nb);
      if (pq_b[0] < 0)
        continue;
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto pq_k = PqCC(tk, iket, nk);
        if (pq_k[0] < 0)
          continue;
        bChi(ibra, iket) = PandyaChiAdcb(Chi, Z.modelspace, lambda, pq_b[0],
                                         pq_b[1], pq_k[0], pq_k[1], Jb, Jk);
        bO(ibra, iket) =
            PandyaOmegaAdcb(Eta, pq_b[0], pq_b[1], pq_k[0], pq_k[1], Jb, Jk);
      }
    }
  }

  std::deque<arma::mat> barG(n_cc);
  for (int ch = 0; ch < n_cc; ++ch) {
    TwoBodyChannel_CC &tbc = Z.modelspace->GetTwoBodyChannel_CC(ch);
    const int nK = tbc.GetNumberKets();
    barG[ch] = arma::mat(2 * nK, 2 * nK, arma::fill::zeros);
  }
  for (int ch_il = 0; ch_il < n_cc; ++ch_il) {
    TwoBodyChannel_CC &til = Z.modelspace->GetTwoBodyChannel_CC(ch_il);
    if (til.GetNumberKets() < 1)
      continue;
    // mid overall (−1)^{Jp}/Ĵp  (tts_ring Path B / test_z_ring)
    const double wL =
        Z.modelspace->phase(til.J) / std::sqrt(2.0 * til.J + 1.0);
    for (int ch_ab = 0; ch_ab < n_cc; ++ch_ab) {
      TwoBodyChannel_CC &tab = Z.modelspace->GetTwoBodyChannel_CC(ch_ab);
      if (tab.GetNumberKets() < 1)
        continue;
      if (not AngMom::Triangle(til.J, tab.J, lambda))
        continue;
      auto itC = barChi.find({ch_il, ch_ab});
      auto itO = barO.find({ch_ab, ch_il});
      if (itC == barChi.end() or itO == barO.end())
        continue;
      const double fac =
          wL * hat_lam_inv * Z.modelspace->phase(tab.J + lambda);
      barG[ch_il] += fac * (itC->second * itO->second);
    }
  }
  barChi.clear();
  barO.clear();

  auto inv_red = [&](index_t i, index_t j, index_t k, index_t l,
                     int J0) -> double {
    // ring X_pqsr with (p,q,s,r)=(i,j,k,l); barG on (i,l)×(k,j).
    // Corrected inv: +Ĵ0 Σ Ĵp sixj barG  (drop AMC-sample minus).
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    if (not AngMom::Triangle(oi.j2 * 0.5, oj.j2 * 0.5, (double)J0) or
        not AngMom::Triangle(ok.j2 * 0.5, ol.j2 * 0.5, (double)J0))
      return 0.0;
    double tot = 0.0;
    const int parity_cc = (oi.l + ol.l) % 2;
    const int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
    const int Jpmin =
        std::max(std::abs(oi.j2 - ol.j2), std::abs(ok.j2 - oj.j2)) / 2;
    const int Jpmax = std::min(oi.j2 + ol.j2, ok.j2 + oj.j2) / 2;
    for (int Jp = Jpmin; Jp <= Jpmax; ++Jp) {
      const double sixj = Z.modelspace->GetSixJ(
          ol.j2 * 0.5, ok.j2 * 0.5, (double)J0, oj.j2 * 0.5, oi.j2 * 0.5,
          Jp);
      if (std::abs(sixj) < 1e-16)
        continue;
      const size_t ch_cc =
          Z.modelspace->GetTwoBodyChannelIndex(Jp, parity_cc, Tz_cc);
      if ((int)ch_cc >= n_cc)
        continue;
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      const int nkets_cc = tbc_cc.GetNumberKets();
      if (nkets_cc < 1)
        continue;
      int indx_il = tbc_cc.GetLocalIndex(std::min((int)i, (int)l),
                                         std::max((int)i, (int)l));
      int indx_kj = tbc_cc.GetLocalIndex(std::min((int)k, (int)j),
                                         std::max((int)k, (int)j));
      if (indx_il < 0 or indx_kj < 0)
        continue;
      indx_il += (i > l ? nkets_cc : 0);
      indx_kj += (k > j ? nkets_cc : 0);
      if (indx_il >= (int)barG[ch_cc].n_rows or
          indx_kj >= (int)barG[ch_cc].n_cols)
        continue;
      tot += hat(Jp) * sixj * barG[ch_cc](indx_il, indx_kj);
    }
    return hat(J0) * tot;
  };

  auto fold_red = [&](index_t i, index_t j, index_t k, index_t l,
                      int J0) -> double {
    const double x = inv_red(i, j, k, l, J0);
    if (not givc_fold_as)
      return x;
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    const double pkl = AngMom::phase((ok.j2 + ol.j2) / 2 - J0);
    const double pij = AngMom::phase((oi.j2 + oj.j2) / 2 - J0);
    return 0.5 * (x - pkl * inv_red(i, j, l, k, J0) -
                  pij * inv_red(j, i, k, l, J0) +
                  pij * pkl * inv_red(j, i, l, k, J0));
  };

  auto &Z2p = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nchp = (int)ch_bra_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nchp; ++ch) {
    const size_t ch_bra = ch_bra_list[ch];
    const size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_ket.J != tbc_bra.J)
      continue;
    const int J0 = tbc_bra.J;
    const double hat_J0 = hat(J0);
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    for (size_t ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t p = bra.p, g = bra.q;
      const size_t ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (size_t iket = ketmin; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t q = ket.p, h = ket.q;
        double z = fold_red(p, g, q, h, J0) / hat_J0;
        if (p == g)
          z /= PhysConst::SQRT2;
        if (q == h)
          z /= PhysConst::SQRT2;
        Z2p.AddToTBME(ch_bra, ch_ket, ibra, iket, z);
      }
    }
  }

  Z.profiler.timer["comm223_232_GIVc_pathB"] += omp_get_wtime() - t_start;
}


////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVc(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  // Γ^{IV_c}: χ^λ DGEMM + leftover (tts_ring Pandya→DGEMM→inv).
  comm223_232_GIVc_pathB(Eta, Gamma, Z);
}

////////////////////////////////////////////////////////////////////////////
/// Factorized 223_132 tts: tensor ω × scalar H → scalar leftover 2b.
/// Pieces from AMC seeds in learn/amc_tts/comm_tts/output/comm223_132tts_*_seed.tex
/// Leftover packaging: unreduced Z, convention-2 λ̂^{-1}.
////////////////////////////////////////////////////////////////////////////
namespace {

double occ_ab(double na, double nb) {
  return (1.0 - na) * nb - na * (1.0 - nb);
}

struct Pack132 {
  std::vector<index_t> allorb;
  int n_orb = 0, n_alloc = 0, max_J = 0, lambda = 0;
  std::vector<int> compact, oj2;
  std::vector<double> ojh, nocc, hatJ;
  double hat_lam_inv = 1.0;
  arma::mat W_ba; // occ(a,b) * η₁(b,a)  — ladder / 1bB
  arma::mat W_ab; // occ(a,b) * η₁(a,b)  — cross
  struct WPair {
    int ia, ib;
    double w;
  };
  std::vector<WPair> wpairs_ab;
};

Pack132 MakePack132(const Operator &Eta, const Operator &Z) {
  Pack132 P;
  P.lambda = Eta.GetJRank();
  P.hat_lam_inv = 1.0 / std::sqrt(2.0 * P.lambda + 1.0);
  P.allorb.assign(Z.modelspace->all_orbits.begin(),
                  Z.modelspace->all_orbits.end());
  P.n_orb = (int)P.allorb.size();
  int max_j2 = 0;
  for (auto o : P.allorb) {
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(o).j2);
    P.n_alloc = std::max(P.n_alloc, (int)o + 1);
  }
  P.max_J = max_j2;
  P.compact.assign(P.n_alloc, -1);
  P.oj2.resize(P.n_orb);
  P.ojh.resize(P.n_orb);
  P.nocc.resize(P.n_orb);
  P.hatJ.resize(P.max_J + 1);
  for (int J = 0; J <= P.max_J; ++J)
    P.hatJ[J] = std::sqrt(2.0 * J + 1.0);
  for (int t = 0; t < P.n_orb; ++t) {
    P.compact[(int)P.allorb[t]] = t;
    Orbit &o = Z.modelspace->GetOrbit(P.allorb[t]);
    P.oj2[t] = o.j2;
    P.ojh[t] = 0.5 * o.j2;
    P.nocc[t] = o.occ;
  }
  P.W_ba.zeros(P.n_orb, P.n_orb);
  P.W_ab.zeros(P.n_orb, P.n_orb);
  for (int ia = 0; ia < P.n_orb; ++ia) {
    const int a = (int)P.allorb[ia];
    for (int ib = 0; ib < P.n_orb; ++ib) {
      if (not AngMom::Triangle(P.oj2[ia], P.oj2[ib], 2 * P.lambda))
        continue;
      const double occ = occ_ab(P.nocc[ia], P.nocc[ib]);
      if (std::abs(occ) < 1e-12)
        continue;
      const int b = (int)P.allorb[ib];
      const double eba = Eta.OneBody(b, a);
      const double eab = Eta.OneBody(a, b);
      if (std::abs(eba) > 1e-14)
        P.W_ba(ia, ib) = occ * eba;
      if (std::abs(eab) > 1e-14) {
        P.W_ab(ia, ib) = occ * eab;
        P.wpairs_ab.push_back({ia, ib, occ * eab});
      }
    }
  }
  return P;
}

double seed_cross_132(const Pack132 &P, const Operator &Eta,
                      const Operator &Gamma, ModelSpace *ms, int i, int j,
                      int k, int l, int J0) {
  const int ci = P.compact[i], cj = P.compact[j], ck = P.compact[k],
            cl = P.compact[l];
  const double ji = P.ojh[ci], jj = P.ojh[cj], jk = P.ojh[ck], jl = P.ojh[cl];
  const int ji2 = P.oj2[ci], jj2 = P.oj2[cj], jk2 = P.oj2[ck];
  const int n_orb = P.n_orb;
  double s = 0.0;
  for (const auto &pr : P.wpairs_ab) {
    const int ia = pr.ia, ib = pr.ib;
    const double w = pr.w;
    const double ja = P.ojh[ia];
    const int a = (int)P.allorb[ia];
    const int ja2 = P.oj2[ia];
    const int J3min_ka = std::abs(jk2 - ja2) / 2;
    const int J3max_ka = (jk2 + ja2) / 2;
    const double jb = P.ojh[ib];
    const int b = (int)P.allorb[ib];
    const int jb2 = P.oj2[ib];
    const int J2min_ib = std::abs(ji2 - jb2) / 2;
    const int J2max_ib = (ji2 + jb2) / 2;
    const int J4min_bj = std::abs(jb2 - jj2) / 2;
    const int J4max_bj = (jb2 + jj2) / 2;
    for (int ic = 0; ic < n_orb; ++ic) {
        const int c = (int)P.allorb[ic];
        const int jc2 = P.oj2[ic];
        const double jc = P.ojh[ic];
        const int J2min_ic = std::abs(ji2 - jc2) / 2;
        const int J2max_ic = (ji2 + jc2) / 2;
        const int J3min_kc = std::abs(jk2 - jc2) / 2;
        const int J3max_kc = (jk2 + jc2) / 2;
        const int J4min_cj = std::abs(jc2 - jj2) / 2;
        const int J4max_cj = (jc2 + jj2) / 2;
        // Term1: η_icka^{J2 J3} Γ_bjcl^{J4} — J2=(i,c), J3=(k,a), J4=(b,j)
        for (int J2 = J2min_ic; J2 <= J2max_ic; ++J2) {
          for (int J3 = J3min_ka; J3 <= J3max_ka; ++J3) {
            if (not AngMom::Triangle(J2, J3, P.lambda))
              continue;
            const double sixj1 = ms->GetSixJ((double)J2, (double)J3,
                                             (double)P.lambda, ja, jb, jk);
            if (std::abs(sixj1) < 1e-12)
              continue;
            const double eta_icka = Eta.TwoBody.GetTBME_J(J2, J3, i, c, k, a);
            if (std::abs(eta_icka) < 1e-14)
              continue;
            for (int J4 = J4min_bj; J4 <= J4max_bj; ++J4) {
              const double g_bjcl = Gamma.TwoBody.GetTBME_J(J4, J4, b, j, c, l);
              if (std::abs(g_bjcl) < 1e-14)
                continue;
              const double n9 = ms->GetNineJ(ji, jc, (double)J2, jj, (double)J4,
                                             jb, (double)J0, jl, jk);
              if (std::abs(n9) < 1e-12)
                continue;
              const int ph =
                  ms->phase(J0 + jk2 / 2) * ms->phase(J4 + jc2 / 2);
              s += ph * w * P.hatJ[J2] * P.hatJ[J3] * (2.0 * J4 + 1.0) *
                   P.hat_lam_inv * sixj1 * n9 * eta_icka * g_bjcl;
            }
          }
        }
        // Term2: η_ibkc^{J2 J3} Γ_cjal^{J4} — J2=(i,b), J3=(k,c), J4=(c,j)
        for (int J2 = J2min_ib; J2 <= J2max_ib; ++J2) {
          for (int J3 = J3min_kc; J3 <= J3max_kc; ++J3) {
            if (not AngMom::Triangle(J2, J3, P.lambda))
              continue;
            const double sixj2 = ms->GetSixJ((double)J3, (double)J2,
                                             (double)P.lambda, jb, ja, ji);
            if (std::abs(sixj2) < 1e-12)
              continue;
            const double eta_ibkc = Eta.TwoBody.GetTBME_J(J2, J3, i, b, k, c);
            if (std::abs(eta_ibkc) < 1e-14)
              continue;
            for (int J4 = J4min_cj; J4 <= J4max_cj; ++J4) {
              const double g_cjal = Gamma.TwoBody.GetTBME_J(J4, J4, c, j, a, l);
              if (std::abs(g_cjal) < 1e-14)
                continue;
              const double n9 = ms->GetNineJ(ji, ja, (double)J3, jj, (double)J4,
                                             jc, (double)J0, jl, jk);
              if (std::abs(n9) < 1e-12)
                continue;
              const int ph =
                  ms->phase(J0 + ji2 / 2) *
                  ms->phase(J2 + J3 + J4 + (ja2 + jb2 + jc2) / 2);
              s += ph * w * P.hatJ[J2] * P.hatJ[J3] * (2.0 * J4 + 1.0) *
                   P.hat_lam_inv * sixj2 * n9 * eta_ibkc * g_cjal;
            }
          }
        }
      }
  }
  return s;
}

} // namespace

void comm223_132_tts_ladder(const Operator &Eta, const Operator &Gamma,
                            Operator &Z) {
  // AMC seed: learn/amc_tts/comm_tts/output/comm223_132tts_ladder_seed.tex
  //   T1: Ĵ0^{-1} (−1)^{J2+ja+jc} Ĵ2 λ̂^{-1} {J0 J2 λ; ja jb jc} η_ba η_cakl Γ_ijcb
  //   T2: −(−1)^{J0} Ĵ0^{-1} (−1)^{ja+jc} Ĵ2 λ̂^{-1} {J0 J2 λ; jb ja jc} η_ba η_ijcb Γ_cakl
  // Convention 2 leftover: AMC λ̂^{-1}. No extra (1−P).
  // Per channel: for each dummy c, Γ tables G(ij,b) and Gg(kl,a) then (a,b) DGEMM in J2.
  Z.modelspace->PreCalculateSixJ();
  const Pack132 P = MakePack132(Eta, Z);
  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  const int n_orb = P.n_orb;

  for (int ich = 0; ich < nch; ich++) {
    size_t ch_bra = ch_bra_list[ich];
    size_t ch_ket = ch_ket_list[ich];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int J0 = tbc_bra.J;
    const double hatJ0inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;

    std::vector<int> bra_i(nbras), bra_j(nbras), ket_k(nkets), ket_l(nkets);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      bra_i[ibra] = (int)bra.p;
      bra_j[ibra] = (int)bra.q;
    }
    for (int iket = 0; iket < nkets; ++iket) {
      Ket &ket = tbc_ket.GetKet(iket);
      ket_k[iket] = (int)ket.p;
      ket_l[iket] = (int)ket.q;
    }

    arma::mat Zacc(nbras, nkets, arma::fill::zeros);
#pragma omp parallel
    {
      arma::mat loc(nbras, nkets, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ic = 0; ic < n_orb; ++ic) {
        const int c = (int)P.allorb[ic];
        const int jc2 = P.oj2[ic];
        const double jc = P.ojh[ic];

        arma::mat G(nbras, n_orb, arma::fill::zeros);
        arma::mat Gg(nkets, n_orb, arma::fill::zeros);
        for (int ibra = 0; ibra < nbras; ++ibra) {
          for (int ib = 0; ib < n_orb; ++ib)
            G(ibra, ib) = Gamma.TwoBody.GetTBME_J(
                J0, J0, bra_i[ibra], bra_j[ibra], c, (int)P.allorb[ib]);
        }
        for (int iket = 0; iket < nkets; ++iket) {
          for (int ia = 0; ia < n_orb; ++ia)
            Gg(iket, ia) = Gamma.TwoBody.GetTBME_J(
                J0, J0, c, (int)P.allorb[ia], ket_k[iket], ket_l[iket]);
        }

        for (int J2 = 0; J2 <= P.max_J; ++J2) {
          if (not AngMom::Triangle((double)J0, (double)J2, (double)P.lambda))
            continue;
          const double pref0 = hatJ0inv * P.hatJ[J2] * P.hat_lam_inv;

          arma::mat PM(n_orb, n_orb, arma::fill::zeros);
          bool anyP = false;
          for (int ia = 0; ia < n_orb; ++ia) {
            if (not AngMom::Triangle(jc2, P.oj2[ia], 2 * J2))
              continue;
            const int ph_a = Z.modelspace->phase(J2 + (P.oj2[ia] + jc2) / 2);
            const double ja = P.ojh[ia];
            for (int ib = 0; ib < n_orb; ++ib) {
              const double w = P.W_ba(ia, ib);
              if (std::abs(w) < 1e-16)
                continue;
              const double sixj = AngMom::SixJ((double)J0, (double)J2,
                                               (double)P.lambda, ja, P.ojh[ib], jc);
              if (std::abs(sixj) < 1e-12)
                continue;
              PM(ia, ib) = w * pref0 * ph_a * sixj;
              anyP = true;
            }
          }
          if (anyP) {
            arma::mat E(nkets, n_orb, arma::fill::zeros);
            for (int iket = 0; iket < nkets; ++iket) {
              for (int ia = 0; ia < n_orb; ++ia) {
                if (not AngMom::Triangle(jc2, P.oj2[ia], 2 * J2))
                  continue;
                E(iket, ia) = Eta.TwoBody.GetTBME_J(
                    J2, J0, c, (int)P.allorb[ia], ket_k[iket], ket_l[iket]);
              }
            }
            loc += G * PM.t() * E.t();
          }

          arma::mat P2(n_orb, n_orb, arma::fill::zeros);
          bool anyP2 = false;
          for (int ib = 0; ib < n_orb; ++ib) {
            if (not AngMom::Triangle(jc2, P.oj2[ib], 2 * J2))
              continue;
            const double jb = P.ojh[ib];
            for (int ia = 0; ia < n_orb; ++ia) {
              const double w = P.W_ba(ia, ib);
              if (std::abs(w) < 1e-16)
                continue;
              const double sixj = AngMom::SixJ((double)J0, (double)J2,
                                               (double)P.lambda, jb, P.ojh[ia], jc);
              if (std::abs(sixj) < 1e-12)
                continue;
              const int ph = Z.modelspace->phase(J0 + (P.oj2[ia] + jc2) / 2);
              P2(ia, ib) = w * pref0 * ph * sixj;
              anyP2 = true;
            }
          }
          if (anyP2) {
            arma::mat Om(nbras, n_orb, arma::fill::zeros);
            for (int ibra = 0; ibra < nbras; ++ibra) {
              for (int ib = 0; ib < n_orb; ++ib) {
                if (not AngMom::Triangle(jc2, P.oj2[ib], 2 * J2))
                  continue;
                Om(ibra, ib) = Eta.TwoBody.GetTBME_J(
                    J0, J2, bra_i[ibra], bra_j[ibra], c, (int)P.allorb[ib]);
              }
            }
            loc -= Om * P2.t() * Gg.t();
          }
        }
      }
#pragma omp critical
      Zacc += loc;
    }

    for (int ibra = 0; ibra < nbras; ++ibra) {
      const int i = bra_i[ibra], j = bra_j[ibra];
      const int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < nkets; ++iket) {
        double zijkl = Zacc(ibra, iket);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (ket_k[iket] == ket_l[iket])
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

void comm223_132_tts_onebody(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {
  const Pack132 P = MakePack132(Eta, Z);
  const int lambda = P.lambda;
  auto &Z2 = Z.TwoBody;
  const int n_orb = P.n_orb;

  // χ^η scalar 1b from AMC 1bA (η2 × η1), then (1-P) insert into scalar Γ.
  arma::mat CHI = Gamma.OneBody * 0;
  std::vector<index_t> porb(Z.modelspace->all_orbits.begin(),
                            Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < (int)porb.size(); ++ip) {
    auto p = porb[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);
    const double jp = 0.5 * op.j2;
    const int ph_p = Z.modelspace->phase(op.j2 / 2);
    const double diag = 1.0 / (op.j2 + 1.0);
    for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
      Orbit &oq = Z.modelspace->GetOrbit(q);
      double chi = 0.0;
      for (int ia = 0; ia < n_orb; ++ia) {
        const int a = (int)P.allorb[ia];
        const double ja = P.ojh[ia];
        const int J2min = std::abs(op.j2 - P.oj2[ia]) / 2;
        const int J2max = (op.j2 + P.oj2[ia]) / 2;
        for (int ib = 0; ib < n_orb; ++ib) {
          const double w = P.W_ba(ia, ib);
          if (std::abs(w) < 1e-16)
            continue;
          const int b = (int)P.allorb[ib];
          const double jb = P.ojh[ib];
          const int J3min = std::abs(oq.j2 - P.oj2[ib]) / 2;
          const int J3max = (oq.j2 + P.oj2[ib]) / 2;
          for (int J2 = J2min; J2 <= J2max; J2++) {
            for (int J3 = J3min; J3 <= J3max; J3++) {
              if (not AngMom::Triangle(J3, lambda, J2))
                continue;
              const double sixj = AngMom::SixJ(
                  (double)J3, (double)lambda, (double)J2, ja, jp, jb);
              if (std::abs(sixj) < 1e-12)
                continue;
              const int ph = Z.modelspace->phase(J2 + P.oj2[ia] / 2);
              const double eta2 = Eta.TwoBody.GetTBME_J(J2, J3, p, a, q, b);
              chi -= ph_p * diag * w * ph * P.hatJ[J2] * P.hatJ[J3] *
                     P.hat_lam_inv * sixj * eta2;
            }
          }
        }
      }
      CHI(p, q) = chi;
    }
  }

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

  for (int ich = 0; ich < nch; ich++) {
    size_t ch_bra = ch_bra_list[ich];
    size_t ch_ket = ch_ket_list[ich];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int J0 = tbc_bra.J;
    const double hatJ0inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;

    std::vector<int> bra_i(nbras), bra_j(nbras), ket_k(nkets), ket_l(nkets);
    std::vector<int> bra_i2(nbras), bra_j2(nbras), ket_k2(nkets), ket_l2(nkets);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      bra_i[ibra] = (int)bra.p;
      bra_j[ibra] = (int)bra.q;
      bra_i2[ibra] = Z.modelspace->GetOrbit(bra.p).j2;
      bra_j2[ibra] = Z.modelspace->GetOrbit(bra.q).j2;
    }
    for (int iket = 0; iket < nkets; ++iket) {
      Ket &ket = tbc_ket.GetKet(iket);
      ket_k[iket] = (int)ket.p;
      ket_l[iket] = (int)ket.q;
      ket_k2[iket] = Z.modelspace->GetOrbit(ket.p).j2;
      ket_l2[iket] = Z.modelspace->GetOrbit(ket.q).j2;
    }

    arma::mat ZB(nbras, nkets, arma::fill::zeros);
    arma::mat ZBkl(nbras, nkets, arma::fill::zeros);
    arma::mat ZB2(nbras, nkets, arma::fill::zeros);
    arma::mat ZB2ij(nbras, nkets, arma::fill::zeros);

    // 1bB: F(k,l,c,J2) independent of bras; ZB.col(iket) += F * η(ij; k,c)
#pragma omp parallel
    {
      arma::mat locB(nbras, nkets, arma::fill::zeros);
      arma::mat locBkl(nbras, nkets, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int iket = 0; iket < nkets; ++iket) {
        const int k = ket_k[iket], l = ket_l[iket];
        const int ck = P.compact[k], cl = P.compact[l];
        const double jk = P.ojh[ck], jl = P.ojh[cl];
        const int jk2 = P.oj2[ck], jl2 = P.oj2[cl];
        const int ph_k = Z.modelspace->phase(J0 + jk2 / 2);
        const int ph_l = Z.modelspace->phase(J0 + jl2 / 2);
        for (int ic = 0; ic < n_orb; ++ic) {
          const int c = (int)P.allorb[ic];
          const int jc2 = P.oj2[ic];
          const double jc = P.ojh[ic];
          auto accum_B = [&](int kk, int ll, double jkk, double jll, int jkk2,
                             int ph_outer, arma::mat &loc) {
            const int J2min = std::abs(jkk2 - jc2) / 2;
            const int J2max = (jkk2 + jc2) / 2;
            for (int J2 = J2min; J2 <= J2max; ++J2) {
              if (not AngMom::Triangle((double)J0, (double)J2, (double)lambda))
                continue;
              const double sixj2 = AngMom::SixJ((double)J2, (double)lambda,
                                                (double)J0, jll, jkk, jc);
              if (std::abs(sixj2) < 1e-12)
                continue;
              double F = 0.0;
              for (int ia = 0; ia < n_orb; ++ia) {
                const int J3min = std::abs(jc2 - P.oj2[ia]) / 2;
                const int J3max = (jc2 + P.oj2[ia]) / 2;
                const double ja = P.ojh[ia];
                const int a = (int)P.allorb[ia];
                for (int ib = 0; ib < n_orb; ++ib) {
                  const double w = P.W_ba(ia, ib);
                  if (std::abs(w) < 1e-16)
                    continue;
                  const double jb = P.ojh[ib];
                  const int b = (int)P.allorb[ib];
                  for (int J3 = J3min; J3 <= J3max; ++J3) {
                    const double sixj3 = AngMom::SixJ(
                        jll, jc, (double)lambda, ja, jb, (double)J3);
                    if (std::abs(sixj3) < 1e-12)
                      continue;
                    const int ph = Z.modelspace->phase(J3 + P.oj2[ib] / 2);
                    F += w * ph * (2.0 * J3 + 1.0) * sixj3 *
                         Gamma.TwoBody.GetTBME_J(J3, J3, c, a, ll, b);
                  }
                }
              }
              const double scale = -ph_outer * hatJ0inv * P.hatJ[J2] *
                                   P.hat_lam_inv * sixj2 * F;
              if (std::abs(scale) < 1e-16)
                continue;
              for (int ibra = 0; ibra < nbras; ++ibra)
                loc(ibra, iket) +=
                    scale * Eta.TwoBody.GetTBME_J(J0, J2, bra_i[ibra],
                                                  bra_j[ibra], kk, c);
            }
          };
          accum_B(k, l, jk, jl, jk2, ph_k, locB);
          accum_B(l, k, jl, jk, jl2, ph_l, locBkl);
        }
      }
#pragma omp critical
      {
        ZB += locB;
        ZBkl += locBkl;
      }
    }

    // 1bB2: inner(a,b,J2) independent of kets; ZB2.row(ibra) += G * η(i,c; kl)
#pragma omp parallel
    {
      arma::mat loc2(nbras, nkets, arma::fill::zeros);
      arma::mat loc2ij(nbras, nkets, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ibra = 0; ibra < nbras; ++ibra) {
        const int i = bra_i[ibra], j = bra_j[ibra];
        const int ci = P.compact[i], cj = P.compact[j];
        const double ji = P.ojh[ci], jj = P.ojh[cj];
        const int ji2 = P.oj2[ci], jj2 = P.oj2[cj];
        const int ph_i = Z.modelspace->phase(ji2 / 2);
        const int ph_j = Z.modelspace->phase(jj2 / 2);
        auto accum_B2 = [&](int ii, int jjv, double jii, double jjj, int jii2,
                            int jjv2, int ph_outer, arma::mat &loc) {
          for (int ic = 0; ic < n_orb; ++ic) {
            const int c = (int)P.allorb[ic];
            const int jc2 = P.oj2[ic];
            const double jc = P.ojh[ic];
            const int J3min = std::abs(jii2 - jc2) / 2;
            const int J3max = (jii2 + jc2) / 2;
            for (int J3 = J3min; J3 <= J3max; ++J3) {
              if (not AngMom::Triangle((double)J0, (double)J3, (double)lambda))
                continue;
              const double sixj3 = AngMom::SixJ(
                  (double)J0, (double)J3, (double)lambda, jc, jjj, jii);
              if (std::abs(sixj3) < 1e-12)
                continue;
              double inner = 0.0;
              for (int ia = 0; ia < n_orb; ++ia) {
                const int J2min = std::abs(jjv2 - P.oj2[ia]) / 2;
                const int J2max = (jjv2 + P.oj2[ia]) / 2;
                const double ja = P.ojh[ia];
                const int a = (int)P.allorb[ia];
                for (int ib = 0; ib < n_orb; ++ib) {
                  const double w = P.W_ba(ia, ib);
                  if (std::abs(w) < 1e-16)
                    continue;
                  const double jb = P.ojh[ib];
                  const int b = (int)P.allorb[ib];
                  for (int J2 = J2min; J2 <= J2max; ++J2) {
                    const double sixj2 = AngMom::SixJ(
                        jc, jjj, (double)lambda, ja, jb, (double)J2);
                    if (std::abs(sixj2) < 1e-12)
                      continue;
                    const int ph =
                        Z.modelspace->phase(J2 + J3 + P.oj2[ib] / 2);
                    inner += w * ph * (2.0 * J2 + 1.0) * sixj2 *
                             Gamma.TwoBody.GetTBME_J(J2, J2, jjv, a, c, b);
                  }
                }
              }
              const double scale = -ph_outer * hatJ0inv * P.hatJ[J3] *
                                   P.hat_lam_inv * sixj3 * inner;
              if (std::abs(scale) < 1e-16)
                continue;
              for (int iket = 0; iket < nkets; ++iket)
                loc(ibra, iket) +=
                    scale * Eta.TwoBody.GetTBME_J(J3, J0, ii, c, ket_k[iket],
                                                  ket_l[iket]);
            }
          }
        };
        accum_B2(i, j, ji, jj, ji2, jj2, ph_i, loc2);
        accum_B2(j, i, jj, ji, jj2, ji2, ph_j, loc2ij);
      }
#pragma omp critical
      {
        ZB2 += loc2;
        ZB2ij += loc2ij;
      }
    }

    for (int ibra = 0; ibra < nbras; ibra++) {
      const int i = bra_i[ibra], j = bra_j[ibra];
      const int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < nkets; iket++) {
        const int k = ket_k[iket], l = ket_l[iket];
        double zijkl = 0.0;
        for (int ia = 0; ia < n_orb; ++ia) {
          const int a = (int)P.allorb[ia];
          const double xajkl = Gamma.TwoBody.GetTBME_J(J0, J0, a, j, k, l);
          const double xiakl = Gamma.TwoBody.GetTBME_J(J0, J0, i, a, k, l);
          const double xijal = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, a, l);
          const double xijka = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, k, a);
          zijkl += CHI(i, a) * xajkl + CHI(j, a) * xiakl - xijal * CHI(a, k) -
                   xijka * CHI(a, l);
        }
        const int ph_kl =
            Z.modelspace->phase((ket_k2[iket] + ket_l2[iket]) / 2 - J0);
        const int ph_ij =
            Z.modelspace->phase((bra_i2[ibra] + bra_j2[ibra]) / 2 - J0);
        // 1bA 122 already matches m χ^η unpack. The overall minus was locked
        // at emax=1 where 1bA vanishes; it belongs only on 1bB/1bB2.
        zijkl -= (ZB(ibra, iket) - ph_kl * ZBkl(ibra, iket)) -
                 (ZB2(ibra, iket) - ph_ij * ZB2ij(ibra, iket));
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

void comm223_132_tts_cross(const Operator &Eta, const Operator &Gamma,
                            Operator &Z) {
  // AMC seed (Path A): learn/amc_tts/comm_tts/output/comm223_132tts_cross_seed.tex
  //   Term1: −(−1)^{J0+jk}(−1)^{J4+jc} Ĵ2 Ĵ3 Ĵ4² λ̂^{-1} {J2 J3 λ; ja jb jk}
  //          {ji jc J2; jj J4 jb; J0 jl jk} η_ab η_icka Γ_bjcl
  //   Term2: −(−1)^{J0+ji}(−1)^{J2+J3+J4+ja+jb+jc} Ĵ2 Ĵ3 Ĵ4² λ̂^{-1} {J3 J2 λ; jb ja ji}
  //          {ji ja J3; jj J4 jc; J0 jl jk} η_ab η_ibkc Γ_cjal
  // Restore (1−Pij)(1−Pkl). Convention 2: AMC λ̂^{-1}.
  // AMC overall minus on both terms dropped: seed is (icka Γ − ibkc Γ).
  // Skip order kept: 6j/η then 9j. W_ab = occ·η₁.
  Z.modelspace->PreCalculateSixJ();
  const Pack132 P = MakePack132(Eta, Z);
  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ich++) {
    size_t ch_bra = ch_bra_list[ich];
    size_t ch_ket = ch_ket_list[ich];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.J != tbc_ket.J)
      continue;
    const int J0 = tbc_bra.J;
    const int nbras = tbc_bra.GetNumberKets();
    const int nkets = tbc_ket.GetNumberKets();
    for (int ibra = 0; ibra < nbras; ibra++) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < nkets; iket++) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int ph_ij = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J0);
        const int ph_kl = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J0);
        const int ph_ijkl =
            Z.modelspace->phase((oi.j2 + oj.j2 + ok.j2 + ol.j2) / 2);
        double zijkl =
            seed_cross_132(P, Eta, Gamma, Z.modelspace, i, j, k, l, J0) -
            ph_ij * seed_cross_132(P, Eta, Gamma, Z.modelspace, j, i, k, l, J0) -
            ph_kl * seed_cross_132(P, Eta, Gamma, Z.modelspace, i, j, l, k, J0) +
            ph_ijkl * seed_cross_132(P, Eta, Gamma, Z.modelspace, j, i, l, k, J0);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

void comm223_132_tts(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  const double t_start = omp_get_wtime();
  if (Z.GetJRank() != 0)
    return;

  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  double t0 = omp_get_wtime();
  comm223_132_tts_ladder(Eta, Gamma, Z);
  Z.profiler.timer["_132_ladder"] += omp_get_wtime() - t0;
  t0 = omp_get_wtime();
  comm223_132_tts_onebody(Eta, Gamma, Z);
  Z.profiler.timer["_132_onebody"] += omp_get_wtime() - t0;
  t0 = omp_get_wtime();
  comm223_132_tts_cross(Eta, Gamma, Z);
  Z.profiler.timer["_132_cross"] += omp_get_wtime() - t0;

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

} // namespace FactorizedDoubleCommutator_eths

} // namespace Commutator
