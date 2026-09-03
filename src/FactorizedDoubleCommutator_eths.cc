
#include "FactorizedDoubleCommutator_eths.hh"
#include "Commutator.hh"
#include "ReferenceImplementations.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <map>
#include <unordered_map>
#include <cstdint>
#include <array>
#include <iostream>
#include <iomanip>
#include <deque>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
namespace Commutator {

namespace FactorizedDoubleCommutator_eths {

void comm223_231_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z);
void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z);

bool use_goose_tank_1b = true; // always include
bool use_goose_tank_2b = true; // always include

bool use_1b_intermediates = true;
bool use_2b_intermediates = true;
bool use_goose_tank_only_1b = false; // only calculate Goose Tanks
bool use_goose_tank_only_2b = false; // only calculate Goose Tanks
bool use_TypeI_1b = true;
bool use_TypeII_1b = true;
bool use_TypeIII_1b = true;
bool use_TypeIIIa_1b = true; // Pandya+DGEMM f^III_a (paper B4c/B5c)
bool use_TypeIIIa_slow = false; // f^III_a W1/W2 oracle (binding stub)
bool use_TypeGI_2b = true;  // Gamma^I / chi^epsilon
bool use_TypeGII_2b = true; // Gamma^II / chi^zeta (universal AMC Path B)
bool use_TypeGIIIa_2b = true; // Gamma^III_a Path B χ^η (2n×2n) → Chi_AS×Γ
bool use_TypeGIIIb_2b = true; // Gamma^III_b Path B Fac Pandya χ → RC → DGEMM
bool use_TypeGIIIc_2b = true; // Gamma^III_c Factorized IIe/IIf (λ=0)
bool use_TypeGIVa_2b = true; // Gamma^IV_a / chi^kappa Factorized DGEMM (λ=0)
bool use_TypeGIVb_2b = true;  // Gamma^IV_b / chi^iota Factorized DGEMM (λ=0)
bool use_TypeGIVc_2b = true;  // Gamma^IV_c / chi^lambda (enable diagram)
bool use_TypeII_2b = true;
bool use_TypeIII_2b = true;
bool use_GT_TypeI_2b = true;
bool use_GT_TypeIV_2b = true;
//  bool SlowVersion = false;

void SetUse_GooseTank_1b(bool tf) { use_goose_tank_1b = tf; }

void SetUse_GooseTank_2b(bool tf) { use_goose_tank_2b = tf; }

void SetUse_1b_Intermediates(bool tf) { use_1b_intermediates = tf; }
void SetUse_2b_Intermediates(bool tf) { use_2b_intermediates = tf; }

void SetUse_GooseTank_only_1b(bool tf) { use_goose_tank_only_1b = tf; }

void SetUse_GooseTank_only_2b(bool tf) { use_goose_tank_only_2b = tf; }

void SetUse_TypeI_1b(bool tf) { use_TypeI_1b = tf; }

void SetUse_TypeII_1b(bool tf) { use_TypeII_1b = tf; }

void SetUse_TypeIII_1b(bool tf) { use_TypeIII_1b = tf; }

void SetUse_TypeIIIa_1b(bool tf) { use_TypeIIIa_1b = tf; }
void SetUse_TypeIIIa_slow(bool tf) { use_TypeIIIa_slow = tf; }

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
void SetUse_TypeII_2b(bool tf) { use_TypeII_2b = tf; }

void SetUse_TypeIII_2b(bool tf) { use_TypeIII_2b = tf; }

void SetUse_GT_TypeI_2b(bool tf) { use_GT_TypeI_2b = tf; }

void SetUse_GT_TypeIV_2b(bool tf) { use_GT_TypeIV_2b = tf; }

//  void UseSlowVersion(bool tf)
//  {
//    SlowVersion = tf;
//  }

// factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
// Eta is tensor (reduced); Gamma and Z are scalar.
void comm223_231_st(const Operator &Eta, const Operator &Gamma, Operator &Z) {

  const Operator *Etap = &Eta;
  const Operator *Gammap = &Gamma;
  Operator Gammatmp;

  // Tensor operators stay reduced (MakeNotReduced refuses rank_J>0).
  // Scalar Gamma may be un-reduced for consistency with scalar factorized code.
  if (Gamma.IsReduced() and Gamma.GetJRank() == 0) {
    Gammatmp = Gamma;
    Gammatmp.MakeNotReduced();
    Gammap = &Gammatmp;
  }

  bool z_was_reduced = Z.IsReduced();
  if (z_was_reduced and Z.GetJRank() == 0)
    Z.MakeNotReduced();

  if (use_1b_intermediates)
    comm223_231_chi1b_tensor(*Etap, *Gammap, Z);
  if (use_2b_intermediates)
    comm223_231_chi2b_tensor(*Etap, *Gammap, Z);

  if (z_was_reduced and Z.GetJRank() == 0)
    Z.MakeReduced();

  return;
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
  //
  // Gold chain (locked): m-scheme ≡ AMC direct ≡ AMC Path B (corrected).
  //   Docs: learn/amc_tts/factored_fIIIa/OMEGA_TT_TO_SCALAR.md
  //         learn/amc_tts/REDUCED_UNREDUCED.md
  //
  // NOT using neithan.tex TT→0 / ph formulas here. Neithan A≡B internally,
  // but Neithan ≠ m / AMC (overall factor ~ λ̂ with J-dependent sign). Prefer
  // AMC printed equations (direct 5×6j; Path B Pandya→RME→inv with
  // chi=barChi, i.e. drop AMC-sample overall minus on inv).
  //
  // Production: Pandya(Ω_red) → RME+DGEMM χ̄ → inv Pandya (AMC Path B)
  //   → MakeNotReduced(χ) → f^III_a ladder.
  // ==================================================================
  if (use_TypeIIIa_1b) {
    // ------------------------------------------------------------------
    // Production Path B: AMC Pandya→RME→inv (NOT Neithan TT→0).
    // Neithan ph formulas disagree with m/AMC by ~λ̂; do not use as gold.
    // Target: learn/amc_tts/factored_fIIIa/output/chi_gamma_via_pandya_*.tex
    //   with chi=barChi (drop AMC-sample overall minus on inv).
    //
    // Pandya(Ω_red) → RME+DGEMM χ̄^γ → inv Pandya
    //   → MakeNotReduced(χ) → f^III_a ladder with Γ_unred.
    //
    // 1) Forward Pandya = IMSRG DoPandya / DoTensorPandya (legs adcb):
    //      λ=0:  bar = −Σ(2J'+1) 6j Ω^{J'}(a,d,c,b)   [comm222_phss]
    //      λ≠0:  bar = −Σ hats (−1)^{jb+jd+Jket+J2} 9j
    //                 Ω^{(J1 J2)λ}(a,d,c,b)             [comm222_phst]
    // 2) χ̄ non-Hermitian: RME [Ω̄⊗Ω̄]^(0) via DGEMM; occ partner (−occ).
    // 3) Inverse Pandya = IMSRG AddInversePandyaTransformation:
    //      χ ← −Σ(2Jp+1) 6j χ̄   (same channel packaging as phss)
    //    Then ForceScalarMakeNotReduced: treat inv output as reduced χ,
    //    convert to unreduced before folding with unreduced Γ.
    // 4) f^III_a: ĵ^{-2} Σ Ĵ^2 (Γ_cpab χ_abcq − χ_pcab Γ_abqc), Γ unreduced.
    // ------------------------------------------------------------------
    Z.modelspace->PreCalculateNineJ();
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
    int nch_ord = Z.modelspace->GetNumberTwoBodyChannels();
    std::deque<arma::mat> barCHI(n_nonzero);

    // IMSRG Pandya of Ω (DoPandyaTransformation / DoTensorPandyaTransformation).
    // ME legs: (a,d,c,b). Tensor Ω is reduced; λ=0 branch uses phss (2J'+1).
    auto pandya_eta = [&](int a, int b, int c, int d, int Jbra,
                          int Jket) -> double {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      Orbit &oc = Z.modelspace->GetOrbit(c);
      Orbit &od = Z.modelspace->GetOrbit(d);
      double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
      double Xbar = 0.0;
      if (lambda == 0) {
        // Commutator::DoPandyaTransformation_SingleChannel
        if (Jbra != Jket)
          return 0.0;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        for (int J_std = jmin; J_std <= jmax; ++J_std) {
          double sixj1 = AngMom::SixJ(ja, jb, Jbra, jc, jd, J_std);
          if (std::abs(sixj1) < 1e-8)
            continue;
          Xbar -= (2 * J_std + 1) * sixj1 *
                  Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
        }
        return Xbar;
      }
      // TensorCommutators::DoTensorPandyaTransformation
      if (not AngMom::Triangle(Jbra, Jket, lambda))
        return 0.0;
      int j1min = std::abs(oa.j2 - od.j2) / 2;
      int j1max = (oa.j2 + od.j2) / 2;
      for (int J1 = j1min; J1 <= j1max; ++J1) {
        int j2min =
            std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
        int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
        for (int J2 = j2min; J2 <= j2max; ++J2) {
          double ninej = Z.modelspace->GetNineJ(ja, jd, J1, jb, jc, J2, Jbra,
                                                Jket, lambda);
          if (std::abs(ninej) < 1e-10)
            continue;
          double hat = std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) *
                                 (2.0 * Jbra + 1.0) * (2.0 * Jket + 1.0));
          Xbar -= hat * Z.modelspace->phase((ob.j2 + od.j2) / 2 + Jket + J2) *
                  ninej * Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
        }
      }
      return Xbar;
    };

    // ---- Build non-Hermitian χ̄^γ in Pandya via RME+DGEMM ----
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      if (nKets_cc < 1)
        continue;
      int J_cc = tbc_cc.J;
      int n2 = nKets_cc * 2;
      int parity_cc = tbc_cc.parity;
      int Tz_cc = tbc_cc.Tz;
      barCHI[ch_cc] = arma::mat(n2, n2, arma::fill::zeros);

      if (lambda == 0) {
        arma::mat bar_Eta(n2, n2, arma::fill::zeros);
        arma::mat nnnbar(n2, n2, arma::fill::zeros);
        for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc) {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          int a = bra_cc.p, b = bra_cc.q;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_a = oa.occ, nbar_a = 1.0 - n_a;
          double n_b = ob.occ, nbar_b = 1.0 - n_b;
          for (int iket_cc = ibra_cc; iket_cc < n2; ++iket_cc) {
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
            double n_c = oc.occ, nbar_c = 1.0 - n_c;
            double n_d = od.occ, nbar_d = 1.0 - n_d;
            double occ_factor =
                nbar_c * nbar_b * n_a * n_d - n_c * n_b * nbar_a * nbar_d;
            if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Eta.GetTRank())
              continue;
            double Xbar = pandya_eta(a, b, c, d, J_cc, J_cc);
            double flip_phase =
                Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);
            if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d)) {
              bar_Eta(ibra_cc, iket_cc) = Xbar;
              nnnbar(ibra_cc, iket_cc) = Xbar * occ_factor;
              if (iket_cc != ibra_cc) {
                bar_Eta(iket_cc, ibra_cc) = hEta * Xbar;
                // NON-HERMITIAN χ packaging: minus on occ for partner
                nnnbar(iket_cc, ibra_cc) = hEta * Xbar * (-occ_factor);
              }
            }
            if (a != b) {
              bar_Eta(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % n2) =
                  Xbar * flip_phase * hEta;
              nnnbar(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % n2) =
                  Xbar * flip_phase * hEta * (-occ_factor);
            }
            if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d)) {
              bar_Eta((iket_cc + nKets_cc) % n2, ibra_cc + nKets_cc) =
                  Xbar * flip_phase;
              nnnbar((iket_cc + nKets_cc) % n2, ibra_cc + nKets_cc) =
                  Xbar * flip_phase * occ_factor;
            }
          }
        }
        barCHI[ch_cc] = bar_Eta * nnnbar; // RME+DGEMM [Ω̄⊗Ω̄]^(0)
      } else {
        // Mid-J RME: χ̄^J = Σ_{J'} λ̂^{-1}(-1)^{J'+λ} Ω̄^{J J'} (occ Ω̄^{J' J})
        for (int Jmid = 0; Jmid <= max_J; ++Jmid) {
          if (not AngMom::Triangle(J_cc, Jmid, lambda))
            continue;
          int ch_mid =
              Z.modelspace->GetTwoBodyChannelIndex(Jmid, parity_cc, Tz_cc);
          if (ch_mid < 0 or ch_mid >= n_nonzero)
            continue;
          TwoBodyChannel_CC &tbc_mid =
              Z.modelspace->GetTwoBodyChannel_CC(ch_mid);
          int nKm = tbc_mid.GetNumberKets();
          if (nKm < 1)
            continue;
          int n2m = nKm * 2;
          arma::mat Om_L(n2, n2m, arma::fill::zeros);
          arma::mat Om_R(n2m, n2, arma::fill::zeros);

          for (int ibra = 0; ibra < n2; ++ibra) {
            int a, b;
            if (ibra < nKets_cc) {
              Ket &bra = tbc_cc.GetKet(ibra);
              a = bra.p;
              b = bra.q;
            } else {
              Ket &bra = tbc_cc.GetKet(ibra - nKets_cc);
              b = bra.p;
              a = bra.q;
            }
            if (ibra >= nKets_cc and a == b)
              continue;
            Orbit &oa = Z.modelspace->GetOrbit(a);
            Orbit &ob = Z.modelspace->GetOrbit(b);
            for (int iket = 0; iket < n2m; ++iket) {
              int c, d;
              if (iket < nKm) {
                Ket &ket = tbc_mid.GetKet(iket);
                c = ket.p;
                d = ket.q;
              } else {
                Ket &ket = tbc_mid.GetKet(iket - nKm);
                d = ket.p;
                c = ket.q;
              }
              if (iket >= nKm and c == d)
                continue;
              Orbit &oc = Z.modelspace->GetOrbit(c);
              Orbit &od = Z.modelspace->GetOrbit(d);
              if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) !=
                  Eta.GetTRank())
                continue;
              Om_L(ibra, iket) = pandya_eta(a, b, c, d, J_cc, Jmid);
            }
          }

          for (int ibra = 0; ibra < n2m; ++ibra) {
            int a, b;
            if (ibra < nKm) {
              Ket &bra = tbc_mid.GetKet(ibra);
              a = bra.p;
              b = bra.q;
            } else {
              Ket &bra = tbc_mid.GetKet(ibra - nKm);
              b = bra.p;
              a = bra.q;
            }
            if (ibra >= nKm and a == b)
              continue;
            Orbit &oa = Z.modelspace->GetOrbit(a);
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double n_a = oa.occ, nbar_a = 1.0 - n_a;
            double n_b = ob.occ, nbar_b = 1.0 - n_b;
            for (int iket = 0; iket < n2; ++iket) {
              int c, d;
              if (iket < nKets_cc) {
                Ket &ket = tbc_cc.GetKet(iket);
                c = ket.p;
                d = ket.q;
              } else {
                Ket &ket = tbc_cc.GetKet(iket - nKets_cc);
                d = ket.p;
                c = ket.q;
              }
              if (iket >= nKets_cc and c == d)
                continue;
              Orbit &oc = Z.modelspace->GetOrbit(c);
              Orbit &od = Z.modelspace->GetOrbit(d);
              double n_c = oc.occ, nbar_c = 1.0 - n_c;
              double n_d = od.occ, nbar_d = 1.0 - n_d;
              double occ_factor =
                  nbar_c * nbar_b * n_a * n_d - n_c * n_b * nbar_a * nbar_d;
              if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) !=
                  Eta.GetTRank())
                continue;
              Om_R(ibra, iket) =
                  pandya_eta(a, b, c, d, Jmid, J_cc) * occ_factor;
            }
          }
          double pref =
              hat_lambda_inv * Z.modelspace->phase(Jmid + lambda);
          barCHI[ch_cc] += pref * Om_L * Om_R;
        }
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fIIIa_pandya"] +=
          omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // ---- Inverse Pandya → scalar χ (still reduced) ----
    Operator ChiOp(*(Z.modelspace), 0, 0, 0, 2);
    ChiOp.SetNonHermitian();
    ChiOp.Erase();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch_ord; ++ch) {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra) {
        Ket &bra = tbc.GetKet(ibra);
        size_t i = bra.p, j = bra.q;
        Orbit &oi = *(bra.op);
        Orbit &oj = *(bra.oq);
        int ji = oi.j2, jj = oj.j2;
        for (int iket = 0; iket < nKets * 2; ++iket) {
          size_t k, l;
          if (iket < nKets) {
            Ket &ket = tbc.GetKet(iket);
            k = ket.p;
            l = ket.q;
          } else {
            Ket &ket = tbc.GetKet(iket - nKets);
            l = ket.p;
            k = ket.q;
          }
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk = ok.j2, jl = ol.j2;
          double commij = 0.0, commji = 0.0;

          int parity_cc = (oi.l + ol.l) % 2;
          int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
          int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
          int Jpmax = std::min(ji + jl, jj + jk) / 2;
          for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
            double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0,
                                                jk * 0.5, jl * 0.5, Jprime);
            if (std::abs(sixj) < 1e-8)
              continue;
            int ch_cc =
                Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
            if (ch_cc < 0 or ch_cc >= n_nonzero)
              continue;
            if (barCHI[ch_cc].n_rows == 0)
              continue;
            TwoBodyChannel_CC &tbc_cc =
                Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int nkets_cc = tbc_cc.GetNumberKets();
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
            if (indx_il < 0 or indx_kj < 0)
              continue;
            indx_il += (i > l ? nkets_cc : 0);
            indx_kj += (k > j ? nkets_cc : 0);
            if (indx_il >= (int)barCHI[ch_cc].n_rows or
                indx_kj >= (int)barCHI[ch_cc].n_cols)
              continue;
            commij -=
                (2 * Jprime + 1) * sixj * barCHI[ch_cc](indx_il, indx_kj);
          }

          parity_cc = (oi.l + ok.l) % 2;
          Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
          Jpmin = std::max(std::abs(jj - jl), std::abs(jk - ji)) / 2;
          Jpmax = std::min(jj + jl, jk + ji) / 2;
          for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
            double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0,
                                                jk * 0.5, jl * 0.5, Jprime);
            if (std::abs(sixj) < 1e-8)
              continue;
            int ch_cc =
                Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
            if (ch_cc < 0 or ch_cc >= n_nonzero)
              continue;
            if (barCHI[ch_cc].n_rows == 0)
              continue;
            TwoBodyChannel_CC &tbc_cc =
                Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int nkets_cc = tbc_cc.GetNumberKets();
            int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
            int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
            if (indx_ik < 0 or indx_lj < 0)
              continue;
            indx_ik += (k > i ? nkets_cc : 0);
            indx_lj += (j > l ? nkets_cc : 0);
            if (indx_lj >= (int)barCHI[ch_cc].n_rows or
                indx_ik >= (int)barCHI[ch_cc].n_cols)
              continue;
            commji -=
                (2 * Jprime + 1) * sixj * barCHI[ch_cc](indx_lj, indx_ik);
          }

          double zijkl =
              (commij - Z.modelspace->phase((ji + jj) / 2 - J0) * commji);
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;
          if (iket < nKets)
            ChiOp.TwoBody.GetMatrix(ch, ch)(ibra, iket) += zijkl;
          if (iket >= nKets)
            ChiOp.TwoBody.GetMatrix(ch, ch)(ibra, iket % nKets) -=
                zijkl * Z.modelspace->phase((jk + jl) / 2 - J0);
        }
      }
    }

    // Inv Pandya wrote χ in *reduced* packaging (Ω was reduced tensor).
    // Γ is scalar unreduced → must MakeNotReduced(χ) before the ladder.
    ForceScalarMakeNotReduced(ChiOp);

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chigamma"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // ---- f^III_a ladder: unreduced χ × unreduced Γ ----
    // Dense table avoids antisym GetTBME index-order pitfalls for NH χ.
    const int n_orb = norbits;
    const int n_J = max_J + 1;
    const size_t chi_stride_l = (size_t)n_J;
    const size_t chi_stride_k = (size_t)n_orb * chi_stride_l;
    const size_t chi_stride_j = (size_t)n_orb * chi_stride_k;
    const size_t chi_stride_i = (size_t)n_orb * chi_stride_j;
    std::vector<double> chi_tab((size_t)n_orb * chi_stride_i, 0.0);
    auto chi_index = [&](index_t i, index_t j, index_t k, index_t l,
                         int J0) -> size_t {
      return (size_t)i * chi_stride_i + (size_t)j * chi_stride_j +
             (size_t)k * chi_stride_k + (size_t)l * chi_stride_l + (size_t)J0;
    };
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < n_orb; ++ii) {
      auto i = allorb_vec[ii];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (auto j : Z.modelspace->all_orbits) {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        for (auto k : Z.modelspace->all_orbits) {
          Orbit &ok = Z.modelspace->GetOrbit(k);
          for (auto l : Z.modelspace->all_orbits) {
            Orbit &ol = Z.modelspace->GetOrbit(l);
            for (int J0 = 0; J0 <= max_J; ++J0) {
              if (not AngMom::Triangle(oi.j2 * 0.5, oj.j2 * 0.5, (double)J0))
                continue;
              if (not AngMom::Triangle(ok.j2 * 0.5, ol.j2 * 0.5, (double)J0))
                continue;
              chi_tab[chi_index(i, j, k, l, J0)] =
                  ChiOp.TwoBody.GetTBME_J(J0, J0, i, j, k, l);
            }
          }
        }
      }
    }

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
        for (auto &c : Z.modelspace->all_orbits)
          for (auto &a : Z.modelspace->all_orbits)
            for (auto &b : Z.modelspace->all_orbits)
              for (int J0 = 0; J0 <= max_J; ++J0) {
                double g1 = Gamma.TwoBody.GetTBME_J(J0, J0, c, p, a, b);
                double g2 = Gamma.TwoBody.GetTBME_J(J0, J0, a, b, q, c);
                if (std::abs(g1) < 1e-14 && std::abs(g2) < 1e-14)
                  continue;
                zpq += (2.0 * J0 + 1.0) *
                       (g1 * chi_tab[chi_index(a, b, c, q, J0)] -
                        chi_tab[chi_index(p, c, a, b, J0)] * g2);
              }
        Z.OneBody(p, q) += zpq / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * zpq / (op.j2 + 1.0);
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_fIIIa"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
  }


  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z) {

  // Scalar Ω/Γ/Z: un-reduce for consistency with the scalar Factorized path.
  // Tensor Ω (rank_J>0) stays reduced — MakeNotReduced refuses rank_J>0
  // (same convention as comm223_231_st).
  const Operator *Etap = &Eta;
  const Operator *Gammap = &Gamma;
  Operator Etatmp, Gammatmp;

  if (Eta.IsReduced() and Eta.GetJRank() == 0) {
    Etatmp = Eta;
    Etatmp.MakeNotReduced();
    Etap = &Etatmp;
  }
  if (Gamma.IsReduced() and Gamma.GetJRank() == 0) {
    Gammatmp = Gamma;
    Gammatmp.MakeNotReduced();
    Gammap = &Gammatmp;
  }

  bool z_was_reduced = Z.IsReduced();
  if (z_was_reduced and Z.GetJRank() == 0)
    Z.MakeNotReduced();

  if (use_1b_intermediates) {
    comm223_232_chi1b_tensor(*Etap, *Gammap,
                             Z); // topology with 1-body intermediate (fast)
  }
  if (use_2b_intermediates) {
    comm223_232_chi2b(*Etap, *Gammap,
                      Z); // topology with 2-body intermediate (slow)
  }

  if (z_was_reduced and Z.GetJRank() == 0)
    Z.MakeReduced();

  return;
}

////////////////////////////////////////////////////////////////////////////
/// factorized 223_232 double commutator with 1b intermediate
///
/// Gamma^I  : chi^epsilon (Omega x Omega -> scalar 1b) x Gamma  -- any lambda
/// Gamma^II : chi^zeta (Gamma x Omega -> tensor 1b) x Omega -> scalar
///            Universal AMC Path B RME (any lambda). Gold: m == AMC direct == Path B.
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
  bool tensor_case =
      (Eta.GetJRank() != 0) && (Gamma.GetJRank() == 0) && Z_is_scalar;
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

          if (tensor_case) {
            if (tbc_ket.J != J)
              continue;
            int phase_pg = bra.Phase(J);
            int phase_qh = ket.Phase(J);
            for (auto d : Z.modelspace->all_orbits) {
              Orbit &od = Z.modelspace->GetOrbit(d);
              if (od.j2 == op.j2)
                zpqrs += CHI_I(p, d) *
                         Gamma.TwoBody.GetTBME_J(J, J, d, q, r, s);
              if (od.j2 == oq.j2)
                zpqrs += phase_pg * CHI_I(q, d) *
                         Gamma.TwoBody.GetTBME_J(J, J, d, p, r, s);
              if (od.j2 == os.j2)
                zpqrs += CHI_I(d, s) *
                         Gamma.TwoBody.GetTBME_J(J, J, p, q, r, d);
              if (od.j2 == oR.j2)
                zpqrs += phase_qh * CHI_I(d, r) *
                         Gamma.TwoBody.GetTBME_J(J, J, p, q, d, s);
            }
            if (p == q)
              zpqrs /= PhysConst::SQRT2;
            if (r == s)
              zpqrs /= PhysConst::SQRT2;
          } else {
            for (auto b : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
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
            for (auto b : Eta.OneBodyChannels.at({oq.l, oq.j2, oq.tz2})) {
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
            for (auto b : Eta.OneBodyChannels.at({oR.l, oR.j2, oR.tz2})) {
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
            for (auto b : Eta.OneBodyChannels.at({os.l, os.j2, os.tz2})) {
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
  // Gamma^II -- chi^zeta (reduced tensor 1b, any lambda) x Omega
  //
  // M-scheme (analyze):
  //   Γ^II = Σ_a [ (1-P_ij) χ_aj Ω_iakl  -  (1-P_kl) χ_ak Ω_ijal ]
  // Bra W and ket V are not Hermitian alone; W - V is.
  // Do NOT use W + h W^T (χ^ζ is non-Hermitian; ≠ W-V).
  //
  // AMC: strip P → G2_Wbra_noperm / G2_Wket_noperm → restore (1-P) by hand.
  // Docs: learn/amc_tts/factored_GII/NOTES.md
  // ######################################################################
  if (do_GII) {
    int max_j2 = 0;
    for (auto x : Z.modelspace->all_orbits)
      max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
    const int max_J = max_j2;
    std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                                Z.modelspace->all_orbits.end());
    const int norb = (int)allorb.size();

    // chi^zeta_ij^lambda = 1/2 * (-1)^{j_j+lambda} * sum w * (-1)^{J0+j_a}
    //   * hatJ0 * hatJ1 * SixJ(lambda,J1,J0; ja,ji,jj) * Gamma_aibc * Omega_bcaj
    // w = n_a n_b nbar_c + nbar_a nbar_b n_c.  Reduced tensor 1b.
    arma::mat Chi_zeta(Gamma.OneBody.n_rows, Gamma.OneBody.n_cols,
                       arma::fill::zeros);

#pragma omp parallel
    {
      arma::mat Chi_loc(Chi_zeta.n_rows, Chi_zeta.n_cols, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ii = 0; ii < norb; ++ii) {
        const index_t i = allorb[ii];
        Orbit &oi = Z.modelspace->GetOrbit(i);
        const double ji = oi.j2 * 0.5;
        for (auto j : allorb) {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          const double jj = oj.j2 * 0.5;
          if (not AngMom::Triangle(ji, jj, (double)lambda))
            continue;
          double sm = 0.0;
          for (auto a : allorb) {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            const double ja = oa.j2 * 0.5;
            const double n_a = oa.occ, nbar_a = 1.0 - n_a;
            for (auto b : allorb) {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              const double n_b = ob.occ, nbar_b = 1.0 - n_b;
              for (auto c : allorb) {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                const double n_c = oc.occ, nbar_c = 1.0 - n_c;
                const double w =
                    n_a * n_b * nbar_c + nbar_a * nbar_b * n_c;
                if (std::abs(w) < 1e-12)
                  continue;
                for (int J0 = 0; J0 <= max_J; ++J0) {
                  if (not AngMom::Triangle(ja, ji, (double)J0))
                    continue;
                  if (not AngMom::Triangle(ob.j2 * 0.5, oc.j2 * 0.5,
                                           (double)J0))
                    continue;
                  const double g =
                      Gamma.TwoBody.GetTBME_J(J0, J0, a, i, b, c);
                  if (std::abs(g) < 1e-16)
                    continue;
                  for (int J1 = 0; J1 <= max_J; ++J1) {
                    if (not AngMom::Triangle(J0, J1, lambda))
                      continue;
                    if (not AngMom::Triangle(ja, jj, (double)J1))
                      continue;
                    const double om =
                        Eta.TwoBody.GetTBME_J(J0, J1, b, c, a, j);
                    if (std::abs(om) < 1e-16)
                      continue;
                    const double sixj =
                        AngMom::SixJ(lambda, J1, J0, ja, ji, jj);
                    if (std::abs(sixj) < 1e-16)
                      continue;
                    const double ph = Z.modelspace->phase(
                        (oj.j2 + oa.j2) / 2 + lambda + J0);
                    sm += ph * hat(J0) * hat(J1) * sixj * w * g * om;
                  }
                }
              }
            }
          }
          Chi_loc(i, j) += 0.5 * sm;
        }
      }
#pragma omp critical
      { Chi_zeta += Chi_loc; }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_232_eths_chizeta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // Fold: Γ^II = W - V with (1-P) restored on each topology.
    // W: AMC G2_Wbra_noperm + bra exchange; V: G2_Wket_noperm + ket exchange.
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

          double W = 0.0; // bra (1-P_ij) χ_ja Ω_iakl
          double V = 0.0; // ket (1-P_kl) χ_ak Ω_ijal
          for (auto a : allorb) {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            const double ja = oa.j2 * 0.5;
            for (int J2 = 0; J2 <= max_J; ++J2) {
              if (not AngMom::Triangle(J, J2, lambda))
                continue;
              const double pref = hatJ_inv * hat(J2) * hat_lambda_inv;

              // W bare (AMC Wbra) + bra exchange
              {
                const double c = Chi_zeta(j, a);
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
                const double c = Chi_zeta(i, a);
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
              // V bare (AMC Wket) + ket exchange
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

          double z = W - V; // Γ^II
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
  // global variables
  double t_start = omp_get_wtime();
  double t_internal = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  int norbits = Z.modelspace->all_orbits.size();
  // Two Body channels
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  size_t nch = ch_bra_list.size();
  int nch_eta = Eta.modelspace->GetNumberTwoBodyChannels();
  // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB
  // channels
  int n_nonzero =
      Eta.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
  auto &Z2 = Z.TwoBody;

  bool Z_is_scalar = true;
  if (Z.TwoBody.rank_T != 0) {
    Z_is_scalar = false;
  }
  // determine symmetry
  int hEta = Eta.IsHermitian() ? 1 : -1;
  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  // int hZ = Z.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  // χ^η / χ^θ rule (TTS): tensor×tensor → scalar (λ_χ=0). Downstream
  // Pandya reverse / DGEMM / RC is always the scalar Factorized path.
  // Ω may still have λ≠0; only the coupled product is rank 0.
  bool tensor_eta_case = (Eta.GetJRank() != 0) && (Gamma.GetJRank() == 0) && Z_is_scalar;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv =
      tensor_eta_case ? 1.0 / std::sqrt(2.0 * lambda + 1.0) : 0.0;

  // Tensor Pandya recoupling for barred Omega with adcb convention:
  // barred indices (a,j,k,b) map to direct TBME legs (a,b;j,k).
  auto barred_eta_tensor = [&](index_t a, index_t j, index_t k, index_t b,
                               int J0, int J1) -> double {
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ob = Z.modelspace->GetOrbit(b);

    if (not AngMom::Triangle(oa.j2 / 2, ob.j2 / 2, J0))
      return 0.0;
    if (not AngMom::Triangle(oj.j2 / 2, ok.j2 / 2, J1))
      return 0.0;
    if (not AngMom::Triangle(J1, J0, lambda))
      return 0.0;

    int pref_phase = Z.modelspace->phase(J0 + (oa.j2 + ok.j2) / 2 + lambda);
    double pref_hat = std::sqrt((2.0 * J0 + 1.0) * (2.0 * J1 + 1.0));

    int J2min = std::abs(oa.j2 - ob.j2) / 2;
    int J2max = (oa.j2 + ob.j2) / 2;
    int J3min = std::abs(oj.j2 - ok.j2) / 2;
    int J3max = (oj.j2 + ok.j2) / 2;

    double sum = 0.0;
    for (int J2 = J2min; J2 <= J2max; ++J2) {
      for (int J3 = J3min; J3 <= J3max; ++J3) {
        if (not AngMom::Triangle(J2, J3, lambda))
          continue;

        int j0min = std::max({std::abs(oa.j2 - 2 * J3),
                              std::abs(ok.j2 - 2 * J0),
                              std::abs(oj.j2 - 2 * lambda)}) /
                    2;
        int j0max = std::min({oa.j2 + 2 * J3,
                              ok.j2 + 2 * J0,
                              oj.j2 + 2 * lambda}) /
                    2;

        for (int j0 = j0min; j0 <= j0max; ++j0) {
          double sixj1 = Z.modelspace->GetSixJ(lambda, J3, J2, oa.j2 / 2.0,
                                               oj.j2 / 2.0, j0);
          double sixj2 = Z.modelspace->GetSixJ(oa.j2 / 2.0, ob.j2 / 2.0, J0,
                                               ok.j2 / 2.0, j0, J3);
          double sixj3 = Z.modelspace->GetSixJ(J1, J0, lambda, j0,
                                               oj.j2 / 2.0, ok.j2 / 2.0);
          if (std::abs(sixj1) < 1e-12 || std::abs(sixj2) < 1e-12 ||
              std::abs(sixj3) < 1e-12)
            continue;

          double omega = Eta.TwoBody.GetTBME_J(J2, J3, a, b, j, k);
          if (std::abs(omega) < 1e-12)
            continue;

          double hats = std::sqrt((2.0 * J2 + 1.0) * (2.0 * J3 + 1.0)) *
                        (2.0 * j0 + 1.0);
          sum += Z.modelspace->phase(J2) * hats * sixj1 * sixj2 * sixj3 * omega;
        }
      }
    }

    return pref_phase * pref_hat * sum;
  };

  // ***********************************************************************************
  // //
  //                             Diagram II and III //
  // ***********************************************************************************
  // //

  //______________________________________________________________________
  // global array
  std::deque<arma::mat> bar_Eta(n_nonzero);   // released
  std::deque<arma::mat> bar_Gamma(n_nonzero); // released
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    // because the restriction a<b in the bar and ket vector, if we want to
    // store the full Pandya transformed matrix, we twice the size of matrix
    bar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    bar_Gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
  }

  std::deque<arma::mat> barCHI_III(n_nonzero);    //  released
  std::deque<arma::mat> bar_CHI_V(n_nonzero);     // released
  std::deque<arma::mat> bar_CHI_VI(n_nonzero);    //  released
  std::deque<arma::mat> bar_CHI_VI_II(n_nonzero); //  released

  /// Pandya transformation
  /// construct bar_Gamma, bar_Eta, nnnbar_Eta, nnnbar_Eta_d, bar_CHI_VI
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    if (nKets_cc < 1) {
      continue;
    }

    arma::mat nnnbar_Eta =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    arma::mat nnnbar_Eta_d =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc) {
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
      double ja = oa.j2 * 0.5;
      double n_a = oa.occ;
      double nbar_a = 1.0 - n_a;

      Orbit &ob = Z.modelspace->GetOrbit(b);
      double jb = ob.j2 * 0.5;
      double n_b = ob.occ;
      double nbar_b = 1.0 - n_b;

      // loop over cross-coupled kets |cd> in this channel
      for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
        if ((iket_cc % nKets_cc) < ibra_cc)
          continue; // We'll get these from symmetry
        int c, d;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
          c = ket_cc_cd.p;
          d = ket_cc_cd.q;
        } else {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc_cd.p;
          c = ket_cc_cd.q;
        }
        // if (iket_cc >= nKets_cc and c == d)
        //  continue;

        Orbit &oc = Z.modelspace->GetOrbit(c);
        double jc = oc.j2 * 0.5;
        double n_c = oc.occ;
        double nbar_c = 1.0 - n_c;

        Orbit &od = Z.modelspace->GetOrbit(d);
        double jd = od.j2 * 0.5;
        double n_d = od.occ;
        double nbar_d = 1.0 - n_d;

        double occ_AbarBC = (nbar_a * n_b * n_c + n_a * nbar_b * nbar_c);
        double occ_ABbarD = (n_a * nbar_b * n_d + nbar_a * n_b * nbar_d);

        double occ_BCDbar = (n_b * n_c * nbar_d + nbar_b * nbar_c * n_d);
        double occ_ACbarD = (n_a * nbar_c * n_d + nbar_a * n_c * nbar_d);

        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        double Etabar = 0;
        double Gammabar = 0;
        int dJ_std = 1;
        if ((a == d or b == c)) {
          dJ_std = 2;
          jmin += jmin % 2;
        }
        for (int J_std = jmin; J_std <= jmax; J_std += dJ_std) {
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8) {
            Etabar -= (2 * J_std + 1) * sixj1 *
                      Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
            Gammabar -= (2 * J_std + 1) * sixj1 *
                        Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
          }
        }

        double flip_phase =
            Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);
        if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d)) {
          // direct term
          bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
          bar_Eta[ch_cc](ibra_cc, iket_cc) = Etabar;
          nnnbar_Eta(ibra_cc, iket_cc) = Etabar * occ_AbarBC;
          nnnbar_Eta_d(ibra_cc, iket_cc) = Etabar * occ_ABbarD;

          if (iket_cc != ibra_cc) {
            // Hermiticity: Xbar_cdab = hX * Xbar_abcd.
            bar_Gamma[ch_cc](iket_cc, ibra_cc) = hGamma * Gammabar;
            bar_Eta[ch_cc](iket_cc, ibra_cc) = hEta * Etabar;
            nnnbar_Eta(iket_cc, ibra_cc) = hEta * Etabar * occ_ACbarD;
            nnnbar_Eta_d(iket_cc, ibra_cc) = hEta * Etabar * occ_BCDbar;
          }
        }

        if (a != b) {
          // By exchange symmetry Xbar_badc = phase * hX * Xbar_abcd.
          bar_Gamma[ch_cc](ibra_cc + nKets_cc,
                           (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Gammabar * flip_phase * hGamma;
          bar_Eta[ch_cc](ibra_cc + nKets_cc,
                         (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Etabar * flip_phase * hEta;
          nnnbar_Eta(ibra_cc + nKets_cc,
                     (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Etabar * flip_phase * hEta * occ_ABbarD;
          nnnbar_Eta_d(ibra_cc + nKets_cc,
                       (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Etabar * flip_phase * hEta * occ_AbarBC;
        }

        if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d)) {
          // Combined exchange symmetry and hermiticity
          // Xbar_dcba = phase * Xbar_abcd
          bar_Gamma[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc),
                           ibra_cc + nKets_cc) = Gammabar * flip_phase;
          bar_Eta[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc),
                         ibra_cc + nKets_cc) = Etabar * flip_phase;
          nnnbar_Eta((iket_cc + nKets_cc) % (2 * nKets_cc),
                     ibra_cc + nKets_cc) = Etabar * flip_phase * occ_BCDbar;
          nnnbar_Eta_d((iket_cc + nKets_cc) % (2 * nKets_cc),
                       ibra_cc + nKets_cc) = Etabar * flip_phase * occ_ACbarD;
        }
      }
      //-------------------
    }
    // chi^\eta
    if (not tensor_eta_case) {
      barCHI_III[ch_cc] = bar_Eta[ch_cc] * nnnbar_Eta;
    } else {
      barCHI_III[ch_cc].zeros();
      for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc) {
        index_t i, l;
        if (ibra_cc < nKets_cc) {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          i = bra_cc.p;
          l = bra_cc.q;
        } else {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
          l = bra_cc.p;
          i = bra_cc.q;
        }

        for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
          index_t k, j;
          if (iket_cc < nKets_cc) {
            Ket &ket_cc = tbc_cc.GetKet(iket_cc);
            k = ket_cc.p;
            j = ket_cc.q;
          } else {
            Ket &ket_cc = tbc_cc.GetKet(iket_cc - nKets_cc);
            j = ket_cc.p;
            k = ket_cc.q;
          }

          Orbit &ok = Z.modelspace->GetOrbit(k);
          double n_k = ok.occ;
          double nbar_k = 1.0 - n_k;

          double sum_eta = 0.0;
          for (auto &a : Z.modelspace->all_orbits) {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits) {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              double occ = nbar_a * n_b * nbar_k + n_a * nbar_b * n_k;
              if (std::abs(occ) < 1e-12)
                continue;

              int J2min = std::max(std::abs(ob.j2 - oa.j2),
                                   std::abs(oa.j2 - ob.j2)) /
                          2;
              int J2max = std::min(ob.j2 + oa.j2, oa.j2 + ob.j2) / 2;
              for (int J2 = J2min; J2 <= J2max; ++J2) {
                if (not AngMom::Triangle(J_cc, J2, lambda))
                  continue;

                double bo1 = barred_eta_tensor(i, b, a, j, J_cc, J2);
                double bo2 = barred_eta_tensor(a, l, k, b, J2, J_cc);
                if (std::abs(bo1) < 1e-12 || std::abs(bo2) < 1e-12)
                  continue;

                sum_eta += occ * Z.modelspace->phase(J2 + lambda) *
                           hat_lambda_inv * bo1 * bo2;
              }
            }
          }

          barCHI_III[ch_cc](ibra_cc, iket_cc) =
              Z.modelspace->phase(J_cc) * sum_eta / (2.0 * J_cc + 1.0);
        }
      }
    }
    bar_CHI_V[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta;
    bar_CHI_VI[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta_d;
    // The following code is stable;
     arma::mat nnnbar_Eta_d_t = nnnbar_Eta_d.t();
     bar_CHI_VI_II[ch_cc] = hEta * (nnnbar_Eta_d_t * bar_Gamma[ch_cc]);
     // The following original code have unexpected behavior, the result is not stable
    //bar_CHI_VI_II[ch_cc] = hEta * (nnnbar_Eta_d).t() * bar_Gamma[ch_cc];
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  //-------------------------------------------------------------------------------
  // intermediate operator for diagram IIa  and IIc
  // Theintermediate two body operator
  //  Chi_III :
  //            eta |
  //           _____|
  //          /\    |
  //   |     (  )
  //   |_____ \/
  //   | eta
  TwoBodyME Chi_III_Op = Z.TwoBody; // scalar: [Ω×Ω]^0 → χ^η (never copy Eta rank)
  Chi_III_Op.Erase();
  // Inverse Pandya (scalar): χ is always λ=0 after tensor×tensor → scalar.
  //  X^J_ijkl  = - ( 1- P_ij )  sum_J' (2J'+1)  { i j J }  \bar{X}^J'_il`kj`
  //                                             { k l J'}
#pragma omp parallel for
  for (int ch = 0; ch < nch_eta; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J0 = tbc.J;
    int nKets = tbc.GetNumberKets();
    for (int ibra = 0; ibra < nKets; ++ibra) {
      Ket &bra = tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;

      for (int iket = 0; iket < nKets * 2; ++iket) {
        size_t k, l;
        if (iket < nKets) {
          Ket &ket = tbc.GetKet(iket);
          k = ket.p;
          l = ket.q;
        } else {
          Ket &ket = tbc.GetKet(iket - nKets);
          l = ket.p;
          k = ket.q;
        }

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commij = 0;
        double commji = 0;

        // ijkl
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
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
          if (indx_il < 0 or indx_kj < 0)
            continue;
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);

          double me1 = barCHI_III[ch_cc](indx_il, indx_kj);
          commij -= (2 * Jprime + 1) * sixj * me1;
        }

        // jikl, exchange i and j
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
          int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
          int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));

          if (indx_ik < 0 or indx_lj < 0)
            continue;
          indx_ik += (k > i ? nkets_cc : 0);
          indx_lj += (j > l ? nkets_cc : 0);
          double me1 = barCHI_III[ch_cc](indx_lj, indx_ik);
          commji -= (2 * Jprime + 1) * sixj * me1;
        }

        double zijkl =
            (commij - Z.modelspace->phase((ji + jj) / 2 - J0) * commji);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        if (iket < nKets)
          Chi_III_Op.GetMatrix(ch, ch)(ibra, iket) += zijkl;
        if (iket >= nKets)
          Chi_III_Op.GetMatrix(ch, ch)(ibra, iket % nKets) -=
              zijkl * Z.modelspace->phase((jk + jl) / 2 - J0);
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  //------------------------------------------------------------------------------
  // intermediate operator for diagram IIIc and IIId
  //------------------------------------------------------------------------------
  //  The intermediate two body operator
  //  Chi_VI :
  //            eta |
  //           _____|
  //          /\    |
  //         (  )     |
  //          \/~~~~~~|
  //            gamma |
  //
  //  Chi_VI_cdqh = \sum_{ab} (nbar_a * n_b * n_c - nbar_a * nbar_b * n_c )
  //                 ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
  //
  //                { J3 J4 J0 } { J3 J4 J0 }
  //                { jc jd ja } { jq jh jb }
  //
  //                \bar{Eta}_bq`ac`  Gamma_dahb
  //-------------------------------------------------------------------------------
  TwoBodyME Chi_VI_Op = Gamma.TwoBody;
  Chi_VI_Op.Erase();

  TwoBodyME Chi_VI_II_Op = Gamma.TwoBody;
  Chi_VI_II_Op.Erase();
  // BUILD CHI_VI
  // Inverse Pandya transformation
#pragma omp parallel for
  for (int ch = 0; ch < nch; ++ch) {
    size_t ch_bra = ch_bra_list[ch];
    size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();
    if (nbras == 0 or nkets == 0)
      continue;
    int J0 = tbc_bra.J;
    for (int ibra = 0; ibra < nbras * 2; ++ibra) {
      size_t i, j;
      if (ibra < nbras) {
        Ket &bra = tbc_bra.GetKet(ibra);
        i = bra.p;
        j = bra.q;
      } else {
        Ket &bra = tbc_bra.GetKet(ibra - nbras);
        i = bra.q;
        j = bra.p;
      }
      // if (ibra >= nbras and i == j)
      //   continue;

      Orbit &oi = Z.modelspace->GetOrbit(i);
      int ji = oi.j2;
      Orbit &oj = Z.modelspace->GetOrbit(j);
      int jj = oj.j2;

      for (int iket = 0; iket < nkets * 2; ++iket) {
        size_t k, l;
        if (iket < nkets) {
          Ket &ket = tbc_ket.GetKet(iket);
          k = ket.p;
          l = ket.q;
        } else {
          Ket &ket = tbc_ket.GetKet(iket - nkets);
          k = ket.q;
          l = ket.p;
        }
        // if (iket >= nkets and k == l)
        //   continue;

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commijlk = 0;

        double commijkld = 0;
        double commjikld = 0;

        // ijkl, direct term        -->  il kj
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
          int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_kj < 0)
            continue;

          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          double me1 = bar_CHI_VI[ch_cc](indx_il, indx_kj);
          commijkl -= (2 * Jprime + 1) * sixj * me1;
          double me12 = bar_CHI_VI_II[ch_cc](indx_il, indx_kj);
          commijkld -= (2 * Jprime + 1) * sixj * me12;
        }

        // ijlk,  exchange k and l -->  ik lj
        parity_cc = (oi.l + ok.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jk), std::abs(jj - jl)) / 2;
        Jpmax = std::min(ji + jk, jj + jl) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5,
                                              jk * 0.5, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
          int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
          if (indx_ik < 0 or indx_lj < 0)
            continue;

          int indx_ki = indx_ik;
          int indx_jl = indx_lj;

          // exchange k and l
          indx_ik += (i > k ? nkets_cc : 0);
          indx_lj += (l > j ? nkets_cc : 0);
          double me2 = bar_CHI_VI[ch_cc](indx_ik, indx_lj);
          commijlk -= (2 * Jprime + 1) * sixj * me2;

          indx_ki += (k > i ? nkets_cc : 0);
          indx_jl += (j > l ? nkets_cc : 0);
          double me21 = bar_CHI_VI_II[ch_cc](indx_jl, indx_ki);
          commjikld -= (2 * Jprime + 1) * sixj * me21;
        }

        double zijkl =
            (commijkl - Z.modelspace->phase((jk + jl) / 2 - J0) * commijlk);
        double zijkl_II =
            (commijkld - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikld);

        if (i == j) {
          zijkl /= PhysConst::SQRT2;
          zijkl_II /= PhysConst::SQRT2;
        }
        if (k == l) {
          zijkl /= PhysConst::SQRT2;
          zijkl_II /= PhysConst::SQRT2;
        }

        if (iket < nkets) {
          if (ibra < nbras) {
            Chi_VI_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket) += zijkl;
          }
          if (ibra >= nbras) {
            Chi_VI_Op.GetMatrix(ch_bra, ch_ket)(ibra % nbras, iket) -=
                zijkl * Z.modelspace->phase((ji + jj) / 2 - J0);
          }
        }

        if (ibra < nbras) {
          if (iket < nkets) {
            Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket) += zijkl_II;
          }
          if (iket >= nkets) {
            Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket % nkets) -=
                zijkl_II * Z.modelspace->phase((jk + jl) / 2 - J0);
          }
        }
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // release memory
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    bar_CHI_VI[ch_cc].clear();
    bar_CHI_VI_II[ch_cc].clear();
  }
  bar_CHI_VI.clear();
  bar_CHI_VI_II.clear();

// Diagram IIa and Diagram IIc
// Diagram IIIc and Diagram IIId
#pragma omp parallel for
  for (int ch = 0; ch < nch; ++ch) {
    size_t ch_bra = ch_bra_list[ch];
    size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();

    if (nbras < 1 or nkets < 1)
      continue;
    // Diagram IIa and IIc
    arma::mat Multi_matirx = Chi_III_Op.GetMatrix(ch_bra, ch_bra) *
                             Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);
    Multi_matirx += hZ * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket) *
                    (Chi_III_Op.GetMatrix(ch_ket, ch_ket).t());
    // Diagram IIIc and Diagram IIId
    Multi_matirx +=
        -Eta.TwoBody.GetMatrix(ch_bra) * Chi_VI_Op.GetMatrix(ch_bra, ch_ket) -
        (Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket) *
         Eta.TwoBody.GetMatrix(ch_ket));
    Z2.GetMatrix(ch_bra, ch_ket) += Multi_matirx;
  } // J0 channel

  // release memory
  Chi_III_Op.Deallocate();
  Chi_VI_Op.Deallocate();
  Chi_VI_II_Op.Deallocate();

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // Diagram IIb and Diagram IId
  std::deque<arma::mat> barCHI_III_RC(
      n_nonzero); // released Recoupled bar CHI_III
  std::deque<arma::mat> bar_CHI_V_RC(n_nonzero); // released
  /// build intermediate bar operator
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    // because the restriction a<b in the bar and ket vector, if we want to
    // store the full Pandya transformed matrix, we twice the size of matrix
    barCHI_III_RC[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    bar_CHI_V_RC[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
  }

  //------------------------------------------------------------------------------
  // intermediate operator for diagram IIb and IId    IIIa and IIIb
  //------------------------------------------------------------------------------
  /// Pandya transformation only recouple the angula momentum
  /// IIb and IId
  /// diagram IIIa - diagram IIIb
  //  \bar{X}^J_ijkl  = sum_J' (2J'+1)  { i j J }  (-)^(j+k+J')
  //  \bar{X}^J'_il`jk`
  //                                    { k l J'}
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc) {
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
      double ja = oa.j2 * 0.5;
      double jb = ob.j2 * 0.5;

      // loop over cross-coupled kets |cd> in this channel
      for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
        int c, d;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
          c = ket_cc_cd.p;
          d = ket_cc_cd.q;
        } else {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc_cd.p;
          c = ket_cc_cd.q;
        }
        if (iket_cc >= nKets_cc and c == d)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        double jc = oc.j2 * 0.5;
        double jd = od.j2 * 0.5;

        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        double XbarIIbd = 0;
        double XbarIIIab = 0;

        for (int J_std = jmin; J_std <= jmax; J_std++) {
          // int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8) {
            int parity_cc = (oa.l + od.l) % 2;
            int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
            int ch_cc_old =
                Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);

            TwoBodyChannel_CC &tbc_cc_old =
                Z.modelspace->GetTwoBodyChannel_CC(ch_cc_old);
            int nkets = tbc_cc_old.GetNumberKets();
            int indx_ad = tbc_cc_old.GetLocalIndex(std::min(int(a), int(d)),
                                                   std::max(int(a), int(d)));
            int indx_bc = tbc_cc_old.GetLocalIndex(std::min(int(b), int(c)),
                                                   std::max(int(b), int(c)));
            if (indx_ad >= 0 and indx_bc >= 0) {
              int indx_cb = indx_bc;
              if (a > d)
                indx_ad += nkets;
              if (b > c)
                indx_bc += nkets;
              if (c > b)
                indx_cb += nkets;
              XbarIIbd -= Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) *
                          (2 * J_std + 1) * sixj1 *
                          (barCHI_III[ch_cc_old](indx_bc, indx_ad) +
                           barCHI_III[ch_cc_old](indx_ad, indx_bc));
              XbarIIIab += Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) *
                           (2 * J_std + 1) * sixj1 *
                           (bar_CHI_V[ch_cc_old](indx_ad, indx_bc) -
                            hZ * bar_CHI_V[ch_cc_old](indx_bc, indx_ad));
            }
          }
        }
        barCHI_III_RC[ch_cc](ibra_cc, iket_cc) = XbarIIbd;
        bar_CHI_V_RC[ch_cc](ibra_cc, iket_cc) = XbarIIIab;
      }
      //-------------------
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  /// release memory
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    barCHI_III[ch_cc].clear();
    bar_CHI_V[ch_cc].clear();
  }
  barCHI_III.clear();
  bar_CHI_V.clear();

  // ##########################################
  //      diagram IIIa - diagram IIIb
  // ##########################################
  std::deque<arma::mat> CHI_V_final(n_nonzero);
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    CHI_V_final[ch_cc] = bar_Eta[ch_cc] * bar_CHI_V_RC[ch_cc];
  }
  /// release memory
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    bar_CHI_V_RC[ch_cc].clear();
  }
  bar_CHI_V_RC.clear();

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  //  Inverse Pandya transformation
  //  diagram IIIa - diagram IIIb
  //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
  //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
  //                                    { k l J'}
#pragma omp parallel for
  for (int ch = 0; ch < nch; ++ch) {
    int ch_bra = ch_bra_list[ch];
    int ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();
    int J0 = tbc_bra.J;

    if (nbras == 0 or nkets == 0)
      continue;

    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;
      int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

      int ketmin = 0;
      if (ch_bra == ch_ket)
        ketmin = ibra;
      for (int iket = ketmin; iket < nkets; ++iket) {
        size_t k, l;
        Ket &ket = tbc_ket.GetKet(iket);
        k = ket.p;
        l = ket.q;

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commjikl = 0;
        double commijlk = 0;
        double commjilk = 0;

        // jikl, direct term        -->  jl  ki
        // ijlk, exchange ij and kl -->  lj  ik
        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);

          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;

          int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_jl < 0 or indx_ik < 0)
            continue;

          int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
          // direct term
          indx_jl += (j > l ? nkets_cc : 0);
          indx_ik += (i > k ? nkets_cc : 0);
          double me1 = CHI_V_final[ch_cc](indx_jl, indx_ik);
          commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

          int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);
          // exchange ij and kl
          double me2 = CHI_V_final[ch_cc](indx_ik, indx_jl);
          commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        // ijkl,  exchange i and j -->  il  kj
        // jilk,  exchange k and l -->  jk li
        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;

          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_jk < 0)
            continue;

          int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
          // exchange k and l
          indx_il += (i > l ? nkets_cc : 0);
          indx_jk += (j > k ? nkets_cc : 0);
          double me1 = CHI_V_final[ch_cc](indx_jk, indx_il);
          commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

          int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);
          // exchange i and j
          double me2 = CHI_V_final[ch_cc](indx_il, indx_jk);
          commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        double zijkl =
            (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;

        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    CHI_V_final[ch_cc].clear();
  }
  CHI_V_final.clear();

  // ######################################
  // declare CHI_IV
  //------------------------------------------------------------------------------
  // intermediate operator for diagram IIe and IIf
  //------------------------------------------------------------------------------
  // The intermediate two body operator
  //  Chi_IV :
  //        q  | eta |  b
  //           |_____|
  //        a  |_____|  c
  //           | eta |
  //        p  |     |  d
  /////////////////////////////////////////////////////////////////////////////////
  std::deque<arma::mat> CHI_IV(nch_eta); // released
  for (int ch = 0; ch < nch_eta; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int nKets = tbc.GetNumberKets();
    // Not symmetric
    CHI_IV[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
  }

  //------------------------------------------------------------------------------
  //                      Factorization of IIIe and IIIf
  //------------------------------------------------------------------------------
  //
  // The intermediate two body operator
  //  CHI_VII :
  //        g  |     |  c
  //           |~~~~~|
  //        a  |     |  b
  //           |_____|
  //        h  |     |  d
  //------------------------------------------------------------------------------
  std::deque<arma::mat> CHI_VII(nch); // released
  for (int ch = 0; ch < nch; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int nKets = tbc.GetNumberKets();
    // Not symmetric
    CHI_VII[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // this loop appears to be broken.
  // full matrix
#pragma omp parallel for
  for (int ch = 0; ch < nch_eta; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J0 = tbc.J;
    int nKets = tbc.GetNumberKets();
    arma::mat Eta_matrix(2 * nKets, 2 * nKets);
    arma::mat Eta_matrix_c(2 * nKets, 2 * nKets);
    arma::mat Eta_matrix_d(2 * nKets, 2 * nKets);
    arma::mat Gamma_matrix(2 * nKets, 2 * nKets);

    for (int ibra = 0; ibra < nKets; ++ibra) {
      Ket &bra = tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;
      double n_i = oi.occ;
      double bar_n_i = 1. - n_i;
      double n_j = oj.occ;
      double bar_n_j = 1. - n_j;

      for (int iket = 0; iket < nKets; ++iket) {
        size_t k, l;
        Ket &ket = tbc.GetKet(iket);
        k = ket.p;
        l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double n_k = ok.occ;
        double bar_n_k = 1. - n_k;
        double n_l = ol.occ;
        double bar_n_l = 1. - n_l;

        double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
        double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

        double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);
        double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);

        Eta_matrix(ibra, iket) = EtaME;
        Eta_matrix_c(ibra, iket) = occfactor_k * EtaME;
        Eta_matrix_d(ibra, iket) = occfactor_l * EtaME;
        if (Z_is_scalar)
          Gamma_matrix(ibra, iket) = GammaME;
        if (i != j) {
          int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
          Eta_matrix(ibra + nKets, iket) = phase * EtaME;
          Eta_matrix_c(ibra + nKets, iket) = occfactor_k * phase * EtaME;
          Eta_matrix_d(ibra + nKets, iket) = occfactor_l * phase * EtaME;
          if (Z_is_scalar)
            Gamma_matrix(ibra + nKets, iket) = phase * GammaME;
          if (k != l) {
            phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
            Eta_matrix(ibra + nKets, iket + nKets) = phase * EtaME;
            Eta_matrix_c(ibra + nKets, iket + nKets) =
                occfactor_l * phase * EtaME;
            Eta_matrix_d(ibra + nKets, iket + nKets) =
                occfactor_k * phase * EtaME;
            if (Z_is_scalar)
              Gamma_matrix(ibra + nKets, iket + nKets) = phase * GammaME;

            phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
            Eta_matrix(ibra, iket + nKets) = phase * EtaME;
            Eta_matrix_c(ibra, iket + nKets) = occfactor_l * phase * EtaME;
            Eta_matrix_d(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            if (Z_is_scalar)
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
          }
        } else {
          if (k != l) {
            int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
            Eta_matrix(ibra, iket + nKets) = phase * EtaME;
            Eta_matrix_c(ibra, iket + nKets) = occfactor_l * phase * EtaME;
            Eta_matrix_d(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            if (Z_is_scalar)
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
          }
        }
      }
    }
    // TODO: We can use symmetry here so that we don't have to use the full
    // matrix.
    if (not tensor_eta_case) {
      CHI_IV[ch] = Eta_matrix * Eta_matrix_c;
      CHI_IV[ch] += (Eta_matrix * Eta_matrix_d).t();
    } else {
      // Tensor chi^theta from two Eta operators; result remains scalar and
      // can reuse the existing scalar chi^theta * Gamma contraction path.
      CHI_IV[ch].zeros();
      for (int ibra_ext = 0; ibra_ext < 2 * nKets; ++ibra_ext) {
        index_t i, j;
        if (ibra_ext < nKets) {
          Ket &bra_ext = tbc.GetKet(ibra_ext);
          i = bra_ext.p;
          j = bra_ext.q;
        } else {
          Ket &bra_ext = tbc.GetKet(ibra_ext - nKets);
          j = bra_ext.p;
          i = bra_ext.q;
        }

        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double n_j = oj.occ;
        double nbar_j = 1.0 - n_j;

        for (int iket_ext = 0; iket_ext < 2 * nKets; ++iket_ext) {
          index_t k, l;
          if (iket_ext < nKets) {
            Ket &ket_ext = tbc.GetKet(iket_ext);
            k = ket_ext.p;
            l = ket_ext.q;
          } else {
            Ket &ket_ext = tbc.GetKet(iket_ext - nKets);
            l = ket_ext.p;
            k = ket_ext.q;
          }

          Orbit &ok = Z.modelspace->GetOrbit(k);
          double n_k = ok.occ;
          double nbar_k = 1.0 - n_k;

          double chi_theta = 0.0;
          for (auto &a : Z.modelspace->all_orbits) {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits) {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              double occ = n_a * n_b * nbar_k + nbar_a * nbar_b * n_k +
                           n_a * n_b * nbar_j + nbar_a * nbar_b * n_j;
              if (std::abs(occ) < 1e-12)
                continue;

              int J2min = std::max(std::abs(oi.j2 - oj.j2),
                                   std::abs(oa.j2 - ob.j2)) /
                          2;
              int J2max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J2 = J2min; J2 <= J2max; ++J2) {
                if (not AngMom::Triangle(J0, J2, lambda))
                  continue;

                double eta1 = Eta.TwoBody.GetTBME_J(J0, J2, i, j, a, b);
                double eta2 = Eta.TwoBody.GetTBME_J(J2, J0, a, b, k, l);
                if (std::abs(eta1) < 1e-12 || std::abs(eta2) < 1e-12)
                  continue;

                chi_theta += occ * Z.modelspace->phase(J2 + lambda) *
                             hat_lambda_inv * eta1 * eta2;
              }
            }
          }

          CHI_IV[ch](ibra_ext, iket_ext) =
              Z.modelspace->phase(J0) * chi_theta / (2.0 * J0 + 1.0);
        }
      }
    }
    if (Z_is_scalar)
      CHI_VII[ch] =
          Gamma_matrix * Eta_matrix_d + hEta * Eta_matrix_d.t() * Gamma_matrix;
  }

  //// Full Isospin tensor CHI_VII
  //// diagram IIIe and IIIf
  if (not Z_is_scalar) {
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch) {
      int ch_bra = ch_bra_list[ch];
      int ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nKets = tbc_ket.GetNumberKets();
      if (nbras == 0 or nKets == 0)
        continue;
      int J0 = tbc_bra.J;

      arma::mat Eta_matrix_bra(2 * nbras, 2 * nbras);
      arma::mat Eta_matrix_ket(2 * nKets, 2 * nKets);
      arma::mat Gamma_matrix(2 * nbras, 2 * nKets);

      // full Gamma_matrix Eta_matrix_bra
      for (int ibra = 0; ibra < nbras; ++ibra) {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = *(bra.op);
        Orbit &oj = *(bra.oq);
        int ji = oi.j2;
        int jj = oj.j2;
        double n_i = oi.occ;
        double bar_n_i = 1. - n_i;
        double n_j = oj.occ;
        double bar_n_j = 1. - n_j;

        // full gamma
        for (int iket = 0; iket < nKets; ++iket) {
          size_t k, l;
          Ket &ket = tbc_ket.GetKet(iket);
          k = ket.p;
          l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk = ok.j2;
          int jl = ol.j2;
          double n_k = ok.occ;
          double bar_n_k = 1. - n_k;
          double n_l = ol.occ;
          double bar_n_l = 1. - n_l;

          double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);
          Gamma_matrix(ibra, iket) = GammaME;
          if (i != j) {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Gamma_matrix(ibra + nbras, iket) = phase * GammaME;
            if (k != l) {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Gamma_matrix(ibra + nbras, iket + nKets) = phase * GammaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
            }
          } else {
            if (k != l) {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
            }
          }
        }

        // full Eta_matrix_bra
        for (int iket = 0; iket < nbras; ++iket) {
          size_t k, l;
          Ket &ket = tbc_bra.GetKet(iket);
          k = ket.p;
          l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk = ok.j2;
          int jl = ol.j2;
          double n_k = ok.occ;
          double bar_n_k = 1. - n_k;
          double n_l = ol.occ;
          double bar_n_l = 1. - n_l;

          double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
          double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

          double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);

          Eta_matrix_bra(ibra, iket) = occfactor_l * EtaME;
          if (i != j) {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix_bra(ibra + nbras, iket) = occfactor_l * phase * EtaME;
            if (k != l) {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix_bra(ibra + nbras, iket + nbras) =
                  occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_bra(ibra, iket + nbras) = occfactor_k * phase * EtaME;
            }
          } else {
            if (k != l) {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_bra(ibra, iket + nbras) = occfactor_k * phase * EtaME;
            }
          }
        }
      }

      // full Eta_matrix_ket
      for (int ibra = 0; ibra < nKets; ++ibra) {
        Ket &bra = tbc_ket.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = *(bra.op);
        Orbit &oj = *(bra.oq);
        int ji = oi.j2;
        int jj = oj.j2;
        double n_i = oi.occ;
        double bar_n_i = 1. - n_i;
        double n_j = oj.occ;
        double bar_n_j = 1. - n_j;

        for (int iket = 0; iket < nKets; ++iket) {
          size_t k, l;
          Ket &ket = tbc_ket.GetKet(iket);
          k = ket.p;
          l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk = ok.j2;
          int jl = ol.j2;
          double n_k = ok.occ;
          double bar_n_k = 1. - n_k;
          double n_l = ol.occ;
          double bar_n_l = 1. - n_l;
          double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
          double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

          double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);

          Eta_matrix_ket(ibra, iket) = occfactor_l * EtaME;
          if (i != j) {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix_ket(ibra + nKets, iket) = occfactor_l * phase * EtaME;
            if (k != l) {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix_ket(ibra + nKets, iket + nKets) =
                  occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_ket(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          } else {
            if (k != l) {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_ket(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
        }
      }

      CHI_VII[ch] = Gamma_matrix * Eta_matrix_ket +
                    hEta * Eta_matrix_bra.t() * Gamma_matrix;
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  //_______________________________________________________________________________
  std::deque<arma::mat> bar_CHI_IV(n_nonzero);     // released
  std::deque<arma::mat> bar_CHI_VII_CC(n_nonzero); // released
  /// build intermediate bar operator
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    if (nKets_cc < 1)
      continue;

    // because the restriction a<b in the bar and ket vector, if we want to
    // store the full Pandya transformed matrix, we twice the size of matrix
    bar_CHI_IV[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    bar_CHI_VII_CC[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
  }

  /// Pandya transformation only recouple the angula momentum
  /// IIe and IIf                 barCHI_III_RC   bar_CHI_IV
  /// diagram IIIe and IIIf       bar_CHI_VII_CC
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    if (nKets_cc < 1) {
      continue;
    }

    int J_cc = tbc_cc.J;
    for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc) {
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
      double ja = oa.j2 * 0.5;
      double jb = ob.j2 * 0.5;

      // loop over cross-coupled kets |cd> in this channel
      for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
        int c, d;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
          c = ket_cc_cd.p;
          d = ket_cc_cd.q;
        } else {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc_cd.p;
          c = ket_cc_cd.q;
        }
        if (iket_cc >= nKets_cc and c == d)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        double jc = oc.j2 * 0.5;
        double jd = od.j2 * 0.5;

        int Tz_J2_bc = (ob.tz2 + oc.tz2) / 2;
        int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
        int parity_J2 = (ob.l + oc.l) % 2;

        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        double XbarIIef = 0;
        double XbarIIIef = 0;
        for (int J_std = jmin; J_std <= jmax; J_std++) {
          int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8) {
            int ch_J2_bc = Z.modelspace->GetTwoBodyChannelIndex(
                J_std, parity_J2, Tz_J2_bc);
            int ch_J2_ad = Z.modelspace->GetTwoBodyChannelIndex(
                J_std, parity_J2, Tz_J2_ad);

            TwoBodyChannel &tbc_J2_bc =
                Z.modelspace->GetTwoBodyChannel(ch_J2_bc);
            TwoBodyChannel &tbc_J2_ad =
                Z.modelspace->GetTwoBodyChannel(ch_J2_ad);
            int nkets_bc = tbc_J2_bc.GetNumberKets();
            int nkets_ad = tbc_J2_ad.GetNumberKets();
            if (nkets_bc < 1 or nkets_ad < 1)
              continue;

            int indx_bc = tbc_J2_bc.GetLocalIndex(std::min(int(b), int(c)),
                                                  std::max(int(b), int(c)));
            int indx_ad = tbc_J2_ad.GetLocalIndex(std::min(int(a), int(d)),
                                                  std::max(int(a), int(d)));

            if (indx_ad < 0 or indx_bc < 0)
              continue;
            if (a > d)
              indx_ad += nkets_ad;

            int indx_cb = indx_bc;
            if (b > c)
              indx_bc += nkets_bc;
            if (c > b)
              indx_cb += nkets_bc;

            int index_ch = ch_J2_ad;
            if (not Z_is_scalar) {
              if (ch_J2_ad > ch_J2_bc) {
                index_ch = ch_J2_ad;
                ch_J2_ad = ch_J2_bc;
                ch_J2_bc = index_ch;
              }
              index_ch = -1;
              for (size_t i = 0; i < nch; i++) {
                if (ch_bra_list[i] == ch_J2_ad and ch_ket_list[i] == ch_J2_bc) {
                  index_ch = i;
                }
              }
            }

            XbarIIef -=
                (2 * J_std + 1) * sixj1 * CHI_IV[ch_J2_bc](indx_ad, indx_cb);
            if (std::abs(Tz_J2_ad - Tz_J2_bc) == Z.TwoBody.rank_T)
              XbarIIIef += phaseFactor * (2 * J_std + 1) * sixj1 *
                           CHI_VII[index_ch](indx_ad, indx_bc);
          }
        }
        bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = XbarIIef;
        bar_CHI_VII_CC[ch_cc](ibra_cc, iket_cc) = XbarIIIef;
      }
      //-------------------
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  /// release memory
  for (int ch = 0; ch < nch_eta; ++ch) {
    CHI_IV[ch].clear();
  }
  for (int ch = 0; ch < nch; ++ch) {
    CHI_VII[ch].clear();
  }
  CHI_IV.clear();
  CHI_VII.clear();

  /////////////////////////////////////////////
  //     diagram    IIe and IIf
  /////////////////////////////////////////////
  std::deque<arma::mat> bar_CHI_gamma(n_nonzero); // released
  /// initial bar_CHI_V
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    // because the restriction a<b in the bar and ket vector, if we want to
    // store the full Pandya transformed matrix, we twice the size of matrix
    bar_CHI_gamma[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
  }

  // calculate bat_chi_IV * bar_gamma
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    bar_CHI_gamma[ch_cc] = bar_CHI_IV[ch_cc] * bar_Gamma[ch_cc];
  }
  // release memroy
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    bar_CHI_IV[ch_cc].clear();
    // bar_Gamma[ch_cc].clear();
  }
  bar_CHI_IV.clear();
  // bar_Gamma.clear();

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  //  Inverse Pandya transformation
  //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
  //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
  //                                    { k l J'}
  //  Diagram e
  //  II(e)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }
  //  \bar{bar_CHI_gamma}^J'_il`kj`
  //                                                                { k l J'}
  //  Diagram f
  //  II(f)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }
  //  \bar{bar_CHI_gamma_II}^J'_il`kj`
  //
#pragma omp parallel for
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
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;

      int ketmin = 0;
      if (ch_bra == ch_ket)
        ketmin = ibra;
      for (int iket = ketmin; iket < nKets; ++iket) {
        size_t k, l;
        Ket &ket = tbc_ket.GetKet(iket);
        k = ket.p;
        l = ket.q;

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commjikl = 0;
        double commijlk = 0;
        double commjilk = 0;

        // ijkl direct term
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
          // jilk, exchange i and j, k and l     ->  jk   li
          int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
          int indx_li = indx_il + (l > i ? nkets_cc : 0);
          double me1 = bar_CHI_gamma[ch_cc](indx_jk, indx_li);
          commjilk -= (2 * Jprime + 1) * sixj * me1;

          // ijkl direct term
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          me1 = bar_CHI_gamma[ch_cc](indx_il, indx_kj);
          commijkl -= (2 * Jprime + 1) * sixj * me1;
        }

        // jikl, exchange i and j    ->  jl ki
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

          // ijlk, exchange k and l     ->  ik lj
          int indx_ik = indx_ki + (i > k ? nkets_cc : 0);
          int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
          double me1 = bar_CHI_gamma[ch_cc](indx_ik, indx_lj);
          commijlk -= (2 * Jprime + 1) * sixj * me1;

          // jikl, exchange i and j    ->  jl ki
          indx_ki += (k > i ? nkets_cc : 0);
          indx_jl += (j > l ? nkets_cc : 0);
          me1 = bar_CHI_gamma[ch_cc](indx_jl, indx_ki);
          commjikl -= (2 * Jprime + 1) * sixj * me1;
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

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)

  {
    bar_CHI_gamma[ch_cc].clear();
  }
  bar_CHI_gamma.clear();

  // ######################################################
  //                 Diagram IIb and IId
  // ######################################################
  std::deque<arma::mat> CHI_III_final(n_nonzero);
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc_bra = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nbras = tbc_cc_bra.GetNumberKets();
    if (nbras < 1)
      continue;
    // Not symmetric
    CHI_III_final[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
  }

#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    CHI_III_final[ch_cc] = bar_Gamma[ch_cc] * barCHI_III_RC[ch_cc];
  }
  /// release memory
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    barCHI_III_RC[ch_cc].clear();
    bar_Gamma[ch_cc].clear();
  }
  barCHI_III_RC.clear();
  bar_Gamma.clear();

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }
  //  Inverse Pandya transformation
  //  diagram IIb and IId
  //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
  //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
  //                                    { k l J'}
#pragma omp parallel for
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
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;
      int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

      int ketmin = 0;
      if (ch_bra == ch_ket)
        ketmin = ibra;
      for (int iket = ketmin; iket < nKets; ++iket) {
        size_t k, l;
        Ket &ket = tbc_ket.GetKet(iket);
        k = ket.p;
        l = ket.q;

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commjikl = 0;
        double commijlk = 0;
        double commjilk = 0;

        // jikl, direct term        -->  jl  ki
        // ijlk, exchange ij and kl -->  lj  ik
        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;

        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);

          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;

          int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_jl < 0 or indx_ik < 0)
            continue;

          int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
          // direct term
          indx_jl += (j > l ? nkets_cc : 0);
          indx_ik += (i > k ? nkets_cc : 0);
          double me1 = CHI_III_final[ch_cc](indx_jl, indx_ik);
          commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

          int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);
          // exchange ij and kl
          double me2 = CHI_III_final[ch_cc](indx_ik, indx_jl);
          commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        // ijkl,  exchange i and j -->  il  kj
        // jilk,  exchange k and l -->  jk li
        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;

        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;

          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_jk < 0)
            continue;

          int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
          // exchange k and l
          indx_il += (i > l ? nkets_cc : 0);
          indx_jk += (j > k ? nkets_cc : 0);
          double me1 = CHI_III_final[ch_cc](indx_jk, indx_il);
          commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

          int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);
          // exchange i and j
          double me2 = CHI_III_final[ch_cc](indx_il, indx_jk);
          commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        double zijkl =
            (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;

        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  for (int ch = 0; ch < n_nonzero; ++ch) {
    CHI_III_final[ch].clear();
  }
  CHI_III_final.clear();

  /////////////////////////////////////////////
  //     diagram    IIIe and IIIf
  /////////////////////////////////////////////
  std::deque<arma::mat> bar_CHI_VII_CC_ef(n_nonzero); // released
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    // because the restriction a<b in the bar and ket vector, if we want to
    // store the full Pandya transformed matrix, we twice the size of matrix
    bar_CHI_VII_CC_ef[ch_cc] =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
  }

  // bar_CHI_VII_CC_ef = bar_CHI_VII_CC * bar_Eta
#pragma omp parallel for
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets = tbc_cc.GetNumberKets();
    if (nKets < 1)
      continue;
    bar_CHI_VII_CC_ef[ch_cc] = bar_CHI_VII_CC[ch_cc] * bar_Eta[ch_cc];
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // release bar_Eta
  // release bar_CHI_VII_CC
  for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++) {
    bar_Eta[ch_cc].clear();
    bar_CHI_VII_CC[ch_cc].clear();
  }
  bar_Eta.clear();
  bar_CHI_VII_CC.clear();

  //  Inverse Pandya transformation
  //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
  //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
  //                                    { k l J'}
#pragma omp parallel for
  for (int ch = 0; ch < nch; ++ch) {
    int ch_bra = ch_bra_list[ch];
    int ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();
    int J0 = tbc_bra.J;
    if (nbras == 0 or nkets == 0)
      continue;

    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = *(bra.op);
      Orbit &oj = *(bra.oq);
      int ji = oi.j2;
      int jj = oj.j2;
      int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

      int ketmin = 0;
      if (ch_bra == ch_ket)
        ketmin = ibra;
      for (int iket = ketmin; iket < nkets; ++iket) {
        size_t k, l;
        Ket &ket = tbc_ket.GetKet(iket);
        k = ket.p;
        l = ket.q;

        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commjikl = 0;
        double commijlk = 0;
        double commjilk = 0;

        // jikl, direct term        -->  jl  ki
        // ijlk, exchange ij and kl -->  lj  ik
        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5,
                                              jl * 0.5, Jprime);
          double sixj2 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5,
                                               jk * 0.5, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;

          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;

          int indx_lj = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_lj < 0 or indx_ik < 0)
            continue;

          // direct term
          int indx_jl = indx_lj + (j > l ? nkets_cc : 0);
          int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
          double me1 = bar_CHI_VII_CC_ef[ch_cc](indx_jl, indx_ki);
          commjikl -= (2 * Jprime + 1) * sixj * me1;

          // exchange ij and kl
          indx_ik += (i > k ? nkets_cc : 0);
          indx_lj += (l > j ? nkets_cc : 0);
          double me2 = bar_CHI_VII_CC_ef[ch_cc](indx_ik, indx_lj);
          commijlk -= (2 * Jprime + 1) * sixj2 * me2;
        }

        // ijkl,  exchange i and j -->  il  kj
        // jilk,  exchange k and l -->  jk li
        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5,
                                              jl * 0.5, Jprime);
          double sixj2 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jl * 0.5,
                                               jk * 0.5, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;

          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);

          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 0)
            continue;

          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_kj < 0)
            continue;

          // exchange k and l
          int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
          int indx_li = indx_il + (l > i ? nkets_cc : 0);
          double me2 = bar_CHI_VII_CC_ef[ch_cc](indx_jk, indx_li);
          commjilk -= (2 * Jprime + 1) * sixj2 * me2;

          // exchange i and j
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          double me1 = bar_CHI_VII_CC_ef[ch_cc](indx_il, indx_kj);
          commijkl -= (2 * Jprime + 1) * sixj * me1;
        }

        double zijkl =
            (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;

        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * 0.5 * zijkl);
      }
    }
  }

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    bar_CHI_VII_CC_ef[ch_cc].clear();
  }
  bar_CHI_VII_CC_ef.clear();

  // Timer
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return;
} //  comm223_232_chi2b

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




/// Pack AMC χ^θ into one IMSRG scalar Operator (equal-J channels only).
/// MakeNotReduced(χ_k) + MakeNotReduced(χ_j) (same slots, no transpose).
/// ChiTab stores bare S = Σ (−1)^{J0+J2+λ} λ̂^{-1} w ΩΩ.
/// AMC: χ_red = S/Ĵ, χ_unred = S/Ĵ². Pack S/Ĵ as reduced then MakeNotReduced
/// → χ_unred for Pandya / unreduced Ĵ²Ĵ² fold (matches m-scheme).
Operator ChiThetaToScalarOperator(ModelSpace &ms, const ChiTab &chi_k,
                                  const ChiTab &chi_j) {
  auto pack_one = [&](const ChiTab &chi) {
    Operator Chi(ms, 0, 0, 0, 2);
    Chi.SetNonHermitian();
    Chi.Erase();
    const int nch = ms.GetNumberTwoBodyChannels();
    for (int ch = 0; ch < nch; ++ch) {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J0 = tbc.J;
      const double hat_J0 = std::sqrt(2.0 * J0 + 1.0);
      const int nk = tbc.GetNumberKets();
      if (nk < 1)
        continue;
      arma::mat &M = Chi.TwoBody.GetMatrix(ch, ch);
      for (int ib = 0; ib < nk; ++ib) {
        Ket &bra = tbc.GetKet(ib);
        const index_t i = bra.p, j = bra.q;
        for (int ik = 0; ik < nk; ++ik) {
          Ket &ket = tbc.GetKet(ik);
          const index_t k = ket.p, l = ket.q;
          double nrm = 1.0;
          if (i == j)
            nrm *= PhysConst::SQRT2;
          if (k == l)
            nrm *= PhysConst::SQRT2;
          // χ_red = S/Ĵ then MakeNotReduced → χ_unred = S/Ĵ²
          M(ib, ik) = chi(i, j, k, l, J0) / (nrm * hat_J0);
        }
      }
    }
    Chi.is_reduced = true;
    Chi.MakeNotReduced();
    return Chi;
  };

  // two unreduced Ops, add (no M+=M.t())
  Operator Chi = pack_one(chi_k);
  Chi += pack_one(chi_j);
  return Chi;
}

/// Factorized 2n layout from non-reduced scalar χ^θ (GetTBME).
void PackChiThetaFactLayout(const Operator &Chi, Operator &Z, int ch, int J0,
                            arma::mat &Out) {
  TwoBodyChannel &tb = Z.modelspace->GetTwoBodyChannel(ch);
  const int nb = tb.GetNumberKets();
  Out = arma::mat(2 * nb, 2 * nb, arma::fill::zeros);
  if (nb < 1)
    return;
  for (int ibra = 0; ibra < nb; ++ibra) {
    Ket &bra = tb.GetKet(ibra);
    const size_t i = bra.p, j = bra.q;
    const int ji = bra.op->j2, jj = bra.oq->j2;
    for (int iket = 0; iket < nb; ++iket) {
      Ket &ket = tb.GetKet(iket);
      const size_t k = ket.p, l = ket.q;
      const int jk = ket.op->j2, jl = ket.oq->j2;
      const double me = Chi.TwoBody.GetTBME_J(J0, i, j, k, l);
      Out(ibra, iket) = me;
      if (i != j) {
        int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
        Out(ibra + nb, iket) = phase * me;
        if (k != l) {
          phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
          Out(ibra + nb, iket + nb) = phase * me;
          phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
          Out(ibra, iket + nb) = phase * me;
        }
      } else if (k != l) {
        int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
        Out(ibra, iket + nb) = phase * me;
      }
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

  // Serial: GetNineJ must not populate its cache under OMP.
  std::map<std::array<int, 2>, arma::mat> bar_Omega;
  for (int ch_b = 0; ch_b < n_cc; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    const int nb = tb.GetNumberKets();
    if (nb < 1)
      continue;
    const int Jb = tb.J;
    for (int ch_k = 0; ch_k < n_cc; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      const int nk = tk.GetNumberKets();
      if (nk < 1)
        continue;
      if (not AngMom::Triangle(Jb, tk.J, lambda))
        continue;
      if ((tb.parity + tk.parity + Eta.GetParity()) % 2 != 0 or tb.Tz != tk.Tz)
        continue;
      const int Jk = tk.J;
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
      bar_Omega[{ch_b, ch_k}] = std::move(Om);
    }
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
/// Gamma^III_b Path B: locked χ^η (normal) → Fac Pandya → code RC → Γ̄ DGEMM
/// → inv Pandya + (1−P_ij)(1−P_kl). χ is scalar and not AS / not Hermitian.
/// Gold: test_G3b_normal_to_RC.py, test_tts_GIIIb.py. Not tts_GIIIb.
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;

  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  auto &Z2 = Z.TwoBody;

  std::deque<arma::mat> Chi2_unused;
  std::deque<arma::mat> barCHI_eta;
  BuildChiEtaPathB(Eta, Z, Chi2_unused, &barCHI_eta, false);

  auto chi_unnorm = [&](int J, int i, int j, int k, int l) -> double {
    return std::sqrt(2.0 * J + 1.0) *
           InvChiEtaRed(Z, barCHI_eta, i, j, k, l, J);
  };

  std::deque<arma::mat> bar_Gamma(n_nonzero);
  std::deque<arma::mat> barCHI_III(n_nonzero);

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    bar_Gamma[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    barCHI_III[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    if (nKets_cc < 1)
      continue;

    // Scalar Pandya of Γ (hermiticity OK).
    for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc) {
      Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
      int a = bra_cc.p, b = bra_cc.q;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;

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
        double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;

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
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8) {
            Gammabar -= (2 * J_std + 1) * sixj1 *
                        Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
          }
        }

        double flip_phase =
            Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);
        if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d)) {
          bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
          if (iket_cc != ibra_cc)
            bar_Gamma[ch_cc](iket_cc, ibra_cc) = hGamma * Gammabar;
        }
        if (a != b) {
          bar_Gamma[ch_cc](ibra_cc + nKets_cc,
                           (iket_cc + nKets_cc) % (2 * nKets_cc)) =
              Gammabar * flip_phase * hGamma;
        }
        if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d)) {
          bar_Gamma[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc),
                           ibra_cc + nKets_cc) = Gammabar * flip_phase;
        }
      }
    }

    // Plain Fac Pandya of χ: all 2n orientations, NO hermiticity, NO √2.
    // χ̄_{ab,cd} ← χ_unnorm_{adcb}. Gold: test_G3b_normal_to_RC.py.
    for (int ibra_cc = 0; ibra_cc < 2 * nKets_cc; ++ibra_cc) {
      int a, b;
      if (ibra_cc < nKets_cc) {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
        a = bra_cc.p;
        b = bra_cc.q;
      } else {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
        a = bra_cc.q;
        b = bra_cc.p;
      }
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
          c = ket_cc.q;
          d = ket_cc.p;
        }
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);
        double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
        double X = 0.0;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        for (int Jp = jmin; Jp <= jmax; ++Jp) {
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, Jp);
          if (std::abs(sixj1) < 1e-8)
            continue;
          X -= (2 * Jp + 1) * sixj1 * chi_unnorm(Jp, a, d, c, b);
        }
        barCHI_III[ch_cc](ibra_cc, iket_cc) = X;
      }
    }
  }

  // RC of barCHI_III (Factorized IIb / code L1828).
  std::deque<arma::mat> barCHI_III_RC(n_nonzero);
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    int nK = Z.modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
    barCHI_III_RC[ch_cc] = arma::mat(2 * nK, 2 * nK, arma::fill::zeros);
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    if (nKets_cc < 1)
      continue;
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
        double XbarIIbd = 0.0;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        for (int J_std = jmin; J_std <= jmax; ++J_std) {
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int parity_cc = (oa.l + od.l) % 2;
          int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
          int ch_old =
              Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_old =
              Z.modelspace->GetTwoBodyChannel_CC(ch_old);
          int nkets = tbc_old.GetNumberKets();
          int indx_ad = tbc_old.GetLocalIndex(std::min(a, d), std::max(a, d));
          int indx_bc = tbc_old.GetLocalIndex(std::min(b, c), std::max(b, c));
          if (indx_ad < 0 or indx_bc < 0)
            continue;
          if (a > d)
            indx_ad += nkets;
          if (b > c)
            indx_bc += nkets;
          XbarIIbd -= Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) *
                      (2 * J_std + 1) * sixj1 *
                      (barCHI_III[ch_old](indx_bc, indx_ad) +
                       barCHI_III[ch_old](indx_ad, indx_bc));
        }
        barCHI_III_RC[ch_cc](ibra_cc, iket_cc) = XbarIIbd;
      }
    }
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    barCHI_III[ch_cc].clear();
  barCHI_III.clear();
  barCHI_eta.clear();

  // DGEMM: CHI_III_final = bar_Gamma * barCHI_III_RC
  std::deque<arma::mat> CHI_III_final(n_nonzero);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    if (bar_Gamma[ch_cc].n_rows < 1)
      continue;
    CHI_III_final[ch_cc] = bar_Gamma[ch_cc] * barCHI_III_RC[ch_cc];
  }
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    bar_Gamma[ch_cc].clear();
    barCHI_III_RC[ch_cc].clear();
  }
  bar_Gamma.clear();
  barCHI_III_RC.clear();

  // Inverse Pandya + (1−P_ij)(1−P_kl) (Factorized IIb/IId)
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    int ch_bra = (int)ch_bra_list[ch];
    int ch_ket = (int)ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();
    int J0 = tbc_bra.J;
    if (nbras == 0 or nkets == 0)
      continue;

    for (int ibra = 0; ibra < (int)nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p, j = bra.q;
      Orbit &oi = *(bra.op), &oj = *(bra.oq);
      int ji = oi.j2, jj = oj.j2;
      int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);
      int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < (int)nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p, l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2, jl = ol.j2;
        double commijkl = 0, commjikl = 0, commijlk = 0, commjilk = 0;

        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_jl < 0 or indx_ik < 0)
            continue;
          int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
          indx_jl += (j > l ? nkets_cc : 0);
          indx_ik += (i > k ? nkets_cc : 0);
          commjikl -= phase1 * (2 * Jprime + 1) * sixj1 *
                      CHI_III_final[ch_cc](indx_jl, indx_ik);
          int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);
          commijlk -= phase2 * (2 * Jprime + 1) * sixj1 *
                      CHI_III_final[ch_cc](indx_ik, indx_jl);
        }

        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5,
                                               jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          int ch_cc =
              Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_jk < 0)
            continue;
          int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
          indx_il += (i > l ? nkets_cc : 0);
          indx_jk += (j > k ? nkets_cc : 0);
          commjilk -= phase1 * (2 * Jprime + 1) * sixj1 *
                      CHI_III_final[ch_cc](indx_jk, indx_il);
          int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);
          commijkl -= phase2 * (2 * Jprime + 1) * sixj1 *
                      CHI_III_final[ch_cc](indx_il, indx_jk);
        }

        double zijkl =
            (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
      }
    }
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    CHI_III_final[ch_cc].clear();
  CHI_III_final.clear();
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
  // χ^θ = MakeNotReduced(χ_k) + MakeNotReduced(χ_j) (same slots)
  Operator Chi_theta = ChiThetaToScalarOperator(*Z.modelspace, chi_k, chi_j);

  // Gamma non-reduced for scalar Pandya
  const Operator *Gp = &Gamma;
  Operator Gtmp;
  if (Gamma.IsReduced() and Gamma.GetJRank() == 0) {
    Gtmp = Gamma;
    Gtmp.MakeNotReduced();
    Gp = &Gtmp;
  }

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  const int nch_eta = Z.modelspace->GetNumberTwoBodyChannels();
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  auto &Z2 = Z.TwoBody;

  // CHI_IV from non-reduced χ^θ (scalar Factorized IIe layout)
  std::deque<arma::mat> CHI_IV(nch_eta);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch_eta; ++ch) {
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    PackChiThetaFactLayout(Chi_theta, Z, ch, tbc.J, CHI_IV[ch]);
  }

  // Pandya Γ̄ and χ̄^θ
  std::deque<arma::mat> bar_Gamma(n_nonzero);
  std::deque<arma::mat> bar_CHI_IV(n_nonzero);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;
    bar_Gamma[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    bar_CHI_IV[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
    if (nKets_cc < 1)
      continue;

    // bar_Gamma (same Pandya as Factorized)
    for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc) {
      Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
      int a = bra_cc.p, b = bra_cc.q;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
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
        double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
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
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8)
            Gammabar -= (2 * J_std + 1) * sixj1 *
                        Gp->TwoBody.GetTBME_J(J_std, a, d, c, b);
        }
        double flip_phase =
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
          int ch_J2_bc = Z.modelspace->GetTwoBodyChannelIndex(
              J_std, parity_J2, Tz_J2_bc);
          TwoBodyChannel &tbc_J2_bc = Z.modelspace->GetTwoBodyChannel(ch_J2_bc);
          int nkets_bc = tbc_J2_bc.GetNumberKets();
          if (nkets_bc < 1)
            continue;
          int indx_bc = tbc_J2_bc.GetLocalIndex(std::min(b, c), std::max(b, c));
          int indx_ad = tbc_J2_bc.GetLocalIndex(std::min(a, d), std::max(a, d));
          // Factorized uses ch_J2_ad for indx_ad — for scalar Z, Tz_ad may
          // differ; match Factorized: separate ch for ad.
          int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
          int ch_J2_ad = Z.modelspace->GetTwoBodyChannelIndex(
              J_std, parity_J2, Tz_J2_ad);
          TwoBodyChannel &tbc_J2_ad = Z.modelspace->GetTwoBodyChannel(ch_J2_ad);
          int nkets_ad = tbc_J2_ad.GetNumberKets();
          if (nkets_ad < 1)
            continue;
          indx_ad = tbc_J2_ad.GetLocalIndex(std::min(a, d), std::max(a, d));
          indx_bc = tbc_J2_bc.GetLocalIndex(std::min(b, c), std::max(b, c));
          if (indx_ad < 0 or indx_bc < 0)
            continue;
          if (a > d)
            indx_ad += nkets_ad;
          int indx_cb = indx_bc;
          if (b > c)
            indx_bc += nkets_bc;
          if (c > b)
            indx_cb += nkets_bc;
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
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    bar_CHI_IV[ch_cc].clear();
    bar_Gamma[ch_cc].clear();
  }
  bar_CHI_IV.clear();
  bar_Gamma.clear();

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


/// Round-trip test: Ω → Pandya → inv Pandya ≟ Ω
/// (A) TensorCommutators::DoTensorPandya + AddInverseTensorPandya
/// (B) ethS fIIIa NineJ forward (legs ad,bc) + matching NineJ inverse
void DebugTensorPandyaRoundTrip(const Operator &Omega) {
  ModelSpace &ms = *Omega.modelspace;
  ms.PreCalculateSixJ();
  ms.PreCalculateNineJ();
  const int lambda = Omega.GetJRank();
  const int hO = Omega.IsHermitian() ? 1 : -1;
  const double n0 = Omega.TwoBodyNorm();

  auto report = [&](const char *tag, double n_rec, double n_diff, double max_abs) {
    double rel = (n0 > 1e-14) ? n_diff / n0 : n_diff;
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  [" << tag << "] ||Ω||=" << n0 << " ||Ω_rec||=" << n_rec
              << " ||Δ||=" << n_diff << " rel=" << rel << " max|Δ|=" << max_abs
              << (rel < 1e-8 ? "  PASS" : "  FAIL") << "\n";
  };

  auto diff_ops = [&](const Operator &Orec) {
    double n_rec = Orec.TwoBodyNorm();
    double n_diff = 0.0, max_abs = 0.0;
    for (auto &it : Omega.TwoBody.MatEl) {
      auto jt = Orec.TwoBody.MatEl.find(it.first);
      if (jt == Orec.TwoBody.MatEl.end())
        continue;
      arma::mat d = it.second - jt->second;
      n_diff += arma::accu(d % d);
      max_abs = std::max(max_abs, arma::abs(d).max());
    }
    return std::make_tuple(n_rec, std::sqrt(n_diff), max_abs);
  };

  std::cout << "==== DebugTensorPandyaRoundTrip (λ=" << lambda
            << ", reduced=" << Omega.IsReduced() << ", h=" << hO << ") ====\n";

  // ------------------------------------------------------------------
  // (A) TensorCommutators Pandya is *ph-restricted* (hh+ph bras only).
  //     A full-Ω round-trip is not defined for that API — skip AddInverse.
  // ------------------------------------------------------------------
  {
    std::map<std::array<index_t, 2>, arma::mat> Zbar;
    DoTensorPandyaTransformation(Omega, Zbar);
    size_t nbar = 0;
    for (auto &kv : Zbar)
      nbar += kv.second.n_elem;
    std::cout << "  [TensorCommutators] fwd OK (ph-restricted): blocks="
              << Zbar.size() << " elems=" << nbar
              << "  (no full-Omega inv round-trip)\n";
  }

  // ------------------------------------------------------------------
  // (B) ethS fIIIa NineJ (legs Ω_{ad,bc}): forward then exact inverse
  //     Forward:  bar(a,b;c,d)^{Jp Jq} = Σ_{J1 J2} hat 9j phase Ω^{J1 J2}(a,d,b,c)
  //     Inverse:  Ω_rec^{J1 J2}(a,d,b,c) = Σ_{Jp Jq} hat 9j phase bar(a,b;c,d)
  //     (same hat/phase both ways; 9j orthogonality ⇒ identity)
  // ------------------------------------------------------------------
  {
    int n_cc = ms.GetNumberTwoBodyChannels_CC();
    std::vector<std::vector<arma::mat>> bar(n_cc, std::vector<arma::mat>(n_cc));

    auto fwd_me = [&](int a, int b, int c, int d, int Jbra, int Jket) -> double {
      Orbit &oa = ms.GetOrbit(a);
      Orbit &ob = ms.GetOrbit(b);
      Orbit &oc = ms.GetOrbit(c);
      Orbit &od = ms.GetOrbit(d);
      double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
      double Xbar = 0.0;
      if (lambda == 0) {
        if (Jbra != Jket)
          return 0.0;
        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
        int dJ = 1;
        if ((a == d or b == c)) {
          dJ = 2;
          jmin += jmin % 2;
        }
        for (int J1 = jmin; J1 <= jmax; J1 += dJ) {
          double sixj = AngMom::SixJ(ja, jb, Jbra, jc, jd, J1);
          if (std::abs(sixj) < 1e-10)
            continue;
          double phase = ms.phase((ob.j2 + oc.j2) / 2 + J1);
          Xbar += (2 * J1 + 1) * sixj * phase *
                  Omega.TwoBody.GetTBME_J(J1, a, d, b, c);
        }
        return Xbar;
      }
      if (not AngMom::Triangle(Jbra, Jket, lambda))
        return 0.0;
      int j1min = std::abs(oa.j2 - od.j2) / 2;
      int j1max = (oa.j2 + od.j2) / 2;
      for (int J1 = j1min; J1 <= j1max; ++J1) {
        int j2min =
            std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
        int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
        for (int J2 = j2min; J2 <= j2max; ++J2) {
          double ninej =
              ms.GetNineJ(ja, jd, J1, jb, jc, J2, Jbra, Jket, lambda);
          if (std::abs(ninej) < 1e-10)
            continue;
          double hat = std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) *
                                 (2.0 * Jbra + 1.0) * (2.0 * Jket + 1.0));
          double phase =
              ms.phase((ob.j2 + oc.j2) / 2 + J2 + Jket - Jbra);
          Xbar += hat * ninej * phase *
                  Omega.TwoBody.GetTBME_J(J1, J2, a, d, b, c);
        }
      }
      return Xbar;
    };

    // Fill bar for all CC channel pairs
    for (int chb = 0; chb < n_cc; ++chb) {
      TwoBodyChannel_CC &tb = ms.GetTwoBodyChannel_CC(chb);
      int nb = tb.GetNumberKets();
      if (nb < 1)
        continue;
      int Jbra = tb.J;
      int n2b = nb * 2;
      for (int chk = 0; chk < n_cc; ++chk) {
        TwoBodyChannel_CC &tk = ms.GetTwoBodyChannel_CC(chk);
        int nk = tk.GetNumberKets();
        if (nk < 1)
          continue;
        int Jket = tk.J;
        if (lambda == 0 and Jbra != Jket)
          continue;
        if (lambda != 0 and not AngMom::Triangle(Jbra, Jket, lambda))
          continue;
        if ((tb.parity + tk.parity + Omega.GetParity()) % 2)
          continue;
        int n2k = nk * 2;
        arma::mat M(n2b, n2k, arma::fill::zeros);
        for (int ib = 0; ib < n2b; ++ib) {
          int a, b;
          if (ib < nb) {
            a = tb.GetKet(ib).p;
            b = tb.GetKet(ib).q;
          } else {
            b = tb.GetKet(ib - nb).p;
            a = tb.GetKet(ib - nb).q;
          }
          if (ib >= nb and a == b)
            continue;
          for (int ik = 0; ik < n2k; ++ik) {
            int c, d;
            if (ik < nk) {
              c = tk.GetKet(ik).p;
              d = tk.GetKet(ik).q;
            } else {
              d = tk.GetKet(ik - nk).p;
              c = tk.GetKet(ik - nk).q;
            }
            if (ik >= nk and c == d)
              continue;
            Orbit &oa = ms.GetOrbit(a);
            Orbit &ob = ms.GetOrbit(b);
            Orbit &oc = ms.GetOrbit(c);
            Orbit &od = ms.GetOrbit(d);
            if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Omega.GetTRank())
              continue;
            M(ib, ik) = fwd_me(a, b, c, d, Jbra, Jket);
          }
        }
        bar[chb][chk] = std::move(M);
      }
    }

    // Compare via GetTBME on all orbits (avoids bra=(a,d) confusion)
    double n_diff = 0.0, max_abs = 0.0, n_rec2 = 0.0;
    int max_J = ms.GetTwoBodyJmax();
    for (auto i : ms.all_orbits)
      for (auto j : ms.all_orbits)
        for (auto k : ms.all_orbits)
          for (auto l : ms.all_orbits)
            for (int J1 = 0; J1 <= max_J; ++J1)
              for (int J2 = 0; J2 <= max_J; ++J2) {
                if (lambda != 0 and not AngMom::Triangle(J1, J2, lambda))
                  continue;
                if (lambda == 0 and J1 != J2)
                  continue;
                double o0 = Omega.TwoBody.GetTBME_J(J1, J2, i, j, k, l);
                // Reconstruct from bar(i,k,l,j) corresponding to Ω(i,j,k,l)
                // with forward legs (a,d,b,c)=(i,j,k,l) → bar(a,b,c,d)=bar(i,k,l,j)
                Orbit &oi = ms.GetOrbit(i);
                Orbit &oj = ms.GetOrbit(j);
                Orbit &ok = ms.GetOrbit(k);
                Orbit &ol = ms.GetOrbit(l);
                double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
                double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
                double sm = 0.0;
                int parity_b = (oi.l + ok.l) % 2;
                int parity_k = (ol.l + oj.l) % 2;
                int Tz_b = std::abs(oi.tz2 - ok.tz2) / 2;
                int Tz_k = std::abs(ol.tz2 - oj.tz2) / 2;
                int Jpmin = std::abs(oi.j2 - ok.j2) / 2;
                int Jpmax = (oi.j2 + ok.j2) / 2;
                for (int Jp = Jpmin; Jp <= Jpmax; ++Jp) {
                  int chb = ms.GetTwoBodyChannelIndex(Jp, parity_b, Tz_b);
                  if (chb < 0 or chb >= n_cc)
                    continue;
                  TwoBodyChannel_CC &tb = ms.GetTwoBodyChannel_CC(chb);
                  int nb = tb.GetNumberKets();
                  int idx_ab =
                      tb.GetLocalIndex(std::min(i, k), std::max(i, k));
                  if (idx_ab < 0)
                    continue;
                  idx_ab += (i > k ? nb : 0);
                  int Jqmin =
                      std::max(std::abs(ol.j2 - oj.j2) / 2, std::abs(Jp - lambda));
                  int Jqmax = std::min((ol.j2 + oj.j2) / 2, Jp + lambda);
                  if (lambda == 0) {
                    Jqmin = Jp;
                    Jqmax = Jp;
                  }
                  for (int Jq = Jqmin; Jq <= Jqmax; ++Jq) {
                    int chk = ms.GetTwoBodyChannelIndex(Jq, parity_k, Tz_k);
                    if (chk < 0 or chk >= n_cc)
                      continue;
                    if (bar[chb][chk].n_rows == 0)
                      continue;
                    TwoBodyChannel_CC &tk = ms.GetTwoBodyChannel_CC(chk);
                    int nk = tk.GetNumberKets();
                    int idx_cd =
                        tk.GetLocalIndex(std::min(l, j), std::max(l, j));
                    if (idx_cd < 0)
                      continue;
                    idx_cd += (l > j ? nk : 0);
                    if (idx_ab >= (int)bar[chb][chk].n_rows or
                        idx_cd >= (int)bar[chb][chk].n_cols)
                      continue;
                    double tbme = bar[chb][chk](idx_ab, idx_cd);
                    if (lambda == 0) {
                      double sixj =
                          AngMom::SixJ(ji, jk, Jp, jl, jj, J1);
                      if (std::abs(sixj) < 1e-10)
                        continue;
                      double phase =
                          ms.phase((ok.j2 + ol.j2) / 2 + J1);
                      sm += (2 * Jp + 1) * sixj * phase * tbme;
                    } else {
                      double ninej = ms.GetNineJ(ji, jj, J1, jk, jl, J2, Jp,
                                                   Jq, lambda);
                      if (std::abs(ninej) < 1e-10)
                        continue;
                      double hat =
                          std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) *
                                    (2.0 * Jp + 1.0) * (2.0 * Jq + 1.0));
                      double phase =
                          ms.phase((ok.j2 + ol.j2) / 2 + J2 + Jq - Jp);
                      sm += hat * ninej * phase * tbme;
                    }
                  }
                }
                double d = o0 - sm;
                n_diff += d * d;
                n_rec2 += sm * sm;
                max_abs = std::max(max_abs, std::abs(d));
              }
    report("ethS NineJ ad,bc (exact inv)", std::sqrt(n_rec2), std::sqrt(n_diff),
           max_abs);
  }
}

/// Ordinary-channel hermiticity of χ vs Pandya (CC) matrix hermiticity.
/// Hypothesis: χ_k+χ_j is Hermitian in pp J-scheme, but χ̄ may not be
/// Hermitian as a 2n CC matrix (Path B never fills with h_χ).
void DebugChiPandyaHermiticity(const Operator &Eta) {
  ModelSpace &ms = *Eta.modelspace;
  ms.PreCalculateSixJ();
  const int lambda = Eta.GetJRank();
  const double hat_lambda_inv =
      1.0 / std::sqrt(2.0 * std::max(lambda, 0) + 1.0);
  int max_j2 = 0;
  for (auto x : ms.all_orbits)
    max_j2 = std::max(max_j2, ms.GetOrbit(x).j2);
  std::vector<index_t> allorb(ms.all_orbits.begin(), ms.all_orbits.end());
  const int n_orb = (int)allorb.size();

  Operator Zdummy(ms, 0, 0, 0, 2);
  Zdummy.SetHermitian();
  ChiTab chi_k, chi_j;
  chi_k.allocate(n_orb, max_j2);
  chi_j.allocate(n_orb, max_j2);
  FillChiThetaG3c_DGEMM(Eta, Zdummy, chi_k, chi_j, hat_lambda_inv, lambda,
                        allorb);

  auto ordinary_herm = [&](const Operator &Chi, const char *name) {
    double num = 0.0, den = 0.0;
    const int nch = ms.GetNumberTwoBodyChannels();
    for (int ch = 0; ch < nch; ++ch) {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J0 = tbc.J;
      const int nk = tbc.GetNumberKets();
      for (int ib = 0; ib < nk; ++ib) {
        Ket &bra = tbc.GetKet(ib);
        for (int ik = 0; ik < nk; ++ik) {
          Ket &ket = tbc.GetKet(ik);
          const double v =
              Chi.TwoBody.GetTBME_J(J0, bra.p, bra.q, ket.p, ket.q);
          const double w =
              Chi.TwoBody.GetTBME_J(J0, ket.p, ket.q, bra.p, bra.q);
          num += (v - w) * (v - w);
          den += v * v;
        }
      }
    }
    const double r = (den > 0) ? std::sqrt(num / den) : 0.0;
    std::cout << "  ordinary " << name << ": ‖χ-χᵀ‖/‖χ‖=" << r
              << (r < 1e-8 ? "  Hermitian" : "  NOT Hermitian") << std::endl;
    return r;
  };

  auto pandya_herm = [&](const Operator &Chi, const char *name) {
    const int nch_eta = ms.GetNumberTwoBodyChannels();
    const int n_nonzero = ms.GetNumberTwoBodyChannels_CC();
    std::deque<arma::mat> CHI_IV(nch_eta);
    for (int ch = 0; ch < nch_eta; ++ch) {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      PackChiThetaFactLayout(Chi, Zdummy, ch, tbc.J, CHI_IV[ch]);
    }
    double num = 0.0, den = 0.0;
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
      TwoBodyChannel_CC &tbc_cc = ms.GetTwoBodyChannel_CC(ch_cc);
      const int nKets_cc = tbc_cc.GetNumberKets();
      const int J_cc = tbc_cc.J;
      if (nKets_cc < 1)
        continue;
      arma::mat bar(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
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
        Orbit &oa = ms.GetOrbit(a);
        Orbit &ob = ms.GetOrbit(b);
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
          Orbit &oc = ms.GetOrbit(c);
          Orbit &od = ms.GetOrbit(d);
          double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
          int Tz_J2_bc = (ob.tz2 + oc.tz2) / 2;
          int parity_J2 = (ob.l + oc.l) % 2;
          int jmin =
              std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
          int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
          double Xbar = 0.0;
          for (int J_std = jmin; J_std <= jmax; ++J_std) {
            double sixj1 = ms.GetSixJ(ja, jb, J_cc, jc, jd, J_std);
            if (std::abs(sixj1) < 1e-8)
              continue;
            int ch_J2_bc =
                ms.GetTwoBodyChannelIndex(J_std, parity_J2, Tz_J2_bc);
            int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
            int ch_J2_ad =
                ms.GetTwoBodyChannelIndex(J_std, parity_J2, Tz_J2_ad);
            TwoBodyChannel &tbc_J2_bc = ms.GetTwoBodyChannel(ch_J2_bc);
            TwoBodyChannel &tbc_J2_ad = ms.GetTwoBodyChannel(ch_J2_ad);
            int nkets_bc = tbc_J2_bc.GetNumberKets();
            int nkets_ad = tbc_J2_ad.GetNumberKets();
            if (nkets_bc < 1 or nkets_ad < 1)
              continue;
            int indx_ad =
                tbc_J2_ad.GetLocalIndex(std::min(a, d), std::max(a, d));
            int indx_bc =
                tbc_J2_bc.GetLocalIndex(std::min(b, c), std::max(b, c));
            if (indx_ad < 0 or indx_bc < 0)
              continue;
            if (a > d)
              indx_ad += nkets_ad;
            int indx_cb = indx_bc;
            if (b > c)
              indx_bc += nkets_bc;
            if (c > b)
              indx_cb += nkets_bc;
            Xbar -=
                (2 * J_std + 1) * sixj1 * CHI_IV[ch_J2_bc](indx_ad, indx_cb);
          }
          bar(ibra_cc, iket_cc) = Xbar;
        }
      }
      arma::mat asym = bar - bar.t();
      num += arma::accu(asym % asym);
      den += arma::accu(bar % bar);
    }
    const double r = (den > 0) ? std::sqrt(num / den) : 0.0;
    std::cout << "  Pandya CC " << name << ": ‖χ̄-χ̄ᵀ‖/‖χ̄‖=" << r
              << (r < 1e-8 ? "  Hermitian in CC" : "  NOT Hermitian in CC")
              << std::endl;
    return r;
  };

  std::cout << "==== DebugChiPandyaHermiticity (λ_Ω=" << lambda
            << ") ====" << std::endl;
  Operator Cs = ChiThetaToScalarOperator(ms, chi_k, chi_j);
  ordinary_herm(Cs, "χ_k+χ_j");
  pandya_herm(Cs, "χ_k+χ_j");
  std::cout << "Note: Path B fills χ̄ elementwise with NO h_χ "
               "(Γ̄ uses h_Γ). Naive CC matrix hermiticity is not assumed."
            << std::endl;
  std::cout << "===========================================" << std::endl;
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
void comm223_232_GIVa(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const double hat_lam_inv =
      (lambda == 0) ? 1.0 : 1.0 / std::sqrt(2.0 * lambda + 1.0);
  auto hatJ = [](int J) { return std::sqrt(2.0 * J + 1.0); };
  auto &Z2 = Z.TwoBody;
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  const int nch_tb = Z.modelspace->GetNumberTwoBodyChannels();

  auto pq_cc = [](TwoBodyChannel_CC &tbc, int idx, int nK) -> std::array<int, 2> {
    Ket &ket = tbc.GetKet(idx % nK);
    if (idx < nK)
      return {(int)ket.p, (int)ket.q};
    if (ket.p == ket.q)
      return {-1, -1};
    return {(int)ket.q, (int)ket.p};
  };
  auto pq_tb = [](TwoBodyChannel &tbc, int idx, int nK) -> std::array<int, 2> {
    Ket &ket = tbc.GetKet(idx % nK);
    if (idx < nK)
      return {(int)ket.p, (int)ket.q};
    if (ket.p == ket.q)
      return {-1, -1};
    return {(int)ket.q, (int)ket.p};
  };

  // IMSRG DoTensorPandya (adcb) — hats × NineJ, no extra scale. λ=0 ⇒ Jbra=Jket.
  auto pandya_eta = [&](int a, int b, int c, int d, int Jbra,
                        int Jket) -> double {
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &ob = Z.modelspace->GetOrbit(b);
    Orbit &oc = Z.modelspace->GetOrbit(c);
    Orbit &od = Z.modelspace->GetOrbit(d);
    const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
    const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
    if (not AngMom::Triangle(Jbra, Jket, lambda))
      return 0.0;
    double Xbar = 0.0;
    const int j1min = std::abs(oa.j2 - od.j2) / 2;
    const int j1max = (oa.j2 + od.j2) / 2;
    for (int J1 = j1min; J1 <= j1max; ++J1) {
      const int j2min =
          std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
      const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
      for (int J2 = j2min; J2 <= j2max; ++J2) {
        const double ninej = Z.modelspace->GetNineJ(
            ja, jd, J1, jb, jc, J2, Jbra, Jket, lambda);
        if (std::abs(ninej) < 1e-14)
          continue;
        const double hats = hatJ(J1) * hatJ(J2) * hatJ(Jbra) * hatJ(Jket);
        const double tbme = Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
        Xbar -= hats * Z.modelspace->phase((ob.j2 + od.j2) / 2 + Jket + J2) *
                ninej * tbme;
      }
    }
    return Xbar;
  };

  // Scalar Pandya Γ per CC channel (Factorized 2×nKets layout).
  std::deque<arma::mat> bar_Gamma(n_nonzero);
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

  // Rectangular CC pairs: χ̄(il;kj)^{J0,J1} = hΩ (−1)^{J0+J1} Ω̄_occ^T · Γ̄^{J1}
  std::vector<std::array<int, 2>> cc_pairs;
  for (int ch_b = 0; ch_b < n_nonzero; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    if (tb.GetNumberKets() < 1)
      continue;
    for (int ch_k = 0; ch_k < n_nonzero; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      if (tk.GetNumberKets() < 1)
        continue;
      if ((tb.parity + tk.parity) % 2 != 0)
        continue;
      if (tb.Tz != tk.Tz)
        continue;
      if (not AngMom::Triangle(tb.J, tk.J, lambda))
        continue;
      cc_pairs.push_back({ch_b, ch_k});
    }
  }

  std::map<std::array<int, 2>, arma::mat> bar_CHI_VI_II;
  for (auto pk : cc_pairs) {
    const int ch_ab = pk[0], ch_il = pk[1];
    TwoBodyChannel_CC &tab = Z.modelspace->GetTwoBodyChannel_CC(ch_ab);
    TwoBodyChannel_CC &til = Z.modelspace->GetTwoBodyChannel_CC(ch_il);
    const int nab = tab.GetNumberKets(), nil = til.GetNumberKets();
    const int J1 = tab.J, J0 = til.J;
    arma::mat Omega_occ(2 * nab, 2 * nil, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nab; ++ibra) {
      auto ab = pq_cc(tab, ibra, nab);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double n_a = oa.occ, nbar_a = 1.0 - n_a;
      const double n_b = ob.occ, nbar_b = 1.0 - n_b;
      for (int iket = 0; iket < 2 * nil; ++iket) {
        auto il = pq_cc(til, iket, nil);
        if (il[0] < 0)
          continue;
        Orbit &od = Z.modelspace->GetOrbit(il[1]);
        const double n_d = od.occ, nbar_d = 1.0 - n_d;
        const double occ_ABbarD = n_a * nbar_b * n_d + nbar_a * n_b * nbar_d;
        if (std::abs(occ_ABbarD) < 1e-12)
          continue;
        Omega_occ(ibra, iket) =
            pandya_eta(ab[0], ab[1], il[0], il[1], J1, J0) * occ_ABbarD;
      }
    }
    const double phaseJJ = Z.modelspace->phase(J0 + J1);
    bar_CHI_VI_II[{ch_il, ch_ab}] =
        hEta * phaseJJ * (Omega_occ.t() * bar_Gamma[ch_ab]);
  }

  auto bar_chi = [&](int i, int j, int k, int l, int J2, int J3) -> double {
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    if (not AngMom::Triangle(J2, J3, lambda))
      return 0.0;
    if (not AngMom::Triangle(oi.j2 * 0.5, ol.j2 * 0.5, (double)J2) or
        not AngMom::Triangle(ok.j2 * 0.5, oj.j2 * 0.5, (double)J3))
      return 0.0;
    const int ch2 = Z.modelspace->GetTwoBodyChannelIndex(
        J2, (oi.l + ol.l) % 2, std::abs(oi.tz2 - ol.tz2) / 2);
    const int ch3 = Z.modelspace->GetTwoBodyChannelIndex(
        J3, (ok.l + oj.l) % 2, std::abs(ok.tz2 - oj.tz2) / 2);
    if (ch2 < 0 or ch3 < 0 or ch2 >= n_nonzero or ch3 >= n_nonzero)
      return 0.0;
    auto it = bar_CHI_VI_II.find({ch2, ch3});
    if (it == bar_CHI_VI_II.end())
      return 0.0;
    TwoBodyChannel_CC &t2 = Z.modelspace->GetTwoBodyChannel_CC(ch2);
    TwoBodyChannel_CC &t3 = Z.modelspace->GetTwoBodyChannel_CC(ch3);
    const int n2 = t2.GetNumberKets(), n3 = t3.GetNumberKets();
    int indx_il = t2.GetLocalIndex(std::min(i, l), std::max(i, l));
    int indx_kj = t3.GetLocalIndex(std::min(k, j), std::max(k, j));
    if (indx_il < 0 or indx_kj < 0)
      return 0.0;
    indx_il += (i > l ? n2 : 0);
    indx_kj += (k > j ? n3 : 0);
    if (indx_il >= (int)it->second.n_rows or indx_kj >= (int)it->second.n_cols)
      return 0.0;
    return it->second(indx_il, indx_kj);
  };

  // AMC tensor inv WITHOUT leading minus (Path B gold / IMSRG via Eq4).
  auto inv_plus = [&](int i, int j, int k, int l, int J0, int J1) -> double {
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    if (not AngMom::Triangle(J0, J1, lambda))
      return 0.0;
    if (not AngMom::Triangle(oi.j2 * 0.5, oj.j2 * 0.5, (double)J0) or
        not AngMom::Triangle(ok.j2 * 0.5, ol.j2 * 0.5, (double)J1))
      return 0.0;
    double tot = 0.0;
    const int max_Jloc =
        std::max({oi.j2, oj.j2, ok.j2, ol.j2, 2 * lambda}) + lambda + 2;
    for (int J2 = 0; J2 <= max_Jloc; ++J2) {
      for (int J3 = 0; J3 <= max_Jloc; ++J3) {
        if (not AngMom::Triangle(J2, J3, lambda))
          continue;
        const double bc = bar_chi(i, j, k, l, J2, J3);
        if (std::abs(bc) < 1e-16)
          continue;
        const double nj = Z.modelspace->GetNineJ(
            lambda, J0, J1, J3, oj.j2 * 0.5, ok.j2 * 0.5, J2, oi.j2 * 0.5,
            ol.j2 * 0.5);
        if (std::abs(nj) < 1e-16)
          continue;
        tot += Z.modelspace->phase(J2) * hatJ(J2) * hatJ(J3) * nj * bc;
      }
    }
    return Z.modelspace->phase(J0 + (oi.j2 + ok.j2) / 2 + lambda) * hatJ(J0) *
           hatJ(J1) * tot;
  };

  // Ordinary pair channels: χ^{J0 J2}_{ij,bd} × (−1)^{jb+jd+λ}  (2n × 2n)
  std::vector<std::array<int, 2>> tb_pairs;
  for (int c0 = 0; c0 < nch_tb; ++c0) {
    TwoBodyChannel &t0 = Z.modelspace->GetTwoBodyChannel(c0);
    if (t0.GetNumberKets() < 1)
      continue;
    for (int c2 = 0; c2 < nch_tb; ++c2) {
      TwoBodyChannel &t2 = Z.modelspace->GetTwoBodyChannel(c2);
      if (t2.GetNumberKets() < 1)
        continue;
      if (not AngMom::Triangle(t0.J, t2.J, lambda))
        continue;
      tb_pairs.push_back({c0, c2});
    }
  }

  std::map<std::array<int, 2>, arma::mat> Chi;
  for (auto pk : tb_pairs) {
    const int c0 = pk[0], c2 = pk[1];
    TwoBodyChannel &t0 = Z.modelspace->GetTwoBodyChannel(c0);
    TwoBodyChannel &t2 = Z.modelspace->GetTwoBodyChannel(c2);
    const int n0 = t0.GetNumberKets(), n2 = t2.GetNumberKets();
    const int J0 = t0.J, J2 = t2.J;
    arma::mat C(2 * n0, 2 * n2, arma::fill::zeros);
    for (int ib = 0; ib < 2 * n0; ++ib) {
      auto ij = pq_tb(t0, ib, n0);
      if (ij[0] < 0)
        continue;
      for (int ik = 0; ik < 2 * n2; ++ik) {
        auto bd = pq_tb(t2, ik, n2);
        if (bd[0] < 0)
          continue;
        const double v = inv_plus(ij[0], ij[1], bd[0], bd[1], J0, J2);
        if (std::abs(v) < 1e-16)
          continue;
        Orbit &ob = Z.modelspace->GetOrbit(bd[0]);
        Orbit &od = Z.modelspace->GetOrbit(bd[1]);
        C(ib, ik) = v * Z.modelspace->phase((ob.j2 + od.j2) / 2 + lambda);
      }
    }
    Chi[pk] = std::move(C);
  }

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

  // W_red = −(−1)^{J0}/Ĵ0 · λ̂^{-1} Σ_{J2} Chi^{J0 J2} · Ω^{J2 J0}(db;kl)
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
    if (nbras < 1 or nkets < 1)
      continue;

    arma::mat W2(2 * nbras, 2 * nkets, arma::fill::zeros);
    arma::mat W2T(2 * nkets, 2 * nbras, arma::fill::zeros); // W_klij layout
    for (int c2 = 0; c2 < nch_tb; ++c2) {
      TwoBodyChannel &t2 = Z.modelspace->GetTwoBodyChannel(c2);
      const int n2 = t2.GetNumberKets();
      const int J2 = t2.J;
      if (n2 < 1 or not AngMom::Triangle(J0, J2, lambda))
        continue;
      auto itC = Chi.find({(int)ch_bra, c2});
      auto itCT = Chi.find({(int)ch_ket, c2});
      arma::mat Om(2 * n2, 2 * nkets, arma::fill::zeros);
      arma::mat OmT(2 * n2, 2 * nbras, arma::fill::zeros);
      for (int ik = 0; ik < 2 * n2; ++ik) {
        auto bd = pq_tb(t2, ik, n2);
        if (bd[0] < 0)
          continue;
        for (int il = 0; il < 2 * nkets; ++il) {
          auto kl = pq_tb(tbc_ket, il, nkets);
          if (kl[0] < 0)
            continue;
          Om(ik, il) = Eta.TwoBody.GetTBME_J(J2, J0, bd[1], bd[0], kl[0], kl[1]);
        }
        for (int il = 0; il < 2 * nbras; ++il) {
          auto ij = pq_tb(tbc_bra, il, nbras);
          if (ij[0] < 0)
            continue;
          OmT(ik, il) =
              Eta.TwoBody.GetTBME_J(J2, J0, bd[1], bd[0], ij[0], ij[1]);
        }
      }
      if (itC != Chi.end())
        W2 += itC->second * Om;
      if (itCT != Chi.end())
        W2T += itCT->second * OmT;
    }
    const double pref = -Z.modelspace->phase(J0) / hatJ(J0) * hat_lam_inv;
    W2 *= pref;
    W2T *= pref;

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
        const double W = W2(ibra, iket);
        const double W_ji = W2(ibra_s, iket);
        const double W_kl = W2T(iket, ibra);
        const double W_lk = W2T(iket_s, ibra);
        const double Zred = (W - Pij * W_ji) + (W_kl - Pkl * W_lk);
        double z = Zred / hatJ(J0);
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

////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_b / chi^iota — Path B for all λ (λ=0 = equal-J, not a fork).
/// Scalar Factorized flow: Pandya χ̄ → RC pack → Ω̄·RC → Inv (1−P)^2.
/// Tensor: GIVa 9j Pandya Ω (no extra scale); NineJ RC/Inv that λ→0 is 6j.
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVb(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const int hZ = hGamma;
  auto hatJ = [](int J) { return std::sqrt(2.0 * J + 1.0); };
  auto &Z2 = Z.TwoBody;
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();

  auto pq_cc = [](TwoBodyChannel_CC &tbc, int idx, int nK) -> std::array<int, 2> {
    Ket &ket = tbc.GetKet(idx % nK);
    if (idx < nK)
      return {(int)ket.p, (int)ket.q};
    if (ket.p == ket.q)
      return {-1, -1};
    return {(int)ket.q, (int)ket.p};
  };
  auto idx2n = [](TwoBodyChannel_CC &tbc, int p, int q, int nK) -> int {
    int loc = tbc.GetLocalIndex(std::min(p, q), std::max(p, q));
    if (loc < 0)
      return -1;
    if (p > q)
      loc += nK;
    return loc;
  };

  // IMSRG DoTensorPandya (adcb) — hats × NineJ, no extra scale. λ=0 ⇒ Jbra=Jket.
  auto pandya_eta = [&](int a, int b, int c, int d, int Jbra,
                        int Jket) -> double {
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &ob = Z.modelspace->GetOrbit(b);
    Orbit &oc = Z.modelspace->GetOrbit(c);
    Orbit &od = Z.modelspace->GetOrbit(d);
    const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
    const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
    if (not AngMom::Triangle(Jbra, Jket, lambda))
      return 0.0;
    double Xbar = 0.0;
    const int j1min = std::abs(oa.j2 - od.j2) / 2;
    const int j1max = (oa.j2 + od.j2) / 2;
    for (int J1 = j1min; J1 <= j1max; ++J1) {
      const int j2min =
          std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
      const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
      for (int J2 = j2min; J2 <= j2max; ++J2) {
        const double ninej = Z.modelspace->GetNineJ(
            ja, jd, J1, jb, jc, J2, Jbra, Jket, lambda);
        if (std::abs(ninej) < 1e-14)
          continue;
        const double hats = hatJ(J1) * hatJ(J2) * hatJ(Jbra) * hatJ(Jket);
        const double tbme = Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
        Xbar -= hats * Z.modelspace->phase((ob.j2 + od.j2) / 2 + Jket + J2) *
                ninej * tbme;
      }
    }
    return Xbar;
  };

  // Scalar Pandya Γ per CC channel (Factorized 2×nKets layout).
  std::deque<arma::mat> bar_Gamma(n_nonzero);
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

  // Rectangular CC pairs. Triangle(Jb,Jk,λ) — λ=0 ⇒ Jb=Jk.
  std::vector<std::array<int, 2>> cc_pairs;
  for (int ch_b = 0; ch_b < n_nonzero; ++ch_b) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    if (tb.GetNumberKets() < 1)
      continue;
    for (int ch_k = 0; ch_k < n_nonzero; ++ch_k) {
      TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
      if (tk.GetNumberKets() < 1)
        continue;
      if ((tb.parity + tk.parity) % 2 != 0)
        continue;
      if (tb.Tz != tk.Tz)
        continue;
      if (not AngMom::Triangle(tb.J, tk.J, lambda))
        continue;
      cc_pairs.push_back({ch_b, ch_k});
    }
  }

  // Ω̄^{Jb Jk}(ab;kj) and χ̄ = Γ̄·(occ_AbarBC ⊙ Ω̄). No hΩ, no (−1)^{J0+J1}.
  std::map<std::array<int, 2>, arma::mat> bar_Omega;
  std::map<std::array<int, 2>, arma::mat> bar_CHI_V;
  for (auto pk : cc_pairs) {
    const int ch_b = pk[0], ch_k = pk[1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    const int Jb = tb.J, Jk = tk.J;
    arma::mat Om(2 * nb, 2 * nk, arma::fill::zeros);
    arma::mat Occ(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto ab = pq_cc(tb, ibra, nb);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double n_a = oa.occ, nbar_a = 1.0 - n_a;
      const double n_b = ob.occ, nbar_b = 1.0 - n_b;
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto kj = pq_cc(tk, iket, nk);
        if (kj[0] < 0)
          continue;
        Orbit &ok = Z.modelspace->GetOrbit(kj[0]);
        const double n_k = ok.occ, nbar_k = 1.0 - n_k;
        const double occ_AbarBC = nbar_a * n_b * n_k + n_a * nbar_b * nbar_k;
        const double bo = pandya_eta(ab[0], ab[1], kj[0], kj[1], Jb, Jk);
        Om(ibra, iket) = bo;
        Occ(ibra, iket) = bo * occ_AbarBC;
      }
    }
    bar_Omega[pk] = std::move(Om);
    bar_CHI_V[pk] = bar_Gamma[ch_b] * Occ;
  }

  auto lookup_chi = [&](int a, int d, int b, int c, int J3, int J4) -> double {
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &od = Z.modelspace->GetOrbit(d);
    Orbit &ob = Z.modelspace->GetOrbit(b);
    Orbit &oc = Z.modelspace->GetOrbit(c);
    const int ch_ad = Z.modelspace->GetTwoBodyChannelIndex(
        J3, (oa.l + od.l) % 2, std::abs(oa.tz2 - od.tz2) / 2);
    const int ch_bc = Z.modelspace->GetTwoBodyChannelIndex(
        J4, (ob.l + oc.l) % 2, std::abs(ob.tz2 - oc.tz2) / 2);
    if (ch_ad < 0 or ch_bc < 0 or ch_ad >= n_nonzero or ch_bc >= n_nonzero)
      return 0.0;
    auto it = bar_CHI_V.find({ch_ad, ch_bc});
    if (it == bar_CHI_V.end())
      return 0.0;
    TwoBodyChannel_CC &tad = Z.modelspace->GetTwoBodyChannel_CC(ch_ad);
    TwoBodyChannel_CC &tbc = Z.modelspace->GetTwoBodyChannel_CC(ch_bc);
    const int nad = tad.GetNumberKets(), nbc = tbc.GetNumberKets();
    const int iad = idx2n(tad, a, d, nad);
    const int ibc = idx2n(tbc, b, c, nbc);
    if (iad < 0 or ibc < 0)
      return 0.0;
    if (iad >= (int)it->second.n_rows or ibc >= (int)it->second.n_cols)
      return 0.0;
    return it->second(iad, ibc);
  };

  // RC[χ]: NineJ middle-row (jd,jc) + Φ so λ→0 is Factorized (2J'+1)6j.
  // Pack χ̄_adbc − hZ χ̄_bcad (χ is not AS — two fold terms).
  std::map<std::array<int, 2>, arma::mat> bar_CHI_V_RC;
  for (auto pk : cc_pairs) {
    const int ch_b = pk[0], ch_k = pk[1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    const int J0 = tb.J, J1 = tk.J;
    arma::mat RC(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nb; ++ibra) {
      auto ab = pq_cc(tb, ibra, nb);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
      for (int iket = 0; iket < 2 * nk; ++iket) {
        auto cd = pq_cc(tk, iket, nk);
        if (cd[0] < 0)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(cd[0]);
        Orbit &od = Z.modelspace->GetOrbit(cd[1]);
        const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
        const int a = ab[0], b = ab[1], c = cd[0], d = cd[1];
        double sm = 0.0;
        const int j3min = std::abs(oa.j2 - od.j2) / 2;
        const int j3max = (oa.j2 + od.j2) / 2;
        for (int J3 = j3min; J3 <= j3max; ++J3) {
          const int j4min =
              std::max(std::abs(ob.j2 - oc.j2) / 2, std::abs(J3 - lambda));
          const int j4max = std::min((ob.j2 + oc.j2) / 2, J3 + lambda);
          for (int J4 = j4min; J4 <= j4max; ++J4) {
            const double pack =
                lookup_chi(a, d, b, c, J3, J4) -
                hZ * lookup_chi(b, c, a, d, J4, J3);
            if (std::abs(pack) < 1e-16)
              continue;
            const double ninej = Z.modelspace->GetNineJ(ja, jb, J0, jd, jc, J1,
                                                         J3, J4, lambda);
            if (std::abs(ninej) < 1e-16)
              continue;
            const double phi =
                Z.modelspace->phase((ob.j2 + od.j2) / 2 + J0 + J3);
            const double w = phi * hatJ(J0) * hatJ(J1) * hatJ(J3) * hatJ(J4) *
                             (hatJ(J3) / hatJ(J0)) * ninej;
            sm += Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J4) * w * pack;
          }
        }
        RC(ibra, iket) = sm;
      }
    }
    bar_CHI_V_RC[pk] = std::move(RC);
  }
  bar_CHI_V.clear();

  // W^{Jb Jk}(jl;ki) = Σ_{ab,Jm} Ω̄^{Jb Jm}(jl;ab) RC^{Jm Jk}(ab;ki)
  std::map<std::array<int, 2>, arma::mat> CHI_V_final;
  for (auto pk : cc_pairs) {
    const int ch_b = pk[0], ch_k = pk[1];
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(ch_b);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(ch_k);
    const int nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    arma::mat W(2 * nb, 2 * nk, arma::fill::zeros);
    for (int ch_m = 0; ch_m < n_nonzero; ++ch_m) {
      TwoBodyChannel_CC &tm = Z.modelspace->GetTwoBodyChannel_CC(ch_m);
      if (tm.GetNumberKets() < 1)
        continue;
      if (not AngMom::Triangle(tb.J, tm.J, lambda) or
          not AngMom::Triangle(tm.J, tk.J, lambda))
        continue;
      if ((tb.parity + tm.parity) % 2 != 0 or (tm.parity + tk.parity) % 2 != 0)
        continue;
      if (tb.Tz != tm.Tz or tm.Tz != tk.Tz)
        continue;
      auto itO = bar_Omega.find({ch_b, ch_m});
      auto itR = bar_CHI_V_RC.find({ch_m, ch_k});
      if (itO == bar_Omega.end() or itR == bar_CHI_V_RC.end())
        continue;
      W += itO->second * itR->second;
    }
    CHI_V_final[pk] = std::move(W);
  }
  bar_Omega.clear();
  bar_CHI_V_RC.clear();

  auto lookup_W = [&](int p, int q, int r, int s, int J3, int J4) -> double {
    Orbit &op = Z.modelspace->GetOrbit(p);
    Orbit &oq = Z.modelspace->GetOrbit(q);
    Orbit &or_ = Z.modelspace->GetOrbit(r);
    Orbit &os = Z.modelspace->GetOrbit(s);
    const int ch_pq = Z.modelspace->GetTwoBodyChannelIndex(
        J3, (op.l + oq.l) % 2, std::abs(op.tz2 - oq.tz2) / 2);
    const int ch_rs = Z.modelspace->GetTwoBodyChannelIndex(
        J4, (or_.l + os.l) % 2, std::abs(or_.tz2 - os.tz2) / 2);
    if (ch_pq < 0 or ch_rs < 0 or ch_pq >= n_nonzero or ch_rs >= n_nonzero)
      return 0.0;
    auto it = CHI_V_final.find({ch_pq, ch_rs});
    if (it == CHI_V_final.end())
      return 0.0;
    TwoBodyChannel_CC &tp = Z.modelspace->GetTwoBodyChannel_CC(ch_pq);
    TwoBodyChannel_CC &tr = Z.modelspace->GetTwoBodyChannel_CC(ch_rs);
    const int np = tp.GetNumberKets(), nr = tr.GetNumberKets();
    const int ip = idx2n(tp, p, q, np);
    const int ir = idx2n(tr, r, s, nr);
    if (ip < 0 or ir < 0)
      return 0.0;
    if (ip >= (int)it->second.n_rows or ir >= (int)it->second.n_cols)
      return 0.0;
    return it->second(ip, ir);
  };

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

  // Inv: NineJ continuous of Factorized L1884–2013. (1−P)^2 only here.
#pragma omp parallel for schedule(dynamic, 1)
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

    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = bra.p, j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ji = oi.j2, jj = oj.j2;
      const int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);
      const int ketmin = (ch_bra == ch_ket) ? ibra : 0;
      for (int iket = ketmin; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = ket.p, l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int jk = ok.j2, jl = ol.j2;
        double commijkl = 0, commjikl = 0, commijlk = 0, commjilk = 0;

        // Sector 1: W(jl;ik) / W(ik;jl)
        const int j3min = std::abs(jj - jl) / 2;
        const int j3max = (jj + jl) / 2;
        for (int J3 = j3min; J3 <= j3max; ++J3) {
          const int j4min =
              std::max(std::abs(ji - jk) / 2, std::abs(J3 - lambda));
          const int j4max = std::min((ji + jk) / 2, J3 + lambda);
          for (int J4 = j4min; J4 <= j4max; ++J4) {
            const double ninej = Z.modelspace->GetNineJ(
                jj * 0.5, ji * 0.5, J0, jl * 0.5, jk * 0.5, J0, J3, J4, lambda);
            if (std::abs(ninej) < 1e-16)
              continue;
            const double phi =
                Z.modelspace->phase((ji + jl) / 2 + J0 + J3);
            const double w = phi * hatJ(J0) * hatJ(J0) * hatJ(J3) * hatJ(J4) *
                             (hatJ(J3) / hatJ(J0)) * ninej;
            const double me1 = lookup_W(j, l, i, k, J3, J4);
            const double me2 = lookup_W(i, k, j, l, J4, J3);
            commjikl -= Z.modelspace->phase(J3 + (ji + jk) / 2) * w * me1;
            commijlk -= Z.modelspace->phase(J3 + (jj + jl) / 2) * w * me2;
          }
        }

        // Sector 2: W(jk;il) / W(il;jk)
        const int j3min2 = std::abs(jj - jk) / 2;
        const int j3max2 = (jj + jk) / 2;
        for (int J3 = j3min2; J3 <= j3max2; ++J3) {
          const int j4min =
              std::max(std::abs(ji - jl) / 2, std::abs(J3 - lambda));
          const int j4max = std::min((ji + jl) / 2, J3 + lambda);
          for (int J4 = j4min; J4 <= j4max; ++J4) {
            const double ninej = Z.modelspace->GetNineJ(
                ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, J0, J3, J4, lambda);
            if (std::abs(ninej) < 1e-16)
              continue;
            const double phi =
                Z.modelspace->phase((jj + jl) / 2 + J0 + J3);
            const double w = phi * hatJ(J0) * hatJ(J0) * hatJ(J3) * hatJ(J4) *
                             (hatJ(J3) / hatJ(J0)) * ninej;
            const double me1 = lookup_W(j, k, i, l, J3, J4);
            const double me2 = lookup_W(i, l, j, k, J4, J3);
            commjilk -= Z.modelspace->phase(J3 + (ji + jl) / 2) * w * me1;
            commijkl -= Z.modelspace->phase(J3 + (jj + jk) / 2) * w * me2;
          }
        }

        double zijkl =
            (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
      }
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_c / chi^lambda — Pandya→DGEMM→inv ring (≡ m ≡ AMC direct).
/// Rank: T×S / S×T → tensor χ^λ.

////////////////////////////////////////////////////////////////////////////
/// GIVc Path B — χ^λ → Pandya → mid-J DGEMM → inv → fermionic AS.
/// λ≠0 only (scalar λ=0 Factorized CHI_VII is separate — do not retune).
///
/// Gold: m ≡ AMC direct (tts_ring) ≡ this Path B (drop AMC-sample inv minus).
/// IMSRG DoTensorPandya(adcb) with index map:
///   AMC barχ(p,b,a,r) = IMSRG(p,r,a,b);  AMC barΩ(a,q,s,b) = IMSRG(a,b,s,q).
/// Then Z = ½(1−P)(1−P) X with J-phases. Store Z_unred = Z_red/Ĵ.
////////////////////////////////////////////////////////////////////////////
static void comm223_232_GIVc_pathB(const Operator &Eta, const Operator &Gamma,
                                   Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  const int lambda = Eta.GetJRank();
  if (lambda == 0) {
    // Tensor Path B only (λ≠0). λ=0: no-op here.
    Z.profiler.timer["comm223_232_GIVc_pathB"] += omp_get_wtime() - t_start;
    return;
  }
  const double hat_lambda =
      std::sqrt(2.0 * std::max(lambda, 0) + 1.0);
  const double hat_lambda_inv = 1.0 / hat_lambda;
  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };

  // Unreduced Γ for AMC chi_lambda (no 6j)
  const Operator *Gp = &Gamma;
  Operator Gunred;
  if (Gamma.IsReduced() and Gamma.GetJRank() == 0) {
    Gunred = Gamma;
    Gunred.MakeNotReduced();
    Gp = &Gunred;
  }

  int max_j2 = 0;
  for (auto x : Z.modelspace->all_orbits)
    max_j2 = std::max(max_j2, Z.modelspace->GetOrbit(x).j2);
  const int max_J = max_j2;

  std::vector<index_t> allorb(Z.modelspace->all_orbits.begin(),
                              Z.modelspace->all_orbits.end());
  int n_orb_alloc = 0;
  for (auto o : allorb)
    n_orb_alloc = std::max(n_orb_alloc, (int)o + 1);
  const int n_orb = (int)allorb.size();

  // χ^λ table: χ(i,j,k,l; J0,J1) in GetTBME_J units (AMC chi_lambda.tex)
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

  ChiTabJJ Chi;
  Chi.allocate(n_orb_alloc, max_J);

#pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < n_orb; ++ii) {
    index_t i = allorb[ii];
    Orbit &oi = Z.modelspace->GetOrbit(i);
    const double ji = oi.j2 / 2.0;
    for (auto j : allorb) {
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const double jj = oj.j2 / 2.0;
      const double n_j = oj.occ, nbar_j = 1.0 - n_j;
      for (auto k : allorb) {
        Orbit &ok = Z.modelspace->GetOrbit(k);
        const double jk = ok.j2 / 2.0;
        for (auto l : allorb) {
          Orbit &ol = Z.modelspace->GetOrbit(l);
          const double jl = ol.j2 / 2.0;
          const double n_l = ol.occ, nbar_l = 1.0 - n_l;
          for (int J0 = 0; J0 <= max_J; ++J0) {
            if (not AngMom::Triangle(ji, jj, (double)J0))
              continue;
            for (int J1 = 0; J1 <= max_J; ++J1) {
              if (not AngMom::Triangle(jk, jl, (double)J1))
                continue;
              if (not AngMom::Triangle(J0, J1, lambda))
                continue;
              double sm = 0.0;
              for (auto a : allorb) {
                Orbit &oa = Z.modelspace->GetOrbit(a);
                const double ja = oa.j2 / 2.0;
                const double n_a = oa.occ, nbar_a = 1.0 - n_a;
                for (auto b : allorb) {
                  Orbit &ob = Z.modelspace->GetOrbit(b);
                  const double jb = ob.j2 / 2.0;
                  const double n_b = ob.occ, nbar_b = 1.0 - n_b;
                  const double w_l =
                      nbar_a * nbar_b * n_l + n_a * n_b * nbar_l;
                  const double w_j =
                      nbar_a * nbar_b * n_j + n_a * n_b * nbar_j;
                  if (std::abs(w_l) > 1e-12 and
                      AngMom::Triangle(ja, jb, (double)J0)) {
                    sm += w_l * Gp->TwoBody.GetTBME_J(J0, J0, i, j, a, b) *
                          Eta.TwoBody.GetTBME_J(J0, J1, a, b, k, l);
                  }
                  if (std::abs(w_j) > 1e-12 and
                      AngMom::Triangle(ja, jb, (double)J1)) {
                    sm += w_j * Eta.TwoBody.GetTBME_J(J0, J1, i, j, a, b) *
                          Gp->TwoBody.GetTBME_J(J1, J1, a, b, k, l);
                  }
                }
              }
              Chi.at(i, j, k, l, J0, J1) = sm;
            }
          }
        }
      }
    }
  }

  // ---- IMSRG tensor Pandya (adcb). Chi from table; Ω from Eta. ----
  auto pandya_tensor = [&](bool use_chi, int a, int b, int c, int d, int Jbra,
                           int Jket) -> double {
    if (not AngMom::Triangle(Jbra, Jket, lambda))
      return 0.0;
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &ob = Z.modelspace->GetOrbit(b);
    Orbit &oc = Z.modelspace->GetOrbit(c);
    Orbit &od = Z.modelspace->GetOrbit(d);
    const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5;
    const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
    double Xbar = 0.0;
    const int j1min = std::abs(oa.j2 - od.j2) / 2;
    const int j1max = (oa.j2 + od.j2) / 2;
    for (int J1 = j1min; J1 <= j1max; ++J1) {
      const int j2min =
          std::max(std::abs(oc.j2 - ob.j2) / 2, std::abs(J1 - lambda));
      const int j2max = std::min((oc.j2 + ob.j2) / 2, J1 + lambda);
      for (int J2 = j2min; J2 <= j2max; ++J2) {
        const double ninej = AngMom::NineJ(ja, jd, (double)J1, jb, jc,
                                           (double)J2, (double)Jbra,
                                           (double)Jket, (double)lambda);
        if (std::abs(ninej) < 1e-10)
          continue;
        const double hats = std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) *
                                      (2.0 * Jbra + 1.0) * (2.0 * Jket + 1.0));
        const double me =
            use_chi ? Chi(a, d, c, b, J1, J2)
                    : Eta.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
        Xbar -= hats * AngMom::phase((ob.j2 + od.j2) / 2 + Jket + J2) * ninej *
                me;
      }
    }
    return Xbar;
  };

  // ---- Mid-J DGEMM: barProd[J]((p,r),(s,q)) = Σ_Jmid pref barχ barΩ ----
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  std::deque<arma::mat> barProd(n_cc);

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    const int nKets_cc = tbc_cc.GetNumberKets();
    if (nKets_cc < 1)
      continue;
    const int J_cc = tbc_cc.J;
    const int n2 = nKets_cc * 2;
    const int parity_cc = tbc_cc.parity;
    const int Tz_cc = tbc_cc.Tz;
    barProd[ch_cc] = arma::mat(n2, n2, arma::fill::zeros);

    for (int Jmid = 0; Jmid <= max_J; ++Jmid) {
      if (not AngMom::Triangle(J_cc, Jmid, lambda))
        continue;
      const int ch_mid =
          Z.modelspace->GetTwoBodyChannelIndex(Jmid, parity_cc, Tz_cc);
      if (ch_mid < 0 or ch_mid >= n_cc)
        continue;
      TwoBodyChannel_CC &tbc_mid = Z.modelspace->GetTwoBodyChannel_CC(ch_mid);
      const int nKm = tbc_mid.GetNumberKets();
      if (nKm < 1)
        continue;
      const int n2m = nKm * 2;
      arma::mat Chi_L(n2, n2m, arma::fill::zeros);
      arma::mat Om_R(n2m, n2, arma::fill::zeros);

      for (int ibra = 0; ibra < n2; ++ibra) {
        int p, r;
        if (ibra < nKets_cc) {
          Ket &bra = tbc_cc.GetKet(ibra);
          p = bra.p;
          r = bra.q;
        } else {
          Ket &bra = tbc_cc.GetKet(ibra - nKets_cc);
          r = bra.p;
          p = bra.q;
        }
        if (ibra >= nKets_cc and p == r)
          continue;
        Orbit &op = Z.modelspace->GetOrbit(p);
        Orbit &orr = Z.modelspace->GetOrbit(r);
        for (int iket = 0; iket < n2m; ++iket) {
          int a, b;
          if (iket < nKm) {
            Ket &ket = tbc_mid.GetKet(iket);
            a = ket.p;
            b = ket.q;
          } else {
            Ket &ket = tbc_mid.GetKet(iket - nKm);
            b = ket.p;
            a = ket.q;
          }
          if (iket >= nKm and a == b)
            continue;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          if (std::abs(op.tz2 + ob.tz2 - orr.tz2 - oa.tz2) != Eta.GetTRank())
            continue;
          // AMC barχ(p,b,a,r) = IMSRG(p,r,a,b)
          Chi_L(ibra, iket) = pandya_tensor(true, p, r, a, b, J_cc, Jmid);
        }
      }

      for (int ibra = 0; ibra < n2m; ++ibra) {
        int a, b;
        if (ibra < nKm) {
          Ket &bra = tbc_mid.GetKet(ibra);
          a = bra.p;
          b = bra.q;
        } else {
          Ket &bra = tbc_mid.GetKet(ibra - nKm);
          b = bra.p;
          a = bra.q;
        }
        if (ibra >= nKm and a == b)
          continue;
        Orbit &oa = Z.modelspace->GetOrbit(a);
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (int iket = 0; iket < n2; ++iket) {
          int s, q;
          if (iket < nKets_cc) {
            Ket &ket = tbc_cc.GetKet(iket);
            s = ket.p;
            q = ket.q;
          } else {
            Ket &ket = tbc_cc.GetKet(iket - nKets_cc);
            q = ket.p;
            s = ket.q;
          }
          if (iket >= nKets_cc and s == q)
            continue;
          Orbit &os = Z.modelspace->GetOrbit(s);
          Orbit &oq = Z.modelspace->GetOrbit(q);
          if (std::abs(oa.tz2 + oq.tz2 - ob.tz2 - os.tz2) != Eta.GetTRank())
            continue;
          // AMC barΩ(a,q,s,b) = IMSRG(a,b,s,q)
          Om_R(ibra, iket) = pandya_tensor(false, a, b, s, q, Jmid, J_cc);
        }
      }

      const double pref =
          hat_lambda_inv * AngMom::phase(Jmid + lambda);
      barProd[ch_cc] += pref * Chi_L * Om_R;
    }
  }

  // Bare X_pqsr from corrected inv Pandya (no sample minus):
  //   X = +Ĵ0 Σ_Jp Ĵp {jr js J0; jq jp Jp} · [(-1)^{Jp}/Ĵp barProd_Jp]
  auto ring_X = [&](index_t p, index_t q, index_t s, index_t r,
                    int J0) -> double {
    Orbit &op = Z.modelspace->GetOrbit(p);
    Orbit &oq = Z.modelspace->GetOrbit(q);
    Orbit &os = Z.modelspace->GetOrbit(s);
    Orbit &orr = Z.modelspace->GetOrbit(r);
    const double jp = op.j2 / 2.0, jq = oq.j2 / 2.0;
    const double js = os.j2 / 2.0, jr = orr.j2 / 2.0;
    if (not AngMom::Triangle(jp, jq, (double)J0) or
        not AngMom::Triangle(js, jr, (double)J0))
      return 0.0;

    double sm = 0.0;
    const int parity_cc = (op.l + orr.l) % 2;
    const int Tz_cc = std::abs(op.tz2 - orr.tz2) / 2;
    const int Jpmin = std::max(std::abs(op.j2 - orr.j2), std::abs(oq.j2 - os.j2)) / 2;
    const int Jpmax = std::min(op.j2 + orr.j2, oq.j2 + os.j2) / 2;
    for (int Jp = Jpmin; Jp <= Jpmax; ++Jp) {
      const double six =
          AngMom::SixJ(jr, js, (double)J0, jq, jp, (double)Jp);
      if (std::abs(six) < 1e-16)
        continue;
      const int ch_cc =
          Z.modelspace->GetTwoBodyChannelIndex(Jp, parity_cc, Tz_cc);
      if (ch_cc < 0 or ch_cc >= n_cc)
        continue;
      if (barProd[ch_cc].n_rows == 0)
        continue;
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      const int nk = tbc_cc.GetNumberKets();
      int indx_pr = tbc_cc.GetLocalIndex(std::min(p, r), std::max(p, r));
      int indx_sq = tbc_cc.GetLocalIndex(std::min(s, q), std::max(s, q));
      if (indx_pr < 0 or indx_sq < 0)
        continue;
      indx_pr += (p > r ? nk : 0);
      indx_sq += (s > q ? nk : 0);
      // Also try ket as (q,s) ordering used in Om_R fill: ket (s,q) with
      // exchange block when q>s stored as swapped.
      // Om_R iket: if iket < nK: (s,q)=(ket.p,ket.q); else (s,q)=(ket.q,ket.p)
      // so (s,q) with s>q uses exchange block: local(min,max)+nk when first>second
      // For (s,q): if s<=q use local(s,q); if s>q use local(q,s)+nk — but we
      // used min(s,q),max(s,q) then += (s>q ? nk : 0). When s>q, min=q,max=s,
      // +nk → exchange of (q,s) which is (s,q) with s>q. Good.
      // For (p,r): same.
      if (indx_pr >= (int)barProd[ch_cc].n_rows or
          indx_sq >= (int)barProd[ch_cc].n_cols)
        continue;
      const double bar = barProd[ch_cc](indx_pr, indx_sq);
      const double mid = AngMom::phase(Jp) / hat(Jp) * bar;
      sm += hat(Jp) * six * mid;
    }
    return hat(J0) * sm; // corrected Path B (drop AMC-sample overall minus)
  };

  auto fold_red = [&](index_t i, index_t j, index_t k, index_t l,
                      int J0) -> double {
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ol = Z.modelspace->GetOrbit(l);
    const double x = ring_X(i, j, k, l, J0);
    const double xkl = ring_X(i, j, l, k, J0);
    const double xij = ring_X(j, i, k, l, J0);
    const double xijkl = ring_X(j, i, l, k, J0);
    const double pkl =
        AngMom::phase((ok.j2 + ol.j2) / 2) * AngMom::phase(J0);
    const double pij =
        AngMom::phase((oi.j2 + oj.j2) / 2) * AngMom::phase(J0);
    return 0.5 * (x - pkl * xkl - pij * xij + pij * pkl * xijkl);
  };

  auto &Z2 = Z.TwoBody;
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

        const double z_red = fold_red(p, g, q, h, J0);
        double z = z_red / hat_J0;
        if (p == g)
          z /= PhysConst::SQRT2;
        if (q == h)
          z /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, z);
      }
    }
  }

  Z.profiler.timer["comm223_232_GIVc_pathB"] += omp_get_wtime() - t_start;
}


////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVc(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  // Γ^{IV_c} / χ^λ ethS entry: Pandya → mid-J DGEMM → inv → AS.
  // Locked gold: m ≡ AMC direct (tts_ring) ≡ this Path B.
  comm223_232_GIVc_pathB(Eta, Gamma, Z);
}


void comm223_132(const Operator &Eta, const Operator &Gamma, Operator &Z) {

  double t_internal = omp_get_wtime(); // timer
  double t_start = omp_get_wtime();    // timer

  Z.modelspace->PreCalculateSixJ();
  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }

  // determine symmetry
  int hEta = Eta.IsHermitian() ? 1 : -1;
  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  // int hZ = Z.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  int nch = Gamma.modelspace->GetNumberTwoBodyChannels();
  {
// full matrix
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          double zijkl = 0;

          for(int iket2 = 0; iket2 < nkets*2; iket2++)  // cution for tensor here, only work for scaler
         {
          size_t c = 0;
          size_t a = 0;
          if(iket2 <  nkets){
             Ket &ket2 =tbc_ket.GetKet(iket2);
           c = ket2.p;
           a = ket2.q;}
          else{
             Ket &ket2 =tbc_ket.GetKet(iket2-nkets);
           c = ket2.q;
           a = ket2.p;
            if(a==c)continue;
          }
          Orbit& oa = Z.modelspace->GetOrbit(a);
          double na=oa.occ;
    for (auto &b : Z.GetOneBodyChannel(oa.l, oa.j2, oa.tz2)) // delta_jd je
          {
            Orbit& ob = Z.modelspace->GetOrbit(b);
            double nb=ob.occ;

            double xijcb = Gamma.TwoBody.GetTBME_J(J, J, i,j,c,b);
            double yijcb = Eta.TwoBody.GetTBME_J(J, J, i,j,c,b);
            double xcakl = Gamma.TwoBody.GetTBME_J(J, J, c,a,k,l);
            double ycakl = Eta.TwoBody.GetTBME_J(J, J, c,a,k,l);

            zijkl += xijcb*ycakl*Eta.OneBody(b,a)*((1.-na)*nb-na*(1.-nb));
            zijkl -= yijcb*xcakl*Eta.OneBody(b,a)*((1.-na)*nb-na*(1.-nb));

          } // b
          }
          //if(std::abs(zijkl)>0.000001)std::cout<< zijkl<< " here"<<std::endl;
          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

   Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch
  }


  comm223_132_cross(Eta, Gamma, Z);


  {
  // now we take care of eta contract to either eta or gamma.

 arma::mat CHI_I = Gamma.OneBody * 0;
  arma::mat CHI_II = Gamma.OneBody * 0;

   std::vector<index_t> p_list, q_list;
  for (auto p : Z.modelspace->all_orbits) {
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
      p_list.push_back(p);
      q_list.push_back(q);
    }
  }
  size_t ipq_max = p_list.size();
/// Build the intermediate one-body operators
  for (size_t ipq = 0; ipq < ipq_max; ipq++) {
    index_t p = p_list[ipq];
    index_t q = q_list[ipq];
    Orbit &op = Z.modelspace->GetOrbit(p);
    Orbit &oq = Z.modelspace->GetOrbit(q);

    double chi_pq = 0;
    double chiY_pq = 0;


      for (auto i : Z.modelspace->all_orbits) {
        Orbit &oi = Z.modelspace->GetOrbit(i);
        double n_i = oi.occ;

        for (auto j : Z.OneBodyChannels.at({oi.l, oi.j2, oi.tz2})){

            Orbit &oj = Z.modelspace->GetOrbit(j);

          double n_j = oj.occ;

          double occfactor = (1.-n_i)*n_j-n_i*(1.0 - n_j);
          if (abs(occfactor) < 1.e-7) continue;

          int J2min =
              std::max(std::abs(oi.j2 - op.j2), std::abs(oj.j2 - oq.j2)) / 2;
          int J2max = std::min(oi.j2 + op.j2, oj.j2 + oq.j2) / 2;

          for (int J2 = J2min; J2 <= J2max; J2++) {
            double xpiqj = Eta.TwoBody.GetTBME_J(J2,  p, i, q, j);
            double ypiqj = Gamma.TwoBody.GetTBME_J(J2,  p, i, q, j);
            chi_pq +=
               occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xpiqj*Eta.OneBody(j, i);
            chiY_pq +=
               occfactor * (2 * J2 + 1) / (oq.j2 + 1) * ypiqj*Eta.OneBody(j, i);
         //   std::cout << p<<" "<<i<< " "<<q<<" "<<j<<" "<< xpiqj<<" "<<Eta.OneBody(j, i)<<occfactor<<std::endl;
            //if(std::abs(xpiqj*Eta.OneBody(j, i))>0.000001)std::cout << Eta.TwoBody.GetTBME_J(J2,  p, i, q, j)<<" "<<Eta.OneBody(j, i)<<" "<<std::endl;
            //if(std::abs(xpiqj*Eta.OneBody(j, i))>0.000001)std::cout << Eta.TwoBody.GetTBME_J(J2,  q, j, p, i)<<" "<<Eta.OneBody(i, j)<<std::endl;
          } // for j2
        } // for j
     } // for i
    //
      //std::cout<<"pq="<< p<<" " <<q<<" " << chi_pq<<std::endl;
     CHI_I(p, q) = chi_pq;
    CHI_II(p, q) = chiY_pq;
    } // for pq






    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          double zijkl = 0;
          for (size_t a : Z.modelspace->all_orbits)
          {
            //           Orbit& oa = Z.modelspace->GetOrbit(a);
            double xajkl = Gamma.TwoBody.GetTBME_J(J, J, a, j, k, l);
            double xiakl = Gamma.TwoBody.GetTBME_J(J, J, i, a, k, l);
            double xijal = Gamma.TwoBody.GetTBME_J(J, J, i, j, a, l);
            double xijka = Gamma.TwoBody.GetTBME_J(J, J, i, j, k, a);

            double yajkl =   Eta.TwoBody.GetTBME_J(J, J, a, j, k, l);
            double yiakl =   Eta.TwoBody.GetTBME_J(J, J, i, a, k, l);
            double yijal =   Eta.TwoBody.GetTBME_J(J, J, i, j, a, l);
            double yijka =   Eta.TwoBody.GetTBME_J(J, J, i, j, k, a);

            zijkl += CHI_I(i, a) * xajkl + CHI_I(j, a) * xiakl - xijal * CHI_I(a, k) - xijka * CHI_I(a, l);
            zijkl -= CHI_II(i, a) * yajkl + CHI_II(j, a) * yiakl - yijal * CHI_II(a, k) - yijka * CHI_II(a, l);
          } // a
          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;
          //if(abs(zijkl)>0.000001)std::cout<< i<< " "<<j<<" "<<k<<" "<<l<<" "<<zijkl<<std::endl;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch



  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  CHI_I.clear();
  CHI_II.clear();
  }

  if (Commutator::verbose) {
    Z.profiler.timer["_" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

}

void comm223_132_cross(const Operator &Eta, const Operator &Gamma, Operator &Z) {
  double t_start = omp_get_wtime();

  Z.modelspace->PreCalculateSixJ();
  auto &Z2 = Z.TwoBody;

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }

  int nch = Gamma.modelspace->GetNumberTwoBodyChannels();
  int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();

  std::deque<arma::mat> IntermediateTwobody(n_nonzero);
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    int J_cc = tbc_cc.J;

    arma::mat Eta_bar =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    arma::mat Eta_bar_nnnn =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    arma::mat Gamma_bar =
        arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);

    for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc) {
      int a = 0;
      int b = 0;
      if (ibra_cc < nKets_cc) {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
        a = bra_cc.p;
        b = bra_cc.q;
      } else {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
        a = bra_cc.q;
        b = bra_cc.p;
        if (a == b)
          continue;
      }

      Orbit &oa = Z.modelspace->GetOrbit(a);
      double ja = oa.j2 * 0.5;
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double jb = ob.j2 * 0.5;

      for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
        int c = 0;
        int d = 0;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
          c = ket_cc_cd.p;
          d = ket_cc_cd.q;
        } else {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc_cd.p;
          c = ket_cc_cd.q;
          if (c == d)
            continue;
        }

        Orbit &oc = Z.modelspace->GetOrbit(c);
        double jc = oc.j2 * 0.5;
        Orbit &od = Z.modelspace->GetOrbit(d);
        double jd = od.j2 * 0.5;

        int jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
        int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;

        double Xbar = 0;
        double Ybar = 0;

        int dJ_std = 1;
        if ((a == d or b == c)) {
          dJ_std = 2;
          jmin += jmin % 2;
        }
        for (int J_std = jmin; J_std <= jmax; J_std += dJ_std) {
          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
          if (std::abs(sixj1) > 1e-8) {
            int phase = Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std);
            Xbar += (2 * J_std + 1) * sixj1 *
                    Eta.TwoBody.GetTBME_J(J_std, a, d, b, c) * phase;
            Ybar += (2 * J_std + 1) * sixj1 *
                    Gamma.TwoBody.GetTBME_J(J_std, a, d, b, c) * phase;
          }
        }

        Gamma_bar(ibra_cc, iket_cc) = Ybar;
        Eta_bar(ibra_cc, iket_cc) = Xbar;
      }
    }

    for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc) {
      int a = 0;
      int b = 0;
      if (ibra_cc < nKets_cc) {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
        a = bra_cc.p;
        b = bra_cc.q;
      } else {
        Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
        b = bra_cc.p;
        a = bra_cc.q;
        if (a == b)
          continue;
      }

      Orbit &oa = Z.modelspace->GetOrbit(a);
      double n_a = oa.occ;
      double nbar_a = 1 - n_a;
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double n_b = ob.occ;
      double nbar_b = 1 - n_b;

      for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc) {
        int c = 0;
        int d = 0;
        if (iket_cc < nKets_cc) {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
          c = ket_cc_cd.p;
          d = ket_cc_cd.q;
        } else {
          Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
          d = ket_cc_cd.p;
          c = ket_cc_cd.q;
          if (c == d)
            continue;
        }

        Orbit &oc = Z.modelspace->GetOrbit(c);
        double n_c = oc.occ;
        double nbar_c = 1 - n_c;
        Orbit &od = Z.modelspace->GetOrbit(d);
        double n_d = od.occ;
        double nbar_d = 1 - n_d;

        if (a == c) {
          Eta_bar_nnnn(ibra_cc, iket_cc) +=
              Eta.OneBody(d, b) * (nbar_b * n_d - n_b * nbar_d);
        }
        if (b == d) {
          Eta_bar_nnnn(ibra_cc, iket_cc) +=
              Eta.OneBody(a, c) * (nbar_a * n_c - n_a * nbar_c);
        }
      }
    }

    IntermediateTwobody[ch_cc] = Eta_bar * Eta_bar_nnnn * Gamma_bar;
  }

  for (int ch = 0; ch < nch; ++ch) {
    int ch_bra = ch_bra_list[ch];
    int ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();
    int J0 = tbc_bra.J;

    if (nbras == 0 or nkets == 0)
      continue;

    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      int ji = oi.j2;
      int jj = oj.j2;

      int ketmin = 0;
      if (ch_bra == ch_ket)
        ketmin = ibra;
      for (int iket = ketmin; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int jk = ok.j2;
        int jl = ol.j2;
        double commijkl = 0;
        double commjikl = 0;
        double commijlk = 0;
        double commjilk = 0;

        int parity_cc = (oi.l + ol.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        int Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0,
                                               jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
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
          double me1 = IntermediateTwobody[ch_cc](indx_jk, indx_li);
          commjilk -= (2 * Jprime + 1) * sixj1 * me1;

          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          me1 = IntermediateTwobody[ch_cc](indx_il, indx_kj);
          commijkl -= (2 * Jprime + 1) * sixj1 * me1;
        }

        parity_cc = (oi.l + ok.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        Jpmin = std::max(std::abs(int(jj - jl)), std::abs(int(jk - ji))) / 2;
        Jpmax = std::min(int(jj + jl), int(jk + ji)) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0,
                                               jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
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
          double me1 = IntermediateTwobody[ch_cc](indx_ik, indx_lj);
          commijlk -= (2 * Jprime + 1) * sixj1 * me1;

          indx_ki += (k > i ? nkets_cc : 0);
          indx_jl += (j > l ? nkets_cc : 0);
          me1 = IntermediateTwobody[ch_cc](indx_jl, indx_ki);
          commjikl -= (2 * Jprime + 1) * sixj1 * me1;
        }

        double zijkl =
            (commijkl - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikl);
        zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commijlk +
                  Z.modelspace->phase((jk + jl + ji + jj) / 2) * commjilk);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;

        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }

  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    IntermediateTwobody[ch_cc].clear();
  }
  IntermediateTwobody.clear();

  if (Commutator::verbose) {
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }
}

} // namespace FactorizedDoubleCommutator_eths

} // namespace Commutator
