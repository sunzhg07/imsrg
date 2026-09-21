
#include "FactorizedDoubleCommutator_eths.hh"
#include "Commutator.hh"
#include "ReferenceImplementations.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <omp.h>
#include <map>
#include <set>
#include <cstdint>
#include <array>
#include <deque>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
namespace Commutator {

namespace FactorizedDoubleCommutator_eths {

// ethS takes Ω as a reduced RME for every λ, including λ=0. The IMSRG scalar
// convention (RandomOp, EOM η at J=0, IMSRG generators) is unreduced; convert
// on entry so λ=0 reproduces the scalar FactorizedDoubleCommutator exactly
// (gold: scalar comm223ss→231ss+232ss+132ss, agreement ~1e-11).
// No copy is made when Ω is already reduced or has λ>0.
namespace {
struct ReducedEtaView {
  std::unique_ptr<Operator> copy;
  const Operator &ref;
  explicit ReducedEtaView(const Operator &Eta) : ref(Pick(Eta, copy)) {}

private:
  static const Operator &Pick(const Operator &Eta,
                              std::unique_ptr<Operator> &copy) {
    if (Eta.GetJRank() != 0 or Eta.IsReduced())
      return Eta;
    copy = std::make_unique<Operator>(Eta);
    copy->MakeReduced();
    return *copy;
  }
};
} // namespace

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

// FDC.cc comm223_232_chi2b ~1227–1393, tensor Ω version with the same operator
// names: one Pandya of Ω (bar_Eta; rank-λ → rectangular CC pair, 9j) and of the
// scalar Γ (bar_Gamma; 6j, cached), occupancy nnnbar_Eta / nnnbar_Eta_d (not a
// second transform), then χ̄ by DGEMM inside the same CC block.
//   bar_CHI_V     = bar_Gamma · nnnbar_Eta                    χ^ι  → GIVb
//   bar_CHI_VI_II = hΩ (−1)^{J+J'} nnnbar_Eta_dᵀ · bar_Gamma   χ^κ  → GIVa
// χ^η (GIIIa) and χ^λ (GIVc) stay diagram-local (AMC Path B seeds).
namespace {
using CcPair = std::array<int, 2>;
using CcMatMap = std::map<CcPair, arma::mat>;
struct PandyaBars232 {
  std::vector<CcPair> cc_pairs; // {ch_bra_cc, ch_ket_cc} coupled by Ω^{λ,T,π}
  CcMatMap bar_Eta;             // Ω̄^{J J'}(ab;cd), IMSRG adcb tensor Pandya
  CcMatMap bar_CHI_V;           // χ̄^ι(il;kj)  keyed {ch_il, ch_kj}
  CcMatMap bar_CHI_VI_II;       // χ̄^κ(il;kj)  keyed {ch_il, ch_kj}
  const std::deque<arma::mat> *bar_Gamma = nullptr; // Γ̄^J(ab;cd), square CC
};
void FillPandyaBars232(const Operator &Eta, const Operator &Gamma, Operator &Z,
                       PandyaBars232 &B);
const std::deque<arma::mat> &CachedBarGammaScalarCC(const Operator &Gamma,
                                                    Operator &Z, int hGamma);
void comm223_232_GIVa_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const PandyaBars232 &B);
void comm223_232_GIVb_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const PandyaBars232 &B);
// Shared 2n ordered-pair channel helpers (defined with TensorChannelPairs).
std::array<int, 2> PqTB(TwoBodyChannel &tbc, int idx, int nK);
int Index2n(TwoBodyChannel &tbc, int p, int q);
arma::mat Omega2n(const Operator &Op, ModelSpace *ms, int ch0, int ch1);
void OccWeights2n(ModelSpace *ms, int ch, arma::vec &w_pp, arma::vec &w_hh);
void FillOmegaOccOmega2n(const Operator &Eta, ModelSpace *ms,
                         std::deque<arma::mat> &T_pp,
                         std::deque<arma::mat> &T_hh);
void TensorChannelPairs(const Operator &Eta, ModelSpace *ms,
                        std::vector<size_t> &ch_bra_list,
                        std::vector<size_t> &ch_ket_list);
} // namespace
// factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
// Eta is tensor (reduced). Gamma and Z are scalar unreduced.
void comm223_231_st(const Operator &Eta_in, const Operator &Gamma,
                    Operator &Z) {
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
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

  ModelSpace *ms = Z.modelspace;
  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  const int nch = ms->GetNumberTwoBodyChannels();
  const int norb1b = (int)Z.OneBody.n_rows;

  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };
  auto same_1b_channel = [](Orbit &oa, Orbit &ob) {
    return oa.j2 == ob.j2 and oa.l == ob.l and oa.tz2 == ob.tz2;
  };

  std::vector<index_t> allorb_vec(ms->all_orbits.begin(),
                                  ms->all_orbits.end());

  // ==================================================================
  // f^I via scalar chi^alpha   (FDC.cc Chi_221_a, tensor Ω)
  // ==================================================================
  if (use_TypeI_1b) {
    // χ^α_de = 1/ĵ_d² Σ_{c J0} [ n_c n_d T_pp(cd;ce) − n̄_c n̄_d T_hh(cd;ce) ]
    //   T = Σ_{J1} (−1)^{J0+J1+λ} λ̂^{-1} Ω^{J0J1} W Ω^{J1J0}  (channel DGEMM,
    //   FillOmegaOccOmega2n). occ: n̄_a n̄_b n_c n_d − n_a n_b n̄_c n̄_d.
    std::deque<arma::mat> T_pp, T_hh;
    FillOmegaOccOmega2n(Eta, ms, T_pp, T_hh);

    arma::mat Chi_alpha(norb1b, norb1b, arma::fill::zeros);
#pragma omp parallel
    {
      arma::mat Chi_loc(norb1b, norb1b, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ch = 0; ch < nch; ++ch) {
        TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
        const int nK = tbc.GetNumberKets();
        for (int ibra = 0; ibra < 2 * nK; ++ibra) {
          auto cd = PqTB(tbc, ibra, nK);
          if (cd[0] < 0)
            continue;
          const int c = cd[0], d = cd[1];
          Orbit &oc = ms->GetOrbit(c);
          Orbit &od = ms->GetOrbit(d);
          const double n_c = oc.occ, nbar_c = 1.0 - n_c;
          const double n_d = od.occ, nbar_d = 1.0 - n_d;
          const double hatj2_inv = 1.0 / (od.j2 + 1.0);
          for (auto e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) {
            const int iket = Index2n(tbc, c, (int)e);
            if (iket < 0)
              continue;
            Chi_loc(d, e) += hatj2_inv * (n_c * n_d * T_pp[ch](ibra, iket) -
                                          nbar_c * nbar_d * T_hh[ch](ibra, iket));
          }
        }
      }
#pragma omp critical
      { Chi_alpha += Chi_loc; }
    }
    T_pp.clear();
    T_hh.clear();

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chialpha"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // f^I (scalar final, FDC.cc A2a):
    //   Z_pq = 1/(2 ĵ_p²) Σ_{J de} (2J+1) (χ_de + χ_ed) Γ^J(ep;dq)
    // on the ordered pairs of each Γ channel block (no per-element lookups).
    // Γ Hermitian ⇒ the (p,q) sum is already hZ-symmetric.
    const arma::mat Chi_sym = Chi_alpha + Chi_alpha.t();
    arma::mat Z1(norb1b, norb1b, arma::fill::zeros);
#pragma omp parallel
    {
      arma::mat Z_loc(norb1b, norb1b, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ch = 0; ch < nch; ++ch) {
        TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
        const int nK = tbc.GetNumberKets();
        if (nK < 1)
          continue;
        const arma::mat Gam2 = Omega2n(Gamma, ms, ch, ch);
        const double twoJp1 = 2.0 * tbc.J + 1.0;
        for (int ibra = 0; ibra < 2 * nK; ++ibra) {
          auto ep = PqTB(tbc, ibra, nK);
          if (ep[0] < 0)
            continue;
          Orbit &oe = ms->GetOrbit(ep[0]);
          Orbit &op = ms->GetOrbit(ep[1]);
          for (int iket = 0; iket < 2 * nK; ++iket) {
            auto dq = PqTB(tbc, iket, nK);
            if (dq[0] < 0)
              continue;
            Orbit &od = ms->GetOrbit(dq[0]);
            Orbit &oq = ms->GetOrbit(dq[1]);
            if (not same_1b_channel(op, oq) or not same_1b_channel(od, oe))
              continue;
            Z_loc(ep[1], dq[1]) +=
                twoJp1 * Chi_sym(dq[0], ep[0]) * Gam2(ibra, iket);
          }
        }
      }
#pragma omp critical
      { Z1 += Z_loc; }
    }
    for (auto p : allorb_vec) {
      Orbit &op = ms->GetOrbit(p);
      for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
        Z.OneBody(p, q) += 0.5 * Z1(p, q) / (op.j2 + 1.0);
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
    std::vector<size_t> ch_bra_list, ch_ket_list;
    TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
    const int npairs = (int)ch_bra_list.size();

    // χ^β: rectangular ladder DGEMM  Γ^{J0} W Ω^{J0 J1 λ}  then 6j + trace c
    //   χ_de^λ = 1/2 Σ_{c J0 J1} (−1)^{je+λ+J0+jc} Ĵ0 Ĵ1
    //            {λ J1 J0; jc jd je} (n̄_c n̄_e T_hh − n_c n_e T_pp)
    //   T_hh = Γ^{J0}(cd;ab) n_a n_b Ω^{J0J1}(ab;ce),  T_pp with n̄_a n̄_b
    // on the {ch0, ch1} blocks of TensorChannelPairs. No ÷ĵ_d²: reduced 1b.
    arma::mat Chi_beta(norb1b, norb1b, arma::fill::zeros);
#pragma omp parallel
    {
      arma::mat Chi_loc(norb1b, norb1b, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ip = 0; ip < npairs; ++ip) {
        const int ch0 = (int)ch_bra_list[ip], ch1 = (int)ch_ket_list[ip];
        TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
        TwoBodyChannel &tbc1 = ms->GetTwoBodyChannel(ch1);
        const int n0 = tbc0.GetNumberKets(), n1 = tbc1.GetNumberKets();
        if (n0 < 1 or n1 < 1)
          continue;
        const int J0 = tbc0.J, J1 = tbc1.J;
        arma::vec w_pp, w_hh;
        OccWeights2n(ms, ch0, w_pp, w_hh);
        const arma::mat Gam2 = Omega2n(Gamma, ms, ch0, ch0); // Γ^{J0}(cd;ab)
        const arma::mat Om = Omega2n(Eta, ms, ch0, ch1);     // Ω^{J0J1}(ab;ce)
        const double hats01 = hat(J0) * hat(J1);
        const arma::mat T_hh = hats01 * (Gam2 * arma::diagmat(w_hh) * Om);
        const arma::mat T_pp = hats01 * (Gam2 * arma::diagmat(w_pp) * Om);

        for (int ibra = 0; ibra < 2 * n0; ++ibra) {
          auto cd = PqTB(tbc0, ibra, n0);
          if (cd[0] < 0)
            continue;
          const int c = cd[0], d = cd[1];
          Orbit &oc = ms->GetOrbit(c);
          Orbit &od = ms->GetOrbit(d);
          const double jc = oc.j2 * 0.5, jd = od.j2 * 0.5;
          const double n_c = oc.occ, nbar_c = 1.0 - n_c;
          for (auto e : allorb_vec) {
            Orbit &oe = ms->GetOrbit(e);
            if (not AngMom::Triangle(jd, oe.j2 * 0.5, (double)lambda))
              continue;
            const int iket = Index2n(tbc1, c, (int)e);
            if (iket < 0)
              continue;
            const double n_e = oe.occ, nbar_e = 1.0 - n_e;
            // occ: n_a n_b n̄_c n̄_e − n̄_a n̄_b n_c n_e
            const double block =
                nbar_c * nbar_e * T_hh(ibra, iket) - n_c * n_e * T_pp(ibra, iket);
            if (std::abs(block) < 1e-16)
              continue;
            const double sixj =
                ms->GetSixJ(lambda, J1, J0, jc, jd, oe.j2 * 0.5);
            if (std::abs(sixj) < 1e-16)
              continue;
            const double phase = ms->phase((oe.j2 + oc.j2) / 2 + lambda + J0);
            Chi_loc(d, e) += 0.5 * phase * sixj * block;
          }
        }
      }
#pragma omp critical
      { Chi_beta += Chi_loc; }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_231_eths_chibeta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // f^II (scalar final): tensor χ^β × tensor Ω → scalar via RME [χ⊗Ω]^(0)
    //   Z_pq = 1/ĵ_p² Σ_{de J3 J4} (−1)^{jp+J4+jd} Ĵ3 Ĵ4 / λ̂
    //          {J4 λ J3; je jp jd} χ_de (Ω^{J3J4}(ep;dq) + hΓ Ω^{J3J4}(eq;dp))
    //        = (Z1 + hΓ Z1ᵀ)/ĵ_p²,  Z1 = the first string over ordered (p,q),
    // walked on the Ω^{J3 J4} blocks of TensorChannelPairs.
    // AMC f2a: combined (−1)^{j_i+J_4+j_d} ≡ phase((jp2+jd2)/2+J4). Do NOT
    // multiply an extra (−1)^{j_i} (see factorized_code_analyze.tex).
    arma::mat Z1(norb1b, norb1b, arma::fill::zeros);
#pragma omp parallel
    {
      arma::mat Z_loc(norb1b, norb1b, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ip = 0; ip < npairs; ++ip) {
        const int ch3 = (int)ch_bra_list[ip], ch4 = (int)ch_ket_list[ip];
        TwoBodyChannel &tbc3 = ms->GetTwoBodyChannel(ch3);
        TwoBodyChannel &tbc4 = ms->GetTwoBodyChannel(ch4);
        const int n3 = tbc3.GetNumberKets(), n4 = tbc4.GetNumberKets();
        if (n3 < 1 or n4 < 1)
          continue;
        const int J3 = tbc3.J, J4 = tbc4.J;
        const arma::mat Om = Omega2n(Eta, ms, ch3, ch4); // Ω^{J3J4}(ep;dq)
        const double pref = hat(J3) * hat(J4) * hat_lambda_inv;
        for (int ibra = 0; ibra < 2 * n3; ++ibra) {
          auto ep = PqTB(tbc3, ibra, n3);
          if (ep[0] < 0)
            continue;
          Orbit &oe = ms->GetOrbit(ep[0]);
          Orbit &op = ms->GetOrbit(ep[1]);
          const double je = oe.j2 * 0.5, jp = op.j2 * 0.5;
          for (int iket = 0; iket < 2 * n4; ++iket) {
            auto dq = PqTB(tbc4, iket, n4);
            if (dq[0] < 0)
              continue;
            Orbit &od = ms->GetOrbit(dq[0]);
            Orbit &oq = ms->GetOrbit(dq[1]);
            if (not same_1b_channel(op, oq))
              continue;
            const double chi = Chi_beta(dq[0], ep[0]);
            if (std::abs(chi) < 1e-16)
              continue;
            const double sixj =
                ms->GetSixJ(J4, lambda, J3, je, jp, od.j2 * 0.5);
            if (std::abs(sixj) < 1e-16)
              continue;
            const double phase = ms->phase((op.j2 + od.j2) / 2 + J4);
            Z_loc(ep[1], dq[1]) += phase * pref * sixj * chi * Om(ibra, iket);
          }
        }
      }
#pragma omp critical
      { Z1 += Z_loc; }
    }
    for (auto p : allorb_vec) {
      Orbit &op = ms->GetOrbit(p);
      for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
        Z.OneBody(p, q) +=
            (Z1(p, q) + hGamma * Z1(q, p)) / (op.j2 + 1.0); // unreduced 1b
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
// Defined after ChiOp2n / BuildChiGammaPathB (needs the no-Pauli layout).
static void comm223_231_fIIIa_leftover_dgemm(const Operator &Eta,
                                            const Operator &Gamma,
                                            Operator &Z);

void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {

  double t_internal = omp_get_wtime();
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();

  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;

  int norbits = Z.modelspace->all_orbits.size();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),
                                  Z.modelspace->all_orbits.end());

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

    // T_pp/T_hh on the ordered pairs of every channel (shared with χ^α, χ^θ);
    // χ^δ only needs the stored-ket block [0,n)×[0,n).
    std::deque<arma::mat> T_pp_all, T_hh_all;
    FillOmegaOccOmega2n(Eta, Z.modelspace, T_pp_all, T_hh_all);

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch; ++ch) {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      if (nKets < 1)
        continue;
      arma::mat &Chi_mat = Chi_delta.GetMatrix(ch, ch);
      const arma::mat T_pp = T_pp_all[ch].submat(0, 0, nKets - 1, nKets - 1);
      const arma::mat T_hh = T_hh_all[ch].submat(0, 0, nKets - 1, nKets - 1);
      T_pp_all[ch].clear();
      T_hh_all[ch].clear();

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


void comm223_232(const Operator &Eta_in, const Operator &Gamma, Operator &Z) {
  // Ω is reduced RME (including λ=0). Gamma and Z are scalar unreduced.
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
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
          // Trace on the spectator c: walk the kets (c,p) of ch0 in both
          // orderings, q in the 1b channel of p (normalized T → GetTBME units).
          for (int icp = 0; icp < n0; ++icp) {
            Ket &kcp = tbc0.GetKet(icp);
            for (int swap = 0; swap < 2; ++swap) {
              if (swap and kcp.p == kcp.q)
                break;
              const index_t c = swap ? kcp.q : kcp.p;
              const index_t p = swap ? kcp.p : kcp.q;
              Orbit &oc = Z.modelspace->GetOrbit(c);
              Orbit &op = Z.modelspace->GetOrbit(p);
              const double n_c = oc.occ, nbar_c = 1.0 - n_c;
              const double ph_cp = (c > p) ? kcp.Phase(J0) : 1.0;
              const double N_cp = (c == p) ? PhysConst::SQRT2 : 1.0;
              const double hatj2 = op.j2 + 1.0; // j_p == j_q
              for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) {
                const int icq = static_cast<int>(
                    tbc0.GetLocalIndex(std::min(c, q), std::max(c, q)));
                if (icq < 0 or icq >= n0)
                  continue;
                Ket &kcq = tbc0.GetKet(icq);
                const double ph_cq = (c > q) ? kcq.Phase(J0) : 1.0;
                const double N_cq = (c == q) ? PhysConst::SQRT2 : 1.0;
                const double me =
                    n_c * T_pp(icp, icq) + nbar_c * T_hh(icp, icq);
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
    ModelSpace *ms = Z.modelspace;
    std::vector<index_t> allorb(ms->all_orbits.begin(), ms->all_orbits.end());
    const int norb1b = (int)Gamma.OneBody.n_rows;

    std::vector<size_t> pair_bra, pair_ket; // {ch0, ch1} coupled by Ω
    TensorChannelPairs(Eta, ms, pair_bra, pair_ket);
    const int npairs = (int)pair_bra.size();
    const int nch_tb = ms->GetNumberTwoBodyChannels();
    // partners_of[ch]   = {ch2 : Ω^{J2 J}(ch2;ch) ≠ 0}   (bra-side fold)
    // partners_from[ch] = {ch2 : Ω^{J J2}(ch;ch2) ≠ 0}   (ket-side fold)
    std::vector<std::vector<int>> partners_of(nch_tb), partners_from(nch_tb);
    for (int ip = 0; ip < npairs; ++ip) {
      partners_of[pair_ket[ip]].push_back((int)pair_bra[ip]);
      partners_from[pair_bra[ip]].push_back((int)pair_ket[ip]);
    }

    arma::mat Chi_zeta(norb1b, norb1b, arma::fill::zeros);
    arma::mat Chi_OG(norb1b, norb1b, arma::fill::zeros);

    // χ^ζ / χ^{ΩΓ}: same dummy (a,b) channel DGEMM as χ^β (Γ W Ω), then 6j on
    // the spectator c. Rows (c,i)@J0 in ch0, cols (c,j)@J1 in ch1:
    //   χ^ζ_ij   = 1/2 Σ_c (−)^{jj+jc+λ+J0} {λ J1 J0; jc ji jj}
    //              [ n̄_c TZ_hh + n_c TZ_pp ],   TZ  = Ĵ0Ĵ1 Γ^{J0} W0 Ω^{J0J1}
    //   χ^{ΩΓ}_ij = 1/2 Σ_c (−)^{jj+jc+λ+J0} {J0 J1 λ; jj ji jc}
    //              [ n̄_c TOG_hh + n_c TOG_pp ], TOG = Ĵ0Ĵ1 Ω^{J0J1} W1 Γ^{J1}
#pragma omp parallel
    {
      arma::mat Zeta_loc(norb1b, norb1b, arma::fill::zeros);
      arma::mat OG_loc(norb1b, norb1b, arma::fill::zeros);
#pragma omp for schedule(dynamic, 1)
      for (int ip = 0; ip < npairs; ++ip) {
        const int ch0 = (int)pair_bra[ip], ch1 = (int)pair_ket[ip];
        TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
        TwoBodyChannel &tbc1 = ms->GetTwoBodyChannel(ch1);
        const int n0 = tbc0.GetNumberKets(), n1 = tbc1.GetNumberKets();
        if (n0 < 1 or n1 < 1)
          continue;
        const int J0 = tbc0.J, J1 = tbc1.J;
        const double hats01 = hat(J0) * hat(J1);
        arma::vec w_pp0, w_hh0, w_pp1, w_hh1;
        OccWeights2n(ms, ch0, w_pp0, w_hh0);
        OccWeights2n(ms, ch1, w_pp1, w_hh1);
        const arma::mat Om = Omega2n(Eta, ms, ch0, ch1);     // Ω^{J0J1}
        const arma::mat G0 = Omega2n(Gamma, ms, ch0, ch0);   // Γ^{J0}
        const arma::mat G1 = Omega2n(Gamma, ms, ch1, ch1);   // Γ^{J1}
        const arma::mat TZ_hh = hats01 * (G0 * arma::diagmat(w_hh0) * Om);
        const arma::mat TZ_pp = hats01 * (G0 * arma::diagmat(w_pp0) * Om);
        const arma::mat TOG_hh = hats01 * (Om * arma::diagmat(w_hh1) * G1);
        const arma::mat TOG_pp = hats01 * (Om * arma::diagmat(w_pp1) * G1);

        for (int ibra = 0; ibra < 2 * n0; ++ibra) {
          auto ci = PqTB(tbc0, ibra, n0);
          if (ci[0] < 0)
            continue;
          const int c = ci[0], i = ci[1];
          Orbit &oc = ms->GetOrbit(c);
          Orbit &oi = ms->GetOrbit(i);
          const double jc = oc.j2 * 0.5, ji = oi.j2 * 0.5;
          const double n_c = oc.occ, nbar_c = 1.0 - n_c;
          for (auto j : allorb) {
            Orbit &oj = ms->GetOrbit(j);
            const double jj = oj.j2 * 0.5;
            if (not AngMom::Triangle(ji, jj, (double)lambda))
              continue;
            const int iket = Index2n(tbc1, c, (int)j);
            if (iket < 0)
              continue;
            const double ph = ms->phase((oj.j2 + oc.j2) / 2 + lambda + J0);
            const double blkZ =
                nbar_c * TZ_hh(ibra, iket) + n_c * TZ_pp(ibra, iket);
            if (std::abs(blkZ) > 1e-16) {
              const double sixj = ms->GetSixJ(lambda, J1, J0, jc, ji, jj);
              Zeta_loc(i, j) += 0.5 * ph * sixj * blkZ;
            }
            const double blkOG =
                nbar_c * TOG_hh(ibra, iket) + n_c * TOG_pp(ibra, iket);
            if (std::abs(blkOG) > 1e-16) {
              const double sixj = ms->GetSixJ(J0, J1, lambda, jj, ji, jc);
              OG_loc(i, j) += 0.5 * ph * sixj * blkOG;
            }
          }
        }
      }
#pragma omp critical
      {
        Chi_zeta += Zeta_loc;
        Chi_OG += OG_loc;
      }
    }

    if (Commutator::verbose) {
      Z.profiler.timer["_232_eths_chizeta"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }

    // Fold: Z = −W_OG − V_ζ. W: AMC Wbra with χ^{ΩΓ}; V: Wket with χ^ζ.
    //   W(ij;kl) = Σ_{J2 a} Ĵ2/(Ĵ λ̂) [ (−)^{ji+ja+J2}   {J J2 λ; ja jj ji} χ^{ΩΓ}_ja Ω^{J2J}(ia;kl)
    //                                + (−)^{J+ji+ja+J2} {J J2 λ; ja ji jj} χ^{ΩΓ}_ia Ω^{J2J}(ja;kl) ]
    //   V(ij;kl) = Σ_{J2 a} Ĵ2/(Ĵ λ̂) [ (−)^{jl+ja+J2}   {J J2 λ; ja jk jl} χ^ζ_ak Ω^{JJ2}(ij;al)
    //                                + (−)^{J+jl+ja+J2} {J J2 λ; ja jl jk} χ^ζ_al Ω^{JJ2}(ij;ak) ]
    // as DGEMM per channel: W = Σ_{ch2} C_W(bras × 2n2) · Ω(ch2,ch)[:, kets],
    //                       V = Σ_{ch2} Ω(ch,ch2)[bras, :] · D_V(2n2 × kets).
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nch; ++ich) {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = ms->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = ms->GetTwoBodyChannel(ch_ket);
      const int J = tbc_bra.J;
      if (tbc_ket.J != J)
        continue;
      const int nbras = tbc_bra.GetNumberKets();
      const int nkets = tbc_ket.GetNumberKets();
      if (nbras < 1 or nkets < 1)
        continue;
      const double hatJ_inv = 1.0 / hat(J);

      arma::mat W(nbras, nkets, arma::fill::zeros);
      arma::mat V(nbras, nkets, arma::fill::zeros);

      // bra side: Ω^{J2 J}((i,a);(k,l)), ch2 → ch_ket
      for (int ch2 : partners_of[ch_ket]) {
        TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
        const int n2 = tbc2.GetNumberKets();
        if (n2 < 1)
          continue;
        const int J2 = tbc2.J;
        const double pref = hatJ_inv * hat(J2) * hat_lambda_inv;
        arma::mat C_W(nbras, 2 * n2, arma::fill::zeros);
        for (int ibra = 0; ibra < nbras; ++ibra) {
          Ket &bra = tbc_bra.GetKet(ibra);
          const int i = (int)bra.p, j = (int)bra.q;
          Orbit &oi = ms->GetOrbit(i);
          Orbit &oj = ms->GetOrbit(j);
          const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
          for (auto a : allorb) {
            Orbit &oa = ms->GetOrbit(a);
            const double ja = oa.j2 * 0.5;
            const int r_ia = Index2n(tbc2, i, (int)a);
            if (r_ia >= 0 and std::abs(Chi_OG(j, a)) > 1e-16)
              C_W(ibra, r_ia) += ms->phase((oi.j2 + oa.j2) / 2 + J2) * pref *
                                 ms->GetSixJ(J, J2, lambda, ja, jj, ji) *
                                 Chi_OG(j, a);
            const int r_ja = Index2n(tbc2, j, (int)a);
            if (r_ja >= 0 and std::abs(Chi_OG(i, a)) > 1e-16)
              C_W(ibra, r_ja) +=
                  ms->phase(J + (oi.j2 + oa.j2) / 2 + J2) * pref *
                  ms->GetSixJ(J, J2, lambda, ja, ji, jj) * Chi_OG(i, a);
          }
        }
        const arma::mat Om = Omega2n(Eta, ms, ch2, (int)ch_ket);
        W += C_W * Om.cols(0, nkets - 1);
      }

      // ket side: Ω^{J J2}((i,j);(a,l)), ch_bra → ch2
      for (int ch2 : partners_from[ch_bra]) {
        TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
        const int n2 = tbc2.GetNumberKets();
        if (n2 < 1)
          continue;
        const int J2 = tbc2.J;
        const double pref = hatJ_inv * hat(J2) * hat_lambda_inv;
        arma::mat D_V(2 * n2, nkets, arma::fill::zeros);
        for (int iket = 0; iket < nkets; ++iket) {
          Ket &ket = tbc_ket.GetKet(iket);
          const int k = (int)ket.p, l = (int)ket.q;
          Orbit &ok = ms->GetOrbit(k);
          Orbit &ol = ms->GetOrbit(l);
          const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
          for (auto a : allorb) {
            Orbit &oa = ms->GetOrbit(a);
            const double ja = oa.j2 * 0.5;
            const int c_al = Index2n(tbc2, (int)a, l);
            if (c_al >= 0 and std::abs(Chi_zeta(a, k)) > 1e-16)
              D_V(c_al, iket) += ms->phase((ol.j2 + oa.j2) / 2 + J2) * pref *
                                 ms->GetSixJ(J, J2, lambda, ja, jk, jl) *
                                 Chi_zeta(a, k);
            const int c_ak = Index2n(tbc2, (int)a, k);
            if (c_ak >= 0 and std::abs(Chi_zeta(a, l)) > 1e-16)
              D_V(c_ak, iket) +=
                  ms->phase(J + (ol.j2 + oa.j2) / 2 + J2) * pref *
                  ms->GetSixJ(J, J2, lambda, ja, jl, jk) * Chi_zeta(a, l);
          }
        }
        const arma::mat Om = Omega2n(Eta, ms, (int)ch_bra, ch2);
        V += Om.rows(0, nbras - 1) * D_V;
      }

      for (int ibra = 0; ibra < nbras; ++ibra) {
        Ket &bra = tbc_bra.GetKet(ibra);
        const int ketmin = (ch_bra == ch_ket) ? ibra : 0;
        for (int iket = ketmin; iket < nkets; ++iket) {
          Ket &ket = tbc_ket.GetKet(iket);
          double z = -W(ibra, iket) - V(ibra, iket); // Γ^II = −W_OG − V_ζ
          if (bra.p == bra.q)
            z /= PhysConst::SQRT2;
          if (ket.p == ket.q)
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
  //   ~1227 allocate bar_Eta, bar_Gamma
  //   ~1243 combined Pandya + nnnbar_Eta / nnnbar_Eta_d occupancy
  //   ~1385 DGEMM bar_CHI_V, bar_CHI_VI_II  (η is AMC Path B; θ,λ diagram-local)
  //   ~1413 GIIIa  χ^η ladder
  //   ~1528 GIVa   χ^κ   inverse Pandya on ch_bra×ch_ket, ladder DGEMM
  //   ~1906 GIVb   χ^ι   inverse Pandya, CC-channel cross DGEMM leftover
  //   ~2393 GIIIc  χ^θ
  //   ~2739 GIIIb  (ethS gold is leftover IIb/IId, not FDC recouple of χ^η)
  //   ~2916 GIVc   χ^λ  (Pandya of χ^λ and Ω together, same 9j)
  PandyaBars232 bars;
  if (use_TypeGIVa_2b or use_TypeGIVb_2b)
    FillPandyaBars232(Eta, Gamma, Z, bars);

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
/// Shared channel helpers for the 232 leftover diagrams.
///
/// Rule: [Ω×Ω]^0 → scalar χ (λ_χ=0), but the tensor χ are not antisymmetric,
/// so they live in 2n ordered-pair blocks (PqTB / PqCC / ChiOp2n), never in a
/// ket-ordered TwoBodyME. GIIIa: χ^η ladder DGEMM. GIIIb: χ̄^η → RC → Γ̄·RC →
/// inverse Pandya (FDC.cc chain). GIIIc: χ^θ channel DGEMM → Pandya → Γ̄ → inv.
////////////////////////////////////////////////////////////////////////////
namespace {

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

double HatJ(int J) { return std::sqrt(2.0 * J + 1.0); }

// Ordinary-channel pairs {ch_bra, ch_ket} coupled by the rank-(λ,T,π) Ω, in
// both orderings (TwoBodyME::Allocate keeps ch_bra ≤ ch_ket only, but the χ
// intermediates are not Hermitian). FDC.cc ch_bra_list / ch_ket_list.
void TensorChannelPairs(const Operator &Eta, ModelSpace *ms,
                        std::vector<size_t> &ch_bra_list,
                        std::vector<size_t> &ch_ket_list) {
  const int nch = (int)ms->GetNumberTwoBodyChannels();
  const int lambda = Eta.GetJRank();
  const int rank_T = Eta.GetTRank();
  const int parity = Eta.GetParity();
  ch_bra_list.clear();
  ch_ket_list.clear();
  for (int ch_bra = 0; ch_bra < nch; ++ch_bra) {
    TwoBodyChannel &tbc_bra = ms->GetTwoBodyChannel(ch_bra);
    for (int ch_ket = 0; ch_ket < nch; ++ch_ket) {
      TwoBodyChannel &tbc_ket = ms->GetTwoBodyChannel(ch_ket);
      if (not AngMom::Triangle(tbc_bra.J, tbc_ket.J, lambda))
        continue;
      if (std::abs(tbc_bra.Tz - tbc_ket.Tz) != rank_T)
        continue;
      if ((tbc_bra.parity + tbc_ket.parity + parity) % 2 != 0)
        continue;
      ch_bra_list.push_back((size_t)ch_bra);
      ch_ket_list.push_back((size_t)ch_ket);
    }
  }
}

// 2n index of the ordered pair (p,q) in an ordinary channel (PqTB layout:
// [0,n) stored kets p ≤ q, [n,2n) the swapped pair); −1 if absent.
int Index2n(TwoBodyChannel &tbc, int p, int q) {
  const int nK = tbc.GetNumberKets();
  const int loc = tbc.GetLocalIndex(std::min(p, q), std::max(p, q));
  if (loc < 0 or loc >= nK)
    return -1;
  return (p > q) ? loc + nK : loc;
}

// Same 2n index on a CC channel (PqCC). Looks up the unordered ket first so a
// missing pair is −1, not NumberKets−1 (GetLocalIndex(q,p) adds nK to −1).
int IndexCC(TwoBodyChannel_CC &tbc, int p, int q) {
  const int nK = tbc.GetNumberKets();
  const size_t loc = tbc.GetLocalIndex(std::min(p, q), std::max(p, q));
  if (loc >= (size_t)nK)
    return -1;
  return (p > q) ? (int)loc + nK : (int)loc;
}

// X^{J0 J1}(pq;rs) on ordered pairs (2n0 × 2n1) in GetTBME_J units, expanded
// from the stored normalized block (√2 for p==q, Ket phase for p>q; partner
// block ch0>ch1 is h (−1)^{J0−J1} Mᵀ exactly as TwoBodyME::GetTBME). Works for
// the tensor Ω (any ch0,ch1 of TensorChannelPairs) and the scalar Γ (ch0==ch1).
arma::mat Omega2n(const Operator &Op, ModelSpace *ms, int ch0, int ch1) {
  TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
  TwoBodyChannel &tbc1 = ms->GetTwoBodyChannel(ch1);
  const int n0 = tbc0.GetNumberKets(), n1 = tbc1.GetNumberKets();
  arma::mat X(2 * n0, 2 * n1, arma::fill::zeros);
  if (n0 < 1 or n1 < 1)
    return X;
  const size_t clo = (size_t)std::min(ch0, ch1), chi = (size_t)std::max(ch0, ch1);
  auto it = Op.TwoBody.MatEl.find({clo, chi});
  if (it == Op.TwoBody.MatEl.end())
    return X;
  arma::mat M;
  if (ch0 <= ch1)
    M = it->second;
  else
    M = (Op.IsHermitian() ? 1.0 : -1.0) * ms->phase(tbc0.J - tbc1.J) *
        it->second.t();
  if ((int)M.n_rows != n0 or (int)M.n_cols != n1)
    return X;
  arma::vec f0_lo(n0), f0_hi(n0), f1_lo(n1), f1_hi(n1);
  for (int i = 0; i < n0; ++i) {
    Ket &k = tbc0.GetKet(i);
    f0_lo(i) = (k.p == k.q) ? PhysConst::SQRT2 : 1.0;
    f0_hi(i) = (k.p == k.q) ? 0.0 : k.Phase(tbc0.J);
  }
  for (int i = 0; i < n1; ++i) {
    Ket &k = tbc1.GetKet(i);
    f1_lo(i) = (k.p == k.q) ? PhysConst::SQRT2 : 1.0;
    f1_hi(i) = (k.p == k.q) ? 0.0 : k.Phase(tbc1.J);
  }
  const arma::mat Mlo = M * arma::diagmat(f1_lo);
  const arma::mat Mhi = M * arma::diagmat(f1_hi);
  X.submat(0, 0, n0 - 1, n1 - 1) = arma::diagmat(f0_lo) * Mlo;
  X.submat(0, n1, n0 - 1, 2 * n1 - 1) = arma::diagmat(f0_lo) * Mhi;
  X.submat(n0, 0, 2 * n0 - 1, n1 - 1) = arma::diagmat(f0_hi) * Mlo;
  X.submat(n0, n1, 2 * n0 - 1, 2 * n1 - 1) = arma::diagmat(f0_hi) * Mhi;
  return X;
}

// Occupancy weights on the ordered pairs of an ordinary channel:
//   w_pp(ab) = n̄_a n̄_b,   w_hh(ab) = n_a n_b
void OccWeights2n(ModelSpace *ms, int ch, arma::vec &w_pp, arma::vec &w_hh) {
  TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
  const int nK = tbc.GetNumberKets();
  w_pp.zeros(2 * nK);
  w_hh.zeros(2 * nK);
  for (int m = 0; m < 2 * nK; ++m) {
    auto ab = PqTB(tbc, m, nK);
    if (ab[0] < 0)
      continue;
    const double na = ms->GetOrbit(ab[0]).occ;
    const double nb = ms->GetOrbit(ab[1]).occ;
    w_pp(m) = (1.0 - na) * (1.0 - nb);
    w_hh(m) = na * nb;
  }
}

// Scalar [Ω⊗Ω]^(0) ladder core shared by χ^α (f^I), χ^δ (f^{III_b}) and
// χ^θ (G^{III_c}); per ordinary channel, 2n layout, GetTBME_J units:
//   T^{J0}_{pp|hh}(ij;kl) = Σ_{J1} (−1)^{J0+J1+λ} λ̂^{-1}
//                            Ω^{J0J1}(ij;ab) w_{pp|hh}(ab) Ω^{J1J0}(ab;kl)
// The diagrams differ only in which occupancies multiply T afterwards.
void FillOmegaOccOmega2n(const Operator &Eta, ModelSpace *ms,
                         std::deque<arma::mat> &T_pp,
                         std::deque<arma::mat> &T_hh) {
  const int nch = ms->GetNumberTwoBodyChannels();
  const int lambda = Eta.GetJRank();
  const double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  std::vector<size_t> ch_bra_list, ch_ket_list;
  TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
  std::vector<std::vector<int>> partners(nch);
  for (size_t ip = 0; ip < ch_bra_list.size(); ++ip)
    partners[ch_bra_list[ip]].push_back((int)ch_ket_list[ip]);

  T_pp.assign(nch, arma::mat());
  T_hh.assign(nch, arma::mat());
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch0 = 0; ch0 < nch; ++ch0) {
    TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
    const int n0 = tbc0.GetNumberKets();
    const int J0 = tbc0.J;
    T_pp[ch0] = arma::mat(2 * n0, 2 * n0, arma::fill::zeros);
    T_hh[ch0] = arma::mat(2 * n0, 2 * n0, arma::fill::zeros);
    if (n0 < 1)
      continue;
    for (int ch1 : partners[ch0]) {
      TwoBodyChannel &tbc1 = ms->GetTwoBodyChannel(ch1);
      if (tbc1.GetNumberKets() < 1)
        continue;
      arma::vec w_pp, w_hh;
      OccWeights2n(ms, ch1, w_pp, w_hh);
      const arma::mat Om = Omega2n(Eta, ms, ch0, ch1);   // Ω^{J0 J1}(ij;ab)
      const arma::mat Om_r = Omega2n(Eta, ms, ch1, ch0); // Ω^{J1 J0}(ab;kl)
      const double ang = ms->phase(J0 + tbc1.J + lambda) * hat_lambda_inv;
      T_pp[ch0] += ang * (Om * arma::diagmat(w_pp) * Om_r);
      T_hh[ch0] += ang * (Om * arma::diagmat(w_hh) * Om_r);
    }
  }
}

// ModelSpace::PreCalculateNineJ only covers the LS-jj 9j. The tensor Pandya
// and its inverse need two more shapes (j half-integer, J integer):
//   A  { ja jd J1 }        B  { λ  J0 J1 }
//      { jb jc J2 }           { J3 jj jk }
//      { Jb Jk λ  }           { J2 ji jl }
// GetNineJ on a cache miss calls GetSixJ with shapes that are not cached,
// which aborts inside a parallel region. Reserve the keys serially and fill
// them in parallel, exactly like ModelSpace::PreCalculateNineJ.
void PreCalculateNineJTensorPandya(ModelSpace *ms, int lambda) {
  static std::set<std::array<int, 2>> done; // {max_j2, λ}; NineJList is static
  int max_j2 = 0;
  for (auto o : ms->all_orbits)
    max_j2 = std::max(max_j2, ms->GetOrbit(o).j2);
  if (done.count({max_j2, lambda}))
    return;
  const double t_start = omp_get_wtime();
  std::vector<uint64_t> KEYS;
  auto reserve = [&](uint64_t key) {
    if (ms->NineJList.count(key) == 0) {
      KEYS.push_back(key);
      ms->NineJList[key] = 0.0;
    }
  };
  auto tri = [](int j2a, int j2b, int &Jmin, int &Jmax) {
    Jmin = std::abs(j2a - j2b) / 2;
    Jmax = (j2a + j2b) / 2;
  };
  for (int j2a = 1; j2a <= max_j2; j2a += 2)
    for (int j2b = 1; j2b <= max_j2; j2b += 2)
      for (int j2c = 1; j2c <= max_j2; j2c += 2)
        for (int j2d = 1; j2d <= max_j2; j2d += 2) {
          const double ja = 0.5 * j2a, jb = 0.5 * j2b, jc = 0.5 * j2c,
                       jd = 0.5 * j2d;
          int J1min, J1max, J2min, J2max, Jbmin, Jbmax, Jkmin, Jkmax;
          // A: rows (ja jd J1) (jb jc J2) (Jb Jk λ); cols (ja jb Jb) (jd jc Jk)
          tri(j2a, j2d, J1min, J1max);
          tri(j2b, j2c, J2min, J2max);
          tri(j2a, j2b, Jbmin, Jbmax);
          tri(j2d, j2c, Jkmin, Jkmax);
          for (int J1 = J1min; J1 <= J1max; ++J1)
            for (int J2 = std::max(J2min, std::abs(J1 - lambda));
                 J2 <= std::min(J2max, J1 + lambda); ++J2)
              for (int Jb = Jbmin; Jb <= Jbmax; ++Jb)
                for (int Jk = std::max(Jkmin, std::abs(Jb - lambda));
                     Jk <= std::min(Jkmax, Jb + lambda); ++Jk)
                  reserve(ms->NineJHash(ja, jd, J1, jb, jc, J2, Jb, Jk, lambda));
          // B: rows (λ J0 J1) (J3 jb jc) (J2 ja jd); cols (λ J3 J2) (J0 jb ja) (J1 jc jd)
          tri(j2a, j2b, J1min, J1max); // J0 ∈ tri(ja,jb)
          tri(j2c, j2d, J2min, J2max); // J1 ∈ tri(jc,jd)
          tri(j2a, j2d, Jbmin, Jbmax); // J2 ∈ tri(ja,jd)
          tri(j2b, j2c, Jkmin, Jkmax); // J3 ∈ tri(jb,jc)
          for (int J0 = J1min; J0 <= J1max; ++J0)
            for (int J1 = std::max(J2min, std::abs(J0 - lambda));
                 J1 <= std::min(J2max, J0 + lambda); ++J1)
              for (int J2 = Jbmin; J2 <= Jbmax; ++J2)
                for (int J3 = std::max(Jkmin, std::abs(J2 - lambda));
                     J3 <= std::min(Jkmax, J2 + lambda); ++J3)
                  reserve(ms->NineJHash(lambda, J0, J1, J3, jb, jc, J2, ja, jd));
        }
  const int nkeys = (int)KEYS.size();
#pragma omp parallel for schedule(dynamic, 64)
  for (int i = 0; i < nkeys; ++i) {
    uint64_t k[9];
    ms->NineJUnHash(KEYS[i], k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7],
                    k[8]);
    ms->NineJList[KEYS[i]] =
        AngMom::NineJ(0.5 * k[0], 0.5 * k[1], 0.5 * k[2], 0.5 * k[3],
                      0.5 * k[4], 0.5 * k[5], 0.5 * k[6], 0.5 * k[7],
                      0.5 * k[8]);
  }
  done.insert({max_j2, lambda});
  ms->profiler.timer["PreCalculateNineJTensorPandya"] +=
      omp_get_wtime() - t_start;
}

// FDC.cc keeps χ^κ, χ^ι after the inverse Pandya in `TwoBodyME Chi_VI_Op`
// etc. The tensor χ are not antisymmetric (χ(pq;..) ≠ −χ(qp;..), and χ(pp;..)
// exists for odd J), so TwoBodyME cannot hold them. ChiOp2n keeps the same
// (J,π,Tz) channel index and MatEl {ch_bra, ch_ket}, but rows/cols run over
// ordered pairs: [0,n) = (p,q) with p ≤ q (no Pauli), [n,2n) = (q,p) (row of
// p == q left zero, same as PqTB). Values in GetTBME_J units.
struct ChiOp2n {
  ModelSpace *ms = nullptr;
  int nch = 0, n_alloc = 0;
  std::vector<std::vector<std::array<int, 2>>> kets; // per channel, p ≤ q
  std::vector<std::vector<int>> local;               // p*n_alloc+q → idx
  CcMatMap MatEl;
  std::vector<arma::mat *> block; // nch*nch → &MatEl[{ch_bra,ch_ket}]

  void Init(ModelSpace *m) {
    ms = m;
    nch = (int)ms->GetNumberTwoBodyChannels();
    n_alloc = 0;
    for (auto o : ms->all_orbits)
      n_alloc = std::max(n_alloc, (int)o + 1);
    kets.assign((size_t)nch, {});
    local.assign((size_t)nch, {});
    std::vector<index_t> allorb(ms->all_orbits.begin(), ms->all_orbits.end());
    for (int ch = 0; ch < nch; ++ch) {
      TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
      local[ch].assign((size_t)n_alloc * n_alloc, -1);
      for (auto p : allorb) {
        Orbit &op = ms->GetOrbit(p);
        for (auto q : allorb) {
          if (q < p)
            continue;
          Orbit &oq = ms->GetOrbit(q);
          if ((op.l + oq.l) % 2 != tbc.parity)
            continue;
          if ((op.tz2 + oq.tz2) != 2 * tbc.Tz)
            continue;
          if (not AngMom::Triangle(op.j2 * 0.5, oq.j2 * 0.5, (double)tbc.J))
            continue;
          local[ch][(size_t)p * n_alloc + q] = (int)kets[ch].size();
          kets[ch].push_back({(int)p, (int)q});
        }
      }
    }
  }
  void Allocate(const std::vector<size_t> &ch_bra_list,
                const std::vector<size_t> &ch_ket_list) {
    MatEl.clear();
    block.assign((size_t)nch * nch, nullptr);
    for (size_t ich = 0; ich < ch_bra_list.size(); ++ich) {
      const int ch_bra = (int)ch_bra_list[ich], ch_ket = (int)ch_ket_list[ich];
      const int nbras = NumberKets(ch_bra), nkets = NumberKets(ch_ket);
      if (nbras < 1 or nkets < 1)
        continue;
      MatEl[{ch_bra, ch_ket}] =
          arma::mat(2 * nbras, 2 * nkets, arma::fill::zeros);
    }
    for (auto &kv : MatEl)
      block[(size_t)kv.first[0] * nch + kv.first[1]] = &kv.second;
  }
  int NumberKets(int ch) const { return (int)kets[ch].size(); }
  // Same convention as PqTB / PqCC.
  std::array<int, 2> Pq(int ch, int idx) const {
    const int n = NumberKets(ch);
    const auto &pq = kets[ch][idx % n];
    if (idx < n)
      return pq;
    if (pq[0] == pq[1])
      return {-1, -1};
    return {pq[1], pq[0]};
  }
  // 2n index of the ordered pair (p,q) in channel ch; −1 if absent.
  int Index(int ch, int p, int q) const {
    const int i = local[ch][(size_t)std::min(p, q) * n_alloc + std::max(p, q)];
    if (i < 0)
      return -1;
    if (p > q)
      return i + NumberKets(ch);
    return i;
  }
  int ChannelIndex(int J, int p, int q) const {
    Orbit &op = ms->GetOrbit(p);
    Orbit &oq = ms->GetOrbit(q);
    if (not AngMom::Triangle(op.j2 * 0.5, oq.j2 * 0.5, (double)J))
      return -1;
    const size_t ch = ms->GetTwoBodyChannelIndex(J, (op.l + oq.l) % 2,
                                                 (op.tz2 + oq.tz2) / 2);
    return ((int)ch < nch) ? (int)ch : -1;
  }
  arma::mat *GetBlock(int ch_bra, int ch_ket) const {
    if (ch_bra < 0 or ch_ket < 0)
      return nullptr;
    return block[(size_t)ch_bra * nch + ch_ket];
  }
  double GetTBME_J(int J0, int J1, int i, int j, int k, int l) const {
    const int ch_bra = ChannelIndex(J0, i, j);
    const int ch_ket = ChannelIndex(J1, k, l);
    arma::mat *M = GetBlock(ch_bra, ch_ket);
    if (M == nullptr)
      return 0.0;
    const int ibra = Index(ch_bra, i, j);
    const int iket = Index(ch_ket, k, l);
    if (ibra < 0 or iket < 0)
      return 0.0;
    return (*M)(ibra, iket);
  }
  // ChiTab-style accessor used by the GIVc Pandya of χ^λ.
  double operator()(index_t p, index_t q, index_t r, index_t s, int J0,
                    int J1) const {
    return GetTBME_J(J0, J1, (int)p, (int)q, (int)r, (int)s);
  }
};

// FDC.cc ~1556–1710 (Chi_VI_Op: inverse Pandya over ch_bra×ch_ket), tensor:
//   χ^{J0 J1}(ij;kl) = (−1)^{J0+ji+jk+λ} Ĵ0 Ĵ1 Σ_{J2 J3} (−1)^{J2} Ĵ2 Ĵ3
//                      { λ  J0 J1 }
//                      { J3 jj jk }  χ̄^{J2 J3}(il;kj)
//                      { J2 ji jl }
// (AMC Path B inverse without the printed leading minus; gold
// test_chi_eta_mscheme / GIVa+GIVb benches). No (1−P): χ is stored non-AS.
// bar_CHI blocks are keyed {ch_il_cc, ch_kj_cc}.
void InvPandyaTensor2n(const CcMatMap &bar_CHI, const Operator &Eta,
                       ChiOp2n &Chi_Op) {
  ModelSpace *ms = Chi_Op.ms;
  const int lambda = Eta.GetJRank();
  ms->PreCalculateSixJ();
  PreCalculateNineJTensorPandya(ms, lambda); // shape B, serial reserve
  const int n_cc = (int)ms->GetNumberTwoBodyChannels_CC();
  std::vector<const arma::mat *> bar_ptr((size_t)n_cc * n_cc, nullptr);
  for (const auto &kv : bar_CHI)
    bar_ptr[(size_t)kv.first[0] * n_cc + kv.first[1]] = &kv.second;
  std::vector<CcPair> keys;
  for (const auto &kv : Chi_Op.MatEl)
    keys.push_back(kv.first);
  const int nblocks = (int)keys.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ib = 0; ib < nblocks; ++ib) {
    const int ch_bra = keys[ib][0], ch_ket = keys[ib][1];
    arma::mat &Chi = *Chi_Op.GetBlock(ch_bra, ch_ket);
    const int J0 = ms->GetTwoBodyChannel(ch_bra).J;
    const int J1 = ms->GetTwoBodyChannel(ch_ket).J;
    const int nbras = Chi_Op.NumberKets(ch_bra);
    const int nkets = Chi_Op.NumberKets(ch_ket);
    const double pref01 = HatJ(J0) * HatJ(J1);
    for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
      auto ij = Chi_Op.Pq(ch_bra, ibra);
      if (ij[0] < 0)
        continue;
      const int i = ij[0], j = ij[1];
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = Chi_Op.Pq(ch_ket, iket);
        if (kl[0] < 0)
          continue;
        const int k = kl[0], l = kl[1];
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
        const int parity_il = (oi.l + ol.l) % 2;
        const int Tz_il = std::abs(oi.tz2 - ol.tz2) / 2;
        const int parity_kj = (ok.l + oj.l) % 2;
        const int Tz_kj = std::abs(ok.tz2 - oj.tz2) / 2;
        const int J2min = std::abs(oi.j2 - ol.j2) / 2;
        const int J2max = (oi.j2 + ol.j2) / 2;
        double zijkl = 0.0;
        for (int J2 = J2min; J2 <= J2max; ++J2) {
          const size_t ch_il = ms->GetTwoBodyChannelIndex(J2, parity_il, Tz_il);
          if ((int)ch_il >= n_cc)
            continue;
          TwoBodyChannel_CC &tcc_il = ms->GetTwoBodyChannel_CC(ch_il);
          if (tcc_il.GetNumberKets() < 1)
            continue;
          const int indx_il = (int)tcc_il.GetLocalIndex(i, l);
          const int J3min = std::max(std::abs(ok.j2 - oj.j2) / 2,
                                     std::abs(J2 - lambda));
          const int J3max = std::min((ok.j2 + oj.j2) / 2, J2 + lambda);
          for (int J3 = J3min; J3 <= J3max; ++J3) {
            const size_t ch_kj =
                ms->GetTwoBodyChannelIndex(J3, parity_kj, Tz_kj);
            if ((int)ch_kj >= n_cc)
              continue;
            const arma::mat *bar = bar_ptr[ch_il * n_cc + ch_kj];
            if (bar == nullptr)
              continue;
            TwoBodyChannel_CC &tcc_kj = ms->GetTwoBodyChannel_CC(ch_kj);
            const int indx_kj = (int)tcc_kj.GetLocalIndex(k, j);
            if (indx_il >= (int)bar->n_rows or indx_kj >= (int)bar->n_cols)
              continue;
            const double me1 = (*bar)(indx_il, indx_kj);
            if (std::abs(me1) < 1e-16)
              continue;
            const double ninej =
                ms->GetNineJ(lambda, J0, J1, J3, jj, jk, J2, ji, jl);
            if (std::abs(ninej) < 1e-16)
              continue;
            zijkl += ms->phase(J2) * HatJ(J2) * HatJ(J3) * ninej * me1;
          }
        }
        Chi(ibra, iket) = ms->phase(J0 + (oi.j2 + ok.j2) / 2 + lambda) *
                          pref01 * zijkl;
      }
    }
  }
}

// (ia)(lb) → (il)(ab): χ̄_{il,ab} from occupancy-AS χ_{ialb}.
// AMC scheme=((1,-3),(2,-4)) on labels (i,a,l,b) — lands on (il)|(ab) so the
// gold product Σ χ_ialb Ω_bjak is the same DGEMM mid as tts_ring.
// Locked: ≡ IMSRG adcb Pandya when T1 χ is AS (even π); ≡ m-gold on odd π T1
// (run/tmp givc_gold_dgemm_lock + learn/.../G4c_gold_il_ab_pandya.txt).
// Chi(i,a,l,b,J,J') is any χ^λ table in GetTBME_J units (ChiOp2n).
template <class ChiT>
double PandyaChiIalb(const ChiT &Chi, ModelSpace *ms, int lambda, int i,
                     int l, int a, int b, int Jil, int Jab) {
  if (not AngMom::Triangle(Jil, Jab, lambda))
    return 0.0;
  Orbit &oi = ms->GetOrbit(i);
  Orbit &ol = ms->GetOrbit(l);
  Orbit &oa = ms->GetOrbit(a);
  Orbit &ob = ms->GetOrbit(b);
  const double ji = oi.j2 * 0.5, jl = ol.j2 * 0.5, ja = oa.j2 * 0.5,
               jb = ob.j2 * 0.5;
  if (not AngMom::Triangle(ji, jl, (double)Jil) or
      not AngMom::Triangle(ja, jb, (double)Jab))
    return 0.0;
  double sm = 0.0;
  const int jiamin = std::abs(oi.j2 - oa.j2) / 2;
  const int jiamax = (oi.j2 + oa.j2) / 2;
  for (int Jia = jiamin; Jia <= jiamax; ++Jia) {
    const int jlbmin =
        std::max(std::abs(ol.j2 - ob.j2) / 2, std::abs(Jia - lambda));
    const int jlbmax = std::min((ol.j2 + ob.j2) / 2, Jia + lambda);
    for (int Jlb = jlbmin; Jlb <= jlbmax; ++Jlb) {
      // ninej{λ, Jil, Jab; Jlb, jl, jb; Jia, ji, ja}
      const double ninej =
          ms->GetNineJ(lambda, Jil, Jab, Jlb, jl, jb, Jia, ji, ja);
      if (std::abs(ninej) < 1e-14)
        continue;
      const double tbme = Chi(i, a, l, b, Jia, Jlb);
      sm += ms->phase(Jia + Jlb) * HatJ(Jia) * HatJ(Jlb) * ninej * tbme;
    }
  }
  const double pref =
      -ms->phase(Jil + Jab + lambda +
                 (oi.j2 + oa.j2 + ol.j2 + ob.j2) / 2) *
      HatJ(Jil) * HatJ(Jab);
  return pref * sm;
}

// (bj)(ak) → (ab)(kj): Ω̄_{ab,kj} from Ω_{bjak}.
// AMC scheme=((3,-1),(4,-2)) on labels (b,j,a,k). Twin of PandyaChiIalb.
double PandyaOmegaBjak(const Operator &Eta, int a, int b, int k, int j, int Jab,
                       int Jkj) {
  const int lambda = Eta.GetJRank();
  if (not AngMom::Triangle(Jab, Jkj, lambda))
    return 0.0;
  ModelSpace *ms = Eta.modelspace;
  Orbit &oa = ms->GetOrbit(a);
  Orbit &ob = ms->GetOrbit(b);
  Orbit &ok = ms->GetOrbit(k);
  Orbit &oj = ms->GetOrbit(j);
  const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5, jk = ok.j2 * 0.5,
               jj = oj.j2 * 0.5;
  if (not AngMom::Triangle(ja, jb, (double)Jab) or
      not AngMom::Triangle(jk, jj, (double)Jkj))
    return 0.0;
  double sm = 0.0;
  const int jbjmin = std::abs(ob.j2 - oj.j2) / 2;
  const int jbjmax = (ob.j2 + oj.j2) / 2;
  for (int Jbj = jbjmin; Jbj <= jbjmax; ++Jbj) {
    const int jakmin =
        std::max(std::abs(oa.j2 - ok.j2) / 2, std::abs(Jbj - lambda));
    const int jakmax = std::min((oa.j2 + ok.j2) / 2, Jbj + lambda);
    for (int Jak = jakmin; Jak <= jakmax; ++Jak) {
      // ninej{λ, Jab, Jkj; Jak, ja, jk; Jbj, jb, jj}
      const double ninej =
          ms->GetNineJ(lambda, Jab, Jkj, Jak, ja, jk, Jbj, jb, jj);
      if (std::abs(ninej) < 1e-14)
        continue;
      const double tbme = Eta.TwoBody.GetTBME_J(Jbj, Jak, b, j, a, k);
      sm += ms->phase(Jbj + Jak) * HatJ(Jbj) * HatJ(Jak) * ninej * tbme;
    }
  }
  return -ms->phase(lambda) * HatJ(Jab) * HatJ(Jkj) * sm;
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

/// Scalar [Ω⊗Ω]^0 intermediates of Pandya topology (χ^η, χ^γ), AMC Path B:
/// same-label tensor Pandya Ω̄ (9j) → occupancy ⊙ Ω̄ → DGEMM in CC blocks →
/// optional scalar inverse Pandya to 2n ordinary blocks (InvChiEtaRed).
///   χ̄^{Jb}(il;kj) = (−1)^{Jb}/Ĵb Σ_{Jk} (−1)^{Jk+λ}/λ̂ Σ_{pq}
///                     Ω̄^{Jb Jk}(il;pq) · w(p,q,k,j) Ω̄^{Jk Jb}(pq;kj)
/// The diagrams differ only in the occupancy weight w on the (pq)×(kj) legs:
///   χ^η  w = n̄_q n_p n̄_k + n_q n̄_p n_k                   GIIIa / GIIIb
///   χ^γ  w = n_p n̄_q n_j n̄_k − n̄_p n_q n̄_j n_k           f^{III_a}
/// χ is scalar and not AS. Store 2n ordinary (not n×n GetMatrix).
template <class WeightF>
void BuildChiPandyaPathB(const Operator &Eta, Operator &Z, WeightF weight,
                         std::deque<arma::mat> &Chi2,
                         std::deque<arma::mat> *barCHI_out, bool fill_chi2) {
  const int lambda = Eta.GetJRank();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();
  PreCalculateNineJTensorPandya(Z.modelspace, lambda); // amc_bar_omega 9j
  const int nch_eta = Z.modelspace->GetNumberTwoBodyChannels();
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  const double hat_lambda = std::sqrt(2.0 * lambda + 1.0);
  double t_start = omp_get_wtime();

  // 9j Pandya of Ω̄ only for ch_b ≤ ch_k. Partner block is reduced transpose
  // hΩ (−1)^{J0−J1} Ω̄^T (LESSONS.md). Occ χ̄ is applied after, on that layout.
  //
  //   Ω̄^{Jb Jk}(i,l; k,j) = −(−1)^{Jb+(ji+jk)/2+λ} Ĵb Ĵk
  //                         Σ_{J2 J3} (−1)^{J2} Ĵ2 Ĵ3
  //                         { λ  Jb Jk }  Ω^{J2 J3}(i,j; k,l)
  //                         { J3 jl jk }
  //                         { J2 ji jj }
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
  std::vector<int> pair_of((size_t)n_cc * n_cc, -1);
  for (int ic = 0; ic < nc; ++ic) {
    TwoBodyChannel_CC &tb =
        Z.modelspace->GetTwoBodyChannel_CC(cc_canon[ic][0]);
    TwoBodyChannel_CC &tk =
        Z.modelspace->GetTwoBodyChannel_CC(cc_canon[ic][1]);
    tmpOm[ic] = arma::mat(2 * tb.GetNumberKets(), 2 * tk.GetNumberKets(),
                          arma::fill::zeros);
    pair_of[(size_t)cc_canon[ic][0] * n_cc + cc_canon[ic][1]] = ic;
  }
  std::vector<size_t> om_bra, om_ket;
  TensorChannelPairs(Eta, Z.modelspace, om_bra, om_ket);
  const int n_om = (int)om_bra.size();
  std::vector<arma::mat> Om2n((size_t)n_om);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < n_om; ++ip)
    Om2n[ip] = Omega2n(Eta, Z.modelspace, (int)om_bra[ip], (int)om_ket[ip]);
#pragma omp parallel for schedule(dynamic, 1)
  for (int iom = 0; iom < n_om; ++iom) {
    TwoBodyChannel &tbc2 = Z.modelspace->GetTwoBodyChannel((int)om_bra[iom]);
    TwoBodyChannel &tbc3 = Z.modelspace->GetTwoBodyChannel((int)om_ket[iom]);
    const int J2 = tbc2.J, J3 = tbc3.J;
    const arma::mat &Om = Om2n[iom];
    const int n2 = tbc2.GetNumberKets(), n3 = tbc3.GetNumberKets();
    const double hats23 = HatJ(J2) * HatJ(J3);
    const double phJ2 = Z.modelspace->phase(J2);
    for (int ibra = 0; ibra < 2 * n2; ++ibra) {
      auto ij = PqTB(tbc2, ibra, n2);
      if (ij[0] < 0)
        continue;
      const int i = ij[0], j = ij[1];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
      for (int iket = 0; iket < 2 * n3; ++iket) {
        auto kl = PqTB(tbc3, iket, n3);
        if (kl[0] < 0)
          continue;
        const double om = Om(ibra, iket);
        if (std::abs(om) < 1e-16)
          continue;
        const int k = kl[0], l = kl[1];
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
        const int parity_il = (oi.l + ol.l) % 2;
        const int Tz_il = std::abs(oi.tz2 - ol.tz2) / 2;
        const int parity_kj = (ok.l + oj.l) % 2;
        const int Tz_kj = std::abs(ok.tz2 - oj.tz2) / 2;
        const int Jbmin = std::abs(oi.j2 - ol.j2) / 2;
        const int Jbmax = (oi.j2 + ol.j2) / 2;
        for (int Jb = Jbmin; Jb <= Jbmax; ++Jb) {
          const size_t ch_il =
              Z.modelspace->GetTwoBodyChannelIndex(Jb, parity_il, Tz_il);
          if ((int)ch_il >= n_cc)
            continue;
          const int row =
              IndexCC(Z.modelspace->GetTwoBodyChannel_CC(ch_il), i, l);
          if (row < 0)
            continue;
          const int Jkmin =
              std::max(std::abs(ok.j2 - oj.j2) / 2, std::abs(Jb - lambda));
          const int Jkmax = std::min((ok.j2 + oj.j2) / 2, Jb + lambda);
          for (int Jk = Jkmin; Jk <= Jkmax; ++Jk) {
            const size_t ch_kj =
                Z.modelspace->GetTwoBodyChannelIndex(Jk, parity_kj, Tz_kj);
            if ((int)ch_kj >= n_cc)
              continue;
            if ((int)ch_il > (int)ch_kj)
              continue;
            const int ip = pair_of[(size_t)ch_il * n_cc + ch_kj];
            if (ip < 0)
              continue;
            const int col =
                IndexCC(Z.modelspace->GetTwoBodyChannel_CC(ch_kj), k, j);
            if (col < 0)
              continue;
            const double n9 = Z.modelspace->GetNineJ(
                (double)lambda, (double)Jb, (double)Jk, (double)J3, jl, jk,
                (double)J2, ji, jj);
            if (std::abs(n9) < 1e-16)
              continue;
            const double add =
                -Z.modelspace->phase(Jb + (oi.j2 + ok.j2) / 2 + lambda) *
                HatJ(Jb) * HatJ(Jk) * phJ2 * hats23 * n9 * om;
            double *dst = &tmpOm[ip](row, col);
#pragma omp atomic
            *dst += add;
          }
        }
      }
    }
  }
  Om2n.clear();
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

  Z.profiler.timer["BuildChiEtaPathB_pandya"] += omp_get_wtime() - t_start;
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
      arma::mat Rocc = itR->second; // Ω̄^{Jk Jb}(pq;kj), then ⊙ w(p,q,k,j)
      for (int ibra = 0; ibra < 2 * nk; ++ibra) {
        auto pq = PqCC(tk, ibra, nk);
        if (pq[0] < 0)
          continue;
        for (int iket = 0; iket < 2 * nb; ++iket) {
          auto kj = PqCC(tb, iket, nb);
          if (kj[0] < 0)
            continue;
          Rocc(ibra, iket) *= weight(pq[0], pq[1], kj[0], kj[1]);
        }
      }
      const double wR = Z.modelspace->phase(Jk + lambda) / hat_lambda;
      barCHI[ch_b] += (wL * wR) * (itL->second * Rocc);
    }
  }
  bar_Omega.clear();

  Z.profiler.timer["BuildChiEtaPathB_dgemm_chi"] += omp_get_wtime() - t_mid;
  t_mid = omp_get_wtime();

  if (fill_chi2) {
    Chi2.assign(nch_eta, arma::mat());
    for (int ch = 0; ch < nch_eta; ++ch) {
      const int nK = Z.modelspace->GetTwoBodyChannel(ch).GetNumberKets();
      Chi2[ch] = arma::mat(2 * nK, 2 * nK, arma::fill::zeros);
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc) {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      const int nkets_cc = tbc_cc.GetNumberKets();
      if (nkets_cc < 1 or barCHI[ch_cc].n_rows < 1)
        continue;
      const int Jp = tbc_cc.J;
      const double hatJp = HatJ(Jp);
      const arma::mat &bar = barCHI[ch_cc];
      for (int ibra = 0; ibra < 2 * nkets_cc; ++ibra) {
        auto il = PqCC(tbc_cc, ibra, nkets_cc);
        if (il[0] < 0)
          continue;
        const int i = il[0], l = il[1];
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const double ji = oi.j2 * 0.5, jl = ol.j2 * 0.5;
        for (int iket = 0; iket < 2 * nkets_cc; ++iket) {
          auto kj = PqCC(tbc_cc, iket, nkets_cc);
          if (kj[0] < 0)
            continue;
          const double me = bar(ibra, iket);
          if (std::abs(me) < 1e-16)
            continue;
          const int k = kj[0], j = kj[1];
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &oj = Z.modelspace->GetOrbit(j);
          const double jk = ok.j2 * 0.5, jj = oj.j2 * 0.5;
          if ((oi.l + oj.l) % 2 != (ok.l + ol.l) % 2)
            continue;
          if (oi.tz2 + oj.tz2 != ok.tz2 + ol.tz2)
            continue;
          const int parity_ij = (oi.l + oj.l) % 2;
          const int Tz_ij = (oi.tz2 + oj.tz2) / 2;
          const int J0min =
              std::max(std::abs(oi.j2 - oj.j2), std::abs(ok.j2 - ol.j2)) / 2;
          const int J0max = std::min(oi.j2 + oj.j2, ok.j2 + ol.j2) / 2;
          for (int J0 = J0min; J0 <= J0max; ++J0) {
            const double sixj =
                Z.modelspace->GetSixJ(jl, jk, (double)J0, jj, ji, Jp);
            if (std::abs(sixj) < 1e-16)
              continue;
            const size_t ch0 = Z.modelspace->GetTwoBodyChannelIndex(
                J0, parity_ij, Tz_ij);
            if ((int)ch0 >= nch_eta)
              continue;
            TwoBodyChannel &tbc0 = Z.modelspace->GetTwoBodyChannel((int)ch0);
            const int row = Index2n(tbc0, i, j);
            const int col = Index2n(tbc0, k, l);
            if (row < 0 or col < 0)
              continue;
            const double add = HatJ(J0) * hatJp * sixj * me;
            double *dst = &Chi2[(int)ch0](row, col);
#pragma omp atomic
            *dst += add;
          }
        }
      }
    }
  }

  Z.profiler.timer["BuildChiEtaPathB_inv"] += omp_get_wtime() - t_mid;

  if (barCHI_out)
    *barCHI_out = std::move(barCHI);
  else {
    for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc)
      barCHI[ch_cc].clear();
    barCHI.clear();
  }
}

/// χ^η (GIIIa / GIIIb). Gold: test_chi_eta_mscheme.py.
void BuildChiEtaPathB(const Operator &Eta, Operator &Z,
                      std::deque<arma::mat> &Chi2,
                      std::deque<arma::mat> *barCHI_out, bool fill_chi2) {
  ModelSpace *ms = Z.modelspace;
  auto w_eta = [ms](int p, int q, int k, int /*j*/) {
    const double np = ms->GetOrbit(p).occ, nq = ms->GetOrbit(q).occ;
    const double nk = ms->GetOrbit(k).occ;
    return (1.0 - nq) * np * (1.0 - nk) + nq * (1.0 - np) * nk;
  };
  BuildChiPandyaPathB(Eta, Z, w_eta, Chi2, barCHI_out, fill_chi2);
}

/// χ^γ (f^{III_a}), unreduced (AMC direct / Ĵ0), in the no-Pauli ChiOp2n
/// (ch,ch) blocks. χ^γ is not AS: χ(pp;ab)^{J0 odd} ≠ 0 and it does
/// contribute to f^{III_a} through Γ(ab;qp) when q ≠ p shares (j,l,tz) with p
/// (first at emax=2: 0s1/2–1s1/2). The Pauli-filtered ordinary 2n layout has
/// no slot for (p,p) at odd J0, so ChiOp2n is required here.
/// Gold: run/test_chi_gamma_pathB_amc.py (m ≡ AMC direct ≡ Path B).
void BuildChiGammaPathB(const Operator &Eta, Operator &Z, ChiOp2n &CHI_gamma) {
  ModelSpace *ms = Z.modelspace;
  auto w_gamma = [ms](int p, int q, int k, int j) {
    const double np = ms->GetOrbit(p).occ, nq = ms->GetOrbit(q).occ;
    const double nk = ms->GetOrbit(k).occ, nj = ms->GetOrbit(j).occ;
    return np * (1.0 - nq) * nj * (1.0 - nk) -
           (1.0 - np) * nq * (1.0 - nj) * nk;
  };
  std::deque<arma::mat> Chi2_unused;
  std::deque<arma::mat> barCHI;
  BuildChiPandyaPathB(Eta, Z, w_gamma, Chi2_unused, &barCHI, false);

  const int nch = ms->GetNumberTwoBodyChannels();
  CHI_gamma.Init(ms);
  std::vector<size_t> ch_list;
  for (int ch = 0; ch < nch; ++ch)
    if (CHI_gamma.NumberKets(ch) > 0)
      ch_list.push_back((size_t)ch);
  CHI_gamma.Allocate(ch_list, ch_list);
  const int nlist = (int)ch_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ic = 0; ic < nlist; ++ic) {
    const int ch = (int)ch_list[ic];
    const int J0 = ms->GetTwoBodyChannel(ch).J;
    const double hatJ0_inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
    const int nK = CHI_gamma.NumberKets(ch);
    arma::mat &M = *CHI_gamma.GetBlock(ch, ch);
    for (int ibra = 0; ibra < 2 * nK; ++ibra) {
      auto ij = CHI_gamma.Pq(ch, ibra);
      if (ij[0] < 0)
        continue;
      for (int iket = 0; iket < 2 * nK; ++iket) {
        auto kl = CHI_gamma.Pq(ch, iket);
        if (kl[0] < 0)
          continue;
        // InvChiEtaRed is Ĵ0·χ_unred
        M(ibra, iket) =
            hatJ0_inv * InvChiEtaRed(Z, barCHI, ij[0], ij[1], kl[0], kl[1], J0);
      }
    }
  }
}

} // namespace

/// f^{III_a} via χ^γ (any λ): AMC Path B χ^γ (BuildChiGammaPathB: tensor
/// Pandya Ω̄ → occ ⊙ Ω̄ → CC DGEMM → inverse Pandya), then the ordinary-channel
/// ladder traced on the spectator c:
///   z_pq = 1/ĵ_p² Σ_{c J0} (2J0+1) [ (Γ χ)^{J0}(cp;cq) − (χ Γ)^{J0}(pc;qc) ]
/// Γ and χ on ordered pairs in the no-Pauli ChiOp2n layout (GetTBME_J units):
/// χ^γ is not AS, so the two products are kept apart and the (c,c) slot at
/// odd J0 is kept (Γ(cc;ab)^{odd}=0 but χ(cc;ab)^{odd}≠0 feeds Γ(ab;qc)).
/// Gold: run/test_chi_gamma_pathB_amc.py; bench TERMS=fIIIa vs Mscheme_fact_fIIIa.
static void comm223_231_fIIIa_leftover_dgemm(const Operator &Eta,
                                            const Operator &Gamma,
                                            Operator &Z) {
  double t_start = omp_get_wtime();
  ModelSpace *ms = Z.modelspace;
  ms->PreCalculateSixJ();
  const int hZ = Gamma.IsHermitian() ? 1 : -1;
  const int nch = ms->GetNumberTwoBodyChannels();

  ChiOp2n CHI_gamma; // χ^γ, no-Pauli 2n (ch,ch) blocks, unreduced
  BuildChiGammaPathB(Eta, Z, CHI_gamma);
  double t_mid = omp_get_wtime();
  Z.profiler.timer["_fIIIa_chi"] += t_mid - t_start;

  // (2J0+1) Γ·χ and (2J0+1) χ·Γ per ordinary channel, ChiOp2n layout.
  std::deque<arma::mat> GamChi(nch), ChiGam(nch);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    const int nK = CHI_gamma.NumberKets(ch);
    if (nK < 1)
      continue;
    const int J0 = ms->GetTwoBodyChannel(ch).J;
    const arma::mat &Chi = *CHI_gamma.GetBlock(ch, ch);
    arma::mat Gam2(2 * nK, 2 * nK, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nK; ++ibra) {
      auto ij = CHI_gamma.Pq(ch, ibra);
      if (ij[0] < 0)
        continue;
      for (int iket = 0; iket < 2 * nK; ++iket) {
        auto kl = CHI_gamma.Pq(ch, iket);
        if (kl[0] < 0)
          continue;
        Gam2(ibra, iket) = // 0 for (p,p) at odd J0
            Gamma.TwoBody.GetTBME_J(J0, J0, ij[0], ij[1], kl[0], kl[1]);
      }
    }
    GamChi[ch] = (2.0 * J0 + 1.0) * (Gam2 * Chi);
    ChiGam[ch] = (2.0 * J0 + 1.0) * (Chi * Gam2);
  }
  CHI_gamma.MatEl.clear();
  Z.profiler.timer["_fIIIa_dgemm"] += omp_get_wtime() - t_mid;
  t_mid = omp_get_wtime();

  // Trace on c.
  std::vector<index_t> allorb(ms->all_orbits.begin(), ms->all_orbits.end());
  const int norb = (int)allorb.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < norb; ++ip) {
    const int p = (int)allorb[ip];
    Orbit &op = ms->GetOrbit(p);
    for (int iq = 0; iq <= ip; ++iq) {
      const int q = (int)allorb[iq];
      Orbit &oq = ms->GetOrbit(q);
      if (oq.j2 != op.j2 or oq.l != op.l or oq.tz2 != op.tz2)
        continue;
      double zpq = 0.0;
      for (auto c : allorb) {
        Orbit &oc = ms->GetOrbit(c);
        const int J0min = std::abs(oc.j2 - op.j2) / 2;
        const int J0max = (oc.j2 + op.j2) / 2;
        for (int J0 = J0min; J0 <= J0max; ++J0) {
          const int ch = CHI_gamma.ChannelIndex(J0, (int)c, p);
          if (ch < 0 or CHI_gamma.NumberKets(ch) < 1)
            continue;
          const int i_cp = CHI_gamma.Index(ch, (int)c, p);
          const int i_cq = CHI_gamma.Index(ch, (int)c, q);
          const int i_pc = CHI_gamma.Index(ch, p, (int)c);
          const int i_qc = CHI_gamma.Index(ch, q, (int)c);
          if (i_cp < 0 or i_cq < 0 or i_pc < 0 or i_qc < 0)
            continue;
          zpq += GamChi[ch](i_cp, i_cq) - ChiGam[ch](i_pc, i_qc);
        }
      }
      Z.OneBody(p, q) += zpq / (op.j2 + 1.0);
      if (p != q)
        Z.OneBody(q, p) += hZ * zpq / (op.j2 + 1.0);
    }
  }
  Z.profiler.timer["_fIIIa_fold"] += omp_get_wtime() - t_mid;
  Z.profiler.timer["comm223_231_tts_fIIIa"] += omp_get_wtime() - t_start;
}

void comm223_232_GIIIa(const Operator &Eta_in, const Operator &Gamma,
                       Operator &Z) {
  // Γ^{III_a}: AMC Path B χ^η (same-label Pandya → occ DGEMM → inv)
  // then ordinary-channel Chi_AS×Γ ladder. χ is scalar and not AS.
  // Gold: test_chi_eta_mscheme.py + test_GIIIa_ladder_mscheme.py
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
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

    const arma::mat Gam2 =
        Omega2n(Gamma, Z.modelspace, (int)ch_bra, (int)ch_ket);

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

  Z.profiler.timer[std::string(__func__) + "_ladder"] +=
      omp_get_wtime() - t_mid;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// Gamma^III_b — IIb and IId. FDC.cc comm223_232_chi2b "Diagram IIb and IId"
/// with the tensor χ̄^η (BuildChiEtaPathB) in place of bar_Eta·nnnbar_Eta.
/// χ^η = [Ω⊗Ω]^0 is scalar for every λ, so after χ̄^η the chain is the scalar
/// FDC.cc one with the same names:
///   barCHI_III     χ̄^η(il;kj)           CC 2n×2n   (BuildChiEtaPathB / Ĵ_cc)
///   barCHI_III_RC  RC[χ̄(bc;ad)+χ̄(ad;bc)]            FDC.cc ~1799
///   CHI_III_final  bar_Gamma · barCHI_III_RC          FDC.cc ~2752
///   Z              (1−Pij)(1−Pkl) inverse Pandya from (j l̄; i k̄)  FDC.cc ~2773
/// Gold chain (learn/amc_tts/factored_GIIIb/NOTES.md): Path B pack → RC → Γ̄·RC
/// → Inv → (1−P)² ≡ fold ≡ 4-index m, λ=0…4 (run/test_G3b_pathB_pack_mscheme.py);
/// bench run/bench_eths_pathB_vs_mscheme.py TERMS=GIIIb vs Mscheme_fact_GIIIb.
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIIIb(const Operator &Eta_in, const Operator &Gamma,
                       Operator &Z) {
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();
  ModelSpace *ms = Z.modelspace;
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const int n_nonzero = ms->GetNumberTwoBodyChannels_CC();
  auto &Z2 = Z.TwoBody;

  // χ̄^η(il;kj) in CC channels. BuildChiEtaPathB keeps Ĵ_cc·χ̄ (GIIIa's
  // InvChiEtaRed absorbs it); the scalar RC wants the plain Pandya χ̄.
  std::deque<arma::mat> Chi2_unused;
  std::deque<arma::mat> barCHI_III;
  BuildChiEtaPathB(Eta, Z, Chi2_unused, &barCHI_III, false);
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
    if (tbc_cc.GetNumberKets() < 1)
      continue;
    barCHI_III[ch_cc] /= std::sqrt(2.0 * tbc_cc.J + 1.0);
  }
  double t_internal = omp_get_wtime();
  Z.profiler.timer["_GIIIb_chi"] += t_internal - t_start;

  // Recoupling, Pandya → Pandya (FDC.cc ~1799):
  //  \dbar{χ}^J_{ab,cd} = −Σ_J' (2J'+1) (−)^{jb+jc+J'} { ja jb J  }
  //                        (χ̄^J'_{bc,ad} + χ̄^J'_{ad,bc})  { jc jd J' }
  std::deque<arma::mat> barCHI_III_RC(n_nonzero);
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    const int nKets_cc = ms->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
    barCHI_III_RC[ch_cc] =
        arma::mat(2 * nKets_cc, 2 * nKets_cc, arma::fill::zeros);
  }
  // Scatter χ̄^{J'}(bc;ad)+χ̄^{J'}(ad;bc) → dest (ab;cd). Same 6j as the
  // dest loop; walk each source CC block once.
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_old = 0; ch_old < n_nonzero; ++ch_old) {
    if (barCHI_III[ch_old].n_rows < 1)
      continue;
    TwoBodyChannel_CC &tbc_old = ms->GetTwoBodyChannel_CC(ch_old);
    const int nold = tbc_old.GetNumberKets();
    if (nold < 1)
      continue;
    const int Jp = tbc_old.J;
    const double hat2 = 2 * Jp + 1;
    const arma::mat &bar = barCHI_III[ch_old];
    for (int ibra = 0; ibra < 2 * nold; ++ibra) {
      auto bc = PqCC(tbc_old, ibra, nold);
      if (bc[0] < 0)
        continue;
      const int b = bc[0], c = bc[1];
      Orbit &ob = ms->GetOrbit(b);
      Orbit &oc = ms->GetOrbit(c);
      const double jb = ob.j2 * 0.5, jc = oc.j2 * 0.5;
      const double ph = ms->phase((ob.j2 + oc.j2) / 2 + Jp);
      for (int iket = 0; iket < 2 * nold; ++iket) {
        auto ad = PqCC(tbc_old, iket, nold);
        if (ad[0] < 0)
          continue;
        const double chi_sum = bar(ibra, iket) + bar(iket, ibra);
        if (std::abs(chi_sum) < 1e-16)
          continue;
        const int a = ad[0], d = ad[1];
        Orbit &oa = ms->GetOrbit(a);
        Orbit &od = ms->GetOrbit(d);
        const double ja = oa.j2 * 0.5, jd = od.j2 * 0.5;
        const int parity_ab = (oa.l + ob.l) % 2;
        const int Tz_ab = std::abs(oa.tz2 - ob.tz2) / 2;
        if ((oc.l + od.l) % 2 != parity_ab)
          continue;
        if (std::abs(oc.tz2 - od.tz2) / 2 != Tz_ab)
          continue;
        const int Jmin =
            std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
        const int Jmax = std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
        for (int J_cc = Jmin; J_cc <= Jmax; ++J_cc) {
          const double sixj = ms->GetSixJ(ja, jb, J_cc, jc, jd, Jp);
          if (std::abs(sixj) < 1e-8)
            continue;
          const int ch_cc =
              ms->GetTwoBodyChannelIndex(J_cc, parity_ab, Tz_ab);
          if (ch_cc < 0 or ch_cc >= n_nonzero)
            continue;
          TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
          const int row = IndexCC(tbc_cc, a, b);
          const int col = IndexCC(tbc_cc, c, d);
          if (row < 0 or col < 0)
            continue;
          const double add = -ph * hat2 * sixj * chi_sum;
          double *dst = &barCHI_III_RC[ch_cc](row, col);
#pragma omp atomic
          *dst += add;
        }
      }
    }
  }
  barCHI_III.clear();
  Z.profiler.timer["_GIIIb_RC"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // CHI_III_final = bar_Gamma · barCHI_III_RC (FDC.cc ~2752)
  const std::deque<arma::mat> &bar_Gamma =
      CachedBarGammaScalarCC(Gamma, Z, hGamma);
  std::deque<arma::mat> CHI_III_final(n_nonzero);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc) {
    if (ms->GetTwoBodyChannel_CC(ch_cc).GetNumberKets() < 1)
      continue;
    CHI_III_final[ch_cc] = bar_Gamma[ch_cc] * barCHI_III_RC[ch_cc];
  }
  barCHI_III_RC.clear();
  Z.profiler.timer["_GIIIb_dgemm"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Inverse Pandya (FDC.cc ~2773):
  //  Z^J_ijkl = (1−Pij)(1−Pkl) (−)^{J+ji+jj} Σ_J' −(2J'+1) (−)^{J'+ji+jk}
  //             { j i J  }  \bar{X}^J'_{j l̄, i k̄}
  //             { k l J' }
  // Every (ibra,iket) is written: one tensor diagram on its own need not be
  // Hermitian, so no ketmin / AddToTBME mirroring.
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch) {
    const size_t ch_bra = ch_bra_list[ch];
    const size_t ch_ket = ch_ket_list[ch];
    TwoBodyChannel &tbc_bra = ms->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = ms->GetTwoBodyChannel(ch_ket);
    const int nbras = tbc_bra.GetNumberKets();
    const int nKets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nKets < 1)
      continue;
    const int J0 = tbc_bra.J;
    if (tbc_ket.J != J0)
      continue;
    arma::mat &Zmat = Z2.GetMatrix(ch_bra, ch_ket);
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const int ji = oi.j2, jj = oj.j2;
      const int phaseFactor = ms->phase(J0 + (ji + jj) / 2);
      for (int iket = 0; iket < nKets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const int jk = ok.j2, jl = ol.j2;
        double commijkl = 0, commjikl = 0, commijlk = 0, commjilk = 0;

        // jikl, direct term        -->  jl  ki
        // ijlk, exchange ij and kl -->  lj  ik
        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          const double sixj1 =
              ms->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          const int ch_cc =
              ms->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          if (ch_cc < 0 or ch_cc >= n_nonzero)
            continue;
          TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
          const int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_jl < 0 or indx_ik < 0)
            continue;
          indx_jl += (j > l ? nkets_cc : 0);
          indx_ik += (i > k ? nkets_cc : 0);
          const int phase1 = ms->phase(Jprime + (ji + jk) / 2);
          const double me1 = CHI_III_final[ch_cc](indx_jl, indx_ik);
          commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * me1;
          const int phase2 = ms->phase(Jprime + (jj + jl) / 2);
          const double me2 = CHI_III_final[ch_cc](indx_ik, indx_jl);
          commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        // ijkl,  exchange i and j -->  il  kj
        // jilk,  exchange k and l -->  jk li
        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime) {
          const double sixj1 =
              ms->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj1) < 1e-8)
            continue;
          const int ch_cc =
              ms->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          if (ch_cc < 0 or ch_cc >= n_nonzero)
            continue;
          TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
          const int nkets_cc = tbc_cc.GetNumberKets();
          if (nkets_cc < 1)
            continue;
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_jk < 0)
            continue;
          indx_il += (i > l ? nkets_cc : 0);
          indx_jk += (j > k ? nkets_cc : 0);
          const int phase1 = ms->phase(Jprime + (ji + jl) / 2);
          const double me1 = CHI_III_final[ch_cc](indx_jk, indx_il);
          commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * me1;
          const int phase2 = ms->phase(Jprime + (jj + jk) / 2);
          const double me2 = CHI_III_final[ch_cc](indx_il, indx_jk);
          commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
        }

        double zijkl = commjikl - ms->phase((ji + jj) / 2 - J0) * commijkl;
        zijkl += -ms->phase((jl + jk) / 2 - J0) * commjilk +
                 ms->phase((jk + jl + ji + jj) / 2) * commijlk;
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Zmat(ibra, iket) += phaseFactor * zijkl;
      }
    }
  }
  CHI_III_final.clear();

  Z.profiler.timer["_GIIIb_fold"] += omp_get_wtime() - t_internal;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

namespace {
// χ^θ (G3c; FDC.cc CHI_IV): scalar [Ω⊗Ω]^(0) by ordinary-channel DGEMM in the
// {ch_bra, ch_ket} blocks of TensorChannelPairs (Ω^{J0 J1} couples the two).
//   T^{J0}_{pp|hh}(ij;kl) = Σ_{J1} (−1)^{J0+J1+λ} λ̂^{-1}
//                            Ω^{J0J1}(ij;ab) w_{pp|hh}(ab) Ω^{J1J0}(ab;kl)
//   CHI_IV^{J0}(ij;kl)   = [ (n_k + n_j) T_pp + (n̄_k + n̄_j) T_hh ] / Ĵ0²
// (χ_k weight on ket k plus χ_j weight on bra j; same T.) Rows/cols are the
// ordered pairs of PqTB (2n, no Pauli): χ^θ has occupancy on j and k, so
// χ(ji;kl) ≠ −phase χ(ij;kl) when n_i ≠ n_j and TwoBodyME cannot hold it.
void FillChiThetaG3c_DGEMM(const Operator &Eta, Operator &Z, int lambda,
                           std::deque<arma::mat> &CHI_IV) {
  ModelSpace *ms = Z.modelspace;
  const int nch = ms->GetNumberTwoBodyChannels();
  (void)lambda; // λ enters through Eta.GetJRank() in FillOmegaOccOmega2n

  std::deque<arma::mat> T_pp_all, T_hh_all;
  FillOmegaOccOmega2n(Eta, ms, T_pp_all, T_hh_all);

  CHI_IV.assign(nch, arma::mat());
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch0 = 0; ch0 < nch; ++ch0) {
    TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
    const int n0 = tbc0.GetNumberKets();
    const int J0 = tbc0.J;
    CHI_IV[ch0] = arma::mat(2 * n0, 2 * n0, arma::fill::zeros);
    if (n0 < 1)
      continue;
    const arma::mat &T_pp = T_pp_all[ch0];
    const arma::mat &T_hh = T_hh_all[ch0];
    const double inv_hat2 = 1.0 / (2.0 * J0 + 1.0);
    for (int ibra = 0; ibra < 2 * n0; ++ibra) {
      auto ij = PqTB(tbc0, ibra, n0);
      if (ij[0] < 0)
        continue;
      const double n_j = ms->GetOrbit(ij[1]).occ;
      for (int iket = 0; iket < 2 * n0; ++iket) {
        auto kl = PqTB(tbc0, iket, n0);
        if (kl[0] < 0)
          continue;
        const double n_k = ms->GetOrbit(kl[0]).occ;
        CHI_IV[ch0](ibra, iket) =
            ((n_k + n_j) * T_pp(ibra, iket) +
             (2.0 - n_k - n_j) * T_hh(ibra, iket)) *
            inv_hat2;
      }
    }
  }
}
} // namespace

void comm223_232_GIIIc(const Operator &Eta_in, const Operator &Gamma,
                       Operator &Z) {
  // χ^θ = T×T→S (FDC.cc CHI_IV / IIe+IIf): ordinary-channel DGEMM χ^θ →
  // Pandya χ̄^θ (bar_CHI_IV) → DGEMM with bar_Gamma → inverse Pandya.
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
  const int lambda = Eta.GetJRank();
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;

  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z.TwoBody.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
  const int nch_eta = Z.modelspace->GetNumberTwoBodyChannels();
  const int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
  auto &Z2 = Z.TwoBody;

  // CHI_IV: χ^θ in ordinary channels, 2n ordered-pair layout.
  std::deque<arma::mat> CHI_IV;
  FillChiThetaG3c_DGEMM(Eta, Z, lambda, CHI_IV);
  if (Commutator::verbose)
    Z.profiler.timer["_GIIIc_chi"] += omp_get_wtime() - t_start;

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
  double norm2b = -1.0; // guards in-place updates of the same Operator
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
  const double norm2b = Gamma.TwoBodyNorm();
  if (g_bar_gamma_cc.gamma == &Gamma && g_bar_gamma_cc.nch == nch &&
      g_bar_gamma_cc.norm2b == norm2b &&
      (int)g_bar_gamma_cc.bar.size() == nch)
    return g_bar_gamma_cc.bar;
  FillBarGammaScalarCC(Gamma, Z, hGamma, g_bar_gamma_cc.bar);
  g_bar_gamma_cc.gamma = &Gamma;
  g_bar_gamma_cc.nch = nch;
  g_bar_gamma_cc.norm2b = norm2b;
  return g_bar_gamma_cc.bar;
}

// FDC.cc ~1243–1393 for a tensor Ω. One ket loop writes bar_Eta (Ω̄, 9j
// DoTensorPandya adcb) and the occupancy copies nnnbar_Eta / nnnbar_Eta_d
// (not a second Pandya); bar_Gamma is the scalar 6j Pandya of Γ (cached).
// Then χ̄ by DGEMM in the same CC block, exactly as FDC.cc ~1385.
//
//   Ω̄^{Jb Jk}(ab;cd) = −Σ_{J1 J2} Ĵ1 Ĵ2 Ĵb Ĵk (−1)^{jb+jd+Jk+J2}
//                        { ja jd J1 }
//                        { jb jc J2 }  Ω^{J1 J2}(ad;cb)
//                        { Jb Jk λ  }
//   partner: Ω̄^{Jk Jb}(cd;ab) = hΩ (−1)^{Jb−Jk} Ω̄^{Jb Jk}(ab;cd)ᵀ
void FillPandyaBars232(const Operator &Eta, const Operator &Gamma, Operator &Z,
                       PandyaBars232 &B) {
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();
  const int lambda = Eta.GetJRank();
  PreCalculateNineJTensorPandya(Z.modelspace, lambda); // shape A (bar_Eta)
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();

  B.bar_Gamma = &CachedBarGammaScalarCC(Gamma, Z, hGamma);
  const std::deque<arma::mat> &bar_Gamma = *B.bar_Gamma;
  B.cc_pairs.clear();
  for (int ch_bra_cc = 0; ch_bra_cc < n_cc; ++ch_bra_cc) {
    TwoBodyChannel_CC &tbc_bra_cc =
        Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
    if (tbc_bra_cc.GetNumberKets() < 1)
      continue;
    for (int ch_ket_cc = 0; ch_ket_cc < n_cc; ++ch_ket_cc) {
      TwoBodyChannel_CC &tbc_ket_cc =
          Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
      if (tbc_ket_cc.GetNumberKets() < 1)
        continue;
      // Odd-π Ω couples opposite-π CC channels (BuildChiEtaPathB / GIVc).
      if ((tbc_bra_cc.parity + tbc_ket_cc.parity + Eta.GetParity()) % 2 != 0)
        continue;
      if (not CcTzCouples(tbc_bra_cc.Tz, tbc_ket_cc.Tz, Eta.GetTRank()))
        continue;
      if (not AngMom::Triangle(tbc_bra_cc.J, tbc_ket_cc.J, lambda))
        continue;
      B.cc_pairs.push_back({ch_bra_cc, ch_ket_cc});
    }
  }

  const int np = (int)B.cc_pairs.size();
  std::vector<arma::mat> bar_Eta(np), nnnbar_Eta(np), nnnbar_Eta_d(np);
  std::vector<int> pair_of((size_t)n_cc * n_cc, -1);
  for (int ip = 0; ip < np; ++ip) {
    TwoBodyChannel_CC &tb = Z.modelspace->GetTwoBodyChannel_CC(B.cc_pairs[ip][0]);
    TwoBodyChannel_CC &tk = Z.modelspace->GetTwoBodyChannel_CC(B.cc_pairs[ip][1]);
    bar_Eta[ip] =
        arma::mat(2 * tb.GetNumberKets(), 2 * tk.GetNumberKets(),
                  arma::fill::zeros);
    pair_of[(size_t)B.cc_pairs[ip][0] * n_cc + B.cc_pairs[ip][1]] = ip;
  }

  // Ω̄^{Jb Jk}(a,b;c,d) from Ω^{J1 J2}(a,d;c,b) — same 9j as the elementwise
  // fill. Walk Omega2n once; write only canonical ch_bra ≤ ch_ket (partner
  // by the reduced CC flip).
  double t_bars = omp_get_wtime();
  std::vector<size_t> om_bra, om_ket;
  TensorChannelPairs(Eta, Z.modelspace, om_bra, om_ket);
  const int n_om = (int)om_bra.size();
  std::vector<arma::mat> Om2n((size_t)n_om);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < n_om; ++ip)
    Om2n[ip] = Omega2n(Eta, Z.modelspace, (int)om_bra[ip], (int)om_ket[ip]);
#pragma omp parallel for schedule(dynamic, 1)
  for (int iom = 0; iom < n_om; ++iom) {
    TwoBodyChannel &tbc1 = Z.modelspace->GetTwoBodyChannel((int)om_bra[iom]);
    TwoBodyChannel &tbc2 = Z.modelspace->GetTwoBodyChannel((int)om_ket[iom]);
    const int J1 = tbc1.J, J2 = tbc2.J;
    const arma::mat &Om = Om2n[iom];
    const int n1 = tbc1.GetNumberKets(), n2 = tbc2.GetNumberKets();
    const double hats12 = HatJ(J1) * HatJ(J2);
    for (int ibra = 0; ibra < 2 * n1; ++ibra) {
      auto ad = PqTB(tbc1, ibra, n1);
      if (ad[0] < 0)
        continue;
      const int a = ad[0], d = ad[1];
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &od = Z.modelspace->GetOrbit(d);
      const double ja = oa.j2 * 0.5, jd = od.j2 * 0.5;
      for (int iket = 0; iket < 2 * n2; ++iket) {
        auto cb = PqTB(tbc2, iket, n2);
        if (cb[0] < 0)
          continue;
        const double om = Om(ibra, iket);
        if (std::abs(om) < 1e-16)
          continue;
        const int c = cb[0], b = cb[1];
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &ob = Z.modelspace->GetOrbit(b);
        const double jc = oc.j2 * 0.5, jb = ob.j2 * 0.5;
        // dest: (a,b)=(a,b), (c,d)=(c,d)
        const int parity_ab = (oa.l + ob.l) % 2;
        const int Tz_ab = std::abs(oa.tz2 - ob.tz2) / 2;
        const int parity_cd = (oc.l + od.l) % 2;
        const int Tz_cd = std::abs(oc.tz2 - od.tz2) / 2;
        const int Jbmin = std::abs(oa.j2 - ob.j2) / 2;
        const int Jbmax = (oa.j2 + ob.j2) / 2;
        for (int Jb = Jbmin; Jb <= Jbmax; ++Jb) {
          const size_t ch_ab =
              Z.modelspace->GetTwoBodyChannelIndex(Jb, parity_ab, Tz_ab);
          if ((int)ch_ab >= n_cc)
            continue;
          TwoBodyChannel_CC &tcc_ab =
              Z.modelspace->GetTwoBodyChannel_CC(ch_ab);
          const int row = IndexCC(tcc_ab, a, b);
          if (row < 0)
            continue;
          const int Jkmin =
              std::max(std::abs(oc.j2 - od.j2) / 2, std::abs(Jb - lambda));
          const int Jkmax = std::min((oc.j2 + od.j2) / 2, Jb + lambda);
          for (int Jk = Jkmin; Jk <= Jkmax; ++Jk) {
            const size_t ch_cd =
                Z.modelspace->GetTwoBodyChannelIndex(Jk, parity_cd, Tz_cd);
            if ((int)ch_cd >= n_cc)
              continue;
            if ((int)ch_ab > (int)ch_cd)
              continue; // partner by CC flip
            const int ip = pair_of[(size_t)ch_ab * n_cc + ch_cd];
            if (ip < 0)
              continue;
            const int col =
                IndexCC(Z.modelspace->GetTwoBodyChannel_CC(ch_cd), c, d);
            if (col < 0)
              continue;
            const double ninej = Z.modelspace->GetNineJ(
                ja, jd, J1, jb, jc, J2, Jb, Jk, lambda);
            if (std::abs(ninej) < 1e-14)
              continue;
            const double hats = hats12 * HatJ(Jb) * HatJ(Jk);
            const double add =
                -hats *
                Z.modelspace->phase((ob.j2 + od.j2) / 2 + Jk + J2) * ninej *
                om;
            double *dst = &bar_Eta[ip](row, col);
#pragma omp atomic
            *dst += add;
          }
        }
      }
    }
  }
  Om2n.clear();
  Z.profiler.timer["_232_bar_Eta_scatter"] += omp_get_wtime() - t_bars;
  t_bars = omp_get_wtime();

  std::map<CcPair, int> ip_of;
  for (int ip = 0; ip < np; ++ip)
    ip_of[B.cc_pairs[ip]] = ip;
  for (int ip = 0; ip < np; ++ip) {
    const int ch_bra_cc = B.cc_pairs[ip][0], ch_ket_cc = B.cc_pairs[ip][1];
    if (ch_bra_cc <= ch_ket_cc)
      continue;
    auto it = ip_of.find({ch_ket_cc, ch_bra_cc});
    if (it == ip_of.end())
      continue;
    TwoBodyChannel_CC &tbc_bra_cc =
        Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
    TwoBodyChannel_CC &tbc_ket_cc =
        Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
    bar_Eta[ip] = TensorBarOmegaPartner(bar_Eta[it->second], tbc_ket_cc.J,
                                        tbc_bra_cc.J, hEta, Z.modelspace);
  }

  // Occupancy on the 2n CC layout (canon from scatter, partner from flip).
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < np; ++ip) {
    TwoBodyChannel_CC &tbc_bra_cc =
        Z.modelspace->GetTwoBodyChannel_CC(B.cc_pairs[ip][0]);
    TwoBodyChannel_CC &tbc_ket_cc =
        Z.modelspace->GetTwoBodyChannel_CC(B.cc_pairs[ip][1]);
    const int nbras_cc = tbc_bra_cc.GetNumberKets();
    const int nkets_cc = tbc_ket_cc.GetNumberKets();
    nnnbar_Eta[ip] = arma::mat(2 * nbras_cc, 2 * nkets_cc, arma::fill::zeros);
    nnnbar_Eta_d[ip] =
        arma::mat(2 * nbras_cc, 2 * nkets_cc, arma::fill::zeros);
    if (bar_Eta[ip].n_rows < 1)
      continue;
    for (int ibra_cc = 0; ibra_cc < 2 * nbras_cc; ++ibra_cc) {
      auto ab = PqCC(tbc_bra_cc, ibra_cc, nbras_cc);
      if (ab[0] < 0)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(ab[0]);
      Orbit &ob = Z.modelspace->GetOrbit(ab[1]);
      const double n_a = oa.occ, nbar_a = 1.0 - n_a;
      const double n_b = ob.occ, nbar_b = 1.0 - n_b;
      for (int iket_cc = 0; iket_cc < 2 * nkets_cc; ++iket_cc) {
        auto cd = PqCC(tbc_ket_cc, iket_cc, nkets_cc);
        if (cd[0] < 0)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(cd[0]);
        Orbit &od = Z.modelspace->GetOrbit(cd[1]);
        const double n_c = oc.occ, nbar_c = 1.0 - n_c;
        const double n_d = od.occ, nbar_d = 1.0 - n_d;
        const double occ_AbarBC = nbar_a * n_b * n_c + n_a * nbar_b * nbar_c;
        const double occ_ABbarD = n_a * nbar_b * n_d + nbar_a * n_b * nbar_d;
        const double Xbar = bar_Eta[ip](ibra_cc, iket_cc);
        if (std::abs(occ_AbarBC) >= 1e-12)
          nnnbar_Eta[ip](ibra_cc, iket_cc) = Xbar * occ_AbarBC;
        if (std::abs(occ_ABbarD) >= 1e-12)
          nnnbar_Eta_d[ip](ibra_cc, iket_cc) = Xbar * occ_ABbarD;
      }
    }
  }
  Z.profiler.timer["_232_bar_Eta_occ"] += omp_get_wtime() - t_bars;
  t_bars = omp_get_wtime();

  // FDC.cc ~1385: χ̄ in the Pandya representation, DGEMM per CC block.
  //   bar_CHI_V[{ch_il,ch_kj}]     = bar_Gamma[ch_il] · nnnbar_Eta[{ch_il,ch_kj}]
  //                                  (χ^ι: Σ_ab occ_AbarBC(a,b,k) Γ̄(il;ab) Ω̄(ab;kj))
  //   bar_CHI_VI_II[{ch_il,ch_kj}] = hΩ (−1)^{J_il+J_kj}
  //                                  nnnbar_Eta_d[{ch_kj,ch_il}]ᵀ · bar_Gamma[ch_kj]
  //                                  (χ^κ: Σ_ab occ_ABbarD(a,b,l) Ω̄(ab;il) Γ̄(ab;kj))
  std::vector<arma::mat> bar_CHI_V(np), bar_CHI_VI_II(np);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < np; ++ip) {
    const int ch_il = B.cc_pairs[ip][0], ch_kj = B.cc_pairs[ip][1];
    if (bar_Gamma[ch_il].n_rows < 1 or bar_Gamma[ch_kj].n_rows < 1)
      continue;
    bar_CHI_V[ip] = bar_Gamma[ch_il] * nnnbar_Eta[ip];
    auto it = ip_of.find({ch_kj, ch_il});
    if (it == ip_of.end())
      continue;
    TwoBodyChannel_CC &tbc_il = Z.modelspace->GetTwoBodyChannel_CC(ch_il);
    TwoBodyChannel_CC &tbc_kj = Z.modelspace->GetTwoBodyChannel_CC(ch_kj);
    const double phaseJJ = Z.modelspace->phase(tbc_il.J + tbc_kj.J);
    arma::mat nnnbar_Eta_d_t = nnnbar_Eta_d[it->second].t();
    bar_CHI_VI_II[ip] = (hEta * phaseJJ) * (nnnbar_Eta_d_t * bar_Gamma[ch_kj]);
  }
  nnnbar_Eta.clear();
  nnnbar_Eta_d.clear();

  B.bar_Eta.clear();
  B.bar_CHI_V.clear();
  B.bar_CHI_VI_II.clear();
  for (int ip = 0; ip < np; ++ip) {
    B.bar_Eta[B.cc_pairs[ip]] = std::move(bar_Eta[ip]);
    if (bar_CHI_V[ip].n_rows > 0)
      B.bar_CHI_V[B.cc_pairs[ip]] = std::move(bar_CHI_V[ip]);
    if (bar_CHI_VI_II[ip].n_rows > 0)
      B.bar_CHI_VI_II[B.cc_pairs[ip]] = std::move(bar_CHI_VI_II[ip]);
  }
  Z.profiler.timer["_232_bar_CHI_dgemm"] += omp_get_wtime() - t_bars;
}
} // namespace

void comm223_232_GIVa(const Operator &Eta_in, const Operator &Gamma,
                      Operator &Z) {
  const ReducedEtaView E(Eta_in);
  PandyaBars232 B;
  comm223_232_GIVa_from_bars(E.ref, Gamma, Z, B);
}

namespace {
////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_a / chi^kappa.  FDC.cc ~1528–1760 (bar_CHI_VI_II → Chi_VI_II_Op
/// → −Chi_VI_II_Op·Eta and its Hermitian twin), tensor Ω:
///   1. bar_CHI_VI_II from the shared Pandya block (DGEMM, FillPandyaBars232)
///   2. inverse Pandya on ch_bra×ch_ket (9j) → Chi_VI_II_Op (2n, non-AS)
///   3. ladder DGEMM per Z channel, then (1−Pij) and the (kl;ij) twin.
/// AMC G4a_Wbra (reduced Ω, unreduced scalar Z):
///   W^{J0}(ij;kl) = −(−1)^{J0}/Ĵ0 λ̂⁻¹ Σ_{J2} Σ_{bd}
///                   χ^{J0 J2}(ij;bd) (−1)^{jb+jd+λ} Ω^{J2 J0}(db;kl)
///   Z(ij;kl) += [W(ij;kl) − Pij W(ji;kl) + W(kl;ij) − Pkl W(lk;ij)] / Ĵ0
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVa_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const PandyaBars232 &B_in) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  PandyaBars232 B_local;
  const PandyaBars232 *Bp = &B_in;
  if (B_in.cc_pairs.empty()) {
    FillPandyaBars232(Eta, Gamma, Z, B_local);
    Bp = &B_local;
  }
  const PandyaBars232 &B = *Bp;
  if (B.bar_CHI_VI_II.empty()) {
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  }

  ModelSpace *ms = Z.modelspace;
  const int lambda = Eta.GetJRank();
  const double hat_lam_inv = 1.0 / HatJ(lambda);
  auto &Z2 = Z.TwoBody;

  // FDC.cc ~1553: Chi_VI_II_Op by inverse Pandya over ch_bra×ch_ket.
  double t_internal = omp_get_wtime();
  std::vector<size_t> ch_bra_list, ch_ket_list;
  TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
  ChiOp2n Chi_VI_II_Op;
  Chi_VI_II_Op.Init(ms);
  Chi_VI_II_Op.Allocate(ch_bra_list, ch_ket_list);
  InvPandyaTensor2n(B.bar_CHI_VI_II, Eta, Chi_VI_II_Op);
  Z.profiler.timer["_GIVa_chi_inv"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // FDC.cc ~1752: Z2 += −Chi_VI_II_Op · Eta (+ twin), per Z channel.
  const int nch_tb = (int)ms->GetNumberTwoBodyChannels();
  std::vector<std::vector<int>> ket_channels(nch_tb);
  for (size_t ich = 0; ich < ch_bra_list.size(); ++ich)
    ket_channels[ch_bra_list[ich]].push_back((int)ch_ket_list[ich]);

  std::vector<size_t> zch_list;
  for (auto &iter : Z2.MatEl)
    if (iter.first[0] == iter.first[1])
      zch_list.push_back(iter.first[0]);
  const int nch = (int)zch_list.size();

#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ++ich) {
    const int ch = (int)zch_list[ich];
    TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
    const int J0 = tbc.J;
    const int nkets = tbc.GetNumberKets();
    const int n0 = Chi_VI_II_Op.NumberKets(ch);
    if (nkets < 1 or n0 < 1)
      continue;

    // W^{J0}(ij;kl) on the 2n layout of channel ch (rows and cols).
    arma::mat W(2 * n0, 2 * n0, arma::fill::zeros);
    for (int ch2 : ket_channels[ch]) {
      const arma::mat *Chi = Chi_VI_II_Op.GetBlock(ch, ch2);
      if (Chi == nullptr)
        continue;
      TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
      const int nK2 = tbc2.GetNumberKets();
      const int n2 = Chi_VI_II_Op.NumberKets(ch2);
      if (nK2 < 1)
        continue;
      // Ω^{J2 J0}(bd;kl) in Pauli 2n, then row-swap to (db) and copy into
      // the ChiOp2n layout. (p,p) at odd J is Pauli-forbidden so Ω=0 there.
      const arma::mat Om = Omega2n(Eta, ms, ch2, ch);
      arma::mat Eta_matrix(2 * n2, 2 * n0, arma::fill::zeros);
      for (int ibd = 0; ibd < 2 * n2; ++ibd) {
        auto bd = Chi_VI_II_Op.Pq(ch2, ibd);
        if (bd[0] < 0)
          continue;
        const int idx_db = Index2n(tbc2, bd[1], bd[0]);
        if (idx_db < 0)
          continue;
        Orbit &ob = ms->GetOrbit(bd[0]);
        Orbit &od = ms->GetOrbit(bd[1]);
        const double ph_bd = ms->phase((ob.j2 + od.j2) / 2 + lambda);
        for (int ikl = 0; ikl < 2 * n0; ++ikl) {
          auto kl = Chi_VI_II_Op.Pq(ch, ikl);
          if (kl[0] < 0)
            continue;
          const int idx_kl = Index2n(tbc, kl[0], kl[1]);
          if (idx_kl < 0)
            continue;
          Eta_matrix(ibd, ikl) = ph_bd * Om(idx_db, idx_kl);
        }
      }
      W += (*Chi) * Eta_matrix;
    }
    W *= -ms->phase(J0) / HatJ(J0) * hat_lam_inv;

    // (1−Pij) W + (kl;ij) twin, then Z normalization.
    arma::mat &Zmat = Z2.GetMatrix(ch, ch);
    for (int ibra = 0; ibra < nkets; ++ibra) {
      Ket &bra = tbc.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const int idx_ij = Chi_VI_II_Op.Index(ch, i, j);
      const int idx_ji = (i == j) ? idx_ij : Chi_VI_II_Op.Index(ch, j, i);
      const double Pij = ms->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const int idx_kl = Chi_VI_II_Op.Index(ch, k, l);
        const int idx_lk = (k == l) ? idx_kl : Chi_VI_II_Op.Index(ch, l, k);
        const double Pkl = ms->phase((ok.j2 + ol.j2) / 2 - J0);
        double z = (W(idx_ij, idx_kl) - Pij * W(idx_ji, idx_kl)) +
                   (W(idx_kl, idx_ij) - Pkl * W(idx_lk, idx_ij));
        z /= HatJ(J0);
        if (i == j)
          z /= PhysConst::SQRT2;
        if (k == l)
          z /= PhysConst::SQRT2;
        Zmat(ibra, iket) += z;
      }
    }
  }
  Z.profiler.timer["_GIVa_leftover"] += omp_get_wtime() - t_internal;

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
void comm223_232_GIVb(const Operator &Eta_in, const Operator &Gamma,
                      Operator &Z) {
  const ReducedEtaView E(Eta_in);
  PandyaBars232 B;
  comm223_232_GIVb_from_bars(E.ref, Gamma, Z, B);
}

namespace {
////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_b / chi^iota.  FDC.cc ~1906 (bar_CHI_V → recouple → DGEMM with
/// bar_Eta → inverse Pandya), tensor Ω:
///   1. bar_CHI_V from the shared Pandya block (DGEMM, FillPandyaBars232)
///   2. inverse Pandya on ch_bra×ch_ket (9j) → Chi_V_Op (2n, non-AS)
///   3. leftover: per CC channel (J6) Pandya-shaped DGEMM over (a,b) for
///      each half-integer j0, then the 6j inverse on Z's ch_bra×ch_ket.
/// AMC G4b_W{1,2}_leftover.tex (unreduced leftover), Z = (1−Pij)(1−Pkl) W:
///   W1 = Σ CG χ_aibk Ω_jbla,   W2 = −hΩ Σ CG χ_akbi Ω_lbja.
/// Do not dagger χ; do not pack χ_adbc − hZ χ_bcad with one Ω̄ (that identity
/// Ω_jalb ≡ hΩ Ω_lbja holds only at λ = 0).
////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVb_from_bars(const Operator &Eta, const Operator &Gamma,
                               Operator &Z, const PandyaBars232 &B_in) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  PandyaBars232 B_local;
  const PandyaBars232 *Bp = &B_in;
  if (B_in.cc_pairs.empty()) {
    FillPandyaBars232(Eta, Gamma, Z, B_local);
    Bp = &B_local;
  }
  const PandyaBars232 &B = *Bp;
  if (B.bar_CHI_V.empty()) {
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  }

  ModelSpace *ms = Z.modelspace;
  const int lambda = Eta.GetJRank();
  const int hEta = Eta.IsHermitian() ? 1 : -1;
  const int rank_T = Eta.GetTRank();
  const int parity_eta = Eta.GetParity();
  const double hat_lam_inv = 1.0 / HatJ(lambda);
  auto &Z2 = Z.TwoBody;

  // Chi_V_Op by inverse Pandya over ch_bra×ch_ket (FDC.cc Chi_*_Op).
  double t_internal = omp_get_wtime();
  std::vector<size_t> ch_bra_list, ch_ket_list;
  TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
  const int nch_pairs = (int)ch_bra_list.size();
  ChiOp2n Chi_V_Op;
  Chi_V_Op.Init(ms);
  Chi_V_Op.Allocate(ch_bra_list, ch_ket_list);
  InvPandyaTensor2n(B.bar_CHI_V, Eta, Chi_V_Op);
  Z.profiler.timer["_GIVb_chi_inv"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Leftover W1/W2 (AMC G4b_W{1,2}_leftover.tex, five 6j). The (a,b)
  // contraction is Pandya-shaped: rows (i,k) and (j,l) both live in the CC
  // channel ch_cc = (J6, π_ik, Tz_ik) of Z's scalar ket; the summed (a,b) run
  // over all ordered pairs of the coupled (π,Tz) class (no J restriction, a
  // couples to λ through the half-integer j0). Per ch_cc and j0:
  //   X1(ik;ab) = (−1)^{ja+jb+λ} Σ_{J2 J3} Ĵ2Ĵ3 {J3 λ J2}{ji jk J6} χ^{J2 J3}(ai;bk)
  //                                            {ja ji j0}{jb j0 J3}
  //   Y1(jl;ab) =                Σ_{J4 J5} Ĵ4Ĵ5 {J4 λ J5}{jl jj J6} Ω^{J4 J5}(jb;la)
  //                                            {ja jl j0}{jb j0 J4}
  //   X2, Y2: i↔k in χ (χ(ak;bi)), j↔l in Ω (Ω(lb;ja)).
  //   CHI_V_final_{1,2}[ch_cc] += (2j0+1)/λ̂ · X·Yᵀ        (FDC CHI_V_final)
  // then the scalar 6j inverse on Z's ch_bra×ch_ket:
  //   W1(ij;kl) = −(−1)^{J0+jj+jk} Σ_{J6} (2J6+1) {jk jl J0} CHI_V_final_1(ik;jl)
  //                                                {jj ji J6}
  //   W2(ij;kl) = −(−1)^{J0+ji+jl} Σ_{J6} (2J6+1) {..} CHI_V_final_2(ik;jl)
  //   Z2 += (1−Pij)(1−Pkl) (W1 − hΩ W2)
  const int n_cc = (int)ms->GetNumberTwoBodyChannels_CC();
  std::vector<index_t> allorb(ms->all_orbits.begin(), ms->all_orbits.end());
  int max_j2 = 0;
  for (auto o : allorb)
    max_j2 = std::max(max_j2, ms->GetOrbit(o).j2);
  const int j0_2max = max_j2 + 2 * lambda;

  // Ordered (a,b) pairs by (π_ab, Tz_ab = |tz_a − tz_b|/2); a == b once.
  std::array<std::array<std::vector<std::array<int, 2>>, 2>, 2> ab_pairs;
  for (auto a : allorb) {
    Orbit &oa = ms->GetOrbit(a);
    for (auto b : allorb) {
      Orbit &ob = ms->GetOrbit(b);
      ab_pairs[(oa.l + ob.l) % 2][std::abs(oa.tz2 - ob.tz2) / 2].push_back(
          {(int)a, (int)b});
    }
  }

  std::vector<CcPair> chi_keys;
  std::vector<const arma::mat *> chi_mats;
  for (const auto &kv : Chi_V_Op.MatEl) {
    chi_keys.push_back(kv.first);
    chi_mats.push_back(&kv.second);
  }
  const int n_chi = (int)chi_keys.size();
  std::vector<arma::mat> Om2n((size_t)nch_pairs);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < nch_pairs; ++ip)
    Om2n[ip] = Omega2n(Eta, ms, (int)ch_bra_list[ip], (int)ch_ket_list[ip]);
  int n_alloc = 0;
  for (auto o : allorb)
    n_alloc = std::max(n_alloc, (int)o + 1);
  const int n_j0 = (j0_2max + 1) / 2;

  // Per-CC metadata. ab_col is 2n×n_ab of that channel; X,Y for one j0 over
  // all channels (~100 MB at emax=5), not n_j0×n_cc.
  std::vector<int> cc_nk(n_cc, 0), cc_nab(n_cc, 0);
  std::vector<std::vector<int>> cc_abcol(n_cc);
  std::deque<arma::mat> CHI_V_final_1(n_cc), CHI_V_final_2(n_cc);
  for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
    const int nkets_cc = tbc_cc.GetNumberKets();
    if (nkets_cc < 1)
      continue;
    const int parity_ab = (tbc_cc.parity + parity_eta) % 2;
    std::vector<std::array<int, 2>> ab_list;
    for (int Tz_ab = 0; Tz_ab < 2; ++Tz_ab)
      if (CcTzCouples(Tz_ab, tbc_cc.Tz, rank_T))
        ab_list.insert(ab_list.end(), ab_pairs[parity_ab][Tz_ab].begin(),
                       ab_pairs[parity_ab][Tz_ab].end());
    if (ab_list.empty())
      continue;
    cc_nk[ch_cc] = nkets_cc;
    cc_nab[ch_cc] = (int)ab_list.size();
    cc_abcol[ch_cc].assign((size_t)n_alloc * n_alloc, -1);
    for (int col = 0; col < cc_nab[ch_cc]; ++col)
      cc_abcol[ch_cc][(size_t)ab_list[col][0] * n_alloc + ab_list[col][1]] =
          col;
    CHI_V_final_1[ch_cc].zeros(2 * nkets_cc, 2 * nkets_cc);
    CHI_V_final_2[ch_cc].zeros(2 * nkets_cc, 2 * nkets_cc);
  }

  auto add_xy = [&](arma::mat &M, int ch, int row, int col, double add) {
    if (ch < 0 or ch >= n_cc or row < 0 or col < 0)
      return;
    if (M.n_rows < 1 or row >= (int)M.n_rows or col >= (int)M.n_cols)
      return;
    if (std::abs(add) < 1e-16)
      return;
    double *dst = &M(row, col);
#pragma omp atomic
    *dst += add;
  };

  // Walk χ/Ω once. Dummy j0 only re-evaluates the two 6j (same formula for
  // X and Y after packing ja,ji,jb,jk). 6j:
  //   { Jx λ Jy ; ja ji j0 } { ji jk J6 ; jb j0 Jx }
  struct GivbRec {
    int ch, row1, row2, col, J6, Jx, Jy, j2a, j2i, j2b, j2k;
    double w;
  };
  const int nthr = std::max(1, omp_get_max_threads());
  std::vector<std::vector<GivbRec>> tlsX((size_t)nthr), tlsY((size_t)nthr);
  const double t_pack = omp_get_wtime();
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    auto &mineX = tlsX[(size_t)tid];
#pragma omp for schedule(dynamic, 1) nowait
    for (int ic = 0; ic < n_chi; ++ic) {
      const int ch2 = chi_keys[ic][0], ch3 = chi_keys[ic][1];
      const int J2 = ms->GetTwoBodyChannel(ch2).J;
      const int J3 = ms->GetTwoBodyChannel(ch3).J;
      const arma::mat &Chi = *chi_mats[ic];
      const int n2 = Chi_V_Op.NumberKets(ch2), n3 = Chi_V_Op.NumberKets(ch3);
      const double hats = HatJ(J2) * HatJ(J3);
      for (int ibra = 0; ibra < 2 * n2; ++ibra) {
        auto pq = Chi_V_Op.Pq(ch2, ibra);
        if (pq[0] < 0)
          continue;
        const int p = pq[0], q = pq[1];
        Orbit &op = ms->GetOrbit(p);
        Orbit &oq = ms->GetOrbit(q);
        for (int iket = 0; iket < 2 * n3; ++iket) {
          auto rs = Chi_V_Op.Pq(ch3, iket);
          if (rs[0] < 0)
            continue;
          const double chi = Chi(ibra, iket);
          if (std::abs(chi) < 1e-16)
            continue;
          const int r = rs[0], s = rs[1];
          Orbit &orr = ms->GetOrbit(r);
          Orbit &os = ms->GetOrbit(s);
          const int parity_qs = (oq.l + os.l) % 2;
          const int Tz_qs = std::abs(oq.tz2 - os.tz2) / 2;
          const double w = ms->phase((op.j2 + orr.j2) / 2 + lambda) * hats * chi;
          const int J6min = std::abs(oq.j2 - os.j2) / 2;
          const int J6max = (oq.j2 + os.j2) / 2;
          for (int J6 = J6min; J6 <= J6max; ++J6) {
            const size_t ch =
                ms->GetTwoBodyChannelIndex(J6, parity_qs, Tz_qs);
            if ((int)ch >= n_cc or cc_nk[ch] < 1)
              continue;
            const int col = cc_abcol[ch][(size_t)p * n_alloc + r];
            if (col < 0)
              continue;
            TwoBodyChannel_CC &tcc = ms->GetTwoBodyChannel_CC(ch);
            const int row1 = IndexCC(tcc, q, s);
            const int row2 = IndexCC(tcc, s, q);
            if (row1 < 0 and row2 < 0)
              continue;
            mineX.push_back({(int)ch, row1, row2, col, J6, J3, J2, op.j2,
                             oq.j2, orr.j2, os.j2, w});
          }
        }
      }
    }
    auto &mineY = tlsY[(size_t)tid];
#pragma omp for schedule(dynamic, 1)
    for (int ip = 0; ip < nch_pairs; ++ip) {
      TwoBodyChannel &tbc4 = ms->GetTwoBodyChannel((int)ch_bra_list[ip]);
      TwoBodyChannel &tbc5 = ms->GetTwoBodyChannel((int)ch_ket_list[ip]);
      const int J4 = tbc4.J, J5 = tbc5.J;
      const arma::mat &Om = Om2n[ip];
      const int n4 = tbc4.GetNumberKets(), n5 = tbc5.GetNumberKets();
      const double hats = HatJ(J4) * HatJ(J5);
      for (int ibra = 0; ibra < 2 * n4; ++ibra) {
        auto pq = PqTB(tbc4, ibra, n4);
        if (pq[0] < 0)
          continue;
        const int p = pq[0], q = pq[1];
        Orbit &op = ms->GetOrbit(p);
        Orbit &oq = ms->GetOrbit(q);
        for (int iket = 0; iket < 2 * n5; ++iket) {
          auto rs = PqTB(tbc5, iket, n5);
          if (rs[0] < 0)
            continue;
          const double om = Om(ibra, iket);
          if (std::abs(om) < 1e-16)
            continue;
          const int r = rs[0], s = rs[1];
          Orbit &orr = ms->GetOrbit(r);
          Orbit &os = ms->GetOrbit(s);
          const int parity_pr = (op.l + orr.l) % 2;
          const int Tz_pr = std::abs(op.tz2 - orr.tz2) / 2;
          const double w = hats * om;
          const int J6min = std::abs(op.j2 - orr.j2) / 2;
          const int J6max = (op.j2 + orr.j2) / 2;
          for (int J6 = J6min; J6 <= J6max; ++J6) {
            const size_t ch =
                ms->GetTwoBodyChannelIndex(J6, parity_pr, Tz_pr);
            if ((int)ch >= n_cc or cc_nk[ch] < 1)
              continue;
            const int col = cc_abcol[ch][(size_t)s * n_alloc + q];
            if (col < 0)
              continue;
            TwoBodyChannel_CC &tcc = ms->GetTwoBodyChannel_CC(ch);
            const int row1 = IndexCC(tcc, p, r);
            const int row2 = IndexCC(tcc, r, p);
            if (row1 < 0 and row2 < 0)
              continue;
            mineY.push_back({(int)ch, row1, row2, col, J6, J4, J5, os.j2,
                             orr.j2, oq.j2, op.j2, w});
          }
        }
      }
    }
  }
  Om2n.clear();
  std::vector<GivbRec> recX, recY;
  size_t nx = 0, ny = 0;
  for (auto &v : tlsX)
    nx += v.size();
  for (auto &v : tlsY)
    ny += v.size();
  recX.reserve(nx);
  recY.reserve(ny);
  for (auto &v : tlsX)
    recX.insert(recX.end(), v.begin(), v.end());
  for (auto &v : tlsY)
    recY.insert(recY.end(), v.begin(), v.end());
  tlsX.clear();
  tlsY.clear();
  Z.profiler.timer["_GIVb_pack"] += omp_get_wtime() - t_pack;
  const int nrx = (int)recX.size(), nry = (int)recY.size();

  for (int ij0 = 0; ij0 < n_j0; ++ij0) {
    const double j0 = 0.5 * (2 * ij0 + 1);
    const double sc = (2.0 * ij0 + 2.0) * hat_lam_inv;
    std::vector<arma::mat> X1(n_cc), X2(n_cc), Y1(n_cc), Y2(n_cc);
    for (int ch = 0; ch < n_cc; ++ch) {
      if (cc_nk[ch] < 1)
        continue;
      X1[ch].zeros(2 * cc_nk[ch], cc_nab[ch]);
      X2[ch].zeros(2 * cc_nk[ch], cc_nab[ch]);
      Y1[ch].zeros(2 * cc_nk[ch], cc_nab[ch]);
      Y2[ch].zeros(2 * cc_nk[ch], cc_nab[ch]);
    }

#pragma omp parallel for schedule(static)
    for (int ir = 0; ir < nrx; ++ir) {
      const GivbRec &R = recX[ir];
      const double ja = 0.5 * R.j2a, ji = 0.5 * R.j2i;
      const double jb = 0.5 * R.j2b, jk = 0.5 * R.j2k;
      if (not AngMom::Triangle(ja, (double)lambda, j0))
        continue;
      const double six =
          ms->GetSixJ(R.Jx, lambda, R.Jy, ja, ji, j0) *
          ms->GetSixJ(ji, jk, R.J6, jb, j0, R.Jx);
      const double add = R.w * six;
      add_xy(X1[R.ch], R.ch, R.row1, R.col, add);
      add_xy(X2[R.ch], R.ch, R.row2, R.col, add);
    }
#pragma omp parallel for schedule(static)
    for (int ir = 0; ir < nry; ++ir) {
      const GivbRec &R = recY[ir];
      const double ja = 0.5 * R.j2a, ji = 0.5 * R.j2i;
      const double jb = 0.5 * R.j2b, jk = 0.5 * R.j2k;
      if (not AngMom::Triangle(ja, (double)lambda, j0))
        continue;
      const double six =
          ms->GetSixJ(R.Jx, lambda, R.Jy, ja, ji, j0) *
          ms->GetSixJ(ji, jk, R.J6, jb, j0, R.Jx);
      const double add = R.w * six;
      add_xy(Y1[R.ch], R.ch, R.row1, R.col, add);
      add_xy(Y2[R.ch], R.ch, R.row2, R.col, add);
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < n_cc; ++ch) {
      if (cc_nk[ch] < 1)
        continue;
      CHI_V_final_1[ch] += sc * (X1[ch] * Y1[ch].t());
      CHI_V_final_2[ch] += sc * (X2[ch] * Y2[ch].t());
    }
  }
  Z.profiler.timer["_GIVb_cc_dgemm"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

  // Inverse (scalar 6j) on Z's channels, then (1−Pij)(1−Pkl).
  std::vector<size_t> zch_list;
  for (auto &iter : Z2.MatEl)
    if (iter.first[0] == iter.first[1])
      zch_list.push_back(iter.first[0]);
  const int nch = (int)zch_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ++ich) {
    const int ch = (int)zch_list[ich];
    TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
    const int J0 = tbc.J;
    const int nkets = tbc.GetNumberKets();
    if (nkets < 1)
      continue;
    arma::mat W2n(2 * nkets, 2 * nkets, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nkets; ++ibra) {
      auto ij = PqTB(tbc, ibra, nkets);
      if (ij[0] < 0)
        continue;
      const int i = ij[0], j = ij[1];
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const double ji = oi.j2 * 0.5, jj = oj.j2 * 0.5;
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = PqTB(tbc, iket, nkets);
        if (kl[0] < 0)
          continue;
        const int k = kl[0], l = kl[1];
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const double jk = ok.j2 * 0.5, jl = ol.j2 * 0.5;
        const int parity_cc = (oi.l + ok.l) % 2;
        const int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        const int J6min =
            std::max(std::abs(oi.j2 - ok.j2), std::abs(oj.j2 - ol.j2)) / 2;
        const int J6max = std::min(oi.j2 + ok.j2, oj.j2 + ol.j2) / 2;
        const double ph1 = -ms->phase(J0 + (oj.j2 + ok.j2) / 2);
        const double ph2 = -ms->phase(J0 + (oi.j2 + ol.j2) / 2);
        double w = 0.0;
        for (int J6 = J6min; J6 <= J6max; ++J6) {
          const size_t ch_cc = ms->GetTwoBodyChannelIndex(J6, parity_cc, Tz_cc);
          if ((int)ch_cc >= n_cc or CHI_V_final_1[ch_cc].n_rows < 1)
            continue;
          const double six5 = ms->GetSixJ(jk, jl, J0, jj, ji, J6);
          if (std::abs(six5) < 1e-16)
            continue;
          TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
          const int indx_ik = (int)tbc_cc.GetLocalIndex(i, k);
          const int indx_jl = (int)tbc_cc.GetLocalIndex(j, l);
          const double pref = (2.0 * J6 + 1.0) * six5;
          const double w1 = ph1 * pref * CHI_V_final_1[ch_cc](indx_ik, indx_jl);
          const double w2 = ph2 * pref * CHI_V_final_2[ch_cc](indx_ik, indx_jl);
          w += w1 - hEta * w2;
        }
        W2n(ibra, iket) = w;
      }
    }

    arma::mat &Zmat = Z2.GetMatrix(ch, ch);
    for (int ibra = 0; ibra < nkets; ++ibra) {
      Ket &bra = tbc.GetKet(ibra);
      const int i = bra.p, j = bra.q;
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const int ibra_s = (i == j) ? ibra : ibra + nkets;
      const double Pij = ms->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = 0; iket < nkets; ++iket) {
        Ket &ket = tbc.GetKet(iket);
        const int k = ket.p, l = ket.q;
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const int iket_s = (k == l) ? iket : iket + nkets;
        const double Pkl = ms->phase((ok.j2 + ol.j2) / 2 - J0);
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

  Z.profiler.timer["_GIVb_leftover"] += omp_get_wtime() - t_internal;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}
} // namespace

////////////////////////////////////////////////////////////////////////////
/// Gamma^IV_c / chi^lambda — χ DGEMM then leftover of χ_{ialb} Ω_{bjak}.
/// Rank: T×S / S×T → tensor χ^λ. Fold T×T→S.
///
/// χ^λ is occupancy-AS (no h_χ). Leftover is gold Path B for all π:
/// scatter χ/Ω 2n blocks with the PandyaChiIalb / PandyaOmegaBjak 9j
/// → mid-J DGEMM → inv 6j, then W2n (1−P)².
/// (Ring adcb ≡ gold only when T1 χ is accidentally AS — even π.)
////////////////////////////////////////////////////////////////////////////
static void comm223_232_GIVc_pathB(const Operator &Eta, const Operator &Gamma,
                                   Operator &Z) {
  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();
  Z.modelspace->PreCalculateNineJ();

  ModelSpace *ms = Z.modelspace;
  const int lambda = Eta.GetJRank();
  PreCalculateNineJTensorPandya(ms, lambda); // PandyaChiIalb / PandyaOmegaBjak
  auto hat = [](double x) { return std::sqrt(2.0 * x + 1.0); };

  // χ^λ (FDC.cc CHI_VII) in ordinary-channel blocks {ch_bra, ch_ket}, 2n
  // layout, GetTBME_J units (AMC chi_lambda.tex). Occupancy-AS, so ChiOp2n.
  //   CHI_VII^{J0 J1} = Γ^{J0} · (w(a,b,l) ⊙ Ω^{J0 J1}) + (w(a,b,j) ⊙ Ω^{J0 J1}) · Γ^{J1}
  //   w(a,b,x) = n̄_a n̄_b n_x + n_a n_b n̄_x
  std::vector<size_t> ch_bra_list, ch_ket_list;
  TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
  ChiOp2n CHI_VII;
  CHI_VII.Init(ms);
  CHI_VII.Allocate(ch_bra_list, ch_ket_list);
  const int nch_pairs = (int)ch_bra_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch_pairs; ++ich) {
    const int ch_bra = (int)ch_bra_list[ich], ch_ket = (int)ch_ket_list[ich];
    arma::mat *Chi = CHI_VII.GetBlock(ch_bra, ch_ket);
    if (Chi == nullptr)
      continue;
    TwoBodyChannel &tbc_bra = ms->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = ms->GetTwoBodyChannel(ch_ket);
    const int nbras = tbc_bra.GetNumberKets(), nkets = tbc_ket.GetNumberKets();
    if (nbras < 1 or nkets < 1)
      continue;
    const arma::mat Gamma_matrix_bra = Omega2n(Gamma, ms, ch_bra, ch_bra);
    const arma::mat Gamma_matrix_ket = Omega2n(Gamma, ms, ch_ket, ch_ket);
    const arma::mat Eta_matrix = Omega2n(Eta, ms, ch_bra, ch_ket);

    arma::mat T(2 * nbras, 2 * nkets, arma::fill::zeros);
    // T1 = Γ^{J0} · (w(a,b,l) ⊙ Ω^{J0 J1})
    if (givc_chi_which != 2) {
      arma::mat Eta_matrix_c = Eta_matrix;
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = PqTB(tbc_ket, iket, nkets);
        if (kl[0] < 0)
          continue;
        const double n_l = ms->GetOrbit(kl[1]).occ;
        for (int ib = 0; ib < 2 * nbras; ++ib) {
          auto ab = PqTB(tbc_bra, ib, nbras);
          if (ab[0] < 0)
            continue;
          const double n_a = ms->GetOrbit(ab[0]).occ;
          const double n_b = ms->GetOrbit(ab[1]).occ;
          Eta_matrix_c(ib, iket) *=
              (1.0 - n_a) * (1.0 - n_b) * n_l + n_a * n_b * (1.0 - n_l);
        }
      }
      T += Gamma_matrix_bra * Eta_matrix_c;
    }
    // T2 = (w(a,b,j) ⊙ Ω^{J0 J1}) · Γ^{J1}
    if (givc_chi_which != 1) {
      arma::mat Eta_matrix_d = Eta_matrix;
      for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
        auto ij = PqTB(tbc_bra, ibra, nbras);
        if (ij[0] < 0)
          continue;
        const double n_j = ms->GetOrbit(ij[1]).occ;
        for (int ik = 0; ik < 2 * nkets; ++ik) {
          auto ab = PqTB(tbc_ket, ik, nkets);
          if (ab[0] < 0)
            continue;
          const double n_a = ms->GetOrbit(ab[0]).occ;
          const double n_b = ms->GetOrbit(ab[1]).occ;
          Eta_matrix_d(ibra, ik) *=
              (1.0 - n_a) * (1.0 - n_b) * n_j + n_a * n_b * (1.0 - n_j);
        }
      }
      T += Eta_matrix_d * Gamma_matrix_ket;
    }

    // TwoBodyChannel 2n → ChiOp2n 2n (superset ket list: (p,p) at odd J).
    for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
      auto ij = PqTB(tbc_bra, ibra, nbras);
      if (ij[0] < 0)
        continue;
      const int row = CHI_VII.Index(ch_bra, ij[0], ij[1]);
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = PqTB(tbc_ket, iket, nkets);
        if (kl[0] < 0)
          continue;
        const int col = CHI_VII.Index(ch_ket, kl[0], kl[1]);
        (*Chi)(row, col) = T(ibra, iket);
      }
    }
  }

  Z.profiler.timer["_GIVc_chi"] += omp_get_wtime() - t_start;

  // ------------------------------------------------------------------
  // Leftover: gold Path B Pandya → DGEMM → inv (all π).
  //
  // m-gold is Σ χ_ialb Ω_bjak. Ring adcb implements χ_ibal Ω_ajkb, which
  // equals gold only when T1 χ is accidentally AS (even π). Odd-π T1 fails.
  //
  // Gold fills (AMC scheme (il)|(ab) / (ab)|(kj)):
  //   PandyaChiIalb  — scheme=((1,-3),(2,-4)) on χ_ialb
  //   PandyaOmegaBjak — scheme=((3,-1),(4,-2)) on Ω_bjak
  // Same ring mid w = (−1)^{Jp}/Ĵp · λ̂^{-1}(−1)^{Jab+λ} and corrected inv.
  // Locked: even ≡ adcb; odd T1 ≡ CG/m (givc_gold_dgemm_lock).
  //
  // T1-only leftover is not Hermitian (T1† ≡ T2); store NonHerm when
  // givc_chi_which != 0 so each fold(p,g,q,h) is kept.
  // ------------------------------------------------------------------
  const double hat_lam_inv =
      1.0 / std::sqrt(2.0 * std::max(lambda, 0) + 1.0);

  // Packaging: χ/Ω WE-reduced; CHI_VII_final/Z reduce=true; store
  // Z_unred = Z_red/Ĵ. 2n unnormalized (GetTBME_J); ÷√2 only on AddToTBME.
  //
  // FDC.cc ~2394 names: bar_CHI_VII_CC (Pandya of χ^λ), bar_Eta_CC (Pandya of
  // Ω in the gold bjak scheme), CHI_VII_final (their CC product).
  const int n_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();
  const int parity_eta = Eta.GetParity();
  const int rank_T = Eta.GetTRank();
  std::vector<CcPair> cc_pairs;
  for (int ch_bra_cc = 0; ch_bra_cc < n_cc; ++ch_bra_cc) {
    TwoBodyChannel_CC &tbc_bra_cc = ms->GetTwoBodyChannel_CC(ch_bra_cc);
    if (tbc_bra_cc.GetNumberKets() < 1)
      continue;
    for (int ch_ket_cc = 0; ch_ket_cc < n_cc; ++ch_ket_cc) {
      TwoBodyChannel_CC &tbc_ket_cc = ms->GetTwoBodyChannel_CC(ch_ket_cc);
      if (tbc_ket_cc.GetNumberKets() < 1)
        continue;
      if ((tbc_bra_cc.parity + tbc_ket_cc.parity) % 2 != parity_eta)
        continue;
      if (not CcTzCouples(tbc_bra_cc.Tz, tbc_ket_cc.Tz, rank_T))
        continue;
      if (not AngMom::Triangle(tbc_bra_cc.J, tbc_ket_cc.J, lambda))
        continue;
      cc_pairs.push_back({ch_bra_cc, ch_ket_cc});
    }
  }
  const int np = (int)cc_pairs.size();
  CcMatMap bar_CHI_VII_CC, bar_Eta_CC;
  std::vector<arma::mat *> chi_ptr(np, nullptr), eta_ptr(np, nullptr);
  std::vector<int> pair_of((size_t)n_cc * n_cc, -1);
  for (int ip = 0; ip < np; ++ip) {
    TwoBodyChannel_CC &tbc_bra_cc = ms->GetTwoBodyChannel_CC(cc_pairs[ip][0]);
    TwoBodyChannel_CC &tbc_ket_cc = ms->GetTwoBodyChannel_CC(cc_pairs[ip][1]);
    const int nbras_cc = tbc_bra_cc.GetNumberKets();
    const int nkets_cc = tbc_ket_cc.GetNumberKets();
    bar_CHI_VII_CC[cc_pairs[ip]] =
        arma::mat(2 * nbras_cc, 2 * nkets_cc, arma::fill::zeros);
    bar_Eta_CC[cc_pairs[ip]] =
        arma::mat(2 * nbras_cc, 2 * nkets_cc, arma::fill::zeros);
    chi_ptr[ip] = &bar_CHI_VII_CC[cc_pairs[ip]];
    eta_ptr[ip] = &bar_Eta_CC[cc_pairs[ip]];
    pair_of[(size_t)cc_pairs[ip][0] * n_cc + cc_pairs[ip][1]] = ip;
  }

  // Gold Pandya by scatter from 2n ordinary blocks (same 9j as
  // PandyaChiIalb / PandyaOmegaBjak). Walk χ/Ω once; write 2n CC bars.
  std::vector<CcPair> chi_keys;
  std::vector<const arma::mat *> chi_mats;
  for (const auto &kv : CHI_VII.MatEl) {
    chi_keys.push_back(kv.first);
    chi_mats.push_back(&kv.second);
  }
  const int n_chi = (int)chi_keys.size();
  const double t_pandya = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ic = 0; ic < n_chi; ++ic) {
    const int ch_ia = chi_keys[ic][0], ch_lb = chi_keys[ic][1];
    const int Jia = ms->GetTwoBodyChannel(ch_ia).J;
    const int Jlb = ms->GetTwoBodyChannel(ch_lb).J;
    const arma::mat &Chi = *chi_mats[ic];
    const int nia = CHI_VII.NumberKets(ch_ia), nlb = CHI_VII.NumberKets(ch_lb);
    const double hats = HatJ(Jia) * HatJ(Jlb);
    const double phJ = ms->phase(Jia + Jlb);
    for (int ibra = 0; ibra < 2 * nia; ++ibra) {
      auto ia = CHI_VII.Pq(ch_ia, ibra);
      if (ia[0] < 0)
        continue;
      const int i = ia[0], a = ia[1];
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oa = ms->GetOrbit(a);
      const double ji = oi.j2 * 0.5, ja = oa.j2 * 0.5;
      for (int iket = 0; iket < 2 * nlb; ++iket) {
        auto lb = CHI_VII.Pq(ch_lb, iket);
        if (lb[0] < 0)
          continue;
        const double chi = Chi(ibra, iket);
        if (std::abs(chi) < 1e-16)
          continue;
        const int l = lb[0], b = lb[1];
        Orbit &ol = ms->GetOrbit(l);
        Orbit &ob = ms->GetOrbit(b);
        const double jl = ol.j2 * 0.5, jb = ob.j2 * 0.5;
        const int parity_il = (oi.l + ol.l) % 2;
        const int Tz_il = std::abs(oi.tz2 - ol.tz2) / 2;
        const int parity_ab = (oa.l + ob.l) % 2;
        const int Tz_ab = std::abs(oa.tz2 - ob.tz2) / 2;
        const int Jilmin = std::abs(oi.j2 - ol.j2) / 2;
        const int Jilmax = (oi.j2 + ol.j2) / 2;
        for (int Jil = Jilmin; Jil <= Jilmax; ++Jil) {
          const size_t ch_il =
              ms->GetTwoBodyChannelIndex(Jil, parity_il, Tz_il);
          if ((int)ch_il >= n_cc)
            continue;
          TwoBodyChannel_CC &tcc_il = ms->GetTwoBodyChannel_CC(ch_il);
          const int row = IndexCC(tcc_il, i, l);
          if (row < 0)
            continue;
          const int Jabmin =
              std::max(std::abs(oa.j2 - ob.j2) / 2, std::abs(Jil - lambda));
          const int Jabmax = std::min((oa.j2 + ob.j2) / 2, Jil + lambda);
          for (int Jab = Jabmin; Jab <= Jabmax; ++Jab) {
            const size_t ch_ab =
                ms->GetTwoBodyChannelIndex(Jab, parity_ab, Tz_ab);
            if ((int)ch_ab >= n_cc)
              continue;
            const int ip = pair_of[(size_t)ch_il * n_cc + ch_ab];
            if (ip < 0)
              continue;
            const int col = IndexCC(ms->GetTwoBodyChannel_CC(ch_ab), a, b);
            if (col < 0)
              continue;
            const double ninej =
                ms->GetNineJ(lambda, Jil, Jab, Jlb, jl, jb, Jia, ji, ja);
            if (std::abs(ninej) < 1e-14)
              continue;
            const double pref =
                -ms->phase(Jil + Jab + lambda +
                           (oi.j2 + oa.j2 + ol.j2 + ob.j2) / 2) *
                HatJ(Jil) * HatJ(Jab);
            const double add = pref * phJ * hats * ninej * chi;
            double *dst = &(*chi_ptr[ip])(row, col);
#pragma omp atomic
            *dst += add;
          }
        }
      }
    }
  }

  std::vector<arma::mat> Om2n((size_t)nch_pairs);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < nch_pairs; ++ip)
    Om2n[ip] = Omega2n(Eta, ms, (int)ch_bra_list[ip], (int)ch_ket_list[ip]);
#pragma omp parallel for schedule(dynamic, 1)
  for (int iom = 0; iom < nch_pairs; ++iom) {
    const int ch_bj = (int)ch_bra_list[iom], ch_ak = (int)ch_ket_list[iom];
    TwoBodyChannel &tbc_bj = ms->GetTwoBodyChannel(ch_bj);
    TwoBodyChannel &tbc_ak = ms->GetTwoBodyChannel(ch_ak);
    const int Jbj = tbc_bj.J, Jak = tbc_ak.J;
    const arma::mat &Om = Om2n[iom];
    const int nbj = tbc_bj.GetNumberKets(), nak = tbc_ak.GetNumberKets();
    const double hats = HatJ(Jbj) * HatJ(Jak);
    const double phJ = ms->phase(Jbj + Jak);
    for (int ibra = 0; ibra < 2 * nbj; ++ibra) {
      auto bj = PqTB(tbc_bj, ibra, nbj);
      if (bj[0] < 0)
        continue;
      const int b = bj[0], j = bj[1];
      Orbit &ob = ms->GetOrbit(b);
      Orbit &oj = ms->GetOrbit(j);
      const double jb = ob.j2 * 0.5, jj = oj.j2 * 0.5;
      for (int iket = 0; iket < 2 * nak; ++iket) {
        auto ak = PqTB(tbc_ak, iket, nak);
        if (ak[0] < 0)
          continue;
        const double om = Om(ibra, iket);
        if (std::abs(om) < 1e-16)
          continue;
        const int a = ak[0], k = ak[1];
        Orbit &oa = ms->GetOrbit(a);
        Orbit &ok = ms->GetOrbit(k);
        const double ja = oa.j2 * 0.5, jk = ok.j2 * 0.5;
        // Ω̄(a,b; k,j) from Ω(b,j; a,k)
        const int parity_ab = (oa.l + ob.l) % 2;
        const int Tz_ab = std::abs(oa.tz2 - ob.tz2) / 2;
        const int parity_kj = (ok.l + oj.l) % 2;
        const int Tz_kj = std::abs(ok.tz2 - oj.tz2) / 2;
        const int Jabmin = std::abs(oa.j2 - ob.j2) / 2;
        const int Jabmax = (oa.j2 + ob.j2) / 2;
        for (int Jab = Jabmin; Jab <= Jabmax; ++Jab) {
          const size_t ch_ab =
              ms->GetTwoBodyChannelIndex(Jab, parity_ab, Tz_ab);
          if ((int)ch_ab >= n_cc)
            continue;
          TwoBodyChannel_CC &tcc_ab = ms->GetTwoBodyChannel_CC(ch_ab);
          const int row = IndexCC(tcc_ab, a, b);
          if (row < 0)
            continue;
          const int Jkjmin =
              std::max(std::abs(ok.j2 - oj.j2) / 2, std::abs(Jab - lambda));
          const int Jkjmax = std::min((ok.j2 + oj.j2) / 2, Jab + lambda);
          for (int Jkj = Jkjmin; Jkj <= Jkjmax; ++Jkj) {
            const size_t ch_kj =
                ms->GetTwoBodyChannelIndex(Jkj, parity_kj, Tz_kj);
            if ((int)ch_kj >= n_cc)
              continue;
            const int ip = pair_of[(size_t)ch_ab * n_cc + ch_kj];
            if (ip < 0)
              continue;
            const int col = IndexCC(ms->GetTwoBodyChannel_CC(ch_kj), k, j);
            if (col < 0)
              continue;
            const double ninej =
                ms->GetNineJ(lambda, Jab, Jkj, Jak, ja, jk, Jbj, jb, jj);
            if (std::abs(ninej) < 1e-14)
              continue;
            const double pref = -ms->phase(lambda) * HatJ(Jab) * HatJ(Jkj);
            const double add = pref * phJ * hats * ninej * om;
            double *dst = &(*eta_ptr[ip])(row, col);
#pragma omp atomic
            *dst += add;
          }
        }
      }
    }
  }
  CHI_VII.MatEl.clear();
  Z.profiler.timer["_GIVc_pandya"] += omp_get_wtime() - t_pandya;

  // CHI_VII_final[ch_il] = Σ_{ch_ab} (−1)^{J_il}/Ĵ_il λ̂⁻¹ (−1)^{J_ab+λ}
  //                        bar_CHI_VII_CC[{ch_il,ch_ab}] · bar_Eta_CC[{ch_ab,ch_il}]
  std::deque<arma::mat> CHI_VII_final(n_cc);
  for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc) {
    const int nkets_cc = ms->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
    CHI_VII_final[ch_cc] =
        arma::mat(2 * nkets_cc, 2 * nkets_cc, arma::fill::zeros);
  }
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_il = 0; ch_il < n_cc; ++ch_il) {
    TwoBodyChannel_CC &tbc_il = ms->GetTwoBodyChannel_CC(ch_il);
    if (tbc_il.GetNumberKets() < 1)
      continue;
    // mid overall (−1)^{Jp}/Ĵp  (same channel factor as tts_ring Path B)
    const double wL = ms->phase(tbc_il.J) / hat(tbc_il.J);
    for (int ch_ab = 0; ch_ab < n_cc; ++ch_ab) {
      TwoBodyChannel_CC &tbc_ab = ms->GetTwoBodyChannel_CC(ch_ab);
      if (tbc_ab.GetNumberKets() < 1)
        continue;
      if (not AngMom::Triangle(tbc_il.J, tbc_ab.J, lambda))
        continue;
      auto itC = bar_CHI_VII_CC.find({ch_il, ch_ab});
      auto itO = bar_Eta_CC.find({ch_ab, ch_il});
      if (itC == bar_CHI_VII_CC.end() or itO == bar_Eta_CC.end())
        continue;
      const double fac = wL * hat_lam_inv * ms->phase(tbc_ab.J + lambda);
      CHI_VII_final[ch_il] += fac * (itC->second * itO->second);
    }
  }
  bar_CHI_VII_CC.clear();
  bar_Eta_CC.clear();

  auto inv_red = [&](index_t i, index_t j, index_t k, index_t l,
                     int J0) -> double {
    // CHI_VII_final on (i,l)×(k,j) from the χ_ialb×Ω_bjak product.
    // Corrected inv: +Ĵ0 Σ Ĵp sixj CHI_VII_final (drop AMC-sample minus).
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
      if (indx_il >= (int)CHI_VII_final[ch_cc].n_rows or
          indx_kj >= (int)CHI_VII_final[ch_cc].n_cols)
        continue;
      tot += hat(Jp) * sixj * CHI_VII_final[ch_cc](indx_il, indx_kj);
    }
    return hat(J0) * tot;
  };

  // T1/T2-only leftover is not Hermitian (T1† ≡ T2). Full-matrix NonHerm
  // keeps fold(p,g,q,h) ≠ fold(q,h,p,g). Both χ terms → Hermitian, upper OK.
  // W2n holds inv_red once; (1−Pij)(1−Pkl) is the GIVb leftover fold (×½).
  const bool nonherm = (givc_chi_which != 0);
  auto &Z2p = Z.TwoBody;
  ch_bra_list.clear();
  ch_ket_list.clear();
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
    const int nbras = (int)tbc_bra.GetNumberKets();
    const int nkets = (int)tbc_ket.GetNumberKets();
    arma::mat W2n(2 * nbras, 2 * nkets, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * nbras; ++ibra) {
      auto ij = PqTB(tbc_bra, ibra, nbras);
      if (ij[0] < 0)
        continue;
      for (int iket = 0; iket < 2 * nkets; ++iket) {
        auto kl = PqTB(tbc_ket, iket, nkets);
        if (kl[0] < 0)
          continue;
        W2n(ibra, iket) = inv_red(ij[0], ij[1], kl[0], kl[1], J0);
      }
    }
    for (int ibra = 0; ibra < nbras; ++ibra) {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t p = bra.p, g = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(p);
      Orbit &oj = Z.modelspace->GetOrbit(g);
      const int ibra_s = (p == g) ? ibra : ibra + nbras;
      const double pij = AngMom::phase((oi.j2 + oj.j2) / 2 - J0);
      const int ketmin = (nonherm or ch_bra != ch_ket) ? 0 : ibra;
      for (int iket = ketmin; iket < nkets; ++iket) {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t q = ket.p, h = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(q);
        Orbit &ol = Z.modelspace->GetOrbit(h);
        const int iket_s = (q == h) ? iket : iket + nkets;
        const double pkl = AngMom::phase((ok.j2 + ol.j2) / 2 - J0);
        double z = W2n(ibra, iket);
        if (givc_fold_as)
          z = 0.5 * (z - pkl * W2n(ibra, iket_s) - pij * W2n(ibra_s, iket) +
                     pij * pkl * W2n(ibra_s, iket_s));
        z /= hat_J0;
        if (p == g)
          z /= PhysConst::SQRT2;
        if (q == h)
          z /= PhysConst::SQRT2;
        if (nonherm)
          Z2p.AddToTBMENonHerm(ch_bra, ch_ket, ibra, iket, z);
        else
          Z2p.AddToTBME(ch_bra, ch_ket, ibra, iket, z);
      }
    }
  }

  Z.profiler.timer["comm223_232_GIVc_pathB"] += omp_get_wtime() - t_start;
}


////////////////////////////////////////////////////////////////////////////
void comm223_232_GIVc(const Operator &Eta_in, const Operator &Gamma,
                      Operator &Z) {
  // Γ^{IV_c}: gold Pandya (ialb/bjak) → DGEMM → inv (all π).
  const ReducedEtaView E(Eta_in);
  comm223_232_GIVc_pathB(E.ref, Gamma, Z);
}

////////////////////////////////////////////////////////////////////////////
/// Factorized 223_132 tts: tensor ω × scalar Γ → scalar leftover 2b.
/// [ω₁, [ω₂, Γ]_3]_{2b}: every diagram is linear in the occupation-weighted
/// one-body ω₁ inserted between the tensor Ω = ω₂ and the scalar Γ,
///   W_ba(a,b) = [(1−n_a)n_b − n_a(1−n_b)] ω₁(b,a),  W_ab likewise with ω₁(a,b).
/// FDC.cc comm223_132 / comm223_132_cross (λ=0) written as channel DGEMMs
/// with W embedded in pair space ("W ⊗ 1" matrices carrying the 6j):
///   ladder : Z  = Γ·M₁·Ω − Ω·M₂·Γ                       (ch0 × ch2 blocks)
///   onebody: χ  (scalar 1b, Ω×W) → Z = [χ⊗1, Γ];  θ (rank-λ 1b, Γ×W) →
///            Z ∋ −Ω·N_B(θ) + N_B2(θ)·Ω  with (1−P_kl), (1−P_ij)
///   cross  : X̂ = Ω·M_X1 + M_X2·Ω (scalar non-AS 2b, ChiOp2n) → Pandya X̄ →
///            X̄·Γ̄ per CC channel (cached Γ̄) → inverse Pandya, (1−P)(1−P)
/// AMC seeds: learn/amc_tts/comm_tts/output/comm223_132tts_*_seed.tex.
/// The 9j of the cross seed is the Pandya chain: {ji jc J2; jj J4 jb; J0 jl jk}
///   = Σ_x (2x+1) {ji jk x; jb jc J2}{jb jc x; jl jj J4}{ji jk x; jl jj J0}.
/// Leftover packaging: unreduced Z, convention-2 λ̂^{-1}.
////////////////////////////////////////////////////////////////////////////
namespace {

double occ_ab(double na, double nb) {
  return (1.0 - na) * nb - na * (1.0 - nb);
}

// Occupation-weighted ω₁ on all orbits, plus the sparsity lists.
struct Eta1Occ132 {
  arma::mat W_ba; // occ(a,b) ω₁(b,a) — ladder / onebody
  arma::mat W_ab; // occ(a,b) ω₁(a,b) — cross
  std::vector<std::vector<int>> b_of_a, a_of_b; // W(a,b) ≠ 0
};

Eta1Occ132 MakeEta1Occ132(const Operator &Eta, ModelSpace *ms) {
  Eta1Occ132 W;
  const int lambda = Eta.GetJRank();
  int n_alloc = 0;
  for (auto o : ms->all_orbits)
    n_alloc = std::max(n_alloc, (int)o + 1);
  W.W_ba.zeros(n_alloc, n_alloc);
  W.W_ab.zeros(n_alloc, n_alloc);
  W.b_of_a.assign(n_alloc, {});
  W.a_of_b.assign(n_alloc, {});
  for (auto a : ms->all_orbits) {
    Orbit &oa = ms->GetOrbit(a);
    for (auto b : ms->all_orbits) {
      Orbit &ob = ms->GetOrbit(b);
      if (not AngMom::Triangle(oa.j2, ob.j2, 2 * lambda))
        continue;
      const double occ = occ_ab(oa.occ, ob.occ);
      if (std::abs(occ) < 1e-12)
        continue;
      const double eba = Eta.OneBody(b, a), eab = Eta.OneBody(a, b);
      if (std::abs(eba) < 1e-14 and std::abs(eab) < 1e-14)
        continue;
      W.W_ba(a, b) = occ * eba;
      W.W_ab(a, b) = occ * eab;
      W.b_of_a[a].push_back((int)b);
      W.a_of_b[b].push_back((int)a);
    }
  }
  return W;
}

// partners[ch] = ordinary channels coupled to ch by Ω^{λ,T,π} (both orders).
std::vector<std::vector<int>> OmegaPartners(const Operator &Eta,
                                            ModelSpace *ms) {
  std::vector<size_t> ch_bra_list, ch_ket_list;
  TensorChannelPairs(Eta, ms, ch_bra_list, ch_ket_list);
  std::vector<std::vector<int>> partners(ms->GetNumberTwoBodyChannels());
  for (size_t ip = 0; ip < ch_bra_list.size(); ++ip)
    partners[ch_bra_list[ip]].push_back((int)ch_ket_list[ip]);
  return partners;
}

// Γ^J(ij;kl) on ordered pairs, every ordinary channel (FDC.cc Gamma blocks).
void FillScalar2n(const Operator &Gamma, ModelSpace *ms,
                  std::deque<arma::mat> &Gam2n) {
  const int nch = ms->GetNumberTwoBodyChannels();
  Gam2n.assign(nch, arma::mat());
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ++ch)
    Gam2n[ch] = Omega2n(Gamma, ms, ch, ch);
}

// Element of a 2n-layout matrix at the swapped ket (l,k) of stored (k,l);
// for k==l the swapped slot does not exist and the value is the same.
inline double SwapKet(const arma::mat &M, int ibra, int iket, int n,
                      bool same) {
  return same ? M(ibra, iket) : M(ibra, iket + n);
}
inline double SwapBra(const arma::mat &M, int ibra, int iket, int n,
                      bool same) {
  return same ? M(ibra, iket) : M(ibra + n, iket);
}

// Small thread-local table for the 6j {J J' λ; ja jb jc} shapes that are not
// in ModelSpace::SixJList (three integer, three half-integer arguments).
struct SixJLocal {
  int max_j2 = 1;
  std::vector<double> tab;
  std::vector<char> have;
  SixJLocal(int mj2) : max_j2(mj2) {
    const size_t n = (size_t)(mj2 + 1);
    tab.assign(n * n * n, 0.0);
    have.assign(n * n * n, 0);
  }
  // {J1 J2 λ; ja jb jc} with (J1,J2,λ) fixed by the caller's channels.
  double get(int J1, int J2, int lambda, int j2a, int j2b, int j2c) {
    const size_t n = (size_t)(max_j2 + 1);
    const size_t key = ((size_t)j2a * n + (size_t)j2b) * n + (size_t)j2c;
    if (!have[key]) {
      tab[key] = AngMom::SixJ((double)J1, (double)J2, (double)lambda,
                              0.5 * j2a, 0.5 * j2b, 0.5 * j2c);
      have[key] = 1;
    }
    return tab[key];
  }
};

int MaxJ2(ModelSpace *ms) {
  int m = 1;
  for (auto o : ms->all_orbits)
    m = std::max(m, ms->GetOrbit(o).j2);
  return m;
}

} // namespace

////////////////////////////////////////////////////////////////////////////
/// 132 ladder.  AMC seed comm223_132tts_ladder_seed.tex:
///   T1:  Ĵ0^{-1} (−1)^{J2+ja+jc} Ĵ2 λ̂^{-1} {J0 J2 λ; ja jb jc} W_ba Ω_cakl^{J2J0} Γ_ijcb^{J0}
///   T2: −Ĵ0^{-1} (−1)^{J0+ja+jc} Ĵ2 λ̂^{-1} {J0 J2 λ; jb ja jc} W_ba Ω_ijcb^{J0J2} Γ_cakl^{J2}
/// Per Z channel ch0 and Ω partner ch2 (FDC.cc comm223_132 "full matrix"):
///   Z(ch0) += Γ(ch0)·M₁(ch0,ch2)·Ω(ch2,ch0) − Ω(ch0,ch2)·M₂(ch0,ch2)·Γ(ch2)
/// with M[(c,b),(c,a)] = the 6j-dressed W_ba(a,b) (one-body ⊗ spectator c).
////////////////////////////////////////////////////////////////////////////
void comm223_132_tts_ladder(const Operator &Eta_in, const Operator &Gamma,
                            Operator &Z) {
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
  double t_start = omp_get_wtime();
  ModelSpace *ms = Z.modelspace;
  ms->PreCalculateSixJ();
  const int lambda = Eta.GetJRank();
  const double hat_lam_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  const Eta1Occ132 W = MakeEta1Occ132(Eta, ms);
  const std::vector<std::vector<int>> partners = OmegaPartners(Eta, ms);
  std::deque<arma::mat> Gam2n;
  FillScalar2n(Gamma, ms, Gam2n);
  const int max_j2 = MaxJ2(ms);

  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ich++) {
    const size_t ch0 = ch_bra_list[ich];
    if (ch_ket_list[ich] != ch0)
      continue;
    TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
    const int n0 = tbc0.GetNumberKets();
    if (n0 < 1)
      continue;
    const int J0 = tbc0.J;
    const double hatJ0inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
    SixJLocal six1(max_j2), six2(max_j2);
    arma::mat Zacc(2 * n0, 2 * n0, arma::fill::zeros);

    for (int ch2 : partners[ch0]) {
      TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
      const int n2 = tbc2.GetNumberKets();
      if (n2 < 1)
        continue;
      const int J2 = tbc2.J;
      const double pref0 = hatJ0inv * HatJ(J2) * hat_lam_inv;
      six1 = SixJLocal(max_j2);
      six2 = SixJLocal(max_j2);
      // T1: Γ^{J0}(ij;cb) M₁[(c,b)∈ch0,(c,a)∈ch2] Ω^{J2J0}(ca;kl)
      arma::mat M1(2 * n0, 2 * n2, arma::fill::zeros);
      bool any1 = false;
      for (int r = 0; r < 2 * n0; ++r) {
        auto cb = PqTB(tbc0, r, n0);
        if (cb[0] < 0)
          continue;
        const int c = cb[0], b = cb[1];
        Orbit &oc = ms->GetOrbit(c);
        Orbit &ob = ms->GetOrbit(b);
        for (int a : W.a_of_b[b]) {
          const double w = W.W_ba(a, b);
          if (std::abs(w) < 1e-16)
            continue;
          const int col = Index2n(tbc2, c, a);
          if (col < 0)
            continue;
          Orbit &oa = ms->GetOrbit(a);
          // (−1)^{J2+ja+jc} {J0 J2 λ; ja jb jc}
          const double s1 = six1.get(J0, J2, lambda, oa.j2, ob.j2, oc.j2);
          if (std::abs(s1) < 1e-12)
            continue;
          M1(r, col) = w * pref0 * ms->phase(J2 + (oa.j2 + oc.j2) / 2) * s1;
          any1 = true;
        }
      }
      if (any1)
        Zacc += Gam2n[ch0] * M1 * Omega2n(Eta, ms, ch2, ch0);

      // T2: Ω^{J0J2}(ij;cb) M₂[(c,b)∈ch2,(c,a)∈ch0] Γ^{J0}(ca;kl)
      arma::mat M2(2 * n2, 2 * n0, arma::fill::zeros);
      bool any2 = false;
      for (int r = 0; r < 2 * n2; ++r) {
        auto cb = PqTB(tbc2, r, n2);
        if (cb[0] < 0)
          continue;
        const int c = cb[0], b = cb[1];
        Orbit &oc = ms->GetOrbit(c);
        Orbit &ob = ms->GetOrbit(b);
        for (int a : W.a_of_b[b]) {
          const double w = W.W_ba(a, b);
          if (std::abs(w) < 1e-16)
            continue;
          const int col = Index2n(tbc0, c, a);
          if (col < 0)
            continue;
          Orbit &oa = ms->GetOrbit(a);
          // (−1)^{J0+ja+jc} {J0 J2 λ; jb ja jc}
          const double s2 = six2.get(J0, J2, lambda, ob.j2, oa.j2, oc.j2);
          if (std::abs(s2) < 1e-12)
            continue;
          M2(r, col) = w * pref0 * ms->phase(J0 + (oa.j2 + oc.j2) / 2) * s2;
          any2 = true;
        }
      }
      if (any2)
        Zacc -= Omega2n(Eta, ms, ch0, ch2) * M2 * Gam2n[ch0];
    }

    for (int ibra = 0; ibra < n0; ++ibra) {
      Ket &bra = tbc0.GetKet(ibra);
      for (int iket = ibra; iket < n0; ++iket) {
        Ket &ket = tbc0.GetKet(iket);
        double zijkl = Zacc(ibra, iket);
        if (bra.p == bra.q)
          zijkl /= PhysConst::SQRT2;
        if (ket.p == ket.q)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch0, ch0, ibra, iket, zijkl);
      }
    }
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// 132 onebody.  AMC seeds comm223_132tts_onebodyA/B/B2_seed.tex.
///  1bA: χ(p,q) = −(−1)^{jp} (2jp+1)^{-1} Σ_{ab J2 J3} W_ba (−1)^{J2+ja} Ĵ2 Ĵ3 λ̂^{-1}
///                 {J3 λ J2; ja jp jb} Ω^{J2J3}(pa;qb)          (scalar 1b)
///       Z = (χ⊗1)·Γ − Γ·(χ⊗1)                                  (FDC.cc CHI_I)
///  θ(p,q) = Σ_{ab J} W_ba (−1)^{J+jb} (2J+1) {jq jp λ; ja jb J} Γ^J(pa;qb)  (rank λ)
///  1bB : Z −= Ω(ch0,ch2)·N_B,  N_B[(k,c),(k,l)] = −(−1)^{J0+jk} Ĵ0^{-1}Ĵ2 λ̂^{-1}
///                                                   {J2 λ J0; jl jk jc} θ(c,l)
///  1bB2: Z += N_B2·Ω(ch3,ch0), N_B2[(i,j),(i,c)] = −(−1)^{ji+J3} Ĵ0^{-1}Ĵ3 λ̂^{-1}
///                                                   {J0 J3 λ; jc jj ji} θ(j,c)
///  then (1−P_kl) on 1bB and (1−P_ij) on 1bB2 (seeds are one-sided).
////////////////////////////////////////////////////////////////////////////
void comm223_132_tts_onebody(const Operator &Eta_in, const Operator &Gamma,
                              Operator &Z) {
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
  double t_start = omp_get_wtime();
  ModelSpace *ms = Z.modelspace;
  ms->PreCalculateSixJ();
  const int lambda = Eta.GetJRank();
  const double hat_lam_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  const Eta1Occ132 W = MakeEta1Occ132(Eta, ms);
  const std::vector<std::vector<int>> partners = OmegaPartners(Eta, ms);
  const int nch_all = ms->GetNumberTwoBodyChannels();
  const int max_j2 = MaxJ2(ms);
  const int n_alloc = (int)W.W_ba.n_rows;
  std::deque<arma::mat> Gam2n;
  FillScalar2n(Gamma, ms, Gam2n);

  // χ (scalar) from Ω ⊗ W over Ω channel pairs; θ (rank λ) from Γ ⊗ W over
  // scalar channels. Both are traces over the ordered-pair slots.
  arma::mat CHI(n_alloc, n_alloc, arma::fill::zeros);
  arma::mat THETA(n_alloc, n_alloc, arma::fill::zeros);
#pragma omp parallel
  {
    arma::mat chi_loc(n_alloc, n_alloc, arma::fill::zeros);
    arma::mat theta_loc(n_alloc, n_alloc, arma::fill::zeros);
    SixJLocal six(max_j2);
#pragma omp for schedule(dynamic, 1)
    for (int ch2 = 0; ch2 < nch_all; ++ch2) {
      TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
      const int n2 = tbc2.GetNumberKets();
      if (n2 < 1)
        continue;
      const int J2 = tbc2.J;
      // χ: Ω^{J2 J3}(pa;qb), p and q in the same one-body channel.
      for (int ch3 : partners[ch2]) {
        TwoBodyChannel &tbc3 = ms->GetTwoBodyChannel(ch3);
        const int n3 = tbc3.GetNumberKets();
        if (n3 < 1)
          continue;
        const int J3 = tbc3.J;
        const arma::mat Om = Omega2n(Eta, ms, ch2, ch3);
        const double pref = HatJ(J2) * HatJ(J3) * hat_lam_inv;
        six = SixJLocal(max_j2); // {J3 λ J2; ja jp jb}
        for (int r = 0; r < 2 * n2; ++r) {
          auto pa = PqTB(tbc2, r, n2);
          if (pa[0] < 0)
            continue;
          const int p = pa[0], a = pa[1];
          Orbit &op = ms->GetOrbit(p);
          Orbit &oa = ms->GetOrbit(a);
          const double fac_p = -ms->phase(op.j2 / 2) / (op.j2 + 1.0) *
                               ms->phase(J2 + oa.j2 / 2) * pref;
          for (int b : W.b_of_a[a]) {
            const double w = W.W_ba(a, b);
            if (std::abs(w) < 1e-16)
              continue;
            Orbit &ob = ms->GetOrbit(b);
            const double sixj = six.get(J3, lambda, J2, oa.j2, op.j2, ob.j2);
            if (std::abs(sixj) < 1e-12)
              continue;
            for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
              const int col = Index2n(tbc3, (int)q, b);
              if (col < 0)
                continue;
              chi_loc(p, q) += fac_p * w * sixj * Om(r, col);
            }
          }
        }
      }
      // θ: Γ^{J2}(pa;qb), (−1)^{J2+jb} (2J2+1) {jq jp λ; ja jb J2} W_ba(a,b)
      const arma::mat &G = Gam2n[ch2];
      for (int r = 0; r < 2 * n2; ++r) {
        auto pa = PqTB(tbc2, r, n2);
        if (pa[0] < 0)
          continue;
        const int p = pa[0], a = pa[1];
        Orbit &op = ms->GetOrbit(p);
        Orbit &oa = ms->GetOrbit(a);
        for (int b : W.b_of_a[a]) {
          const double w = W.W_ba(a, b);
          if (std::abs(w) < 1e-16)
            continue;
          Orbit &ob = ms->GetOrbit(b);
          const double fac = w * ms->phase(J2 + ob.j2 / 2) * (2.0 * J2 + 1.0);
          for (auto q : ms->all_orbits) {
            Orbit &oq = ms->GetOrbit(q);
            if (not AngMom::Triangle(oq.j2, op.j2, 2 * lambda))
              continue;
            const int c = Index2n(tbc2, (int)q, b);
            if (c < 0)
              continue;
            // {jq jp λ; ja jb J2}: cached shape, all four triads hold here.
            const double sixj = ms->GetSixJ(0.5 * oq.j2, 0.5 * op.j2, lambda,
                                            0.5 * oa.j2, 0.5 * ob.j2, J2);
            if (std::abs(sixj) < 1e-12)
              continue;
            theta_loc(p, q) += fac * sixj * G(r, c);
          }
        }
      }
    }
#pragma omp critical
    {
      CHI += chi_loc;
      THETA += theta_loc;
    }
  }
  Z.profiler.timer["_132_onebody_chi"] += omp_get_wtime() - t_start;
  double t_mid = omp_get_wtime();

  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ich++) {
    const size_t ch0 = ch_bra_list[ich];
    if (ch_ket_list[ich] != ch0)
      continue;
    TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
    const int n0 = tbc0.GetNumberKets();
    if (n0 < 1)
      continue;
    const int J0 = tbc0.J;
    const double hatJ0inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
    const arma::mat &G0 = Gam2n[ch0];

    // 1bA: X[(i,j),(a,b)] = δ_jb χ(i,a) + δ_ia χ(j,b);  Z = X·Γ − Γ·X
    arma::mat X(2 * n0, 2 * n0, arma::fill::zeros);
    for (int r = 0; r < 2 * n0; ++r) {
      auto ij = PqTB(tbc0, r, n0);
      if (ij[0] < 0)
        continue;
      const int i = ij[0], j = ij[1];
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      for (auto a : Z.OneBodyChannels.at({oi.l, oi.j2, oi.tz2})) {
        const int col = Index2n(tbc0, (int)a, j);
        if (col >= 0)
          X(r, col) += CHI(i, a);
      }
      for (auto b : Z.OneBodyChannels.at({oj.l, oj.j2, oj.tz2})) {
        const int col = Index2n(tbc0, i, (int)b);
        if (col >= 0)
          X(r, col) += CHI(j, b);
      }
    }
    arma::mat Zacc = X * G0 - G0 * X;

    // 1bB / 1bB2 over Ω partners of ch0.
    arma::mat ZB(2 * n0, 2 * n0, arma::fill::zeros);
    arma::mat ZB2(2 * n0, 2 * n0, arma::fill::zeros);
    SixJLocal sixB(max_j2), sixB2(max_j2);
    for (int ch2 : partners[ch0]) {
      TwoBodyChannel &tbc2 = ms->GetTwoBodyChannel(ch2);
      const int n2 = tbc2.GetNumberKets();
      if (n2 < 1)
        continue;
      const int J2 = tbc2.J;
      const double pref = -hatJ0inv * HatJ(J2) * hat_lam_inv;
      sixB = SixJLocal(max_j2);  // {J2 λ J0; jl jk jc}
      sixB2 = SixJLocal(max_j2); // {J0 J2 λ; jc jj ji}
      // N_B (ch2 rows (k,c), ch0 cols (k,l)):
      //   −(−1)^{J0+jk} Ĵ0^{-1}Ĵ2 λ̂^{-1} {J2 λ J0; jl jk jc} θ(c,l)
      arma::mat NB(2 * n2, 2 * n0, arma::fill::zeros);
      // N_B2 (ch0 rows (i,j), ch2 cols (i,c)):
      //   −(−1)^{ji+J2} Ĵ0^{-1}Ĵ2 λ̂^{-1} {J0 J2 λ; jc jj ji} θ(j,c)
      arma::mat NB2(2 * n0, 2 * n2, arma::fill::zeros);
      bool anyB = false, anyB2 = false;
      for (int r = 0; r < 2 * n0; ++r) {
        auto kl = PqTB(tbc0, r, n0);
        if (kl[0] < 0)
          continue;
        const int k = kl[0], l = kl[1];
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        for (auto c : ms->all_orbits) {
          Orbit &oc = ms->GetOrbit(c);
          // N_B: row (k,c) in ch2, θ(c,l)
          {
            const double th = THETA(c, l);
            const int row = (std::abs(th) > 1e-16) ? Index2n(tbc2, k, (int)c) : -1;
            if (row >= 0) {
              const double sixj = sixB.get(J2, lambda, J0, ol.j2, ok.j2, oc.j2);
              if (std::abs(sixj) > 1e-12) {
                NB(row, r) = pref * ms->phase(J0 + ok.j2 / 2) * sixj * th;
                anyB = true;
              }
            }
          }
          // N_B2: row (i,j)=(k,l) here, col (i,c) in ch2, θ(j,c)
          {
            const double th = THETA(l, c);
            const int col = (std::abs(th) > 1e-16) ? Index2n(tbc2, k, (int)c) : -1;
            if (col >= 0) {
              const double sixj = sixB2.get(J0, J2, lambda, oc.j2, ol.j2, ok.j2);
              if (std::abs(sixj) > 1e-12) {
                NB2(r, col) = pref * ms->phase(ok.j2 / 2 + J2) * sixj * th;
                anyB2 = true;
              }
            }
          }
        }
      }
      if (anyB)
        ZB += Omega2n(Eta, ms, ch0, ch2) * NB;
      if (anyB2)
        ZB2 += NB2 * Omega2n(Eta, ms, ch2, ch0);
    }

    for (int ibra = 0; ibra < n0; ibra++) {
      Ket &bra = tbc0.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const int ph_ij = ms->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = ibra; iket < n0; iket++) {
        Ket &ket = tbc0.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const int ph_kl = ms->phase((ok.j2 + ol.j2) / 2 - J0);
        double zijkl = Zacc(ibra, iket);
        zijkl -= ZB(ibra, iket) - ph_kl * SwapKet(ZB, ibra, iket, n0, k == l);
        zijkl += ZB2(ibra, iket) - ph_ij * SwapBra(ZB2, ibra, iket, n0, i == j);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch0, ch0, ibra, iket, zijkl);
      }
    }
  }
  Z.profiler.timer["_132_onebody_fold"] += omp_get_wtime() - t_mid;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// 132 cross.  AMC seed comm223_132tts_cross_seed.tex (Path A), with the 9j
/// resolved into the Pandya chain (FDC.cc comm223_132_cross for λ=0):
///   X̂^{J2}(ic;kb) = Σ_{a J3} (−1)^{J2+jk+jb} Ĵ2^{-1}Ĵ3 λ̂^{-1} {J2 J3 λ; ja jb jk}
///                    W_ab(a,b) Ω^{J2J3}(ic;ka)                       [T1]
///   X̂^{J3}(ia;kc) −= Σ_{b J2} (−1)^{J2+ji+jb} Ĵ3^{-1}Ĵ2 λ̂^{-1} {J3 J2 λ; jb ja ji}
///                    W_ab(a,b) Ω^{J2J3}(ib;kc)                       [T2]
///   X̄^x(ik;bc) = Σ_J (2J+1) {ji jk x; jb jc J} (−1)^{jk+jb+J} X̂^J(ic;kb)   (Pandya)
///   Ī^x        = X̄^x · Γ̄^x                                (CC DGEMM, cached Γ̄)
///   S(ij;kl)   = −(−1)^{J0+jk−jl} Σ_x (2x+1) {ji jk x; jl jj J0} Ī^x(ik;lj)
///   Z = (1−P_ij)(1−P_kl) S.
/// X̂ is scalar but not antisymmetric (X̂(ia;kc) ≠ 0 for i==a at odd J), so
/// it lives in the no-Pauli ChiOp2n layout.
////////////////////////////////////////////////////////////////////////////
void comm223_132_tts_cross(const Operator &Eta_in, const Operator &Gamma,
                            Operator &Z) {
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;
  double t_start = omp_get_wtime();
  ModelSpace *ms = Z.modelspace;
  ms->PreCalculateSixJ();
  const int lambda = Eta.GetJRank();
  const double hat_lam_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;
  const Eta1Occ132 W = MakeEta1Occ132(Eta, ms);
  const std::vector<std::vector<int>> partners = OmegaPartners(Eta, ms);
  const int nch_all = ms->GetNumberTwoBodyChannels();
  const int n_cc = ms->GetNumberTwoBodyChannels_CC();
  const int max_j2 = MaxJ2(ms);

  // X̂ in ChiOp2n, block diagonal.
  ChiOp2n Xhat;
  Xhat.Init(ms);
  std::vector<size_t> ch_list;
  for (int ch = 0; ch < nch_all; ++ch)
    if (Xhat.NumberKets(ch) > 0)
      ch_list.push_back((size_t)ch);
  Xhat.Allocate(ch_list, ch_list);

#pragma omp parallel for schedule(dynamic, 1)
  for (int chX = 0; chX < nch_all; ++chX) {
    TwoBodyChannel &tbcX = ms->GetTwoBodyChannel(chX);
    const int nX = tbcX.GetNumberKets();
    const int NX = Xhat.NumberKets(chX);
    if (nX < 1 or NX < 1)
      continue;
    const int JX = tbcX.J;
    const double hatJXinv = 1.0 / std::sqrt(2.0 * JX + 1.0);
    arma::mat &XB = *Xhat.GetBlock(chX, chX);
    // Pauli slot → no-Pauli slot of chX
    std::vector<int> np(2 * nX, -1);
    for (int s = 0; s < 2 * nX; ++s) {
      auto pq = PqTB(tbcX, s, nX);
      if (pq[0] >= 0)
        np[s] = Xhat.Index(chX, pq[0], pq[1]);
    }
    SixJLocal six1(max_j2), six2(max_j2);
    for (int chY : partners[chX]) {
      TwoBodyChannel &tbcY = ms->GetTwoBodyChannel(chY);
      const int nY = tbcY.GetNumberKets();
      if (nY < 1)
        continue;
      const int JY = tbcY.J;
      const double pref = hatJXinv * HatJ(JY) * hat_lam_inv;
      six1 = SixJLocal(max_j2); // {JX JY λ; ja jb jk}
      six2 = SixJLocal(max_j2); // {JX JY λ; jb ja ji}

      // T1: X̂(chX)(ic;kb) += Σ Ω^{JX JY}(ic;ka) M1[(k,a),(k,b)]
      arma::mat M1(2 * nY, 2 * NX, arma::fill::zeros);
      bool any1 = false;
      for (int r = 0; r < 2 * nY; ++r) {
        auto ka = PqTB(tbcY, r, nY);
        if (ka[0] < 0)
          continue;
        const int k = ka[0], a = ka[1];
        Orbit &ok = ms->GetOrbit(k);
        Orbit &oa = ms->GetOrbit(a);
        for (int b : W.b_of_a[a]) {
          const double w = W.W_ab(a, b);
          if (std::abs(w) < 1e-16)
            continue;
          const int col = Xhat.Index(chX, k, b);
          if (col < 0)
            continue;
          Orbit &ob = ms->GetOrbit(b);
          const double sixj = six1.get(JX, JY, lambda, oa.j2, ob.j2, ok.j2);
          if (std::abs(sixj) < 1e-12)
            continue;
          M1(r, col) = ms->phase(JX + (ok.j2 + ob.j2) / 2) * pref * sixj * w;
          any1 = true;
        }
      }
      if (any1) {
        const arma::mat R = Omega2n(Eta, ms, chX, chY) * M1; // 2nX × 2NX
        for (int s = 0; s < 2 * nX; ++s)
          if (np[s] >= 0)
            XB.row(np[s]) += R.row(s);
      }

      // T2: X̂(chX)(ia;kc) += Σ M2[(i,a),(i,b)] Ω^{JY JX}(ib;kc)
      arma::mat M2(2 * NX, 2 * nY, arma::fill::zeros);
      bool any2 = false;
      for (int c = 0; c < 2 * nY; ++c) {
        auto ib = PqTB(tbcY, c, nY);
        if (ib[0] < 0)
          continue;
        const int i = ib[0], b = ib[1];
        Orbit &oi = ms->GetOrbit(i);
        Orbit &ob = ms->GetOrbit(b);
        for (int a : W.a_of_b[b]) {
          const double w = W.W_ab(a, b);
          if (std::abs(w) < 1e-16)
            continue;
          const int row = Xhat.Index(chX, i, a);
          if (row < 0)
            continue;
          Orbit &oa = ms->GetOrbit(a);
          const double sixj = six2.get(JX, JY, lambda, ob.j2, oa.j2, oi.j2);
          if (std::abs(sixj) < 1e-12)
            continue;
          M2(row, c) = -ms->phase(JY + (oi.j2 + ob.j2) / 2) * pref * sixj * w;
          any2 = true;
        }
      }
      if (any2) {
        const arma::mat R = M2 * Omega2n(Eta, ms, chY, chX); // 2NX × 2nX
        for (int s = 0; s < 2 * nX; ++s)
          if (np[s] >= 0)
            XB.col(np[s]) += R.col(s);
      }
    }
  }
  Z.profiler.timer["_132_cross_chi"] += omp_get_wtime() - t_start;
  double t_mid = omp_get_wtime();

  // Pandya of X̂ (FDC.cc Eta_bar convention) and DGEMM with the cached Γ̄:
  //   X̄^x(ab;cd) = Σ_J (2J+1) {ja jb x; jc jd J} (−1)^{jb+jc+J} X̂^J(ad;bc)
  const std::deque<arma::mat> &bar_Gamma =
      CachedBarGammaScalarCC(Gamma, Z, hGamma);
  std::deque<arma::mat> bar_I(n_cc);
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch_cc = 0; ch_cc < n_cc; ++ch_cc) {
    TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
    const int ncc = tbc_cc.GetNumberKets();
    if (ncc < 1)
      continue;
    const int x = tbc_cc.J;
    arma::mat bar_X(2 * ncc, 2 * ncc, arma::fill::zeros);
    for (int ibra = 0; ibra < 2 * ncc; ++ibra) {
      auto ab = PqCC(tbc_cc, ibra, ncc);
      if (ab[0] < 0)
        continue;
      const int a = ab[0], b = ab[1];
      Orbit &oa = ms->GetOrbit(a);
      Orbit &ob = ms->GetOrbit(b);
      for (int iket = 0; iket < 2 * ncc; ++iket) {
        auto cd = PqCC(tbc_cc, iket, ncc);
        if (cd[0] < 0)
          continue;
        const int c = cd[0], d = cd[1];
        Orbit &oc = ms->GetOrbit(c);
        Orbit &od = ms->GetOrbit(d);
        if ((oa.l + od.l) % 2 != (ob.l + oc.l) % 2)
          continue;
        if (oa.tz2 + od.tz2 != ob.tz2 + oc.tz2)
          continue;
        const int Jmin =
            std::max(std::abs(oa.j2 - od.j2), std::abs(ob.j2 - oc.j2)) / 2;
        const int Jmax = std::min(oa.j2 + od.j2, ob.j2 + oc.j2) / 2;
        double xbar = 0.0;
        for (int J = Jmin; J <= Jmax; ++J) {
          const int ch = Xhat.ChannelIndex(J, a, d);
          const arma::mat *M = Xhat.GetBlock(ch, ch);
          if (M == nullptr)
            continue;
          const int iad = Xhat.Index(ch, a, d), ibc = Xhat.Index(ch, b, c);
          if (iad < 0 or ibc < 0)
            continue;
          const double sixj =
              ms->GetSixJ(0.5 * oa.j2, 0.5 * ob.j2, x, 0.5 * oc.j2,
                          0.5 * od.j2, J);
          if (std::abs(sixj) < 1e-12)
            continue;
          xbar += (2.0 * J + 1.0) * sixj *
                  ms->phase((ob.j2 + oc.j2) / 2 + J) * (*M)(iad, ibc);
        }
        bar_X(ibra, iket) = xbar;
      }
    }
    bar_I[ch_cc] = bar_X * bar_Gamma[ch_cc];
  }
  Xhat.MatEl.clear();
  Z.profiler.timer["_132_cross_pandya"] += omp_get_wtime() - t_mid;
  t_mid = omp_get_wtime();

  // Inverse: S(ij;kl) = −(−1)^{J0+jk−jl} Σ_x (2x+1) {ji jk x; jl jj J0} Ī^x(ik;lj)
  auto seed = [&](int i, int j, int k, int l, int J0) -> double {
    Orbit &oi = ms->GetOrbit(i);
    Orbit &oj = ms->GetOrbit(j);
    Orbit &ok = ms->GetOrbit(k);
    Orbit &ol = ms->GetOrbit(l);
    const int parity_cc = (oi.l + ok.l) % 2;
    const int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
    const int xmin =
        std::max(std::abs(oi.j2 - ok.j2), std::abs(oj.j2 - ol.j2)) / 2;
    const int xmax = std::min(oi.j2 + ok.j2, oj.j2 + ol.j2) / 2;
    double s = 0.0;
    for (int x = xmin; x <= xmax; ++x) {
      const double sixj = ms->GetSixJ(0.5 * oj.j2, 0.5 * oi.j2, J0, 0.5 * ok.j2,
                                      0.5 * ol.j2, x);
      if (std::abs(sixj) < 1e-12)
        continue;
      const size_t ch_cc = ms->GetTwoBodyChannelIndex(x, parity_cc, Tz_cc);
      if ((int)ch_cc >= n_cc or bar_I[ch_cc].n_rows < 1)
        continue;
      TwoBodyChannel_CC &tbc_cc = ms->GetTwoBodyChannel_CC(ch_cc);
      const int ncc = tbc_cc.GetNumberKets();
      int iik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
      int ilj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
      if (iik < 0 or ilj < 0)
        continue;
      iik += (i > k) ? ncc : 0;
      ilj += (l > j) ? ncc : 0;
      s += (2.0 * x + 1.0) * sixj * bar_I[ch_cc](iik, ilj);
    }
    return -ms->phase(J0 + (ok.j2 - ol.j2) / 2) * s;
  };

  auto &Z2 = Z.TwoBody;
  std::vector<size_t> ch_bra_list, ch_ket_list;
  for (auto &iter : Z2.MatEl) {
    ch_bra_list.push_back(iter.first[0]);
    ch_ket_list.push_back(iter.first[1]);
  }
  const int nch = (int)ch_bra_list.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ich++) {
    const size_t ch0 = ch_bra_list[ich];
    if (ch_ket_list[ich] != ch0)
      continue;
    TwoBodyChannel &tbc0 = ms->GetTwoBodyChannel(ch0);
    const int J0 = tbc0.J;
    const int n0 = tbc0.GetNumberKets();
    for (int ibra = 0; ibra < n0; ibra++) {
      Ket &bra = tbc0.GetKet(ibra);
      const int i = (int)bra.p, j = (int)bra.q;
      Orbit &oi = ms->GetOrbit(i);
      Orbit &oj = ms->GetOrbit(j);
      const int ph_ij = ms->phase((oi.j2 + oj.j2) / 2 - J0);
      for (int iket = ibra; iket < n0; iket++) {
        Ket &ket = tbc0.GetKet(iket);
        const int k = (int)ket.p, l = (int)ket.q;
        Orbit &ok = ms->GetOrbit(k);
        Orbit &ol = ms->GetOrbit(l);
        const int ph_kl = ms->phase((ok.j2 + ol.j2) / 2 - J0);
        const int ph_ijkl = ms->phase((oi.j2 + oj.j2 + ok.j2 + ol.j2) / 2);
        double zijkl = seed(i, j, k, l, J0) - ph_ij * seed(j, i, k, l, J0) -
                       ph_kl * seed(i, j, l, k, J0) +
                       ph_ijkl * seed(j, i, l, k, J0);
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch0, ch0, ibra, iket, zijkl);
      }
    }
  }
  Z.profiler.timer["_132_cross_inv"] += omp_get_wtime() - t_mid;
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_132_tts(const Operator &Eta_in, const Operator &Gamma,
                     Operator &Z) {
  const double t_start = omp_get_wtime();
  if (Z.GetJRank() != 0)
    return;
  const ReducedEtaView E(Eta_in);
  const Operator &Eta = E.ref;

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
