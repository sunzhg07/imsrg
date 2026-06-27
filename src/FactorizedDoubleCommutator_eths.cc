
#include "FactorizedDoubleCommutator_eths.hh"
#include "Commutator.hh"
#include "ReferenceImplementations.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"

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
bool use_TypeII_1b = true;
bool use_TypeIII_1b = true;
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

void SetUse_TypeII_1b(bool tf) { use_TypeII_1b = tf; }

void SetUse_TypeIII_1b(bool tf) { use_TypeIII_1b = tf; }

void SetUse_TypeII_2b(bool tf) { use_TypeII_2b = tf; }

void SetUse_TypeIII_2b(bool tf) { use_TypeIII_2b = tf; }

void SetUse_GT_TypeI_2b(bool tf) { use_GT_TypeI_2b = tf; }

void SetUse_GT_TypeIV_2b(bool tf) { use_GT_TypeIV_2b = tf; }

//  void UseSlowVersion(bool tf)
//  {
//    SlowVersion = tf;
//  }

// factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
void comm223_231_st(const Operator &Eta, const Operator &Gamma, Operator &Z) {

  // Do some extra work to check if the operators being passed in are reduced,
  // because the comm223_231 routines assume not-reduced matrix elements
  const Operator *Etanred =
      &Eta; // Pointer to the non-reduced version of the operator
  const Operator *Gammanred = &Gamma;
  Operator Etatmp,
      Gammatmp;        // Declared, but not yet allocated, for reasons of scope.
  if (Eta.IsReduced()) // comm223_231 doesn't expect reduced operators. Need to
                       // make it not reduced.
  {
    Etatmp = Eta; // Actual copy, not a reference
    Etatmp.MakeNotReduced();
    Etanred = &Etatmp; // Now Etanred points to the not-reduced copy Etatmp
  }
  if (Gamma.IsReduced()) {
    Gammatmp = Gamma;
    Gammatmp.MakeNotReduced();
    Gammanred =
        &Gammatmp; // Now Gammanred points to the not-reduced copy Gammatmp
  }

  bool z_was_reduced = Z.IsReduced();
  if (z_was_reduced)
    Z.MakeNotReduced();

  if (use_1b_intermediates)
    comm223_231_chi1b_tensor(*Etanred, *Gammanred, Z);
  if (use_2b_intermediates)
    comm223_231_chi2b_tensor(*Etanred, *Gammanred, Z);
  if (z_was_reduced)
    Z.MakeReduced();

  return;
} // comm223_231_st

////////////////////////////////////////////////////////////////////////////
/// tensor extension of 223_231 one-body block (function-by-function)
/// Implements diag2 term-I structure: chi^alpha -> f^(I), no explicit 3b.
////////////////////////////////////////////////////////////////////////////
void comm223_231_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {

  double t_internal = omp_get_wtime();
  double t_start = omp_get_wtime();

  Z.modelspace->PreCalculateSixJ();

  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  int lambda = Eta.GetJRank();

  auto Chi_alpha = Z.OneBody;
  auto Chi_beta = Z.OneBody;
  Chi_alpha.zeros();
  Chi_beta.zeros();

  int norbits = Z.modelspace->all_orbits.size();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),
                                  Z.modelspace->all_orbits.end());

  double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);

  // chi^alpha_{de} from diag2_compact (scalar intermediate from Omega*Omega).
#pragma omp parallel for schedule(dynamic, 1)
  for (int indexd = 0; indexd < norbits; ++indexd) {
    auto d = allorb_vec[indexd];
    Orbit &od = Z.modelspace->GetOrbit(d);
    for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) {
      if (e > d)
        continue;

      double chi_de = 0.0;
      for (auto &a : Z.modelspace->all_orbits) {
        Orbit &oa = Z.modelspace->GetOrbit(a);
        double n_a = oa.occ;
        double nbar_a = 1.0 - n_a;

        for (auto &b : Z.modelspace->all_orbits) {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_b = ob.occ;
          double nbar_b = 1.0 - n_b;

          int Jab_min = std::abs(oa.j2 - ob.j2) / 2;
          int Jab_max = (oa.j2 + ob.j2) / 2;

          for (auto &c : Z.modelspace->all_orbits) {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            double occfactor = nbar_a * nbar_b * n_c * n_d -
                               n_a * n_b * nbar_c * nbar_d;
            if (std::abs(occfactor) < 1e-10)
              continue;

            int Jcd_min = std::abs(oc.j2 - od.j2) / 2;
            int Jcd_max = (oc.j2 + od.j2) / 2;
            int Jce_min = std::abs(oc.j2 - Z.modelspace->GetOrbit(e).j2) / 2;
            int Jce_max = (oc.j2 + Z.modelspace->GetOrbit(e).j2) / 2;

            int J0_min = std::max(Jab_min, Jcd_min);
            int J0_max = std::min(Jab_max, Jcd_max);
            int J1_min = std::max(Jab_min, Jce_min);
            int J1_max = std::min(Jab_max, Jce_max);

            for (int J0 = J0_min; J0 <= J0_max; ++J0) {
              for (int J1 = J1_min; J1 <= J1_max; ++J1) {
                if (not AngMom::Triangle(J0, J1, lambda))
                  continue;

                double x_cdab = Eta.TwoBody.GetTBME_J(J0, J1, c, d, a, b);
                double x_abce = Eta.TwoBody.GetTBME_J(J1, J0, a, b, c, e);
                if (std::abs(x_cdab) < 1e-12 || std::abs(x_abce) < 1e-12)
                  continue;

                int phase = Z.modelspace->phase(J0 + J1 + lambda);
                chi_de += 0.5 * occfactor * phase * hat_lambda_inv * x_cdab *
                          x_abce;
              }
            }
          }
        }
      }

      Chi_alpha(d, e) += chi_de / (od.j2 + 1.0);
      if (d != e)
        Chi_alpha(e, d) += chi_de / (od.j2 + 1.0);
    }
  }

  if (Commutator::verbose) {
    Z.profiler.timer["_231_F_tensor_chialpha"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // chi^beta_{de}^lambda from diag2_compact (Gamma*Omega tensor intermediate).
#pragma omp parallel for schedule(dynamic, 1)
  for (int indexd = 0; indexd < norbits; ++indexd) {
    auto d = allorb_vec[indexd];
    Orbit &od = Z.modelspace->GetOrbit(d);
    for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) {
      Orbit &oe = Z.modelspace->GetOrbit(e);

      double chi_de = 0.0;
      int pref_phase = Z.modelspace->phase(oe.j2 / 2 + lambda);

      for (auto &a : Z.modelspace->all_orbits) {
        Orbit &oa = Z.modelspace->GetOrbit(a);
        double n_a = oa.occ;
        double nbar_a = 1.0 - n_a;

        for (auto &b : Z.modelspace->all_orbits) {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_b = ob.occ;
          double nbar_b = 1.0 - n_b;

          int Jab_min = std::abs(oa.j2 - ob.j2) / 2;
          int Jab_max = (oa.j2 + ob.j2) / 2;

          for (auto &c : Z.modelspace->all_orbits) {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            double occfactor = nbar_a * nbar_b * n_c * n_d -
                               n_a * n_b * nbar_c * nbar_d;
            if (std::abs(occfactor) < 1e-10)
              continue;

            int J0_min = std::max(std::abs(oc.j2 - od.j2), Jab_min) / 1;
            int J0_max = std::min(oc.j2 + od.j2, 2 * Jab_max) / 2;
            int J1_min = std::max(Jab_min, std::abs(oc.j2 - oe.j2) / 2);
            int J1_max = std::min(Jab_max, (oc.j2 + oe.j2) / 2);

            for (int J0 = J0_min; J0 <= J0_max; ++J0) {
              for (int J1 = J1_min; J1 <= J1_max; ++J1) {
                if (not AngMom::Triangle(J0, J1, lambda))
                  continue;

                double omega_cdab = Eta.TwoBody.GetTBME_J(J0, J1, c, d, a, b);
                double gamma_abce = Gamma.TwoBody.GetTBME_J(J1, J1, a, b, c, e);
                if (std::abs(omega_cdab) < 1e-12 || std::abs(gamma_abce) < 1e-12)
                  continue;

                double sixj = Z.modelspace->GetSixJ(J0, J1, lambda,
                                                    oe.j2 / 2.0, od.j2 / 2.0,
                                                    oc.j2 / 2.0);
                if (std::abs(sixj) < 1e-12)
                  continue;

                int phase = Z.modelspace->phase(J0 + oc.j2 / 2);
                chi_de += 0.5 * pref_phase * occfactor * phase *
                          std::sqrt((2.0 * J0 + 1.0) * (2.0 * J1 + 1.0)) *
                          sixj * omega_cdab * gamma_abce;
              }
            }
          }
        }
      }
      Chi_beta(d, e) = chi_de;
    }
  }

  if (Commutator::verbose) {
    Z.profiler.timer["_231_F_tensor_chibeta"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  // f^(I): contract chi^alpha with scalar Gamma as in diag2_compact.
#pragma omp parallel for schedule(dynamic, 1)
  for (int indexp = 0; indexp < norbits; ++indexp) {
    auto p = allorb_vec[indexp];
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) {
      if (q > p)
        continue;
      Orbit &oq = Z.modelspace->GetOrbit(q);

      double zI = 0.0;
      double zII = 0.0;
      for (auto &a : Z.modelspace->all_orbits) {
        Orbit &oa = Z.modelspace->GetOrbit(a);
        for (auto &b : Z.modelspace->all_orbits) {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          if (oa.j2 != ob.j2)
            continue;

          int Jmin = std::max(std::abs(ob.j2 - op.j2), std::abs(oa.j2 - oq.j2)) / 2;
          int Jmax = std::min(ob.j2 + op.j2, oa.j2 + oq.j2) / 2;

          for (int J0 = Jmin; J0 <= Jmax; ++J0) {
            double g1 = Gamma.TwoBody.GetTBME_J(J0, J0, b, p, a, q);
            double g2 = Gamma.TwoBody.GetTBME_J(J0, J0, a, p, b, q);
            if (std::abs(g1) < 1e-12 && std::abs(g2) < 1e-12)
              continue;
            zI += (2 * J0 + 1) * Chi_alpha(a, b) * (g1 + g2);

            int J0_II_min = std::abs(ob.j2 - op.j2) / 2;
            int J0_II_max = (ob.j2 + op.j2) / 2;
            int J1_II_min = std::abs(oa.j2 - oq.j2) / 2;
            int J1_II_max = (oa.j2 + oq.j2) / 2;
            for (int J0_II = J0_II_min; J0_II <= J0_II_max; ++J0_II) {
              for (int J1 = J1_II_min; J1 <= J1_II_max; ++J1) {
                if (not AngMom::Triangle(J0_II, J1, lambda))
                  continue;

                double sixj_ab = Z.modelspace->GetSixJ(J1, lambda, J0_II,
                                                     ob.j2 / 2.0, op.j2 / 2.0,
                                                     oa.j2 / 2.0);
                double sixj_ba = Z.modelspace->GetSixJ(J1, lambda, J0_II,
                                                     oa.j2 / 2.0, op.j2 / 2.0,
                                                     ob.j2 / 2.0);
                if (std::abs(sixj_ab) < 1e-12 && std::abs(sixj_ba) < 1e-12)
                  continue;

                double omega_ab = Eta.TwoBody.GetTBME_J(J0_II, J1, b, p, a, q);
                double omega_ba = Eta.TwoBody.GetTBME_J(J0_II, J1, a, p, b, q);

                double pref = std::sqrt((2.0 * J0_II + 1.0) * (2.0 * J1 + 1.0)) *
                            hat_lambda_inv;
                zII += Z.modelspace->phase(J1 + oa.j2 / 2) * pref * sixj_ab *
                     Chi_beta(a, b) * omega_ab;
                zII -= Z.modelspace->phase(J1 + ob.j2 / 2) * pref * sixj_ba *
                     Chi_beta(b, a) * omega_ba;
              }
            }
          }
        }
      }

      double pref_I = 1.0 / (op.j2 + 1.0);
      double pref_II = Z.modelspace->phase(op.j2 / 2) / (op.j2 + 1.0);
      Z.OneBody(p, q) += pref_I * zI + pref_II * zII;
      if (p != q)
        Z.OneBody(q, p) += hZ * (pref_I * zI + pref_II * zII);
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////////////////////
/// tensor extension of 223_231 two-body-intermediate block
/// Currently implements the f^(III_b) topology from diag2 in direct form.
////////////////////////////////////////////////////////////////////////////
void comm223_231_chi2b_tensor(const Operator &Eta, const Operator &Gamma,
                              Operator &Z) {

  double t_start = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();

  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv = 1.0 / std::sqrt(2.0 * lambda + 1.0);

  int norbits = Z.modelspace->all_orbits.size();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(),
                                  Z.modelspace->all_orbits.end());

  auto barred_tbme = [&](const Operator &Op, index_t a, index_t j, index_t k,
                         index_t b, int J0, int J1) -> double {
    const int lam = Op.GetJRank();
    Orbit &oa = Z.modelspace->GetOrbit(a);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    Orbit &ob = Z.modelspace->GetOrbit(b);

    int pref_phase = Z.modelspace->phase(J0 + (oa.j2 + ok.j2) / 2 + lam);
    double pref_hat = std::sqrt((2.0 * J0 + 1.0) * (2.0 * J1 + 1.0));

    double sum = 0.0;
    if (not AngMom::Triangle(oa.j2 / 2, ob.j2 / 2, J0))
      return 0.0;
    if (not AngMom::Triangle(oj.j2 / 2, ok.j2 / 2, J1))
      return 0.0;
    if (not AngMom::Triangle(J1, J0, lam))
      return 0.0;

    // Use cross-coupled ordering consistent with scalar Pandya coding:
    // barred indices (a,j,k,b) are mapped to underlying TBME legs (a,b;j,k).
    int J2min = std::abs(oa.j2 - ob.j2) / 2;
    int J2max = (oa.j2 + ob.j2) / 2;
    int J3min = std::abs(oj.j2 - ok.j2) / 2;
    int J3max = (oj.j2 + ok.j2) / 2;

    for (int J2 = J2min; J2 <= J2max; ++J2) {
      for (int J3 = J3min; J3 <= J3max; ++J3) {
        if (not AngMom::Triangle(J2, J3, lam))
          continue;

        int j0min = std::max({std::abs(oa.j2 - 2 * J3),
                              std::abs(ok.j2 - 2 * J0),
                              std::abs(oj.j2 - 2 * lam)}) /
                    2;
        int j0max = std::min({oa.j2 + 2 * J3,
                              ok.j2 + 2 * J0,
                              oj.j2 + 2 * lam}) /
                    2;

        for (int j0 = j0min; j0 <= j0max; ++j0) {
          double sixj1 = Z.modelspace->GetSixJ(lam, J3, J2, oa.j2 / 2.0,
                                               oj.j2 / 2.0, j0);
          double sixj2 = Z.modelspace->GetSixJ(oa.j2 / 2.0, ob.j2 / 2.0, J0,
                                               ok.j2 / 2.0, j0, J3);
          double sixj3 = Z.modelspace->GetSixJ(J1, J0, lam, j0,
                                               oj.j2 / 2.0, ok.j2 / 2.0);
          if (std::abs(sixj1) < 1e-12 || std::abs(sixj2) < 1e-12 ||
              std::abs(sixj3) < 1e-12)
            continue;

          double omega = Op.TwoBody.GetTBME_J(J2, J3, a, b, j, k);
          if (std::abs(omega) < 1e-12)
            continue;

          double hats = std::sqrt((2.0 * J2 + 1.0) * (2.0 * J3 + 1.0)) *
                        (2.0 * j0 + 1.0);
          sum += Z.modelspace->phase(J2) * hats * sixj1 * sixj2 * sixj3 *
                 omega;
        }
      }
    }

    return pref_phase * pref_hat * sum;
  };

  auto bar_chi_gamma = [&](index_t i, index_t l, index_t k, index_t j,
                           int J0) -> double {
    Orbit &ol = Z.modelspace->GetOrbit(l);
    Orbit &ok = Z.modelspace->GetOrbit(k);
    double n_l = ol.occ;
    double n_k = ok.occ;

    double sum = 0.0;
    for (auto &a : Z.modelspace->all_orbits) {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      double n_a = oa.occ;
      double nbar_a = 1.0 - n_a;
      for (auto &b : Z.modelspace->all_orbits) {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        double n_b = ob.occ;
        double nbar_b = 1.0 - n_b;

        int J2min = std::max(std::abs(oa.j2 - Z.modelspace->GetOrbit(j).j2),
                             std::abs(oa.j2 - ol.j2)) /
                    2;
        int J2max = std::min(oa.j2 + Z.modelspace->GetOrbit(j).j2,
                             oa.j2 + ol.j2) /
                    2;
        for (int J2 = J2min; J2 <= J2max; ++J2) {
          if (not AngMom::Triangle(J0, J2, lambda))
            continue;

          double occ = n_a * nbar_b * n_l * (1.0 - n_k) -
                       nbar_a * n_b * (1.0 - n_l) * n_k;
          if (std::abs(occ) < 1e-12)
            continue;

          double bo1 = barred_tbme(Eta, i, b, a, j, J0, J2);
          double bo2 = barred_tbme(Eta, a, l, k, b, J2, J0);
          if (std::abs(bo1) < 1e-12 || std::abs(bo2) < 1e-12)
            continue;

          sum += occ * Z.modelspace->phase(J2 + lambda) * hat_lambda_inv *
                 bo1 * bo2;
        }
      }
    }

    return Z.modelspace->phase(J0) * sum / (2.0 * J0 + 1.0);
  };

  auto chi_delta = [&](index_t i, index_t j, index_t k, index_t l,
                       int J0) -> double {
    Orbit &oi = Z.modelspace->GetOrbit(i);
    Orbit &oj = Z.modelspace->GetOrbit(j);
    double n_i = oi.occ;
    double n_j = oj.occ;

    double val = 0.0;
    for (auto &m : Z.modelspace->all_orbits) {
      Orbit &om = Z.modelspace->GetOrbit(m);
      double n_m = om.occ;
      double nbar_m = 1.0 - n_m;
      for (auto &n : Z.modelspace->all_orbits) {
        Orbit &on = Z.modelspace->GetOrbit(n);
        double n_n = on.occ;
        double nbar_n = 1.0 - n_n;

        double occ = (1.0 - n_i) * (1.0 - n_j) * n_m * n_n -
                     n_i * n_j * nbar_m * nbar_n;
        if (std::abs(occ) < 1e-10)
          continue;

        int J2min = std::max(std::abs(oi.j2 - oj.j2), std::abs(om.j2 - on.j2)) / 2;
        int J2max = std::min(oi.j2 + oj.j2, om.j2 + on.j2) / 2;
        for (int J2 = J2min; J2 <= J2max; ++J2) {
          if (not AngMom::Triangle(J0, J2, lambda))
            continue;

          double x1 = Eta.TwoBody.GetTBME_J(J0, J2, i, j, m, n);
          double x2 = Eta.TwoBody.GetTBME_J(J2, J0, m, n, k, l);
          if (std::abs(x1) < 1e-12 || std::abs(x2) < 1e-12)
            continue;

          int phase = Z.modelspace->phase(J2 + lambda);
          val += occ * phase * hat_lambda_inv * x1 * x2;
        }
      }
    }

    return 0.25 * Z.modelspace->phase(J0) * val / (2.0 * J0 + 1.0);
  };

#pragma omp parallel for schedule(dynamic, 1)
  for (int ip = 0; ip < norbits; ++ip) {
    auto p = allorb_vec[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);

    for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) {
      if (q > p)
        continue;

      double zIIIb = 0.0;
      double zIIIa = 0.0;
      for (auto &c : Z.modelspace->all_orbits) {
        Orbit &oc = Z.modelspace->GetOrbit(c);
        for (auto &a : Z.modelspace->all_orbits) {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (auto &b : Z.modelspace->all_orbits) {
            Orbit &ob = Z.modelspace->GetOrbit(b);

            int J0min = std::max({std::abs(oc.j2 - op.j2), std::abs(oa.j2 - ob.j2),
                                  std::abs(oc.j2 - Z.modelspace->GetOrbit(q).j2)}) /
                        2;
            int J0max = std::min({oc.j2 + op.j2, oa.j2 + ob.j2,
                                  oc.j2 + Z.modelspace->GetOrbit(q).j2}) /
                        2;
            for (int J0 = J0min; J0 <= J0max; ++J0) {
              double bg_i_cab = bar_chi_gamma(p, c, a, b, J0);
              double bg_c_jab = bar_chi_gamma(c, q, a, b, J0);
              double bG_abjc = barred_tbme(Gamma, a, b, q, c, J0, J0);
              double bG_abci = barred_tbme(Gamma, a, b, c, p, J0, J0);
              zIIIa += (2.0 * J0 + 1.0) *
                       (bg_i_cab * bG_abjc - bg_c_jab * bG_abci);

              double chi_ciab = chi_delta(c, p, a, b, J0);
              double chi_abcq = chi_delta(a, b, c, q, J0);
              double g_abcq = Gamma.TwoBody.GetTBME_J(J0, J0, a, b, c, q);
              double g_cpab = Gamma.TwoBody.GetTBME_J(J0, J0, c, p, a, b);
              if (std::abs(chi_ciab) < 1e-12 && std::abs(chi_abcq) < 1e-12)
                continue;
              zIIIb += (2.0 * J0 + 1.0) * (chi_ciab * g_abcq - g_cpab * chi_abcq);
            }
          }
        }
      }

      double pref = 1.0 / (op.j2 + 1.0);
      Z.OneBody(p, q) += pref * (zIIIa + zIIIb);
      if (p != q)
        Z.OneBody(q, p) += hZ * pref * (zIIIa + zIIIb);
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z) {

  // Do some extra work to check if the operators being passed in are reduced,
  // because the comm223_232 routines assume not-reduced matrix elements
  const Operator *Etanred =
      &Eta; // Pointer to the non-reduced version of the operator
  const Operator *Gammanred = &Gamma;
  Operator Etatmp,
      Gammatmp;        // Declared, but not yet allocated, for reasons of scope.
  if (Eta.IsReduced()) // comm223_231 doesn't expect reduced operators. Need to
                       // make it not reduced.
  {
    Etatmp = Eta; // Actual copy, not a reference
    Etatmp.MakeNotReduced();
    Etanred = &Etatmp; // Now Etanred points to the not-reduced copy Etatmp
  }
  if (Gamma.IsReduced()) {
    Gammatmp = Gamma;
    Gammatmp.MakeNotReduced();
    Gammanred =
        &Gammatmp; // Now Gammanred points to the not-reduced copy Gammatmp
  }

  bool z_was_reduced = Z.IsReduced();
  if (z_was_reduced)
    Z.MakeNotReduced();

  if (use_1b_intermediates) {
    comm223_232_chi1b_tensor(*Etanred, *Gammanred,
                             Z); // topology with 1-body intermediate (fast)
  }
  if (use_2b_intermediates) {
    comm223_232_chi2b(*Etanred, *Gammanred,
                      Z); // topology with 2-body intermediate (slow)
  }

  if (z_was_reduced)
    Z.MakeReduced(); // put it back the way we found it...

  return;
}

////////////////////////////////////////////////////////////////////////////
/// factorized 223_232 double commutator with 1b intermediate
////////////////////////////////////////////////////////////////////////////
void comm223_232_chi1b_tensor(const Operator &Eta, const Operator &Gamma,
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
  // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB
  // channels int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC(); //
  // number of CC channels
  auto &Z2 = Z.TwoBody;

  bool Z_is_scalar = true;
  if (Z.TwoBody.rank_T != 0) {
    Z_is_scalar = false;
  }
  // determine symmetry
  int hEta = Eta.IsHermitian() ? 1 : -1;
  int hGamma = Gamma.IsHermitian() ? 1 : -1;
  int hZ = hGamma;
  bool tensor_case = (Eta.GetJRank() != 0) && (Gamma.GetJRank() == 0) && Z_is_scalar;
  int lambda = Eta.GetJRank();
  double hat_lambda_inv =
      tensor_case ? 1.0 / std::sqrt(2.0 * lambda + 1.0) : 0.0;
  // ####################################################################################
  //                      Factorization of Ia, Ib, IVa and IVb
  // ####################################################################################

  arma::mat CHI_I = Gamma.OneBody * 0;
  arma::mat CHI_II = Gamma.OneBody * 0;

  // The intermidate one body operator
  //  CHI_I :                            //  CHI_II :
  //          eta | p                    //          eta | p
  //         _____|                      //         _____|
  //       /\     |                      //       /\     |
  //   a  (  ) b  | c                    //   a  (  ) b  | c
  //       \/_____|                      //       \/~~~~~|
  //          eta |                      //        gamma |
  //              | q                    //              | q
  //-------------------------------------------------------------------------------------
  // CHI_I_pq  = 1/2 \sum_abcJ2 \hat(J_2) ( \bar{n}_a \bar{n}_c n_b - \bar{n}_c
  // n_a n_c )
  //             eta^J2_bpac eta^J2_acbq
  //
  // CHI_II_pq = 1/2 \sum_abcJ2 \hat(J_2) ( \bar{n}_b \bar{n}_c n_a - \bar{n}_a
  // n_b n_c )
  //             eta^J2_bcaq gamma^J2_apbc
  //-------------------------------------------------------------------------------------

  // Rolling the two loops into one helps with load balancing
  std::vector<index_t> p_list, q_list;
  for (auto p : Z.modelspace->all_orbits) {
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
      p_list.push_back(p);
      q_list.push_back(q);
    }
  }
  size_t ipq_max = p_list.size();
/// Build one-body intermediates (scalar original path + tensor chi^epsilon/chi^zeta path)
#pragma omp parallel for schedule(dynamic)
    for (size_t ipq = 0; ipq < ipq_max; ipq++) {
      index_t p = p_list[ipq];
      index_t q = q_list[ipq];
      //      for (size_t p = 0; p < norbits; p++)
      //      {
      Orbit &op = Z.modelspace->GetOrbit(p);
      //        for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
      //        {
      Orbit &oq = Z.modelspace->GetOrbit(q);

      double chi_pq = 0;
      double chiY_pq = 0;

      if (!tensor_case) {
        for (auto a : Z.modelspace->all_orbits) {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;
          if (nbar_a < 1e-6)
            continue;

          for (auto i : Z.modelspace->holes) {
            Orbit &oi = Z.modelspace->GetOrbit(i);
            double n_i = oi.occ;

            for (auto j : Z.modelspace->holes) {
              Orbit &oj = Z.modelspace->GetOrbit(j);
              double n_j = oj.occ;

              double occfactor = nbar_a * n_i * n_j;
              if (occfactor < 1.e-7) {
                continue;
              }

              int J2min =
                  std::max(std::abs(oa.j2 - oq.j2), std::abs(oi.j2 - oj.j2)) / 2;
              int J2max = std::min(oa.j2 + oq.j2, oi.j2 + oj.j2) / 2;

              for (int J2 = J2min; J2 <= J2max; J2++) {
                double xijaq = Eta.TwoBody.GetTBME_J(J2, J2, i, j, a, q);
                double xapij, yapij;
                if (Z_is_scalar) {
                  Eta.TwoBody.GetTBME_J_twoOps(Gamma.TwoBody, J2, J2, a, p, i, j,
                                               xapij, yapij);
                } else {
                  xapij = Eta.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                  yapij = Gamma.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                }
                chi_pq +=
                    0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xapij * xijaq;
                chiY_pq +=
                    0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yapij * xijaq;
              }
            } // for j

            for (auto b : Z.modelspace->all_orbits) {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double occfactor = nbar_a * nbar_b * n_i;

              if (std::abs(occfactor) < 1e-7)
                continue;

              int J2min =
                  std::max({std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2),
                            std::abs(oi.j2 - op.j2)}) /
                  2;
              int J2max =
                  std::min({oa.j2 + ob.j2, oi.j2 + oq.j2, oi.j2 + op.j2}) / 2;

              for (int J2 = J2min; J2 <= J2max; J2++) {
                double xabiq = Eta.TwoBody.GetTBME_J(J2, J2, a, b, i, q);
                double xipab, yipab;

                if (Z_is_scalar) {
                  Eta.TwoBody.GetTBME_J_twoOps(Gamma.TwoBody, J2, J2, i, p, a, b,
                                               xipab, yipab);
                } else {
                  xipab = Eta.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                  yipab = Gamma.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                }

                chi_pq +=
                    0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xipab * xabiq;
                chiY_pq +=
                    0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yipab * xabiq;
              }
            } // for b

          } // for i
        } // for a
        CHI_I(p, q) = chi_pq;
        CHI_II(p, q) = chiY_pq;
      } else {
        int zeta_pref_phase = Z.modelspace->phase(oq.j2 / 2 + lambda);
        for (auto &a : Z.modelspace->all_orbits) {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits) {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;

            int Jab_min = std::abs(oa.j2 - ob.j2) / 2;
            int Jab_max = (oa.j2 + ob.j2) / 2;

            for (auto &c : Z.modelspace->all_orbits) {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              double occ_eps_zeta = nbar_a * nbar_b * n_c + n_a * n_b * nbar_c;
              if (std::abs(occ_eps_zeta) < 1e-12)
                continue;

              int Jcp_min = std::abs(oc.j2 - op.j2) / 2;
              int Jcp_max = (oc.j2 + op.j2) / 2;
              int Jcq_min = std::abs(oc.j2 - oq.j2) / 2;
              int Jcq_max = (oc.j2 + oq.j2) / 2;
              int Jad_min = std::abs(oa.j2 - op.j2) / 2;
              int Jad_max = (oa.j2 + op.j2) / 2;
              int Jaq_min = std::abs(oa.j2 - oq.j2) / 2;
              int Jaq_max = (oa.j2 + oq.j2) / 2;
              int Jbc_min = std::abs(ob.j2 - oc.j2) / 2;
              int Jbc_max = (ob.j2 + oc.j2) / 2;

              int J0_eps_min = std::max(Jab_min, Jcp_min);
              int J0_eps_max = std::min(Jab_max, Jcp_max);
              int J1_eps_min = std::max(Jab_min, Jcq_min);
              int J1_eps_max = std::min(Jab_max, Jcq_max);
              for (int J0 = J0_eps_min; J0 <= J0_eps_max; ++J0) {
                for (int J1 = J1_eps_min; J1 <= J1_eps_max; ++J1) {
                  if (not AngMom::Triangle(J0, J1, lambda))
                    continue;
                  double o1 = Eta.TwoBody.GetTBME_J(J0, J1, c, p, a, b);
                  double o2 = Eta.TwoBody.GetTBME_J(J1, J0, a, b, c, q);
                  if (std::abs(o1) < 1e-12 || std::abs(o2) < 1e-12)
                    continue;
                  chi_pq += 0.5 * occ_eps_zeta *
                            Z.modelspace->phase(J0 + J1 + lambda) *
                            hat_lambda_inv * o1 * o2;
                }
              }

              int J0_zeta_min = std::max(Jbc_min, Jad_min);
              int J0_zeta_max = std::min(Jbc_max, Jad_max);
              int J1_zeta_min = Jaq_min;
              int J1_zeta_max = Jaq_max;
              for (int J0 = J0_zeta_min; J0 <= J0_zeta_max; ++J0) {
                for (int J1 = J1_zeta_min; J1 <= J1_zeta_max; ++J1) {
                  if (not AngMom::Triangle(J0, J1, lambda))
                    continue;

                  double sixj = Z.modelspace->GetSixJ(
                      lambda, J1, J0, oa.j2 / 2.0, op.j2 / 2.0, oq.j2 / 2.0);
                  if (std::abs(sixj) < 1e-12)
                    continue;

                  double gamma_apbc = Gamma.TwoBody.GetTBME_J(J0, J0, a, p, b, c);
                  double omega_bcaq = Eta.TwoBody.GetTBME_J(J0, J1, b, c, a, q);
                  if (std::abs(gamma_apbc) < 1e-12 || std::abs(omega_bcaq) < 1e-12)
                    continue;

                  chiY_pq += 0.5 * zeta_pref_phase * occ_eps_zeta *
                             Z.modelspace->phase(J0 + oa.j2 / 2) *
                             std::sqrt((2.0 * J0 + 1.0) * (2.0 * J1 + 1.0)) *
                             sixj * gamma_apbc * omega_bcaq;
                }
              }
            }
          }
        }

        CHI_I(p, q) = chi_pq / (op.j2 + 1.0);
        CHI_II(p, q) = chiY_pq;
      }
      //        } // for q
    } // for p

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

/// Now use the intermediate to form the double commutator
#pragma omp parallel for schedule(dynamic, 1)
  for (int ich = 0; ich < nch; ich++) {
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
      int phasepq = bra.Phase(J);

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
        double gamma2b, eta2b;

        // SRS Modified so that we don't have to re-look up the channel and
        // bra/ket indices. This requires a bit more work to get the
        // normalization and phases right but at emax=6 is speeds this up by a
        // factor of 2.
        for (auto b : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2})) {
          auto ibra_bq = tbc_bra.GetLocalIndex(std::min(b, q), std::max(b, q));
          if (ibra_bq < 0 or ibra_bq > nbras)
            continue;
          double norm = (b == q ? PhysConst::SQRT2 : 1) *
                        (p == q ? 1 / PhysConst::SQRT2 : 1);
          if (b > q)
            norm *= bra.Phase(tbc_bra.J);
          zpqrs += norm * CHI_I(p, b) *
                   Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq, iket);
          if (Z_is_scalar and not tensor_case)
            zpqrs += norm * hZ * CHI_II(b, p) *
                     Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq, iket); // add to ket side? not sure
          // zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
          // zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r,
          // s);
        }
        for (auto b : Eta.OneBodyChannels.at({oq.l, oq.j2, oq.tz2})) {
          auto ibra_pb = tbc_bra.GetLocalIndex(std::min(p, b), std::max(p, b));
          if (ibra_pb < 0 or ibra_pb > nbras)
            continue;
          double norm = (b == p ? PhysConst::SQRT2 : 1) *
                        (p == q ? 1 / PhysConst::SQRT2 : 1);
          if (p > b)
            norm *= bra.Phase(tbc_bra.J);
          zpqrs += norm * CHI_I(q, b) *
                   Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb, iket);
          if (Z_is_scalar and not tensor_case)
            zpqrs += norm * hZ * CHI_II(b, q) *
                     Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb, iket);//  why
          // zpqrs += CHI_I(q, b) *     Gamma.TwoBody.GetTBME_J(J, J, p, b, r,
          // s); zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b,
          // r, s);
        }
        for (auto b : Eta.OneBodyChannels.at({oR.l, oR.j2, oR.tz2})) {
          auto iket_bs = tbc_ket.GetLocalIndex(std::min(b, s), std::max(b, s));
          if (iket_bs < 0 or iket_bs > nkets)
            continue;
          double norm = (b == s ? PhysConst::SQRT2 : 1) *
                        (r == s ? 1 / PhysConst::SQRT2 : 1);
          if (b > s)
            norm *= ket.Phase(tbc_ket.J);
          zpqrs += norm *
                   Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) *
                   CHI_I(b, r);
          if (Z_is_scalar and not tensor_case)
            zpqrs -= norm *
                     Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) *
                     CHI_II(b, r);
          // zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_I(b, r);
          // zpqrs -=   Eta.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_II(b, r);
        }
        for (auto b : Eta.OneBodyChannels.at({os.l, os.j2, os.tz2})) {
          auto iket_rb = tbc_ket.GetLocalIndex(std::min(r, b), std::max(r, b));
          if (iket_rb < 0 or iket_rb > nkets)
            continue;
          double norm = (b == r ? PhysConst::SQRT2 : 1) *
                        (r == s ? 1 / PhysConst::SQRT2 : 1);
          if (r > b)
            norm *= ket.Phase(tbc_ket.J);
          zpqrs += norm *
                   Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) *
                   CHI_I(b, s);
          if (Z_is_scalar and not tensor_case)
            zpqrs -= norm *
                     Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) *
                     CHI_II(b, s);
          // zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_I(b, s);
          // zpqrs -=   Eta.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_II(b, s);
        }

        //            if (p == q)
        //              zpqrs /= PhysConst::SQRT2;
        //            if (r == s)
        //              zpqrs /= PhysConst::SQRT2;
        if (tensor_case) {
          int J0 = J;
          double gammaII = 0.0;
          double hatJ0_inv = 1.0 / std::sqrt(2.0 * J0 + 1.0);
          int phase_i = Z.modelspace->phase(op.j2 / 2);
          int phase_l = Z.modelspace->phase(os.j2 / 2);

          for (auto a : Z.modelspace->all_orbits) {
            Orbit &oa = Z.modelspace->GetOrbit(a);

            int J2min_ia_kl = std::max(std::abs(op.j2 - oa.j2), std::abs(oR.j2 - os.j2)) / 2;
            int J2max_ia_kl = std::min(op.j2 + oa.j2, oR.j2 + os.j2) / 2;
            int J2min_ij_al = std::max(std::abs(oa.j2 - os.j2), std::abs(op.j2 - oq.j2)) / 2;
            int J2max_ij_al = std::min(oa.j2 + os.j2, op.j2 + oq.j2) / 2;
            int J2min = std::min(J2min_ia_kl, J2min_ij_al);
            int J2max = std::max(J2max_ia_kl, J2max_ij_al);

            for (int J2 = J2min; J2 <= J2max; ++J2) {
              if (not AngMom::Triangle(J0, J2, lambda))
                continue;

              double pref = Z.modelspace->phase(J2 + oa.j2 / 2) *
                            std::sqrt(2.0 * J2 + 1.0) * hat_lambda_inv;

              if (J2 >= J2min_ia_kl && J2 <= J2max_ia_kl) {
                double sixj_1 = Z.modelspace->GetSixJ(J0, J2, lambda,
                                                      oa.j2 / 2.0, oq.j2 / 2.0,
                                                      op.j2 / 2.0);
                if (std::abs(sixj_1) > 1e-12) {
                  double zeta_ja = CHI_II(q, a);
                  double omega_iakl = Eta.TwoBody.GetTBME_J(J2, J0, p, a, r, s);
                  gammaII += phase_i * hatJ0_inv * pref * sixj_1 * zeta_ja * omega_iakl;
                }
              }

              if (J2 >= J2min_ij_al && J2 <= J2max_ij_al) {
                double sixj_2 = Z.modelspace->GetSixJ(J0, J2, lambda,
                                                      oa.j2 / 2.0, oR.j2 / 2.0,
                                                      os.j2 / 2.0);
                if (std::abs(sixj_2) > 1e-12) {
                  double zeta_ak = CHI_II(a, r);
                  double omega_ijal = Eta.TwoBody.GetTBME_J(J0, J2, p, q, a, s);
                  gammaII -= phase_l * hatJ0_inv * pref * sixj_2 * zeta_ak * omega_ijal;
                }
              }
            }
          }
          zpqrs += gammaII;
        }
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
      } // for iket
    } // for ibra

  } // for itmat

  if (Commutator::verbose) {
    Z.profiler
        .timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] +=
        omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
  }

  CHI_I.clear();
  CHI_II.clear();

  // Timer
  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return;
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
  TwoBodyME Chi_III_Op = Eta.TwoBody;
  Chi_III_Op.Erase();
  // Inverse Pandya transformation
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
          if (not tensor_eta_case) {
            commij -= (2 * Jprime + 1) * sixj * me1;
          } else {
            // Tensor reverse Pandya recoupling (adcb convention in barred element).
            // Here (i,j,k,l) -> barred access (i,l ; k,j).
            double rec = 0.0;
            int j0min = std::max({std::abs(ji - jj), std::abs(ji - 2 * Jprime),
                                  std::abs(jk - 2 * J0)}) /
                        2;
            int j0max = std::min({ji + jj, ji + 2 * Jprime, jk + 2 * J0}) / 2;
            for (int j0 = j0min; j0 <= j0max; ++j0) {
              double sixj1_t = Z.modelspace->GetSixJ(lambda, Jprime, Jprime,
                                                     ji * 0.5, jj * 0.5, j0);
              double sixj2_t = Z.modelspace->GetSixJ(ji * 0.5, jl * 0.5, J0,
                                                     jk * 0.5, j0, Jprime);
              double sixj3_t = Z.modelspace->GetSixJ(J0, J0, lambda, j0,
                                                     jj * 0.5, jk * 0.5);
              if (std::abs(sixj1_t) < 1e-12 || std::abs(sixj2_t) < 1e-12 ||
                  std::abs(sixj3_t) < 1e-12)
                continue;
              rec += Z.modelspace->phase(Jprime) * (2.0 * Jprime + 1.0) *
                     (2.0 * j0 + 1.0) * sixj1_t * sixj2_t * sixj3_t;
            }
            double pref = Z.modelspace->phase(J0 + ji / 2 + jk / 2 + lambda) *
                          (2.0 * J0 + 1.0);
            commij += pref * rec * me1;
          }
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
          if (not tensor_eta_case) {
            commji -= (2 * Jprime + 1) * sixj * me1;
          } else {
            // Tensor reverse Pandya recoupling for exchanged ij leg.
            double rec = 0.0;
            int j0min = std::max({std::abs(jj - ji), std::abs(jj - 2 * Jprime),
                                  std::abs(jk - 2 * J0)}) /
                        2;
            int j0max = std::min({jj + ji, jj + 2 * Jprime, jk + 2 * J0}) / 2;
            for (int j0 = j0min; j0 <= j0max; ++j0) {
              double sixj1_t = Z.modelspace->GetSixJ(lambda, Jprime, Jprime,
                                                     jj * 0.5, ji * 0.5, j0);
              double sixj2_t = Z.modelspace->GetSixJ(jj * 0.5, jl * 0.5, J0,
                                                     jk * 0.5, j0, Jprime);
              double sixj3_t = Z.modelspace->GetSixJ(J0, J0, lambda, j0,
                                                     ji * 0.5, jk * 0.5);
              if (std::abs(sixj1_t) < 1e-12 || std::abs(sixj2_t) < 1e-12 ||
                  std::abs(sixj3_t) < 1e-12)
                continue;
              rec += Z.modelspace->phase(Jprime) * (2.0 * Jprime + 1.0) *
                     (2.0 * j0 + 1.0) * sixj1_t * sixj2_t * sixj3_t;
            }
            double pref = Z.modelspace->phase(J0 + jj / 2 + jk / 2 + lambda) *
                          (2.0 * J0 + 1.0);
            commji += pref * rec * me1;
          }
        }

        double zijkl =
            (commij - Z.modelspace->phase((ji + jj) / 2 - J0) * commji);
        //          Chi_III[ch](ibra, iket) += zijkl;
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
