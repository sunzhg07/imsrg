#include "evc.hh"
#include "AngMom.hh"
#include "BCH.hh"
#include "PhysicalConstants.hh"
#include "ReferenceImplementations.hh"

#include <algorithm>
#include <cmath>
#include <deque>
#include <omp.h>
#include <utility>
#include <vector>

namespace
{
  constexpr double kHalf = 1.0 / 2.0;
  constexpr double kSixth = 1.0 / 6.0;
  constexpr double kTwentyFourth = 1.0 / 24.0;

  bool SameOneBodyChannel(const Orbit &oa, const Orbit &ob)
  {
    return oa.j2 == ob.j2 and oa.tz2 == ob.tz2 and ((oa.l + ob.l) % 2 == 0);
  }

  int TwoBodyJmin(const Orbit &oa, const Orbit &ob, const Orbit &oc, const Orbit &od)
  {
    return std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
  }

  int TwoBodyJmax(const Orbit &oa, const Orbit &ob, const Orbit &oc, const Orbit &od)
  {
    return std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
  }

  // AMC unreduced χ_hp_ck = ĵ_c^{-2} Σ_{bj J} Ĵ^{2} t^J_{bcjk} z_bj
  // ĵ = √(2j+1) so ĵ_c^{-2} = 1/(2j_c+1)
  arma::mat BuildChiHP(const Operator &t, const Operator &z,
                       const std::vector<index_t> &holes,
                       const std::vector<index_t> &particles)
  {
    const size_t norb = t.modelspace->GetNumberOrbits();
    arma::mat chi(norb, norb, arma::fill::zeros);

    std::vector<std::pair<index_t, index_t>> pairs;
    for (auto k : holes)
    {
      const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
      for (auto c : particles)
      {
        const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
        if (SameOneBodyChannel(oc, ok))
        {
          pairs.emplace_back(c, k);
        }
      }
    }

    const int npairs = static_cast<int>(pairs.size());
#pragma omp parallel for schedule(dynamic)
    for (int ip = 0; ip < npairs; ++ip)
    {
      const index_t c = pairs[ip].first;
      const index_t k = pairs[ip].second;
      const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
      const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
      double sum = 0.0;
      for (auto j : holes)
      {
        const Orbit &oj = t.modelspace->GetOrbit(static_cast<int>(j));
        for (auto b : particles)
        {
          const Orbit &ob = t.modelspace->GetOrbit(static_cast<int>(b));
          if (not SameOneBodyChannel(ob, oj))
          {
            continue;
          }
          const double z1 = z.OneBody(b, j);
          if (std::abs(z1) < 1e-16)
          {
            continue;
          }
          const int Jmin = TwoBodyJmin(ob, oc, oj, ok);
          const int Jmax = TwoBodyJmax(ob, oc, oj, ok);
          for (int J = Jmin; J <= Jmax; ++J)
          {
            sum += (2 * J + 1) * t.TwoBody.GetTBME_J(J, static_cast<int>(b), static_cast<int>(c), static_cast<int>(j), static_cast<int>(k)) * z1;
          }
        }
      }
      // AMC unreduced: ĵ_c^{-2} = 1/(2j_c+1)
      chi(c, k) = sum / (oc.j2 + 1.0);
    }
    return chi;
  }

  // AMC unreduced χ_hh_ji (term 5)
  arma::mat BuildChiHH(const Operator &t, const Operator &z,
                       const std::vector<index_t> &holes,
                       const std::vector<index_t> &particles)
  {
    const size_t norb = t.modelspace->GetNumberOrbits();
    arma::mat chi(norb, norb, arma::fill::zeros);

    std::vector<std::pair<index_t, index_t>> pairs;
    for (auto i : holes)
    {
      const Orbit &oi = t.modelspace->GetOrbit(static_cast<int>(i));
      for (auto j : holes)
      {
        const Orbit &oj = t.modelspace->GetOrbit(static_cast<int>(j));
        if (SameOneBodyChannel(oj, oi))
        {
          pairs.emplace_back(j, i);
        }
      }
    }

    const int npairs = static_cast<int>(pairs.size());
#pragma omp parallel for schedule(dynamic)
    for (int ip = 0; ip < npairs; ++ip)
    {
      const index_t j = pairs[ip].first;
      const index_t i = pairs[ip].second;
      const Orbit &oj = t.modelspace->GetOrbit(static_cast<int>(j));
      const Orbit &oi = t.modelspace->GetOrbit(static_cast<int>(i));
      double sum = 0.0;
      for (auto k : holes)
      {
        const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
        for (auto b : particles)
        {
          const Orbit &ob = t.modelspace->GetOrbit(static_cast<int>(b));
          for (auto c : particles)
          {
            const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
            const int Jmin = std::max(TwoBodyJmin(ob, oc, oj, ok), TwoBodyJmin(ob, oc, ok, oi));
            const int Jmax = std::min(TwoBodyJmax(ob, oc, oj, ok), TwoBodyJmax(ob, oc, ok, oi));
            for (int J = Jmin; J <= Jmax; ++J)
            {
              const int phase = AngMom::phase((oi.j2 + ok.j2) / 2 + J);
              sum += phase * (2 * J + 1)
                     * t.TwoBody.GetTBME_J(J, static_cast<int>(b), static_cast<int>(c), static_cast<int>(j), static_cast<int>(k))
                     * z.TwoBody.GetTBME_J(J, static_cast<int>(b), static_cast<int>(c), static_cast<int>(k), static_cast<int>(i));
            }
          }
        }
      }
      // AMC unreduced: ĵ_i^{-2} = 1/(2j_i+1)
      chi(j, i) = 0.5 * sum / (oi.j2 + 1.0);
    }
    return chi;
  }

  // AMC unreduced χ_pp_ba (term 4)
  arma::mat BuildChiPP(const Operator &t, const Operator &z,
                       const std::vector<index_t> &holes,
                       const std::vector<index_t> &particles)
  {
    const size_t norb = t.modelspace->GetNumberOrbits();
    arma::mat chi(norb, norb, arma::fill::zeros);

    std::vector<std::pair<index_t, index_t>> pairs;
    for (auto a : particles)
    {
      const Orbit &oa = t.modelspace->GetOrbit(static_cast<int>(a));
      for (auto b : particles)
      {
        const Orbit &ob = t.modelspace->GetOrbit(static_cast<int>(b));
        if (SameOneBodyChannel(ob, oa))
        {
          pairs.emplace_back(b, a);
        }
      }
    }

    const int npairs = static_cast<int>(pairs.size());
#pragma omp parallel for schedule(dynamic)
    for (int ip = 0; ip < npairs; ++ip)
    {
      const index_t b = pairs[ip].first;
      const index_t a = pairs[ip].second;
      const Orbit &ob = t.modelspace->GetOrbit(static_cast<int>(b));
      const Orbit &oa = t.modelspace->GetOrbit(static_cast<int>(a));
      double sum = 0.0;
      for (auto c : particles)
      {
        const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
        for (auto j : holes)
        {
          const Orbit &oj = t.modelspace->GetOrbit(static_cast<int>(j));
          for (auto k : holes)
          {
            const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
            const int Jmin = std::max(TwoBodyJmin(ob, oc, oj, ok), TwoBodyJmin(oc, oa, oj, ok));
            const int Jmax = std::min(TwoBodyJmax(ob, oc, oj, ok), TwoBodyJmax(oc, oa, oj, ok));
            for (int J = Jmin; J <= Jmax; ++J)
            {
              const int phase = AngMom::phase((oa.j2 + oc.j2) / 2 + J);
              sum += phase * (2 * J + 1)
                     * t.TwoBody.GetTBME_J(J, static_cast<int>(b), static_cast<int>(c), static_cast<int>(j), static_cast<int>(k))
                     * z.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(a), static_cast<int>(j), static_cast<int>(k));
            }
          }
        }
      }
      // AMC unreduced: ĵ_a^{-2} = 1/(2j_a+1)
      chi(b, a) = 0.5 * sum / (oa.j2 + 1.0);
    }
    return chi;
  }

  Operator BuildWFromZ1(const Operator &z,
                        const std::vector<index_t> &holes,
                        const std::vector<index_t> &particles)
  {
    Operator w = z;
    w.Erase();
    w.SetNonHermitian();
    ModelSpace &ms = *z.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const Ket &bra = tbc.GetKet(static_cast<int>(pp(ib)));
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const int phase = AngMom::phase(J + (oa.j2 + ob.j2) / 2);
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const Ket &ket = tbc.GetKet(static_cast<int>(hh(ik)));
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          double me = 0.0;
          if (SameOneBodyChannel(oa, oi) and SameOneBodyChannel(ob, oj))
          {
            me += z.OneBody(a, i) * z.OneBody(b, j);
          }
          if (SameOneBodyChannel(oa, oj) and SameOneBodyChannel(ob, oi))
          {
            me -= phase * z.OneBody(a, j) * z.OneBody(b, i);
          }
          if (std::abs(me) < 1e-16)
          {
            continue;
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          w.TwoBody.SetTBME_J(J, a, b, i, j, me / nrm);
        }
      }
    }
    (void)holes;
    (void)particles;
    return w;
  }

  arma::mat BuildI6(const Operator &t, const Operator &z, const arma::mat &chi_pp,
                    const std::vector<index_t> &holes,
                    const std::vector<index_t> &particles)
  {
    const size_t norb = t.modelspace->GetNumberOrbits();
    arma::mat i6(norb, norb, arma::fill::zeros);
    for (auto a : particles)
    {
      const Orbit &oa = t.modelspace->GetOrbit(static_cast<int>(a));
      for (auto c : particles)
      {
        const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
        if (not SameOneBodyChannel(oa, oc))
        {
          continue;
        }
        double s = 0.0;
        for (auto k : holes)
        {
          const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
          if (SameOneBodyChannel(ok, oa))
          {
            s += t.OneBody(c, k) * z.OneBody(a, k);
          }
        }
        i6(a, c) = s;
      }
    }
    (void)chi_pp;
    return i6;
  }

  arma::mat BuildI7(const Operator &t, const Operator &z, const arma::mat &chi_hh,
                    const std::vector<index_t> &holes,
                    const std::vector<index_t> &particles)
  {
    const size_t norb = t.modelspace->GetNumberOrbits();
    arma::mat i7(norb, norb, arma::fill::zeros);
    for (auto i : holes)
    {
      const Orbit &oi = t.modelspace->GetOrbit(static_cast<int>(i));
      for (auto k : holes)
      {
        const Orbit &ok = t.modelspace->GetOrbit(static_cast<int>(k));
        if (not SameOneBodyChannel(ok, oi))
        {
          continue;
        }
        double s = 0.0;
        for (auto c : particles)
        {
          const Orbit &oc = t.modelspace->GetOrbit(static_cast<int>(c));
          if (SameOneBodyChannel(oc, oi))
          {
            s += t.OneBody(c, k) * z.OneBody(c, i);
          }
        }
        i7(k, i) = s;
      }
    }
    (void)chi_hh;
    return i7;
  }

  void FoldI6I7(const arma::mat &i6, const arma::mat &i7, const Operator &z, Operator &dz)
  {
    ModelSpace &ms = *z.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const int ph_ab = AngMom::phase(J + (oa.j2 + ob.j2) / 2);
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          const int ph_ij = AngMom::phase(J + (oi.j2 + oj.j2) / 2);
          double r = 0.0;
          for (auto c : ms.particles)
          {
            const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
            if (SameOneBodyChannel(oc, oa))
            {
              r -= i6(a, c) * z.TwoBody.GetTBME_J(J, static_cast<int>(c), b, i, j);
            }
            if (SameOneBodyChannel(oc, ob))
            {
              r += ph_ab * i6(b, c) * z.TwoBody.GetTBME_J(J, static_cast<int>(c), a, i, j);
            }
          }
          for (auto k : ms.holes)
          {
            const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
            if (SameOneBodyChannel(ok, oi))
            {
              r -= i7(k, i) * z.TwoBody.GetTBME_J(J, a, b, static_cast<int>(k), j);
            }
            if (SameOneBodyChannel(ok, oj))
            {
              r += ph_ij * i7(k, j) * z.TwoBody.GetTBME_J(J, a, b, static_cast<int>(k), i);
            }
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          dmat(ibra, iket) += r / nrm;
        }
      }
    }
  }

  int HalfIntPhase(int two_j_a, int two_j_b, int J)
  {
    return AngMom::phase(J) * AngMom::phase(std::abs((two_j_a - two_j_b) / 2));
  }

  // AMC z2_terms.tex eq 7–8 (unreduced). T2† Z1² Z2; not the naive I6/I7 R67 fold.
  void FoldGH(const Operator &t, const Operator &z, Operator &dz)
  {
    ModelSpace &ms = *t.modelspace;
    const size_t norb = ms.GetNumberOrbits();
    arma::mat i6h(norb, norb, arma::fill::zeros);
    arma::mat i7g(norb, norb, arma::fill::zeros);

    for (auto i : ms.holes)
    {
      const Orbit &oi = ms.GetOrbit(static_cast<int>(i));
      for (auto l : ms.holes)
      {
        const Orbit &ol = ms.GetOrbit(static_cast<int>(l));
        if (not SameOneBodyChannel(ol, oi))
        {
          continue;
        }
        double s = 0.0;
        for (auto c : ms.particles)
        {
          const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
          if (not SameOneBodyChannel(oc, oi))
          {
            continue;
          }
          const double z_ci = z.OneBody(static_cast<int>(c), static_cast<int>(i));
          if (std::abs(z_ci) < 1e-16)
          {
            continue;
          }
          for (auto d : ms.particles)
          {
            const Orbit &od = ms.GetOrbit(static_cast<int>(d));
            for (auto k : ms.holes)
            {
              const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
              if (not SameOneBodyChannel(ok, od))
              {
                continue;
              }
              const double z_dk = z.OneBody(static_cast<int>(d), static_cast<int>(k));
              if (std::abs(z_dk) < 1e-16)
              {
                continue;
              }
              const int Jmin = TwoBodyJmin(oc, od, ok, ol);
              const int Jmax = TwoBodyJmax(oc, od, ok, ol);
              for (int J = Jmin; J <= Jmax; ++J)
              {
                const int ph = AngMom::phase(J + (oi.j2 + od.j2) / 2);
                s += ph * (2 * J + 1)
                     * t.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(d), static_cast<int>(k), static_cast<int>(l))
                     * z_ci * z_dk;
              }
            }
          }
        }
        i7g(l, i) = s / (oi.j2 + 1.0);
      }
    }

    for (auto a : ms.particles)
    {
      const Orbit &oa = ms.GetOrbit(static_cast<int>(a));
      for (auto d : ms.particles)
      {
        const Orbit &od = ms.GetOrbit(static_cast<int>(d));
        if (not SameOneBodyChannel(od, oa))
        {
          continue;
        }
        double s = 0.0;
        for (auto k : ms.holes)
        {
          const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
          if (not SameOneBodyChannel(ok, oa))
          {
            continue;
          }
          const double z_ak = z.OneBody(static_cast<int>(a), static_cast<int>(k));
          if (std::abs(z_ak) < 1e-16)
          {
            continue;
          }
          for (auto c : ms.particles)
          {
            const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
            for (auto l : ms.holes)
            {
              const Orbit &ol = ms.GetOrbit(static_cast<int>(l));
              if (not SameOneBodyChannel(ol, oc))
              {
                continue;
              }
              const double z_cl = z.OneBody(static_cast<int>(c), static_cast<int>(l));
              if (std::abs(z_cl) < 1e-16)
              {
                continue;
              }
              const int Jmin = TwoBodyJmin(oc, od, ok, ol);
              const int Jmax = TwoBodyJmax(oc, od, ok, ol);
              for (int J = Jmin; J <= Jmax; ++J)
              {
                const int ph = AngMom::phase(J + (oa.j2 + oc.j2) / 2);
                s += ph * (2 * J + 1)
                     * t.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(d), static_cast<int>(k), static_cast<int>(l))
                     * z_ak * z_cl;
              }
            }
          }
        }
        i6h(d, a) = s / (oa.j2 + 1.0);
      }
    }

    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          double r = 0.0;
          for (auto l : ms.holes)
          {
            const int li = static_cast<int>(l);
            const Orbit &ol = ms.GetOrbit(li);
            if (SameOneBodyChannel(ol, oi))
            {
              r += i7g(li, i) * z.TwoBody.GetTBME_J(J, a, b, li, j);
            }
            if (SameOneBodyChannel(ol, oj))
            {
              r += HalfIntPhase(oi.j2, oj.j2, J) * i7g(li, j) * z.TwoBody.GetTBME_J(J, a, b, li, i);
            }
          }
          for (auto d : ms.particles)
          {
            const int di = static_cast<int>(d);
            const Orbit &od = ms.GetOrbit(di);
            if (SameOneBodyChannel(od, oa))
            {
              r += i6h(di, a) * z.TwoBody.GetTBME_J(J, di, b, i, j);
            }
            if (SameOneBodyChannel(od, ob))
            {
              r += HalfIntPhase(oa.j2, ob.j2, J) * i6h(di, b) * z.TwoBody.GetTBME_J(J, di, a, i, j);
            }
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          dmat(ibra, iket) += r / nrm;
        }
      }
    }
  }

  // AMC C,D (z2_terms.tex eq 3–4): 1-body intermediates with χ_hh-style phases, then fold.
  void FoldCD(const Operator &t, const Operator &z, Operator &dz)
  {
    ModelSpace &ms = *t.modelspace;
    const size_t norb = ms.GetNumberOrbits();
    arma::mat i6d(norb, norb, arma::fill::zeros);
    arma::mat i7c(norb, norb, arma::fill::zeros);

    for (auto b : ms.particles)
    {
      const Orbit &ob = ms.GetOrbit(static_cast<int>(b));
      for (auto c : ms.particles)
      {
        const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
        if (not SameOneBodyChannel(ob, oc))
        {
          continue;
        }
        double s2 = 0.0;
        for (auto d : ms.particles)
        {
          const Orbit &od = ms.GetOrbit(static_cast<int>(d));
          for (auto k : ms.holes)
          {
            const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
            for (auto l : ms.holes)
            {
              const Orbit &ol = ms.GetOrbit(static_cast<int>(l));
              const int Jmin = TwoBodyJmin(oc, od, ok, ol);
              const int Jmax = TwoBodyJmax(oc, od, ok, ol);
              for (int J = Jmin; J <= Jmax; ++J)
              {
                const int ph = AngMom::phase(J + (ob.j2 + od.j2) / 2);
                s2 += ph * (2 * J + 1)
                      * t.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(d), static_cast<int>(k), static_cast<int>(l))
                      * z.TwoBody.GetTBME_J(J, static_cast<int>(d), static_cast<int>(b), static_cast<int>(k), static_cast<int>(l));
              }
            }
          }
        }
        i6d(c, b) = 0.5 * s2 / (ob.j2 + 1.0);
      }
    }

    for (auto i : ms.holes)
    {
      const Orbit &oi = ms.GetOrbit(static_cast<int>(i));
      for (auto k : ms.holes)
      {
        const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
        if (not SameOneBodyChannel(ok, oi))
        {
          continue;
        }
        double s2 = 0.0;
        for (auto c : ms.particles)
        {
          const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
          for (auto d : ms.particles)
          {
            const Orbit &od = ms.GetOrbit(static_cast<int>(d));
            for (auto l : ms.holes)
            {
              const Orbit &ol = ms.GetOrbit(static_cast<int>(l));
              const int Jmin = TwoBodyJmin(oc, od, ok, ol);
              const int Jmax = TwoBodyJmax(oc, od, ok, ol);
              for (int J = Jmin; J <= Jmax; ++J)
              {
                const int ph = AngMom::phase(J + (oi.j2 + ol.j2) / 2);
                s2 += ph * (2 * J + 1)
                      * t.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(d), static_cast<int>(k), static_cast<int>(l))
                      * z.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(d), static_cast<int>(l), static_cast<int>(i));
              }
            }
          }
        }
        i7c(k, i) = 0.5 * s2 / (oi.j2 + 1.0);
      }
    }

    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          double r = 0.0;
          for (auto k : ms.holes)
          {
            const Orbit &ok = ms.GetOrbit(static_cast<int>(k));
            if (SameOneBodyChannel(ok, oj))
            {
              r += i7c(k, j) * z.TwoBody.GetTBME_J(J, a, b, i, static_cast<int>(k));
            }
            if (SameOneBodyChannel(ok, oi))
            {
              r += HalfIntPhase(oj.j2, oi.j2, J) * i7c(k, i) * z.TwoBody.GetTBME_J(J, a, b, j, static_cast<int>(k));
            }
          }
          for (auto c : ms.particles)
          {
            const Orbit &oc = ms.GetOrbit(static_cast<int>(c));
            if (SameOneBodyChannel(oc, ob))
            {
              r += i6d(c, b) * z.TwoBody.GetTBME_J(J, a, static_cast<int>(c), i, j);
            }
            if (SameOneBodyChannel(oc, oa))
            {
              r += HalfIntPhase(ob.j2, oa.j2, J) * i6d(c, a) * z.TwoBody.GetTBME_J(J, b, static_cast<int>(c), i, j);
            }
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          dmat(ibra, iket) += r / nrm;
        }
      }
    }
  }

  // AMC: I9 = 1/2 sum_{cd} t τ, R = 1/2 sum_{kl} τ I9 on unreduced GetTBME_J.
  // Stored (normalized) kets already absorb those 1/2's: T_st^T Tau_st = I9_st.
  void FoldI9(const Operator &t, const Operator &tau, Operator &dz)
  {
    ModelSpace &ms = *t.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      if (pp.n_rows == 0 or hh.n_rows == 0)
      {
        continue;
      }
      const arma::mat Tph = t.TwoBody.GetMatrix(ch, ch).submat(pp, hh);
      const arma::mat Tau = tau.TwoBody.GetMatrix(ch, ch).submat(pp, hh);
      dz.TwoBody.GetMatrix(ch, ch).submat(pp, hh) += Tau * (Tph.t() * Tau);
    }
  }

  void PandyaFill2n(const Operator &op, std::deque<arma::mat> &bar)
  {
    ModelSpace &ms = *op.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels_CC();
    bar.assign(nch, arma::mat());
    for (int ch = 0; ch < nch; ++ch)
    {
      const int n = ms.GetTwoBodyChannel_CC(ch).GetNumberKets();
      bar[ch].zeros(2 * n, 2 * n);
    }
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel_CC &tbc = ms.GetTwoBodyChannel_CC(ch);
      const int n = tbc.GetNumberKets();
      const int Jcc = tbc.J;
      for (int ibra = 0; ibra < 2 * n; ++ibra)
      {
        int a, b;
        if (ibra < n)
        {
          const Ket &bra = tbc.GetKet(ibra);
          a = bra.p;
          b = bra.q;
        }
        else
        {
          const Ket &bra = tbc.GetKet(ibra - n);
          b = bra.p;
          a = bra.q;
        }
        if (ibra >= n and a == b)
        {
          continue;
        }
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const double ja = 0.5 * oa.j2;
        const double jb = 0.5 * ob.j2;
        for (int iket = 0; iket < 2 * n; ++iket)
        {
          int c, d;
          if (iket < n)
          {
            const Ket &ket = tbc.GetKet(iket);
            c = ket.p;
            d = ket.q;
          }
          else
          {
            const Ket &ket = tbc.GetKet(iket - n);
            d = ket.p;
            c = ket.q;
          }
          if (iket >= n and c == d)
          {
            continue;
          }
          const Orbit &oc = ms.GetOrbit(c);
          const Orbit &od = ms.GetOrbit(d);
          const double jc = 0.5 * oc.j2;
          const double jd = 0.5 * od.j2;
          int Jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
          int Jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
          double xbar = 0.0;
          for (int J = Jmin; J <= Jmax; ++J)
          {
            const double sixj = ms.GetSixJ(ja, jb, Jcc, jc, jd, J);
            if (std::abs(sixj) < 1e-12)
            {
              continue;
            }
            xbar -= (2 * J + 1) * sixj * op.TwoBody.GetTBME_J(J, a, d, c, b);
          }
          bar[ch](ibra, iket) = xbar;
        }
      }
    }
  }

  void InversePandyaAdd2n(const std::deque<arma::mat> &bar, Operator &dz)
  {
    ModelSpace &ms = *dz.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const double ja = 0.5 * oa.j2;
        const double jb = 0.5 * ob.j2;
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          const double ji = 0.5 * oi.j2;
          const double jj = 0.5 * oj.j2;

          double commij = 0.0;
          int parity_cc = (oa.l + oj.l) % 2;
          int Tz_cc = std::abs(oa.tz2 - oj.tz2) / 2;
          int Jpmin = std::max(std::abs(oa.j2 - oj.j2), std::abs(oi.j2 - ob.j2)) / 2;
          int Jpmax = std::min(oa.j2 + oj.j2, oi.j2 + ob.j2) / 2;
          for (int Jp = Jpmin; Jp <= Jpmax; ++Jp)
          {
            const double sixj = ms.GetSixJ(ja, jb, J, ji, jj, Jp);
            if (std::abs(sixj) < 1e-12)
            {
              continue;
            }
            const int ch_cc = ms.GetTwoBodyChannelIndex(Jp, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = ms.GetTwoBodyChannel_CC(ch_cc);
            const int n = tbc_cc.GetNumberKets();
            const int loc_aj = static_cast<int>(tbc_cc.GetLocalIndex(std::min(a, j), std::max(a, j)));
            const int loc_ib = static_cast<int>(tbc_cc.GetLocalIndex(std::min(i, b), std::max(i, b)));
            if (loc_aj < 0 or loc_aj >= n or loc_ib < 0 or loc_ib >= n)
            {
              continue;
            }
            const int indx_aj = loc_aj + (a > j ? n : 0);
            const int indx_ib = loc_ib + (i > b ? n : 0);
            commij -= (2 * Jp + 1) * sixj * bar[ch_cc](indx_aj, indx_ib);
          }

          double commji = 0.0;
          parity_cc = (oa.l + oi.l) % 2;
          Tz_cc = std::abs(oa.tz2 - oi.tz2) / 2;
          Jpmin = std::max(std::abs(ob.j2 - oj.j2), std::abs(oi.j2 - oa.j2)) / 2;
          Jpmax = std::min(ob.j2 + oj.j2, oi.j2 + oa.j2) / 2;
          for (int Jp = Jpmin; Jp <= Jpmax; ++Jp)
          {
            const double sixj = ms.GetSixJ(jb, ja, J, ji, jj, Jp);
            if (std::abs(sixj) < 1e-12)
            {
              continue;
            }
            const int ch_cc = ms.GetTwoBodyChannelIndex(Jp, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = ms.GetTwoBodyChannel_CC(ch_cc);
            const int n = tbc_cc.GetNumberKets();
            const int loc_bi = static_cast<int>(tbc_cc.GetLocalIndex(std::min(b, i), std::max(b, i)));
            const int loc_ja = static_cast<int>(tbc_cc.GetLocalIndex(std::min(j, a), std::max(j, a)));
            if (loc_bi < 0 or loc_bi >= n or loc_ja < 0 or loc_ja >= n)
            {
              continue;
            }
            const int indx_bi = loc_bi + (b > i ? n : 0);
            const int indx_ja = loc_ja + (j > a ? n : 0);
            commji -= (2 * Jp + 1) * sixj * bar[ch_cc](indx_bi, indx_ja);
          }

          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          const double zijkl = (commij - ms.phase((oi.j2 + oj.j2) / 2 - J) * commji) / nrm;
          dmat(ibra, iket) += zijkl;
        }
      }
    }
  }

  void FoldPH(const Operator &t, const Operator &z, const Operator &w, Operator &dz)
  {
    std::deque<arma::mat> tbar, zbar, wbar;
    PandyaFill2n(t, tbar);
    PandyaFill2n(z, zbar);
    PandyaFill2n(w, wbar);
    const int nch = static_cast<int>(tbar.size());
    std::deque<arma::mat> rbar(nch);
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      // Non-Hermitian: full 2n×2n, no transpose fill.
      // 1/2 Z T Z (term E) + Z T W (term J).
      rbar[ch] = zbar[ch] * tbar[ch].t() * zbar[ch] + zbar[ch] * tbar[ch].t() * wbar[ch]
                 + wbar[ch] * tbar[ch].t() * zbar[ch];
    }
    InversePandyaAdd2n(rbar, dz);
  }

  // AMC z2_terms.tex eq 5 (unreduced). Four P(ab)P(ij) pieces. Ĵ² = 2J+1.
  void FoldE(const Operator &t, const Operator &z, Operator &dz)
  {
    ModelSpace &ms = *t.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J0 = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const double ja = 0.5 * oa.j2;
        const double jb = 0.5 * ob.j2;
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          const double ji = 0.5 * oi.j2;
          const double jj = 0.5 * oj.j2;
          double r = 0.0;
          for (auto c : ms.particles)
          {
            const int ci = static_cast<int>(c);
            const Orbit &oc = ms.GetOrbit(ci);
            const double jc = 0.5 * oc.j2;
            for (auto d : ms.particles)
            {
              const int di = static_cast<int>(d);
              const Orbit &od = ms.GetOrbit(di);
              const double jd = 0.5 * od.j2;
              for (auto k : ms.holes)
              {
                const int ki = static_cast<int>(k);
                const Orbit &ok = ms.GetOrbit(ki);
                const double jk = 0.5 * ok.j2;
                for (auto l : ms.holes)
                {
                  const int li = static_cast<int>(l);
                  const Orbit &ol = ms.GetOrbit(li);
                  const double jl = 0.5 * ol.j2;
                  const int jsum = (oc.j2 + od.j2 + ok.j2 + ol.j2) / 2;
                  const int J2min = TwoBodyJmin(oc, od, ok, ol);
                  const int J2max = TwoBodyJmax(oc, od, ok, ol);
                  for (int J2 = J2min; J2 <= J2max; ++J2)
                  {
                    const double t2 = t.TwoBody.GetTBME_J(J2, ci, di, ki, li);
                    if (std::abs(t2) < 1e-16)
                    {
                      continue;
                    }
                    const int hat2 = 2 * J2 + 1;

                    auto add_term = [&](int a3, int i3, int b4, int j4, const Orbit &oa3, const Orbit &oi3,
                                        const Orbit &ob4, const Orbit &oj4, double ja3, double ji3, double jb4,
                                        double jj4, double sh1, double sh2, double sp1, double sp2, int sign,
                                        int ph_outer) {
                      const int J3min = TwoBodyJmin(oa3, oc, oi3, ok);
                      const int J3max = TwoBodyJmax(oa3, oc, oi3, ok);
                      const int J4min = TwoBodyJmin(od, ob4, ol, oj4);
                      const int J4max = TwoBodyJmax(od, ob4, ol, oj4);
                      for (int J3 = J3min; J3 <= J3max; ++J3)
                      {
                        const double z3 = z.TwoBody.GetTBME_J(J3, a3, ci, i3, ki);
                        if (std::abs(z3) < 1e-16)
                        {
                          continue;
                        }
                        for (int J4 = J4min; J4 <= J4max; ++J4)
                        {
                          const double z4 = z.TwoBody.GetTBME_J(J4, di, b4, li, j4);
                          if (std::abs(z4) < 1e-16)
                          {
                            continue;
                          }
                          const int J5min = std::max(std::abs(oi3.j2 - oa3.j2), std::abs(oj4.j2 - ob4.j2)) / 2;
                          const int J5max = std::min(oi3.j2 + oa3.j2, oj4.j2 + ob4.j2) / 2;
                          for (int J5 = J5min; J5 <= J5max; ++J5)
                          {
                            const double s1 = ms.GetSixJ(jl, jd, J5, jc, jk, J2);
                            const double s2 = ms.GetSixJ(ji3, ja3, J5, jc, jk, J3);
                            const double s3 = ms.GetSixJ(jj4, jb4, J5, jd, jl, J4);
                            const double s4 = ms.GetSixJ(sh1, sh2, J0, sp1, sp2, J5);
                            const double six = s1 * s2 * s3 * s4;
                            if (std::abs(six) < 1e-16)
                            {
                              continue;
                            }
                            const int ph = AngMom::phase(ph_outer + J2 + J3 + J4 + jsum);
                            r += sign * ph * hat2 * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * six * t2 * z3 * z4;
                          }
                        }
                      }
                    };

                    add_term(a, i, b, j, oa, oi, ob, oj, ja, ji, jb, jj, ji, jj, jb, ja, +1, J0);
                    add_term(a, j, b, i, oa, oj, ob, oi, ja, jj, jb, ji, jj, ji, jb, ja, -1, (oi.j2 + oj.j2) / 2);
                    add_term(b, i, a, j, ob, oi, oa, oj, jb, ji, ja, jj, ji, jj, ja, jb, -1, (oa.j2 + ob.j2) / 2);
                    add_term(b, j, a, i, ob, oj, oa, oi, jb, jj, ja, ji, jj, ji, ja, jb, +1,
                             J0 + (oa.j2 + ob.j2 + oi.j2 + oj.j2) / 2);
                  }
                }
              }
            }
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          dmat(ibra, iket) += 0.5 * r / nrm;
        }
      }
    }
  }

  // AMC z2_terms.tex eq 10 (unreduced). No 1/2. Vacuous if Z1=0.
  void FoldJ(const Operator &t, const Operator &z, Operator &dz)
  {
    ModelSpace &ms = *t.modelspace;
    const int nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
      const int J0 = tbc.J;
      const auto &pp = tbc.GetKetIndex_pp();
      const auto &hh = tbc.GetKetIndex_hh();
      arma::mat &dmat = dz.TwoBody.GetMatrix(ch, ch);
      for (arma::uword ib = 0; ib < pp.n_rows; ++ib)
      {
        const int ibra = static_cast<int>(pp(ib));
        const Ket &bra = tbc.GetKet(ibra);
        const int a = bra.p;
        const int b = bra.q;
        const Orbit &oa = ms.GetOrbit(a);
        const Orbit &ob = ms.GetOrbit(b);
        const double ja = 0.5 * oa.j2;
        const double jb = 0.5 * ob.j2;
        for (arma::uword ik = 0; ik < hh.n_rows; ++ik)
        {
          const int iket = static_cast<int>(hh(ik));
          const Ket &ket = tbc.GetKet(iket);
          const int i = ket.p;
          const int j = ket.q;
          const Orbit &oi = ms.GetOrbit(i);
          const Orbit &oj = ms.GetOrbit(j);
          const double ji = 0.5 * oi.j2;
          const double jj = 0.5 * oj.j2;
          double r = 0.0;
          for (auto c : ms.particles)
          {
            const int ci = static_cast<int>(c);
            const Orbit &oc = ms.GetOrbit(ci);
            for (auto d : ms.particles)
            {
              const int di = static_cast<int>(d);
              const Orbit &od = ms.GetOrbit(di);
              const double jd = 0.5 * od.j2;
              for (auto k : ms.holes)
              {
                const int ki = static_cast<int>(k);
                const Orbit &ok = ms.GetOrbit(ki);
                for (auto l : ms.holes)
                {
                  const int li = static_cast<int>(l);
                  const Orbit &ol = ms.GetOrbit(li);
                  const double jl = 0.5 * ol.j2;
                  const int J2min = TwoBodyJmin(oc, od, ok, ol);
                  const int J2max = TwoBodyJmax(oc, od, ok, ol);
                  for (int J2 = J2min; J2 <= J2max; ++J2)
                  {
                    const double t2 = t.TwoBody.GetTBME_J(J2, ci, di, ki, li);
                    if (std::abs(t2) < 1e-16)
                    {
                      continue;
                    }

                    auto add_term = [&](int p_z1, int h_z1, int k_z1, int a_z1, int b4, int j4,
                                        const Orbit &ob4, const Orbit &oj4, double ja1, double ji1, double jb4,
                                        double jj4, double sh1, double sh2, double sp1, double sp2, int sign,
                                        int ph_outer) {
                      if (not SameOneBodyChannel(ms.GetOrbit(p_z1), ms.GetOrbit(h_z1)))
                      {
                        return;
                      }
                      if (not SameOneBodyChannel(ms.GetOrbit(a_z1), ms.GetOrbit(k_z1)))
                      {
                        return;
                      }
                      const double z_ci = z.OneBody(p_z1, h_z1);
                      const double z_ak = z.OneBody(a_z1, k_z1);
                      if (std::abs(z_ci * z_ak) < 1e-16)
                      {
                        return;
                      }
                      const int J3min = TwoBodyJmin(od, ob4, ol, oj4);
                      const int J3max = TwoBodyJmax(od, ob4, ol, oj4);
                      const int ja2 = static_cast<int>(std::lround(2.0 * ja1));
                      const int ji2 = static_cast<int>(std::lround(2.0 * ji1));
                      const int J4min = std::max(std::max(std::abs(ja2 - ji2), std::abs(oj4.j2 - ob4.j2)),
                                                 std::abs(od.j2 - ol.j2))
                                        / 2;
                      const int J4max = std::min(std::min(ja2 + ji2, oj4.j2 + ob4.j2), od.j2 + ol.j2) / 2;
                      for (int J3 = J3min; J3 <= J3max; ++J3)
                      {
                        const double z3 = z.TwoBody.GetTBME_J(J3, di, b4, li, j4);
                        if (std::abs(z3) < 1e-16)
                        {
                          continue;
                        }
                        for (int J4 = J4min; J4 <= J4max; ++J4)
                        {
                          const double s1 = ms.GetSixJ(ja1, ji1, J4, jd, jl, J2);
                          const double s2 = ms.GetSixJ(jj4, jb4, J4, jd, jl, J3);
                          const double s3 = ms.GetSixJ(sh1, sh2, J0, sp1, sp2, J4);
                          const double six = s1 * s2 * s3;
                          if (std::abs(six) < 1e-16)
                          {
                            continue;
                          }
                          const int ph = AngMom::phase(ph_outer + J2 + J3 + (od.j2 + ol.j2) / 2);
                          r += sign * ph * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * six * t2 * z_ci * z_ak * z3;
                        }
                      }
                    };

                    add_term(ci, i, ki, a, b, j, ob, oj, ja, ji, jb, jj, ji, jj, jb, ja, -1, J0);
                    add_term(ci, j, ki, a, b, i, ob, oi, ja, jj, jb, ji, jj, ji, jb, ja, +1, (oi.j2 + oj.j2) / 2);
                    add_term(ci, i, ki, b, a, j, oa, oj, jb, ji, ja, jj, ji, jj, ja, jb, +1, (oa.j2 + ob.j2) / 2);
                    add_term(ci, j, ki, b, a, i, oa, oi, jb, jj, ja, ji, jj, ji, ja, jb, -1,
                             J0 + (oa.j2 + ob.j2 + oi.j2 + oj.j2) / 2);
                  }
                }
              }
            }
          }
          const double nrm = (bra.delta_pq() ? PhysConst::SQRT2 : 1.0) * (ket.delta_pq() ? PhysConst::SQRT2 : 1.0);
          dmat(ibra, iket) += r / nrm;
        }
      }
    }
  }
}

EVC::EVC(ModelSpace &ms)
    : modelspace_(&ms)
{
}

Operator EVC::MakeExcitationOperatorLike(const Operator &op) const
{
  Operator out(*modelspace_, op.GetJRank(), op.GetTRank(), op.GetParity(), std::max(2, op.GetParticleRank()));
  out *= 0.0;
  out.SetNonHermitian();
  return out;
}

void EVC::PromoteToThreeBody(Operator &op) const
{
  if (op.GetParticleRank() < 3)
  {
    op.SetParticleRank(3);
  }
  op.ThreeBody.SetMode("pn");
}

void EVC::KeepExcitationBlocks(Operator &op, bool keep_zero_body) const
{
  double z0 = keep_zero_body ? op.ZeroBody : 0.0;
  Operator projected = MakeExcitationOperatorLike(op);
  projected.ZeroBody = z0;

  for (auto i : modelspace_->holes)
  {
    for (auto a : modelspace_->particles)
    {
      projected.OneBody(a, i) = op.OneBody(a, i);
    }
  }

  for (auto &iter : op.TwoBody.MatEl)
  {
    const size_t ch_bra = iter.first[0];
    const size_t ch_ket = iter.first[1];
    const TwoBodyChannel &tbc_bra = modelspace_->GetTwoBodyChannel(ch_bra);
    const TwoBodyChannel &tbc_ket = modelspace_->GetTwoBodyChannel(ch_ket);
    arma::mat &dst = projected.TwoBody.GetMatrix(ch_bra, ch_ket);
    const arma::mat &src = iter.second;
    for (auto ibra : tbc_bra.GetKetIndex_pp())
    {
      for (auto iket : tbc_ket.GetKetIndex_hh())
      {
        dst(ibra, iket) = src(ibra, iket);
      }
    }
  }

  op = std::move(projected);
}

void EVC::KeepDeexcitationBlocks(Operator &op) const
{
  Operator projected = MakeExcitationOperatorLike(op);

  for (auto i : modelspace_->holes)
  {
    for (auto a : modelspace_->particles)
    {
      projected.OneBody(i, a) = op.OneBody(i, a);
    }
  }

  for (auto &iter : op.TwoBody.MatEl)
  {
    const size_t ch_bra = iter.first[0];
    const size_t ch_ket = iter.first[1];
    const TwoBodyChannel &tbc_bra = modelspace_->GetTwoBodyChannel(ch_bra);
    const TwoBodyChannel &tbc_ket = modelspace_->GetTwoBodyChannel(ch_ket);
    arma::mat &dst = projected.TwoBody.GetMatrix(ch_bra, ch_ket);
    const arma::mat &src = iter.second;
    for (auto ibra : tbc_bra.GetKetIndex_hh())
    {
      for (auto iket : tbc_ket.GetKetIndex_pp())
      {
        dst(ibra, iket) = src(ibra, iket);
      }
    }
  }

  op = std::move(projected);
}

namespace
{
  Operator OnlyOneBody(const Operator &op)
  {
    Operator out(op);
    out.EraseZeroBody();
    out.EraseTwoBody();
    out.SetNonHermitian();
    return out;
  }

  Operator OnlyTwoBody(const Operator &op)
  {
    Operator out(op);
    out.EraseZeroBody();
    out.EraseOneBody();
    out.SetNonHermitian();
    return out;
  }

  Operator BuildExplicitCommutator(const Operator &X, const Operator &Y)
  {
    ModelSpace &ms = *X.GetModelSpace();
    int jrank = X.GetJRank() + Y.GetJRank();
    int trank = X.GetTRank() + Y.GetTRank();
    int parity = (X.GetParity() + Y.GetParity()) % 2;
    Operator Z(ms, jrank, trank, parity, 2);
    Z.SetNonHermitian();
    Z.Erase();

    Operator X1 = OnlyOneBody(X);
    Operator X2 = OnlyTwoBody(X);
    Operator Y1 = OnlyOneBody(Y);
    Operator Y2 = OnlyTwoBody(Y);

    ReferenceImplementations::comm110ss(X1, Y1, Z);
    ReferenceImplementations::comm220ss(X2, Y2, Z);

    ReferenceImplementations::comm111ss(X1, Y1, Z);
    ReferenceImplementations::comm121ss(X1, Y2, Z);

    Operator tmp(ms, jrank, trank, parity, 2);
    tmp.SetNonHermitian();
    tmp.Erase();
    ReferenceImplementations::comm121ss(Y1, X2, tmp);
    Z -= tmp;

    ReferenceImplementations::comm221ss(X2, Y2, Z);

    ReferenceImplementations::comm122ss(X1, Y2, Z);
    tmp.Erase();
    ReferenceImplementations::comm122ss(Y1, X2, tmp);
    Z -= tmp;

    ReferenceImplementations::comm222_pp_hhss(X2, Y2, Z);
    ReferenceImplementations::comm222_phss(X2, Y2, Z);

    return Z;
  }
}

Operator EVC::ExtractExcitationPart(const Operator &op) const
{
  Operator out = op;
  out.SetNonHermitian();
  KeepExcitationBlocks(out, false);
  return out;
}

Operator EVC::ExtractDeexcitationPart(const Operator &op) const
{
  Operator out = op;
  out.SetNonHermitian();
  KeepDeexcitationBlocks(out);
  return out;
}

Operator EVC::ExtractTFromOmega(const Operator &omega) const
{
  return ExtractExcitationPart(omega);
}

Operator EVC::ExtractTDaggerFromOmega(const Operator &omega) const
{
  Operator tdagger = ExtractDeexcitationPart(omega);
  tdagger *= -1.0;
  return tdagger;
}

Operator EVC::BuildClusterOperator(const ClusterAmplitudes &amps) const
{
  Operator z = amps.z;
  z.ZeroBody = amps.z0;
  z.SetNonHermitian();
  KeepExcitationBlocks(z, true);
  return z;
}

Operator EVC::BuildZ1RHS(const Operator &t_in, const Operator &z_in) const
{
  double t_start = omp_get_wtime();
  Operator t = t_in;
  Operator z = z_in;
  t.SetNonHermitian();
  z.SetNonHermitian();
  KeepExcitationBlocks(t, false);
  KeepExcitationBlocks(z, false);

  Operator dz = MakeExcitationOperatorLike(z);
  dz.SetNonHermitian();
  dz.Erase();

  const std::vector<index_t> holes(modelspace_->holes.begin(), modelspace_->holes.end());
  const std::vector<index_t> particles(modelspace_->particles.begin(), modelspace_->particles.end());

  // CCM-style 1-body intermediates. AMC: run/evc/amc_eqs/output/z1_intermediates.tex
  const arma::mat chi_hp = BuildChiHP(t, z, holes, particles);
  const arma::mat chi_hh = BuildChiHH(t, z, holes, particles);
  const arma::mat chi_pp = BuildChiPP(t, z, holes, particles);

  std::vector<std::pair<index_t, index_t>> pairs;
  for (auto i : holes)
  {
    const Orbit &oi = modelspace_->GetOrbit(static_cast<int>(i));
    for (auto a : particles)
    {
      const Orbit &oa = modelspace_->GetOrbit(static_cast<int>(a));
      if (SameOneBodyChannel(oa, oi))
      {
        pairs.emplace_back(a, i);
      }
    }
  }

  const int npairs = static_cast<int>(pairs.size());
#pragma omp parallel for schedule(dynamic)
  for (int ip = 0; ip < npairs; ++ip)
  {
    const index_t a = pairs[ip].first;
    const index_t i = pairs[ip].second;
    const Orbit &oa = modelspace_->GetOrbit(static_cast<int>(a));
    const Orbit &oi = modelspace_->GetOrbit(static_cast<int>(i));
    const double hatj2 = oa.j2 + 1.0; // AMC ĵ_a^{-2} = 1/(2j_a+1)
    double r = 0.0;

    // Term 1: R_ai = ĵ_a^{-2} Σ_{bj J} Ĵ^{2} t_bj z^J_abij
    for (auto j : holes)
    {
      const Orbit &oj = modelspace_->GetOrbit(static_cast<int>(j));
      for (auto b : particles)
      {
        const Orbit &ob = modelspace_->GetOrbit(static_cast<int>(b));
        if (not SameOneBodyChannel(ob, oj))
        {
          continue;
        }
        const double t1 = t.OneBody(b, j);
        if (std::abs(t1) < 1e-16)
        {
          continue;
        }
        const int Jmin = TwoBodyJmin(oa, ob, oi, oj);
        const int Jmax = TwoBodyJmax(oa, ob, oi, oj);
        for (int J = Jmin; J <= Jmax; ++J)
        {
          r += (2 * J + 1) * t1 * z.TwoBody.GetTBME_J(J, static_cast<int>(a), static_cast<int>(b), static_cast<int>(i), static_cast<int>(j));
        }
      }
    }
    r /= hatj2;

    // Terms 2+6: I5_bj = t_bj + χ_hp_bj, then - I5_bj z_bi z_aj
    for (auto j : holes)
    {
      const Orbit &oj = modelspace_->GetOrbit(static_cast<int>(j));
      if (not SameOneBodyChannel(oj, oa))
      {
        continue;
      }
      for (auto b : particles)
      {
        const Orbit &ob = modelspace_->GetOrbit(static_cast<int>(b));
        if (not SameOneBodyChannel(ob, oa))
        {
          continue;
        }
        r -= (t.OneBody(b, j) + chi_hp(b, j)) * z.OneBody(b, i) * z.OneBody(a, j);
      }
    }

    // Term 3: ĵ_a^{-2} Σ_{ck J} Ĵ^{2} χ_hp_ck z^J_caki
    double r3 = 0.0;
    for (auto k : holes)
    {
      const Orbit &ok = modelspace_->GetOrbit(static_cast<int>(k));
      for (auto c : particles)
      {
        const Orbit &oc = modelspace_->GetOrbit(static_cast<int>(c));
        if (not SameOneBodyChannel(oc, ok))
        {
          continue;
        }
        const double ch = chi_hp(c, k);
        if (std::abs(ch) < 1e-16)
        {
          continue;
        }
        const int Jmin = TwoBodyJmin(oc, oa, ok, oi);
        const int Jmax = TwoBodyJmax(oc, oa, ok, oi);
        for (int J = Jmin; J <= Jmax; ++J)
        {
          r3 += (2 * J + 1) * ch * z.TwoBody.GetTBME_J(J, static_cast<int>(c), static_cast<int>(a), static_cast<int>(k), static_cast<int>(i));
        }
      }
    }
    r += r3 / hatj2;

    // Term 4: Σ_b χ_pp_ba z_bi
    for (auto b : particles)
    {
      const Orbit &ob = modelspace_->GetOrbit(static_cast<int>(b));
      if (not SameOneBodyChannel(ob, oa))
      {
        continue;
      }
      r += chi_pp(b, a) * z.OneBody(b, i);
    }

    // Term 5: Σ_j z_aj χ_hh_ji
    for (auto j : holes)
    {
      const Orbit &oj = modelspace_->GetOrbit(static_cast<int>(j));
      if (not SameOneBodyChannel(oj, oa))
      {
        continue;
      }
      r += z.OneBody(a, j) * chi_hh(j, i);
    }

    dz.OneBody(a, i) = -r;
  }

  t_in.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return dz;
}

Operator EVC::BuildZ2RHS(const Operator &t_in, const Operator &z_in) const
{
  double t_start = omp_get_wtime();
  Operator t = t_in;
  Operator z = z_in;
  t.SetNonHermitian();
  z.SetNonHermitian();
  KeepExcitationBlocks(t, false);
  KeepExcitationBlocks(z, false);

  Operator dz = MakeExcitationOperatorLike(z);
  dz.SetNonHermitian();
  dz.Erase();

  const std::vector<index_t> holes(modelspace_->holes.begin(), modelspace_->holes.end());
  const std::vector<index_t> particles(modelspace_->particles.begin(), modelspace_->particles.end());

  const arma::mat chi_pp = BuildChiPP(t, z, holes, particles);
  const arma::mat chi_hh = BuildChiHH(t, z, holes, particles);
  const arma::mat i6 = BuildI6(t, z, chi_pp, holes, particles);
  const arma::mat i7 = BuildI7(t, z, chi_hh, holes, particles);
  FoldI6I7(i6, i7, z, dz); // A,B
  FoldGH(t, z, dz);        // G,H
  FoldCD(t, z, dz);        // C,D

  Operator w = BuildWFromZ1(z, holes, particles);
  Operator tau = z;
  tau.TwoBody += w.TwoBody;
  FoldI9(t, tau, dz); // F,I,K,L
  FoldE(t, z, dz);    // E
  FoldJ(t, z, dz);    // J

  dz *= -1.0;
  t_in.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return dz;
}

double EVC::ConnectedOV_OOVV(const Operator &op, const Operator &z) const
{
  double r = 0.0;
  for (auto i : modelspace_->holes)
  {
    const Orbit &oi = modelspace_->GetOrbit(static_cast<int>(i));
    for (auto a : modelspace_->particles)
    {
      const Orbit &oa = modelspace_->GetOrbit(static_cast<int>(a));
      if (not SameOneBodyChannel(oa, oi))
      {
        continue;
      }
      r += (oa.j2 + 1.0) * op.OneBody(a, i) * z.OneBody(a, i);
    }
  }

  for (auto a : modelspace_->particles)
  {
    const Orbit &oa = modelspace_->GetOrbit(static_cast<int>(a));
    for (auto b : modelspace_->particles)
    {
      const Orbit &ob = modelspace_->GetOrbit(static_cast<int>(b));
      for (auto i : modelspace_->holes)
      {
        const Orbit &oi = modelspace_->GetOrbit(static_cast<int>(i));
        for (auto j : modelspace_->holes)
        {
          const Orbit &oj = modelspace_->GetOrbit(static_cast<int>(j));
          const int Jmin = TwoBodyJmin(oa, ob, oi, oj);
          const int Jmax = TwoBodyJmax(oa, ob, oi, oj);
          for (int J = Jmin; J <= Jmax; ++J)
          {
            const double ojme = op.TwoBody.GetTBME_J(J, static_cast<int>(a), static_cast<int>(b), static_cast<int>(i), static_cast<int>(j));
            const double hatJ2 = 2 * J + 1.0;
            r += 0.25 * hatJ2 * ojme
                 * z.TwoBody.GetTBME_J(J, static_cast<int>(a), static_cast<int>(b), static_cast<int>(i), static_cast<int>(j));
            if (SameOneBodyChannel(oa, oi) and SameOneBodyChannel(ob, oj))
            {
              r += 0.5 * hatJ2 * ojme * z.OneBody(a, i) * z.OneBody(b, j);
            }
          }
        }
      }
    }
  }
  return r;
}

double EVC::BuildZ0RHS(const Operator &t_in, const Operator &z_in) const
{
  double t_start = omp_get_wtime();
  Operator t = t_in;
  Operator z = z_in;
  t.SetNonHermitian();
  z.SetNonHermitian();
  KeepExcitationBlocks(t, false);
  KeepExcitationBlocks(z, false);
  const double dz0 = -ConnectedOV_OOVV(t, z);
  t_in.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return dz0;
}

double EVC::CoupledClusterEnergy(const Operator &h, const Operator &z_in) const
{
  // Additive H_0 from the unitary BCH only. Z_0 is not used here.
  const double h0 = h.ZeroBody;
  Operator z = z_in;
  z.SetNonHermitian();
  z.EraseZeroBody();
  KeepExcitationBlocks(z, false);
  return h0 + ConnectedOV_OOVV(h, z);
}

double EVC::HamiltonianKernel(const Operator &h, const ClusterAmplitudes &amps) const
{
  // <0|H e^Z|0> = e^{Z_0} (H_0 + Δ_CC(H, Z_1, Z_2)). Never H_0+Z_0 or BCH(H,Z).
  return std::exp(amps.z0) * CoupledClusterEnergy(h, amps.z);
}

double EVC::NormKernel(const ClusterAmplitudes &amps) const
{
  return std::exp(amps.z0);
}

std::array<double, 2> EVC::EvaluateKernels(const Operator &h0, const Operator &omega_i, const Operator &omega_k) const
{
  Operator wi = omega_i;
  Operator wk = omega_k;
  // Magnus Ω_0 ≡ 0. Junk here would be a fake e^{Ω_0} on top of the ODE's e^{Z_0}.
  wi.EraseZeroBody();
  wk.EraseZeroBody();
  Operator neg_i = wi;
  neg_i *= -1.0;
  Operator h_i = BCH::BCH_Transform(h0, neg_i); // keep h_i.ZeroBody = H^{(i)}_0
  Operator omega_d = BCH::BCH_Product(neg_i, wk);
  omega_d.EraseZeroBody();
  ClusterAmplitudes amps = SolveFromOmega(omega_d);
  const double n_ik = NormKernel(amps);
  const double h_ik = HamiltonianKernel(h_i, amps);
  return {h_ik, n_ik};
}

Operator EVC::SimilarityTransformCCSD(const Operator &tdagger_in, const Operator &z_in) const
{
  Operator tdagger = tdagger_in;
  Operator z = z_in;
  tdagger.SetNonHermitian();
  z.SetNonHermitian();
  KeepDeexcitationBlocks(tdagger);
  KeepExcitationBlocks(z, true);

  Operator transformed = tdagger;
  transformed.SetNonHermitian();

  Operator c1 = BuildExplicitCommutator(tdagger, z);
  c1.SetNonHermitian();
  transformed += c1;

  Operator c2 = BuildExplicitCommutator(c1, z);
  c2.SetNonHermitian();
  transformed += kHalf * c2;

  Operator c3 = BuildExplicitCommutator(c2, z);
  c3.SetNonHermitian();
  transformed += kSixth * c3;

  Operator c4 = BuildExplicitCommutator(c3, z);
  c4.SetNonHermitian();
  transformed += kTwentyFourth * c4;

  return transformed;
}

EVC::RHSDiagrams EVC::BuildRHSDiagrams(const Operator &tdagger_in, const Operator &z_in) const
{
  RHSDiagrams rhs(*modelspace_);

  Operator tdagger = tdagger_in;
  Operator z = z_in;
  tdagger.SetNonHermitian();
  z.SetNonHermitian();
  KeepDeexcitationBlocks(tdagger);
  KeepExcitationBlocks(z, true);
  rhs.c1 = BuildExplicitCommutator(tdagger, z);
  rhs.c1.SetNonHermitian();
  rhs.c2 = BuildExplicitCommutator(rhs.c1, z);
  rhs.c2.SetNonHermitian();
  rhs.c3 = BuildExplicitCommutator(rhs.c2, z);
  rhs.c3.SetNonHermitian();
  rhs.c4 = BuildExplicitCommutator(rhs.c3, z);
  rhs.c4.SetNonHermitian();

  rhs.dz0 = -(tdagger.ZeroBody + rhs.c1.ZeroBody + kHalf * rhs.c2.ZeroBody + kSixth * rhs.c3.ZeroBody + kTwentyFourth * rhs.c4.ZeroBody);

  KeepExcitationBlocks(rhs.c1, false);
  KeepExcitationBlocks(rhs.c2, false);
  KeepExcitationBlocks(rhs.c3, false);
  KeepExcitationBlocks(rhs.c4, false);

  rhs.c1 *= -1.0;
  rhs.c2 *= -kHalf;
  rhs.c3 *= -kSixth;
  rhs.c4 *= -kTwentyFourth;
  return rhs;
}

EVC::RHS EVC::BuildRHS(const Operator &t, const ClusterAmplitudes &amps) const
{
  RHS rhs(*modelspace_);
  Operator z = amps.z;
  z.SetNonHermitian();
  KeepExcitationBlocks(z, false);

  rhs.dz = BuildZ1RHS(t, z);
  Operator dz2 = BuildZ2RHS(t, z);
  rhs.dz.TwoBody += dz2.TwoBody;
  rhs.dz.SetNonHermitian();
  KeepExcitationBlocks(rhs.dz, false);

  rhs.dz0 = BuildZ0RHS(t, z);
  rhs.transformed = rhs.dz;
  rhs.transformed *= -1.0;
  return rhs;
}

void EVC::AddRHS(ClusterAmplitudes &amps, const RHS &rhs, double scale) const
{
  amps.z0 += scale * rhs.dz0;
  amps.z += scale * rhs.dz;
  amps.z.SetNonHermitian();
  KeepExcitationBlocks(amps.z, false);
}

EVC::ClusterAmplitudes EVC::Integrate(const Operator &t_in, double lambda_max, bool rk4) const
{
  double t_start = omp_get_wtime();
  Operator t = ExtractExcitationPart(t_in);

  ClusterAmplitudes amps(*modelspace_);
  amps.z = t;
  amps.z0 = 0.0;
  amps.z.ZeroBody = 0.0;

  const std::size_t nsteps = euler_steps_;
  if (nsteps == 0 or std::abs(lambda_max) < 1e-16)
  {
    t_in.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return amps;
  }

  const double h = lambda_max / static_cast<double>(nsteps);
  for (std::size_t istep = 0; istep < nsteps; ++istep)
  {
    if (not rk4)
    {
      RHS k = BuildRHS(t, amps);
      AddRHS(amps, k, h);
      continue;
    }
    RHS k1 = BuildRHS(t, amps);

    ClusterAmplitudes z2 = amps;
    AddRHS(z2, k1, 0.5 * h);
    RHS k2 = BuildRHS(t, z2);

    ClusterAmplitudes z3 = amps;
    AddRHS(z3, k2, 0.5 * h);
    RHS k3 = BuildRHS(t, z3);

    ClusterAmplitudes z4 = amps;
    AddRHS(z4, k3, h);
    RHS k4 = BuildRHS(t, z4);

    AddRHS(amps, k1, h / 6.0);
    AddRHS(amps, k2, h / 3.0);
    AddRHS(amps, k3, h / 3.0);
    AddRHS(amps, k4, h / 6.0);
  }

  // Z_0 stays in amps.z0. Do not copy it onto z.ZeroBody (that slot is BCH/H only).
  amps.z.EraseZeroBody();
  t_in.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return amps;
}

EVC::ClusterAmplitudes EVC::Solve(const Operator &t, double lambda_max) const
{
  return Integrate(t, lambda_max, use_rk4_);
}

EVC::ClusterAmplitudes EVC::SolveEuler(const Operator &t, const Operator &tdagger, double lambda_max) const
{
  (void)tdagger;
  return Integrate(t, lambda_max, false);
}

EVC::ClusterAmplitudes EVC::SolveFromOmega(const Operator &omega, double lambda_max) const
{
  return Solve(ExtractTFromOmega(omega), lambda_max);
}

EVC::ClusterAmplitudes EVC::SolveEulerFromOmega(const Operator &omega, double lambda_max) const
{
  Operator t = ExtractTFromOmega(omega);
  Operator tdagger = ExtractTDaggerFromOmega(omega);
  return SolveEuler(t, tdagger, lambda_max);
}