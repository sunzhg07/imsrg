#include "FactorizedDoubleCommutator.hh"
#include "Commutator.hh"
#include "TensorCommutators.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Commutator
{

namespace FactorizedDoubleCommutator
{

namespace
{

using ChannelPair = std::array<index_t, 2>;
using TwoBodyChannelPair = std::array<size_t, 2>;
using PandyaMap = std::map<ChannelPair, arma::mat>;
using TensorPandyaFullMap = std::map<ChannelPair, arma::mat>;

double Hat(int J)
{
  return std::sqrt(2.0 * J + 1.0);
}

double HatOrbit(const Orbit &orb)
{
  return std::sqrt(orb.j2 + 1.0);
}

int Phase(int exponent)
{
  return AngMom::phase(exponent);
}

int HalfIntegerPhaseExponent(std::initializer_list<int> two_j_values)
{
  int sum = 0;
  for (int two_j : two_j_values)
    sum += two_j;
  return sum / 2;
}

double ScalarCoupledTensorProductFactor(int Jout, int Jint, int K)
{
  if ((Jout + Jint < K) or (std::abs(Jout - Jint) > K))
    return 0.0;
  return Phase(Jout + Jint + K) / (Hat(Jout) * Hat(K));
}

bool AllowExperimentalTensor231()
{
  const char *flag = std::getenv("IMSRG_ALLOW_EXPERIMENTAL_TENSOR_223231");
  if (flag == nullptr)
    return false;
  const std::string value(flag);
  return value == "1" or value == "true" or value == "TRUE" or value == "on" or value == "ON";
}

double TensorOneBodyTensorTwoBodyToScalarFactor(const Operator &Z, int Jout, int Jint, int K,
                                               const Orbit &changed, const Orbit &replacement, const Orbit &spectator)
{
  if ((Jout + Jint < K) or (std::abs(Jout - Jint) > K))
    return 0.0;
  double sixj = Z.modelspace->GetSixJ(Jint, Jout, K,
                                      changed.j2 * 0.5, replacement.j2 * 0.5, spectator.j2 * 0.5);
  return Phase(Jout + K + HalfIntegerPhaseExponent({changed.j2, spectator.j2})) * Hat(Jout) * Hat(Jint) * sixj / Hat(K);
}

int JMin(const Orbit &a, const Orbit &b, const Orbit &c, const Orbit &d)
{
  return std::max(std::abs(a.j2 - b.j2), std::abs(c.j2 - d.j2)) / 2;
}

int JMax(const Orbit &a, const Orbit &b, const Orbit &c, const Orbit &d)
{
  return std::min(a.j2 + b.j2, c.j2 + d.j2) / 2;
}

void CheckTensorEtaInputs(const Operator &Eta, const Operator &Gamma, const Operator &Z, const std::string &where)
{
  if (Eta.GetJRank() == 0)
  {
    throw std::invalid_argument(where + " expects tensor Eta with Jrank > 0. Use the scalar factorized routines for scalar Eta.");
  }
  if (Gamma.GetJRank() != 0)
  {
    throw std::invalid_argument(where + " expects scalar Gamma with Jrank == 0.");
  }
  if (Z.GetJRank() != 0)
  {
    throw std::invalid_argument(where + " currently codes the scalar-coupled tensor-eta double commutator, so Z must have Jrank == 0.");
  }
  if (not Eta.IsReduced())
  {
    throw std::invalid_argument(where + " expects Eta in reduced tensor storage.");
  }
}

bool TensorChannelPairAllowed(const Operator &Op, const TwoBodyChannel_CC &bra, const TwoBodyChannel_CC &ket)
{
  if ((bra.J + ket.J < Op.GetJRank()) or (std::abs(bra.J - ket.J) > Op.GetJRank()))
    return false;
  if ((bra.parity + ket.parity + Op.GetParity()) % 2 > 0)
    return false;
  if (not((bra.Tz + ket.Tz == Op.GetTRank()) or (std::abs(bra.Tz - ket.Tz) == Op.GetTRank())))
    return false;
  return true;
}

PandyaMap AllocateTensorPandyaLike(const Operator &Op)
{
  PandyaMap mats;
  for (auto ch_bra_cc : Op.modelspace->SortedTwoBodyChannels_CC)
  {
    const auto &tbc_bra_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
    arma::uvec bras_ph = arma::join_cols(tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph());
    for (auto ch_ket_cc : Op.modelspace->SortedTwoBodyChannels_CC)
    {
      const auto &tbc_ket_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
      if (not TensorChannelPairAllowed(Op, tbc_bra_cc, tbc_ket_cc))
        continue;
      mats[{ch_bra_cc, ch_ket_cc}] = arma::mat(2 * bras_ph.n_rows, tbc_ket_cc.GetNumberKets(), arma::fill::zeros);
    }
  }
  return mats;
}

std::vector<TwoBodyChannelPair> StoredTwoBodyChannelPairs(const Operator &Op)
{
  std::vector<TwoBodyChannelPair> channel_pairs;
  channel_pairs.reserve(Op.TwoBody.MatEl.size());
  for (const auto &itmat : Op.TwoBody.MatEl)
    channel_pairs.push_back(itmat.first);
  return channel_pairs;
}

void BuildTensorPandya(const Operator &Op, PandyaMap &barOp)
{
  barOp.clear();
  DoTensorPandyaTransformation(Op, barOp);
}

void BuildFullSquareTensorPandya(const Operator &Op, TensorPandyaFullMap &barOp)
{
  barOp.clear();
  const int K = Op.GetJRank();

  for (auto ch_bra_cc : Op.modelspace->SortedTwoBodyChannels_CC)
  {
    TwoBodyChannel_CC &tbc_bra_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
    const int Jbra_cc = tbc_bra_cc.J;
    const size_t nbras_cc = tbc_bra_cc.GetNumberKets();

    for (auto ch_ket_cc : Op.modelspace->SortedTwoBodyChannels_CC)
    {
      TwoBodyChannel_CC &tbc_ket_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
      const int Jket_cc = tbc_ket_cc.J;
      const size_t nkets_cc = tbc_ket_cc.GetNumberKets();

      if ((Jbra_cc + Jket_cc < K) or (std::abs(Jbra_cc - Jket_cc) > K))
        continue;
      if ((tbc_bra_cc.parity + tbc_ket_cc.parity + Op.GetParity()) % 2 > 0)
        continue;
      if (nbras_cc == 0 or nkets_cc == 0)
        continue;

      arma::mat &bar = barOp[{ch_bra_cc, ch_ket_cc}];
      bar.zeros(2 * nbras_cc, 2 * nkets_cc);

      for (size_t ibra_cc = 0; ibra_cc < 2 * nbras_cc; ++ibra_cc)
      {
        Ket &bra_cc = tbc_bra_cc.GetKet(ibra_cc % nbras_cc);
        const index_t a = ibra_cc < nbras_cc ? bra_cc.p : bra_cc.q;
        const index_t b = ibra_cc < nbras_cc ? bra_cc.q : bra_cc.p;
        Orbit &oa = Op.modelspace->GetOrbit(a);
        Orbit &ob = Op.modelspace->GetOrbit(b);
        const double ja = oa.j2 * 0.5;
        const double jb = ob.j2 * 0.5;

        for (size_t iket_cc = 0; iket_cc < 2 * nkets_cc; ++iket_cc)
        {
          Ket &ket_cc = tbc_ket_cc.GetKet(iket_cc % nkets_cc);
          const index_t c = iket_cc < nkets_cc ? ket_cc.p : ket_cc.q;
          const index_t d = iket_cc < nkets_cc ? ket_cc.q : ket_cc.p;
          Orbit &oc = Op.modelspace->GetOrbit(c);
          Orbit &od = Op.modelspace->GetOrbit(d);
          const double jc = oc.j2 * 0.5;
          const double jd = od.j2 * 0.5;

          double sum = 0.0;
          int J1min = std::abs(ja - jd);
          int J1max = ja + jd;
          for (int J1 = J1min; J1 <= J1max; ++J1)
          {
            int J2min = std::max(int(std::abs(jc - jb)), std::abs(J1 - K));
            int J2max = std::min(int(jc + jb), J1 + K);
            for (int J2 = J2min; J2 <= J2max; ++J2)
            {
              double recoupling = 0.0;
              if (K == 0)
              {
                recoupling = AngMom::phase(jb + jd + J1 + Jbra_cc) * Op.modelspace->GetSixJ(ja, jb, Jbra_cc, jc, jd, J1) / std::sqrt((2 * J2 + 1) * (2 * Jbra_cc + 1));
              }
              else
              {
                recoupling = Op.modelspace->GetNineJ(ja, jd, J1, jb, jc, J2, Jbra_cc, Jket_cc, K);
              }
              if (std::abs(recoupling) < 1e-12)
                continue;

              double hatfactor = std::sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * Jbra_cc + 1) * (2 * Jket_cc + 1));
              double tbme = Op.TwoBody.GetTBME_J(J1, J2, a, d, c, b);
              sum -= hatfactor * AngMom::phase(jb + jd + Jket_cc + J2) * recoupling * tbme;
            }
          }
          bar(ibra_cc, iket_cc) = sum;
        }
      }
    }
  }
}

void BuildFullSquareTensorPandyaOccWeighted(const Operator &Op, int occ_type, TensorPandyaFullMap &barOp)
{
  BuildFullSquareTensorPandya(Op, barOp);

  for (auto &itmat : barOp)
  {
    TwoBodyChannel_CC &tbc_bra_cc = Op.modelspace->GetTwoBodyChannel_CC(itmat.first[0]);
    TwoBodyChannel_CC &tbc_ket_cc = Op.modelspace->GetTwoBodyChannel_CC(itmat.first[1]);
    const size_t nbras_cc = tbc_bra_cc.GetNumberKets();
    const size_t nkets_cc = tbc_ket_cc.GetNumberKets();
    arma::mat &bar = itmat.second;

    for (size_t ibra_cc = 0; ibra_cc < 2 * nbras_cc; ++ibra_cc)
    {
      Ket &bra_cc = tbc_bra_cc.GetKet(ibra_cc % nbras_cc);
      const index_t a = ibra_cc < nbras_cc ? bra_cc.p : bra_cc.q;
      const index_t b = ibra_cc < nbras_cc ? bra_cc.q : bra_cc.p;
      Orbit &oa = Op.modelspace->GetOrbit(a);
      Orbit &ob = Op.modelspace->GetOrbit(b);
      const double n_a = oa.occ;
      const double n_b = ob.occ;
      const double nbar_a = 1.0 - n_a;
      const double nbar_b = 1.0 - n_b;

      for (size_t iket_cc = 0; iket_cc < 2 * nkets_cc; ++iket_cc)
      {
        Ket &ket_cc = tbc_ket_cc.GetKet(iket_cc % nkets_cc);
        const index_t c = iket_cc < nkets_cc ? ket_cc.p : ket_cc.q;
        const index_t d = iket_cc < nkets_cc ? ket_cc.q : ket_cc.p;
        Orbit &oc = Op.modelspace->GetOrbit(c);
        Orbit &od = Op.modelspace->GetOrbit(d);
        const double n_c = oc.occ;
        const double n_d = od.occ;
        const double nbar_c = 1.0 - n_c;
        const double nbar_d = 1.0 - n_d;
        double occfactor = 0.0;
        if (occ_type == 0)
          occfactor = nbar_a * n_b * n_c + n_a * nbar_b * nbar_c;
        else
          occfactor = n_a * nbar_b * n_d + nbar_a * n_b * nbar_d;
        bar(ibra_cc, iket_cc) *= occfactor;
      }
    }
  }
}

void AddScalarCoupledTensorProductFullSquare(const Operator &Eta,
                                             const TensorPandyaFullMap &barLeft,
                                             const TensorPandyaFullMap &barRight,
                                             TensorPandyaFullMap &barZ,
                                             double scale)
{
  const int K = Eta.GetJRank();
  for (auto ch0_cc : Eta.modelspace->SortedTwoBodyChannels_CC)
  {
    TwoBodyChannel_CC &tbc0 = Eta.modelspace->GetTwoBodyChannel_CC(ch0_cc);
    const int J0 = tbc0.J;
    const size_t n0 = tbc0.GetNumberKets();
    if (n0 == 0)
      continue;
    arma::mat &accum = barZ[{ch0_cc, ch0_cc}];
    if (accum.n_rows == 0)
      accum.zeros(2 * n0, 2 * n0);

    for (auto ch2_cc : Eta.modelspace->SortedTwoBodyChannels_CC)
    {
      TwoBodyChannel_CC &tbc2 = Eta.modelspace->GetTwoBodyChannel_CC(ch2_cc);
      const int J2 = tbc2.J;
      if ((J0 + J2 < K) or (std::abs(J0 - J2) > K))
        continue;
      auto left_iter = barLeft.find({ch0_cc, ch2_cc});
      auto right_iter = barRight.find({ch2_cc, ch0_cc});
      if (left_iter == barLeft.end() or right_iter == barRight.end())
        continue;
      const arma::mat &left = left_iter->second;
      const arma::mat &right = right_iter->second;
      if (left.n_cols != right.n_rows)
        throw std::runtime_error("Full-square tensor Pandya product has incompatible dimensions.");
      double factor = scale * Phase(J0 + J2 + K) / (Hat(J0) * Hat(K));
      accum += factor * (left * right);
    }
  }
}

void AddToFullSquareMap(TensorPandyaFullMap &out, const ChannelPair &key, const arma::mat &term)
{
  arma::mat &accum = out[key];
  if (accum.n_rows == 0 and accum.n_cols == 0)
    accum.zeros(term.n_rows, term.n_cols);
  if (accum.n_rows != term.n_rows or accum.n_cols != term.n_cols)
    throw std::runtime_error("Full-square Pandya map accumulation has incompatible dimensions.");
  accum += term;
}

void BuildFullSquareScalarPandya(const Operator &Op, TensorPandyaFullMap &barOp)
{
  if (Op.GetJRank() != 0)
    throw std::invalid_argument("BuildFullSquareScalarPandya expects a scalar operator.");

  BuildFullSquareTensorPandya(Op, barOp);
  for (auto it = barOp.begin(); it != barOp.end();)
  {
    if (it->first[0] != it->first[1])
      it = barOp.erase(it);
    else
      ++it;
  }
}

void AddLeftScalarFullSquareProduct(const TensorPandyaFullMap &barScalar,
                                    const TensorPandyaFullMap &barTensor,
                                    TensorPandyaFullMap &barOut,
                                    double scale)
{
  for (const auto &tensor_it : barTensor)
  {
    const ChannelPair key = tensor_it.first;
    auto scalar_it = barScalar.find({key[0], key[0]});
    if (scalar_it == barScalar.end())
      continue;
    if (scalar_it->second.n_cols != tensor_it.second.n_rows)
      throw std::runtime_error("Left scalar full-square Pandya product has incompatible dimensions.");
    AddToFullSquareMap(barOut, key, scale * scalar_it->second * tensor_it.second);
  }
}

void AddRightScalarFullSquareProduct(const TensorPandyaFullMap &barTensor,
                                     const TensorPandyaFullMap &barScalar,
                                     TensorPandyaFullMap &barOut,
                                     double scale)
{
  for (const auto &tensor_it : barTensor)
  {
    const ChannelPair key = tensor_it.first;
    auto scalar_it = barScalar.find({key[1], key[1]});
    if (scalar_it == barScalar.end())
      continue;
    if (tensor_it.second.n_cols != scalar_it->second.n_rows)
      throw std::runtime_error("Right scalar full-square Pandya product has incompatible dimensions.");
    AddToFullSquareMap(barOut, key, scale * tensor_it.second * scalar_it->second);
  }
}

void AddTransposeTensorRightScalarFullSquareProduct(const TensorPandyaFullMap &barTensor,
                                                    const TensorPandyaFullMap &barScalar,
                                                    TensorPandyaFullMap &barOut,
                                                    double scale)
{
  for (const auto &tensor_it : barTensor)
  {
    const ChannelPair in_key = tensor_it.first;
    const ChannelPair out_key = {in_key[1], in_key[0]};
    auto scalar_it = barScalar.find({in_key[0], in_key[0]});
    if (scalar_it == barScalar.end())
      continue;
    if (tensor_it.second.n_rows != scalar_it->second.n_rows)
      throw std::runtime_error("Transposed tensor-right-scalar full-square Pandya product has incompatible dimensions.");
    AddToFullSquareMap(barOut, out_key, scale * tensor_it.second.t() * scalar_it->second);
  }
}

void AddInverseScalarPandyaFullDiagonalToOperator(const TensorPandyaFullMap &barOp, Operator &Z)
{
  std::vector<TwoBodyChannelPair> channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch = channel_pairs[ich][0];
    if (ch != channel_pairs[ich][1])
      continue;
    TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
    const int J0 = tbc.J;
    const size_t nKets = tbc.GetNumberKets();
    arma::mat &Zmat = Z.TwoBody.GetMatrix(ch, ch);

    for (size_t ibra = 0; ibra < nKets; ++ibra)
    {
      Ket &bra = tbc.GetKet(ibra);
      const index_t i = bra.p;
      const index_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ji = oi.j2;
      const int jj = oj.j2;

      for (size_t iket = ibra; iket < nKets; ++iket)
      {
        Ket &ket = tbc.GetKet(iket);
        const index_t k = ket.p;
        const index_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int jk = ok.j2;
        const int jl = ol.j2;
        double commij = 0.0;
        double commji = 0.0;

        int parity_cc = (oi.l + ol.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        int Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_kj < 0)
            continue;
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          commij -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_il, indx_kj);
        }

        parity_cc = (oj.l + ok.l) % 2;
        Tz_cc = std::abs(oj.tz2 - ok.tz2) / 2;
        Jpmin = std::max(std::abs(jj - jl), std::abs(jk - ji)) / 2;
        Jpmax = std::min(jj + jl, jk + ji) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
          int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
          if (indx_ik < 0 or indx_lj < 0)
            continue;
          indx_ik += (k > i ? nkets_cc : 0);
          indx_lj += (j > l ? nkets_cc : 0);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          commji -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_lj, indx_ik);
        }

        double zijkl = commij - Z.modelspace->phase((ji + jj) / 2 - J0) * commji;
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Zmat(ibra, iket) += zijkl;
        if (iket != ibra)
          Zmat(iket, ibra) += Z.modelspace->phase((ji + jj + jk + jl) / 2) * zijkl;
      }
    }
  }
}

void RecoupleScalarChiIII232(const TensorPandyaFullMap &barChiIII, TensorPandyaFullMap &barChiIII_RC, const Operator &Z)
{
  barChiIII_RC.clear();
  for (auto ch_cc : Z.modelspace->SortedTwoBodyChannels_CC)
  {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    const size_t nkets_cc = tbc_cc.GetNumberKets();
    if (nkets_cc == 0)
      continue;
    arma::mat &out = barChiIII_RC[{ch_cc, ch_cc}];
    out.zeros(2 * nkets_cc, 2 * nkets_cc);
    const int J_cc = tbc_cc.J;

    for (size_t ibra_cc = 0; ibra_cc < 2 * nkets_cc; ++ibra_cc)
    {
      Ket &bra_cc = tbc_cc.GetKet(ibra_cc % nkets_cc);
      const index_t a = ibra_cc < nkets_cc ? bra_cc.p : bra_cc.q;
      const index_t b = ibra_cc < nkets_cc ? bra_cc.q : bra_cc.p;
      if (ibra_cc >= nkets_cc and a == b)
        continue;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);

      for (size_t iket_cc = 0; iket_cc < 2 * nkets_cc; ++iket_cc)
      {
        Ket &ket_cc = tbc_cc.GetKet(iket_cc % nkets_cc);
        const index_t c = iket_cc < nkets_cc ? ket_cc.p : ket_cc.q;
        const index_t d = iket_cc < nkets_cc ? ket_cc.q : ket_cc.p;
        if (iket_cc >= nkets_cc and c == d)
          continue;
        Orbit &oc = Z.modelspace->GetOrbit(c);
        Orbit &od = Z.modelspace->GetOrbit(d);

        const int Jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(ob.j2 - oc.j2)) / 2;
        const int Jmax = std::min(oa.j2 + od.j2, ob.j2 + oc.j2) / 2;
        double xbar = 0.0;
        for (int J_std = Jmin; J_std <= Jmax; ++J_std)
        {
          const double sixj = Z.modelspace->GetSixJ(oa.j2 * 0.5, ob.j2 * 0.5, J_cc, oc.j2 * 0.5, od.j2 * 0.5, J_std);
          if (std::abs(sixj) < 1e-12)
            continue;
          const int parity_cc = (oa.l + od.l) % 2;
          const int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
          const index_t old_ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);
          auto old_iter = barChiIII.find({old_ch_cc, old_ch_cc});
          if (old_iter == barChiIII.end())
            continue;
          TwoBodyChannel_CC &old_tbc = Z.modelspace->GetTwoBodyChannel_CC(old_ch_cc);
          const size_t old_nkets = old_tbc.GetNumberKets();
          int indx_ad = old_tbc.GetLocalIndex(std::min(a, d), std::max(a, d));
          int indx_bc = old_tbc.GetLocalIndex(std::min(b, c), std::max(b, c));
          if (indx_ad < 0 or indx_bc < 0)
            continue;
          indx_ad += (a > d ? old_nkets : 0);
          indx_bc += (b > c ? old_nkets : 0);
          xbar -= Phase((ob.j2 + oc.j2) / 2 + J_std) * (2 * J_std + 1) * sixj *
                  (old_iter->second(indx_bc, indx_ad) + old_iter->second(indx_ad, indx_bc));
        }
        out(ibra_cc, iket_cc) = xbar;
      }
    }
  }
}

void AddInverseScalarPandyaFullPhasedToOperator(const TensorPandyaFullMap &barOp, Operator &Z, double scale)
{
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    if (nbras == 0 or nkets == 0)
      continue;
    const int J0 = tbc_bra.J;

    for (size_t ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t i = bra.p;
      const index_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ji = oi.j2;
      const int jj = oj.j2;
      const int phaseFactor = Phase(J0 + (ji + jj) / 2);

      const size_t ketmin = ch_bra == ch_ket ? ibra : 0;
      for (size_t iket = ketmin; iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t k = ket.p;
        const index_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int jk = ok.j2;
        const int jl = ol.j2;
        double commijkl = 0.0;
        double commjikl = 0.0;
        double commijlk = 0.0;
        double commjilk = 0.0;

        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_jl < 0 or indx_ik < 0)
            continue;
          indx_jl += (j > l ? nkets_cc : 0);
          indx_ik += (i > k ? nkets_cc : 0);
          commjikl -= Phase(Jprime + (ji + jk) / 2) * (2 * Jprime + 1) * sixj * bar_iter->second(indx_jl, indx_ik);
          commijlk -= Phase(Jprime + (jj + jl) / 2) * (2 * Jprime + 1) * sixj * bar_iter->second(indx_ik, indx_jl);
        }

        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_jk < 0)
            continue;
          indx_il += (i > l ? nkets_cc : 0);
          indx_jk += (j > k ? nkets_cc : 0);
          commjilk -= Phase(Jprime + (ji + jl) / 2) * (2 * Jprime + 1) * sixj * bar_iter->second(indx_jk, indx_il);
          commijkl -= Phase(Jprime + (jj + jk) / 2) * (2 * Jprime + 1) * sixj * bar_iter->second(indx_il, indx_jk);
        }

        double zijkl = commjikl - Phase((ji + jj) / 2 - J0) * commijkl;
        zijkl += -Phase((jl + jk) / 2 - J0) * commjilk + Phase((jk + jl + ji + jj) / 2) * commijlk;
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, scale * phaseFactor * zijkl);
      }
    }
  }
}

void BuildStandardOccWeightedTensorTwoBody(const Operator &Eta, int occ_type, Operator &Weighted)
{
  Weighted.Erase();
  for (auto &mat_it : Weighted.TwoBody.MatEl)
  {
    const size_t ch_bra = mat_it.first[0];
    const size_t ch_ket = mat_it.first[1];
    TwoBodyChannel &tbc_bra = Eta.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Eta.modelspace->GetTwoBodyChannel(ch_ket);
    arma::mat &weighted_mat = mat_it.second;
    const arma::mat &eta_mat = Eta.TwoBody.GetMatrix(ch_bra, ch_ket);
    for (size_t ibra = 0; ibra < tbc_bra.GetNumberKets(); ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      Orbit &oi = Eta.modelspace->GetOrbit(bra.p);
      Orbit &oj = Eta.modelspace->GetOrbit(bra.q);
      const double n_i = oi.occ;
      const double n_j = oj.occ;
      const double nbar_i = 1.0 - n_i;
      const double nbar_j = 1.0 - n_j;
      for (size_t iket = 0; iket < tbc_ket.GetNumberKets(); ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        Orbit &ok = Eta.modelspace->GetOrbit(ket.p);
        Orbit &ol = Eta.modelspace->GetOrbit(ket.q);
        const double n_k = ok.occ;
        const double n_l = ol.occ;
        const double occfactor = (occ_type == 0)
                                     ? (nbar_i * nbar_j * n_k + n_i * n_j * (1.0 - n_k))
                                     : (nbar_i * nbar_j * n_l + n_i * n_j * (1.0 - n_l));
        weighted_mat(ibra, iket) = occfactor * eta_mat(ibra, iket);
      }
    }
  }
}

void AddTransposeTwoBodyMatrices(const Operator &Source, Operator &Target, double scale)
{
  for (auto &target_it : Target.TwoBody.MatEl)
  {
    const size_t ch_bra = target_it.first[0];
    const size_t ch_ket = target_it.first[1];
    auto source_it = Source.TwoBody.MatEl.find({ch_ket, ch_bra});
    if (source_it == Source.TwoBody.MatEl.end())
      continue;
    target_it.second += scale * source_it->second.t();
  }
}

void AddInverseScalarPandyaFullUnphasedToOperator(const TensorPandyaFullMap &barOp, Operator &Z, double scale)
{
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    if (nbras == 0 or nkets == 0)
      continue;
    const int J0 = tbc_bra.J;

    for (size_t ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t i = bra.p;
      const index_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ji = oi.j2;
      const int jj = oj.j2;
      const size_t ketmin = ch_bra == ch_ket ? ibra : 0;

      for (size_t iket = ketmin; iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t k = ket.p;
        const index_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int jk = ok.j2;
        const int jl = ol.j2;
        double commijkl = 0.0;
        double commjikl = 0.0;
        double commijlk = 0.0;
        double commjilk = 0.0;

        int parity_cc = (oi.l + ol.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        int Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
          if (indx_il < 0 or indx_kj < 0)
            continue;
          int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
          int indx_li = indx_il + (l > i ? nkets_cc : 0);
          commjilk -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_jk, indx_li);
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          commijkl -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_il, indx_kj);
        }

        parity_cc = (oi.l + ok.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        Jpmin = std::max(std::abs(jj - jl), std::abs(jk - ji)) / 2;
        Jpmax = std::min(jj + jl, jk + ji) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_ki = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
          int indx_jl = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
          if (indx_ki < 0 or indx_jl < 0)
            continue;
          int indx_ik = indx_ki + (i > k ? nkets_cc : 0);
          int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
          commijlk -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_ik, indx_lj);
          indx_ki += (k > i ? nkets_cc : 0);
          indx_jl += (j > l ? nkets_cc : 0);
          commjikl -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_jl, indx_ki);
        }

        double zijkl = commijkl - Phase((ji + jj) / 2 - J0) * commjikl;
        zijkl += -Phase((jl + jk) / 2 - J0) * commijlk + Phase((jk + jl + ji + jj) / 2) * commjilk;
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, scale * zijkl);
      }
    }
  }
}

void AddInverseScalarPandyaFullVIIToOperator(const TensorPandyaFullMap &barOp, Operator &Z, double scale)
{
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    if (nbras == 0 or nkets == 0)
      continue;
    const int J0 = tbc_bra.J;

    for (size_t ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t i = bra.p;
      const index_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      const int ji = oi.j2;
      const int jj = oj.j2;
      const int phaseFactor = Phase(J0 + (ji + jj) / 2);
      const size_t ketmin = ch_bra == ch_ket ? ibra : 0;

      for (size_t iket = ketmin; iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t k = ket.p;
        const index_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        const int jk = ok.j2;
        const int jl = ol.j2;
        double commijkl = 0.0;
        double commjikl = 0.0;
        double commijlk = 0.0;
        double commjilk = 0.0;

        int parity_cc = (oi.l + ok.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
        int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
        int Jpmax = std::min(jj + jl, ji + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          double sixj2 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12 and std::abs(sixj2) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_lj = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
          int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
          if (indx_lj < 0 or indx_ik < 0)
            continue;
          int indx_jl = indx_lj + (j > l ? nkets_cc : 0);
          int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
          commjikl -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_jl, indx_ki);
          indx_ik += (i > k ? nkets_cc : 0);
          indx_lj += (l > j ? nkets_cc : 0);
          commijlk -= (2 * Jprime + 1) * sixj2 * bar_iter->second(indx_ik, indx_lj);
        }

        parity_cc = (oi.l + ol.l) % 2;
        Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
        Jpmax = std::min(ji + jl, jj + jk) / 2;
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
          double sixj2 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
          if (std::abs(sixj) < 1e-12 and std::abs(sixj2) < 1e-12)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          auto bar_iter = barOp.find({static_cast<index_t>(ch_cc), static_cast<index_t>(ch_cc)});
          if (bar_iter == barOp.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          const size_t nkets_cc = tbc_cc.GetNumberKets();
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
          if (indx_il < 0 or indx_kj < 0)
            continue;
          int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
          int indx_li = indx_il + (l > i ? nkets_cc : 0);
          commjilk -= (2 * Jprime + 1) * sixj2 * bar_iter->second(indx_jk, indx_li);
          indx_il += (i > l ? nkets_cc : 0);
          indx_kj += (k > j ? nkets_cc : 0);
          commijkl -= (2 * Jprime + 1) * sixj * bar_iter->second(indx_il, indx_kj);
        }

        double zijkl = commjikl - Phase((ji + jj) / 2 - J0) * commijkl;
        zijkl += -Phase((jl + jk) / 2 - J0) * commjilk + Phase((jk + jl + ji + jj) / 2) * commijlk;
        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, scale * phaseFactor * zijkl);
      }
    }
  }
}

void BuildTensorChiIII232(const Operator &Eta, Operator &ChiIII)
{
  TensorPandyaFullMap barEta;
  TensorPandyaFullMap nnnbarEta;
  TensorPandyaFullMap barChiIII;
  BuildFullSquareTensorPandya(Eta, barEta);
  BuildFullSquareTensorPandyaOccWeighted(Eta, 0, nnnbarEta);
  AddScalarCoupledTensorProductFullSquare(Eta, barEta, nnnbarEta, barChiIII, 1.0);
  AddInverseScalarPandyaFullDiagonalToOperator(barChiIII, ChiIII);
}

void AddTensorChiIII232DirectToZ(const Operator &Gamma, const Operator &ChiIII, Operator &Z)
{
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);
  const int hGamma = Gamma.IsHermitian() ? 1 : -1;

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    if (tbc_bra.GetNumberKets() == 0 or tbc_ket.GetNumberKets() == 0)
      continue;

    arma::mat update = ChiIII.TwoBody.GetMatrix(ch_bra, ch_bra) * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);
    update += hGamma * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket) * ChiIII.TwoBody.GetMatrix(ch_ket, ch_ket).t();
    Z.TwoBody.GetMatrix(ch_bra, ch_ket) += update;
  }
}

void AddInverseFullSquareTensorPandyaToOperator(const TensorPandyaFullMap &barOp, Operator &Op)
{
  const int K = Op.GetJRank();
  const int hOp = Op.IsHermitian() ? 1 : -1;
  auto channel_pairs = StoredTwoBodyChannelPairs(Op);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Op.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Op.modelspace->GetTwoBodyChannel(ch_ket);
    const int Jbra = tbc_bra.J;
    const int Jket = tbc_ket.J;
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    arma::mat &OpMat = Op.TwoBody.GetMatrix(ch_bra, ch_ket);

    for (size_t ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      const index_t i = bra.p;
      const index_t j = bra.q;
      Orbit &oi = Op.modelspace->GetOrbit(i);
      Orbit &oj = Op.modelspace->GetOrbit(j);
      const double ji = oi.j2 * 0.5;
      const double jj = oj.j2 * 0.5;

      const size_t ketmin = ch_bra == ch_ket ? ibra : 0;
      for (size_t iket = ketmin; iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        const index_t k = ket.p;
        const index_t l = ket.q;
        Orbit &ok = Op.modelspace->GetOrbit(k);
        Orbit &ol = Op.modelspace->GetOrbit(l);
        const double jk = ok.j2 * 0.5;
        const double jl = ol.j2 * 0.5;

        double commij = 0.0;
        double commji = 0.0;

        int parity_bra_cc = (oi.l + ol.l) % 2;
        int parity_ket_cc = (ok.l + oj.l) % 2;
        int Tz_bra_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Tz_ket_cc = std::abs(ok.tz2 - oj.tz2) / 2;
        int J3min = std::abs(int(ji - jl));
        int J3max = ji + jl;
        for (int J3 = J3min; J3 <= J3max; ++J3)
        {
          index_t ch_bra_cc = Op.modelspace->GetTwoBodyChannelIndex(J3, parity_bra_cc, Tz_bra_cc);
          TwoBodyChannel_CC &tbc_bra_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
          const size_t nbras_cc = tbc_bra_cc.GetNumberKets();
          int indx_il = tbc_bra_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
          if (indx_il < 0)
            continue;
          indx_il += (i > l ? nbras_cc : 0);

          int J4min = std::max(std::abs(int(jk - jj)), std::abs(J3 - K));
          int J4max = std::min(int(jk + jj), J3 + K);
          for (int J4 = J4min; J4 <= J4max; ++J4)
          {
            index_t ch_ket_cc = Op.modelspace->GetTwoBodyChannelIndex(J4, parity_ket_cc, Tz_ket_cc);
            TwoBodyChannel_CC &tbc_ket_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
            const size_t nkets_cc = tbc_ket_cc.GetNumberKets();
            int indx_kj = tbc_ket_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_kj < 0)
              continue;
            indx_kj += (k > j ? nkets_cc : 0);

            auto bar_iter = barOp.find({ch_bra_cc, ch_ket_cc});
            if (bar_iter == barOp.end())
              continue;

            double recoupling = 0.0;
            if (K == 0)
              recoupling = Phase(int(jj + jl) + Jbra + J3) * Op.modelspace->GetSixJ(ji, jj, Jbra, jk, jl, J3) / std::sqrt((2 * Jket + 1) * (2 * J4 + 1));
            else
              recoupling = Op.modelspace->GetNineJ(ji, jl, J3, jj, jk, J4, Jbra, Jket, K);
            if (std::abs(recoupling) < 1e-12)
              continue;

            double hatfactor = std::sqrt((2 * Jbra + 1) * (2 * Jket + 1) * (2 * J3 + 1) * (2 * J4 + 1));
            double tbme = bar_iter->second(indx_il, indx_kj);
            commij += hatfactor * Phase(int(jj + jl) + Jket + J4) * recoupling * tbme;
          }
        }

        if (i == j)
        {
          commji = commij;
        }
        else
        {
          parity_bra_cc = (oj.l + ol.l) % 2;
          parity_ket_cc = (ok.l + oi.l) % 2;
          Tz_bra_cc = std::abs(oj.tz2 - ol.tz2) / 2;
          Tz_ket_cc = std::abs(ok.tz2 - oi.tz2) / 2;
          J3min = std::abs(int(jj - jl));
          J3max = jj + jl;
          for (int J3 = J3min; J3 <= J3max; ++J3)
          {
            index_t ch_bra_cc = Op.modelspace->GetTwoBodyChannelIndex(J3, parity_bra_cc, Tz_bra_cc);
            TwoBodyChannel_CC &tbc_bra_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
            const size_t nbras_cc = tbc_bra_cc.GetNumberKets();
            int indx_jl = tbc_bra_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            if (indx_jl < 0)
              continue;
            indx_jl += (j > l ? nbras_cc : 0);

            int J4min = std::max(std::abs(int(jk - ji)), std::abs(J3 - K));
            int J4max = std::min(int(jk + ji), J3 + K);
            for (int J4 = J4min; J4 <= J4max; ++J4)
            {
              index_t ch_ket_cc = Op.modelspace->GetTwoBodyChannelIndex(J4, parity_ket_cc, Tz_ket_cc);
              TwoBodyChannel_CC &tbc_ket_cc = Op.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
              const size_t nkets_cc = tbc_ket_cc.GetNumberKets();
              int indx_ki = tbc_ket_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
              if (indx_ki < 0)
                continue;
              indx_ki += (k > i ? nkets_cc : 0);

              auto bar_iter = barOp.find({ch_bra_cc, ch_ket_cc});
              if (bar_iter == barOp.end())
                continue;

              double recoupling = 0.0;
              if (K == 0)
                recoupling = Phase(int(ji + jl) + Jbra + J3) * Op.modelspace->GetSixJ(jj, ji, Jbra, jk, jl, J3) / std::sqrt((2 * Jket + 1) * (2 * J4 + 1));
              else
                recoupling = Op.modelspace->GetNineJ(jj, jl, J3, ji, jk, J4, Jbra, Jket, K);
              if (std::abs(recoupling) < 1e-12)
                continue;

              double hatfactor = std::sqrt((2 * Jbra + 1) * (2 * Jket + 1) * (2 * J3 + 1) * (2 * J4 + 1));
              double tbme = bar_iter->second(indx_jl, indx_ki);
              commji += hatfactor * Phase(int(ji + jl) + Jket + J4) * recoupling * tbme;
            }
          }
        }

        double norm = bra.delta_pq() == ket.delta_pq() ? 1.0 + bra.delta_pq() : PhysConst::SQRT2;
        OpMat(ibra, iket) += (commij - Phase(int(ji + jj) - Jbra) * commji) / norm;
        if (ch_bra == ch_ket and iket != ibra)
          OpMat(iket, ibra) = hOp * OpMat(ibra, iket);
      }
    }
  }
}

void AddScalarCoupledTwoBodyProductToOperator(const Operator &Left,
                                              const Operator &Right,
                                              Operator &Z,
                                              double scale)
{
  const int K = Left.GetJRank();
  if (Right.GetJRank() != K or Z.GetJRank() != 0)
    throw std::invalid_argument("AddScalarCoupledTwoBodyProductToOperator expects two rank-K operators coupled to scalar Z.");

  auto channel_pairs = StoredTwoBodyChannelPairs(Z);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    const int Jbra = tbc_bra.J;
    const int Jket = tbc_ket.J;
    arma::mat update(Z.TwoBody.GetMatrix(ch_bra, ch_ket).n_rows,
                     Z.TwoBody.GetMatrix(ch_bra, ch_ket).n_cols,
                     arma::fill::zeros);

    for (const auto &mid_it : Left.TwoBody.MatEl)
    {
      if (mid_it.first[0] != ch_bra)
        continue;
      const size_t ch_mid = mid_it.first[1];
      auto right_it = Right.TwoBody.MatEl.find({ch_mid, ch_ket});
      if (right_it == Right.TwoBody.MatEl.end())
        continue;

      TwoBodyChannel &tbc_mid = Z.modelspace->GetTwoBodyChannel(ch_mid);
      const int Jmid = tbc_mid.J;
      if ((Jbra + Jmid < K) or (std::abs(Jbra - Jmid) > K) or
          (Jmid + Jket < K) or (std::abs(Jmid - Jket) > K))
        continue;

      const double factor = scale * ScalarCoupledTensorProductFactor(Jbra, Jmid, K);
      if (std::abs(factor) < 1e-12)
        continue;
      update += factor * mid_it.second * right_it->second;
    }
    Z.TwoBody.GetMatrix(ch_bra, ch_ket) += update;
  }
}

void BuildStandardDoubleOccWeightedTensorTwoBody231(const Operator &Eta, Operator &EtaOcc)
{
  EtaOcc.Erase();
  auto channel_pairs = StoredTwoBodyChannelPairs(EtaOcc);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Eta.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Eta.modelspace->GetTwoBodyChannel(ch_ket);
    const arma::mat &eta_mat = Eta.TwoBody.GetMatrix(ch_bra, ch_ket);
    arma::mat &occ_mat = EtaOcc.TwoBody.GetMatrix(ch_bra, ch_ket);

    for (size_t ibra = 0; ibra < tbc_bra.GetNumberKets(); ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      const double ni = bra.op->occ;
      const double nj = bra.oq->occ;
      for (size_t iket = 0; iket < tbc_ket.GetNumberKets(); ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        const double nk = ket.op->occ;
        const double nl = ket.oq->occ;
        const double occ = ni * nj * (1.0 - nk) * (1.0 - nl)
                         - (1.0 - ni) * (1.0 - nj) * nk * nl;
        occ_mat(ibra, iket) = occ * eta_mat(ibra, iket);
      }
    }
  }
}

void BuildScalarChi221a231(const Operator &Eta, Operator &ChiTwoBody, arma::mat &Chi221a)
{
  Operator EtaOcc = Eta;
  EtaOcc.Erase();
  ChiTwoBody.Erase();
  BuildStandardDoubleOccWeightedTensorTwoBody231(Eta, EtaOcc);
  AddScalarCoupledTwoBodyProductToOperator(EtaOcc, Eta, ChiTwoBody, 2.0);
  AddTransposeTwoBodyMatrices(ChiTwoBody, ChiTwoBody, 1.0);

  std::vector<index_t> allorb_vec(Eta.modelspace->all_orbits.begin(), Eta.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t id = 0; id < allorb_vec.size(); ++id)
  {
    index_t d = allorb_vec[id];
    Orbit &od = Eta.modelspace->GetOrbit(d);
    for (auto e : Eta.modelspace->all_orbits)
    {
      Orbit &oe = Eta.modelspace->GetOrbit(e);
      if (e > d or od.l != oe.l or od.j2 != oe.j2 or od.tz2 != oe.tz2)
        continue;
      double sum = 0.0;
      for (auto b : Eta.modelspace->all_orbits)
      {
        Orbit &ob = Eta.modelspace->GetOrbit(b);
        int Jmin = std::abs(od.j2 - ob.j2) / 2;
        int Jmax = (od.j2 + ob.j2) / 2;
        for (int J = Jmin; J <= Jmax; ++J)
          sum += ChiTwoBody.TwoBody.GetTBME_J(J, J, b, d, b, e);
      }
      Chi221a(d, e) += sum / (od.j2 + 1.0);
      if (d != e)
        Chi221a(e, d) += sum / (od.j2 + 1.0);
    }
  }
}

void AddDiagramI231ToOneBody(const arma::mat &Chi221a, const Operator &Gamma, Operator &Z, double scale)
{
  const int hZ = Z.IsHermitian() ? 1 : -1;
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < allorb_vec.size(); ++ip)
  {
    index_t p = allorb_vec[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      if (q > p)
        continue;
      Orbit &oq = Z.modelspace->GetOrbit(q);
      double zpq = 0.0;
      for (auto d : Z.modelspace->all_orbits)
      {
        Orbit &od = Z.modelspace->GetOrbit(d);
        for (auto e : Z.modelspace->all_orbits)
        {
          Orbit &oe = Z.modelspace->GetOrbit(e);
          if (std::abs(Chi221a(d, e)) < 1e-14)
            continue;
          int Jmin = std::abs(od.j2 - oq.j2) / 2;
          int Jmax = (od.j2 + oq.j2) / 2;
          for (int J = Jmin; J <= Jmax; ++J)
            zpq += (2 * J + 1) * Chi221a(d, e) * Gamma.TwoBody.GetTBME_J(J, J, e, p, d, q);
        }
      }
      Z.OneBody(p, q) += scale * zpq / (op.j2 + 1.0);
      if (p != q)
        Z.OneBody(q, p) += scale * hZ * zpq / (op.j2 + 1.0);
    }
  }
}

void BuildTensorChi221b231(const Operator &Eta, const Operator &Gamma, arma::mat &Chi221b)
{
  const int K = Eta.GetJRank();
  std::vector<index_t> allorb_vec(Eta.modelspace->all_orbits.begin(), Eta.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t id = 0; id < allorb_vec.size(); ++id)
  {
    index_t d = allorb_vec[id];
    Orbit &od = Eta.modelspace->GetOrbit(d);
    for (auto e : Eta.modelspace->all_orbits)
    {
      Orbit &oe = Eta.modelspace->GetOrbit(e);
      if ((std::abs(od.j2 - oe.j2) > 2 * K) or (od.j2 + oe.j2 < 2 * K))
        continue;
      double chi = 0.0;
      for (auto a : Eta.modelspace->all_orbits)
      {
        Orbit &oa = Eta.modelspace->GetOrbit(a);
        for (auto b : Eta.modelspace->all_orbits)
        {
          Orbit &ob = Eta.modelspace->GetOrbit(b);
          for (auto c : Eta.modelspace->all_orbits)
          {
            Orbit &oc = Eta.modelspace->GetOrbit(c);
            double occ = (1.0 - oa.occ) * (1.0 - oe.occ) * ob.occ * oc.occ
                       - (1.0 - ob.occ) * (1.0 - oc.occ) * oa.occ * oe.occ;
            if (std::abs(occ) < 1e-12)
              continue;
            int Jgamma_min = JMin(oa, od, ob, oc);
            int Jgamma_max = JMax(oa, od, ob, oc);
            for (int Jgamma = Jgamma_min; Jgamma <= Jgamma_max; ++Jgamma)
            {
              double gamma = Gamma.TwoBody.GetTBME_J(Jgamma, Jgamma, a, d, b, c);
              if (std::abs(gamma) < 1e-14)
                continue;
              int Jeta_min = std::max(std::abs(ob.j2 - oc.j2) / 2, std::abs(Jgamma - K));
              int Jeta_max = std::min((ob.j2 + oc.j2) / 2, Jgamma + K);
              for (int Jeta = Jeta_min; Jeta <= Jeta_max; ++Jeta)
              {
                double factor = ScalarCoupledTensorProductFactor(Jgamma, Jeta, K);
                if (std::abs(factor) < 1e-14)
                  continue;
                chi += occ * (2 * Jgamma + 1) * factor * Eta.TwoBody.GetTBME_J(Jeta, Jgamma, b, c, a, e) * gamma;
              }
            }
          }
        }
      }
      Chi221b(d, e) += chi / (od.j2 + 1.0);
    }
  }
}

void AddDiagramIII231ToOneBody(const arma::mat &Chi221b, const Operator &Eta, Operator &Z, double scale)
{
  const int hZ = Z.IsHermitian() ? 1 : -1;
  const int K = Eta.GetJRank();
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < allorb_vec.size(); ++ip)
  {
    index_t p = allorb_vec[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      if (q > p)
        continue;
      double zpq_a = 0.0;
      double zpq_b = 0.0;
      for (auto d : Z.modelspace->all_orbits)
      {
        Orbit &od = Z.modelspace->GetOrbit(d);
        for (auto e : Z.modelspace->all_orbits)
        {
          Orbit &oe = Z.modelspace->GetOrbit(e);
          if (std::abs(Chi221b(d, e)) < 1e-14 and std::abs(Chi221b(e, d)) < 1e-14)
            continue;
          int Jext_min = std::abs(oe.j2 - op.j2) / 2;
          int Jext_max = (oe.j2 + op.j2) / 2;
          for (int Jext = Jext_min; Jext <= Jext_max; ++Jext)
          {
            int Jint_min = std::max(std::abs(od.j2 - Z.modelspace->GetOrbit(q).j2) / 2, std::abs(Jext - K));
            int Jint_max = std::min((od.j2 + Z.modelspace->GetOrbit(q).j2) / 2, Jext + K);
            for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
            {
              double factor = ScalarCoupledTensorProductFactor(Jext, Jint, K);
              if (std::abs(factor) < 1e-14)
                continue;
              double eta = Eta.TwoBody.GetTBME_J(Jext, Jint, e, p, d, q);
              zpq_a += factor * Chi221b(d, e) * eta;
              zpq_b += factor * Chi221b(e, d) * eta;
            }
          }
        }
      }
      double update = scale * (zpq_a - hZ * zpq_b) / (op.j2 + 1.0);
      Z.OneBody(p, q) += update;
      if (p != q)
        Z.OneBody(q, p) += hZ * update;
    }
  }
}

void BuildScalarChi222bStandard231(const Operator &Eta, const Operator &Gamma, Operator &Chi222b)
{
  Operator EtaOcc = Eta;
  Operator EtaEtaOcc = Chi222b;
  Operator EtaOccEta = Chi222b;
  EtaOcc.Erase();
  EtaEtaOcc.Erase();
  EtaOccEta.Erase();
  Chi222b.Erase();
  BuildStandardDoubleOccWeightedTensorTwoBody231(Eta, EtaOcc);
  AddScalarCoupledTwoBodyProductToOperator(Eta, EtaOcc, EtaEtaOcc, 1.0);
  AddScalarCoupledTwoBodyProductToOperator(EtaOcc, Eta, EtaOccEta, 1.0);

  auto channel_pairs = StoredTwoBodyChannelPairs(Chi222b);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Chi222b.modelspace->GetTwoBodyChannel(ch_bra);
    int J = tbc_bra.J;
    arma::mat &out = Chi222b.TwoBody.GetMatrix(ch_bra, ch_ket);
    out += 4.0 * (2 * J + 1) * EtaEtaOcc.TwoBody.GetMatrix(ch_bra, ch_bra) * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);
    out -= 4.0 * (2 * J + 1) * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket) * EtaOccEta.TwoBody.GetMatrix(ch_ket, ch_ket);
  }
}

void AddStandardTrace231ToOneBody(const Operator &Chi222b, Operator &Z, double scale)
{
  const int hZ = Z.IsHermitian() ? 1 : -1;
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < allorb_vec.size(); ++ip)
  {
    index_t p = allorb_vec[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      if (q > p)
        continue;
      double zpq = 0.0;
      for (auto c : Z.modelspace->all_orbits)
      {
        Orbit &oc = Z.modelspace->GetOrbit(c);
        int Jmin = std::abs(oc.j2 - op.j2) / 2;
        int Jmax = (oc.j2 + op.j2) / 2;
        for (int J = Jmin; J <= Jmax; ++J)
          zpq += Chi222b.TwoBody.GetTBME_J(J, J, c, p, c, q);
      }
      Z.OneBody(p, q) += scale * zpq / (op.j2 + 1.0);
      if (p != q)
        Z.OneBody(q, p) += scale * hZ * zpq / (op.j2 + 1.0);
    }
  }
}

void BuildFullSquareTensorPandyaDoubleOccWeighted231(const Operator &Eta, TensorPandyaFullMap &barEtaOcc)
{
  BuildFullSquareTensorPandya(Eta, barEtaOcc);
  for (auto &itmat : barEtaOcc)
  {
    TwoBodyChannel_CC &tbc_bra = Eta.modelspace->GetTwoBodyChannel_CC(itmat.first[0]);
    TwoBodyChannel_CC &tbc_ket = Eta.modelspace->GetTwoBodyChannel_CC(itmat.first[1]);
    const size_t nbras = tbc_bra.GetNumberKets();
    const size_t nkets = tbc_ket.GetNumberKets();
    arma::mat &bar = itmat.second;
    for (size_t ibra = 0; ibra < 2 * nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra % nbras);
      const index_t a = ibra < nbras ? bra.p : bra.q;
      const index_t b = ibra < nbras ? bra.q : bra.p;
      Orbit &oa = Eta.modelspace->GetOrbit(a);
      Orbit &ob = Eta.modelspace->GetOrbit(b);
      for (size_t iket = 0; iket < 2 * nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket % nkets);
        const index_t c = iket < nkets ? ket.p : ket.q;
        const index_t d = iket < nkets ? ket.q : ket.p;
        Orbit &oc = Eta.modelspace->GetOrbit(c);
        Orbit &od = Eta.modelspace->GetOrbit(d);
        double occ = (1.0 - oc.occ) * (1.0 - ob.occ) * oa.occ * od.occ
                   - oc.occ * ob.occ * (1.0 - oa.occ) * (1.0 - od.occ);
        bar(ibra, iket) *= occ;
      }
    }
  }
}

void BuildScalarChi222aPandya231(const Operator &Eta, const Operator &Gamma, TensorPandyaFullMap &barChi222a)
{
  TensorPandyaFullMap barEta;
  TensorPandyaFullMap barEtaOcc;
  TensorPandyaFullMap barGamma;
  TensorPandyaFullMap barEtaEtaOcc;
  BuildFullSquareTensorPandya(Eta, barEta);
  BuildFullSquareTensorPandyaDoubleOccWeighted231(Eta, barEtaOcc);
  BuildFullSquareScalarPandya(Gamma, barGamma);
  AddScalarCoupledTensorProductFullSquare(Eta, barEta, barEtaOcc, barEtaEtaOcc, 1.0);
  AddRightScalarFullSquareProduct(barEtaEtaOcc, barGamma, barChi222a, 1.0);
}

void AddPandyaTrace231ToOneBody(const TensorPandyaFullMap &barChi222a, Operator &Z, double scale)
{
  const int hZ = Z.IsHermitian() ? 1 : -1;
  std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < allorb_vec.size(); ++ip)
  {
    index_t p = allorb_vec[ip];
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      if (q > p)
        continue;
      double zpq = 0.0;
      for (auto e : Z.modelspace->all_orbits)
      {
        Orbit &oe = Z.modelspace->GetOrbit(e);
        int Jmin = std::abs(op.j2 - oe.j2) / 2;
        int Jmax = (op.j2 + oe.j2) / 2;
        int parity_cc = (op.l + oe.l) % 2;
        int Tz_cc = std::abs(op.tz2 - oe.tz2) / 2;
        for (int J = Jmin; J <= Jmax; ++J)
        {
          index_t ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J, parity_cc, Tz_cc);
          auto iter = barChi222a.find({ch_cc, ch_cc});
          if (iter == barChi222a.end())
            continue;
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          size_t n = tbc_cc.GetNumberKets();
          int ind_pe = tbc_cc.GetLocalIndex(std::min(p, e), std::max(p, e));
          int ind_qe = tbc_cc.GetLocalIndex(std::min(q, e), std::max(q, e));
          int ind_eq = tbc_cc.GetLocalIndex(std::min(e, q), std::max(e, q));
          int ind_ep = tbc_cc.GetLocalIndex(std::min(e, p), std::max(e, p));
          if (ind_pe < 0 or ind_qe < 0 or ind_eq < 0 or ind_ep < 0)
            continue;
          ind_pe += (p > e ? n : 0);
          ind_qe += (q > e ? n : 0);
          ind_eq += (e > q ? n : 0);
          ind_ep += (e > p ? n : 0);
          zpq += iter->second(ind_pe, ind_qe);
          zpq -= iter->second(ind_eq, ind_ep);
        }
      }
      Z.OneBody(p, q) += scale * zpq / (op.j2 + 1.0);
      if (p != q)
        Z.OneBody(q, p) += scale * hZ * zpq / (op.j2 + 1.0);
    }
  }
}

void AddScalarCoupledTensorProduct(const Operator &Eta,
                                   const PandyaMap &barLeft,
                                   const PandyaMap &barRight,
                                   PandyaMap &barZ,
                                   double scale)
{
  const int K = Eta.GetJRank();

  for (auto ch0_cc : Eta.modelspace->SortedTwoBodyChannels_CC)
  {
    const auto &tbc0 = Eta.modelspace->GetTwoBodyChannel_CC(ch0_cc);
    const int J0 = tbc0.J;
    const double prefactor_J0 = scale * AngMom::phase(J0) / Hat(J0);

    arma::uvec z_rows_ph = arma::join_cols(tbc0.GetKetIndex_hh(), tbc0.GetKetIndex_ph());
    arma::mat accum(2 * z_rows_ph.n_rows, tbc0.GetNumberKets(), arma::fill::zeros);

    for (auto ch2_cc : Eta.modelspace->SortedTwoBodyChannels_CC)
    {
      const auto &tbc2 = Eta.modelspace->GetTwoBodyChannel_CC(ch2_cc);
      const int J2 = tbc2.J;
      if ((J0 + J2 < K) or (std::abs(J0 - J2) > K))
        continue;

      auto left_iter = barLeft.find({ch0_cc, ch2_cc});
      auto right_iter = barRight.find({ch2_cc, ch0_cc});
      if ((left_iter == barLeft.end()) or (right_iter == barRight.end()))
        continue;

      const arma::mat &left = left_iter->second;
      const arma::mat &right = right_iter->second;

      if (left.n_cols != right.n_rows)
      {
        throw std::runtime_error("Tensor Pandya product has incompatible dimensions; doubled ph row/column layout still needs the comm222_phst block embedding.");
      }

      double factor = prefactor_J0 * AngMom::phase(J2 + K) / Hat(K);
      accum += factor * (left * right);
    }

    if (arma::norm(accum, "fro") > 0)
      barZ[{ch0_cc, ch0_cc}] += accum;
  }
}

void AddScalarCoupledTensorCommutatorProduct(const Operator &Eta,
                                             const PandyaMap &barEta,
                                             const PandyaMap &barY,
                                             PandyaMap &barZ,
                                             double scale)
{
  AddScalarCoupledTensorProduct(Eta, barEta, barY, barZ, scale);
  AddScalarCoupledTensorProduct(Eta, barY, barEta, barZ, -scale);
}

void ThrowDiagramBodyMissing(const std::string &where)
{
  throw std::runtime_error(where + " needs a true tensor implementation: loop over stored (bra channel, ket channel) pairs and keep Jbra/Jket throughout the occupation-weighted intermediate. Do not call the scalar factorized kernel or collapse tensor Eta to GetTBME_J(J,J,...).");
}

} // namespace

void comm223_231st(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  (void)Eta;
  (void)Gamma;
  (void)Z;
  throw std::runtime_error(std::string(__func__) + " is disabled: tensor factorized public path is turned off while validating reference implementations.");
}

void comm223_231st_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  double t_start = omp_get_wtime();

  arma::mat Chi221a = Z.OneBody;
  arma::mat Chi221b = Eta.OneBody;
  Chi221a.zeros();
  Chi221b.zeros();

  Operator ChiTwoBody = Z;
  ChiTwoBody.Erase();
  BuildScalarChi221a231(Eta, ChiTwoBody, Chi221a);
  AddDiagramI231ToOneBody(Chi221a, Gamma, Z, 0.5);

  if (not use_goose_tank_only_1b)
  {
    BuildTensorChi221b231(Eta, Gamma, Chi221b);
    AddDiagramIII231ToOneBody(Chi221b, Eta, Z, 0.5);
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_231st_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  double t_start = omp_get_wtime();

  Operator Chi222b = Z;
  Chi222b.Erase();
  BuildScalarChi222bStandard231(Eta, Gamma, Chi222b);
  AddStandardTrace231ToOneBody(Chi222b, Z, 0.25);

  TensorPandyaFullMap barChi222a;
  BuildScalarChi222aPandya231(Eta, Gamma, barChi222a);
  AddPandyaTrace231ToOneBody(barChi222a, Z, 1.0);

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_232st(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  (void)Eta;
  (void)Gamma;
  (void)Z;
  throw std::runtime_error(std::string(__func__) + " is disabled: tensor factorized public path is turned off while validating reference implementations.");
}

void comm223_232st_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  double t_start = omp_get_wtime();
  const int K = Eta.GetJRank();
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

  arma::mat CHI_I = Z.OneBody;
  arma::mat CHI_II = Z.OneBody;
  CHI_I.zeros();
  CHI_II.zeros();

  std::vector<index_t> scalar_p_list;
  std::vector<index_t> scalar_q_list;
  for (auto p : Z.modelspace->all_orbits)
  {
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      scalar_p_list.push_back(p);
      scalar_q_list.push_back(q);
    }
  }

  std::vector<index_t> tensor_p_list;
  std::vector<index_t> tensor_q_list;
  for (auto p : Z.modelspace->all_orbits)
  {
    Orbit &op = Z.modelspace->GetOrbit(p);
    for (auto q : Eta.GetOneBodyChannel(op.l, op.j2, op.tz2))
    {
      tensor_p_list.push_back(p);
      tensor_q_list.push_back(q);
    }
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ipq = 0; ipq < scalar_p_list.size(); ++ipq)
  {
    index_t p = scalar_p_list[ipq];
    index_t q = scalar_q_list[ipq];
    Orbit &op = Z.modelspace->GetOrbit(p);
    Orbit &oq = Z.modelspace->GetOrbit(q);
    double chi_pq = 0.0;

    for (auto a : Z.modelspace->all_orbits)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      double nbar_a = 1.0 - oa.occ;
      if (std::abs(nbar_a) < 1e-12)
        continue;

      for (auto i : Z.modelspace->holes)
      {
        Orbit &oi = Z.modelspace->GetOrbit(i);
        double n_i = oi.occ;
        if (std::abs(n_i) < 1e-12)
          continue;

        for (auto j : Z.modelspace->holes)
        {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          double occfactor = nbar_a * n_i * oj.occ;
          if (std::abs(occfactor) < 1e-12)
            continue;

          int Jext_min = JMin(oa, op, oi, oj);
          int Jext_max = JMax(oa, op, oi, oj);
          for (int Jext = Jext_min; Jext <= Jext_max; ++Jext)
          {
            int Jint_min = std::max(std::abs(oi.j2 - oj.j2) / 2, std::abs(Jext - K));
            int Jint_max = std::min((oi.j2 + oj.j2) / 2, Jext + K);
            for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
            {
              double factor = ScalarCoupledTensorProductFactor(Jext, Jint, K);
              if (std::abs(factor) < 1e-12)
                continue;
              double eta_apij = Eta.TwoBody.GetTBME_J(Jext, Jint, a, p, i, j);
              double eta_ijaq = Eta.TwoBody.GetTBME_J(Jint, Jext, i, j, a, q);
              chi_pq += 0.5 * occfactor * factor * eta_apij * eta_ijaq / (oq.j2 + 1.0);
            }
          }
        }

        for (auto b : Z.modelspace->all_orbits)
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double occfactor = nbar_a * (1.0 - ob.occ) * n_i;
          if (std::abs(occfactor) < 1e-12)
            continue;

          int Jext_min = JMin(oi, op, oa, ob);
          int Jext_max = JMax(oi, op, oa, ob);
          for (int Jext = Jext_min; Jext <= Jext_max; ++Jext)
          {
            int Jint_min = std::max(std::abs(oa.j2 - ob.j2) / 2, std::abs(Jext - K));
            int Jint_max = std::min((oa.j2 + ob.j2) / 2, Jext + K);
            for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
            {
              double factor = ScalarCoupledTensorProductFactor(Jext, Jint, K);
              if (std::abs(factor) < 1e-12)
                continue;
              double eta_ipab = Eta.TwoBody.GetTBME_J(Jext, Jint, i, p, a, b);
              double eta_abiq = Eta.TwoBody.GetTBME_J(Jint, Jext, a, b, i, q);
              chi_pq += 0.5 * occfactor * factor * eta_ipab * eta_abiq / (oq.j2 + 1.0);
            }
          }
        }
      }
    }

    CHI_I(p, q) = chi_pq;
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ipq = 0; ipq < tensor_p_list.size(); ++ipq)
  {
    index_t p = tensor_p_list[ipq];
    index_t q = tensor_q_list[ipq];
    Orbit &op = Z.modelspace->GetOrbit(p);
    Orbit &oq = Z.modelspace->GetOrbit(q);
    double chi_pq = 0.0;

    for (auto a : Z.modelspace->all_orbits)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      double nbar_a = 1.0 - oa.occ;
      if (std::abs(nbar_a) < 1e-12)
        continue;

      for (auto i : Z.modelspace->holes)
      {
        Orbit &oi = Z.modelspace->GetOrbit(i);
        double n_i = oi.occ;
        if (std::abs(n_i) < 1e-12)
          continue;

        for (auto j : Z.modelspace->holes)
        {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          double occfactor = nbar_a * n_i * oj.occ;
          if (std::abs(occfactor) < 1e-12)
            continue;

          int Jgamma_min = JMin(oa, op, oi, oj);
          int Jgamma_max = JMax(oa, op, oi, oj);
          for (int Jgamma = Jgamma_min; Jgamma <= Jgamma_max; ++Jgamma)
          {
            double gamma_apij = Gamma.TwoBody.GetTBME_J(Jgamma, Jgamma, a, p, i, j);
            if (std::abs(gamma_apij) < 1e-12)
              continue;
            int Jeta_min = std::max(std::abs(oi.j2 - oj.j2) / 2, std::abs(Jgamma - K));
            int Jeta_max = std::min((oi.j2 + oj.j2) / 2, Jgamma + K);
            for (int Jeta = Jeta_min; Jeta <= Jeta_max; ++Jeta)
            {
              double factor = ScalarCoupledTensorProductFactor(Jgamma, Jeta, K);
              if (std::abs(factor) < 1e-12)
                continue;
              double eta_ijaq = Eta.TwoBody.GetTBME_J(Jeta, Jgamma, i, j, a, q);
              chi_pq += 0.5 * occfactor * factor * gamma_apij * eta_ijaq / (oq.j2 + 1.0);
            }
          }
        }

        for (auto b : Z.modelspace->all_orbits)
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double occfactor = nbar_a * (1.0 - ob.occ) * n_i;
          if (std::abs(occfactor) < 1e-12)
            continue;

          int Jgamma_min = JMin(oi, op, oa, ob);
          int Jgamma_max = JMax(oi, op, oa, ob);
          for (int Jgamma = Jgamma_min; Jgamma <= Jgamma_max; ++Jgamma)
          {
            double gamma_ipab = Gamma.TwoBody.GetTBME_J(Jgamma, Jgamma, i, p, a, b);
            if (std::abs(gamma_ipab) < 1e-12)
              continue;
            int Jeta_min = std::max(std::abs(oa.j2 - ob.j2) / 2, std::abs(Jgamma - K));
            int Jeta_max = std::min((oa.j2 + ob.j2) / 2, Jgamma + K);
            for (int Jeta = Jeta_min; Jeta <= Jeta_max; ++Jeta)
            {
              double factor = ScalarCoupledTensorProductFactor(Jgamma, Jeta, K);
              if (std::abs(factor) < 1e-12)
                continue;
              double eta_abiq = Eta.TwoBody.GetTBME_J(Jeta, Jgamma, a, b, i, q);
              chi_pq += 0.5 * occfactor * factor * gamma_ipab * eta_abiq / (oq.j2 + 1.0);
            }
          }
        }
      }
    }

    CHI_II(p, q) = chi_pq;
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    size_t nbras = tbc_bra.GetNumberKets();
    size_t nkets = tbc_ket.GetNumberKets();

    for (size_t ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      index_t p = bra.p;
      index_t q = bra.q;
      Orbit &op = Z.modelspace->GetOrbit(p);
      Orbit &oq = Z.modelspace->GetOrbit(q);
      size_t ketmin = (ch_bra == ch_ket) ? ibra : 0;

      for (size_t iket = ketmin; iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        index_t r = ket.p;
        index_t s = ket.q;
        Orbit &oR = Z.modelspace->GetOrbit(r);
        Orbit &os = Z.modelspace->GetOrbit(s);
        double zpqrs = 0.0;

        for (auto b : Z.GetOneBodyChannel(op.l, op.j2, op.tz2))
        {
          auto ibra_bq = tbc_bra.GetLocalIndex(std::min(b, q), std::max(b, q));
          if (ibra_bq < 0 or static_cast<size_t>(ibra_bq) >= nbras)
            continue;
          double norm = (b == q ? PhysConst::SQRT2 : 1.0) * (p == q ? 1.0 / PhysConst::SQRT2 : 1.0);
          if (b > q)
            norm *= bra.Phase(tbc_bra.J);
          zpqrs += norm * CHI_I(p, b) * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq, iket);
        }

        for (auto b : Eta.GetOneBodyChannel(op.l, op.j2, op.tz2))
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          int Jint_min = std::abs(ob.j2 - oq.j2) / 2;
          int Jint_max = (ob.j2 + oq.j2) / 2;
          for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
          {
            double factor = TensorOneBodyTensorTwoBodyToScalarFactor(Z, tbc_bra.J, Jint, K, op, ob, oq);
            if (std::abs(factor) < 1e-12)
              continue;
            zpqrs += factor * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(Jint, tbc_ket.J, b, q, r, s);
          }
        }

        for (auto b : Z.GetOneBodyChannel(oq.l, oq.j2, oq.tz2))
        {
          auto ibra_pb = tbc_bra.GetLocalIndex(std::min(p, b), std::max(p, b));
          if (ibra_pb < 0 or static_cast<size_t>(ibra_pb) >= nbras)
            continue;
          double norm = (b == p ? PhysConst::SQRT2 : 1.0) * (p == q ? 1.0 / PhysConst::SQRT2 : 1.0);
          if (p > b)
            norm *= bra.Phase(tbc_bra.J);
          zpqrs += norm * CHI_I(q, b) * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb, iket);
        }

        for (auto b : Eta.GetOneBodyChannel(oq.l, oq.j2, oq.tz2))
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          int Jint_min = std::abs(op.j2 - ob.j2) / 2;
          int Jint_max = (op.j2 + ob.j2) / 2;
          for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
          {
            double factor = TensorOneBodyTensorTwoBodyToScalarFactor(Z, tbc_bra.J, Jint, K, oq, ob, op);
            if (std::abs(factor) < 1e-12)
              continue;
            zpqrs += factor * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(Jint, tbc_ket.J, p, b, r, s);
          }
        }

        for (auto b : Z.GetOneBodyChannel(oR.l, oR.j2, oR.tz2))
        {
          auto iket_bs = tbc_ket.GetLocalIndex(std::min(b, s), std::max(b, s));
          if (iket_bs < 0 or static_cast<size_t>(iket_bs) >= nkets)
            continue;
          double norm = (b == s ? PhysConst::SQRT2 : 1.0) * (r == s ? 1.0 / PhysConst::SQRT2 : 1.0);
          if (b > s)
            norm *= ket.Phase(tbc_ket.J);
          zpqrs += norm * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) * CHI_I(b, r);
        }

        for (auto b : Eta.GetOneBodyChannel(oR.l, oR.j2, oR.tz2))
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          int Jint_min = std::abs(ob.j2 - os.j2) / 2;
          int Jint_max = (ob.j2 + os.j2) / 2;
          for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
          {
            double factor = TensorOneBodyTensorTwoBodyToScalarFactor(Z, tbc_ket.J, Jint, K, oR, ob, os);
            if (std::abs(factor) < 1e-12)
              continue;
            zpqrs -= factor * Eta.TwoBody.GetTBME_J(tbc_bra.J, Jint, p, q, b, s) * CHI_II(b, r);
          }
        }

        for (auto b : Z.GetOneBodyChannel(os.l, os.j2, os.tz2))
        {
          auto iket_rb = tbc_ket.GetLocalIndex(std::min(r, b), std::max(r, b));
          if (iket_rb < 0 or static_cast<size_t>(iket_rb) >= nkets)
            continue;
          double norm = (b == r ? PhysConst::SQRT2 : 1.0) * (r == s ? 1.0 / PhysConst::SQRT2 : 1.0);
          if (r > b)
            norm *= ket.Phase(tbc_ket.J);
          zpqrs += norm * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) * CHI_I(b, s);
        }

        for (auto b : Eta.GetOneBodyChannel(os.l, os.j2, os.tz2))
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          int Jint_min = std::abs(oR.j2 - ob.j2) / 2;
          int Jint_max = (oR.j2 + ob.j2) / 2;
          for (int Jint = Jint_min; Jint <= Jint_max; ++Jint)
          {
            double factor = TensorOneBodyTensorTwoBodyToScalarFactor(Z, tbc_ket.J, Jint, K, os, ob, oR);
            if (std::abs(factor) < 1e-12)
              continue;
            zpqrs -= factor * Eta.TwoBody.GetTBME_J(tbc_bra.J, Jint, p, q, r, b) * CHI_II(b, s);
          }
        }

        Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
      }
    }
  }

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_232st_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  double t_start = omp_get_wtime();
  const int hEta = Eta.IsHermitian() ? 1 : -1;

  Operator ChiIII = Z;
  ChiIII.Erase();
  BuildTensorChiIII232(Eta, ChiIII);
  AddTensorChiIII232DirectToZ(Gamma, ChiIII, Z);

  Operator ChiVI = Eta;
  Operator ChiVI_II = Eta;
  Operator EtaC = Eta;
  Operator EtaD = Eta;
  Operator ChiIV = Z;
  Operator ChiIV_D = Z;
  Operator ChiV = Eta;
  Operator ChiVII = Eta;
  ChiVI.Erase();
  ChiVI_II.Erase();
  EtaC.Erase();
  EtaD.Erase();
  ChiIV.Erase();
  ChiIV_D.Erase();
  ChiV.Erase();
  ChiVII.Erase();

  TensorPandyaFullMap barEta;
  TensorPandyaFullMap barGamma;
  TensorPandyaFullMap nnnbarEta;
  TensorPandyaFullMap nnnbarEtaD;
  TensorPandyaFullMap barChiIII;
  TensorPandyaFullMap barChiIII_RC;
  TensorPandyaFullMap chiIIIFinal;
  TensorPandyaFullMap barChiIV;
  TensorPandyaFullMap barChiGamma;
  TensorPandyaFullMap barChiV_RC;
  TensorPandyaFullMap chiVFinal;
  TensorPandyaFullMap barChiVII;
  TensorPandyaFullMap chiVIIFinal;
  TensorPandyaFullMap barChiV;
  TensorPandyaFullMap barChiVI;
  TensorPandyaFullMap barChiVI_II;
  BuildFullSquareTensorPandya(Eta, barEta);
  BuildFullSquareScalarPandya(Gamma, barGamma);
  BuildFullSquareTensorPandyaOccWeighted(Eta, 0, nnnbarEta);
  BuildFullSquareTensorPandyaOccWeighted(Eta, 1, nnnbarEtaD);
  AddScalarCoupledTensorProductFullSquare(Eta, barEta, nnnbarEta, barChiIII, 1.0);
  AddLeftScalarFullSquareProduct(barGamma, nnnbarEta, barChiV, 1.0);
  AddLeftScalarFullSquareProduct(barGamma, nnnbarEtaD, barChiVI, 1.0);
  AddTransposeTensorRightScalarFullSquareProduct(nnnbarEtaD, barGamma, barChiVI_II, Eta.IsHermitian() ? 1.0 : -1.0);
  AddInverseFullSquareTensorPandyaToOperator(barChiVI, ChiVI);
  AddInverseFullSquareTensorPandyaToOperator(barChiVI_II, ChiVI_II);
  AddScalarCoupledTwoBodyProductToOperator(Eta, ChiVI, Z, -1.0);
  AddScalarCoupledTwoBodyProductToOperator(ChiVI_II, Eta, Z, -1.0);

  RecoupleScalarChiIII232(barChiIII, barChiIII_RC, Z);
  AddLeftScalarFullSquareProduct(barGamma, barChiIII_RC, chiIIIFinal, 1.0);
  AddInverseScalarPandyaFullPhasedToOperator(chiIIIFinal, Z, 1.0);

  AddInverseFullSquareTensorPandyaToOperator(barChiV, ChiV);
  BuildFullSquareTensorPandya(ChiV, barChiV_RC);
  AddScalarCoupledTensorProductFullSquare(Eta, barEta, barChiV_RC, chiVFinal, 1.0);
  AddInverseScalarPandyaFullPhasedToOperator(chiVFinal, Z, 1.0);

  BuildStandardOccWeightedTensorTwoBody(Eta, 0, EtaC);
  BuildStandardOccWeightedTensorTwoBody(Eta, 1, EtaD);
  AddScalarCoupledTwoBodyProductToOperator(Eta, EtaC, ChiIV, 1.0);
  AddScalarCoupledTwoBodyProductToOperator(Eta, EtaD, ChiIV_D, 1.0);
  AddTransposeTwoBodyMatrices(ChiIV_D, ChiIV, 1.0);
  BuildFullSquareScalarPandya(ChiIV, barChiIV);
  AddRightScalarFullSquareProduct(barChiIV, barGamma, barChiGamma, 1.0);
  AddInverseScalarPandyaFullUnphasedToOperator(barChiGamma, Z, 0.5);

  auto channel_pairs = StoredTwoBodyChannelPairs(ChiVII);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    const size_t ch_bra = channel_pairs[ich][0];
    const size_t ch_ket = channel_pairs[ich][1];
    auto gamma_bra_it = Gamma.TwoBody.MatEl.find({ch_bra, ch_bra});
    auto gamma_ket_it = Gamma.TwoBody.MatEl.find({ch_ket, ch_ket});
    arma::mat &chi_mat = ChiVII.TwoBody.GetMatrix(ch_bra, ch_ket);
    if (gamma_bra_it != Gamma.TwoBody.MatEl.end())
      chi_mat += gamma_bra_it->second * EtaD.TwoBody.GetMatrix(ch_bra, ch_ket);
    auto eta_reverse_it = EtaD.TwoBody.MatEl.find({ch_ket, ch_bra});
    if (eta_reverse_it != EtaD.TwoBody.MatEl.end() and gamma_ket_it != Gamma.TwoBody.MatEl.end())
      chi_mat += hEta * eta_reverse_it->second.t() * gamma_ket_it->second;
  }
  BuildFullSquareTensorPandya(ChiVII, barChiVII);
  AddScalarCoupledTensorProductFullSquare(Eta, barChiVII, barEta, chiVIIFinal, 1.0);
  AddInverseScalarPandyaFullVIIToOperator(chiVIIFinal, Z, 0.5);

  Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void comm223_132st(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  (void)Eta;
  (void)Gamma;
  (void)Z;
  throw std::runtime_error(std::string(__func__) + " is disabled: tensor factorized public path is turned off while validating reference implementations.");
}

void comm223_132st_ladder(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  const int K = Eta.GetJRank();
  auto &Z2 = Z.TwoBody;
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    int J0 = tbc_bra.J;
    int J1 = tbc_ket.J;
    if (J0 != J1)
      continue;

    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
    for (int ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      for (int iket = (ch_bra == ch_ket ? ibra : 0); iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        double zijkl = 0.0;

        for (auto a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          for (auto b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double occfactor = (1.0 - n_a) * ob.occ - n_a * (1.0 - ob.occ);
            if (std::abs(occfactor) < 1e-12)
              continue;

            double eta_ba = Eta.OneBody(b, a);
            if (std::abs(eta_ba) < 1e-12)
              continue;

            for (auto c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              int J2min = JMin(oc, oa, Z.modelspace->GetOrbit(k), Z.modelspace->GetOrbit(l));
              int J2max = JMax(oc, oa, Z.modelspace->GetOrbit(k), Z.modelspace->GetOrbit(l));
              for (int J2 = J2min; J2 <= J2max; ++J2)
              {
                double eta_cakl = Eta.TwoBody.GetTBME_J(J2, J0, c, a, k, l);
                double gamma_ijcb = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, c, b);
                if (std::abs(eta_cakl * gamma_ijcb) > 1e-24)
                {
                  double sixj = Z.modelspace->GetSixJ(J0, J2, K, oa.j2 * 0.5, ob.j2 * 0.5, oc.j2 * 0.5);
                  double factor = occfactor * Phase(J2 + HalfIntegerPhaseExponent({oa.j2, oc.j2})) * Hat(J2) / Hat(K) * sixj;
                  zijkl += factor * eta_ba * eta_cakl * gamma_ijcb;
                }

                double gamma_cakl = Gamma.TwoBody.GetTBME_J(J0, J0, c, a, k, l);
                double eta_ijcb = Eta.TwoBody.GetTBME_J(J0, J2, i, j, c, b);
                if (std::abs(gamma_cakl * eta_ijcb) > 1e-24)
                {
                  double sixj = Z.modelspace->GetSixJ(J2, K, J0, oa.j2 * 0.5, oc.j2 * 0.5, ob.j2 * 0.5);
                  double factor = -Phase(J0 + HalfIntegerPhaseExponent({oa.j2, oc.j2})) * occfactor * Hat(J2) / Hat(K) * sixj;
                  zijkl += factor * eta_ba * gamma_cakl * eta_ijcb;
                }
              }
            }
          }
        }

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

void comm223_132st_onebody(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  const int K = Eta.GetJRank();
  auto &Z2 = Z.TwoBody;
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    int J0 = tbc_bra.J;
    int J1 = tbc_ket.J;
    if (J0 != J1)
      continue;

    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
    for (int ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      for (int iket = (ch_bra == ch_ket ? ibra : 0); iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        double zijkl = 0.0;

        for (auto a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          for (auto b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double occfactor = (1.0 - n_a) * ob.occ - n_a * (1.0 - ob.occ);
            double eta_ba = Eta.OneBody(b, a);
            if (std::abs(occfactor * eta_ba) < 1e-12)
              continue;

            for (auto c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              if (oc.j2 == oi.j2)
              {
                double gamma_cjkl = Gamma.TwoBody.GetTBME_J(J0, J0, c, j, k, l);
                for (int J2 = JMin(oi, oa, oc, ob); J2 <= JMax(oi, oa, oc, ob); ++J2)
                {
                  for (int J3 = JMin(oi, oa, oc, ob); J3 <= JMax(oi, oa, oc, ob); ++J3)
                  {
                    double eta_iacb = Eta.TwoBody.GetTBME_J(J2, J3, i, a, c, b);
                    double sixj = Z.modelspace->GetSixJ(J3, K, J2, oa.j2 * 0.5, oi.j2 * 0.5, ob.j2 * 0.5);
                    double factor = Phase(J2 + HalfIntegerPhaseExponent({oi.j2, oa.j2})) * Hat(J0) * Hat(J2) * Hat(J3) / ((oi.j2 + 1.0) * Hat(K)) * sixj;
                    zijkl += factor * occfactor * eta_ba * eta_iacb * gamma_cjkl;
                  }
                }
              }
              if (oc.j2 == oj.j2)
              {
                double gamma_iclk = Gamma.TwoBody.GetTBME_J(J0, J0, i, c, k, l);
                for (int J2 = JMin(oj, oa, oc, ob); J2 <= JMax(oj, oa, oc, ob); ++J2)
                {
                  for (int J3 = JMin(oj, oa, oc, ob); J3 <= JMax(oj, oa, oc, ob); ++J3)
                  {
                    double eta_jacb = Eta.TwoBody.GetTBME_J(J2, J3, j, a, c, b);
                    double sixj = Z.modelspace->GetSixJ(J3, K, J2, oa.j2 * 0.5, oj.j2 * 0.5, ob.j2 * 0.5);
                    double factor = Phase(J2 + HalfIntegerPhaseExponent({oj.j2, oa.j2})) * Hat(J0) * Hat(J2) * Hat(J3) / ((oj.j2 + 1.0) * Hat(K)) * sixj;
                    zijkl += factor * occfactor * eta_ba * eta_jacb * gamma_iclk;
                  }
                }
              }
              if (oc.j2 == ok.j2)
              {
                double gamma_ijcl = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, c, l);
                for (int J2 = JMin(oc, oa, ok, ob); J2 <= JMax(oc, oa, ok, ob); ++J2)
                {
                  for (int J3 = JMin(oc, oa, ok, ob); J3 <= JMax(oc, oa, ok, ob); ++J3)
                  {
                    double eta_cakb = Eta.TwoBody.GetTBME_J(J2, J3, c, a, k, b);
                    double sixj = Z.modelspace->GetSixJ(J3, K, J2, oa.j2 * 0.5, ok.j2 * 0.5, ob.j2 * 0.5);
                    double factor = -Phase(J2 + HalfIntegerPhaseExponent({ok.j2, oa.j2})) * Hat(J0) * Hat(J2) * Hat(J3) / ((ok.j2 + 1.0) * Hat(K)) * sixj;
                    zijkl += factor * occfactor * eta_ba * gamma_ijcl * eta_cakb;
                  }
                }
              }
              if (oc.j2 == ol.j2)
              {
                double gamma_ijkc = Gamma.TwoBody.GetTBME_J(J0, J0, i, j, k, c);
                for (int J2 = JMin(oc, oa, ol, ob); J2 <= JMax(oc, oa, ol, ob); ++J2)
                {
                  for (int J3 = JMin(oc, oa, ol, ob); J3 <= JMax(oc, oa, ol, ob); ++J3)
                  {
                    double eta_calb = Eta.TwoBody.GetTBME_J(J2, J3, c, a, l, b);
                    double sixj = Z.modelspace->GetSixJ(J3, K, J2, oa.j2 * 0.5, ol.j2 * 0.5, ob.j2 * 0.5);
                    double factor = -Phase(J2 + HalfIntegerPhaseExponent({ol.j2, oa.j2})) * Hat(J0) * Hat(J2) * Hat(J3) / ((ol.j2 + 1.0) * Hat(K)) * sixj;
                    zijkl += factor * occfactor * eta_ba * gamma_ijkc * eta_calb;
                  }
                }

                for (int J2 = JMin(oi, oa, oc, ob); J2 <= JMax(oi, oa, oc, ob); ++J2)
                {
                  double gamma_iacb = Gamma.TwoBody.GetTBME_J(J2, J2, i, a, c, b);
                  double sixj_left = Z.modelspace->GetSixJ(oc.j2 * 0.5, oi.j2 * 0.5, K, oa.j2 * 0.5, ob.j2 * 0.5, J2);
                  if (std::abs(gamma_iacb * sixj_left) < 1e-24)
                    continue;
                  for (int J3 = JMin(oc, oj, Z.modelspace->GetOrbit(k), ol); J3 <= JMax(oc, oj, Z.modelspace->GetOrbit(k), ol); ++J3)
                  {
                    double eta_cjkl = Eta.TwoBody.GetTBME_J(J3, J0, c, j, k, l);
                    double sixj_right = Z.modelspace->GetSixJ(J0, J3, K, oc.j2 * 0.5, oi.j2 * 0.5, oj.j2 * 0.5);
                    double factor = -Phase(J0 + J2 + HalfIntegerPhaseExponent({oi.j2, oj.j2, ob.j2, oc.j2})) * Hat(J2) * Hat(J2) * Hat(J3) / Hat(K) * sixj_left * sixj_right;
                    zijkl += factor * occfactor * eta_ba * gamma_iacb * eta_cjkl;
                  }
                }

                for (int J2 = JMin(oj, oa, oc, ob); J2 <= JMax(oj, oa, oc, ob); ++J2)
                {
                  double gamma_jacb = Gamma.TwoBody.GetTBME_J(J2, J2, j, a, c, b);
                  double sixj_left = Z.modelspace->GetSixJ(oc.j2 * 0.5, oj.j2 * 0.5, K, oa.j2 * 0.5, ob.j2 * 0.5, J2);
                  if (std::abs(gamma_jacb * sixj_left) < 1e-24)
                    continue;
                  for (int J3 = JMin(oi, oc, Z.modelspace->GetOrbit(k), ol); J3 <= JMax(oi, oc, Z.modelspace->GetOrbit(k), ol); ++J3)
                  {
                    double eta_ickl = Eta.TwoBody.GetTBME_J(J3, J0, i, c, k, l);
                    double sixj_right = Z.modelspace->GetSixJ(J0, J3, K, oc.j2 * 0.5, oj.j2 * 0.5, oi.j2 * 0.5);
                    double factor = Phase(J2 + J3 + HalfIntegerPhaseExponent({oi.j2, ob.j2})) * Hat(J2) * Hat(J2) * Hat(J3) / Hat(K) * sixj_left * sixj_right;
                    zijkl += factor * occfactor * eta_ba * gamma_jacb * eta_ickl;
                  }
                }
              }

              if (oc.j2 == ol.j2)
              {
                for (int J2 = JMin(oi, oj, oc, ol); J2 <= JMax(oi, oj, oc, ol); ++J2)
                {
                  double eta_ijcl = Eta.TwoBody.GetTBME_J(J0, J2, i, j, c, l);
                  if (std::abs(eta_ijcl) < 1e-12)
                    continue;
                  for (int J3 = JMin(oc, oa, ok, ob); J3 <= JMax(oc, oa, ok, ob); ++J3)
                  {
                    double gamma_cakb = Gamma.TwoBody.GetTBME_J(J3, J3, c, a, k, b);
                    double sixj_left = Z.modelspace->GetSixJ(ok.j2 * 0.5, oc.j2 * 0.5, K, oa.j2 * 0.5, ob.j2 * 0.5, J3);
                    double sixj_right = Z.modelspace->GetSixJ(J2, K, J0, ok.j2 * 0.5, ol.j2 * 0.5, oc.j2 * 0.5);
                    double factor = Phase(J2 + J3 + HalfIntegerPhaseExponent({ok.j2, ol.j2, ob.j2, oc.j2})) * Hat(J2) * Hat(J3) * Hat(J3) / Hat(K) * sixj_left * sixj_right;
                    zijkl += factor * occfactor * eta_ba * eta_ijcl * gamma_cakb;
                  }
                }

                for (int J2 = JMin(oi, oj, ok, oc); J2 <= JMax(oi, oj, ok, oc); ++J2)
                {
                  double eta_ijkc = Eta.TwoBody.GetTBME_J(J0, J2, i, j, k, c);
                  if (std::abs(eta_ijkc) < 1e-12)
                    continue;
                  for (int J3 = JMin(oc, oa, ol, ob); J3 <= JMax(oc, oa, ol, ob); ++J3)
                  {
                    double gamma_calb = Gamma.TwoBody.GetTBME_J(J3, J3, c, a, l, b);
                    double sixj_left = Z.modelspace->GetSixJ(ol.j2 * 0.5, oc.j2 * 0.5, K, oa.j2 * 0.5, ob.j2 * 0.5, J3);
                    double sixj_right = Z.modelspace->GetSixJ(J2, K, J0, ol.j2 * 0.5, ok.j2 * 0.5, oc.j2 * 0.5);
                    double factor = -Phase(J0 + J3 + HalfIntegerPhaseExponent({ok.j2, ob.j2})) * Hat(J2) * Hat(J3) * Hat(J3) / Hat(K) * sixj_left * sixj_right;
                    zijkl += factor * occfactor * eta_ba * eta_ijkc * gamma_calb;
                  }
                }
              }
            }
          }
        }

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

void comm223_132st_cross(const Operator &Eta, const Operator &Gamma, Operator &Z)
{
  CheckTensorEtaInputs(Eta, Gamma, Z, __func__);
  const int K = Eta.GetJRank();
  auto &Z2 = Z.TwoBody;
  auto channel_pairs = StoredTwoBodyChannelPairs(Z);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < channel_pairs.size(); ++ich)
  {
    size_t ch_bra = channel_pairs[ich][0];
    size_t ch_ket = channel_pairs[ich][1];
    TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
    int J0 = tbc_bra.J;
    int J1 = tbc_ket.J;
    if (J0 != J1)
      continue;

    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
    for (int ibra = 0; ibra < nbras; ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      for (int iket = (ch_bra == ch_ket ? ibra : 0); iket < nkets; ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        double zijkl = 0.0;

        auto add_cross = [&](int sign, int phase_base, size_t sixj_a, size_t sixj_b, size_t sixj_c,
                             size_t ninej_1, size_t ninej_2, int ninej_3_is_J2, size_t ninej_4, int ninej_5_is_J4, size_t ninej_6,
                             size_t ninej_8, size_t ninej_9,
                             size_t eta1, size_t eta2, size_t eta3, size_t eta4,
                             size_t gamma1, size_t gamma2, size_t gamma3, size_t gamma4,
                             double eta_ab, double occfactor) {
          Orbit &osixj_a = Z.modelspace->GetOrbit(sixj_a);
          Orbit &osixj_b = Z.modelspace->GetOrbit(sixj_b);
          Orbit &osixj_c = Z.modelspace->GetOrbit(sixj_c);
          Orbit &on1 = Z.modelspace->GetOrbit(ninej_1);
          Orbit &on2 = Z.modelspace->GetOrbit(ninej_2);
          Orbit &on4 = Z.modelspace->GetOrbit(ninej_4);
          Orbit &on6 = Z.modelspace->GetOrbit(ninej_6);
          Orbit &on8 = Z.modelspace->GetOrbit(ninej_8);
          Orbit &on9 = Z.modelspace->GetOrbit(ninej_9);
          Orbit &oe1 = Z.modelspace->GetOrbit(eta1);
          Orbit &oe2 = Z.modelspace->GetOrbit(eta2);
          Orbit &oe3 = Z.modelspace->GetOrbit(eta3);
          Orbit &oe4 = Z.modelspace->GetOrbit(eta4);
          Orbit &og1 = Z.modelspace->GetOrbit(gamma1);
          Orbit &og2 = Z.modelspace->GetOrbit(gamma2);
          Orbit &og3 = Z.modelspace->GetOrbit(gamma3);
          Orbit &og4 = Z.modelspace->GetOrbit(gamma4);

          for (int J2 = JMin(oe1, oe2, oe3, oe4); J2 <= JMax(oe1, oe2, oe3, oe4); ++J2)
          {
            for (int J3 = JMin(oe1, oe2, oe3, oe4); J3 <= JMax(oe1, oe2, oe3, oe4); ++J3)
            {
              double eta_me = Eta.TwoBody.GetTBME_J(J2, J3, eta1, eta2, eta3, eta4);
              if (std::abs(eta_me) < 1e-12)
                continue;
              double sixj = ninej_3_is_J2 ? Z.modelspace->GetSixJ(J2, J3, K, osixj_a.j2 * 0.5, osixj_b.j2 * 0.5, osixj_c.j2 * 0.5)
                                          : Z.modelspace->GetSixJ(J3, J2, K, osixj_b.j2 * 0.5, osixj_a.j2 * 0.5, osixj_c.j2 * 0.5);
              if (std::abs(sixj) < 1e-12)
                continue;
              for (int J4 = JMin(og1, og2, og3, og4); J4 <= JMax(og1, og2, og3, og4); ++J4)
              {
                double gamma_me = Gamma.TwoBody.GetTBME_J(J4, J4, gamma1, gamma2, gamma3, gamma4);
                if (std::abs(gamma_me) < 1e-12)
                  continue;
                double ninej = Z.modelspace->GetNineJ(on1.j2 * 0.5, on2.j2 * 0.5, ninej_3_is_J2 ? J2 : J3,
                                                      on4.j2 * 0.5, J4, on6.j2 * 0.5,
                                                      J0, on8.j2 * 0.5, on9.j2 * 0.5);
                if (std::abs(ninej) < 1e-12)
                  continue;
                int phase_exp = phase_base + (ninej_5_is_J4 ? J4 : J2 + J3 + J4);
                double factor = sign * Phase(phase_exp) * Hat(J0) * Hat(J2) * Hat(J3) * Hat(J4) * Hat(J4) / Hat(K);
                zijkl += factor * occfactor * eta_ab * sixj * ninej * eta_me * gamma_me;
              }
            }
          }
        };

        for (auto a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          for (auto b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double occfactor = (1.0 - n_a) * ob.occ - n_a * (1.0 - ob.occ);
            double eta_ab = Eta.OneBody(a, b);
            if (std::abs(occfactor * eta_ab) < 1e-12)
              continue;
            for (auto c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              add_cross(-1, J0 + HalfIntegerPhaseExponent({ok.j2, oc.j2}), a, b, k, i, c, 1, j, 1, b, l, k, i, c, k, a, b, j, c, l, eta_ab, occfactor);
              add_cross(-1, J0 + HalfIntegerPhaseExponent({oi.j2, oa.j2, ob.j2, oc.j2}), b, a, i, i, a, 0, j, 1, c, l, k, i, b, k, c, c, j, a, l, eta_ab, occfactor);
              add_cross(-1, HalfIntegerPhaseExponent({ok.j2, oc.j2}), a, b, l, i, c, 1, j, 1, b, k, l, i, c, l, a, b, j, c, k, eta_ab, occfactor);
              add_cross(+1, HalfIntegerPhaseExponent({oi.j2, ok.j2, ol.j2, oa.j2, ob.j2, oc.j2}), b, a, i, i, a, 0, j, 1, c, k, l, i, b, l, c, c, j, a, k, eta_ab, occfactor);
              add_cross(+1, HalfIntegerPhaseExponent({oi.j2, oj.j2, ok.j2, oc.j2}), a, b, k, j, c, 1, i, 1, b, l, k, j, c, k, a, b, i, c, l, eta_ab, occfactor);
              add_cross(-1, HalfIntegerPhaseExponent({oi.j2, oa.j2, ob.j2, oc.j2}), b, a, j, j, a, 0, i, 1, c, l, k, j, b, k, c, c, i, a, l, eta_ab, occfactor);
              add_cross(+1, J0 + HalfIntegerPhaseExponent({oi.j2, oj.j2, ok.j2, oc.j2}), a, b, l, j, c, 1, i, 1, b, k, l, j, c, l, a, b, i, c, k, eta_ab, occfactor);
              add_cross(+1, J0 + HalfIntegerPhaseExponent({oi.j2, ok.j2, ol.j2, oa.j2, ob.j2, oc.j2}), b, a, j, j, a, 0, i, 1, c, k, l, j, b, l, c, c, i, a, k, eta_ab, occfactor);
            }
          }
        }

        if (i == j)
          zijkl /= PhysConst::SQRT2;
        if (k == l)
          zijkl /= PhysConst::SQRT2;
        Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
      }
    }
  }
}

} // namespace FactorizedDoubleCommutator

} // namespace Commutator