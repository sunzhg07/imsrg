#include "EOM.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"
#include "Commutator.hh"
#include "IMSRG3Commutators.hh"
#include "FactorizedDoubleCommutator.hh"
#include "FactorizedDoubleCommutator_eths.hh"
#include "ReadWrite.hh"
#include "UnitTest.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using PhysConst::SQRT2;

// we have two constructor for the EOM, w/o the rdm for multi-reference and
// single reference
EOM::EOM(Operator &Hs, Operator &rdm, int J2,  int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), rdm(rdm), J2(J2), parity(parity),
      itz(itz), is_multiref(true) {
  eom_dims = 0;
  qv_dim = 0;
  ph_dim = 0;
  ppvv_dim = 0;
  pphv_dim = 0;
  pphh_dim = 0;
  rdm_ms = this->rdm.modelspace;
  BuildOrbMap();
};

EOM::EOM(Operator &Hs, const std::string &tdm_file, int J2,  int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), J2(J2), parity(parity),
      itz(itz), is_multiref(true) {
  eom_dims = 0;
  qv_dim = 0;
  ph_dim = 0;
  ppvv_dim = 0;
  pphv_dim = 0;
  pphh_dim = 0;
  rdm = ReadTdm(tdm_file);
  BuildOrbMap();
};

EOM::EOM(Operator &Hs, int J2,  int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), J2(J2), parity(parity),
      itz(itz), is_multiref(false) {
  rdm_ms = modelspace;
  BuildOrbMap();
};

EOM::EOM(Operator &Hs, int J2, int parity, int itz, SREOMMode mode)
    : modelspace(Hs.modelspace), Hs(Hs), J2(J2), parity(parity),
      itz(itz), is_multiref(false), sr_mode(mode) {
  rdm_ms = modelspace;
  BuildOrbMap();
};

void EOM::BuildOrbMap() {
  size_t norbits = modelspace->norbits;
  hs_to_rdm_orb.assign(norbits, -1);
  for (size_t i = 0; i < norbits; ++i) {
    Orbit &oi = modelspace->GetOrbit(i);
    size_t idx = rdm_ms->GetOrbitIndex(oi.n, oi.l, oi.j2, oi.tz2);
    if ((int)idx != ModelSpace::NOT_AN_ORBIT)
      hs_to_rdm_orb[i] = (int)idx;
  }
}

double EOM::RdmOB(size_t i_hs, size_t j_hs) const {
  // Only valence orbits in Hs have a corresponding rdm entry
  if (modelspace->GetOrbit(i_hs).cvq != 1) return 0.0;
  if (modelspace->GetOrbit(j_hs).cvq != 1) return 0.0;
  int i_rdm = (i_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[i_hs] : -1;
  int j_rdm = (j_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[j_hs] : -1;
  if (i_rdm < 0 || j_rdm < 0) return 0.0;
  return rdm.OneBody(i_rdm, j_rdm);
}

double EOM::RdmTB_J(double J, size_t a_hs, size_t b_hs, size_t c_hs, size_t d_hs) const {
  // Only valence orbits in Hs have a corresponding rdm entry
  if (modelspace->GetOrbit(a_hs).cvq != 1) return 0.0;
  if (modelspace->GetOrbit(b_hs).cvq != 1) return 0.0;
  if (modelspace->GetOrbit(c_hs).cvq != 1) return 0.0;
  if (modelspace->GetOrbit(d_hs).cvq != 1) return 0.0;
  int a = (a_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[a_hs] : -1;
  int b = (b_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[b_hs] : -1;
  int c = (c_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[c_hs] : -1;
  int d = (d_hs < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[d_hs] : -1;
  if (a < 0 || b < 0 || c < 0 || d < 0) return 0.0;
  // Guard: check that the channel J actually exists in rdm_modelspace before
  // calling GetTBME_J_norm.  GetTwoBodyChannelIndex is a pure formula (no
  // bounds check), so passing a J that is larger than any channel in the
  // (small) rdm modelspace would cause an out-of-bounds access inside
  // GetTBME_norm -> segfault.
  {
    Orbit &oa = rdm_ms->GetOrbit(a), &ob = rdm_ms->GetOrbit(b);
    Orbit &oc = rdm_ms->GetOrbit(c), &od = rdm_ms->GetOrbit(d);
    int parity_bra = (oa.l + ob.l) % 2;
    int parity_ket = (oc.l + od.l) % 2;
    int Tz_bra = (oa.tz2 + ob.tz2) / 2;
    int Tz_ket = (oc.tz2 + od.tz2) / 2;
    int Ji = (int)J;
    // Keep as size_t: casting to int breaks the comparison when
    // GetTwoBodyChannelIndex returns a large wrapping value for
    // an invalid (J, parity, Tz) combination in the small rdm modelspace.
    size_t ch_bra = rdm_ms->GetTwoBodyChannelIndex(Ji, parity_bra, Tz_bra);
    size_t ch_ket = rdm_ms->GetTwoBodyChannelIndex(Ji, parity_ket, Tz_ket);
    size_t nch    = rdm_ms->GetNumberTwoBodyChannels();
    if (ch_bra >= nch || ch_ket >= nch) return 0.0;
  }
  return rdm.TwoBody.GetTBME_J_norm((int)J, (int)J, a, b, c, d);
}

double EOM::RdmThreeBody_J(int Jab, size_t a_hs, size_t b_hs, size_t c_hs,
                           int Jed, size_t d_hs, size_t e_hs, size_t f_hs, int twoJ) const {
  if (!use_rdm3 || !rdm.ThreeBody.IsAllocated())
    return 0.0;
  size_t norbits_hs = modelspace->GetNumberOrbits();
  if (a_hs >= norbits_hs || b_hs >= norbits_hs || c_hs >= norbits_hs ||
      d_hs >= norbits_hs || e_hs >= norbits_hs || f_hs >= norbits_hs) return 0.0;

  // All six orbits must exist in the rdm model space.
  for (size_t idx : {a_hs, b_hs, c_hs, d_hs, e_hs, f_hs}) {
    if (modelspace->GetOrbit(idx).cvq != 1) return 0.0;
  }

  // Translate Hs orbit indices to rdm model-space orbit indices.
  auto map = [&](size_t hs_idx) -> int {
    return (hs_idx < hs_to_rdm_orb.size()) ? hs_to_rdm_orb[hs_idx] : -1;
  };
  int a = map(a_hs), b = map(b_hs), c = map(c_hs);
  int d = map(d_hs), e = map(e_hs), f = map(f_hs);
  if (a < 0 || b < 0 || c < 0 || d < 0 || e < 0 || f < 0) return 0.0;
  size_t norbits_rdm = rdm_ms->GetNumberOrbits();
  if (static_cast<size_t>(a) >= norbits_rdm || static_cast<size_t>(b) >= norbits_rdm ||
      static_cast<size_t>(c) >= norbits_rdm || static_cast<size_t>(d) >= norbits_rdm ||
      static_cast<size_t>(e) >= norbits_rdm || static_cast<size_t>(f) >= norbits_rdm) return 0.0;

  // Guard: verify the 3-body channel (twoJ, parity, twoTz) exists in rdm_ms.
  {
    Orbit& oa = rdm_ms->GetOrbit(a); Orbit& ob = rdm_ms->GetOrbit(b); Orbit& oc = rdm_ms->GetOrbit(c);
    Orbit& od = rdm_ms->GetOrbit(d); Orbit& oe = rdm_ms->GetOrbit(e); Orbit& of_ = rdm_ms->GetOrbit(f);
    int par_bra  = (oa.l + ob.l + oc.l) % 2;
    int twoTz_bra = oa.tz2 + ob.tz2 + oc.tz2;
    int par_ket  = (od.l + oe.l + of_.l) % 2;
    int twoTz_ket = od.tz2 + oe.tz2 + of_.tz2;
    size_t ch_bra = rdm_ms->GetThreeBodyChannelIndex(twoJ, par_bra,  twoTz_bra);
    size_t ch_ket = rdm_ms->GetThreeBodyChannelIndex(twoJ, par_ket,  twoTz_ket);
    size_t nch    = rdm_ms->GetNumberThreeBodyChannels();
    if (ch_bra >= nch || ch_ket >= nch) return 0.0;
  }

  // GetME_pn handles SortOrbits + recoupling internally.
  return rdm.ThreeBody.GetME_pn(Jab, Jed, twoJ, a, b, c, d, e, f);
}

///  In case we want to construct the A matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void EOM::SetIncludeConfigs(bool qv, bool ph, bool ppvv, bool pphv, bool pphh) {
  include_qv = qv;
  include_ph = ph;
  include_ppvv = ppvv;
  include_pphv = pphv;
  include_pphh = pphh;
  PrintIncludeConfigs();
}

void EOM::PrintIncludeConfigs() const {
  std::cout << "EOM include configs:"
            << " qv=" << (include_qv ? 1 : 0)
            << " ph=" << (include_ph ? 1 : 0)
            << " ppvv=" << (include_ppvv ? 1 : 0)
            << " pphv=" << (include_pphv ? 1 : 0)
            << " pphh=" << (include_pphh ? 1 : 0) << std::endl;
}

void EOM::ClearConfigBookkeeping() {
  eom_confs.clear();
  eom_dims = 0;
  qv_dim = ph_dim = ppvv_dim = pphv_dim = pphh_dim = 0;
  qv_start = qv_end = 0;
  ph_start = ph_end = 0;
  ppvv_start = ppvv_end = 0;
  pphv_start = pphv_end = 0;
  pphh_start = pphh_end = 0;
}

bool EOM::TensorOneBodyAllowed(const Orbit &oa, const Orbit &oi) const {
  if ((oa.l + oi.l + parity) % 2)
    return false;
  if (std::abs(oa.tz2 - oi.tz2) != 2 * itz)
    return false;
  return AngMom::Triangle(oa.j2 * 0.5, oi.j2 * 0.5, (double)J2);
}

bool EOM::TensorTwoBodyAllowed(const TwoBodyChannel &bra,
                               const TwoBodyChannel &ket) const {
  if (!AngMom::Triangle(bra.J, ket.J, J2))
    return false;
  if (std::abs(bra.Tz - ket.Tz) != itz)
    return false;
  if ((bra.parity + ket.parity + parity) % 2)
    return false;
  return true;
}

void EOM::ConstructConfigs() {
  // Generate configuration for fock space EOM (scalar: bra/ket share JπTz).
  // Blocks gated by SetIncludeConfigs (default: all on).
  ClearConfigBookkeeping();
  eom_tensor_configs = false;

  std::cout << "Constructing EOM configurations for J2=" << J2
            << " parity=" << parity << " itz=" << itz << std::endl;
  PrintIncludeConfigs();

  int norbits = modelspace->norbits;

  // qv: 1b excluded ← valence
  qv_start = 0;
  if (include_qv) {
    for (index_t i_orb = 0; i_orb < norbits; i_orb++) {
      Orbit &oi = modelspace->GetOrbit(i_orb);
      if (oi.cvq != 2)
        continue;
      for (index_t j_orb = 0; j_orb < norbits; j_orb++) {
        Orbit &oj = modelspace->GetOrbit(j_orb);
        if (oj.cvq != 1)
          continue;
        if (oj.l != oi.l)
          continue;
        if (oj.j2 != oi.j2)
          continue;
        if (oj.tz2 != oi.tz2)
          continue;
        eom_confs.push_back({i_orb, j_orb, 0, eom_dims});
        eom_dims += 1;
        qv_dim += 1;
      }
    }
  }
  if (qv_dim > 0)
    qv_end = qv_start + qv_dim - 1;
  std::cout << "dimension EOM qv: " << qv_start << " " << qv_end
            << " (dim=" << qv_dim << ")" << std::endl;

  // ph: 1b (v∪q) ← core
  ph_start = eom_confs.size();
  if (include_ph) {
    for (index_t i_orb = 0; i_orb < norbits; i_orb++) {
      Orbit &oi = modelspace->GetOrbit(i_orb);
      if (oi.cvq == 0)
        continue;
      for (index_t j_orb = 0; j_orb < norbits; j_orb++) {
        Orbit &oj = modelspace->GetOrbit(j_orb);
        if (oj.cvq != 0)
          continue;
        if (oj.l != oi.l)
          continue;
        if (oj.j2 != oi.j2)
          continue;
        if (oj.tz2 != oi.tz2)
          continue;
        eom_confs.push_back({i_orb, j_orb, 0, eom_dims});
        eom_dims += 1;
        ph_dim += 1;
      }
    }
  }
  if (ph_dim > 0)
    ph_end = ph_start + ph_dim - 1;
  std::cout << "dimension EOM ph: " << ph_start << " " << ph_end
            << " (dim=" << ph_dim << ")" << std::endl;

  // ppvv / qqvv: bra qq∪qv, ket vv
  ppvv_start = eom_confs.size();
  size_t number_channels = modelspace->GetNumberTwoBodyChannels();
  if (include_ppvv) {
    for (index_t ich = 0; ich < number_channels; ich++) {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      for (auto &ibra :
           VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv())) {
        for (auto &iket : tbc.GetKetIndex_vv()) {
          eom_confs.push_back({ibra, iket, ich, ich});
          eom_dims += 1;
          ppvv_dim += 1;
        }
      }
    }
  }
  if (ppvv_dim > 0)
    ppvv_end = ppvv_start + ppvv_dim - 1;
  std::cout << "dimension EOM ppvv: " << ppvv_start << " " << ppvv_end
            << " (dim=" << ppvv_dim << ")" << std::endl;

  // pphv / qqhv: bra qq∪qv∪vv, ket vc
  pphv_start = eom_confs.size();
  if (include_pphv) {
    for (index_t ich = 0; ich < number_channels; ich++) {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      for (auto &ibra :
           VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(),
                       tbc.GetKetIndex_vv())) {
        for (auto &iket : tbc.GetKetIndex_vc()) {
          eom_confs.push_back({ibra, iket, ich, ich});
          eom_dims += 1;
          pphv_dim += 1;
        }
      }
    }
  }
  if (pphv_dim > 0)
    pphv_end = pphv_start + pphv_dim - 1;
  std::cout << "dimension EOM pphv: " << pphv_start << " " << pphv_end
            << " (dim=" << pphv_dim << ")" << std::endl;

  // pphh: bra qq∪qv∪vv, ket cc
  pphh_start = eom_confs.size();
  if (include_pphh) {
    for (index_t ich = 0; ich < number_channels; ich++) {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      for (auto &ibra :
           VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(),
                       tbc.GetKetIndex_vv())) {
        for (auto &iket : tbc.GetKetIndex_cc()) {
          eom_confs.push_back({ibra, iket, ich, ich});
          eom_dims += 1;
          pphh_dim += 1;
        }
      }
    }
  }
  if (pphh_dim > 0)
    pphh_end = pphh_start + pphh_dim - 1;

  std::cout << "dimension EOM pphh: " << pphh_start << " " << pphh_end
            << " (dim=" << pphh_dim << ")" << std::endl;

  std::cout << "Total dimension of EOM: " << eom_confs.size() << std::endl;
}

void EOM::ConstructConfigs_tensor() {
  // Tensor EOM: bra/ket JπTz couple to (J2, parity, itz), they do not match.
  // 2b stored as {ibra, iket, ch_bra, ch_ket} with ch_bra <= ch_ket
  // (same triangle as TwoBodyME::Allocate). Off-diagonal also keeps the
  // swapped (ket-type, bra-type) block, which is independent in the stored
  // half-matrix (hermitian conjugate lives in the unstored ch_ket, ch_bra).
  ClearConfigBookkeeping();
  eom_tensor_configs = true;

  std::cout << "Constructing tensor EOM configurations for J=" << J2
            << " parity=" << parity << " Tz=" << itz << std::endl;
  PrintIncludeConfigs();

  int norbits = modelspace->norbits;

  auto add_1b = [&](index_t i_orb, index_t j_orb, index_t &dim) {
    eom_confs.push_back({i_orb, j_orb, 0, eom_dims});
    eom_dims += 1;
    dim += 1;
  };

  // qv: ⟨q|v⟩
  qv_start = 0;
  if (include_qv) {
    for (index_t i_orb = 0; i_orb < norbits; i_orb++) {
      Orbit &oi = modelspace->GetOrbit(i_orb);
      if (oi.cvq != 2)
        continue;
      for (index_t j_orb = 0; j_orb < norbits; j_orb++) {
        Orbit &oj = modelspace->GetOrbit(j_orb);
        if (oj.cvq != 1)
          continue;
        if (!TensorOneBodyAllowed(oi, oj))
          continue;
        add_1b(i_orb, j_orb, qv_dim);
      }
    }
  }
  if (qv_dim > 0)
    qv_end = qv_start + qv_dim - 1;
  std::cout << "dimension EOM qv: " << qv_start << " " << qv_end
            << " (dim=" << qv_dim << ")" << std::endl;

  // ph: ⟨p|h⟩  p = v∪q, h = core
  ph_start = eom_confs.size();
  if (include_ph) {
    for (index_t i_orb = 0; i_orb < norbits; i_orb++) {
      Orbit &oi = modelspace->GetOrbit(i_orb);
      if (oi.cvq == 0)
        continue;
      for (index_t j_orb = 0; j_orb < norbits; j_orb++) {
        Orbit &oj = modelspace->GetOrbit(j_orb);
        if (oj.cvq != 0)
          continue;
        if (!TensorOneBodyAllowed(oi, oj))
          continue;
        add_1b(i_orb, j_orb, ph_dim);
      }
    }
  }
  if (ph_dim > 0)
    ph_end = ph_start + ph_dim - 1;
  std::cout << "dimension EOM ph: " << ph_start << " " << ph_end
            << " (dim=" << ph_dim << ")" << std::endl;

  const size_t nch = modelspace->GetNumberTwoBodyChannels();

  auto add_2b_pair = [&](index_t ibra, index_t iket, size_t ch_bra,
                          size_t ch_ket, index_t &dim) {
    eom_confs.push_back(
        {ibra, iket, (index_t)ch_bra, (index_t)ch_ket});
    eom_dims += 1;
    dim += 1;
  };

  auto add_2b_block = [&](auto get_bras, auto get_kets, index_t &dim) {
    for (size_t ch_bra = 0; ch_bra < nch; ++ch_bra) {
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      for (size_t ch_ket = ch_bra; ch_ket < nch; ++ch_ket) {
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
        if (!TensorTwoBodyAllowed(tbc_bra, tbc_ket))
          continue;
        for (auto &ibra : get_bras(tbc_bra)) {
          for (auto &iket : get_kets(tbc_ket))
            add_2b_pair(ibra, iket, ch_bra, ch_ket, dim);
        }
        if (ch_bra == ch_ket)
          continue;
        for (auto &ibra : get_kets(tbc_bra)) {
          for (auto &iket : get_bras(tbc_ket))
            add_2b_pair(ibra, iket, ch_bra, ch_ket, dim);
        }
      }
    }
  };

  auto bras_qq_qv = [](TwoBodyChannel &tbc) {
    return VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv());
  };
  auto bras_pp = [](TwoBodyChannel &tbc) {
    return VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(),
                       tbc.GetKetIndex_vv());
  };
  auto kets_vv = [](TwoBodyChannel &tbc) { return tbc.GetKetIndex_vv(); };
  auto kets_vc = [](TwoBodyChannel &tbc) { return tbc.GetKetIndex_vc(); };
  auto kets_cc = [](TwoBodyChannel &tbc) { return tbc.GetKetIndex_cc(); };

  // ppvv: ⟨qq,qv|vv⟩
  ppvv_start = eom_confs.size();
  if (include_ppvv)
    add_2b_block(bras_qq_qv, kets_vv, ppvv_dim);
  if (ppvv_dim > 0)
    ppvv_end = ppvv_start + ppvv_dim - 1;
  std::cout << "dimension EOM ppvv: " << ppvv_start << " " << ppvv_end
            << " (dim=" << ppvv_dim << ")" << std::endl;

  // pphv: ⟨qq,qv,vv|hv⟩
  pphv_start = eom_confs.size();
  if (include_pphv)
    add_2b_block(bras_pp, kets_vc, pphv_dim);
  if (pphv_dim > 0)
    pphv_end = pphv_start + pphv_dim - 1;
  std::cout << "dimension EOM pphv: " << pphv_start << " " << pphv_end
            << " (dim=" << pphv_dim << ")" << std::endl;

  // pphh: ⟨pp|hh⟩
  pphh_start = eom_confs.size();
  if (include_pphh)
    add_2b_block(bras_pp, kets_cc, pphh_dim);
  if (pphh_dim > 0)
    pphh_end = pphh_start + pphh_dim - 1;
  std::cout << "dimension EOM pphh: " << pphh_start << " " << pphh_end
            << " (dim=" << pphh_dim << ")" << std::endl;

  std::cout << "Total dimension of tensor EOM: " << eom_confs.size()
            << std::endl;
}

arma::vec EOM::GetEnergies() { return Energies; }

void EOM::ShowModel() {
  std::cout << "Modelspace single-particle orbits:" << std::endl;
  std::cout << "idx n l j t cvq" << std::endl;

  for (size_t i = 0; i < modelspace->GetNumberOrbits(); ++i) {
    Orbit &o = modelspace->GetOrbit(i);
    double j = 0.5 * static_cast<double>(o.j2);
    double t = 0.5 * static_cast<double>(o.tz2);
    std::cout << i << " "
              << o.n << " "
              << o.l << " "
              << j << " "
              << t << " "
              << o.cvq << std::endl;
  }
}

void EOM::ConstructNormMatrix() {
  if (eom_tensor_configs)
    throw std::runtime_error(
        "ConstructNormMatrix: scalar-only; use ConstructNormMatrix_tensor");

  Nkernel.set_size(eom_dims, eom_dims);

  Nkernel.zeros();
  // B4, becnhmarked with srg
  if (qv_dim != 0) {
    for (index_t i = qv_start; i <= qv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);

      for (index_t j = qv_start; j <= qv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        if (cf_bra[0] != cf_ket[0])
          continue;
        Orbit &obra = modelspace->GetOrbit(cf_bra[1]);
        Nkernel(i, j) += RdmOB(cf_bra[1], cf_ket[1]) * sqrt(obra.j2 + 1.);
      }
    }
  }

  // A1, B5, becnhmarked with srg
  if (ph_dim != 0) {
    for (index_t i = ph_start; i <= ph_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      Orbit &obra = modelspace->GetOrbit(cf_bra[1]);
      for (index_t j = ph_start; j <= ph_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);

        if (i == j) {
          Nkernel(i, j) += 1. * (obra.j2 + 1.);
        }
        if (cf_bra[1] != cf_ket[1])
          continue;
        Nkernel(i, j) -= RdmOB(cf_ket[0], cf_bra[0]) * sqrt(obra.j2 + 1.);
      }
    }
  }
  // C1, benchmarked, maybe, it too many
  if (ppvv_dim != 0) {

    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);

      for (index_t j = ppvv_start; j <= ppvv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);

        if (cf_bra[2] != cf_ket[2])
          continue;  // same channel
        if (cf_bra[0] != cf_ket[0])
          continue; 
        TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
        Ket &kbra = tbc_bra.GetKet(cf_bra[1]);
        Ket &kket = tbc_bra.GetKet(cf_ket[1]);
        double val =
            RdmTB_J(tbc_bra.J, kbra.p, kbra.q, kket.p, kket.q) *
            sqrt(2 * tbc_bra.J + 1.);
        Nkernel(i, j) += val;
      }
    }
  }
  // B1, benchmarked
  if (pphv_dim != 0) {

    for (index_t i = pphv_start; i <= pphv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      for (index_t j = pphv_start; j <= pphv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        if (cf_bra[2] != cf_ket[2])
          continue; // must be in the same channel
        if (cf_bra[0] != cf_ket[0])
          continue; // pp=pp

        TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
        Ket &dbra = tbc_bra.GetKet(cf_bra[0]);
        Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
        Ket &dket2 = tbc_bra.GetKet(cf_ket[1]);
        size_t e1 = dket1.p;
        size_t c1 = dket1.q;
        size_t e2 = dket2.p;
        size_t c2 = dket2.q;
        if (e1 != e2)
          continue;
        Orbit &oc1 = modelspace->GetOrbit(c1);
        Orbit &oc2 = modelspace->GetOrbit(c2);
        if (oc1.l != oc2.l)
          continue;
        if (oc1.j2 != oc2.j2)
          continue;
        if (oc1.tz2 != oc2.tz2)
          continue;
        if (oc1.cvq != 1)
          continue; // need be valence
        if (oc2.cvq != 1)
          continue;
        double j1 = tbc_bra.J;

        Nkernel(i, j) +=
            RdmOB(c1, c2) * (2 * tbc_bra.J + 1.) / sqrt(oc1.j2 + 1.);
       
      }
    }
  }
  // C4
  if (pphv_dim != 0) {

    for (index_t i = pphv_start; i <= pphv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a1 = dbra1.p;
      size_t b1 = dbra1.q; // permute ec to make c the valence
      size_t c1 = dket1.p;
      size_t d1 = dket1.q;
      Orbit &oa1 = modelspace->GetOrbit(a1);
      Orbit &oc1 = modelspace->GetOrbit(c1);
      Orbit &ob1 = modelspace->GetOrbit(b1);
      Orbit &od1 = modelspace->GetOrbit(d1);
      if (oa1.cvq != 1 && ob1.cvq != 1)
        continue;
      if (od1.cvq != 1)
        continue;
      double j1 = tbc_bra.J;

      double norm_fact1 = 1.;
      if (a1 == b1)
        norm_fact1 = sqrt(2);

      for (index_t j = pphv_start; j <= pphv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
        Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
        Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
        size_t a2 = dbra2.p;
        size_t b2 = dbra2.q; // permute to make b the valence
        size_t c2 = dket2.p;
        size_t d2 = dket2.q;
        Orbit &oa2 = modelspace->GetOrbit(a2);
        Orbit &oc2 = modelspace->GetOrbit(c2);
        Orbit &ob2 = modelspace->GetOrbit(b2);
        Orbit &od2 = modelspace->GetOrbit(d2);

        if (oa2.cvq != 1 && ob2.cvq != 1)
          continue;
        if (od2.cvq != 1)
          continue;
        double j2 = tbc_ket.J;
        double norm_fact2 = 1.;
        if (a2 == b2)
          norm_fact2 = sqrt(2);
        double norm_fact = norm_fact1 * norm_fact2;
        if (c1 != c2)
          continue;
        double val = 0;
        if (b1 == a2 && oa1.cvq == 1 && ob2.cvq == 1) {
          val += norm_fact * Core_Diagram(d1, b2, a1, d2, a2, c1, j1, j2) *
                 dket1.Phase(j1);
        }
        if (a1 == a2 && ob1.cvq == 1 && ob2.cvq == 1 && a1 != b1) {
          val += norm_fact * Core_Diagram(d1, b2, b1, d2, a2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra1.Phase(j1);
        }
        if (b1 == b2 && oa1.cvq == 1 && oa2.cvq == 1 && a2 != b2) {
          val += norm_fact * Core_Diagram(d1, a2, a1, d2, b2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra2.Phase(j2);
        }
        if (a1 == b2 && ob1.cvq == 1 && oa2.cvq == 1 && a1 != b1 && a2 != b2) {
          val += norm_fact * Core_Diagram(d1, a2, b1, d2, b2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra1.Phase(j1) * dbra2.Phase(j2);
        }
        Nkernel(i, j) += val;
      }
    }
  }
  // A2 C3 benchmarked
  if (pphh_dim != 0) {
    // first We compute A2
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Nkernel(i, i) += (2 * tbc_bra.J + 1.);
    }

    // C3
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a = dbra1.p;
      size_t b = dbra1.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      if (oa.cvq != 1)
        continue; // need be valence
      if (ob.cvq != 1)
        continue;

      for (index_t j = pphh_start; j <= pphh_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        if (cf_ket[2] != cf_bra[2])
          continue;
        if (cf_ket[1] != cf_bra[1])
          continue;
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
        Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
        Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
        size_t c = dbra2.p;
        size_t d = dbra2.q;
        Orbit &oc = modelspace->GetOrbit(c);
        Orbit &od = modelspace->GetOrbit(d);
        if (oc.cvq != 1)
          continue; // need be valence
        if (od.cvq != 1)
          continue;

        double val2 =
            RdmTB_J(tbc_bra.J, a, b, c, d) *
            sqrt(2 * tbc_bra.J + 1.);
        Nkernel(i, j) += val2;
      }
    }
  }

  // B2, benchmarked
  for (index_t i = pphh_start; i <= pphh_end; i++) {
    std::array<index_t, 4> &cf_bra = eom_confs.at(i);
    TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
    Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
    Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
    size_t a = dbra1.p;
    size_t b = dbra1.q;
    Orbit &oa = modelspace->GetOrbit(a);
    Orbit &ob = modelspace->GetOrbit(b);
    if (oa.cvq != 1 && ob.cvq != 1)
      continue; // at least one need be valence
    double norm_fact1 = 1.;
    if (a == b)
      norm_fact1 = sqrt(2.);
    for (index_t j = pphh_start; j <= pphh_end; j++) {
      std::array<index_t, 4> &cf_ket = eom_confs.at(j);
      if (cf_ket[2] != cf_bra[2])
        continue; // same channel
      if (cf_ket[1] != cf_bra[1])
        continue; // same holes
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
      Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
      Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
      size_t c = dbra2.p;
      size_t d = dbra2.q;
      double norm_fact2 = 1.;
      if (c == d)
        norm_fact2 = sqrt(2.);

      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &od = modelspace->GetOrbit(d);

      // std::cout<<" start "<< std::endl;
      double norm_fact = norm_fact1 * norm_fact2;
      if (b == d) {
        Nkernel(i, j) -= norm_fact * RdmOB(c, a) * (2 * tbc_bra.J + 1.) /
                         sqrt(oc.j2 + 1.);
      }

      if (b != a && ob.cvq == 1 && a == d) {
        Nkernel(i, j) -= norm_fact * RdmOB(c, b) * (2 * tbc_bra.J + 1.) *
                         dbra1.Phase(tbc_bra.J) / sqrt(oc.j2 + 1.);
      }
      if (c != d && od.cvq == 1 && b == c) {
        Nkernel(i, j) -= norm_fact * RdmOB(d, a) * (2 * tbc_bra.J + 1.) *
                         dbra2.Phase(tbc_ket.J) / sqrt(od.j2 + 1.);
      }
      if (b != a && c != d && od.cvq == 1 && ob.cvq == 1 && a == c) {
        Nkernel(i, j) -= norm_fact * RdmOB(d, b) * (2 * tbc_bra.J + 1.) *
                         dbra2.Phase(tbc_ket.J) * dbra1.Phase(tbc_bra.J) /
                         sqrt(od.j2 + 1.);
      }
    }
  }

  // now off diagonal B3 benchmarked
  if (pphv_dim != 0 && ph_dim != 0) {

    for (index_t i = pphv_start; i <= pphv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a = dbra1.p;
      size_t b = dbra1.q; // permute ec to make c the valence
      size_t c = dket1.p;
      size_t d = dket1.q;
      double norm_fact = 1;
      if (a == b)
        norm_fact = sqrt(2.);

      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &od = modelspace->GetOrbit(d);
      double j1 = tbc_bra.J;
      if (oa.cvq != 1 && ob.cvq != 1)
        continue;
      if (od.cvq != 1)
        continue;
      for (index_t j = ph_start; j <= ph_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        size_t c1 = cf_ket[0];
        size_t b1 = cf_ket[1];
        if (b1 != c)
          continue;

        if (c1 == b) {
          Nkernel(i, j) += norm_fact * dbra1.Phase(j1) * RdmOB(a, d) *
                           (2 * j1 + 1.) / sqrt(oa.j2 + 1.);
          Nkernel(j, i) += norm_fact * dbra1.Phase(j1) * RdmOB(a, d) *
                           (2 * j1 + 1.) / sqrt(oa.j2 + 1.);
        }
        if (c1 == a && a != b) {
          Nkernel(i, j) +=
              norm_fact * RdmOB(b, d) * (2 * j1 + 1.) / sqrt(ob.j2 + 1.);
          Nkernel(j, i) +=
              norm_fact * RdmOB(b, d) * (2 * j1 + 1.) / sqrt(ob.j2 + 1.);
        }
      }
    }
  }
  //
  //
  // C5
  if (pphv_dim != 0 && ph_dim != 0) {

    for (index_t i = pphv_start; i <= pphv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a = dbra1.p;
      size_t b = dbra1.q; // permute ec to make c the valence
      size_t c = dket1.p;
      size_t d = dket1.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &od = modelspace->GetOrbit(d);

      double j1 = tbc_bra.J;
      if (oa.cvq != 1 || ob.cvq != 1 || od.cvq != 1)
        continue;
      for (index_t j = ph_start; j <= ph_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        size_t c1 = cf_ket[0];
        size_t b1 = cf_ket[1];
        if (b1 == c) {
          double norm_fact = 1.;
          if (c1 == d)
            norm_fact = sqrt(2.);
          double val2 =
              RdmTB_J(tbc_bra.J, a, b, c1, d) *
              sqrt(j1 * 2. + 1.);
          Nkernel(i, j) -= val2 * norm_fact;
          Nkernel(j, i) -= val2 * norm_fact;
        }
      }
    }
  }


  // C1 threebody diagram (valence ρ₃)
  if (use_rdm3 && rdm.ThreeBody.IsAllocated() && ppvv_dim != 0) {
    std::cout << "start three body norm " << std::endl;
    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a1 = dbra1.p;
      size_t b1 = dbra1.q; 
      size_t c1 = dket1.p;
      size_t d1 = dket1.q;
      Orbit &oa1 = modelspace->GetOrbit(a1);

      for (index_t j = ppvv_start; j <= ppvv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
      Ket &dbra1 = tbc_ket.GetKet(cf_ket[0]);
      Ket &dket1 = tbc_ket.GetKet(cf_ket[1]);
      size_t a2 = dbra1.p;
      size_t b2 = dbra1.q; 
      size_t c2 = dket1.p;
      size_t d2 = dket1.q;
      Orbit &oa2 = modelspace->GetOrbit(a2);
      
      
      if(b1==b2 && oa1.cvq==1 && oa2.cvq==1){ // both valence not contributing to norm.
      // here we use ThreeBody_Diagram_Entries_Internal, to represent a contraction, instead of ThreeBody_Diagram that mimic comm223ss
      double val = ThreeBody_Diagram_Entries_Internal(c1, d1, a2, a1, d2, c2, b1,
                                     tbc_bra.J, tbc_ket.J);
       Nkernel(i, j) += val;
        //std::cout<<"Threebody: "<< a1<<" "<<b1<<" "<<c1<<" "<<d1<<" "<<a2<<" "<<b2<<" "<<c2<<" "<<d2<<" "<< tbc_bra.J<<" "<<tbc_ket.J<<" " <<val<<std::endl;
    }
        
      }
    }
  } else if (!use_rdm3) {
    std::cout << "skip three body norm (SetUseRdm3 false)" << std::endl;
  }


  // c2
  if (ppvv_dim != 0 && qv_dim != 0) {

    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a = dbra1.p;
      size_t b = dbra1.q; // permute ec to make c the valence
      size_t c = dket1.p;
      size_t d = dket1.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &od = modelspace->GetOrbit(d);

      double j1 = tbc_bra.J;
      if (oa.cvq != 1 || oc.cvq != 1 || od.cvq != 1)
        continue;
      for (index_t j = qv_start; j <= qv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        size_t c1 = cf_ket[0];
        size_t b1 = cf_ket[1];
        if (b == c1) {
          double norm_fact = 1.;
          if (a == b1)
            norm_fact = sqrt(2);
          double val2 =
              RdmTB_J(tbc_bra.J, a, b1, c, d) *
              sqrt(2 * j1 + 1.);
          Nkernel(i, j) += val2 * norm_fact;
          Nkernel(j, i) += val2 * norm_fact;
        }
      }
    }



  }
}



void EOM::ConstructNormMatrix_tensor() {
  // Tensor N-kernel: same Wick strings as ConstructNormMatrix, AMC with
  // X,Y rank λ (Y index-reversed so λ×λ→0). Formulas:
  //   learn/amc_tts/eom_norm/output/*_tensor.tex
  // ρ is always scalar 0+.
  // Empirical hats vs ½[Q⁺,Q⁻] (He4-core, TTS leftover × GetVSEOM_Overlap):
  //   110/220 (A1, A2): λ̂^{-2}  (AMC had λ̂^{-1}; 0-body m-trace is λ̂^{-2})
  //   leftover 2b×ρ C1 (comm222 × ρ·Ĵ): Ĵ_vv^{-1} λ̂^{-1}
  //     (not λ̂^{-2}; that only matches when J_vv=λ)
  //   leftover 2b×ρ C3 (comm222 × ρ·Ĵ): Ĵ_pp^{-1} λ̂^{-1}
  //     (same leftover packaging as C1; λ̂^{-2} only matches when J_pp=λ)
  //   C2/C5 (comm122 leftover × ρ·Ĵ): Ĵ_other λ̂^{-1}
  //   C4 is the unreduced XY scalar, λ̂^{-1}, not the XY−YX commutator
  //   leftover 1b×ρ (B4, B1, B2, B5, B3): λ̂^{-1} / ĵ of the contracted orbit
  // C4 is the scalar four-perm pipeline on Core_Diagram_tensor
  // (X daggered: X_cdab ρ_dfae Y_bfce).
  // Y_ia / Y_ijab folded into stored Y_ai / Y_abij via IMSRG tensor flip:
  //   1b: O_ji = (-1)^{j_i-j_j} O_ij   (Hermitian)
  //   2b: ⟨J1||O||J0⟩ = (-1)^{J1-J0} ⟨J0||O||J1⟩
  // ρ3 not wired yet.
  if (!eom_tensor_configs)
    throw std::runtime_error(
        "ConstructNormMatrix_tensor: call ConstructConfigs_tensor first");

  Nkernel.set_size(eom_dims, eom_dims);
  Nkernel.zeros();

  const int lam = J2;
  const double lamhat = std::sqrt(2.0 * lam + 1.0);
  const double lamhatinv = 1.0 / lamhat;
  const double lamhatinv2 = lamhatinv * lamhatinv;
  auto hat = [](int J) { return std::sqrt(2.0 * J + 1.0); };
  auto jhat = [](const Orbit &o) { return std::sqrt(o.j2 + 1.0); };
  auto ph = [](int n) { return (double)AngMom::phase(n); };

  // Fold Y_ia into stored Y_ai (SetOneBody).
  auto we1 = [&](const Orbit &oa, const Orbit &oi) {
    return ph((oa.j2 - oi.j2) / 2);
  };
  // Fold Y_ijab^{J1 J0} into stored Y_abij^{J0 J1} (GetTwoBody).
  auto we2 = [&](int J0, int J1) { return ph(J1 - J0); };

  // Tensor pphv stores both ⟨pp||χ||hv⟩ and the independent ⟨hv||χ||pp⟩.
  // C4/B1 Wick strings are X_abcd with (ab)=pp, (cd)=hv. Swapped-type
  // configs must not be labeled as pp-hv (scalar hermiticity did that).
  auto is_pp_ket = [&](const Ket &k) {
    return modelspace->GetOrbit(k.p).cvq != 0 &&
           modelspace->GetOrbit(k.q).cvq != 0;
  };
  auto is_hv_ket = [&](const Ket &k) {
    int ca = modelspace->GetOrbit(k.p).cvq;
    int cb = modelspace->GetOrbit(k.q).cvq;
    return (ca + cb == 1);
  };
  auto pphv_is_pp_hv = [&](index_t idx) {
    auto &cf = eom_confs.at(idx);
    Ket &kbra = modelspace->GetTwoBodyChannel(cf[2]).GetKet(cf[0]);
    Ket &kket = modelspace->GetTwoBodyChannel(cf[3]).GetKet(cf[1]);
    return is_pp_ket(kbra) && is_hv_ket(kket);
  };

  // B4 qv-qv: 110 identity (Fermi n_v(1-n_q), same hats as A1) plus leftover 1b×ρ
  if (qv_dim != 0) {
    for (index_t i = qv_start; i <= qv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      Orbit &oa = modelspace->GetOrbit(cf_bra[0]);
      Orbit &oi = modelspace->GetOrbit(cf_bra[1]);
      const double phase = ph((oa.j2 + oi.j2) / 2 + lam);
      const double wocc = oi.occ * (1.0 - oa.occ);
      if (std::abs(wocc) > 1e-8)
        Nkernel(i, i) += -phase * lamhatinv2 * we1(oa, oi) * wocc;
      for (index_t j = qv_start; j <= qv_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        if (cf_bra[0] != cf_ket[0])
          continue;
        Orbit &oj = modelspace->GetOrbit(cf_ket[1]);
        if (oj.j2 != oi.j2)
          continue;
        Nkernel(i, j) +=
            -phase * lamhatinv / jhat(oi) *
            we1(oa, oj) * RdmOB(cf_bra[1], cf_ket[1]);
      }
    }
  }

  // A1 (110): λ̂^{-2}. B5 leftover 1b×ρ: λ̂^{-1} / ĵ_a
  if (ph_dim != 0) {
    for (index_t i = ph_start; i <= ph_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      Orbit &oa = modelspace->GetOrbit(cf_bra[0]);
      Orbit &oi = modelspace->GetOrbit(cf_bra[1]);
      const double phase = ph((oa.j2 + oi.j2) / 2 + lam);
      for (index_t j = ph_start; j <= ph_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        if (cf_bra[1] != cf_ket[1])
          continue;
        Orbit &ob = modelspace->GetOrbit(cf_ket[0]);
        if (i == j)
          Nkernel(i, j) += -phase * lamhatinv2 * we1(oa, oi);
        if (ob.j2 != oa.j2)
          continue;
        Nkernel(i, j) +=
            phase * lamhatinv / jhat(oa) * we1(ob, oi) *
            RdmOB(cf_ket[0], cf_bra[0]);
      }
    }
  }

  // C1 ppvv leftover 2b×ρ. Plus 220 identity on ⟨qq||vv⟩
  // (ladder drops swapped ⟨vv||qq⟩, so do not put identity there).
  // comm222 leftover at J_vv has Ĵ_vv^{-2} λ̂^{-1}; overlap supplies
  // ρ·Ĵ_vv, so N keeps Ĵ_vv^{-1} λ̂^{-1}. Using λ̂^{-2} only matches
  // when J_vv=λ (the P=1 random-op case).
  if (ppvv_dim != 0) {
    auto ket_qq_qv = [&](const Ket &k) {
      int ca = modelspace->GetOrbit(k.p).cvq;
      int cb = modelspace->GetOrbit(k.q).cvq;
      return ca >= 2 || cb >= 2;
    };
    auto ket_vv = [&](const Ket &k) {
      return modelspace->GetOrbit(k.p).cvq == 1 &&
             modelspace->GetOrbit(k.q).cvq == 1;
    };
    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &k0 = tbX.GetKet(cf_bra[0]);
      Ket &k1 = tkX.GetKet(cf_bra[1]);
      const int J0 = tbX.J;
      const int J1 = tkX.J;
      if (ket_qq_qv(k0) && ket_vv(k1)) {
        Orbit &oq1 = modelspace->GetOrbit(k0.p);
        Orbit &oq2 = modelspace->GetOrbit(k0.q);
        Orbit &ov1 = modelspace->GetOrbit(k1.p);
        Orbit &ov2 = modelspace->GetOrbit(k1.q);
        const double wocc =
            ov1.occ * ov2.occ * (1.0 - oq1.occ) * (1.0 - oq2.occ);
        if (std::abs(wocc) > 1e-8)
          Nkernel(i, i) +=
              ph(J0 + J1 + lam) * lamhatinv2 * we2(J0, J1) * wocc;
      }
      if (!(ket_qq_qv(k0) && ket_vv(k1)))
        continue;
      Ket &kvvX = k1;
      for (index_t j = ppvv_start; j <= ppvv_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        if (cf_bra[0] != cf_ket[0] || cf_bra[2] != cf_ket[2])
          continue;
        TwoBodyChannel &tkY = modelspace->GetTwoBodyChannel(cf_ket[3]);
        if (tkY.J != J1)
          continue;
        Ket &k0Y = modelspace->GetTwoBodyChannel(cf_ket[2]).GetKet(cf_ket[0]);
        Ket &k1Y = tkY.GetKet(cf_ket[1]);
        if (!(ket_qq_qv(k0Y) && ket_vv(k1Y)))
          continue;
        Ket &kvvY = k1Y;
        double val = RdmTB_J(J1, kvvX.p, kvvX.q, kvvY.p, kvvY.q);
        // qv bra (exactly one occupied leg): leftover 2b flips vs qq||vv.
        const int nocc_bra =
            (modelspace->GetOrbit(k0.p).occ > 1e-8 ? 1 : 0) +
            (modelspace->GetOrbit(k0.q).occ > 1e-8 ? 1 : 0);
        const double qv_sign = (nocc_bra == 1) ? -1.0 : 1.0;
        Nkernel(i, j) += qv_sign * ph(J0 + J1 + lam) * (lamhatinv / hat(J1)) *
                         we2(J0, J1) * val;
      }
    }
    // ppvv leftover 1b×ρ: B2 on valence, spectator is the other valence.
    // comm221 Nocc = n̄_q n̄_q n_spectator → keep occupied spectator
    // (opposite of pphh, where the pair is hh and Nocc needs n̄_spectator).
    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kqq = tbX.GetKet(cf_bra[0]);
      Ket &kvv = tkX.GetKet(cf_bra[1]);
      if (!(ket_qq_qv(kqq) && ket_vv(kvv)))
        continue;
      // comm221 pair is the qq ket: Nocc ∝ n̄_q1 n̄_q2. A qv bra has n̄_v=0.
      if ((1.0 - modelspace->GetOrbit(kqq.p).occ) < 1e-8 ||
          (1.0 - modelspace->GetOrbit(kqq.q).occ) < 1e-8)
        continue;
      size_t a = kvv.p, b = kvv.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      const int J0 = tbX.J;
      const int J1 = tkX.J;
      const double ph_vv = kvv.Phase(J1);
      const double w0 =
          ph(J0 + J1 + lam) * lamhatinv * we2(J0, J1);
      const double norm_fact1 = (a == b) ? std::sqrt(2.) : 1.;
      for (index_t j = ppvv_start; j <= ppvv_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        if (cf_bra[0] != cf_ket[0] || cf_bra[2] != cf_ket[2])
          continue;
        TwoBodyChannel &tkY = modelspace->GetTwoBodyChannel(cf_ket[3]);
        if (tkY.J != J1)
          continue;
        Ket &k0Y = modelspace->GetTwoBodyChannel(cf_ket[2]).GetKet(cf_ket[0]);
        Ket &kvvY = tkY.GetKet(cf_ket[1]);
        if (!(ket_qq_qv(k0Y) && ket_vv(kvvY)))
          continue;
        size_t c = kvvY.p, d = kvvY.q;
        Orbit &oc = modelspace->GetOrbit(c);
        Orbit &od = modelspace->GetOrbit(d);
        const double ph_vvY = kvvY.Phase(J1);
        const double nf = norm_fact1 * ((c == d) ? std::sqrt(2.) : 1.);
        const bool occ_b = ob.occ > 1e-8;
        const bool occ_a = oa.occ > 1e-8;
        if (occ_b && b == d && oc.j2 == oa.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(c, a) / jhat(oa);
        if (occ_a && b != a && a == d && oc.j2 == ob.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(c, b) * ph_vv / jhat(ob);
        if (occ_b && c != d && b == c && od.j2 == oa.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(d, a) * ph_vvY / jhat(oa);
        if (occ_a && b != a && c != d && a == c && od.j2 == ob.j2)
          Nkernel(i, j) +=
              w0 * nf * RdmOB(d, b) * ph_vv * ph_vvY / jhat(ob);
      }
    }
  }

  // B1 pphv leftover 1b×ρ: λ̂^{-1} / ĵ_c  (+ 220 identity)
  if (pphv_dim != 0 && include_norm_B1) {
    for (index_t i = pphv_start; i <= pphv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kpp = tbX.GetKet(cf_bra[0]);
      Ket &kvcX = tkX.GetKet(cf_bra[1]);
      if (pphv_is_pp_hv(i)) {
        Orbit &oa = modelspace->GetOrbit(kpp.p);
        Orbit &ob = modelspace->GetOrbit(kpp.q);
        Orbit &oh = modelspace->GetOrbit(kvcX.p);
        Orbit &ov = modelspace->GetOrbit(kvcX.q);
        const double wocc =
            oh.occ * ov.occ * (1.0 - oa.occ) * (1.0 - ob.occ);
        if (std::abs(wocc) > 1e-8)
          Nkernel(i, i) += ph(tbX.J + tkX.J + lam) * lamhatinv2 *
                           we2(tbX.J, tkX.J) * wocc;
      }
      size_t e1 = kvcX.p;
      size_t c1 = kvcX.q;
      Orbit &oc1 = modelspace->GetOrbit(c1);
      if (oc1.cvq != 1)
        continue;
      if (!pphv_is_pp_hv(i))
        continue;
      const int J1 = tkX.J;
      for (index_t j = pphv_start; j <= pphv_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        if (cf_bra[0] != cf_ket[0] || cf_bra[2] != cf_ket[2])
          continue;
        TwoBodyChannel &tkY = modelspace->GetTwoBodyChannel(cf_ket[3]);
        if (tkY.J != J1)
          continue;
        Ket &kvcY = tkY.GetKet(cf_ket[1]);
        if (kvcY.p != e1)
          continue;
        size_t c2 = kvcY.q;
        Orbit &oc2 = modelspace->GetOrbit(c2);
        if (oc2.cvq != 1 || oc2.j2 != oc1.j2 || oc1.l != oc2.l ||
            oc1.tz2 != oc2.tz2)
          continue;
        // Dagger ⟨ab|hv⟩† → ⟨hv|ab⟩ gives (-1)^{J_ab-J_hv}.
        // AMC of X_eiab ρ Y_abej: J_0 couples (e,i), J_1 couples (a,b),
        // phase (-1)^{J_0+J_1+λ} λ̂^{-1}/ĵ_c.
        const int Jab = tbX.J;
        const int Jhv = J1;
        Nkernel(i, j) +=
            ph(Jab - Jhv) * ph(Jhv + Jab + lam) * lamhatinv / jhat(oc1) *
            RdmOB(c1, c2);
      }
    }
  }

  // C4 pphv-pphv. Same four permutations as the scalar loop.
  // Each perm is one (e,f) of ⟨af|ce⟩^{j1 j2} ⟨eb|fd⟩^{j3 j4}.
  // Dagger phase (−1)^{J_ab−J_hv} and Ket::Phase sit in the caller.
  if (pphv_dim != 0 && include_norm_C4) {
    for (index_t i = pphv_start; i <= pphv_end; i++) {
      if (!pphv_is_pp_hv(i))
        continue;
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kpp1 = tbX.GetKet(cf_bra[0]);
      Ket &khv1 = tkX.GetKet(cf_bra[1]);
      const size_t a1 = kpp1.p, b1 = kpp1.q, h1 = khv1.p, d1 = khv1.q;
      Orbit &oa1 = modelspace->GetOrbit(a1);
      Orbit &ob1 = modelspace->GetOrbit(b1);
      Orbit &od1 = modelspace->GetOrbit(d1);
      if (oa1.cvq != 1 && ob1.cvq != 1)
        continue;
      if (od1.cvq != 1)
        continue;
      const int JabX = tbX.J;
      const int JhvX = tkX.J;
      const double ph_dag = ph(JabX - JhvX);
      const double ph_hv = (double)khv1.Phase(JhvX);
      const double norm1 = (a1 == b1) ? std::sqrt(2.) : 1.;
      for (index_t j = pphv_start; j <= pphv_end; j++) {
        if (!pphv_is_pp_hv(j))
          continue;
        auto &cf_ket = eom_confs.at(j);
        TwoBodyChannel &tbY = modelspace->GetTwoBodyChannel(cf_ket[2]);
        TwoBodyChannel &tkY = modelspace->GetTwoBodyChannel(cf_ket[3]);
        Ket &kpp2 = tbY.GetKet(cf_ket[0]);
        Ket &khv2 = tkY.GetKet(cf_ket[1]);
        const size_t a2 = kpp2.p, b2 = kpp2.q, h2 = khv2.p, d2 = khv2.q;
        Orbit &oa2 = modelspace->GetOrbit(a2);
        Orbit &ob2 = modelspace->GetOrbit(b2);
        Orbit &od2 = modelspace->GetOrbit(d2);
        if (oa2.cvq != 1 && ob2.cvq != 1)
          continue;
        if (od2.cvq != 1)
          continue;
        if (h1 != h2)
          continue;
        const int JabY = tbY.J;
        const int JhvY = tkY.J;
        const double norm =
            norm1 * ((a2 == b2) ? std::sqrt(2.) : 1.);
        const double base = ph_dag * ph_hv * norm;
        double val = 0.;
        // Same orbit slots as scalar Core_Diagram: ⟨af|ce⟩^{j1 j2} ⟨eb|fd⟩^{j3 j4}.
        // j1 couples (a,f)=(d1,h), j2 couples (c,e), j3 couples (e,b), j4 couples (f,d).
        if (b1 == a2 && oa1.cvq == 1 && ob2.cvq == 1)
          val += base * Core_Diagram_tensor(d1, b2, a1, d2, a2, h1, JhvX,
                                            JabX, JabY, JhvY, lam);
        if (a1 == a2 && ob1.cvq == 1 && ob2.cvq == 1 && a1 != b1)
          val += base * kpp1.Phase(JabX) *
                 Core_Diagram_tensor(d1, b2, b1, d2, a2, h1, JhvX, JabX,
                                     JabY, JhvY, lam);
        if (b1 == b2 && oa1.cvq == 1 && oa2.cvq == 1 && a2 != b2)
          val += base * kpp2.Phase(JabY) *
                 Core_Diagram_tensor(d1, a2, a1, d2, b2, h1, JhvX, JabX,
                                     JabY, JhvY, lam);
        if (a1 == b2 && ob1.cvq == 1 && oa2.cvq == 1 && a1 != b1 &&
            a2 != b2)
          val += base * kpp1.Phase(JabX) * kpp2.Phase(JabY) *
                 Core_Diagram_tensor(d1, a2, b1, d2, b2, h1, JhvX, JabX,
                                     JabY, JhvY, lam);
        Nkernel(i, j) += val;
      }
    }
  }

  // Logical pp/hh legs. Tensor storage also keeps the swapped ⟨hh||χ||pp⟩
  // block (ch_bra ≤ ch_ket); those must not be read as if cf[0] were particles.
  auto ket_noncore = [&](const Ket &k) {
    return modelspace->GetOrbit(k.p).cvq != 0 &&
           modelspace->GetOrbit(k.q).cvq != 0;
  };
  auto ket_hh = [&](const Ket &k) {
    return modelspace->GetOrbit(k.p).cvq == 0 &&
           modelspace->GetOrbit(k.q).cvq == 0;
  };
  auto pphh_legs = [&](index_t idx, size_t &a, size_t &b, size_t &h1,
                       size_t &h2, int &Jpp, int &Jhh, int &Jst0, int &Jst1,
                       double &ph_pp, double &ph_hh_unused) -> bool {
    (void)ph_hh_unused;
    auto &cf = eom_confs.at(idx);
    TwoBodyChannel &tb = modelspace->GetTwoBodyChannel(cf[2]);
    TwoBodyChannel &tk = modelspace->GetTwoBodyChannel(cf[3]);
    Ket &k0 = tb.GetKet(cf[0]);
    Ket &k1 = tk.GetKet(cf[1]);
    Jst0 = tb.J;
    Jst1 = tk.J;
    if (ket_noncore(k0) && ket_hh(k1)) {
      a = k0.p;
      b = k0.q;
      h1 = k1.p;
      h2 = k1.q;
      Jpp = tb.J;
      Jhh = tk.J;
      ph_pp = k0.Phase(Jpp);
      return true;
    }
    if (ket_hh(k0) && ket_noncore(k1)) {
      a = k1.p;
      b = k1.q;
      h1 = k0.p;
      h2 = k0.q;
      Jpp = tk.J;
      Jhh = tb.J;
      ph_pp = k1.Phase(Jpp);
      return true;
    }
    return false;
  };

  // A2 pphh identity (220): λ̂^{-2}, Fermi n_h n_h nbar_p nbar_p
  if (pphh_dim != 0) {
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      size_t a, b, h1, h2;
      int Jpp, Jhh, Jst0, Jst1;
      double ph_pp, dummy = 0.0;
      if (!pphh_legs(i, a, b, h1, h2, Jpp, Jhh, Jst0, Jst1, ph_pp, dummy))
        continue;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &oh1 = modelspace->GetOrbit(h1);
      Orbit &oh2 = modelspace->GetOrbit(h2);
      const double wocc = oh1.occ * oh2.occ * (1.0 - oa.occ) * (1.0 - ob.occ);
      if (std::abs(wocc) < 1e-8)
        continue;
      Nkernel(i, i) +=
          ph(Jst0 + Jst1 + lam) * lamhatinv2 * we2(Jst0, Jst1) * wocc;
    }
  }

  // C3 pphh leftover 2b×ρ. Same leftover packaging as C1: comm222 at J_pp
  // has Ĵ_pp^{-2} λ̂^{-1}; overlap supplies ρ·Ĵ_pp, so N keeps
  // Ĵ_pp^{-1} λ̂^{-1}. Using λ̂^{-2} only matches when J_pp=λ (T=0
  // even-P cases where the only firing C3 is J_pp=λ). Even-P T=1 C3
  // lives on p-shell vv (even pair parity) at J_pp=2 ≠ λ=3.
  if (pphh_dim != 0) {
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      size_t a, b, h1, h2;
      int Jpp, Jhh, Jst0, Jst1;
      double ph_pp, dummy = 0.0;
      if (!pphh_legs(i, a, b, h1, h2, Jpp, Jhh, Jst0, Jst1, ph_pp, dummy))
        continue;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      if (oa.cvq != 1 || ob.cvq != 1)
        continue;
      for (index_t j = pphh_start; j <= pphh_end; j++) {
        size_t c, d, h1b, h2b;
        int JppY, JhhY, Jst0Y, Jst1Y;
        double ph_ppY;
        if (!pphh_legs(j, c, d, h1b, h2b, JppY, JhhY, Jst0Y, Jst1Y, ph_ppY,
                       dummy))
          continue;
        if (h1b != h1 || h2b != h2 || JhhY != Jhh || JppY != Jpp)
          continue;
        Orbit &oc = modelspace->GetOrbit(c);
        Orbit &od = modelspace->GetOrbit(d);
        if (oc.cvq != 1 || od.cvq != 1)
          continue;
        double val = RdmTB_J(Jpp, a, b, c, d);
        Nkernel(i, j) +=
            ph(Jst0 + Jst1 + lam) * (lamhatinv / hat(Jpp)) * we2(Jst0, Jst1) *
            val;
      }
    }
  }

  // B2 pphh leftover 1b×ρ: λ̂^{-1} / ĵ of the contracted particle
  if (pphh_dim != 0) {
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      size_t a, b, h1, h2;
      int Jpp, Jhh, Jst0, Jst1;
      double ph_pp, dummy = 0.0;
      if (!pphh_legs(i, a, b, h1, h2, Jpp, Jhh, Jst0, Jst1, ph_pp, dummy))
        continue;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      if (oa.cvq != 1 && ob.cvq != 1)
        continue;
      double norm_fact1 = (a == b) ? std::sqrt(2.) : 1.;
      const double w0 =
          -ph(Jst0 + Jst1 + lam) * lamhatinv * we2(Jst0, Jst1);
      for (index_t j = pphh_start; j <= pphh_end; j++) {
        size_t c, d, h1b, h2b;
        int JppY, JhhY, Jst0Y, Jst1Y;
        double ph_ppY;
        if (!pphh_legs(j, c, d, h1b, h2b, JppY, JhhY, Jst0Y, Jst1Y, ph_ppY,
                       dummy))
          continue;
        if (h1b != h1 || h2b != h2 || JhhY != Jhh || JppY != Jpp)
          continue;
        double norm_fact2 = (c == d) ? std::sqrt(2.) : 1.;
        Orbit &oc = modelspace->GetOrbit(c);
        Orbit &od = modelspace->GetOrbit(d);
        double nf = norm_fact1 * norm_fact2;
        // comm221 Nocc needs a Fermi-empty spectator on the uncontracted
        // particle. Occupied-occupied (He8 valence) has Nocc=0.
        const bool nbar_b = (1.0 - ob.occ) > 1e-8;
        const bool nbar_a = (1.0 - oa.occ) > 1e-8;
        if (nbar_b && b == d && oc.j2 == oa.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(c, a) / jhat(oa);
        if (nbar_a && b != a && ob.cvq == 1 && a == d && oc.j2 == ob.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(c, b) * ph_pp / jhat(ob);
        if (nbar_b && c != d && od.cvq == 1 && b == c && od.j2 == oa.j2)
          Nkernel(i, j) += w0 * nf * RdmOB(d, a) * ph_ppY / jhat(oa);
        if (nbar_a && b != a && c != d && od.cvq == 1 && ob.cvq == 1 &&
            a == c && od.j2 == ob.j2)
          Nkernel(i, j) +=
              w0 * nf * RdmOB(d, b) * ph_pp * ph_ppY / jhat(ob);
      }
    }
  }

  // B3 pphv-ph leftover 1b×ρ.
  // Bra is stored ⟨ab|cd⟩ (pp-hv). Dagger → X_cdab, phase (−1)^{J_ab−J_hv}.
  // Ket is stored ph Y_ac. AMC: X_cdab ρ_bd Y_ac,
  //   δ_{jd jb} (−1)^{J_ab+ja+jb} Ĵ_ab Ĵ_hv λ̂^{-1}
  //   6j{λ J_ab J_hv; jb jc ja} / ĵ_b
  // Other particle: Y_bc, ρ_ad, Ket::Phase on the pp swap.
  if (pphv_dim != 0 && ph_dim != 0) {
    for (index_t i = pphv_start; i <= pphv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kpp = tbX.GetKet(cf_bra[0]);
      Ket &kvc = tkX.GetKet(cf_bra[1]);
      size_t a = kpp.p, b = kpp.q, c = kvc.p, d = kvc.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &od = modelspace->GetOrbit(d);
      if (oa.cvq != 1 && ob.cvq != 1)
        continue;
      if (od.cvq != 1)
        continue;
      if (!pphv_is_pp_hv(i))
        continue;
      const int Jab = tbX.J;
      const int Jhv = tkX.J;
      const double ja = oa.j2 * 0.5, jb = ob.j2 * 0.5, jc = oc.j2 * 0.5;
      const double ph_dag = ph(Jab - Jhv);
      const double ph_amc = ph(Jab + (oa.j2 + ob.j2) / 2);
      const double hats = hat(Jab) * hat(Jhv) * lamhatinv;
      const double norm_fact = (a == b) ? std::sqrt(2.) : 1.;
      for (index_t j = ph_start; j <= ph_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        size_t c1 = cf_ket[0];
        size_t b1 = cf_ket[1];
        if (b1 != c)
          continue;
        if (c1 == a && od.j2 == ob.j2) {
          double six =
              AngMom::SixJ((double)lam, (double)Jab, (double)Jhv, jb, jc, ja);
          double val = ph_dag * ph_amc * hats / jhat(ob) * six * norm_fact *
                       RdmOB(b, d);
          Nkernel(i, j) += val;
          Nkernel(j, i) += val;
        }
        if (c1 == b && a != b && od.j2 == oa.j2) {
          double six =
              AngMom::SixJ((double)lam, (double)Jab, (double)Jhv, ja, jc, jb);
          double val = ph_dag * ph_amc * hats / jhat(oa) * six * norm_fact *
                       kpp.Phase(Jab) * RdmOB(a, d);
          Nkernel(i, j) += val;
          Nkernel(j, i) += val;
        }
      }
    }
  }

  // C5 pphv-ph leftover 2b×ρ.
  // Bra is stored ⟨ab|cd⟩. Dagger → X_cdab, phase (−1)^{J_ab−J_hv}.
  // Ket is stored ph Y_ec. AMC of −X_cdab ρ_abed Y_ec:
  //   −(−1)^{J_ab+jd+je} Ĵ_hv Ĵ_ab λ̂^{-1} 6j{J_hv J_ab λ; je jc jd}
  // comm122 leftover at J_ab has Ĵ_ab^{-1} Ĵ_hv λ̂^{-1}; overlap supplies
  // ρ·Ĵ_ab, so N keeps Ĵ_hv λ̂^{-1} (same packaging as C2).
  if (pphv_dim != 0 && ph_dim != 0) {
    for (index_t i = pphv_start; i <= pphv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kpp = tbX.GetKet(cf_bra[0]);
      Ket &kvc = tkX.GetKet(cf_bra[1]);
      size_t a = kpp.p, b = kpp.q, c = kvc.p, d = kvc.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &oc = modelspace->GetOrbit(c);
      Orbit &od = modelspace->GetOrbit(d);
      if (oa.cvq != 1 || ob.cvq != 1 || od.cvq != 1)
        continue;
      if (!pphv_is_pp_hv(i))
        continue;
      const int Jab = tbX.J;
      const int Jhv = tkX.J;
      const double ph_dag = ph(Jab - Jhv);
      for (index_t j = ph_start; j <= ph_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        size_t e = cf_ket[0];
        size_t b1 = cf_ket[1];
        if (b1 != c)
          continue;
        Orbit &oe = modelspace->GetOrbit(e);
        double six = AngMom::SixJ((double)Jhv, (double)Jab, (double)lam,
                                   oe.j2 * 0.5, oc.j2 * 0.5, od.j2 * 0.5);
        double nf = (e == d) ? std::sqrt(2.) : 1.;
        double val = -ph_dag * ph(Jab + (od.j2 + oe.j2) / 2) * hat(Jhv) *
                     lamhatinv * six * RdmTB_J(Jab, a, b, e, d) * nf;
        Nkernel(i, j) += val;
        Nkernel(j, i) += val;
      }
    }
  }

  // C2 ppvv-qv leftover 2b×ρ.
  // AMC of X_abij ρ_akij Y_kb prints Ĵ_qv Ĵ_vv λ̂^{-1}. comm122tts leftover
  // at J_vv already has Ĵ_vv^{-1} Ĵ_qv λ̂^{-1}; GetVSEOM_Overlap then multiplies
  // ρ·Ĵ_vv, so the gold is Ĵ_qv λ̂^{-1}. The extra AMC Ĵ_vv is that overlap
  // weight, not a second factor in N. Using Ĵ_qv Ĵ_vv λ̂^{-2} (the C1 leftover-2b
  // rule) only matches when Ĵ_vv=λ̂, i.e. J_vv=λ.
  if (ppvv_dim != 0 && qv_dim != 0) {
    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      auto &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbX = modelspace->GetTwoBodyChannel(cf_bra[2]);
      TwoBodyChannel &tkX = modelspace->GetTwoBodyChannel(cf_bra[3]);
      Ket &kpp = tbX.GetKet(cf_bra[0]);
      Ket &kvv = tkX.GetKet(cf_bra[1]);
      size_t a = kpp.p, b = kpp.q, cv = kvv.p, dv = kvv.q;
      Orbit &oa = modelspace->GetOrbit(a);
      Orbit &ob = modelspace->GetOrbit(b);
      Orbit &oc = modelspace->GetOrbit(cv);
      Orbit &od = modelspace->GetOrbit(dv);
      if (oa.cvq != 1 || oc.cvq != 1 || od.cvq != 1)
        continue;
      const int J0 = tbX.J;
      const int J1 = tkX.J;
      for (index_t j = qv_start; j <= qv_end; j++) {
        auto &cf_ket = eom_confs.at(j);
        size_t q = cf_ket[0];
        size_t k = cf_ket[1];
        if (b != q)
          continue;
        Orbit &ok = modelspace->GetOrbit(k);
        double six = AngMom::SixJ((double)J0, (double)J1, (double)lam,
                                   ok.j2 * 0.5, ob.j2 * 0.5, oa.j2 * 0.5);
        double nf = (a == k) ? std::sqrt(2.) : 1.;
        double val = ph(J0 + (oa.j2 + ob.j2) / 2) * hat(J0) * lamhatinv *
                     six * we1(ob, ok) * RdmTB_J(J1, a, k, cv, dv) * nf;
        Nkernel(i, j) += val;
        Nkernel(j, i) += val;
      }
    }
  }

  std::cout << "ConstructNormMatrix_tensor: λ=" << lam
            << " eom_dims=" << eom_dims
            << " pphv_B1=" << (include_norm_B1 ? 1 : 0)
            << " pphv_C4=" << (include_norm_C4 ? 1 : 0) << std::endl;
}

std::vector<std::tuple<size_t,size_t,size_t,double>> EOM::ThreeBody_Diagram_Entries(size_t a, size_t b, size_t c, size_t d, size_t e,
                          size_t f, size_t g, double j0, double j2) {
  std::vector<std::tuple<size_t,size_t,size_t,double>> result;
  Orbit &oa = modelspace->GetOrbit(a);
  Orbit &ob = modelspace->GetOrbit(b);
  Orbit &oc = modelspace->GetOrbit(c);
  Orbit &od = modelspace->GetOrbit(d);
  Orbit &oe = modelspace->GetOrbit(e);
  Orbit &of = modelspace->GetOrbit(f);
  Orbit &og = modelspace->GetOrbit(g);

  size_t j1_min = abs(od.j2 - oe.j2)/2;
  size_t j1_max = abs(od.j2 + oe.j2)/2;

  size_t j0_min = abs(oa.j2 - ob.j2)/2;
  size_t j0_max = abs(oa.j2 + ob.j2)/2;

  int jabc_max_2j, jabc_min_2j, jdef_max_2j, jdef_min_2j;

  {
    int largest = std::max({oa.j2, ob.j2, oc.j2});
    int sum_others = oa.j2 + ob.j2 + oc.j2 - largest;
    jabc_max_2j = oa.j2 + ob.j2 + oc.j2;
    jabc_min_2j = (largest > sum_others) ? (largest - sum_others) : 1;
  }

  {
    int largest = std::max({od.j2, oe.j2, of.j2});
    int sum_others = od.j2 + oe.j2 + of.j2 - largest;
    jdef_max_2j = od.j2 + oe.j2 + of.j2;
    jdef_min_2j = (largest > sum_others) ? (largest - sum_others) : 1;
  }

  int jmin_2j = std::max(jabc_min_2j, jdef_min_2j);
  int jmax_2j = std::min(jabc_max_2j, jdef_max_2j);
  size_t jmin = (jmin_2j > 0) ? static_cast<size_t>(jmin_2j) : 1;
  size_t jmax = (jmax_2j >= jmin_2j) ? static_cast<size_t>(jmax_2j) : 0;

  size_t j1_len = (j1_max >= j1_min) ? (j1_max - j1_min + 1) : 0;
  size_t j0_len = (j0_max >= j0_min) ? (j0_max - j0_min + 1) : 0;
  arma::mat j1_array(j0_len, j1_len, arma::fill::zeros);

  double norm_fact = 1.;
  if (a == b)
    norm_fact *= sqrt(2.);
  if (g == d)
    norm_fact *= sqrt(2.);
  if (c == g)
    norm_fact *= sqrt(2.);
  if (e == f)
    norm_fact *= sqrt(2.);

  for (auto jtot = jmin; jtot <= jmax; jtot += 2) {
    j1_array.zeros();
    size_t j0_new = j0 - j0_min;

    for (size_t k = 0; k < j1_len; k++) {
      size_t j1 = j1_min + k;
      j1_array(j0_new, k) = sqrt((2. * j0 + 1.) * (2. * j1 + 1.)) * (2 * j2 + 1.) *
            AngMom::SixJ(od.j2 * 0.5, og.j2 * 0.5, j0, oc.j2 * 0.5, jtot * 0.5, j2) *
            AngMom::SixJ(od.j2 * 0.5, oe.j2 * 0.5, j1, of.j2 * 0.5, jtot * 0.5, j2) * norm_fact;
    }

    arma::mat j1_array_bak = j1_array;
    j1_array.zeros();

    if (b == c || a == c) {
      if (b == c) {
        for (size_t k = 0; k < j0_len; k++) {
          size_t j0x = j0_min + k;
          size_t k0 = j0 - j0_min;
          double phase_bc = -pow(-1, ob.j2 * 0.5 + oc.j2 * 0.5 + j0x + j0) * sqrt((2. * j0x + 1.) * (2. * j0 + 1.));
          double angmom_bc = AngMom::SixJ(ob.j2 * 0.5, oa.j2 * 0.5, j0x, oc.j2 * 0.5, jtot * 0.5, j0);
          j1_array.row(k) += j1_array_bak.row(k0) * phase_bc * angmom_bc;
        }
      }
      if (a == c) {
        for (size_t k = 0; k < j0_len; k++) {
          size_t j0x = j0_min + k;
          size_t k0 = j0 - j0_min;
          double phase_ac = sqrt((2. * j0x + 1.) * (2. * j0 + 1.));
          double angmom_ac = AngMom::SixJ(oa.j2 * 0.5, ob.j2 * 0.5, j0x, oc.j2 * 0.5, jtot * 0.5, j0);
          j1_array.row(k) += j1_array_bak.row(k0) * phase_ac * angmom_ac;
        }
      }
    }

    j1_array_bak += j1_array;
    j1_array.zeros();

    if (d == e || d == f) {
      if (d == e) {
        for (size_t k = 0; k < j1_len; k++) {
          size_t j1 = j1_min + k;
          double phase_de = -pow(-1, od.j2 * 0.5 + oe.j2 * 0.5 - static_cast<double>(j1));
          j1_array.col(k) += j1_array_bak.col(k) * phase_de;
        }
      }

      if (d == f) {
        for (size_t k = 0; k < j1_len; k++) {
          size_t j1 = j1_min + k;
          for (size_t m = 0; m < j1_len; m++) {
            size_t j3 = j1_min + m;
            double phase_df = sqrt((2. * j1 + 1.) * (2. * j3 + 1.));
            double angmom_df = AngMom::SixJ(od.j2 * 0.5, oe.j2 * 0.5, j1, of.j2 * 0.5, jtot * 0.5, j3);
            j1_array.col(k) += j1_array_bak.col(m) * phase_df * angmom_df;
          }
        }
      }
    }

    j1_array_bak += j1_array;
    j1_array = j1_array_bak;

    for (size_t r = 0; r < j0_len; r++) {
      for (size_t c1 = 0; c1 < j1_len; c1++) {
        double v = j1_array(r, c1);
        if (std::abs(v) > 1e-14)
          result.emplace_back(j0_min + r, j1_min + c1, static_cast<size_t>(jtot), v);
      }
    }
  }

  return result;
}

double EOM::Norm_abc(size_t p, size_t q, size_t r)
{
  int d_pq = (p == q) ? 1 : 0;
  int d_qr = (q == r) ? 1 : 0;
  int d_pr = (p == r) ? 1 : 0;
  int g_pqr = 1 + d_pq + d_qr + d_pr + 2 * (d_pq * d_qr);
  return std::sqrt(static_cast<double>(g_pqr));
}

double EOM::ThreeBody_Diagram_Entries_Internal(size_t a, size_t b, size_t c,
                                               size_t d, size_t e, size_t f,
                                               size_t g, double j0, double j2)
{
  if (!use_rdm3 || !rdm.ThreeBody.IsAllocated())
    return 0.0;

  double norm_denom = Norm_abc(a, b, c) * Norm_abc(d, e, f);

  double value = 0.0;
  for (const auto &entry : ThreeBody_Diagram_Entries(a, b, c, d, e, f, g, j0, j2)) {
    size_t Jab = std::get<0>(entry);
    size_t Jde = std::get<1>(entry);
    size_t twoJ = std::get<2>(entry);
    double zme = std::get<3>(entry);

    double rdm_me = RdmThreeBody_J(static_cast<int>(Jab), a, b, c,
                                   static_cast<int>(Jde), d, e, f,
                                   static_cast<int>(twoJ));
    if (std::abs(rdm_me) < 1e-14)
      continue;

    value += std::sqrt(double(twoJ) + 1.0) * zme * rdm_me
             / (norm_denom * norm_denom);
  }

  return value;
}

double EOM::ThreeBody_Diagram(size_t a, size_t b, size_t c,
                              size_t d, size_t e, size_t f,
                              size_t g, double j0, double j2)
{
  std::array<size_t, 3> bra = {a, b, c};
  std::array<size_t, 3> ket = {d, e, f};
  std::sort(bra.begin(), bra.end());
  std::sort(ket.begin(), ket.end());

  double convention_factor = (bra == ket) ? 2.0 : 1.0;
  return convention_factor
      * ThreeBody_Diagram_Entries_Internal(a, b, c, d, e, f, g, j0, j2);
}

double EOM::Core_Diagram(size_t a, size_t b, size_t c, size_t d, size_t e,
                         size_t f, double j1, double j2) {
  double val = 0.;
  Orbit &oa = modelspace->GetOrbit(a);
  Orbit &ob = modelspace->GetOrbit(b);
  Orbit &oc = modelspace->GetOrbit(c);
  Orbit &od = modelspace->GetOrbit(d);
  Orbit &oe = modelspace->GetOrbit(e);
  Orbit &of = modelspace->GetOrbit(f);
  size_t jmin = std::max(abs(oa.j2 - ob.j2) / 2, abs(oc.j2 - od.j2) / 2);
  size_t jmax = std::min(abs(oa.j2 + ob.j2) / 2, abs(oc.j2 + od.j2) / 2);
  size_t p_bra = (oa.l + ob.l, 2) % 2;
  size_t p_ket = (oc.l + od.l, 2) % 2;
  if (p_bra != p_ket)
    return (val);
  int itz_bra = (oa.tz2 + ob.tz2) / 2;
  int itz_ket = (oc.tz2 + od.tz2) / 2;
  if (itz_bra != itz_ket)
    return (val);
  double norm_fact = 1.;
  if (a == b)
    norm_fact *= sqrt(2.);
  if (c == d)
    norm_fact *= sqrt(2.);
  // this factor is 2/sqrt{2}, 2 is for permutation, and sqrt2 is from
  // unnormalized to normalized very important
  for (auto jtot = jmin; jtot <= jmax; jtot += 1) {
    double val1 = (2. * j1 + 1.) * (2. * j2 + 1.) *
                  AngMom::NineJ(j1, of.j2 * 0.5, oa.j2 * 0.5, oe.j2 * 0.5, j2,
                                ob.j2 * 0.5, oc.j2 * 0.5, od.j2 * 0.5, jtot);
    double val2 = RdmTB_J((double)jtot, a, b, c, d) *
                  sqrt(2 * jtot + 1.);
    val += norm_fact * val1 * val2 *
           (pow(-1, oa.j2 * 0.5 + ob.j2 * 0.5 + oc.j2 * 0.5 + od.j2 * 0.5));
  }
  return (val);
}

double EOM::Core_Diagram_tensor(size_t a, size_t b, size_t c, size_t d,
                                 size_t e, size_t f, int j1, int j2, int j3,
                                 int j4, int lam) {
  // One orbit pair (e,f). Unreduced scalar Z_abcd^J from
  // ⟨af|ce⟩^{j1 j2 λ} ⟨eb|fd⟩^{j3 j4 λ}. lam is the rank of Q.
  //   −(−1)^{J+j1+j2+j3+j4+ja+jd+je+jf+λ} Ĵ1 Ĵ2 Ĵ3 Ĵ4 ĵ0² λ̂^{-1}
  //   6j{j2 λ j1; jf ja j0} 6j{j3 λ j4; jf jd j0}
  //   9j{jc je j2; jd j3 j0; J jb ja}
  // Times ρ^J_abcd Ĵ and the identical-orbit √2, as in Core_Diagram.
  double val = 0.;
  Orbit &oa = modelspace->GetOrbit(a);
  Orbit &ob = modelspace->GetOrbit(b);
  Orbit &oc = modelspace->GetOrbit(c);
  Orbit &od = modelspace->GetOrbit(d);
  Orbit &oe = modelspace->GetOrbit(e);
  Orbit &of = modelspace->GetOrbit(f);

  if (((oa.l + ob.l) & 1) != ((oc.l + od.l) & 1))
    return val;
  if ((oa.tz2 + ob.tz2) != (oc.tz2 + od.tz2))
    return val;

  const double ja = oa.j2 * 0.5;
  const double jb = ob.j2 * 0.5;
  const double jc = oc.j2 * 0.5;
  const double jd = od.j2 * 0.5;
  const double je = oe.j2 * 0.5;
  const double jf = of.j2 * 0.5;
  const double dlam = (double)lam;
  const double dj1 = (double)j1;
  const double dj2 = (double)j2;
  const double dj3 = (double)j3;
  const double dj4 = (double)j4;
  if (!AngMom::Triangle(dj1, ja, jf) || !AngMom::Triangle(dj2, jc, je) ||
      !AngMom::Triangle(dj3, je, jb) || !AngMom::Triangle(dj4, jf, jd) ||
      !AngMom::Triangle(dj2, dlam, dj1) || !AngMom::Triangle(dj3, dlam, dj4))
    return val;

  double norm_fact = 1.;
  if (a == b)
    norm_fact *= std::sqrt(2.);
  if (c == d)
    norm_fact *= std::sqrt(2.);

  const double hats = std::sqrt(2.0 * j1 + 1.0) * std::sqrt(2.0 * j2 + 1.0) *
                      std::sqrt(2.0 * j3 + 1.0) * std::sqrt(2.0 * j4 + 1.0) /
                      std::sqrt(2.0 * lam + 1.0);
  const int phase_j =
      (oa.j2 + od.j2 + oe.j2 + of.j2) / 2;

  const int Jmin = std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
  const int Jmax = std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
  const int two_j0_min = std::max(std::max(std::abs(2 * j2 - oa.j2),
                                           std::abs(of.j2 - 2 * lam)),
                                  std::abs(2 * j3 - od.j2));
  const int two_j0_max =
      std::min(std::min(2 * j2 + oa.j2, of.j2 + 2 * lam), 2 * j3 + od.j2);

  for (int J = Jmin; J <= Jmax; ++J) {
    if (!AngMom::Triangle(ja, jb, (double)J) ||
        !AngMom::Triangle(jc, jd, (double)J))
      continue;
    double rho = RdmTB_J((double)J, a, b, c, d);
    if (std::abs(rho) < 1e-16)
      continue;
    const double hatJ = std::sqrt(2.0 * J + 1.0);
    const double phaseJ = -(double)AngMom::phase(J + j1 + j2 + j3 + j4 + lam +
                                                 phase_j);
    for (int two_j0 = two_j0_min; two_j0 <= two_j0_max; two_j0 += 2) {
      double j0 = 0.5 * two_j0;
      if (!AngMom::Triangle(dj2, ja, j0) || !AngMom::Triangle(jf, dlam, j0) ||
          !AngMom::Triangle(dj3, jd, j0))
        continue;
      double six1 = AngMom::SixJ(dj2, dlam, dj1, jf, ja, j0);
      double six2 = AngMom::SixJ(dj3, dlam, dj4, jf, jd, j0);
      double nine = AngMom::NineJ(jc, je, dj2, jd, dj3, j0, (double)J, jb, ja);
      val += phaseJ * hats * (2.0 * j0 + 1.0) * six1 * six2 * nine * rho *
             hatJ;
    }
  }
  return val * norm_fact;
}

void EOM::PrintConfigs() {
  // 1-body types: c[0]=a, c[1]=b (orbit indices directly)
  auto print_1b = [&](const std::string &label, index_t start, index_t dim) {
    if (dim == 0) return;
    for (index_t i = 0; i < dim; i++) {
      auto &cfs = eom_confs.at(start + i);
      std::cout << label << " " << (start + i)
                << " " << cfs[0] << " " << cfs[1] << std::endl;
    }
  };
  // 2-body types: c[0]=ibra, c[1]=iket, c[2]=ch_bra, c[3]=ch_ket
  auto print_2b = [&](const std::string &label, index_t start, index_t dim) {
    if (dim == 0) return;
    for (index_t i = 0; i < dim; i++) {
      auto &cfs = eom_confs.at(start + i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cfs[2]);
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cfs[3]);
      Ket &kbra = tbc_bra.GetKet(cfs[0]);
      Ket &kket = tbc_ket.GetKet(cfs[1]);
      std::cout << label << " " << (start + i)
                << " " << kbra.p << " " << kbra.q
                << " " << kket.p << " " << kket.q
                << " " << cfs[2] << " " << cfs[3] << std::endl;
    }
  };
  print_1b("qv",   qv_start,   qv_dim);
  print_1b("ph",   ph_start,   ph_dim);
  print_2b("ppvv", ppvv_start, ppvv_dim);
  print_2b("pphv", pphv_start, pphv_dim);
  print_2b("pphh", pphh_start, pphh_dim);
}

void EOM::ConstructProjectMatrix() {
  Prj_kernel.set_size(eom_dims, eom_dims);
  Prj_kernel.zeros();
  if (eom_dims == 0)
    return;

  // Tensor configs couple different 2b channels (ch_bra ≠ ch_ket). The scalar
  // per-channel blocking is wrong; use one positive-N projector on the full kernel.
  if (eom_tensor_configs) {
    std::vector<int> all_idx;
    all_idx.reserve(eom_dims);
    for (index_t i = 0; i < eom_dims; ++i)
      all_idx.push_back((int)i);
    std::cout << "ConstructProjectMatrix (tensor): full-N positive projector, dim="
              << eom_dims << std::endl;
    block_svd(all_idx);
    return;
  }


  size_t number_channels = modelspace->GetNumberTwoBodyChannels();
  int nconf_pphh_ch[number_channels] = {0};
  int nconf_pphv_ch[number_channels] = {0};
  int nconf_ppvv_ch[number_channels] = {0};

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    index_t ich = eom_confs.at(i)[2];
    nconf_ppvv_ch[ich] += 1;
  }
  for (index_t i = pphh_start; i <= pphh_end; i++) {
    index_t ich = eom_confs.at(i)[2];
    nconf_pphh_ch[ich] += 1;
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    index_t ich = eom_confs.at(i)[2];
    nconf_pphv_ch[ich] += 1;
  }

  // pphh: no coupling to other diagrams, no coupling between channels
  for (index_t ich = 0; ich < number_channels; ich++) {
    if (nconf_pphh_ch[ich] == 0)
      continue;
    std::vector<int> coupled_idx;
    for (index_t i = pphh_start; i <= pphh_end; i++) {
      if (ich == eom_confs.at(i)[2])
        coupled_idx.push_back(i);
    }
    if (coupled_idx.size() == 0)
      continue;
    block_svd(coupled_idx);
  }

  // vqvv coupled with 1-body qv
  std::vector<int> coupled_idx;
  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    index_t ich = eom_confs.at(i)[2];
    size_t ibra = eom_confs.at(i)[0];
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    Ket &dbra = tbc.GetKet(ibra);
    index_t p = dbra.p;
    index_t q = dbra.q;
    Orbit &op = modelspace->GetOrbit(p);
    Orbit &oq = modelspace->GetOrbit(q);
    if (op.cvq != 1)
      continue; // p should be valence
    if (oq.cvq == 1)
      continue; // q should be exclude
    coupled_idx.push_back(i);
  }

  for (index_t i = qv_start; i <= qv_end; i++)
    coupled_idx.push_back(i);

  if (coupled_idx.size() != 0)
    block_svd(coupled_idx);

  // qqvv channels (no coupling between channels)
  for (index_t ich = 0; ich < number_channels; ich++) {

    if (nconf_ppvv_ch[ich] == 0)
      continue;
    coupled_idx.clear();

    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      if (ich != eom_confs.at(i)[2])
        continue;

      size_t ibra = eom_confs.at(i)[0];
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      Ket &dbra = tbc.GetKet(ibra);
      index_t p = dbra.p;
      index_t q = dbra.q;
      Orbit &op = modelspace->GetOrbit(p);
      Orbit &oq = modelspace->GetOrbit(q);
      if (op.cvq == 1)
        continue; // p should be exclude
      if (oq.cvq == 1)
        continue; // q should be exclude
      coupled_idx.push_back(i);
    }

    if (coupled_idx.size() == 0)
      continue;
    block_svd(coupled_idx);
  }

  // pphv coupled with 1-body ph, sorted by hole orbit
  // vvhv and vh
  for (index_t i = 0; i < modelspace->norbits; i++) {
    Orbit &oi = modelspace->GetOrbit(i);
    if (oi.occ != 1)
      continue;
    // find all pphv configs contains oi
    coupled_idx.clear();

    for (index_t j = pphv_start; j <= pphv_end; j++) {
      size_t ibra = eom_confs.at(j)[0];
      size_t iket = eom_confs.at(j)[1];
      size_t ich = eom_confs.at(j)[2];
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      Ket &dket = tbc.GetKet(iket);
      Ket &dbra = tbc.GetKet(ibra);
      index_t p = dbra.p;
      Orbit &op = modelspace->GetOrbit(p);
      index_t q = dbra.q;
      Orbit &oq = modelspace->GetOrbit(q);
      index_t r = dket.p;
      index_t s = dket.q;
      Orbit &os = modelspace->GetOrbit(s);
      if (r != i)
        continue;
      if (op.cvq != 1)
        continue;
      if (oq.cvq != 1)
        continue;
      if (os.cvq != 1)
        continue;
      coupled_idx.push_back(j);
    }
    for (index_t j = ph_start; j <= ph_end; j++) {
      if (eom_confs.at(j)[1] != i)
        continue;
      Orbit &op = modelspace->GetOrbit(eom_confs.at(j)[0]);
      if (op.cvq != 1)
        continue;
      coupled_idx.push_back(j);
    }

    if (coupled_idx.size() == 0)
      continue;
    block_svd(coupled_idx);
  }

  // vphv and qh
  for (index_t i = 0; i < modelspace->norbits; i++) {
    Orbit &oi = modelspace->GetOrbit(i);
    if (oi.occ != 1)
      continue;
    // find all pphv configs contains oi
    coupled_idx.clear();

    for (index_t j = pphv_start; j <= pphv_end; j++) {
      size_t ibra = eom_confs.at(j)[0];
      size_t iket = eom_confs.at(j)[1];
      size_t ich = eom_confs.at(j)[2];
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      Ket &dket = tbc.GetKet(iket);
      Ket &dbra = tbc.GetKet(ibra);
      index_t p = dbra.p;
      Orbit &op = modelspace->GetOrbit(p);
      index_t q = dbra.q;
      Orbit &oq = modelspace->GetOrbit(q);
      index_t r = dket.p;
      index_t s = dket.q;
      Orbit &os = modelspace->GetOrbit(s);
      if (r != i)
        continue;
      if (op.cvq != 1)
        continue;
      if (oq.cvq == 1)
        continue;
      if (os.cvq != 1)
        continue;
      coupled_idx.push_back(j);
    }
    for (index_t j = ph_start; j <= ph_end; j++) {
      if (eom_confs.at(j)[1] != i)
        continue;
      Orbit &op = modelspace->GetOrbit(eom_confs.at(j)[0]);
      if (op.cvq == 1)
        continue;
      coupled_idx.push_back(j);
    }

    if (coupled_idx.size() == 0)
      continue;
    block_svd(coupled_idx);
  }

  // pphv itself, no coupling between channels
  for (index_t ich = 0; ich < number_channels; ich++) {

    if (nconf_pphv_ch[ich] == 0)
      continue;
    coupled_idx.clear();
    for (index_t i = pphv_start; i <= pphv_end; i++) {
      if (ich != eom_confs.at(i)[2])
        continue;
      size_t ibra = eom_confs.at(i)[0];
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      Ket &dbra = tbc.GetKet(ibra);
      index_t p = dbra.p;
      index_t q = dbra.q;
      Orbit &op = modelspace->GetOrbit(p);
      Orbit &oq = modelspace->GetOrbit(q);
      if (op.cvq == 1)
        continue; // p should be exclude
      if (oq.cvq == 1)
        continue; // q should be exclude
      coupled_idx.push_back(i);
    }
    if (coupled_idx.size() == 0)
      continue;

    block_svd(coupled_idx);
  }
}

// True orthogonal projector P = U_r * U_r^T  (idempotent: P*P = P)
// Projects onto the range of the norm matrix, removing the null space.
void EOM::SqrtMat(arma::mat& Amat, size_t n)
{
    arma::mat U, V;
    arma::vec s;
    arma::svd(U, s, V, Amat);

    double s_max = (s.n_elem > 0) ? arma::max(s) : 0.0;

    // If the entire block is zero, the projector is zero.
    if (s_max < 1e-6) {
        Amat.zeros(n, n);
        return;
    }

    // Use a relative threshold to cleanly separate range from null space.
    // Vectors with s_i < 1e-6 * s_max are floating-point noise, not in the range.
    arma::uvec range_idx = arma::find(s >= 1e-6 * s_max);

    Amat.zeros(n, n);
    if (range_idx.n_elem > 0) {
        arma::mat Us = U.cols(range_idx);
        Amat = Us * Us.t();
    }
    return;

}

double EOM::ComputeNorm(Operator &Op1, Operator &Op2) {
  // Flatten Op1 to vector (raw matrix elements, angular momentum factors are in Nkernel)
  arma::vec v1(eom_confs.size());
  v1.fill(0.);

  for (index_t i = qv_start; i <= qv_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    v1(i) = Op1.GetOneBody(p, q);
  }
  for (index_t i = ph_start; i <= ph_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    v1(i) = Op1.GetOneBody(p, q);
  }

  for (index_t i = pphh_start; i <= pphh_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v1(i) = Op1.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v1(i) = Op1.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v1(i) = Op1.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }

  // Flatten Op2 to vector (raw matrix elements, angular momentum factors are in Nkernel)
  arma::vec v2(eom_confs.size());
  v2.fill(0.);

  for (index_t i = qv_start; i <= qv_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    v2(i) = Op2.GetOneBody(p, q);
  }
  for (index_t i = ph_start; i <= ph_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    v2(i) = Op2.GetOneBody(p, q);
  }

  for (index_t i = pphh_start; i <= pphh_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v2(i) = Op2.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v2(i) = Op2.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    v2(i) = Op2.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3], ibra, iket);
  }

  // Compute v1.t() * Nkernel * v2
  arma::vec temp = Nkernel * v2;
  double norm = arma::dot(v1, temp);

  return norm;
}

void EOM::ProjectOprator(Operator &Qin) {
  arma::vec Qflat(eom_confs.size());
  Qflat.fill(0.);

  auto get1 = [&](index_t i) {
    Qflat(i) = Qin.GetOneBody(eom_confs.at(i)[0], eom_confs.at(i)[1]);
  };
  auto get2 = [&](index_t i) {
    Qflat(i) = Qin.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3],
                              eom_confs.at(i)[0], eom_confs.at(i)[1]);
  };
  if (qv_dim > 0)
    for (index_t i = qv_start; i <= qv_end; ++i) get1(i);
  if (ph_dim > 0)
    for (index_t i = ph_start; i <= ph_end; ++i) get1(i);
  if (pphh_dim > 0)
    for (index_t i = pphh_start; i <= pphh_end; ++i) get2(i);
  if (ppvv_dim > 0)
    for (index_t i = ppvv_start; i <= ppvv_end; ++i) get2(i);
  if (pphv_dim > 0)
    for (index_t i = pphv_start; i <= pphv_end; ++i) get2(i);

  Qin *= 0.;

  arma::vec Qout = Prj_kernel * Qflat;

  auto put1 = [&](index_t i) {
    Qin.SetOneBody(eom_confs.at(i)[0], eom_confs.at(i)[1], Qout(i));
  };
  auto put2 = [&](index_t i) {
    Qin.TwoBody.SetTBME(eom_confs.at(i)[2], eom_confs.at(i)[3],
                        eom_confs.at(i)[0], eom_confs.at(i)[1], Qout(i));
  };
  if (qv_dim > 0)
    for (index_t i = qv_start; i <= qv_end; ++i) put1(i);
  if (ph_dim > 0)
    for (index_t i = ph_start; i <= ph_end; ++i) put1(i);
  if (pphh_dim > 0)
    for (index_t i = pphh_start; i <= pphh_end; ++i) put2(i);
  if (ppvv_dim > 0)
    for (index_t i = ppvv_start; i <= ppvv_end; ++i) put2(i);
  if (pphv_dim > 0)
    for (index_t i = pphv_start; i <= pphv_end; ++i) put2(i);
}

void EOM::block_svd(std::vector<int> &coupled_vector) {

  arma::mat Amat;
  int n = coupled_vector.size();
  Amat.set_size(n, n);
  Amat.fill(0.);
  for (index_t i = 0; i < n; i++) {
    index_t ia = coupled_vector.at(i);
    for (index_t j = 0; j < n; j++) {
      index_t ib = coupled_vector.at(j);
      Amat(i, j) = Nkernel(ia, ib);
    }
  }

  SqrtMat(Amat, n);

  for (index_t i = 0; i < n; i++) {
    index_t ia = coupled_vector.at(i);
    for (index_t j = 0; j < n; j++) {
      index_t ib = coupled_vector.at(j);
      Prj_kernel(ia, ib) = Amat(i, j);
    }
  }
}

arma::vec EOM::FlattenOperator(Operator &Op) const
{
  arma::vec v(eom_dims, arma::fill::zeros);
  auto set1 = [&](index_t i) {
    v(i) = Op.GetOneBody(eom_confs.at(i)[0], eom_confs.at(i)[1]);
  };
  auto set2 = [&](index_t i) {
    v(i) = Op.GetTwoBody(eom_confs.at(i)[2], eom_confs.at(i)[3],
                         eom_confs.at(i)[0], eom_confs.at(i)[1]);
  };
  if (qv_dim > 0)
    for (index_t i = qv_start; i <= qv_end; ++i) set1(i);
  if (ph_dim > 0)
    for (index_t i = ph_start; i <= ph_end; ++i) set1(i);
  if (pphh_dim > 0)
    for (index_t i = pphh_start; i <= pphh_end; ++i) set2(i);
  if (ppvv_dim > 0)
    for (index_t i = ppvv_start; i <= ppvv_end; ++i) set2(i);
  if (pphv_dim > 0)
    for (index_t i = pphv_start; i <= pphv_end; ++i) set2(i);
  return v;
}

Operator EOM::ZeroAmplitude() const
{
  Operator op(*modelspace, J2, itz, parity, 2);
  op.SetHermitian();
  op *= 0.0;
  return op;
}

void EOM::UnflattenOperator(Operator &Op, const arma::vec &v) const
{
  Op *= 0.0;
  auto put1 = [&](index_t i) {
    Op.SetOneBody(eom_confs.at(i)[0], eom_confs.at(i)[1], v(i));
  };
  auto put2 = [&](index_t i) {
    Op.TwoBody.SetTBME(eom_confs.at(i)[2], eom_confs.at(i)[3],
                       eom_confs.at(i)[0], eom_confs.at(i)[1], v(i));
  };
  if (qv_dim > 0)
    for (index_t i = qv_start; i <= qv_end; ++i) put1(i);
  if (ph_dim > 0)
    for (index_t i = ph_start; i <= ph_end; ++i) put1(i);
  if (pphh_dim > 0)
    for (index_t i = pphh_start; i <= pphh_end; ++i) put2(i);
  if (ppvv_dim > 0)
    for (index_t i = ppvv_start; i <= ppvv_end; ++i) put2(i);
  if (pphv_dim > 0)
    for (index_t i = pphv_start; i <= pphv_end; ++i) put2(i);
}

int EOM::BuildCanonicalTransform(double eps)
{
  const int N = (int)eom_dims;
  if (N <= 0)
    throw std::runtime_error("BuildCanonicalTransform: empty EOM space");

  // Dense symmetric overlap S ≡ Nkernel
  arma::mat S(N, N, arma::fill::zeros);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j <= i; ++j)
    {
      double sij = Nkernel(i, j);
      S(i, j) = S(j, i) = sij;
    }

  arma::vec s;
  arma::mat U;
  // Symmetric eigendecomposition (S is PSD in exact arithmetic)
  arma::eig_sym(s, U, S);

  double s_max = (s.n_elem > 0) ? arma::max(s) : 0.0;
  double thresh = eps * std::max(s_max, 1.0);
  arma::uvec keep = arma::find(s >= thresh);
  const int M = (int)keep.n_elem;
  std::cout << "Canonical ortho: N=" << N << "  M=" << M
            << "  s_max=" << s_max << "  thresh=" << thresh
            << "  (eps=" << eps << ")" << std::endl;
  if (M == 0)
    throw std::runtime_error("BuildCanonicalTransform: no retained modes");

  // X = U_sub * s_sub^{-1/2}   (N x M)
  arma::mat Usub = U.cols(keep);
  arma::vec ssub = s.elem(keep);
  arma::mat Sinvsqrt = arma::diagmat(1.0 / arma::sqrt(ssub));
  X_canon = Usub * Sinvsqrt;
  have_X_canon = true;

  // Sanity: X^T N X ≈ I
  arma::mat G = X_canon.t() * S * X_canon;
  double max_off = arma::max(arma::max(arma::abs(G - arma::eye(M, M))));
  std::cout << "  max |X^T N X - I| = " << max_off << std::endl;
  return M;
}

EOM::ArnoldiResult
EOM::SolveGEPCanonical(int state_want, double eps, int n_h2_subspace)
{
  using Clock = std::chrono::steady_clock;
  auto t0 = Clock::now();
  auto elapsed = [](Clock::time_point a) {
    return std::chrono::duration<double>(Clock::now() - a).count();
  };

  if (eom_dims == 0)
    throw std::runtime_error("SolveGEPCanonical: call ConstructConfigs first");
  if (Nkernel.n_rows != eom_dims)
    throw std::runtime_error("SolveGEPCanonical: call ConstructNormMatrix first");

  // Full decoupled H (keep vv); eigenvalues are Delta E.
  PrepareHamiltonianForArnoldi();

  const int N = (int)eom_dims;
  const int M = BuildCanonicalTransform(eps);

  // Dense N
  arma::mat S(N, N, arma::fill::zeros);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j <= i; ++j)
    {
      double sij = Nkernel(i, j);
      S(i, j) = S(j, i) = sij;
    }

  // --- Step 1: hermitized H2 in config basis, then H̃₂ = Xᵀ H₂ X ---
  bool save_proj = arnoldi_use_projection;
  arnoldi_use_projection = false; // canonical X replaces Euclidean P

  arma::mat H2mat(N, N, arma::fill::zeros);
  std::cout << "SolveGEPCanonical: building H2 (" << N << " Htc actions)..."
            << std::endl;

  Operator scratch = ZeroAmplitude();
  std::vector<arma::vec> htc_flat(N);
  for (int j = 0; j < N; ++j)
  {
    arma::vec ej(N, arma::fill::zeros);
    ej(j) = 1.0;
    UnflattenOperator(scratch, ej);
    Operator hj = HtcMultiref(Hs, scratch);
    htc_flat[j] = FlattenOperator(hj);
    if ((j + 1) % 50 == 0 || j + 1 == N)
      std::cout << "  Htc column " << j + 1 << " / " << N << std::endl;
  }

  for (int j = 0; j < N; ++j)
  {
    arma::vec Nh = S * htc_flat[j];
    for (int i = 0; i <= j; ++i)
    {
      double hij = Nh(i);
      double hji = arma::dot(htc_flat[i], S.col(j));
      double h2 = 0.5 * (hij + hji);
      H2mat(i, j) = H2mat(j, i) = h2;
    }
  }
  // free Htc cache
  htc_flat.clear();
  htc_flat.shrink_to_fit();

  arma::mat Htilde2 = X_canon.t() * H2mat * X_canon;
  Htilde2 = arma::symmatu(0.5 * (Htilde2 + Htilde2.t()));

  arma::vec evals2;
  arma::mat evecs2;
  arma::eig_sym(evals2, evecs2, Htilde2);
  std::cout << "H2 spectrum (lowest 5):";
  for (int k = 0; k < std::min(5, (int)evals2.n_elem); ++k)
    std::cout << " " << evals2(k);
  std::cout << std::endl;

  arma::vec e;
  arma::mat C; // N x nwant back-transformed eigenvectors
  arma::mat Htilde_final;
  std::string reason = "canonical_gep_h2";

  if (!arnoldi_use_h3)
  {
    const int nwant = std::min(state_want, (int)evals2.n_elem);
    e = evals2.head(nwant);
    C = X_canon * evecs2.cols(0, nwant - 1);
    Htilde_final = Htilde2;
  }
  else
  {
    // --- Step 2: lowest n_sub H2 eigenvectors; build H2+H3 in that subspace ---
    const int n_sub = std::min(n_h2_subspace, (int)evals2.n_elem);
    std::cout << "SolveGEPCanonical: H3 in lowest " << n_sub
              << " H2 eigenvectors (polarization)..." << std::endl;

    // Orthonormal (N-metric) Operators φ_a = X * evecs2.col(a)
    std::vector<Operator> phi;
    phi.reserve(n_sub);
    for (int a = 0; a < n_sub; ++a)
    {
      arma::vec coef = X_canon * evecs2.col(a);
      Operator pa = ZeroAmplitude();
      UnflattenOperator(pa, coef);
      double nn = ComputeNorm(pa, pa);
      if (nn > 0.0)
        pa = pa / std::sqrt(nn);
      phi.push_back(pa);
    }

    std::vector<double> h3_diag(n_sub, 0.0);
    for (int a = 0; a < n_sub; ++a)
    {
      h3_diag[a] = DcomMultiref(Hs, phi[a]);
      if ((a + 1) % 10 == 0 || a + 1 == n_sub)
        std::cout << "  H3 diag " << a + 1 << " / " << n_sub << std::endl;
    }

    // In the (re-normalized) H2 eigenbasis, use E2 on the diagonal for H2.
    // Rebuild off-diagonal H2 from Htc only if needed; for exact H2 evecs
    // of the dense GEP, off-diagonal H2 vanishes in the X-metric. After
    // N-renorm of φ, keep a small Htc rebuild for consistency with Arnoldi.
    arma::mat Hsub(n_sub, n_sub, arma::fill::zeros);
    std::vector<Operator> htc_phi(n_sub);
    for (int a = 0; a < n_sub; ++a)
      htc_phi[a] = HtcMultiref(Hs, phi[a]);

    std::cout << "  a   E2_a        h2_rebuild   H3_aa       nab" << std::endl;
    for (int a = 0; a < std::min(8, n_sub); ++a)
    {
      double nab = ComputeNorm(phi[a], phi[a]);
      double h2r = ComputeNorm(phi[a], htc_phi[a]);
      std::cout << "  " << a << "  " << evals2(a) << "  " << h2r
                << "  " << h3_diag[a] << "  " << nab << std::endl;
    }

    for (int b = 0; b < n_sub; ++b)
    {
      for (int a = 0; a <= b; ++a)
      {
        double h2_ab;
        if (a == b)
        {
          // Prefer the dense-GEP eigenvalue (exact in X basis before renorm)
          h2_ab = evals2(a);
        }
        else
        {
          h2_ab = 0.5 * (ComputeNorm(phi[a], htc_phi[b])
                         + ComputeNorm(phi[b], htc_phi[a]));
        }

        double h3_ab = 0.0;
        if (a == b)
          h3_ab = h3_diag[a];
        else
        {
          Operator psum = phi[a] + phi[b];
          double cross = DcomMultiref(Hs, psum);
          h3_ab = 0.5 * (cross - h3_diag[a] - h3_diag[b]);
        }
        Hsub(a, b) = Hsub(b, a) = h2_ab + h3_ab;
      }
      if ((b + 1) % 10 == 0 || b + 1 == n_sub)
        std::cout << "  H2+H3 column " << b + 1 << " / " << n_sub << std::endl;
    }

    Hsub = arma::symmatu(0.5 * (Hsub + Hsub.t()));

    // First-order PT reference: E2_a + H3_aa (no off-diagonal H3)
    arma::vec fo(n_sub);
    for (int a = 0; a < n_sub; ++a)
      fo(a) = evals2(a) + h3_diag[a];
    arma::uvec fo_ord = arma::sort_index(fo);
    std::cout << "H2+H3 FO-PT (lowest 5):";
    for (int k = 0; k < std::min(5, n_sub); ++k)
      std::cout << " " << fo(fo_ord(k));
    std::cout << std::endl;

    arma::vec evals3;
    arma::mat evecs3;
    arma::eig_sym(evals3, evecs3, Hsub);
    std::cout << "H2+H3 subspace (lowest 5):";
    for (int k = 0; k < std::min(5, (int)evals3.n_elem); ++k)
      std::cout << " " << evals3(k);
    std::cout << std::endl;

    const int nwant = std::min(state_want, (int)evals3.n_elem);
    e = evals3.head(nwant);
    // Back-transform: c = sum_a evecs3(a,k) * φ_a  → flatten
    C.set_size(N, nwant);
    C.zeros();
    for (int k = 0; k < nwant; ++k)
    {
      arma::vec ck(N, arma::fill::zeros);
      for (int a = 0; a < n_sub; ++a)
        ck += evecs3(a, k) * FlattenOperator(phi[a]);
      C.col(k) = ck;
    }
    Htilde_final = Hsub;
    reason = "canonical_gep_h2_then_h3_subspace";
  }

  arnoldi_use_projection = save_proj;

  const int nwant = (int)e.n_elem;
  std::vector<Operator> ritz;
  for (int k = 0; k < nwant; ++k)
  {
    Operator rk = ZeroAmplitude();
    UnflattenOperator(rk, C.col(k));
    double nn = ComputeNorm(rk, rk);
    if (nn > 0.0)
      rk = rk / std::sqrt(nn);
    ritz.push_back(rk);
  }

  ArnoldiResult result;
  result.energies = e;
  result.eigvecs = C;
  result.ritz = ritz;
  result.hall = Htilde_final;
  result.residuals = arma::vec(nwant, arma::fill::zeros);
  result.max_ortho = 0.0;
  result.steps = M;
  result.converged = true;
  result.stop_reason = reason;

  std::cout << "SolveGEPCanonical done in " << elapsed(t0) << " s"
            << "  M=" << M << "  nwant=" << nwant
            << "  reason=" << reason << std::endl;
  const double eref = GetVSEOM_Overlap_multiref(Hs);
  for (int k = 0; k < nwant; ++k)
    std::cout << "  E(" << k << ") = " << e(k)
              << "  absolute=" << e(k) + eref
              << std::endl;
  return result;
}


void EOM::EraseValence(Operator &H) {
  // DISABLED: zeroing vv and shifting by E_val*N is incorrect for MR-EOM.
  // The connected kernel must use the full decoupled H (vv kept); eigenvalues
  // are already Delta E.  This routine is a no-op kept for API compatibility.
  (void)H;
  ClearReferenceEnergyShift();
  std::cerr << "EOM::EraseValence: DISABLED (no-op). "
            << "Use full decoupled H; do not drop vv / E_val*N.\n";
}

void EOM::EraseQspace(Operator &H) {
  // zero out all qspace↔qspace matrix elements, to remove q space
  for (auto &i : H.modelspace->qspace) {
    for (auto &a : H.modelspace->qspace) {
      H.OneBody(a, i) = 0.;
      H.OneBody(i, a) = 0.;
    }
  }

  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

    for (auto &iket : tbc_ket.GetKetIndex_qq()) {
      for (auto &ibra : tbc_bra.GetKetIndex_qq()) {
        H.TwoBody.SetTBME(ch_bra, ch_ket, ibra, iket, 0.);
      }
    }
  }
} 

void EOM::PrepareHamiltonianForArnoldi() {
  force_decouple(Hs);
  // Full decoupled H (vv kept). Eigenvalues of the connected kernel are Delta E;
  // do not EraseValence or shift by E_val*N.
  Hs_full = Hs;
  have_Hs_full = true;
  ClearReferenceEnergyShift();
}

// we enforce the decoupling condition for numerical stability
void EOM::force_decouple(Operator &H) {
  // Enforce decoupling: zero out all ph / qq-vv / qq-vc / qq-vv matrix elements.
  // One-body: kill core↔active and valence↔qspace off-diagonal blocks.
  for (auto &i : H.modelspace->core) {
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      H.OneBody(a, i) = 0.;
      H.OneBody(i, a) = 0.;
    }
  }

  for (auto &i : H.modelspace->valence) {
    for (auto &a : H.modelspace->qspace) {
      H.OneBody(a, i) = 0.;
      H.OneBody(i, a) = 0.;
    }
  }


  // Two body piece only stored half channel, no need to change
  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

    // PPHH
    //
    for (auto &iket : tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes'
                                                // refer to the reference state)
    {
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                       tbc_bra.GetKetIndex_qv())) {
        H.TwoBody.SetTBME(ch_bra, ch_ket, ibra, iket, 0.);
      }
    }

    // PPvq
    //
    for (auto &iket : tbc_ket.GetKetIndex_vc()) // cc means core-core ('holes'
                                                // refer to the reference state)
    {
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                       tbc_bra.GetKetIndex_qv())) {
        H.TwoBody.SetTBME(ch_bra, ch_ket, ibra, iket, 0.);
      }
    }

    // PPvv
    //
    for (auto &iket : tbc_ket.GetKetIndex_vv()) // cc means core-core ('holes'
                                                // refer to the reference state)
    {
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv())) {
        H.TwoBody.SetTBME(ch_bra, ch_ket, ibra, iket, 0.);
      }
    }
  }
}

Operator EOM::GetVSEOM_ladder_single(Operator &H, int herm) {

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of
  // H int herm, 0 for hermit, and 1 for antihermit
  int hZ = H.IsHermitian() ? +1 : -1;

  int herm_phase = 1;
  Operator Hod = 0.0 * H;

  if (herm == 1) {
    Hod.SetHermitian();
    herm_phase = 1;
  }
  if (herm == -1) {
    Hod.SetAntiHermitian();
    herm_phase = -1;
  }

  // One body: ph block (gated by include_ph).
  if (include_ph) {
    for (auto &i : H.modelspace->core) {
      for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
        Hod.SetOneBody(a, i, H.OneBody(a, i));
      }
    }
  }

  // Two body: pphh block (gated by include_pphh).
  if (include_pphh) {
    for (auto &iter : H.TwoBody.MatEl) {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat &H2 = iter.second;

      if (ch_bra == ch_ket) {
        for (auto &iket : tbc_ket.GetKetIndex_cc()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                           tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }
      }

      if (ch_bra != ch_ket) {
        for (auto &iket : tbc_ket.GetKetIndex_cc()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                           tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }

        // conjugate hhpp → pphh with (anti)hermitian phase
        for (auto &ibra : tbc_bra.GetKetIndex_cc()) {
          for (auto &iket :
               VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                           tbc_ket.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket,
                                  H2(ibra, iket) * herm_phase * hZ);
          }
        }
      }
    }
  }

  return Hod;
}

double EOM::GetVSEOM_Overlap_single(Operator &H1, Operator &H2) {
  if (!H1.IsReduced() && !H2.IsReduced()) {
    // scaler-scaler, we use commutator to compute the overlap, which is more efficient and numerically stable than direct summation
    Operator H2d = GetVSEOM_ladder_single(H2, -1);
    Operator nop = H1 * 0.0;
    nop = Commutator::Commutator(H1, H2d);
    return nop.ZeroBody / 2.0;
  }
  double ovlp = 0.;
  // tensor-tensor, we compute the overlap by direct summation of matrix elements,since we dont have tensor-tensor commutator implemented yet. We apply the same phase factor as in GetVSEOM_ladder_single to ensure the consistency.
  int hZ = H2.IsHermitian() ? +1 : -1;
  hZ=1.;
  // One-body contribution: ph and vq blocks
  for (auto &i : H1.modelspace->holes) {
    Orbit &oi = H1.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H1.modelspace->valence, H1.modelspace->qspace)) {
    Orbit &oa = H1.modelspace->GetOrbit(a);
      //double phase_flip= AngMom::phase((oa.j2-oi.j2)/2.)*hZ;
      double phase_flip= AngMom::phase((oa.j2-oi.j2)/2.);
      //ovlp += phase_flip*H1.OneBody(a, i) * H2.OneBody(a, i)/(2.*H1.rank_J+1.);
      ovlp += H1.OneBody(a, i) * H2.OneBody(a, i)/(2.*H1.rank_J+1.);
    }
  }
 

  for (auto &iter : H1.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];

    TwoBodyChannel &tbc_bra = H1.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H1.modelspace->GetTwoBodyChannel(ch_ket);

    arma::mat &Hmat1 = iter.second;
    arma::mat &Hmat2 = H2.TwoBody.GetMatrix(ch_bra, ch_ket);

    // diagonal channel
    if (ch_bra == ch_ket) {
      for (auto &iket :
           tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the
                                     // reference state)
      {
        for (auto &ibra :
             VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                         tbc_bra.GetKetIndex_qv())) {
        double phase_flip= AngMom::phase(tbc_ket.J-tbc_bra.J)*hZ;
          //ovlp += phase_flip*Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
        }
      }
    }

    // off-diagonal channel
    if (ch_bra != ch_ket) {

      for (auto &iket :
           tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the
                                     // reference state)
      {
        for (auto &ibra :
             VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                         tbc_bra.GetKetIndex_qv())) {
        double phase_flip= AngMom::phase(tbc_ket.J-tbc_bra.J)*hZ;
          //ovlp += phase_flip*Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
        }
      }

      for (auto &iket :
           VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                       tbc_ket.GetKetIndex_qv())) {
        for (auto &ibra :
             tbc_bra.GetKetIndex_cc()) // cc means core-core ('holes' refer to
                                       // the reference state)
        {
            // both H1 and H2 times (iph(j1-j2))*hZ
          double phase_flip= AngMom::phase(tbc_ket.J-tbc_bra.J)*hZ;
          //ovlp += phase_flip*Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket)/(2.*H1.rank_J+1);
        }
      }
    }
  }
  //return ovlp*AngMom::phase(H1.rank_J)*hZ;
  return ovlp;
}

Operator EOM::GetVSEOM_ladder_multiref(Operator &H, int herm) {

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of
  // H int herm, 0 for hermit, and 1 for antihermit
  int hZ = H.IsHermitian() ? +1 : -1;
  Operator Hod = 0.0 * H;

  int herm_phase = -1;

  if (herm == 1) {
    Hod.SetHermitian();
    herm_phase = 1;
  }
  if (herm == -1) {
    Hod.SetAntiHermitian();
    herm_phase = -1;
  }

  // One body: ph = (v∪q) ← core
  if (include_ph) {
    for (auto &i : H.modelspace->core) {
      for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
        Hod.SetOneBody(a, i, H.OneBody(a, i));
      }
    }
  }

  // Two body: pphh = (qq∪qv∪vv) ← cc
  if (include_pphh) {
    for (auto &iter : H.TwoBody.MatEl) {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat &H2 = iter.second;

      if (ch_bra == ch_ket) {
        for (auto &iket : tbc_ket.GetKetIndex_cc()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                           tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }
      }

      if (ch_bra != ch_ket) {
        for (auto &iket : tbc_ket.GetKetIndex_cc()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                           tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }

        for (auto &ibra : tbc_bra.GetKetIndex_cc()) {
          for (auto &iket :
               VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                           tbc_ket.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket,
                                  H2(ibra, iket) * hZ * herm_phase);
          }
        }
      }
    }
  }

  // One body: qv = q ← v
  if (include_qv) {
    for (auto &i : H.modelspace->valence) {
      for (auto &a : H.modelspace->qspace) {
        Hod.SetOneBody(a, i, H.OneBody(a, i));
      }
    }
  }

  // Two body: pphv (vc) and ppvv (vv)
  if (include_pphv || include_ppvv) {
    for (auto &iter : H.TwoBody.MatEl) {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat &H2 = iter.second;

      if (include_pphv) {
        for (auto &iket : tbc_ket.GetKetIndex_vc()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                           tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }
      }

      if (include_ppvv) {
        for (auto &iket : tbc_ket.GetKetIndex_vv()) {
          for (auto &ibra :
               VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv())) {
            Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
          }
        }
      }
    }
  }

  return Hod;
}

double EOM::GetVSEOM_Overlap_multiref(Operator &H) {
  double ovlp  = 0;
  double ovlp1 = 0;
  double ovlp2 = 0;
  double ovlp3 = 0;

  ovlp += H.ZeroBody;
  
  for (auto &i : H.modelspace->valence) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &j : H.modelspace->valence) {
      ovlp1 += H.OneBody(i, j) * RdmOB(i, j) * sqrt((oi.j2 + 1));
    }
  }

  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
    arma::mat &H2 = iter.second;

    for (auto &iket : tbc_ket.GetKetIndex_vv())
    {

      Ket &dket = tbc_ket.GetKet(iket);

      for (auto &ibra : tbc_bra.GetKetIndex_vv()) {
        Ket &dbra = tbc_bra.GetKet(ibra);
        if (tbc_bra.J != tbc_ket.J)
          continue;

        double rdm2 = RdmTB_J(tbc_bra.J, dbra.p, dbra.q, dket.p, dket.q);
        // if (std::abs(rdm2) > 1e-12) {
        //   std::cout << "H2: " << H2(ibra, iket)
        //             << " RDM2: " << rdm2 << std::endl;
        // }

        ovlp2 += H2(ibra, iket) * rdm2 * std::sqrt(2.0 * tbc_bra.J + 1.0);
            
      }
    }
  }

  if (use_rdm3 && rdm.ThreeBody.IsAllocated() && H.ThreeBody.IsAllocated() ) {
    for (auto &it : H.ThreeBody.Get_ch_start()) {
      size_t ch_bra = it.first.ch_bra;
      size_t ch_ket = it.first.ch_ket;
      ThreeBodyChannel &tbc_bra = H.modelspace->GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel &tbc_ket = H.modelspace->GetThreeBodyChannel(ch_ket);
      int twoJ = tbc_bra.twoJ;
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      for (size_t ibra = 0; ibra < nbras; ++ibra) {
        Ket3 &bra = tbc_bra.GetKet(ibra);

        for (size_t iket = 0; iket < nkets; ++iket) {
          Ket3 &ket = tbc_ket.GetKet(iket);
          double h3 = H.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket);
          if (std::abs(h3) < 1e-12)
            continue;

          double r3 = RdmThreeBody_J(bra.Jpq, bra.p, bra.q, bra.r,
                                     ket.Jpq, ket.p, ket.q, ket.r, twoJ);
          if (std::abs(r3) < 1e-12)
            continue;

          double norm_denom = Norm_abc(bra.p, bra.q, bra.r)
                            * Norm_abc(ket.p, ket.q, ket.r);
          ovlp3 += h3 * r3 * std::sqrt(double(twoJ) + 1.0)
                  / (norm_denom * norm_denom);
        }
      }
    }
  }
  //  std::cout << "zero-body contribution to norm: " << ovlp << std::endl;
  //  std::cout << "one-body contribution to norm: " << ovlp1 << std::endl;
  //  std::cout << "two-body contribution to norm: " << ovlp2 << std::endl;
  //  std::cout << "three-body contribution to norm: " << ovlp3 << std::endl;

  return (ovlp + ovlp1 + ovlp2 + ovlp3);
}

// ============================================================
//  Lanczos / Arnoldi helpers  (translated from run/lanczos.py)
// ============================================================

// -----------  norm helpers  ---------------------------------

/// Thin wrapper around GetVSEOM_Overlap_single.
double EOM::NormSingle(Operator &T1, Operator &T2)
{
  Operator T2d = GetVSEOM_ladder_single(T2, -1);
  Operator nop = T1 * 0.0;
  nop = Commutator::Commutator(T1, T2d);
  return nop.ZeroBody / 2.0;
}

/// <T1|T2> using the multiref metric:
///   T2d = ladder(T2, -1)   (anti-Hermitian)
///   nop2 = [T1, T2d]       (IMSRG(2) only → 0b+1b+2b)
///   If use_rdm3: t3 = [T1, T2d]_3 via comm223ss (2b×2b → 3b),
///   then contract t3 with valence ρ₃.
///   result = (⟨nop2⟩_ρ + ⟨t3⟩_ρ) / 2
double EOM::NormMultiref(Operator &T1, Operator &T2)
{
  if (T1.GetJRank() != 0 or T2.GetJRank() != 0)
    return NormMultiref_tensor(T1, T2);

  Operator T2d = T1 * 0.0;
  T2d = GetVSEOM_ladder_multiref(T2, -1);
  Operator nop = T1 * 0.0;
  nop.SetHermitian();
  nop = Commutator::Commutator(T1, T2d);
  double rst = GetVSEOM_Overlap_multiref(nop);

  // ρ₃ piece: [χ₂, χ₂]_3 is not produced by the default IMSRG(2) Commutator.
  if (use_rdm3 && rdm.ThreeBody.IsAllocated())
  {
    Operator t3(*modelspace, 0, 0, 0, 3);
    t3.ThreeBody.SetMode("pn");
    t3 *= 0.0;
    // T1 Hermitian, T2d anti-Hermitian → [T1,T2d] Hermitian
    t3.SetHermitian();
    Commutator::comm223ss(T1, T2d, t3);
    rst += GetVSEOM_Overlap_multiref(t3);
  }

  return rst / 2.0;
}

/// Tensor analogue of NormMultiref: leftover [χ_λ, χ^{-}_λ] is a scalar 0+
/// operator. Production Commutator() exits on T×T, so assemble it from TTS.
double EOM::NormMultiref_tensor(Operator &T1, Operator &T2)
{
  Operator T2d = T2 * 0.0;
  T2d = GetVSEOM_ladder_multiref(T2, -1);

  Operator nop(*modelspace, 0, 0, 0, 2);
  nop.SetHermitian();
  nop *= 0.0;

  Commutator::comm110tts(T1, T2d, nop);
  Commutator::comm220tts(T1, T2d, nop);
  Commutator::comm111tts(T1, T2d, nop);
  Commutator::comm121tts(T1, T2d, nop);
  Commutator::comm221tts(T1, T2d, nop);
  Commutator::comm122tts(T1, T2d, nop);
  Commutator::comm222_pp_hhtts(T1, T2d, nop);
  Commutator::comm222_phtts(T1, T2d, nop);

  return GetVSEOM_Overlap_multiref(nop) / 2.0;
}

/// 3-body contribution to <t1|t2>:
///   t3  = comm223ss(t2, haml)
///   nop = comm232ss(t1,t3) + comm231ss(t1,t3) + comm132ss(t1,t3)
///   result = GetVSEOM_Overlap_multiref(nop) / 2
double EOM::Norm3Multiref(Operator &t1, Operator &t2, Operator &haml)
{
  Operator t3(*modelspace, 0, 0, 0, 3);
  t3.ThreeBody.SetMode("pn");
  t3 *= 0.0;
  t3.SetAntiHermitian();

  Operator nop = t1 * 0.0;
  nop.SetHermitian();

  Commutator::comm223ss(haml, t2, t3);
  Commutator::comm232ss(t1, t3, nop);
  Commutator::comm231ss(t1, t3, nop);
  Commutator::comm132ss(t1, t3, nop);

  return GetVSEOM_Overlap_multiref(nop) / 2.0;
}

// -----------  H * v helpers  --------------------------------

/// Single-reference action: [Haml, chi], then apply ladder (herm=0).
Operator EOM::HtcSingle(Operator &haml, Operator &chi)
{
 
  Operator ht_plus = Commutator::Commutator(haml, chi);
  return GetVSEOM_ladder_single(ht_plus, 1);
}

/// Multiref action: [Haml, chi], then apply ladder (herm=1) and project.
Operator EOM::HtcMultiref(Operator &haml, Operator &chi)
{
  Operator ht_minus = Commutator::Commutator(haml, chi);
  Operator heom = GetVSEOM_ladder_multiref(ht_minus, 1);
  if (arnoldi_use_projection)
    ProjectOprator(heom);
  return heom;
}

// -----------  double commutator diagonal  -------------------

/// Diagonal double-commutator contribution H3:
///   [χ, [χ, H]_3] folded to 1b+2b, then ⟨·⟩_ρ / 2.
/// Scalar: FactorizedDoubleCommutator 223_231/232/132.
/// Tensor χ (λ≠0): FactorizedDoubleCommutator_eths 223_231_st / 223_232 /
///   223_132_tts (leftover is a scalar 0+ operator).
double EOM::DcomMultiref(Operator &haml, Operator &chi_in)
{
  using Clock = std::chrono::steady_clock;
  auto elapsed = [](Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
  };

  auto timer = Clock::now();
  Operator chi = GetVSEOM_ladder_multiref(chi_in, -1);
  if (arnoldi_print_timing)
    dcom_time_ladder += elapsed(timer);

  if (IsTensorEOM()) {
    Operator opa(*modelspace, 0, 0, 0, 2);
    opa.SetHermitian();
    opa *= 0.0;

    namespace eth = Commutator::FactorizedDoubleCommutator_eths;
    eth::SetUse_1b_Intermediates(true);
    eth::SetUse_2b_Intermediates(true);

    timer = Clock::now();
    eth::comm223_231_st(chi, haml, opa);
    if (arnoldi_print_timing)
      dcom_time_223_231 += elapsed(timer);

    timer = Clock::now();
    eth::comm223_232(chi, haml, opa);
    if (arnoldi_print_timing)
      dcom_time_223_232 += elapsed(timer);

    timer = Clock::now();
    eth::comm223_132_tts(chi, haml, opa);
    if (arnoldi_print_timing)
      dcom_time_223_132 += elapsed(timer);

    timer = Clock::now();
    double rst = GetVSEOM_Overlap_multiref(opa);
    if (arnoldi_print_timing)
    {
      dcom_time_overlap += elapsed(timer);
      ++dcom_profile_calls;
    }
    return rst / 2.0;
  }

  Operator opa = chi_in * 0.0;
  opa.SetHermitian();

  Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
  Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(true);

  timer = Clock::now();
  Commutator::FactorizedDoubleCommutator::comm223_231(chi, haml, opa);
  if (arnoldi_print_timing)
    dcom_time_223_231 += elapsed(timer);

  timer = Clock::now();
  Commutator::FactorizedDoubleCommutator::comm223_232(chi, haml, opa);
  if (arnoldi_print_timing)
    dcom_time_223_232 += elapsed(timer);

  timer = Clock::now();
  Commutator::FactorizedDoubleCommutator::comm223_132(chi, haml, opa);
  if (arnoldi_print_timing)
    dcom_time_223_132 += elapsed(timer);

  timer = Clock::now();
  double rst = GetVSEOM_Overlap_multiref(opa);
  if (arnoldi_print_timing)
  {
    dcom_time_overlap += elapsed(timer);
    ++dcom_profile_calls;
  }
  return rst / 2.0;
}

// ============================================================
//  LanczosSolve
//  Translated from lanczos_proc() in run/lanczos.py
// ============================================================

std::pair<arma::vec, std::vector<Operator>>
EOM::LanczosSolve(Operator &vi, int max_iter, int state_want)
{
  std::vector<Operator> lanczos_vector;
  arma::mat hall(max_iter, max_iter, arma::fill::zeros);

  // normalize initial vector
  double nn = GetVSEOM_Overlap_single(vi, vi);
  Operator v0 = vi / std::sqrt(nn);
  lanczos_vector.push_back(v0);

  double norm_e_old = -1000.0;
  double norm_e_new = -1000.0;
  arma::vec e(state_want, arma::fill::zeros);
  int j_final = 0;
  double bj = 0.0;
  bool converged = false;

  for (int j = 0; j < max_iter; ++j)
  {
    j_final = j;
    Operator w = HtcSingle(Hs, lanczos_vector[j]);

    double ai = GetVSEOM_Overlap_single(w, lanczos_vector[j]);
    hall(j, j) = ai;

    if (j > 0)
      w = w - ai * lanczos_vector[j] - bj * lanczos_vector[j - 1];
    else
      w = w - ai * lanczos_vector[j];

    // Full double-pass reorthogonalization against ALL previous Lanczos vectors.
    // The three-term recurrence only subtracts j and j-1, but in finite-precision
    // arithmetic orthogonality to earlier vectors drifts away.  Once any converged
    // Ritz value has lost orthogonality the corresponding eigenvalue is reintroduced
    // as a "ghost" — a duplicate copy in the spectrum.  Two passes of Classical
    // Gram-Schmidt eliminate this completely at O(j) extra inner-product evaluations
    // per step.
    for (int pass = 0; pass < 2; ++pass)
    {
      for (int i = 0; i <= j; ++i)
      {
        double cij = GetVSEOM_Overlap_single(w, lanczos_vector[i]);
        w = w - cij * lanczos_vector[i];
      }
    }

    double nm = GetVSEOM_Overlap_single(w, w);
    bj = std::sqrt(nm);

    if (j < max_iter - 1)
    {
      hall(j, j + 1) = bj;
      hall(j + 1, j) = bj;
    }

    if (bj < 1e-10)
    {
      std::cout << "lanczos: exact breakdown at step " << j + 1 << ", stopping" << std::endl;
      break;
    }
    if (bj < 0.01)
      std::cout << "lanczos: bj is small (" << bj << "), continuing" << std::endl;

    lanczos_vector.push_back(w / bj);

    if ((int)lanczos_vector.size() > state_want && j % 4 == 0)
    {
      int dim = std::min((int)lanczos_vector.size() - 1, max_iter); // subspace built so far
      arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, sub);

      int nshow = std::min(state_want, (int)eigval.n_elem);
      std::cout << "lanczos energy @ step " << j << ": ";
      for (int k = 0; k < nshow; ++k)
        std::cout << eigval(k) << " ";
      std::cout << std::endl;

      e = eigval.head(nshow);

      norm_e_new = 0.0;
      for (int k = 0; k < (int)e.n_elem; ++k)
        norm_e_new += e(k) * e(k);

      if (std::abs(norm_e_new - norm_e_old) < 0.01)
      {
        std::cout << "lanczos: energy converged" << std::endl;
        converged = true;
        break;
      }
      norm_e_old = norm_e_new;
    }
  }

  // Final eigensolution only when convergence was not already reached during
  // the iteration — avoids spurious low eigenvalues from extra basis vectors
  // added after the converged subspace.
  // Bug fix: use size()-1 not size(): the last pushed vector (w/bj) is the
  // boundary vector whose diagonal element a_{j+1} was never computed, so
  // including it would add a spurious zero row/col to the Lanczos matrix.
  if (!converged)
  {
    int dim = std::min((int)lanczos_vector.size() - 1, max_iter);
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    arma::vec eigval_f;
    arma::mat eigvec_f;
    arma::eig_sym(eigval_f, eigvec_f, sub);
    int nshow = std::min(state_want, (int)eigval_f.n_elem);
    e = eigval_f.head(nshow);
  }

  std::cout << "Lanczos finished with " << j_final + 1 << " steps" << std::endl;
  for (int k = 0; k < (int)e.n_elem; ++k)
    std::cout << "E(" << k << ") = " << e(k) << std::endl;

  return {e, lanczos_vector};
}

// ============================================================
//  ArnoldiSolve
//  Translated from arnoldi_proc() in run/lanczos.py
// ============================================================

// ------------------------------------------------------------
// ExpectationValue / ExpectationValueFull
// ------------------------------------------------------------
// Both use the Arnoldi H = H1 + H2 split on an S^(+/-) Ritz vector:
//
//   H1 : <S | Htc(H,S) >   (2b-truncated matvec, then metric)
//   H2 : DcomMultiref(H,S) (induced 3b double commutator; already /2
//                           because S is (anti)Hermitian vs chi)
//
// On the full decoupled H (vv kept), the Rayleigh quotient is Delta E.
// ExpectationValue and ExpectationValueFull are therefore equivalent
// when Hs_full has been set; ExpectationValue uses working Hs.
//
double EOM::ExpectationValue(Operator &Psi)
{
  double nrm = ComputeNorm(Psi, Psi);
  if (std::abs(nrm) < 1e-30) return 0.0;

  Operator h1psi = HtcMultiref(Hs, Psi);
  double h1exp   = ComputeNorm(Psi, h1psi);
  double h2exp   = DcomMultiref(Hs, Psi);

  return (h1exp + h2exp) / nrm;
}

double EOM::ExpectationValueFull(Operator &Psi, bool use_rdm_norm)
{
  if (!have_Hs_full)
    throw std::runtime_error(
        "EOM::ExpectationValueFull: Hs_full not set; call "
        "PrepareHamiltonianForArnoldi() or RunMR first");

  double nrm = 0.0;
  double h1exp = 0.0;
  Operator h1psi = HtcMultiref(Hs_full, Psi);
  if (use_rdm_norm)
  {
    nrm   = NormMultiref(Psi, Psi);
    h1exp = NormMultiref(Psi, h1psi);
  }
  else
  {
    nrm   = ComputeNorm(Psi, Psi);
    h1exp = ComputeNorm(Psi, h1psi);
  }
  if (std::abs(nrm) < 1e-30) return 0.0;

  // DcomMultiref already returns <[S^-,[H,S]]>_{3b}/2 (S vs chi factor).
  double h2exp = DcomMultiref(Hs_full, Psi);
  return (h1exp + h2exp) / nrm;
}

EOM::ArnoldiResult
EOM::ArnoldiSolve(Operator &vi, int max_iter, int state_want)
{
  using Clock = std::chrono::steady_clock;
  auto elapsed = [](Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
  };
  auto total_start = Clock::now();
  dcom_time_ladder = 0.0;
  dcom_time_223_231 = 0.0;
  dcom_time_223_232 = 0.0;
  dcom_time_223_132 = 0.0;
  dcom_time_overlap = 0.0;
  dcom_profile_calls = 0;
  double time_h1_action = 0.0;
  double time_hall_h1 = 0.0;
  double time_hall_h2 = 0.0;
  double time_hall_norm = 0.0;
  double time_gs = 0.0;
  double time_h2_diag = 0.0;
  double time_eigensolve = 0.0;
  double time_ritz = 0.0;
  int h2_cross_calls = 0;
  int h2_diag_calls = 0;

  const double bj_tol   = 1e-4;
  const double null_tol = 1e-4;
  const int    min_iter = state_want + 1;
  const int    stable_window = std::max(1, arnoldi_stable_window);

  std::vector<Operator> lanczos_vector;
  std::vector<Operator> h1v_cache;
  std::vector<double>   h2_diag;
  arma::mat hall(max_iter, max_iter, arma::fill::zeros);

  double nn = ComputeNorm(vi, vi);
  Operator v0 = vi / std::sqrt(nn);
  lanczos_vector.push_back(v0);
  auto timer = Clock::now();
  h2_diag.push_back(arnoldi_use_h3 ? DcomMultiref(Hs, v0) : 0.0);
  time_h2_diag += elapsed(timer);
  ++h2_diag_calls;
  const double cn0 = ComputeNorm(v0, v0);

  arma::vec e(state_want, arma::fill::zeros);
  arma::mat vs(1, 1, arma::fill::zeros);
  arma::vec residuals(state_want, arma::fill::zeros);
  std::vector<arma::vec> energy_history;
  double max_ortho_seen = 0.0;
  bool   converged_flag = false;
  std::string stop_reason = "max_iter";
  int       j_final      = 0;

  auto eigensolve_sub = [&hall](int dim, arma::vec &ev, arma::mat &evec)
  {
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    sub = arma::symmatu(0.5 * (sub + sub.t()));
    arma::eig_sym(ev, evec, sub);
  };

  // Physical root pick on full-H Hall (Delta E). Soft E < soft_floor
  // roots of the truncated kernel are discarded when prefer_positive.
  auto pick_physical_roots =
      [this, state_want](const arma::vec &eigval, const arma::mat &eigvec,
                         int j_col, double beta,
                         arma::vec &e_out, arma::mat &vs_out,
                         arma::vec &resid_out)
  {
    std::vector<int> idx;
    idx.reserve((size_t)state_want);
    if (arnoldi_prefer_positive)
    {
      for (int k = 0; k < (int)eigval.n_elem; ++k)
        if (eigval(k) >= arnoldi_soft_floor)
          idx.push_back(k);
      if ((int)idx.size() < state_want)
      {
        for (int k = 0; k < (int)eigval.n_elem && (int)idx.size() < state_want; ++k)
          if (eigval(k) < arnoldi_soft_floor)
            idx.push_back(k);
      }
    }
    else
    {
      for (int k = 0; k < std::min(state_want, (int)eigval.n_elem); ++k)
        idx.push_back(k);
    }
    if ((int)idx.size() > state_want)
      idx.resize((size_t)state_want);

    const int nwant = (int)idx.size();
    e_out.set_size(nwant);
    resid_out.set_size(nwant);
    vs_out = eigvec; // keep full subspace; selection is via column indices below
    // Compact selected eigenvectors into the first nwant columns of vs_out
    arma::mat vs_sel(eigvec.n_rows, nwant, arma::fill::zeros);
    for (int i = 0; i < nwant; ++i)
    {
      const int k = idx[(size_t)i];
      e_out(i) = eigval(k);
      vs_sel.col(i) = eigvec.col(k);
      if (j_col >= 0 && j_col < (int)eigvec.n_rows)
        resid_out(i) = beta * std::abs(eigvec(j_col, k));
      else
        resid_out(i) = 0.0;
    }
    vs_out = std::move(vs_sel);
  };

  auto report_expectation_check =
      [this, &hall, &lanczos_vector](const arma::vec &ev, const arma::mat &evec,
                              int dim, int nshow, const std::string &tag)
  {
    std::cout << tag << std::endl;
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    sub = arma::symmatu(0.5 * (sub + sub.t()));
    for (int k = 0; k < nshow; ++k)
    {
      Operator ritz = lanczos_vector[0] * 0.0;
      for (int m = 0; m < dim; ++m)
        ritz = ritz + evec(m, k) * lanczos_vector[m];

      arma::vec coeff = evec.col(k).head(dim);
      double hall_expect = arma::as_scalar(coeff.t() * sub * coeff);
      double nrm = ComputeNorm(ritz, ritz);
      Operator h1ritz = HtcMultiref(Hs, ritz);
      double h1exp = ComputeNorm(ritz, h1ritz);
      double h2exp = arnoldi_use_h3 ? DcomMultiref(Hs, ritz) : 0.0;
      double direct_num = h1exp + h2exp;
      double direct_quotient = (std::abs(nrm) < 1e-30) ? 0.0 : direct_num / nrm;
      std::cout << "  state " << k
                << ": diag=" << ev(k)
                << "  hall=<c|H|c>=" << hall_expect
                << "  norm=" << nrm
                << "  direct_num=" << direct_num
                << "  direct_E=" << direct_quotient
                << "  delta_num=" << (direct_num - hall_expect)
                << "  delta_E=" << (direct_quotient - ev(k)) << std::endl;
    }
  };

  auto write_ritz_vectors =
      [&lanczos_vector](const arma::mat &evec, int dim, int nshow,
                        const std::string &prefix)
  {
    ReadWrite rw;
    for (int k = 0; k < nshow; ++k)
    {
      Operator ritz = lanczos_vector[0] * 0.0;
      for (int m = 0; m < dim; ++m)
        ritz = ritz + evec(m, k) * lanczos_vector[m];
      rw.WriteOperator(ritz, prefix + std::to_string(k) + ".op");
    }
  };

  auto gram_offdiag_max = [this, &lanczos_vector]() -> double
  {
    double mx = 0.0;
    const int n = (int)lanczos_vector.size();
    for (int i = 0; i < n; ++i)
      for (int jj = i + 1; jj < n; ++jj)
        mx = std::max(mx, std::abs(ComputeNorm(lanczos_vector[i], lanczos_vector[jj])));
    return mx;
  };

  std::cout << "ArnoldiSolve: full H (vv kept), Hall = (Htc_ab+Htc_ba)/2"
            << (arnoldi_use_h3 ? " + Dcom(H3)" : " (H3 off; PT after)")
            << (arnoldi_prefer_positive
                    ? ("  prefer ΔE>=" + std::to_string(arnoldi_soft_floor))
                    : "  (all roots)")
            << "  (no E_val shift)" << std::endl;

  for (int j = 0; j < max_iter - 1; ++j)
  {
    j_final = j;

    timer = Clock::now();
    Operator h1v_j = HtcMultiref(Hs, lanczos_vector[j]);
    time_h1_action += elapsed(timer);
    h1v_cache.push_back(h1v_j);

    for (int i = 0; i <= j; ++i)
    {
      double h1_sym;
      double h2_sym;
      if (i == j)
      {
        timer = Clock::now();
        h1_sym = ComputeNorm(lanczos_vector[j], h1v_j);
        time_hall_h1 += elapsed(timer);
        h2_sym = h2_diag[j];
      }
      else
      {
        timer = Clock::now();
        double h1ij = ComputeNorm(lanczos_vector[i], h1v_j);
        double h1ji = ComputeNorm(lanczos_vector[j], h1v_cache[i]);
        h1_sym      = 0.5 * (h1ij + h1ji);
        time_hall_h1 += elapsed(timer);

        if (arnoldi_use_h3)
        {
          Operator v_sum  = lanczos_vector[i] + lanczos_vector[j];
          timer = Clock::now();
          double h2_cross = DcomMultiref(Hs, v_sum);
          time_hall_h2 += elapsed(timer);
          ++h2_cross_calls;
          h2_sym          = 0.5 * (h2_cross - h2_diag[i] - h2_diag[j]);
        }
        else
        {
          h2_sym = 0.0;
        }
      }

      timer = Clock::now();
      hall(i, j) = hall(j, i) = h1_sym + h2_sym;
      time_hall_norm += elapsed(timer);
    }

    timer = Clock::now();
    Operator w = h1v_j * 1.0;
    if (arnoldi_use_projection)
      ProjectOprator(w);
    auto gs_pass = [&]() {
      for (int i = 0; i <= j; ++i)
      {
        double cij = ComputeNorm(lanczos_vector[i], w);
        w = w - cij * lanczos_vector[i];
      }
      if (arnoldi_use_projection)
        ProjectOprator(w);
    };
    gs_pass();
    gs_pass();

    double bj = ComputeNorm(w, w);
    time_gs += elapsed(timer);
    double bj_kernel = bj;

    if (std::abs(bj_kernel) < null_tol * cn0)
    {
      std::cout << "arnoldi: null vector (breakdown) at step " << j + 1
                << " (ComputeNorm=" << bj_kernel << "), stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      pick_physical_roots(ev, evec, -1, 0.0, e, vs, residuals);
      stop_reason = "null_breakdown";
      break;
    }
    if (std::abs(bj) < bj_tol)
    {
      std::cout << "arnoldi: exact breakdown at step " << j + 1 << ", stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      pick_physical_roots(ev, evec, -1, 0.0, e, vs, residuals);
      stop_reason = "exact_breakdown";
      break;
    }
    if (bj < 0.0)
    {
      std::cout << "arnoldi: bj=" << bj << " < 0 at step " << j + 1
                << " (indefinite metric), stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      pick_physical_roots(ev, evec, -1, 0.0, e, vs, residuals);
      stop_reason = "indefinite_metric";
      break;
    }

    // Hall/Ritz in span(v0..vj); residual ~ beta * |s_{j,k}|.
    const double beta = std::sqrt(bj);
    Operator new_vec = w / beta;

    if (arnoldi_monitor_ortho)
    {
      double max_ov = 0.0;
      for (int i = 0; i <= j; ++i)
        max_ov = std::max(max_ov, std::abs(ComputeNorm(lanczos_vector[i], new_vec)));
      if (max_ov > arnoldi_ortho_warn)
      {
        for (int pass = 0; pass < 2; ++pass)
        {
          for (int i = 0; i <= j; ++i)
          {
            double cij = ComputeNorm(lanczos_vector[i], new_vec);
            new_vec = new_vec - cij * lanczos_vector[i];
          }
          if (arnoldi_use_projection)
            ProjectOprator(new_vec);
        }
        double nn_new = ComputeNorm(new_vec, new_vec);
        if (nn_new <= 0.0)
        {
          std::cout << "arnoldi: reortho produced non-positive norm at step "
                    << j + 1 << ", stopping." << std::endl;
          stop_reason = "reortho_breakdown";
          arma::vec ev; arma::mat evec;
          eigensolve_sub(j + 1, ev, evec);
          pick_physical_roots(ev, evec, -1, 0.0, e, vs, residuals);
          break;
        }
        new_vec = new_vec / std::sqrt(nn_new);
        max_ov = 0.0;
        for (int i = 0; i <= j; ++i)
          max_ov = std::max(max_ov, std::abs(ComputeNorm(lanczos_vector[i], new_vec)));
      }
      max_ortho_seen = std::max(max_ortho_seen, max_ov);
      if (max_ov > arnoldi_ortho_warn)
        std::cout << "arnoldi: ortho warn at step " << j + 1
                  << "  max|<v_i|v_new>|=" << max_ov << std::endl;
      if (max_ov > arnoldi_ortho_fail)
      {
        std::cout << "arnoldi: orthogonality lost at step " << j + 1
                  << "  max|<v_i|v_new>|=" << max_ov
                  << " > fail tol " << arnoldi_ortho_fail
                  << " — stopping (Ritz unreliable)." << std::endl;
        stop_reason = "ortho_lost";
      }
    }

    lanczos_vector.push_back(new_vec);
    timer = Clock::now();
    h2_diag.push_back(arnoldi_use_h3 ? DcomMultiref(Hs, new_vec) : 0.0);
    time_h2_diag += elapsed(timer);
    ++h2_diag_calls;

    if (arnoldi_monitor_ortho && ((j + 1) % 5 == 0))
    {
      double gmax = gram_offdiag_max();
      max_ortho_seen = std::max(max_ortho_seen, gmax);
      std::cout << "arnoldi: Gram offdiag max @ step " << j + 1
                << " = " << gmax << std::endl;
      if (gmax > arnoldi_ortho_fail)
      {
        std::cout << "arnoldi: full Gram orthogonality lost — stopping." << std::endl;
        stop_reason = "ortho_lost";
      }
    }

    if (j + 1 >= min_iter)
    {
      const int hall_dim = j + 1;
      arma::vec eigval; arma::mat eigvec;
      timer = Clock::now();
      eigensolve_sub(hall_dim, eigval, eigvec);
      time_eigensolve += elapsed(timer);

      pick_physical_roots(eigval, eigvec, j, beta, e, vs, residuals);
      const int nwant = (int)e.n_elem;
      double max_resid = 0.0;
      for (int k = 0; k < nwant; ++k)
        max_resid = std::max(max_resid, residuals(k));

      // Soft-root diagnostics (full spectrum before selection)
      if ((j + 1) % 5 == 0 && arnoldi_prefer_positive)
      {
        int nsoft = 0;
        for (int k = 0; k < (int)eigval.n_elem; ++k)
          if (eigval(k) < arnoldi_soft_floor)
            ++nsoft;
        if (nsoft > 0)
          std::cout << "arnoldi @ step " << j + 1 << ": discarded " << nsoft
                    << " soft root(s) < " << arnoldi_soft_floor
                    << " MeV (lowest soft=" << eigval(0) << ")" << std::endl;
      }

      energy_history.push_back(e);
      if ((int)energy_history.size() > stable_window)
        energy_history.erase(energy_history.begin());

      bool energy_stable = false;
      double delta_window = 0.0;
      if ((int)energy_history.size() >= stable_window)
      {
        const arma::vec &e_old = energy_history.front();
        const arma::vec &e_new = energy_history.back();
        if (e_old.n_elem == e_new.n_elem)
        {
          delta_window = arma::max(arma::abs(e_new - e_old));
          energy_stable = (delta_window < arnoldi_energy_tol);
        }
      }
      const bool resid_ok = (max_resid < arnoldi_resid_tol);

      if ((j + 1) % 5 == 0 || (energy_stable && resid_ok))
      {
        std::cout << "arnoldi @ step " << j + 1 << ":";
        for (int k = 0; k < nwant; ++k)
          std::cout << " E" << k << "=" << e(k)
                    << " (r=" << residuals(k) << ")";
        std::cout << "  dE_win=" << delta_window
                  << "  max_ortho=" << max_ortho_seen << std::endl;
        if (arnoldi_check_expectation)
          report_expectation_check(e, vs, hall_dim, nwant,
                                   "arnoldi expectation check:");
        if (j + 1 == 100)
        {
          write_ritz_vectors(vs, hall_dim, nwant,
                             "output/arnoldi_step100_state");
          std::cout << "wrote Ritz wavefunctions at step 100 to output/arnoldi_step100_state*.op"
                    << std::endl;
        }
      }

      if (stop_reason == "ortho_lost")
        break;

      if (energy_stable && resid_ok)
      {
        std::cout << "Arnoldi converged at step " << j + 1
                  << "  (dE_window=" << delta_window
                  << " < " << arnoldi_energy_tol
                  << ", max_resid=" << max_resid
                  << " < " << arnoldi_resid_tol << ")" << std::endl;
        converged_flag = true;
        stop_reason = "converged";
        break;
      }
    }
    else if (stop_reason == "ortho_lost")
    {
      break;
    }
  }

  {
    int nb_cur = (int)h1v_cache.size();
    if (nb_cur == 0)
      nb_cur = std::min((int)lanczos_vector.size(), 1);
    arma::vec eigval_f; arma::mat eigvec_f;
    timer = Clock::now();
    eigensolve_sub(nb_cur, eigval_f, eigvec_f);
    time_eigensolve += elapsed(timer);
    pick_physical_roots(eigval_f, eigvec_f, -1, 0.0, e, vs, residuals);
  }

  std::cout << "Arnoldi finished with " << j_final + 1 << " steps"
            << "  reason=" << stop_reason
            << "  max_ortho=" << max_ortho_seen << std::endl;
  if (arnoldi_check_expectation)
    report_expectation_check(e, vs, (int)vs.n_rows,
                             std::min(state_want, (int)e.n_elem),
                             "final arnoldi expectation check:");
  for (int k = 0; k < (int)e.n_elem; ++k)
  {
    std::cout << "E(" << k << ") = " << e(k);
    if (k < (int)residuals.n_elem)
      std::cout << "  resid~" << residuals(k);
    std::cout << std::endl;
  }

  timer = Clock::now();
  std::vector<Operator> ritz_vecs;
  int nb_ritz = (int)h1v_cache.size();
  if (nb_ritz == 0)
    nb_ritz = std::min((int)lanczos_vector.size(), 1);
  nb_ritz = std::min(nb_ritz, (int)lanczos_vector.size());
  for (int k = 0; k < (int)e.n_elem; ++k)
  {
    Operator vec = lanczos_vector[0] * 0.0;
    for (int m = 0; m < std::min(nb_ritz, (int)vs.n_rows); ++m)
      vec = vec + (double)vs(m, k) * lanczos_vector[m];
    ritz_vecs.push_back(vec);
  }
  time_ritz += elapsed(timer);

  // H3 first-order PT on physical (H2) Ritz vectors — not mixed into Hall.
  if (!arnoldi_use_h3)
  {
    std::cout << "H3 first-order PT on physical Ritz (Hall was H2-only):" << std::endl;
    for (int k = 0; k < (int)ritz_vecs.size(); ++k)
    {
      double nrm = ComputeNorm(ritz_vecs[k], ritz_vecs[k]);
      double h3 = DcomMultiref(Hs, ritz_vecs[k]);
      double h3n = (std::abs(nrm) < 1e-30) ? 0.0 : h3 / nrm;
      std::cout << "  state " << k
                << ": E_H2=" << e(k)
                << "  H3/N=" << h3n
                << "  E_H2+H3=" << (e(k) + h3n) << std::endl;
    }
  }

  if (arnoldi_print_timing)
  {
    double total_time = elapsed(total_start);
    std::cout << "\nArnoldi timing summary (seconds):" << std::endl;
    std::cout << "  total                 " << total_time << std::endl;
    std::cout << "  HtcMultiref actions    " << time_h1_action << std::endl;
    std::cout << "  hall H1 inner products " << time_hall_h1 << std::endl;
    std::cout << "  hall H2 Dcom cross     " << time_hall_h2
              << "  calls=" << h2_cross_calls << std::endl;
    std::cout << "  hall norm shifts       " << time_hall_norm << std::endl;
    std::cout << "  Gram-Schmidt/project   " << time_gs << std::endl;
    std::cout << "  H2 Dcom diagonals      " << time_h2_diag
              << "  calls=" << h2_diag_calls << std::endl;
    std::cout << "  Dcom internal total    "
          << (dcom_time_ladder + dcom_time_223_231 + dcom_time_223_232
            + dcom_time_223_132 + dcom_time_overlap)
          << "  calls=" << dcom_profile_calls << std::endl;
    std::cout << "    ladder_multiref      " << dcom_time_ladder << std::endl;
    std::cout << "    comm223_231          " << dcom_time_223_231 << std::endl;
    std::cout << "    comm223_232          " << dcom_time_223_232 << std::endl;
    std::cout << "    comm223_132          " << dcom_time_223_132 << std::endl;
    std::cout << "    overlap              " << dcom_time_overlap << std::endl;
    std::cout << "  eigensolves            " << time_eigensolve << std::endl;
    std::cout << "  Ritz vector build      " << time_ritz << std::endl;
  }

  ArnoldiResult result;
  result.energies = e;
  result.eigvecs  = vs;
  result.ritz     = ritz_vecs;
  int hall_dim = (int)h1v_cache.size();
  if (hall_dim == 0)
    hall_dim = std::min((int)lanczos_vector.size(), 1);
  result.hall      = hall.submat(0, 0, hall_dim - 1, hall_dim - 1);
  result.residuals = residuals;
  result.max_ortho = max_ortho_seen;
  result.steps = j_final + 1;
  result.converged = converged_flag;
  result.stop_reason = stop_reason;
  return result;
}

// ============================================================
//  ArnoldiSolveH2
//  Arnoldi on the *full* decoupled H (no H_PP / H_QQ split):
//
//  Physical kernel is the double commutator quadratic form
//      Q(χ) = <χ|[H,χ]_2> + <χ|[H,χ]_3>  =  Htc + Dcom.
//  Only the 2b piece yields a new χ′ in the operator space, so:
//      Krylov matvec = HtcMultiref(Hs_full, ·)     // generate subspace
//      Hall_ab = 1/2 (Htc_ab + Htc_ba)              // hermitized H2
//              + Dcom polarization                  // H3 (default on)
//  Rayleigh–Ritz on Hall then matches ExpectationValue = Q/N.
//  No EraseValence, no E_val N shift.  Eigenvalues are Delta E.
// ============================================================
EOM::ArnoldiResult
EOM::ArnoldiSolveH2(Operator &vi, int max_iter, int state_want)
{
  using Clock = std::chrono::steady_clock;
  auto elapsed = [](Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
  };
  auto total_start = Clock::now();
  double time_h1_action = 0.0;
  double time_hall_h1 = 0.0;
  double time_hall_norm = 0.0;
  double time_gs = 0.0;
  double time_eigensolve = 0.0;
  double time_ritz = 0.0;

  // Working Hamiltonian: full decoupled H (vv kept). Do not EraseValence.
  if (!have_Hs_full)
  {
    force_decouple(Hs);
    Hs_full = Hs;
    have_Hs_full = true;
  }
  Operator &Hwork = Hs_full;

  const double bj_tol   = 1e-4;
  const double null_tol = 1e-4;
  const int    min_iter = state_want + 1;
  const int    stable_window = std::max(1, arnoldi_stable_window);

  std::vector<Operator> lanczos_vector;
  std::vector<Operator> h1v_cache; // Htc(v_i), reused for <v_j|Htc(v_i)>
  std::vector<double>   h3_diag;   // Dcom(Hwork, v_i) when arnoldi_use_h3
  arma::mat hall(max_iter, max_iter, arma::fill::zeros);
  double time_hall_h3 = 0.0;
  double time_h3_diag = 0.0;
  int h3_cross_calls = 0;
  int h3_diag_calls = 0;

  double nn = ComputeNorm(vi, vi);
  Operator v0 = vi / std::sqrt(nn);
  lanczos_vector.push_back(v0);
  const double cn0 = ComputeNorm(v0, v0);

  {
    auto timer = Clock::now();
    h3_diag.push_back(arnoldi_use_h3 ? DcomMultiref(Hwork, v0) : 0.0);
    time_h3_diag += elapsed(timer);
    ++h3_diag_calls;
  }

  arma::vec e(state_want, arma::fill::zeros);
  arma::mat vs(1, 1, arma::fill::zeros);
  arma::vec residuals(state_want, arma::fill::zeros);
  std::vector<arma::vec> energy_history;
  double max_ortho_seen = 0.0;
  bool   converged_flag = false;
  std::string stop_reason = "max_iter";
  int       j_final      = 0;

  auto eigensolve_sub = [&hall](int dim, arma::vec &ev, arma::mat &evec)
  {
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    // enforce Hermitian numerically before eig_sym
    sub = arma::symmatu(0.5 * (sub + sub.t()));
    arma::eig_sym(ev, evec, sub);
  };

  auto report_expectation_check =
      [this, &hall, &lanczos_vector, &Hwork](const arma::vec &ev, const arma::mat &evec,
                              int dim, int nshow, const std::string &tag)
  {
    std::cout << tag << std::endl;
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    sub = arma::symmatu(0.5 * (sub + sub.t()));
    for (int k = 0; k < nshow; ++k)
    {
      Operator ritz = lanczos_vector[0] * 0.0;
      for (int m = 0; m < dim; ++m)
        ritz = ritz + evec(m, k) * lanczos_vector[m];

      arma::vec coeff = evec.col(k).head(dim);
      double hall_expect = arma::as_scalar(coeff.t() * sub * coeff);
      double nrm = ComputeNorm(ritz, ritz);
      Operator h1ritz = HtcMultiref(Hwork, ritz);
      double h1exp = ComputeNorm(ritz, h1ritz);
      double h3exp = arnoldi_use_h3 ? DcomMultiref(Hwork, ritz) : 0.0;
      double direct_num = h1exp + h3exp;
      double direct_quotient = (std::abs(nrm) < 1e-30) ? 0.0 : direct_num / nrm;
      std::cout << "  state " << k
                << ": diag=" << ev(k)
                << "  hall=<c|H|c>=" << hall_expect
                << "  norm=" << nrm
                << "  H2=" << h1exp
                << "  H3=" << h3exp
                << "  direct_E=" << direct_quotient
                << "  delta_num=" << (direct_num - hall_expect)
                << "  delta_E=" << (direct_quotient - ev(k)) << std::endl;
    }
  };

  auto gram_offdiag_max = [this, &lanczos_vector]() -> double
  {
    double mx = 0.0;
    const int n = (int)lanczos_vector.size();
    for (int i = 0; i < n; ++i)
      for (int jj = i + 1; jj < n; ++jj)
        mx = std::max(mx, std::abs(ComputeNorm(lanczos_vector[i], lanczos_vector[jj])));
    return mx;
  };

  std::cout << "ArnoldiSolveH2: Krylov=Htc(2b), Hall=(Htc_ab+Htc_ba)/2"
            << (arnoldi_use_h3 ? " + Dcom(H3) [= ExpectationValue kernel]"
                               : " (H3 OFF — debug; ExpectationValue will differ)")
            << "  (no E_val shift)" << std::endl;

  for (int j = 0; j < max_iter - 1; ++j)
  {
    j_final = j;
    auto timer = Clock::now();

    Operator h1v_j = HtcMultiref(Hwork, lanczos_vector[j]);
    time_h1_action += elapsed(timer);
    h1v_cache.push_back(h1v_j);

    // Fill Hall column/row j: hermitized Htc + optional Dcom polarization.
    for (int i = 0; i <= j; ++i)
    {
      timer = Clock::now();
      double h2_sym;
      double h3_sym;
      if (i == j)
      {
        h2_sym = ComputeNorm(lanczos_vector[j], h1v_j);
        h3_sym = h3_diag[j];
      }
      else
      {
        double hij = ComputeNorm(lanczos_vector[i], h1v_j);
        double hji = ComputeNorm(lanczos_vector[j], h1v_cache[i]);
        h2_sym = 0.5 * (hij + hji);
        if (arnoldi_use_h3)
        {
          Operator v_sum = lanczos_vector[i] + lanczos_vector[j];
          auto t3 = Clock::now();
          double h3_cross = DcomMultiref(Hwork, v_sum);
          time_hall_h3 += elapsed(t3);
          ++h3_cross_calls;
          h3_sym = 0.5 * (h3_cross - h3_diag[i] - h3_diag[j]);
        }
        else
        {
          h3_sym = 0.0;
        }
      }
      time_hall_h1 += elapsed(timer);

      // No E_val / E_0 shift — full H already yields Delta E.
      hall(i, j) = hall(j, i) = h2_sym + h3_sym;
    }

    timer = Clock::now();
    Operator w = h1v_j * 1.0;
    if (arnoldi_use_projection)
      ProjectOprator(w);
    auto gs_pass = [&]() {
      for (int i = 0; i <= j; ++i)
      {
        double cij = ComputeNorm(lanczos_vector[i], w);
        w = w - cij * lanczos_vector[i];
      }
      if (arnoldi_use_projection)
        ProjectOprator(w);
    };
    gs_pass();
    gs_pass();

    double bj = ComputeNorm(w, w);
    time_gs += elapsed(timer);
    double bj_kernel = bj;

    if (std::abs(bj_kernel) < null_tol * cn0)
    {
      std::cout << "arnoldiH2: null vector (breakdown) at step " << j + 1
                << " (ComputeNorm=" << bj_kernel << "), stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      stop_reason = "null_breakdown";
      break;
    }
    if (std::abs(bj) < bj_tol)
    {
      std::cout << "arnoldiH2: exact breakdown at step " << j + 1
                << ", stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      stop_reason = "exact_breakdown";
      break;
    }
    if (bj < 0.0)
    {
      std::cout << "arnoldiH2: bj=" << bj << " < 0 at step " << j + 1
                << " (indefinite metric), stopping." << std::endl;
      arma::vec ev; arma::mat evec;
      timer = Clock::now();
      eigensolve_sub(j + 1, ev, evec);
      time_eigensolve += elapsed(timer);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      stop_reason = "indefinite_metric";
      break;
    }

    const double beta = std::sqrt(bj);
    Operator new_vec = w / beta;

    if (arnoldi_monitor_ortho)
    {
      double max_ov = 0.0;
      for (int i = 0; i <= j; ++i)
        max_ov = std::max(max_ov, std::abs(ComputeNorm(lanczos_vector[i], new_vec)));
      if (max_ov > arnoldi_ortho_warn)
      {
        for (int pass = 0; pass < 2; ++pass)
        {
          for (int i = 0; i <= j; ++i)
          {
            double cij = ComputeNorm(lanczos_vector[i], new_vec);
            new_vec = new_vec - cij * lanczos_vector[i];
          }
          if (arnoldi_use_projection)
            ProjectOprator(new_vec);
        }
        double nn_new = ComputeNorm(new_vec, new_vec);
        if (nn_new <= 0.0)
        {
          std::cout << "arnoldiH2: reortho produced non-positive norm at step "
                    << j + 1 << ", stopping." << std::endl;
          stop_reason = "reortho_breakdown";
          arma::vec ev; arma::mat evec;
          eigensolve_sub(j + 1, ev, evec);
          e  = ev.head(std::min(state_want, (int)ev.n_elem));
          vs = evec;
          break;
        }
        new_vec = new_vec / std::sqrt(nn_new);
        max_ov = 0.0;
        for (int i = 0; i <= j; ++i)
          max_ov = std::max(max_ov, std::abs(ComputeNorm(lanczos_vector[i], new_vec)));
      }
      max_ortho_seen = std::max(max_ortho_seen, max_ov);
      if (max_ov > arnoldi_ortho_warn)
        std::cout << "arnoldiH2: ortho warn at step " << j + 1
                  << "  max|<v_i|v_new>|=" << max_ov << std::endl;
      if (max_ov > arnoldi_ortho_fail)
      {
        std::cout << "arnoldiH2: orthogonality lost at step " << j + 1
                  << "  max|<v_i|v_new>|=" << max_ov
                  << " > fail tol " << arnoldi_ortho_fail
                  << " — stopping (Ritz unreliable)." << std::endl;
        stop_reason = "ortho_lost";
      }
    }

    lanczos_vector.push_back(new_vec);
    {
      auto t3 = Clock::now();
      h3_diag.push_back(arnoldi_use_h3 ? DcomMultiref(Hwork, new_vec) : 0.0);
      time_h3_diag += elapsed(t3);
      ++h3_diag_calls;
    }

    if (arnoldi_monitor_ortho && ((j + 1) % 5 == 0))
    {
      double gmax = gram_offdiag_max();
      max_ortho_seen = std::max(max_ortho_seen, gmax);
      std::cout << "arnoldiH2: Gram offdiag max @ step " << j + 1
                << " = " << gmax << std::endl;
      if (gmax > arnoldi_ortho_fail)
      {
        std::cout << "arnoldiH2: full Gram orthogonality lost — stopping."
                  << std::endl;
        stop_reason = "ortho_lost";
      }
    }

    if (j + 1 >= min_iter)
    {
      const int hall_dim = j + 1;
      arma::vec eigval; arma::mat eigvec;
      timer = Clock::now();
      eigensolve_sub(hall_dim, eigval, eigvec);
      time_eigensolve += elapsed(timer);

      const int nwant = std::min(state_want, (int)eigval.n_elem);
      e  = eigval.head(nwant);
      vs = eigvec;

      residuals.set_size(nwant);
      double max_resid = 0.0;
      for (int k = 0; k < nwant; ++k)
      {
        residuals(k) = beta * std::abs(eigvec(j, k));
        max_resid = std::max(max_resid, residuals(k));
      }

      energy_history.push_back(e);
      if ((int)energy_history.size() > stable_window)
        energy_history.erase(energy_history.begin());

      bool energy_stable = false;
      double delta_window = 0.0;
      if ((int)energy_history.size() >= stable_window)
      {
        const arma::vec &e_old = energy_history.front();
        const arma::vec &e_new = energy_history.back();
        if (e_old.n_elem == e_new.n_elem)
        {
          delta_window = arma::max(arma::abs(e_new - e_old));
          energy_stable = (delta_window < arnoldi_energy_tol);
        }
      }
      const bool resid_ok = (max_resid < arnoldi_resid_tol);

      if ((j + 1) % 5 == 0 || (energy_stable && resid_ok))
      {
        std::cout << "arnoldiH2 @ step " << j + 1 << ":";
        for (int k = 0; k < nwant; ++k)
          std::cout << " E" << k << "=" << e(k)
                    << " (r=" << residuals(k) << ")";
        std::cout << "  dE_win=" << delta_window
                  << "  max_ortho=" << max_ortho_seen << std::endl;
        if (arnoldi_check_expectation)
          report_expectation_check(eigval, eigvec, hall_dim, nwant,
                                   "arnoldiH2 expectation check:");
      }

      if (stop_reason == "ortho_lost")
        break;

      if (energy_stable && resid_ok)
      {
        std::cout << "ArnoldiSolveH2 converged at step " << j + 1
                  << "  (dE_window=" << delta_window
                  << " < " << arnoldi_energy_tol
                  << ", max_resid=" << max_resid
                  << " < " << arnoldi_resid_tol << ")" << std::endl;
        converged_flag = true;
        stop_reason = "converged";
        break;
      }
    }
    else if (stop_reason == "ortho_lost")
    {
      break;
    }
  }

  {
    int nb_cur = (int)h1v_cache.size();
    if (nb_cur == 0)
      nb_cur = std::min((int)lanczos_vector.size(), 1);
    arma::vec eigval_f; arma::mat eigvec_f;
    auto timer = Clock::now();
    eigensolve_sub(nb_cur, eigval_f, eigvec_f);
    time_eigensolve += elapsed(timer);
    const int nwant = std::min(state_want, (int)eigval_f.n_elem);
    e  = eigval_f.head(nwant);
    vs = eigvec_f;
    if ((int)residuals.n_elem != nwant)
      residuals = arma::vec(nwant, arma::fill::zeros);
  }

  std::cout << "ArnoldiSolveH2 finished with " << j_final + 1 << " steps"
            << "  reason=" << stop_reason
            << "  max_ortho=" << max_ortho_seen << std::endl;
  if (arnoldi_check_expectation)
    report_expectation_check(e, vs, (int)vs.n_rows,
                             std::min(state_want, (int)e.n_elem),
                             "final arnoldiH2 expectation check:");
  for (int k = 0; k < (int)e.n_elem; ++k)
  {
    std::cout << "E(" << k << ") = " << e(k);
    if (k < (int)residuals.n_elem)
      std::cout << "  resid~" << residuals(k);
    std::cout << std::endl;
  }

  auto timer = Clock::now();
  std::vector<Operator> ritz_vecs;
  int nb_ritz = (int)h1v_cache.size();
  if (nb_ritz == 0)
    nb_ritz = std::min((int)lanczos_vector.size(), 1);
  nb_ritz = std::min(nb_ritz, (int)lanczos_vector.size());
  for (int k = 0; k < (int)e.n_elem; ++k)
  {
    Operator vec = lanczos_vector[0] * 0.0;
    for (int m = 0; m < std::min(nb_ritz, (int)vs.n_rows); ++m)
      vec = vec + (double)vs(m, k) * lanczos_vector[m];
    ritz_vecs.push_back(vec);
  }
  time_ritz += elapsed(timer);

  if (arnoldi_print_timing)
  {
    double total_time = elapsed(total_start);
    std::cout << "\nArnoldiSolveH2 timing summary (seconds):" << std::endl;
    std::cout << "  total                 " << total_time << std::endl;
    std::cout << "  HtcMultiref actions    " << time_h1_action << std::endl;
    std::cout << "  hall Htc products      " << time_hall_h1 << std::endl;
    std::cout << "  hall H3 Dcom cross     " << time_hall_h3
              << "  calls=" << h3_cross_calls << std::endl;
    std::cout << "  H3 Dcom diagonals      " << time_h3_diag
              << "  calls=" << h3_diag_calls << std::endl;
    std::cout << "  Gram-Schmidt/project   " << time_gs << std::endl;
    std::cout << "  eigensolves            " << time_eigensolve << std::endl;
    std::cout << "  Ritz vector build      " << time_ritz << std::endl;
  }

  ArnoldiResult result;
  result.energies = e;
  result.eigvecs  = vs;
  result.ritz     = ritz_vecs;
  int hall_dim = (int)h1v_cache.size();
  if (hall_dim == 0)
    hall_dim = std::min((int)lanczos_vector.size(), 1);
  result.hall      = hall.submat(0, 0, hall_dim - 1, hall_dim - 1);
  result.residuals = residuals;
  result.max_ortho = max_ortho_seen;
  result.steps = j_final + 1;
  result.converged = converged_flag;
  result.stop_reason = stop_reason;
  return result;
}

// ============================================================
//  ArnoldiSolve_old
//  Original implementation kept verbatim for benchmarking.
//  Differences vs ArnoldiSolve:
//    * No i==j short-circuit; the polarization formula is applied on
//      the diagonal too (1/2 (<2v|H2|2v> - 2<v|H2|v>) = <v|H2|v>).
//    * Plain arma::eig_sym (no explicit symmatu symmetrization).
//    * No asymmetry diagnostic, no in-loop variational consistency print.
//  Inner products use ComputeNorm (matching current ArnoldiSolve).
//  In exact arithmetic ArnoldiSolve and ArnoldiSolve_old must produce the
//  same hall matrix; compare result.hall to verify equivalence.
// ============================================================
EOM::ArnoldiResult
EOM::ArnoldiSolve_old(Operator &vi, int max_iter, int state_want)
{
  const double tol      = 1e-8;
  const double bj_tol   = 1e-4;
  const double null_tol = 1e-4;
  const int    min_iter = state_want + 1;

  std::vector<Operator> lanczos_vector;
  std::vector<Operator> h1v_cache;
  std::vector<double>   h2_diag;
  arma::mat hall(max_iter, max_iter, arma::fill::zeros);

  // normalize initial vector
  double nn = ComputeNorm(vi, vi);
  Operator v0 = vi / std::sqrt(nn);
  lanczos_vector.push_back(v0);
  h2_diag.push_back(DcomMultiref(Hs, v0));
  const double cn0 = ComputeNorm(v0, v0);

  arma::vec e(state_want, arma::fill::zeros);
  arma::mat vs(1, 1, arma::fill::zeros);
  arma::vec prev_e;
  bool      prev_e_valid = false;
  int       j_final      = 0;

  for (int j = 0; j < max_iter - 1; ++j)
  {
    j_final = j;

    Operator h1v_j = HtcMultiref(Hs, lanczos_vector[j]);
    h1v_cache.push_back(h1v_j);

    // fill row / column j using full polarization identity (no diagonal shortcut)
    for (int i = 0; i <= j; ++i)
    {
      double h1ij   = ComputeNorm(lanczos_vector[i], h1v_j);
      double h1ji   = ComputeNorm(lanczos_vector[j], h1v_cache[i]);
      double h1_sym = 0.5 * (h1ij + h1ji);

      Operator v_sum  = lanczos_vector[i] + lanczos_vector[j];
      double h2_cross = DcomMultiref(Hs, v_sum);
      double h2_sym   = 0.5 * (h2_cross - h2_diag[i] - h2_diag[j]);

      hall(i, j) = hall(j, i) = h1_sym + h2_sym;
    }

    // double-pass classical Gram-Schmidt orthogonalization
    Operator w = h1v_j * 1.0;
    ProjectOprator(w);
    for (int pass = 0; pass < 2; ++pass)
    {
      for (int i = 0; i <= j; ++i)
      {
        double cij = ComputeNorm(lanczos_vector[i], w);
        w = w - cij * lanczos_vector[i];
      }
      ProjectOprator(w);
    }

    double bj        = ComputeNorm(w, w);
    double bj_kernel = bj;
    //std::cout << bj << " " << bj_kernel << " kernels " << std::endl;

    // null-space breakdown
    if (std::abs(bj_kernel) < null_tol * cn0)
    {
      std::cout << "arnoldi_old: null vector (breakdown) at step " << j + 1
                << " (ComputeNorm=" << bj_kernel << "), stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // exact breakdown
    if (std::abs(bj) < bj_tol)
    {
      std::cout << "arnoldi_old: exact breakdown at step " << j + 1 << ", stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // negative norm
    if (bj < 0.0)
    {
      std::cout << "arnoldi_old: bj=" << bj << " < 0 at step " << j + 1
                << " (indefinite metric), stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // normal step
    Operator new_vec = w / std::sqrt(bj);
    lanczos_vector.push_back(new_vec);
    h2_diag.push_back(DcomMultiref(Hs, new_vec));

    // eigenvalue check / convergence
    if (j + 1 >= min_iter)
    {
      arma::mat sub = hall.submat(0, 0, j, j);
      arma::vec eigval; arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, sub);

      e  = eigval.head(std::min(state_want, (int)eigval.n_elem));
      vs = eigvec;

      if ((j + 1) % 5 == 0)
      {
        std::cout << "arnoldi_old eigenvalues @ step " << j + 1 << ": ";
        for (int k = 0; k < (int)e.n_elem; ++k)
          std::cout << e(k) << " ";
        std::cout << std::endl;
      }

      if (prev_e_valid && prev_e.n_elem == e.n_elem)
      {
        double delta = arma::max(arma::abs(e - prev_e));
        double scale = std::max(1.0, arma::max(arma::abs(e)));
        if (delta < tol || delta / scale < tol)
        {
          std::cout << "Arnoldi_old converged at step " << j + 1 << std::endl;
          break;
        }
      }
      prev_e       = e;
      prev_e_valid = true;
    }
  }

  // Final eigensolution on the filled projected subspace only.
  {
    int nb_cur = (int)h1v_cache.size();
    if (nb_cur == 0)
      nb_cur = std::min((int)lanczos_vector.size(), 1);
    arma::mat sub = hall.submat(0, 0, nb_cur - 1, nb_cur - 1);
    arma::vec eigval_f; arma::mat eigvec_f;
    arma::eig_sym(eigval_f, eigvec_f, sub);
    e  = eigval_f.head(std::min(state_want, (int)eigval_f.n_elem));
    vs = eigvec_f;
  }

  std::cout << "Arnoldi_old finished with " << j_final + 1 << " steps" << std::endl;
  for (int k = 0; k < (int)e.n_elem; ++k)
    std::cout << "E(" << k << ") = " << e(k) << std::endl;

  // Build Ritz vectors.
  std::vector<Operator> ritz_vecs;
  int nb = (int)lanczos_vector.size();
  for (int k = 0; k < (int)e.n_elem; ++k)
  {
    Operator vec = lanczos_vector[0] * 0.0;
    for (int m = 0; m < std::min(nb, (int)vs.n_rows); ++m)
      vec = vec + (double)vs(m, k) * lanczos_vector[m];
    ritz_vecs.push_back(vec);
  }

  ArnoldiResult result;
  result.energies = e;
  result.eigvecs  = vs;
  result.ritz     = ritz_vecs;
  int hall_dim = (int)h1v_cache.size();
  if (hall_dim == 0)
    hall_dim = std::min((int)lanczos_vector.size(), 1);
  result.hall     = hall.submat(0, 0, hall_dim - 1, hall_dim - 1);
  return result;
}

EOM::ArnoldiTraceDiffResult
EOM::CompareArnoldiHallBuild(Operator &vi, int max_iter, double tol)
{
  ArnoldiTraceDiffResult result;
  const double bj_tol   = 1e-4;
  const double null_tol = 1e-4;

  std::vector<Operator> lanczos_new;
  std::vector<Operator> lanczos_old;
  std::vector<Operator> h1v_cache_new;
  std::vector<Operator> h1v_cache_old;
  std::vector<double>   h2_diag_new;
  std::vector<double>   h2_diag_old;

  double nn = ComputeNorm(vi, vi);
  Operator v0_new = vi / std::sqrt(nn);
  Operator v0_old = vi / std::sqrt(nn);
  lanczos_new.push_back(v0_new);
  lanczos_old.push_back(v0_old);
  h2_diag_new.push_back(DcomMultiref(Hs, v0_new));
  h2_diag_old.push_back(DcomMultiref(Hs, v0_old));
  const double cn0 = ComputeNorm(v0_new, v0_new);

  for (int j = 0; j < max_iter - 1; ++j)
  {
    Operator h1v_new = HtcMultiref(Hs, lanczos_new[j]);
    Operator h1v_old = HtcMultiref(Hs, lanczos_old[j]);
    h1v_cache_new.push_back(h1v_new);
    h1v_cache_old.push_back(h1v_old);

    for (int i = 0; i <= j; ++i)
    {
      double h1ij_new   = ComputeNorm(lanczos_new[i], h1v_new);
      double h1ji_new   = ComputeNorm(lanczos_new[j], h1v_cache_new[i]);
      double h1_sym_new = 0.5 * (h1ij_new + h1ji_new);

      double h2_cross_new;
      double h2_sym_new;
      if (i == j)
      {
        Operator v_sum_new = lanczos_new[j] + lanczos_new[j];
        h2_cross_new = DcomMultiref(Hs, v_sum_new);
        h2_sym_new   = 0.5 * (h2_cross_new - h2_diag_new[j] - h2_diag_new[j]);
      }
      else
      {
        Operator v_sum_new = lanczos_new[i] + lanczos_new[j];
        h2_cross_new = DcomMultiref(Hs, v_sum_new);
        h2_sym_new   = 0.5 * (h2_cross_new - h2_diag_new[i] - h2_diag_new[j]);
      }

      double h1ij_old   = ComputeNorm(lanczos_old[i], h1v_old);
      double h1ji_old   = ComputeNorm(lanczos_old[j], h1v_cache_old[i]);
      double h1_sym_old = 0.5 * (h1ij_old + h1ji_old);
      Operator v_sum_old = lanczos_old[i] + lanczos_old[j];
      double h2_cross_old = DcomMultiref(Hs, v_sum_old);
      double h2_sym_old   = 0.5 * (h2_cross_old - h2_diag_old[i] - h2_diag_old[j]);

      double norm_new = ComputeNorm(lanczos_new[i], lanczos_new[j]);
      double norm_old = ComputeNorm(lanczos_old[i], lanczos_old[j]);
      (void)norm_new;
      (void)norm_old;
      double hall_new = h1_sym_new + h2_sym_new;
      double hall_old = h1_sym_old + h2_sym_old;
      double abs_diff = std::abs(hall_new - hall_old);
      double rel_diff = abs_diff / (std::abs(hall_old) + 1e-30);
      result.max_abs_diff = std::max(result.max_abs_diff, abs_diff);
      result.max_rel_diff = std::max(result.max_rel_diff, rel_diff);

      if (abs_diff > tol)
      {
        result.found = true;
        result.step = j;
        result.i = i;
        result.j = j;
        result.hall_new = hall_new;
        result.hall_old = hall_old;
        result.delta_hall = hall_new - hall_old;
        result.h1_sym_new = h1_sym_new;
        result.h1_sym_old = h1_sym_old;
        result.h2_cross_new = h2_cross_new;
        result.h2_cross_old = h2_cross_old;
        result.h2_sym_new = h2_sym_new;
        result.h2_sym_old = h2_sym_old;
        result.h2_diag_i_new = h2_diag_new[i];
        result.h2_diag_j_new = h2_diag_new[j];
        result.h2_diag_i_old = h2_diag_old[i];
        result.h2_diag_j_old = h2_diag_old[j];
        return result;
      }
    }

    Operator w_new = h1v_new * 1.0;
    Operator w_old = h1v_old * 1.0;
    ProjectOprator(w_new);
    ProjectOprator(w_old);
    for (int pass = 0; pass < 2; ++pass)
    {
      for (int i = 0; i <= j; ++i)
      {
        double cij_new = ComputeNorm(lanczos_new[i], w_new);
        double cij_old = ComputeNorm(lanczos_old[i], w_old);
        w_new = w_new - cij_new * lanczos_new[i];
        w_old = w_old - cij_old * lanczos_old[i];
      }
      ProjectOprator(w_new);
      ProjectOprator(w_old);
      
    }

    double bj_new = ComputeNorm(w_new, w_new);
    double bj_old = ComputeNorm(w_old, w_old);
    double bj_kernel_new = bj_new;
    double bj_kernel_old = bj_old;

    bool break_new = (std::abs(bj_kernel_new) < null_tol * cn0)
                  || (std::abs(bj_new) < bj_tol) || (bj_new < 0.0);
    bool break_old = (std::abs(bj_kernel_old) < null_tol * cn0)
                  || (std::abs(bj_old) < bj_tol) || (bj_old < 0.0);

    if (break_new || break_old)
      break;

    Operator new_vec_new = w_new / std::sqrt(bj_new);
    Operator new_vec_old = w_old / std::sqrt(bj_old);
    lanczos_new.push_back(new_vec_new);
    lanczos_old.push_back(new_vec_old);
    h2_diag_new.push_back(DcomMultiref(Hs, new_vec_new));
    h2_diag_old.push_back(DcomMultiref(Hs, new_vec_old));
  }

  return result;
}

// ============================================================
//  ReadTdm
//  Translated from read_tdm() in run/lanczos.py
// ============================================================

/// Read a transition density matrix file and populate a scalar 2-body Operator.
Operator EOM::ReadTdm(const std::string &tdm_file)
{
  std::ifstream infile(tdm_file);
  if (!infile)
    throw std::runtime_error("ReadTdm: cannot open file " + tdm_file);

  std::vector<std::string> lines;
  {
    std::string buf;
    while (std::getline(infile, buf))
      lines.push_back(buf);
  }
  infile.close();

  int lidx = 0;

  // --- line 0: J_total ---
  {
    std::istringstream ss(lines[lidx++]);
    double jtotal;
    ss >> jtotal;
    double factor = std::sqrt(2.0 * jtotal + 1.0);

    // --- line 1: number of single-particle orbits ---
    {
      std::istringstream ss2(lines[lidx++]);
      int norb;
      ss2 >> norb;

      // Pre-scan orbit lines to determine emax BEFORE AddOrbit,
      // because orbits_3body_space_ is populated during AddOrbit
      // only when 2n+l <= emax_3body_, which defaults to 0.
      struct OrbEntry { int idx, l, tz2; };
      struct RawOrb { int n, l, j2, tz2; };
      std::vector<RawOrb> raw_orbits(norb);
      int emax = 0;
      int lidx_save = lidx;
      for (int obs = 0; obs < norb; ++obs)
      {
        std::istringstream sl(lines[lidx++]);
        int dummy;
        sl >> dummy >> raw_orbits[obs].n >> raw_orbits[obs].l
                    >> raw_orbits[obs].j2 >> raw_orbits[obs].tz2;
        emax = std::max(emax, 2 * raw_orbits[obs].n + raw_orbits[obs].l);
      }
      (void)lidx_save;

      // Use default then set emax_3body_ via SetEmax3Body BEFORE AddOrbit,
      // because AddOrbit checks 2n+l <= emax_3body_ to populate orbits_3body_space_.
      rdm_modelspace = ModelSpace();
      rdm_modelspace.SetEmax(emax);          // sets SixJCache and Emax correctly
      rdm_modelspace.SetEmax3Body(emax);     // allow all orbits into orbits_3body_space_
      rdm_modelspace.SetE3max(3 * emax);   // allow all (p,q,r) ket triples in Setup3bKets

      std::vector<OrbEntry> ob_idx(norb);
      for (int obs = 0; obs < norb; ++obs)
      {
        rdm_modelspace.AddOrbit(raw_orbits[obs].n, raw_orbits[obs].l,
                                raw_orbits[obs].j2, raw_orbits[obs].tz2, 0.0, 1);
        ob_idx[obs] = {obs, raw_orbits[obs].l, raw_orbits[obs].tz2};
      }
      rdm_modelspace.FindEFermi();
      rdm_modelspace.SetupKets();
      rdm_modelspace.Setup3bKets();   // must be called separately after SetupKets()
      rdm_ms = &rdm_modelspace;

      Operator ops(*rdm_ms, 0, 0, 0, 3);

      ops.ThreeBody.SetMode("pn");
      
      ops *= 0.0;

      // --- one-body density matrix elements ---
      {
        std::istringstream sn(lines[lidx++]);
        int n_obtd;
        sn >> n_obtd;
        for (int obs = 0; obs < n_obtd; ++obs)
        {
          std::istringstream sl2(lines[lidx++]);
          std::vector<std::string> tok;
          std::string w;
          while (sl2 >> w) tok.push_back(w);

          int aa = ob_idx[std::stoi(tok[1]) - 1].idx;
          int bb = ob_idx[std::stoi(tok[2]) - 1].idx;
          double rd = std::stod(tok.back()) / factor;
          ops.SetOneBody(aa, bb, rd);
        }
      }

      // --- two-body density matrix elements ---
      // Format per line:
      //   label  a b c d  Jab Jcd  [ignored columns...]  value
      // Orbit indices are 1-based in file order. Any trailing metadata columns
      // between Jcd and the final value are ignored.
      {
        std::istringstream sn(lines[lidx++]);
        int n_tbtd;
        sn >> n_tbtd;
        for (int obs = 0; obs < n_tbtd; ++obs)
        {
          std::istringstream sl2(lines[lidx++]);
          std::vector<std::string> tok;
          std::string w;
          while (sl2 >> w) tok.push_back(w);

          if (tok.size() < 8)
            throw std::runtime_error("ReadTdm: malformed TBTD line in " + tdm_file);

          int ia = std::stoi(tok[1]) - 1;
          int ib = std::stoi(tok[2]) - 1;
          int ic = std::stoi(tok[3]) - 1;
          int id = std::stoi(tok[4]) - 1;

          int aa = ob_idx[ia].idx;
          int bb = ob_idx[ib].idx;
          int cc = ob_idx[ic].idx;
          int dd = ob_idx[id].idx;

          int jij = std::stoi(tok[5]);
          int pij = (ob_idx[ia].l + ob_idx[ib].l) % 2;
          int tij = (ob_idx[ia].tz2 + ob_idx[ib].tz2) / 2;

          int jkl = std::stoi(tok[6]);
          int pkl = (ob_idx[ic].l + ob_idx[id].l) % 2;
          int tkl = (ob_idx[ic].tz2 + ob_idx[id].tz2) / 2;

          double rd = std::stod(tok.back()) / factor;
          ops.SetTwoBody(jij, pij, tij, jkl, pkl, tkl, aa, bb, cc, dd, rd);
        }
      }

      // --- three-body density matrix elements ---
      // Format per line (1-based orbit indices, 2*Jab, 2*Jde, 2*Jtot):
      //   label  a b c  d e f  2Jab 2Jde 2Jtot  value
      // Fortran must write (a,b,c) with a<=b<=c and (d,e,f) with d<=e<=f
      // using 1-based indices that match the orbit order in the file header.
      // The C++ canonical order for a ket (p,q,r,Jpq) requires p<=q<=r.
      // Since orbits are added to rdm_modelspace in file order,
      // 1-based index k maps to rdm orbit index k-1, so a<=b<=c (1-based)
      // is equivalent to (a-1)<=(b-1)<=(c-1) which is p<=q<=r canonical.
      // Stored always; contractions gated by use_rdm3 (default false).
      {
        std::istringstream sn(lines[lidx++]);
        int n_3btd;
        sn >> n_3btd;
        for (int obs = 0; obs < n_3btd; ++obs)
        {
          std::istringstream sl2(lines[lidx++]);
          std::vector<std::string> tok;
          std::string w;
          while (sl2 >> w) tok.push_back(w);

          // orbit indices in rdm_modelspace (0-based)
          int aa = ob_idx[std::stoi(tok[1]) - 1].idx;  // bra p
          int bb = ob_idx[std::stoi(tok[2]) - 1].idx;  // bra q
          int cc = ob_idx[std::stoi(tok[3]) - 1].idx;  // bra r
          int dd = ob_idx[std::stoi(tok[4]) - 1].idx;  // ket p
          int ee = ob_idx[std::stoi(tok[5]) - 1].idx;  // ket q
          int ff = ob_idx[std::stoi(tok[6]) - 1].idx;  // ket r

          int two_jab = std::stoi(tok[7]);  // 2*Jab (must be even → integer Jab)
          int two_jde = std::stoi(tok[8]);  // 2*Jde
          int two_tot = std::stoi(tok[9]);  // 2*Jtot

          int jab = two_jab / 2;  // integer Jab (coupling of bra p,q pair)
          int jde = two_jde / 2;  // integer Jde (coupling of ket p,q pair)

          double rd = std::stod(tok.back()) / factor;
          // OSM TRBTD files already include an extra sqrt(2J+1) factor
          // in the stored 3-body value, unlike the C++ writer.
          rd /= std::sqrt(double(two_tot) + 1.0);

          // Identify channels by (2*Jtot, parity, twoTz)
          Orbit& oa  = rdm_ms->GetOrbit(aa);
          Orbit& ob_ = rdm_ms->GetOrbit(bb);
          Orbit& oc  = rdm_ms->GetOrbit(cc);
          Orbit& od  = rdm_ms->GetOrbit(dd);
          Orbit& oe  = rdm_ms->GetOrbit(ee);
          Orbit& of_ = rdm_ms->GetOrbit(ff);

          int par_bra  = (oa.l + ob_.l + oc.l) % 2;
          int twoTz_bra = oa.tz2 + ob_.tz2 + oc.tz2;
          int par_ket  = (od.l + oe.l + of_.l) % 2;
          int twoTz_ket = od.tz2 + oe.tz2 + of_.tz2;

          size_t ch_bra = rdm_ms->GetThreeBodyChannelIndex(two_tot, par_bra, twoTz_bra);
          size_t ch_ket = rdm_ms->GetThreeBodyChannelIndex(two_tot, par_ket, twoTz_ket);

          // Skip if channel not present in this model space
          if (ch_bra == (size_t)-1 || ch_ket == (size_t)-1) continue;

          ThreeBodyChannel& TBC_bra = rdm_ms->GetThreeBodyChannel(ch_bra);
          ThreeBodyChannel& TBC_ket = rdm_ms->GetThreeBodyChannel(ch_ket);

          // Direct lookup: (p,q,r,Jpq) must be in canonical order p<=q<=r
          size_t local_ibra = TBC_bra.GetLocalIndex(aa, bb, cc, jab);
          size_t local_iket = TBC_ket.GetLocalIndex(dd, ee, ff, jde);

          // SetME_pn_ch directly sets the storage element — no recoupling
          ops.ThreeBody.SetME_pn_ch(ch_bra, ch_ket, local_ibra, local_iket,
                                    (ThreeBME_type)rd);
        }
      }
      return ops;
    } // norb scope
  } // jtotal / factor scope
  return Operator(); // unreachable
}

// ============================================================
//  WriteTdm
//  Write op in the exact format ReadTdm reads:
//    line 0  : jtotal (always 0.0 for scalar density)
//    line 1  : norb
//    norb lines: idx n l j2 tz2
//    line    : n_obtd
//    n_obtd lines: OBTD a+1 b+1 value*sqrt(2J+1)
//    line    : n_tbtd
//    n_tbtd lines: TBTD a+1 b+1 c+1 d+1 Jab Jcd value*sqrt(2J+1)
//                  (ReadTdm also accepts extra ignored columns before value)
//    line    : n_3btd
//    n_3btd lines: TRBTD a+1 b+1 c+1 d+1 e+1 f+1 2Jab 2Jde 2Jtot value*sqrt(2J+1)
//  3-body elements are written in native channel/ket memory order (ch_bra<=ch_ket,
//  ibra<=iket when ch_bra==ch_ket).
// ============================================================
void EOM::WriteTdm(const Operator &op, const std::string &filename) const
{
  ModelSpace *ms = op.modelspace;
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("WriteTdm: cannot open " + filename);

  // jtotal = 0 for a scalar density operator
  const double jtotal = 0.0;
  const double factor = std::sqrt(2.0 * jtotal + 1.0); // = 1
  out << std::fixed << std::setprecision(1) << jtotal << "\n";

  // --- orbits ---
  int norb = ms->GetNumberOrbits();
  out << norb << "\n";
  for (int i = 0; i < norb; ++i)
  {
    Orbit &oi = ms->GetOrbit(i);
    out << (i + 1) << " " << oi.n << " " << oi.l << " " << oi.j2 << " " << oi.tz2 << "\n";
  }

  // --- 1-body ---
  // count non-zero elements first
  std::vector<std::tuple<int,int,double>> ob_lines;
  for (int a = 0; a < norb; ++a)
    for (int b = 0; b < norb; ++b)
    {
      double v = op.OneBody(a, b);
      if (v != 0.0) ob_lines.emplace_back(a, b, v);
    }
  out << ob_lines.size() << "\n";
  for (auto &t : ob_lines)
    out << "OBTD " << (std::get<0>(t) + 1) << " " << (std::get<1>(t) + 1) << " " << std::setprecision(10) << std::get<2>(t) * factor << "\n";

  // --- 2-body ---
  std::vector<std::tuple<int,int,int,int,int,double>> tb_lines;
  int nch2 = ms->GetNumberTwoBodyChannels();
  for (int ch = 0; ch < nch2; ++ch)
  {
    TwoBodyChannel &tbc = ms->GetTwoBodyChannel(ch);
    int J2  = tbc.J;
    int nk  = tbc.GetNumberKets();
    for (int ib = 0; ib < nk; ++ib)
    {
      Ket &kb = tbc.GetKet(ib);
      for (int ik = ib; ik < nk; ++ik)
      {
        Ket &kk = tbc.GetKet(ik);
        double v = op.TwoBody.GetTBME_norm(ch, ib, ik);
        if (v != 0.0)
          tb_lines.emplace_back(kb.p, kb.q, kk.p, kk.q, J2, v);
      }
    }
  }
  out << tb_lines.size() << "\n";
  for (auto &t : tb_lines)
    out << "TBTD " << (std::get<0>(t)+1) << " " << (std::get<1>(t)+1) << " " << (std::get<2>(t)+1) << " " << (std::get<3>(t)+1)
        << " " << (std::get<4>(t)*2) << " " << (std::get<4>(t)*2)
        << " " << std::setprecision(10) << std::get<5>(t) * factor << "\n";

  // --- 3-body: iterate in native memory order ---
  size_t nch3 = ms->GetNumberThreeBodyChannels();
  std::vector<std::tuple<int,int,int,int, int,int,int,int, int, double>> lines3b;

  for (size_t ch_bra = 0; ch_bra < nch3; ++ch_bra)
  {
    ThreeBodyChannel &TBC_bra = ms->GetThreeBodyChannel(ch_bra);
    size_t nk_bra = TBC_bra.GetNumber3bKets();
    for (size_t ch_ket = ch_bra; ch_ket < nch3; ++ch_ket)
    {
      ThreeBodyChannel &TBC_ket = ms->GetThreeBodyChannel(ch_ket);
      size_t nk_ket = TBC_ket.GetNumber3bKets();

      // Check this channel pair is actually stored
      auto &ch_start = op.ThreeBody.Get_ch_start();
      if (ch_start.find({ch_bra, ch_ket}) == ch_start.end()) continue;

      int twoJ_bra = TBC_bra.twoJ;
      int twoJ_ket = TBC_ket.twoJ;
      // For a scalar operator ch_bra == ch_ket always (same twoJ, par, Tz)

      for (size_t ibra = 0; ibra < nk_bra; ++ibra)
      {
        Ket3 &kb = TBC_bra.GetKet(ibra);
        size_t iket_start = (ch_bra == ch_ket) ? ibra : 0;
        for (size_t iket = iket_start; iket < nk_ket; ++iket)
        {
          Ket3 &kk = TBC_ket.GetKet(iket);
          double v = op.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket);
          if (v == 0.0) continue;
          lines3b.emplace_back(
            (int)kb.p, (int)kb.q, (int)kb.r, kb.Jpq * 2,
            (int)kk.p, (int)kk.q, (int)kk.r, kk.Jpq * 2,
            twoJ_bra,   // bra==ket channel so same twoJ
            v * factor);
        }
      }
    }
  }

  out << lines3b.size() << "\n";
  for (auto &t : lines3b)
    out << "TRBTD "
        << (std::get<0>(t)+1) << " " << (std::get<1>(t)+1) << " " << (std::get<2>(t)+1) << " "
        << (std::get<3>(t)) << " "
        << (std::get<4>(t)+1) << " " << (std::get<5>(t)+1) << " " << (std::get<6>(t)+1) << " "
        << (std::get<7>(t)) << " " << (std::get<8>(t)) << " "
        << std::setprecision(10) << std::get<9>(t) << "\n";
}

// ============================================================
//  SR dagger EOM (1-PA / 1-PR) — ladder on odd-leg operators
// ============================================================

index_t EOM::PickQOrbitParticle(int J2, int parity, int itz) const
{
  int pick = -1;
  for (auto i : modelspace->all_orbits)
  {
    Orbit &o = modelspace->GetOrbit(i);
    if (o.j2 != J2 or (o.l % 2) != parity or o.tz2 != itz)
      continue;
    if (modelspace->holes.count(i))
      continue;
    if (pick < 0)
      pick = (int)i;
  }
  if (pick < 0)
  {
    std::ostringstream oss;
    oss << "EOM::PickQOrbitParticle: no particle orbit for J2=" << J2
        << " parity=" << parity << " tz2=" << itz
        << " (all_orbits=" << modelspace->all_orbits.size()
        << " holes=" << modelspace->holes.size()
        << " particles=" << modelspace->particles.size() << ")";
    throw std::runtime_error(oss.str());
  }
  return (index_t)pick;
}

index_t EOM::PickQOrbitHole(int J2, int parity, int itz) const
{
  int pick = -1;
  int best_n = 9999;
  for (auto i : modelspace->all_orbits)
  {
    Orbit &o = modelspace->GetOrbit(i);
    if (o.j2 != J2 or (o.l % 2) != parity or o.tz2 != itz)
      continue;
    if (!modelspace->holes.count(i))
      continue;
    if ((int)o.n < best_n)
    {
      best_n = (int)o.n;
      pick = (int)i;
    }
  }
  if (pick < 0)
  {
    std::ostringstream oss;
    oss << "EOM::PickQOrbitHole: no hole orbit for J2=" << J2
        << " parity=" << parity << " tz2=" << itz
        << " (all_orbits=" << modelspace->all_orbits.size()
        << " holes=" << modelspace->holes.size() << ")";
    throw std::runtime_error(oss.str());
  }
  return (index_t)pick;
}

Operator EOM::MakeEmptyDagger(index_t Q) const
{
  Operator R(*modelspace);
  R.SetNumberLegs(3);
  R.SetQSpaceOrbit(Q);
  R.SetNonHermitian();
  R.EraseOneBody();
  R.EraseThreeLeg();
  return R;
}

Operator EOM::LadderPA(const Operator &R) const
{
  index_t Q = R.GetQSpaceOrbit();
  Orbit &oQ = modelspace->GetOrbit(Q);
  Operator out = MakeEmptyDagger(Q);

  for (auto p : R.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (modelspace->holes.count(p))
      continue;
    out.OneBody(p, 0) = R.OneBody(p, 0);
  }

  for (size_t ch = 0; ch < modelspace->GetNumberTwoBodyChannels(); ++ch)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
    auto &pp = tbc.GetKetIndex_pp();
    if (pp.empty())
      continue;
    auto &mat_in = R.ThreeLeg.GetMatrix(ch);
    auto &mat_out = out.ThreeLeg.GetMatrix(ch);
    for (auto ibra : pp)
    {
      for (auto h : modelspace->holes)
        mat_out(ibra, h) = mat_in(ibra, h);
    }
  }
  return out;
}

Operator EOM::LadderPR(const Operator &R) const
{
  index_t Q = R.GetQSpaceOrbit();
  Orbit &oQ = modelspace->GetOrbit(Q);
  Operator out = MakeEmptyDagger(Q);

  for (auto h : R.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (!modelspace->holes.count(h))
      continue;
    out.OneBody(h, 0) = R.OneBody(h, 0);
  }

  for (size_t ch = 0; ch < modelspace->GetNumberTwoBodyChannels(); ++ch)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
    auto &hh = tbc.GetKetIndex_hh();
    if (hh.empty())
      continue;
    auto &mat_in = R.ThreeLeg.GetMatrix(ch);
    auto &mat_out = out.ThreeLeg.GetMatrix(ch);
    for (auto ibra : hh)
    {
      for (auto p : modelspace->particles)
        mat_out(ibra, p) = mat_in(ibra, p);
    }
  }
  return out;
}

double EOM::OverlapPA(const Operator &R1, const Operator &R2) const
{
  index_t Q = R1.GetQSpaceOrbit();
  Orbit &oQ = modelspace->GetOrbit(Q);
  double ov = 0.0;

  for (auto p : R1.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (modelspace->holes.count(p))
      continue;
    ov += R1.OneBody(p, 0) * R2.OneBody(p, 0);
  }

  for (size_t ch = 0; ch < modelspace->GetNumberTwoBodyChannels(); ++ch)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
    auto &pp = tbc.GetKetIndex_pp();
    if (pp.empty())
      continue;
    auto &m1 = R1.ThreeLeg.GetMatrix(ch);
    auto &m2 = R2.ThreeLeg.GetMatrix(ch);
    for (auto ibra : pp)
      for (auto h : modelspace->holes)
        ov += m1(ibra, h) * m2(ibra, h);
  }
  return ov;
}

double EOM::OverlapPR(const Operator &R1, const Operator &R2) const
{
  index_t Q = R1.GetQSpaceOrbit();
  Orbit &oQ = modelspace->GetOrbit(Q);
  double ov = 0.0;

  for (auto h : R1.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (!modelspace->holes.count(h))
      continue;
    ov += R1.OneBody(h, 0) * R2.OneBody(h, 0);
  }

  for (size_t ch = 0; ch < modelspace->GetNumberTwoBodyChannels(); ++ch)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
    auto &hh = tbc.GetKetIndex_hh();
    if (hh.empty())
      continue;
    auto &m1 = R1.ThreeLeg.GetMatrix(ch);
    auto &m2 = R2.ThreeLeg.GetMatrix(ch);
    for (auto ibra : hh)
      for (auto p : modelspace->particles)
        ov += m1(ibra, p) * m2(ibra, p);
  }
  return ov;
}

Operator EOM::HtcPA(const Operator &R) const
{
  return LadderPA(Commutator::Commutator(Hs, R));
}

Operator EOM::HtcPR(const Operator &R) const
{
  Operator Rp = LadderPR(Commutator::Commutator(Hs, R));
  Rp *= -1.0;
  return Rp;
}

Operator EOM::PureAdagSeed(index_t Q) const
{
  Operator R = MakeEmptyDagger(Q);
  R.OneBody(Q, 0) = 1.0;
  return LadderPA(R);
}

std::pair<arma::vec, std::vector<Operator>>
EOM::LanczosDaggerSolve(Operator &vi, int max_iter, int state_want, bool particle_removed)
{
  std::vector<Operator> basis;
  arma::mat hall(max_iter, max_iter, arma::fill::zeros);

  double nn = particle_removed ? OverlapPR(vi, vi) : OverlapPA(vi, vi);
  if (nn <= 0.0)
    throw std::runtime_error("EOM::LanczosDaggerSolve: initial vector has non-positive norm");
  basis.push_back(vi / std::sqrt(nn));

  int j_final = 0;
  double bj = 0.0;

  for (int j = 0; j < max_iter; ++j)
  {
    j_final = j;
    Operator w = particle_removed ? HtcPR(basis[j]) : HtcPA(basis[j]);
    double ai = particle_removed ? OverlapPR(w, basis[j]) : OverlapPA(w, basis[j]);
    hall(j, j) = ai;

    w = w - ai * basis[j];
    if (j > 0)
      w = w - bj * basis[j - 1];

    for (int i = 0; i <= j; ++i)
    {
      double cij = particle_removed ? OverlapPR(w, basis[i]) : OverlapPA(w, basis[i]);
      w = w - cij * basis[i];
    }

    bj = particle_removed ? OverlapPR(w, w) : OverlapPA(w, w);
    if (bj < 0.0 && std::abs(bj) < 1e-14)
      bj = 0.0;
    if (bj < 0.0)
      throw std::runtime_error("EOM::LanczosDaggerSolve: negative beta^2");
    bj = std::sqrt(bj);

    if (bj < 1e-14 || j == max_iter - 1)
      break;

    hall(j, j + 1) = bj;
    hall(j + 1, j) = bj;
    basis.push_back(w / bj);
  }

  int m = j_final + 1;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, hall.submat(0, 0, m - 1, m - 1));
  int nret = std::min(state_want, (int)eigval.n_elem);
  arma::vec e = eigval.head(nret);

  index_t Q = vi.GetQSpaceOrbit();
  std::vector<Operator> ritz;
  for (int k = 0; k < nret; ++k)
  {
    Operator rk = MakeEmptyDagger(Q);
    for (int i = 0; i < m; ++i)
      rk += eigvec(i, k) * basis[i];
    double nr = particle_removed ? OverlapPR(rk, rk) : OverlapPA(rk, rk);
    if (nr > 0.0)
      rk = rk / std::sqrt(nr);
    ritz.push_back(std::move(rk));
  }

  std::cout << "Lanczos (dagger) finished with " << m << " steps" << std::endl;
  return {e, ritz};
}

EOM::RunResult EOM::RunPA(int max_iter, int state_want)
{
  index_t Q = PickQOrbitParticle(J2, parity, itz);
  Orbit &oQ = modelspace->GetOrbit(Q);
  std::cout << "1-PA-EOM  J2=" << J2 << " parity=" << parity << " Tz=" << itz
            << "  Q=" << Q << " (n=" << oQ.n << " l=" << oQ.l << " j=" << oQ.j2 / 2.
            << " tz=" << oQ.tz2 << ")" << std::endl;

  Operator v0 = PureAdagSeed(Q);
  double n0 = OverlapPA(v0, v0);
  std::cout << "  initial N-norm = " << n0 << "  (seed=a†_Q)" << std::endl;

  auto lr = LanczosDaggerSolve(v0, max_iter, state_want, false);
  double eref = Hs.ZeroBody;

  std::cout << "\n  1-PA-EOM eigenvalues (Lanczos):" << std::endl;
  std::cout << "  ZeroBody (eref) = " << eref << " MeV" << std::endl;
  for (int k = 0; k < (int)lr.first.n_elem; ++k)
    std::cout << "    ω(" << k << "): attachment=" << lr.first(k)
              << "  absolute=" << lr.first(k) + eref << " MeV" << std::endl;

  std::cout << "\n  Hs 1b SPE (particles in Q channel):" << std::endl;
  for (auto p : Hs.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (modelspace->holes.count(p))
      continue;
    std::cout << "    orbit " << p << ": n=" << modelspace->GetOrbit(p).n
              << "  ε=" << Hs.OneBody(p, p) << " MeV" << std::endl;
  }

  ArnoldiResult ar;
  ar.energies = lr.first;
  ar.ritz = std::move(lr.second);
  RunResult out{eref, ar, (int)Q, Hs.OneBody(Q, Q), 0.0};
  return out;
}

EOM::RunResult EOM::RunPR(int max_iter, int state_want)
{
  index_t Q = PickQOrbitHole(J2, parity, itz);
  Orbit &oQ = modelspace->GetOrbit(Q);
  std::cout << "1-PR-EOM  J2=" << J2 << " parity=" << parity << " Tz=" << itz
            << "  Q=" << Q << " (n=" << oQ.n << " l=" << oQ.l << " j=" << oQ.j2 / 2.
            << " tz=" << oQ.tz2 << ") [hole]" << std::endl;

  Operator R1 = MakeEmptyDagger(Q);
  R1.OneBody(Q, 0) = 1.0;
  Operator v0 = LadderPR(R1);
  double n0 = OverlapPR(v0, v0);
  std::cout << "  initial N-norm = " << n0 << "  (seed=a†_Q hole)" << std::endl;

  Operator w1 = HtcPR(R1);
  double ray = OverlapPR(w1, R1) / OverlapPR(R1, R1);
  double spe = Hs.OneBody(Q, Q);
  std::cout << "  Koopmans: SPE=" << spe << "  Rayleigh(-[H,a†])=" << ray
            << "  (want ≈ -SPE)" << std::endl;

  auto lr = LanczosDaggerSolve(v0, max_iter, state_want, true);
  double eref = Hs.ZeroBody;

  std::cout << "\n  1-PR-EOM eigenvalues (Lanczos):" << std::endl;
  std::cout << "  ZeroBody (eref) = " << eref << " MeV" << std::endl;
  for (int k = 0; k < (int)lr.first.n_elem; ++k)
    std::cout << "    ω(" << k << "): removal=" << lr.first(k)
              << "  E(A-1)=" << lr.first(k) + eref << " MeV" << std::endl;

  std::cout << "\n  Hs 1b SPE (holes in Q channel):" << std::endl;
  for (auto h : Hs.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
  {
    if (!modelspace->holes.count(h))
      continue;
    std::cout << "    orbit " << h << ": n=" << modelspace->GetOrbit(h).n
              << "  ε=" << Hs.OneBody(h, h) << " MeV" << std::endl;
  }

  ArnoldiResult ar;
  ar.energies = lr.first;
  ar.ritz = std::move(lr.second);
  RunResult out{eref, ar, (int)Q, spe, ray};
  return out;
}

// ============================================================
//  Run
//  Mirrors lines 133-173 of run/mr_eom.py:
//    ConstructConfigs / ConstructNormMatrix / ConstructProjectMatrix
//    compute eref, build random projected initial vector, run ArnoldiSolveH2.
// ============================================================

EOM::RunResult EOM::Run(int max_iter, int state_want)
{
  if (is_multiref)
    return RunMR(max_iter, state_want);
  switch (sr_mode)
  {
  case SREOMMode::ParticleAttached:
    return RunPA(max_iter, state_want);
  case SREOMMode::ParticleRemoved:
    return RunPR(max_iter, state_want);
  default:
    return RunSR(max_iter, state_want);
  }
}

// ---------------------------------------------------------------------------
// SR solve: mirrors run/sr_eom.py
//   force_decouple → random initial vector via GetVSEOM_ladder_single
//   → normalise → LanczosSolve(HtcSingle, NormSingle)
// ---------------------------------------------------------------------------
EOM::RunResult EOM::RunSR(int max_iter, int state_want)
{
//force_decouple(Hs);

  UnitTest unt(*modelspace);
  Operator h_rand = unt.RandomOp(*modelspace, J2, itz,parity, 2, 1);
  Operator chi    = GetVSEOM_ladder_single(h_rand, 1); 

  //double nm = NormSingle(chi, chi);
  double nm = GetVSEOM_Overlap_single(chi, chi);
  if (nm <= 0.0)
    throw std::runtime_error("EOM::RunSR: initial vector has zero norm");
  chi = chi / std::sqrt(nm);

  auto lr = LanczosSolve(chi, max_iter, state_want);
  arma::vec e               = lr.first;
  std::vector<Operator> ritz = std::move(lr.second);

  double eref = Hs.ZeroBody;

  std::cout << "\n  SR EOM eigenvalues (Lanczos):" << std::endl;
  std::cout << "  ZeroBody (eref) = " << eref << " MeV" << std::endl;
  for (int k = 0; k < (int)e.n_elem; ++k)
    std::cout << "    E(" << k << "): excitation=" << e(k)
              << "  absolute=" << e(k) + eref << " MeV" << std::endl;

  ArnoldiResult ar;
  ar.energies = e;
  ar.ritz     = std::move(ritz);
  return RunResult{eref, ar};
}

// ---------------------------------------------------------------------------
// MR solve: ConstructConfigs → NormMatrix → ProjectMatrix
//   → random projected initial vector → ArnoldiSolveH2 (full H)
// No EraseValence / E_val·N shift; eigenvalues are Delta E.
// ---------------------------------------------------------------------------
EOM::RunResult EOM::RunMR(int max_iter, int state_want)
{
  // --- (1) Prepare Hamiltonian (full decoupled H; keep vv) ---
  force_decouple(Hs);
  Hs_full = Hs;
  have_Hs_full = true;
  ClearReferenceEnergyShift();
  double eref = GetVSEOM_Overlap_multiref(Hs);
  std::cout << "  E_ref (<H>_RDM) = " << eref
            << "   ZeroBody = " << Hs.ZeroBody << " MeV" << std::endl;
  std::cout << "  Arnoldi on full H (no EraseValence / E_val·N shift)"
            << std::endl;

  // --- (2) Setup ---
  if (IsTensorEOM()) {
    ConstructConfigs_tensor();
    ConstructNormMatrix_tensor();
    if (arnoldi_use_h3)
      std::cout << "  tensor EOM H3: ethS 223_231_st + 223_232 + 223_132_tts"
                << std::endl;
    else
      std::cout << "  tensor EOM: H3 off (Hall = Htc only)" << std::endl;
  } else {
    ConstructConfigs();
    ConstructNormMatrix();
  }
  ConstructProjectMatrix();

  // --- (3) Random projected initial vector ---
  UnitTest unt(*modelspace);
  Operator h_rand = unt.RandomOp(*modelspace, J2, itz, parity, 2, 1);
  Operator chi_b  = GetVSEOM_ladder_multiref(h_rand, 1);
  if (arnoldi_use_projection)
    ProjectOprator(chi_b);
  double n0 = ComputeNorm(chi_b, chi_b);
  std::cout << "Initial vector N-norm after projection: " << n0 << std::endl;
  if (n0 <= 0.0)
    throw std::runtime_error("EOM::RunMR: initial vector has non-positive N-norm");

  // --- (4) Solve ---
  ArnoldiResult ar = ArnoldiSolveH2(chi_b, max_iter, state_want);

  // --- (5) Print summary ---
  std::cout << "\n  MR EOM eigenvalues (ArnoldiSolveH2):" << std::endl;
  std::cout << "  E_ref = " << eref << " MeV" << std::endl;
  for (int k = 0; k < (int)ar.energies.n_elem; ++k)
    std::cout << "    E(" << k << "): excitation=" << ar.energies(k)
              << "  absolute=" << ar.energies(k) + eref << " MeV" << std::endl;

  return RunResult{eref, ar};
}
