#include "EOM.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"
#include "Commutator.hh"
#include "IMSRG3Commutators.hh"
#include "FactorizedDoubleCommutator.hh"
#include "UnitTest.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
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
void EOM::ConstructConfigs() {
  // Generate configuration for fock space EOM
  // First ppvv
  std::cout << "Constructing EOM configurations for J2=" << J2
            << " parity=" << parity << " itz=" << itz << std::endl;
  // first we do one body
  qv_start = 0;
  qv_end = 0;
  qv_dim = 0;
  int norbits = modelspace->norbits;
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
  if (qv_dim > 0)
    qv_end = qv_start + qv_dim - 1;
std::cout << "dimension EOM qv: " << qv_start << " " << qv_end << std::endl;
  ph_start = eom_confs.size();
  ph_end = 0;
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
  if (ph_dim > 0)
    ph_end = ph_start + ph_dim - 1;
  std::cout << "dimension EOM ph: " << ph_start << " " << ph_end << std::endl;

  ppvv_start = eom_confs.size();
  ppvv_end = 0;
  size_t number_channels = modelspace->GetNumberTwoBodyChannels();
  for (index_t ich = 0; ich < number_channels; ich++) {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    for (auto &ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv())) {
      for (auto &iket : tbc.GetKetIndex_vv()) {
        eom_confs.push_back({ibra, iket, ich, eom_dims});
        eom_dims += 1;
        ppvv_dim += 1;
      }
    }
  }
  if (ppvv_dim > 0)
    ppvv_end = ppvv_start + ppvv_dim - 1;

  std::cout << "dimension EOM ppvv: " << ppvv_start << " " << ppvv_end
            << std::endl;
  // for (index_t i = ppvv_start; i < ppvv_start + ppvv_dim; i++) {
  //   auto &c = eom_confs.at(i);
  //   TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(c[2]);
  //   Ket &kbra = tbc.GetKet(c[0]);
  //   int cvq_p = modelspace->GetOrbit(kbra.p).cvq;
  //   int cvq_q = modelspace->GetOrbit(kbra.q).cvq;
  //   // vpvv: bra ket has one valence (1) and one particle (2)
  //   if ((cvq_p == 2 && cvq_q == 1) || (cvq_p == 1 && cvq_q == 2)) {
  //     Ket &kket = tbc.GetKet(c[1]);
  //     std::cout << "vpvvhere " << i
  //               << " " << kbra.p << " " << kbra.q
  //               << " " << kket.p << " " << kket.q
  //               << " " << tbc.J << std::endl;
  //   }
  // }

  pphv_start = eom_confs.size();
  pphv_end = 0;
  number_channels = modelspace->GetNumberTwoBodyChannels();
  for (index_t ich = 0; ich < number_channels; ich++) {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    for (auto &ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(),
                                  tbc.GetKetIndex_vv())) {
      for (auto &iket : tbc.GetKetIndex_vc()) {
        eom_confs.push_back({ibra, iket, ich, eom_dims});
        eom_dims += 1;
        pphv_dim += 1;
      }
    }
  }
  if (pphv_dim > 0)
    pphv_end = pphv_start + pphv_dim - 1;
  std::cout << "dimension EOM pphv: " << pphv_start << " " << pphv_end
            << std::endl;
  // for (index_t i = pphv_start; i < pphv_start + pphv_dim; i++) {
  //   auto &c = eom_confs.at(i);
  //   TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(c[2]);
  //   Ket &kbra = tbc.GetKet(c[0]);
  //   int cvq_p = modelspace->GetOrbit(kbra.p).cvq;
  //   int cvq_q = modelspace->GetOrbit(kbra.q).cvq;
  //   // vvhv: bra ket has two valence (1,1)
  //   if (cvq_p == 1 && cvq_q == 1) {
  //     Ket &kket = tbc.GetKet(c[1]);
  //     std::cout << "vvhv " << i
  //               << " " << kbra.p << " " << kbra.q
  //               << " " << kket.p << " " << kket.q
  //               << " " << c[2] << std::endl;
  //   }
  // }

  pphh_start = eom_confs.size();
  pphh_end = 0;
  number_channels = modelspace->GetNumberTwoBodyChannels();
  for (index_t ich = 0; ich < number_channels; ich++) {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    for (auto &ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(),
                                  tbc.GetKetIndex_vv())) {
      for (auto &iket : tbc.GetKetIndex_cc()) {
        eom_confs.push_back({ibra, iket, ich, eom_dims});
        eom_dims += 1;
        pphh_dim += 1;
      }
    }
  }
  if (pphh_dim > 0)
    pphh_end = pphh_start + pphh_dim - 1;

  std::cout << "dimension EOM pphh: " << pphh_start << " " << pphh_end
            << std::endl;

  std::cout << "w: " << eom_confs.size() << std::endl;
}

arma::vec EOM::GetEnergies() { return Energies; }

void EOM::ConstructNormMatrix() {

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
        Orbit &oc2 = modelspace->GetOrbit(c1);
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
        // if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<" " << Nkernel(i,j) << std::endl;

        // if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<"
        // "<<dbra.p<<dbra.q<<e1<<c1<<e2<<c2<<" " << Nkernel(i,j) << std::endl;
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


  // C1 threebody diagram
   
  if (ppvv_dim != 0) {

    for (index_t i = ppvv_start; i <= ppvv_end; i++) {
      std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
      Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
      size_t a1 = dbra1.p;
      size_t b1 = dbra1.q; 
      size_t c1 = dket1.p;
      size_t d1 = dket1.q;
      Orbit &ob1 = modelspace->GetOrbit(b1);

      for (index_t j = ppvv_start; j <= ppvv_end; j++) {
        std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
      Ket &dbra1 = tbc_ket.GetKet(cf_ket[0]);
      Ket &dket1 = tbc_ket.GetKet(cf_ket[1]);
      size_t a2 = dbra1.p;
      size_t b2 = dbra1.q; 
      size_t c2 = dket1.p;
      size_t d2 = dket1.q;
      Orbit &ob2 = modelspace->GetOrbit(b2);
      if(b1==b2 && ob1.cvq==1 && ob2.cvq==1){ // both valence not contributing to norm.
      double val = ThreeBody_Diagram(c1, d1, a2, a1, d2, c2, b1,
                                     tbc_bra.J, tbc_ket.J);
        Nkernel(i, j) += val;
    }
        
      }
    }
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

double EOM::ThreeBody_Diagram_Entries_Internal(size_t a, size_t b, size_t c,
                                               size_t d, size_t e, size_t f,
                                               size_t g, double j0, double j2)
{
  if (!rdm.ThreeBody.IsAllocated())
    return 0.0;

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

    value += std::sqrt(double(twoJ) + 1.0) * zme * rdm_me;
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
  // 2-body types: c[0]=ibra, c[1]=iket, c[2]=ich -> unpack to a,b,c,d
  auto print_2b = [&](const std::string &label, index_t start, index_t dim) {
    if (dim == 0) return;
    for (index_t i = 0; i < dim; i++) {
      auto &cfs = eom_confs.at(start + i);
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(cfs[2]);
      Ket &kbra = tbc.GetKet(cfs[0]);
      Ket &kket = tbc.GetKet(cfs[1]);
      std::cout << label << " " << (start + i)
                << " " << kbra.p << " " << kbra.q
                << " " << kket.p << " " << kket.q
                << " " << cfs[2] << std::endl;
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
    if (s_max < 1e-4) {
        Amat.zeros(n, n);
        return;
    }

    // Use a relative threshold to cleanly separate range from null space.
    // Vectors with s_i < 1e-6 * s_max are floating-point noise, not in the range.
    arma::uvec range_idx = arma::find(s >= 1e-3 * s_max);

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
    size_t ich = eom_confs.at(i)[2];
    v1(i) = Op1.GetTwoBody(ich, ich, ibra, iket);
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    v1(i) = Op1.GetTwoBody(ich, ich, ibra, iket);
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    v1(i) = Op1.GetTwoBody(ich, ich, ibra, iket);
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
    size_t ich = eom_confs.at(i)[2];
    v2(i) = Op2.GetTwoBody(ich, ich, ibra, iket);
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    v2(i) = Op2.GetTwoBody(ich, ich, ibra, iket);
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    v2(i) = Op2.GetTwoBody(ich, ich, ibra, iket);
  }

  // Compute v1.t() * Nkernel * v2
  arma::vec temp = Nkernel * v2;
  double norm = arma::dot(v1, temp);

  return norm;
}

void EOM::ProjectOprator(Operator &Qin) {
  arma::vec Qflat(eom_confs.size());
  Qflat.fill(0.);

  for (index_t i = qv_start; i <= qv_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    Qflat(i) = Qin.GetOneBody(p, q);
  }
  for (index_t i = ph_start; i <= ph_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    Qflat(i) = Qin.GetOneBody(p, q);
  }

  for (index_t i = pphh_start; i <= pphh_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qflat(i) = Qin.GetTwoBody(ich, ich, ibra, iket);
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qflat(i) = Qin.GetTwoBody(ich, ich, ibra, iket);
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qflat(i) = Qin.GetTwoBody(ich, ich, ibra, iket);
  }

  Qin *= 0.;

  arma::vec Qout = Prj_kernel * Qflat;

  for (index_t i = qv_start; i <= qv_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    Qin.SetOneBody(p, q, Qout(i));
  }
  for (index_t i = ph_start; i <= ph_end; i++) {
    size_t p = eom_confs.at(i)[0];
    size_t q = eom_confs.at(i)[1];
    Qin.SetOneBody(p, q, Qout(i));
  }

  for (index_t i = pphh_start; i <= pphh_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qin.TwoBody.SetTBME(ich, ich, ibra, iket, Qout(i));
  }

  for (index_t i = ppvv_start; i <= ppvv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qin.TwoBody.SetTBME(ich, ich, ibra, iket, Qout(i));
  }
  for (index_t i = pphv_start; i <= pphv_end; i++) {
    size_t ibra = eom_confs.at(i)[0];
    size_t iket = eom_confs.at(i)[1];
    size_t ich = eom_confs.at(i)[2];
    Qin.TwoBody.SetTBME(ich, ich, ibra, iket, Qout(i));
  }
}

void EOM::block_svd(std::vector<int> &coupled_vector) {

  arma::mat Amat;
  int n = coupled_vector.size();
  Amat.set_size(n, n);
  Amat.fill(0.);
  for (index_t i = 0; i < n; i++) {
    index_t ia = coupled_vector.at(i);
    std::array<index_t, 4> &cfs_bra = eom_confs.at(ia);
    for (index_t j = 0; j < n; j++) {
      index_t ib = coupled_vector.at(j);
      std::array<index_t, 4> &cfs_ket = eom_confs.at(ib);
      Amat(i, j) = Nkernel(cfs_bra[3], cfs_ket[3]);
    }
  }

  SqrtMat(Amat, n);

  for (index_t i = 0; i < n; i++) {
    std::array<index_t, 4> &cfs_bra = eom_confs.at(coupled_vector.at(i));
    for (index_t j = 0; j < n; j++) {
      std::array<index_t, 4> &cfs_ket = eom_confs.at(coupled_vector.at(j));
      Prj_kernel(cfs_bra[3], cfs_ket[3]) = Amat(i, j);
    }
  }
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

  // One body: copy ph and vq blocks, applying (anti-)hermitian phase.
  for (auto &i : H.modelspace->core) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      Orbit &oa = H.modelspace->GetOrbit(a);
      Hod.SetOneBody(a,i,H.OneBody(a, i));
    }
  }

  // Two body: copy pp-hh, pp-vh and pp-vv blocks; (anti-)hermitian partner added automatically.
  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
    arma::mat &H2 = iter.second;

    // diagonal and off-diagonal channels treated the same for pp-hh block
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

    // in a <J1|O^{LM}|J2> block, we have pphh and hhpp, hhpp-> pphh in conjugate channel->hhpp
      for (auto &ibra : tbc_bra.GetKetIndex_cc()) {
        for (auto &iket :
             VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                         tbc_ket.GetKetIndex_qv())) {
          Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket)*herm_phase*hZ);
        }
      }

    }
  }

  return Hod;
}

double EOM::GetVSEOM_Overlap_single(Operator &H1, Operator &H2) {
  double ovlp = 0;
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
  std::cout<< "onebody norm: "<< ovlp<< std::endl;

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

  // One body piece -- eliminate ph bits
  //
  for (auto &i : H.modelspace->core) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      Orbit &oa = H.modelspace->GetOrbit(a);
      Hod.SetOneBody(a,i,H.OneBody(a, i));
    }
  }

  // Two body piece only stored half channel, no need to change
  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

    arma::mat &H2 = iter.second;

    // diagonal channel
    // add <ab|ij> and <ij||ab> will be automatically added
    // with proper factor from (anti)hermitian
    if (ch_bra == ch_ket) {

      for (auto &iket :
           tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the
                                     // reference state)
      {
        for (auto &ibra :
             VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                         tbc_bra.GetKetIndex_qv())) {
          Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
        }
      }
    }

    if (ch_bra != ch_ket) {

      // off-diagonal channel

      for (auto &iket :
           tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the
                                     // reference state)
      {
        for (auto &ibra :
             VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                         tbc_bra.GetKetIndex_qv())) {
          Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
        }
      }

      for (auto &ibra :
           tbc_bra.GetKetIndex_cc()) // cc means core-core ('holes' refer to the
                                     // reference state)
      {
        for (auto &iket :
             VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                         tbc_ket.GetKetIndex_qv())) {
          Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket,
                                H2(ibra, iket) * hZ * herm_phase);
        }
      }
    }
  }

  for (auto &i : H.modelspace->valence) {
    for (auto &a : H.modelspace->qspace) {
      Hod.SetOneBody(a,i,H.OneBody(a, i));
    }
  }


  // Two body piece only stored half channel, no need to change
  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
    arma::mat &H2 = iter.second;

    // PPvc: valence-core pairs
    for (auto &iket : tbc_ket.GetKetIndex_vc())
    {
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                       tbc_bra.GetKetIndex_qv())) {
        Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
      }
    }

    // PPvv: valence-valence pairs
    for (auto &iket : tbc_ket.GetKetIndex_vv())
    {
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv())) {
        Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
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

  auto dabc = [](const Ket3 &ket) -> double {
    int d_ab = (ket.p == ket.q) ? 1 : 0;
    int d_bc = (ket.q == ket.r) ? 1 : 0;
    int d_ac = (ket.p == ket.r) ? 1 : 0;
    int g_abc = 1 + d_ab + d_bc + d_ac + 2 * (d_ab * d_bc);
    return std::sqrt(static_cast<double>(g_abc));
  };

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

        ovlp2 += H2(ibra, iket) *
                 RdmTB_J(tbc_bra.J, dbra.p, dbra.q, dket.p, dket.q) *
                 std::sqrt(2.0 * tbc_bra.J + 1.0);
      }
    }
  }

  if (H.ThreeBody.IsAllocated() && rdm.ThreeBody.IsAllocated()) {
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

          double norm_denom = dabc(bra) * dabc(ket);
          ovlp3 += h3 * r3 * std::sqrt(double(twoJ) + 1.0)
                   / (norm_denom * norm_denom);
        }
      }
    }
  }
  

  return (ovlp + ovlp1 + ovlp2 + ovlp3);
}

// ============================================================
//  Lanczos / Arnoldi helpers  (translated from run/lanczos.py)
// ============================================================

// -----------  norm helpers  ---------------------------------

/// Thin wrapper around GetVSEOM_Overlap_single.
double EOM::NormSingle(Operator &T1, Operator &T2)
{
 if(T1.IsReduced()){// for tensors, the wave function operator is reduced

 return GetVSEOM_Overlap_single(T1, T2);
 }
 else{// scaler is not reduced, then we use commutator
 Operator T1d = EOM::GetVSEOM_ladder_single(T1,1);
 Operator nop=T1*0.0;
 nop.SetHermitian();
 Commutator::comm110ss(T1d,T2,nop);
 Commutator::comm220ss(T1d,T2,nop);
 return(nop.ZeroBody/2.);}

 
}

/// <T1|T2> using the multiref metric:
///   D† = GetVSEOM_ladder_multiref(T1, 0)
///   nop = [D†, T2]
///   result = GetVSEOM_Overlap_multiref(nop) / 2
double EOM::NormMultiref(Operator &T1, Operator &T2)
{
  Operator T1d = GetVSEOM_ladder_multiref(T1, 1);
  Operator nop = T1*0.;
  nop.SetHermitian();
  nop = Commutator::Commutator(T1d, T2);
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
  t3.SetHermitian();

  Operator nop = t1 * 0.0;
  nop.SetHermitian();

  Commutator::comm223ss(t2, haml, t3);
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
  return GetVSEOM_ladder_single(ht_plus, -1);
}

/// Multiref action: [Haml, chi], then apply ladder (herm=1) and project.
Operator EOM::HtcMultiref(Operator &haml, Operator &chi)
{
  Operator ht_minus = Commutator::Commutator(haml, chi);
  Operator heom = GetVSEOM_ladder_multiref(ht_minus, -1);
  ProjectOprator(heom);
  return heom;
}

// -----------  double commutator diagonal  -------------------

/// Diagonal double-commutator contribution:
///   opa += comm223_231 + comm223_232 + comm223_132
///   returns { GetVSEOM_Overlap_multiref(opa)/2, opa }
std::pair<double, Operator> EOM::DcomMultiref(Operator &haml, Operator &chi)
{
  Operator opa = chi * 0.0;
  opa.SetHermitian();

  Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
  Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(true);
  Commutator::FactorizedDoubleCommutator::comm223_231(chi, haml, opa);
  Commutator::FactorizedDoubleCommutator::comm223_232(chi, haml, opa);
  Commutator::FactorizedDoubleCommutator::comm223_132(chi, haml, opa);

  double rst = GetVSEOM_Overlap_multiref(opa);
  return {rst / 2.0, opa};
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
  double nn = NormSingle(vi, vi);
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

    double ai = NormSingle(w, lanczos_vector[j]);
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
        double cij = NormSingle(w, lanczos_vector[i]);
        w = w - cij * lanczos_vector[i];
      }
    }

    double nm = NormSingle(w, w);
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

EOM::ArnoldiResult
EOM::ArnoldiSolve(Operator &vi, int max_iter, int state_want)
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
  double nn = NormMultiref(vi, vi);
  Operator v0 = vi / std::sqrt(nn);
  lanczos_vector.push_back(v0);
  h2_diag.push_back(DcomMultiref(Hs, v0).first);
  const double cn0 = ComputeNorm(v0, v0);

  arma::vec e(state_want, arma::fill::zeros);
  arma::mat vs(1, 1, arma::fill::zeros);
  arma::vec prev_e;
  bool      prev_e_valid = false;
  int       j_final      = 0;

  for (int j = 0; j < max_iter - 1; ++j)
  {
    j_final = j;

    // H1 * v_j  (raw, before projection — used for matrix elements)
    Operator h1v_j = HtcMultiref(Hs, lanczos_vector[j]);
    h1v_cache.push_back(h1v_j);

    // fill row / column j of Hamiltonian matrix using polarization identity
    for (int i = 0; i <= j; ++i)
    {
      double h1ij   = NormMultiref(lanczos_vector[i], h1v_j);
      double h1ji   = NormMultiref(lanczos_vector[j], h1v_cache[i]);
      double h1_sym = 0.5 * (h1ij + h1ji);

      Operator v_sum  = lanczos_vector[i] + lanczos_vector[j];
      double h2_cross = DcomMultiref(Hs, v_sum).first;
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
        double cij = NormMultiref(lanczos_vector[i], w);
        w = w - cij * lanczos_vector[i];
      }
      ProjectOprator(w);
    }

    double bj        = NormMultiref(w, w);
    double bj_kernel = ComputeNorm(w, w);
    std::cout<< bj<<" "<< bj_kernel<<" kernels "<<std::endl;

    // null-space breakdown: w exhausted the physical subspace
    if (std::abs(bj_kernel) < null_tol * cn0)
    {
      std::cout << "arnoldi: null vector (breakdown) at step " << j + 1
                << " (ComputeNorm=" << bj_kernel << "), stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim-1, dim-1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // exact breakdown: w == 0
    if (std::abs(bj) < bj_tol)
    {
      std::cout << "arnoldi: exact breakdown at step " << j + 1 << ", stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim-1, dim-1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // negative norm: indefinite metric hit — solve and stop
    if (bj < 0.0)
    {
      std::cout << "arnoldi: bj=" << bj << " < 0 at step " << j + 1
                << " (indefinite metric), stopping." << std::endl;
      int dim = j + 1;
      arma::mat sub = hall.submat(0, 0, dim-1, dim-1);
      arma::vec ev; arma::mat evec;
      arma::eig_sym(ev, evec, sub);
      e  = ev.head(std::min(state_want, (int)ev.n_elem));
      vs = evec;
      break;
    }

    // normal step: add new basis vector
    Operator new_vec = w / std::sqrt(bj);
    lanczos_vector.push_back(new_vec);
    h2_diag.push_back(DcomMultiref(Hs, new_vec).first);

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
        std::cout << "arnoldi eigenvalues @ step " << j + 1 << ": ";
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
          std::cout << "Arnoldi converged at step " << j + 1 << std::endl;
          break;
        }
      }
      prev_e       = e;
      prev_e_valid = true;
    }
  }

  // Final eigensolution on the full accumulated subspace.
  {
    int nb_cur = (int)lanczos_vector.size();
    arma::mat sub = hall.submat(0, 0, nb_cur - 1, nb_cur - 1);
    arma::vec eigval_f; arma::mat eigvec_f;
    arma::eig_sym(eigval_f, eigvec_f, sub);
    e  = eigval_f.head(std::min(state_want, (int)eigval_f.n_elem));
    vs = eigvec_f;
  }

  std::cout << "Arnoldi finished with " << j_final + 1 << " steps" << std::endl;
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
  result.hall     = hall.submat(0, 0, nb - 1, nb - 1);
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
      std::cout << "hello" << std::endl;
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
//  Run
//  Mirrors lines 133-173 of run/mr_eom.py:
//    ConstructConfigs / ConstructNormMatrix / ConstructProjectMatrix
//    compute eref, build random projected initial vector, run ArnoldiSolve.
// ============================================================

EOM::RunResult EOM::Run(int max_iter, int state_want)
{
  if (is_multiref)
    return RunMR(max_iter, state_want);
  else
    return RunSR(max_iter, state_want);
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
  Operator chi    = GetVSEOM_ladder_single(h_rand, -1); 

  double nm = NormSingle(chi, chi);
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
// MR solve: mirrors run/mr_eom.py lines 133-173
//   ConstructConfigs → NormMatrix → ProjectMatrix
//   → random projected initial vector → ArnoldiSolve (MR)
// ---------------------------------------------------------------------------
EOM::RunResult EOM::RunMR(int max_iter, int state_want)
{
  // --- (1) Setup ---
  ConstructConfigs();
  ConstructNormMatrix();
  ConstructProjectMatrix();

  // --- (2) Reference energy = <Hs> in the valence space ---
  double eref = GetVSEOM_Overlap_multiref(Hs);
  std::cout << "  E_ref (valence) = " << eref-Hs.ZeroBody
            << "   ZeroBody = " << Hs.ZeroBody
            << "   E_ref total = " << eref << " MeV" << std::endl;

  // --- (3) Random projected initial vector ---
  UnitTest unt(*modelspace);
  Operator h_rand = unt.RandomOp(*modelspace, 0, 0, 0, 2, 1); // now we can only deal with scaler
  Operator chi_b  = GetVSEOM_ladder_multiref(h_rand, -1);
  ProjectOprator(chi_b);

  // --- (4) Solve ---
  ArnoldiResult ar = ArnoldiSolve(chi_b, max_iter, state_want);

  // --- (5) Print summary ---
  std::cout << "\n  MR EOM eigenvalues (Arnoldi):" << std::endl;
  std::cout << "  E_ref = " << eref << " MeV" << std::endl;
  for (int k = 0; k < (int)ar.energies.n_elem; ++k)
    std::cout << "    E(" << k << "): excitation=" << ar.energies(k)
              << "  absolute=" << ar.energies(k) + eref << " MeV" << std::endl;

  return RunResult{eref, ar};
}
