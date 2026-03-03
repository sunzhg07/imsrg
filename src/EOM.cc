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
EOM::EOM(Operator &Hs, Operator &rdm, int J2, int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), rdm(rdm), J2(J2), parity(parity),
      itz(itz), is_multiref(true) {
  eom_dims = 0;
  qv_dim = 0;
  ph_dim = 0;
  ppvv_dim = 0;
  pphv_dim = 0;
  pphh_dim = 0;
};

EOM::EOM(Operator &Hs, const std::string &tdm_file, int J2, int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), rdm(*Hs.modelspace), J2(J2), parity(parity),
      itz(itz), is_multiref(true) {
  eom_dims = 0;
  qv_dim = 0;
  ph_dim = 0;
  ppvv_dim = 0;
  pphv_dim = 0;
  pphh_dim = 0;
  rdm = ReadTdm(tdm_file);
};

EOM::EOM(Operator &Hs, int J2, int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), J2(J2), parity(parity),
      itz(itz), is_multiref(false) {};

///  In case we want to construct the A matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void EOM::ConstructConfigs() {
  // Generate configuration for fock space EOM
  // First ppvv

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

  std::cout << "dimension EOM all: " << eom_confs.size() << std::endl;
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
        Nkernel(i, j) += rdm.OneBody(cf_bra[1], cf_ket[1]) * sqrt(obra.j2 + 1.);
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
        Nkernel(i, j) -= rdm.OneBody(cf_ket[0], cf_bra[0]) * sqrt(obra.j2 + 1.);
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
        double val =
            rdm.GetTwoBody(cf_bra[2], cf_ket[2], cf_bra[1], cf_ket[1]) *
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
            rdm.OneBody(c1, c2) * (2 * tbc_bra.J + 1.) / sqrt(oc1.j2 + 1.);
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
            rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J, a, b, c, d) *
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
        Nkernel(i, j) -= norm_fact * rdm.OneBody(c, a) * (2 * tbc_bra.J + 1.) /
                         sqrt(oc.j2 + 1.);
      }

      if (b != a && ob.cvq == 1 && a == d) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(c, b) * (2 * tbc_bra.J + 1.) *
                         dbra1.Phase(tbc_bra.J) / sqrt(oc.j2 + 1.);
      }
      if (c != d && od.cvq == 1 && b == c) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(d, a) * (2 * tbc_bra.J + 1.) *
                         dbra2.Phase(tbc_ket.J) / sqrt(od.j2 + 1.);
      }
      if (b != a && c != d && od.cvq == 1 && ob.cvq == 1 && a == c) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(d, b) * (2 * tbc_bra.J + 1.) *
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
          Nkernel(i, j) += norm_fact * dbra1.Phase(j1) * rdm.OneBody(a, d) *
                           (2 * j1 + 1.) / sqrt(oa.j2 + 1.);
          Nkernel(j, i) += norm_fact * dbra1.Phase(j1) * rdm.OneBody(a, d) *
                           (2 * j1 + 1.) / sqrt(oa.j2 + 1.);
        }
        if (c1 == a && a != b) {
          Nkernel(i, j) +=
              norm_fact * rdm.OneBody(b, d) * (2 * j1 + 1.) / sqrt(ob.j2 + 1.);
          Nkernel(j, i) +=
              norm_fact * rdm.OneBody(b, d) * (2 * j1 + 1.) / sqrt(ob.j2 + 1.);
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
              rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J, a, b, c1, d) *
              sqrt(j1 * 2. + 1.);
          Nkernel(i, j) -= val2 * norm_fact;
          Nkernel(j, i) -= val2 * norm_fact;
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
              rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J, a, b1, c, d) *
              sqrt(2 * j1 + 1.);
          Nkernel(i, j) += val2 * norm_fact;
          Nkernel(j, i) += val2 * norm_fact;
        }
      }
    }
  }
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
    double val2 = rdm.TwoBody.GetTBME_J_norm(jtot, jtot, a, b, c, d) *
                  sqrt(2 * jtot + 1.);
    val += norm_fact * val1 * val2 *
           (pow(-1, oa.j2 * 0.5 + ob.j2 * 0.5 + oc.j2 * 0.5 + od.j2 * 0.5));
  }
  return (val);
}
void EOM::PrintConfigs() {
  for (std::array<index_t, 4> &cfs : eom_confs) {
    std::cout << cfs[0] << " " << cfs[1] << " " << cfs[2] << " " << cfs[3]
              << std::endl;
  }
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
    if (s_max < 1e-10) {
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

  if (herm == 0) {
    Hod.SetHermitian();
    herm_phase = 1;
  }
  if (herm == 1) {
    Hod.SetAntiHermitian();
    herm_phase = -1;
  }

  // One body: copy ph and vq blocks, applying (anti-)hermitian phase.
  for (auto &i : H.modelspace->core) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      Orbit &oa = H.modelspace->GetOrbit(a);
      Hod.OneBody(a, i) = H.OneBody(a, i);
      int phase_ia = H.modelspace->phase((oa.j2 - oi.j2) / 2);
      Hod.OneBody(i, a) = H.OneBody(a, i) * phase_ia * herm_phase;
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

  return Hod;
}

double EOM::GetVSEOM_Overlap_single(Operator &H1, Operator &H2) {
  double ovlp = 0;
  // One-body contribution: ph and vq blocks
  for (auto &i : H1.modelspace->holes) {
    for (auto &a : VectorUnion(H1.modelspace->valence, H1.modelspace->qspace)) {
      ovlp += H1.OneBody(a, i) * H2.OneBody(a, i);
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
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket);
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
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket);
        }
      }

      for (auto &iket :
           VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(),
                       tbc_ket.GetKetIndex_qv())) {
        for (auto &ibra :
             tbc_bra.GetKetIndex_cc()) // cc means core-core ('holes' refer to
                                       // the reference state)
        {
          ovlp += Hmat1(ibra, iket) * Hmat2(ibra, iket);
        }
      }
    }
  }
  return ovlp;
}

Operator EOM::GetVSEOM_ladder_multiref(Operator &H, int herm) {

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of
  // H int herm, 0 for hermit, and 1 for antihermit
  int hZ = H.IsHermitian() ? +1 : -1;
  Operator Hod = 0.0 * H;

  int herm_phase = -1;

  if (herm == 0) {
    Hod.SetHermitian();
    herm_phase = 1;
  }
  if (herm == 1) {
    Hod.SetAntiHermitian();
    herm_phase = -1;
  }

  // One body piece -- eliminate ph bits
  //
  for (auto &i : H.modelspace->core) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      Orbit &oa = H.modelspace->GetOrbit(a);
      Hod.OneBody(a, i) = H.OneBody(a, i);
      int phase_ia = H.modelspace->phase((oa.j2 - oi.j2) / 2);
      Hod.OneBody(i, a) = H.OneBody(a, i) * phase_ia * herm_phase;
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
      Hod.OneBody(a, i) = H.OneBody(a, i);
      Hod.OneBody(i, a) = H.OneBody(a, i) * herm_phase;
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

  ovlp += H.ZeroBody;

  for (auto &i : H.modelspace->valence) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &j : H.modelspace->valence) {
      ovlp1 += H.OneBody(i, j) * rdm.OneBody(i, j) * sqrt((oi.j2 + 1));
    }
  }

  for (auto &iter : H.TwoBody.MatEl) {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1];
    TwoBodyChannel &tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
    arma::mat &H2 = iter.second;
    arma::mat &r2 = rdm.TwoBody.GetMatrix(ch_bra, ch_ket);

    for (auto &iket : tbc_ket.GetKetIndex_vv())
    {

      Ket &dket = tbc_ket.GetKet(iket);

      for (auto &ibra : tbc_bra.GetKetIndex_vv()) {
        Ket &dbra = tbc_bra.GetKet(ibra);

        double norm_fact = 1.0;
        if (dbra.p == dbra.q)
          norm_fact *= sqrt(2.);
        if (dket.p == dket.q)
          norm_fact *= sqrt(2.);

        ovlp2 += H2(ibra, iket) * r2(ibra, iket) * sqrt(2 * tbc_bra.J + 1.);
      }
    }
  }

  return (ovlp + ovlp1 + ovlp2);
}

// ============================================================
//  Lanczos / Arnoldi helpers  (translated from run/lanczos.py)
// ============================================================

// -----------  norm helpers  ---------------------------------

/// Thin wrapper around GetVSEOM_Overlap_single.
double EOM::NormSingle(Operator &T1, Operator &T2)
{
// Operator T1d = EOM::GetVSEOM_ladder_single(T1,0);
// Operator nop=T1*0.0;
// nop.SetHermitian();
// Commutator::comm110ss(T1d,T2,nop);
// Commutator::comm220ss(T1d,T2,nop);
// 
// return(nop.ZeroBody/2.);

 
  return GetVSEOM_Overlap_single(T1, T2);
}

/// <T1|T2> using the multiref metric:
///   D† = GetVSEOM_ladder_multiref(T1, 0)
///   nop = [D†, T2]
///   result = GetVSEOM_Overlap_multiref(nop) / 2
double EOM::NormMultiref(Operator &T1, Operator &T2)
{
  Operator T1d = GetVSEOM_ladder_multiref(T1, 0);
  Operator nop = Commutator::Commutator(T1d, T2);
  nop.SetHermitian();
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
  return GetVSEOM_ladder_single(ht_plus, 1);
}

/// Multiref action: [Haml, chi], then apply ladder (herm=1) and project.
Operator EOM::HtcMultiref(Operator &haml, Operator &chi)
{
  Operator ht_minus = Commutator::Commutator(haml, chi);
  Operator heom = GetVSEOM_ladder_multiref(ht_minus, 1);
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
  if (!converged)
  {
    int dim = std::min((int)lanczos_vector.size(), max_iter);
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
  const double bj_tol   = 1e-10;
  const double null_tol = 1e-6;
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
  // Create a scalar 2-body operator (Jrank=0, Trank=0, parity=0, particle_rank=2)
  Operator ops(*modelspace, 0, 0, 0, 2);
  ops *= 0.0;

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

      // ob_idx[obs] = {orbit_index, l, tz2}
      struct OrbEntry { int idx, l, tz2; };
      std::vector<OrbEntry> ob_idx(norb);

      for (int obs = 0; obs < norb; ++obs)
      {
        std::istringstream sl(lines[lidx++]);
        int dummy, nn, ll, jj, tt;
        sl >> dummy >> nn >> ll >> jj >> tt;
        ob_idx[obs].idx = (int)modelspace->GetOrbitIndex(nn, ll, jj, tt);
        ob_idx[obs].l   = ll;
        ob_idx[obs].tz2 = tt;
      }

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

          int jkl = std::stoi(tok[5]);   // same column as jij in this format
          int pkl = (ob_idx[ic].l + ob_idx[id].l) % 2;
          int tkl = (ob_idx[ic].tz2 + ob_idx[id].tz2) / 2;

          double rd = std::stod(tok.back()) / factor;
          ops.SetTwoBody(jij, pij, tij, jkl, pkl, tkl, aa, bb, cc, dd, rd);
        }
      }

      // --- three-body density matrix elements ---
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

          int aa = ob_idx[std::stoi(tok[1]) - 1].idx;
          int bb = ob_idx[std::stoi(tok[2]) - 1].idx;
          int cc = ob_idx[std::stoi(tok[3]) - 1].idx;
          int ee = ob_idx[std::stoi(tok[4]) - 1].idx;
          int ff = ob_idx[std::stoi(tok[5]) - 1].idx;
          int kk = ob_idx[std::stoi(tok[6]) - 1].idx;

          int jab  = std::stoi(tok[7]) / 2;
          int jef  = std::stoi(tok[8]) / 2;
          int jtot = std::stoi(tok[9]);

          double rd = std::stod(tok.back()) / factor;
          ops.ThreeBody.SetME_pn(jab, jef, jtot, aa, bb, cc, ee, ff, kk, (ThreeBME_type)rd);
        }
      }
    } // norb scope
  } // jtotal / factor scope

  return ops;
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
  Operator h_rand = unt.RandomOp(*modelspace, J2, parity, itz, 2, 1);
  Operator chi    = GetVSEOM_ladder_single(h_rand, 1); 

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
  std::cout << "  E_ref (valence) = " << eref
            << "   ZeroBody = " << Hs.ZeroBody
            << "   E_ref total = " << eref << " MeV" << std::endl;

  // --- (3) Random projected initial vector ---
  UnitTest unt(*modelspace);
  Operator h_rand = unt.RandomOp(*modelspace, 0, 0, 0, 2, 1); // now we can only deal with scaler
  Operator chi_b  = GetVSEOM_ladder_multiref(h_rand, 1);
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
