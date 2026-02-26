#include "EOM.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"
using PhysConst::SQRT2;

// we have two constructor for the EOM, w/o the rdm for multi-reference and
// single reference
EOM::EOM(Operator &Hs, Operator &rdm, int J2, int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), rdm(rdm), J2(J2), parity(parity),
      itz(itz) {
  eom_dims = 0;
  qv_dim = 0;
  ph_dim = 0;
  ppvv_dim = 0;
  pphv_dim = 0;
  pphh_dim = 0;
};

EOM::EOM(Operator &Hs, int J2, int parity, int itz)
    : modelspace(Hs.modelspace), Hs(Hs), J2(J2), parity(parity), itz(itz) {};

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
//std::cout <<"qv conf "<< eom_dims<< " "<< i_orb<<" "<<j_orb<<std::endl;
      eom_dims += 1;
      qv_dim += 1;
    }
  }
  if (qv_dim > 0)
    qv_end = qv_start + qv_dim - 1;
  //std::cout << "dimension EOM qv: " << qv_start << " " << qv_end << std::endl;

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
     // std::cout <<"ph conf "<< eom_dims<< "with i,j = "<<i_orb<<" "<<j_orb<<std::endl;
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
   // std::cout<<"ppvv conf"<< eom_dims<< "ich, bra,ket = "<< ich<< " " << ibra<< " " <<iket<<std::endl;
        eom_dims += 1;
        ppvv_dim += 1;
      }
    }
    // std::cout << "Dimension of ppvv at channel "<< tbc.J<<" "<<tbc.parity <<
    // " "<<tbc.Tz<<" is: "<< ppvv_dim<<std::endl;
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
  //      std::cout <<" Config pphv "<< eom_dims<<"ich bra ket =  "<<ich<< " " <<ibra<<" "<<iket << std::endl;
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
     //   std::cout <<" Config pphh "<< eom_dims<<"ich bra ket =  "<<ich<< " " <<ibra<<" "<<iket << std::endl;
        eom_dims += 1;
        pphh_dim += 1;
      }
    }
    // std::cout<< "dim pphh "<< eom_dims<<std::endl;
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
  std::cout << "done dimension EOM ph: " << ph_start << " " << ph_end
            << std::endl;
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
//        Ket &dbra1 = tbc_bra.GetKet(cf_bra[1]);
//        Ket &dket1 = tbc_bra.GetKet(cf_ket[1]);
//std::cout<< "ppvv norm"<< i<<" "<<j<< " is  "<< dbra1.p<<" "<<dbra1.q<<" "<<dket1.p<<" " <<dket1.q<< " "<< tbc_bra.J<<std::endl;
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
      // for (index_t i =4290; i <= 4290; i++){
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
        // for (index_t j =4049; j <= 4049; j++){
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
          // std::cout<< "a run " <<
          // norm_fact*Core_Diagram(d1,b2,a2,d2,a2,c1,j1,j2)*dket1.Phase(j1) <<
          // std::endl;
        }
        if (a1 == a2 && ob1.cvq == 1 && ob2.cvq == 1 && a1 != b1) {
          val += norm_fact * Core_Diagram(d1, b2, b1, d2, a2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra1.Phase(j1);
          // std::cout<< "b run " <<
          // norm_fact*Core_Diagram(d1,b2,b1,d2,a2,c1,j1,j2)*dket1.Phase(j1)*dbra1.Phase(j1)
          // << std::endl;
        }
        if (b1 == b2 && oa1.cvq == 1 && oa2.cvq == 1 && a2 != b2) {
          val += norm_fact * Core_Diagram(d1, a2, a1, d2, b2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra2.Phase(j2);
          // std::cout<< "c run " <<
          // norm_fact*Core_Diagram(d1,a2,a1,d2,b2,c1,j1,j2)*dket1.Phase(j1)*dbra2.Phase(j2)
          // << std::endl;
        }
        if (a1 == b2 && ob1.cvq == 1 && oa2.cvq == 1 && a1 != b1 && a2 != b2) {
          val += norm_fact * Core_Diagram(d1, a2, b1, d2, b2, c1, j1, j2) *
                 dket1.Phase(j1) * dbra1.Phase(j1) * dbra2.Phase(j2);
          // std::cout<< "d run " <<
          // norm_fact*Core_Diagram(d1,a2,b1,d2,b2,c1,j1,j2)*dket1.Phase(j1)*dbra1.Phase(j1)*dbra2.Phase(j2)
          // << std::endl;
        }
        Nkernel(i, j) += val;
        // std::cout<<a1<<" "<<b1<<" "<<c1<<" "<<d1<<" "<<a2<<" "<<b2<<"
        // "<<c2<<" "<<d2<<" "<<Nkernel(i,j)<<std::endl;
        // if(abs(Nkernel(i,j))>0.0000001)std::cout<<i<<" "<<j<<"
        // "<<Nkernel(i,j)<<std::endl;
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

        // if(abs(Nkernel(i,j))>0.00001)std::cout << i<<" "<<j<<"
        // "<<Nkernel(i,j)<<std::endl;
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
        // std::cout<<" a cont "<< std::endl;
      }

      if (b != a && ob.cvq == 1 && a == d) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(c, b) * (2 * tbc_bra.J + 1.) *
                         dbra1.Phase(tbc_bra.J) / sqrt(oc.j2 + 1.);
        //   std::cout<<" b cont"<< std::endl;
      }
      if (c != d && od.cvq == 1 && b == c) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(d, a) * (2 * tbc_bra.J + 1.) *
                         dbra2.Phase(tbc_ket.J) / sqrt(od.j2 + 1.);
        // std::cout<<" c cont"<< std::endl;
      }
      if (b != a && c != d && od.cvq == 1 && ob.cvq == 1 && a == c) {
        Nkernel(i, j) -= norm_fact * rdm.OneBody(d, b) * (2 * tbc_bra.J + 1.) *
                         dbra2.Phase(tbc_ket.J) * dbra1.Phase(tbc_bra.J) /
                         sqrt(od.j2 + 1.);
        // std::cout<<" d cont"<< std::endl;
      }
      // if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<" " <<
      // Nkernel(i,j) << std::endl;
      //  if(abs(Nkernel(i,j))>0.001)std::cout <<
      //  cf_bra[2]<<cf_bra[0]<<cf_bra[1]<<"
      //  "<<cf_ket[2]<<cf_ket[0]<<cf_ket[1]<<Nkernel(i,j)<<" "<<tbc_bra.J <<"
      //  "<<(2*tbc_bra.J+1.)/(oc.j2+1.)<<std::endl;
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
        //     if(abs(Nkernel(i,j))>0.001)std::cout <<i<<" "<<j << "
        //     "<<Nkernel(i,j)<<std::endl;
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
        //            if(abs(Nkernel(i,j))>0.000001) std::cout <<i<<" "<<j<<" "
        //            << Nkernel(i,j)<<std::endl;
      }
    }
  }

  //// c2
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
        // if(abs(Nkernel(i,j))>0.000001) std::cout <<i<<" "<<j<<" " <<
        // Nkernel(i,j)<<std::endl;
      }
    }
  }
  //  all done

  // std::cout<< "Number of nonzero matrix elements in norm kernel: " <<
  // Nkernel.n_nonzero << std::endl; std::cout << "\nNon-zero elements of A:" <<
  // std::endl;

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

  // first we do pphh, it have no coupling to other diagrams, and no coupling
  // bettwen channels
//  std::cout << " svd for pphh : " << std::endl;

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

  // ppvv
  // pv is coupled to all vqvv, but qqvv is not coupled
  //
  // we do vqvv together with onebody qv first
  std::vector<int> coupled_idx;
 // std::cout << " svd for vqvv and pv : " << std::endl;
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

//  std::cout << " svd for vqvv and pv a : " << qv_start << " " << qv_end
//            << std::endl;

  for (index_t i = qv_start; i <= qv_end; i++)
    coupled_idx.push_back(i);
//  std::cout << " svd for vqvv and pv b : " << std::endl;

  if (coupled_idx.size() != 0)
    block_svd(coupled_idx);

  // now all qqvv channels
  //
 // std::cout << " svd for qqvv : " << std::endl;
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

  // now ph and pphv,
  // vvhv and vh, sort by hole line
//  std::cout << " svd for vvhv and vh : " << std::endl;
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

    // std::cout << " vvhv part 1 dim in New Channel: " << n<< std::endl;
  }

  // vphv and qh
//  std::cout << " svd for vpvh and qh : " << std::endl;
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
    // std::cout << " vphv dim in New Channel: " << n<< std::endl;
  }

  // pphv itself, no coupling between channels
  //
  //
//  std::cout << " svd for qqhv : " << std::endl;
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
    if (s_max < 1e-14) {
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
 // std::cout << "New svd with block size: " << n << std::endl;
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

void EOM::SolveEOM() {}

void EOM::Setup_rdm() {}

// we enforce the decoupling condition for numerical stability
void EOM::force_decouple(Operator &H) {

  // One body piece -- eliminate ph bits
  //  particle hole excitation
  //
  for (auto &i : H.modelspace->core) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : VectorUnion(H.modelspace->valence, H.modelspace->qspace)) {
      Orbit &oa = H.modelspace->GetOrbit(a);
      H.OneBody(a, i) = 0.;
      H.OneBody(i, a) = 0.;
    }
  }

  for (auto &i : H.modelspace->valence) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : H.modelspace->qspace) {
      Orbit &oa = H.modelspace->GetOrbit(a);

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

  int herm_phase = 0;
  Operator Hod = 0.0 * H;

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

    arma::mat &Hod2 = Hod.TwoBody.GetMatrix(ch_bra, ch_ket);

    // diagonal channel
    // diagnal channel, add <ab|ij> and <ij||ab> will be automatically added
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
          // Hod2(ibra,iket)= H2(ibra,iket);
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
          // Hod2(ibra,iket)= H2(ibra,iket)*hZ*herm_phase;
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
  // One body piece -- eliminate ph bits
  for (auto &i : H1.modelspace->holes) {

    Orbit &oi = H1.modelspace->GetOrbit(i);

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
 // std::cout<< "input herm "<< herm<<" and herm phase is: "<< herm_phase<< std::endl;


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

    arma::mat &Hod2 = Hod.TwoBody.GetMatrix(ch_bra, ch_ket);

    // diagonal channel
    // diagnal channel, add <ab|ij> and <ij||ab> will be automatically added
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
          // Hod2(ibra,iket)= H2(ibra,iket);
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
          // Hod2(ibra,iket)= H2(ibra,iket)*hZ*herm_phase;
          Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket,
                                H2(ibra, iket) * hZ * herm_phase);
        }
      }
    }
  }

  for (auto &i : H.modelspace->valence) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : H.modelspace->qspace) {
      Orbit &oa = H.modelspace->GetOrbit(a);

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
    arma::mat &Hod2 = Hod.TwoBody.GetMatrix(ch_bra, ch_ket);

    // PPvc
    //
    for (auto &iket : tbc_ket.GetKetIndex_vc()) // cc means core-core ('holes'
                                                // refer to the reference state)
    {
      Ket &dket = tbc_ket.GetKet(iket);
      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(),
                       tbc_bra.GetKetIndex_qv())) {
        Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
      }
    }

    // PPvv
    //
    for (auto &iket : tbc_ket.GetKetIndex_vv()) // cc means core-core ('holes'
                                                // refer to the reference state)
    {
      Ket &dket = tbc_ket.GetKet(iket);

      for (auto &ibra :
           VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv())) {
        Hod.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, H2(ibra, iket));
      }
    }
  }

  return Hod;
}

double EOM::GetVSEOM_Overlap_multiref(Operator &H) {
  // [S+S', S-S']=[S',S] - [S,S']=2[S',S]
  double ovlp = 0;
  double ovlp1 = 0;
  double ovlp2 = 0;

  ovlp += H.ZeroBody;

  for (auto &i : H.modelspace->valence) {
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &j : H.modelspace->valence) {
      Orbit &oj = H.modelspace->GetOrbit(j);

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

//std::cout << " overlap, 0b, 1b,2b : "<<ovlp <<" " << ovlp1 << " " <<ovlp2<< std::endl;

  return (ovlp + ovlp1 + ovlp2);
}
