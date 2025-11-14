#include "EOM.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"
using PhysConst::SQRT2;


EOM::EOM(ModelSpace& ms, Operator &rdm)
 : modelspace(&ms), rdm(rdm)
{
    eom_dims=0;
qv_dim=0;
ph_dim=0;
ppvv_dim=0;
pphv_dim=0;
pphh_dim=0;
};

EOM::EOM(Operator& h)
 : modelspace(h.modelspace), H(h) 
{
eom_dims=0;
qv_dim=0;
ph_dim=0;
ppvv_dim=0;
pphv_dim=0;
pphh_dim=0;
};


///  In case we want to construct the A matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void EOM::ConstructConfigs()
{
    // Generate configuration for fock space EOM
    // First ppvv
 
// first we do one body
qv_start=0;
qv_end =0;
    int norbits = modelspace->norbits;
   for (index_t i_orb = 0; i_orb < norbits; i_orb++)
   {
      Orbit &oi = modelspace->GetOrbit(i_orb);

       if(oi.cvq !=2 )continue;
     
   for (index_t j_orb = 0; j_orb < norbits; j_orb++){
      Orbit &oj = modelspace->GetOrbit(j_orb);
       if(oj.cvq !=1 )continue;
      if(oj.l != oi.l)continue;
      if(oj.j2 != oi.j2)continue;
      if(oj.tz2 != oi.tz2)continue;
       eom_confs.push_back({i_orb,j_orb,0,eom_dims});
       eom_dims+=1;
       qv_dim+=1;
   }}
qv_end=qv_start+qv_dim-1;
    std::cout << "dimension EOM qv: "<< qv_start <<" "<<qv_end << std::endl;

   

ph_start = eom_confs.size();
   for (index_t i_orb = 0; i_orb < norbits; i_orb++)
   {
      Orbit &oi = modelspace->GetOrbit(i_orb);
       if(oi.cvq ==0 )continue;
   for (index_t j_orb = 0; j_orb < norbits; j_orb++){
      Orbit &oj = modelspace->GetOrbit(j_orb);
       if(oj.cvq !=0 )continue;
      if(oj.l != oi.l)continue;
      if(oj.j2 != oi.j2)continue;
      if(oj.tz2 != oi.tz2)continue;
       eom_confs.push_back({i_orb,j_orb,0,eom_dims});
       eom_dims+=1;
       ph_dim+=1;
       
   }}
   ph_end = ph_start+ph_dim-1;
    std::cout << "dimension EOM ph: "<< ph_start <<" "<<ph_end << std::endl;
   
   ppvv_start=eom_confs.size();
    
    size_t number_channels = modelspace->GetNumberTwoBodyChannels();
    for (index_t ich =0; ich < number_channels; ich ++){
        TwoBodyChannel &tbc = modelspace ->GetTwoBodyChannel(ich);
        for( auto &ibra:  VectorUnion(tbc.GetKetIndex_qq(),tbc.GetKetIndex_qv())){
        for( auto &iket:  tbc.GetKetIndex_vv()){
        eom_confs.push_back({ibra,iket,ich,eom_dims});
       eom_dims+=1;
       ppvv_dim+=1;
    }}}
   ppvv_end=ppvv_start+ppvv_dim-1;

    std::cout << "dimension EOM ppvv: "<< ppvv_start <<" "<<ppvv_end << std::endl;

   pphv_start=eom_confs.size();
    
    number_channels = modelspace->GetNumberTwoBodyChannels();
    for (index_t ich =0; ich < number_channels; ich ++){
        TwoBodyChannel &tbc = modelspace ->GetTwoBodyChannel(ich);
        for( auto &ibra:  VectorUnion(tbc.GetKetIndex_qq(),tbc.GetKetIndex_qv(),tbc.GetKetIndex_vv())){
        for( auto &iket:  tbc.GetKetIndex_vc()){
           Ket& dbra=tbc.GetKet(ibra);
           Ket& dket=tbc.GetKet(iket);
//        std::cout <<" Config pphv "<< eom_dims<<" "<<dbra.p<<dbra.q<<dket.p<<dket.q<<" "<<ibra<<" "<<iket<<" "<<ich<< " "<<tbc.J << std::endl;
        eom_confs.push_back({ibra,iket,ich,eom_dims});
       eom_dims+=1;
       pphv_dim+=1;
    }}}
   pphv_end=pphv_start+pphv_dim-1;
    std::cout << "dimension EOM pphv: "<< pphv_start <<" "<<pphv_end << std::endl;


   pphh_start=eom_confs.size();
    
    number_channels = modelspace->GetNumberTwoBodyChannels();
    for (index_t ich =0; ich < number_channels; ich ++){
        TwoBodyChannel &tbc = modelspace ->GetTwoBodyChannel(ich);
        for( auto &ibra:  VectorUnion(tbc.GetKetIndex_qq(),tbc.GetKetIndex_qv(),tbc.GetKetIndex_vv())){
        for( auto &iket:  tbc.GetKetIndex_cc()){
           Ket& dbra=tbc.GetKet(ibra);
           Ket& dket=tbc.GetKet(iket);
        //std::cout <<" Config pphh "<< eom_dims<<" "<<dbra.p<<dbra.q<<dket.p<<dket.q<<" "<<ibra<<" "<<iket<<" "<<ich<< " "<<tbc.J << std::endl;
        eom_confs.push_back({ibra,iket,ich,eom_dims});
       eom_dims+=1;
       pphh_dim+=1;
    }}}
   pphh_end=pphh_start+pphh_dim-1;

    std::cout << "dimension EOM pphh: "<< pphh_start <<" "<<pphh_end << std::endl;

    std::cout << "dimension EOM all: "<< eom_confs.size() << std::endl;


}



arma::vec EOM::GetEnergies()
{
   return Energies;
}

void EOM::ConstructNormMatrix()
{

    
    Nkernel.set_size(eom_dims,eom_dims);
    Nkernel.zeros();
   // B4, becnhmarked with srg
   if(qv_dim !=0){
       for (index_t i =qv_start; i <= qv_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);

       for (index_t j =qv_start; j <= qv_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
               if(cf_bra[0] != cf_ket[0])continue;
        Orbit &obra = modelspace->GetOrbit(cf_bra[1]) ;  
        Nkernel(i,j)+=rdm.OneBody(cf_bra[1],cf_ket[1])*sqrt(obra.j2+1.);
       }
   }
   }

  //A1, B5, becnhmarked with srg
   if(ph_dim !=0){
       for (index_t i =ph_start; i <= ph_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
            Orbit &obra = modelspace->GetOrbit(cf_bra[1]);
       for (index_t j =ph_start; j <= ph_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
        
     if(i==j){
        Nkernel(i,j) +=1.*(obra.j2+1.);
     } 
     if(cf_bra[1] != cf_ket[1])continue;
        Nkernel(i,j) -=rdm.OneBody(cf_ket[0],cf_bra[0])*sqrt(obra.j2+1.);
       }

   }
   }
  //C1, benchmarked, maybe, it too many
  if(ppvv_dim!=0){

       for (index_t i =ppvv_start; i <= ppvv_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);

       for (index_t j =ppvv_start; j <= ppvv_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);

               if( cf_bra[2] !=cf_ket[2] )continue;
               if(cf_bra[0]!=cf_ket[0]  )continue;
        TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
        double val = rdm.GetTwoBody(cf_bra[2],cf_ket[2],cf_bra[1],cf_ket[1])*sqrt(2*tbc_bra.J+1.);

        Nkernel(i,j) =val;
   //     if(abs(Nkernel(i,j))>0.0000001)std::cout<<i<<" "<<j<<" "<<Nkernel(i,j)<<" cpp"<<std::endl;

       }}
   }
  //B1, benchmarked
  if(pphv_dim!=0){

       for (index_t i =pphv_start; i <= pphv_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
       for (index_t j =pphv_start; j <= pphv_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
           if(cf_bra[2]!=cf_ket[2])continue; // must be in the same channel
           if(cf_bra[0]!=cf_ket[0])continue; // pp=pp

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
            Ket &dbra = tbc_bra.GetKet(cf_bra[0]);
            Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
            Ket &dket2 = tbc_bra.GetKet(cf_ket[1]);
            size_t e1=dket1.p; size_t c1=dket1.q;
            size_t e2=dket2.p; size_t c2=dket2.q;
            if(e1!=e2)continue;
            Orbit &oc1=modelspace->GetOrbit(c1);
            Orbit &oc2=modelspace->GetOrbit(c1);
            if(oc1.l !=oc2.l)continue;
            if(oc1.j2 !=oc2.j2)continue;
            if(oc1.tz2 !=oc2.tz2)continue;
            if(oc1.cvq != 1)continue; // need be valence
            if(oc2.cvq != 1)continue;
            double j1=tbc_bra.J;
            
            Nkernel(i,j) =rdm.OneBody(c1,c2)*(2*tbc_bra.J+1.)/sqrt(oc1.j2+1.);
            //if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<" " << Nkernel(i,j) << std::endl;
    
            //if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<" "<<dbra.p<<dbra.q<<e1<<c1<<e2<<c2<<" " << Nkernel(i,j) << std::endl;
       }}
  }
// C4
  if(pphv_dim!=0){

       for (index_t i =pphv_start; i <= pphv_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
            Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
            Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
            size_t e1=dbra1.q; size_t c=dbra1.p; // permute ec to make c the valence
            size_t f1=dket1.p; size_t a=dket1.q;

            Orbit &oa=modelspace->GetOrbit(a);
            Orbit &oc=modelspace->GetOrbit(c);
            Orbit &oe1=modelspace->GetOrbit(e1);
            Orbit &of1=modelspace->GetOrbit(f1);
            if(oe1.cvq !=1 && oc.cvq !=1)continue;
            double j1=tbc_bra.J;
            double norm_fact1 =1.;
            if(e1==c) norm_fact1 =sqrt(2.);
       for (index_t j =pphv_start; j <= pphv_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
            Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
            Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
            size_t e2=dbra2.q; size_t b=dbra2.p; // permute to make b the valence
            size_t f2=dket2.p; size_t d=dket2.q;
            Orbit &od=modelspace->GetOrbit(d);
            Orbit &ob=modelspace->GetOrbit(b);
            Orbit &oe2=modelspace->GetOrbit(e2);
            Orbit &of2=modelspace->GetOrbit(f2);
            if(oe2.cvq !=1 && ob.cvq !=1)continue;
            double j2=tbc_ket.J;
            double norm_fact2 =1.;
            if(e2==b) norm_fact2 =sqrt(2.);
            double norm_fact = norm_fact1*norm_fact2;
            if(f1 !=f2)continue;

            if(e1==e2){
            size_t jmin = std::max(abs(oa.j2-ob.j2)/2, abs(oc.j2-od.j2)/2);
            size_t jmax = std::min(abs(oa.j2+ob.j2)/2, abs(oc.j2+od.j2)/2);
            size_t p_bra = (oa.l+ob.l,2)%2;
            size_t p_ket =  (oc.l+od.l,2)%2;
            if(p_bra !=p_ket)continue;
            int itz_bra= (oa.tz2 + ob.tz2)/2;
            int itz_ket= (oc.tz2 + od.tz2)/2;
            if(itz_bra!=itz_ket)continue;
            
            for(auto jtot = jmin; jtot <= jmax; jtot+=1){
           // double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1*2,oe1.j2, oc.j2,of1.j2, j2*2,od.j2, oa.j2, ob.j2, 2*jtot);
            double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1,of1.j2*0.5, oa.j2*0.5,oe1.j2*0.5, j2,ob.j2*0.5, oc.j2*0.5, od.j2*0.5, jtot);
            double val2= norm_fact*rdm.TwoBody.GetTBME_J_norm(jtot, jtot,a,b,c,d)*sqrt(2*jtot+1.);
          Nkernel(i,j) +=val*val2*dbra2.Phase(j2)*dket1.Phase(j1)*(pow(-1,oa.j2*0.5+ob.j2*0.5+oc.j2*0.5+od.j2*0.5));
            
            }
        }
        if(e1!=c && c==e2){
            size_t jmin = std::max(abs(oa.j2-ob.j2)/2, abs(oe1.j2-od.j2)/2);
            size_t jmax = std::min(abs(oa.j2+ob.j2)/2, abs(oe1.j2+od.j2)/2);
            size_t p_bra = (oa.l+ob.l,2)%2;
            size_t p_ket =  (oe1.l+od.l,2)%2;
            if(p_bra !=p_ket)continue;
            int itz_bra= (oa.tz2 + ob.tz2)/2;
            int itz_ket= (oe1.tz2 + od.tz2)/2;
            if(itz_bra!=itz_ket)continue;
            
            for(auto jtot = jmin; jtot <= jmax; jtot+=1){
           // double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1*2,oe1.j2, oc.j2,of1.j2, j2*2,od.j2, oa.j2, ob.j2, 2*jtot);
            double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1,of1.j2*0.5, oa.j2*0.5,oc.j2*0.5, j2,ob.j2*0.5, oe1.j2*0.5, od.j2*0.5, jtot);
            double val2= norm_fact*rdm.TwoBody.GetTBME_J_norm(jtot, jtot,a,b,e1,d)*sqrt(2*jtot+1.);
          Nkernel(i,j) +=val*val2*dbra2.Phase(j2)*dket1.Phase(j1)*dbra1.Phase(j1)*(pow(-1,oa.j2*0.5+ob.j2*0.5+oe1.j2*0.5+od.j2*0.5));
            
            }
        }

        if(e1==b && e2!=b){
            size_t jmin = std::max(abs(oa.j2-oe2.j2)/2, abs(oc.j2-od.j2)/2);
            size_t jmax = std::min(abs(oa.j2+oe2.j2)/2, abs(oc.j2+od.j2)/2);
            size_t p_bra = (oa.l+oe2.l,2)%2;
            size_t p_ket =  (oc.l+od.l,2)%2;
            if(p_bra !=p_ket)continue;
            int itz_bra= (oa.tz2 + oe2.tz2)/2;
            int itz_ket= (oc.tz2 + od.tz2)/2;
            if(itz_bra!=itz_ket)continue;
            
            for(auto jtot = jmin; jtot <= jmax; jtot+=1){
           // double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1*2,oe1.j2, oc.j2,of1.j2, j2*2,od.j2, oa.j2, ob.j2, 2*jtot);
            double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1,of1.j2*0.5, oa.j2*0.5,oe1.j2*0.5, j2,oe2.j2*0.5, oc.j2*0.5, od.j2*0.5, jtot);
            double val2= norm_fact*rdm.TwoBody.GetTBME_J_norm(jtot, jtot,a,e2,c,d)*sqrt(2*jtot+1.);
          Nkernel(i,j) +=val*val2*dket1.Phase(j1)*(pow(-1,oa.j2*0.5+oe2.j2*0.5+oc.j2*0.5+od.j2*0.5));
            
            }
        }

        if(c==b && c!=e1 && b!=e2){
            size_t jmin = std::max(abs(oa.j2-oe2.j2)/2, abs(oe1.j2-od.j2)/2);
            size_t jmax = std::min(abs(oa.j2+oe2.j2)/2, abs(oe1.j2+od.j2)/2);
            size_t p_bra = (oa.l+oe2.l,2)%2;
            size_t p_ket =  (oe1.l+od.l,2)%2;
            if(p_bra !=p_ket)continue;
            int itz_bra= (oa.tz2 + oe2.tz2)/2;
            int itz_ket= (oe1.tz2 + od.tz2)/2;
            if(itz_bra!=itz_ket)continue;
            
            for(auto jtot = jmin; jtot <= jmax; jtot+=1){
           // double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1*2,oe1.j2, oc.j2,of1.j2, j2*2,od.j2, oa.j2, ob.j2, 2*jtot);
            double val=(2.*j1+1.)*(2.*j2+1.)*AngMom::NineJ(j1,of1.j2*0.5, oa.j2*0.5,oc.j2*0.5, j2,oe2.j2*0.5, oe1.j2*0.5, od.j2*0.5, jtot);
            double val2= norm_fact*rdm.TwoBody.GetTBME_J_norm(jtot, jtot,a,e2,e1,d)*sqrt(2*jtot+1.);
          Nkernel(i,j) +=val*val2*dbra1.Phase(j1)*dket1.Phase(j1)*(pow(-1,oa.j2*0.5+oe2.j2*0.5+oe1.j2*0.5+od.j2*0.5));
            
            }
        }

//if(abs(Nkernel(i,j))>0.0000001)std::cout<<i<<" "<<j<<" "<<c<<" "<<e1<<" "<<a<<" "<<f1<<" "<<e2<<" "<<b<<" "<<f2<<" "<<d<<" "<<\
//    " "<<Nkernel(i,j)<<std::endl;



            }

       }
   }
// A2 C3 benchmarked
  if(pphh_dim !=0){
    // first We compute A2
    for (index_t i =pphh_start; i <= pphh_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
      Nkernel(i,i) += (2*tbc_bra.J+1.);
  }

    //C3
    for (index_t i =pphh_start; i <= pphh_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
           TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
            Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
            Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
            size_t a=dbra1.p; size_t b=dbra1.q;
            Orbit &oa=modelspace->GetOrbit(a);
            Orbit &ob=modelspace->GetOrbit(b);
            if(oa.cvq != 1)continue; // need be valence
            if(ob.cvq != 1)continue;

    for (index_t j =pphh_start; j <= pphh_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
           if(cf_ket[2] != cf_bra[2])continue;
           if(cf_ket[1] != cf_bra[1])continue;
           TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
            Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
            Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
            size_t c =dbra2.p; size_t d=dbra2.q;
            Orbit &oc=modelspace->GetOrbit(c);
            Orbit &od=modelspace->GetOrbit(d);
            if(oc.cvq != 1)continue; // need be valence
            if(od.cvq != 1)continue;

            double val2= rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J,a,b,c,d)*sqrt(2*tbc_bra.J+1.);
            Nkernel(i,j) += val2;

     // if(abs(Nkernel(i,j))>0.00001)std::cout << i<<" "<<j<<" "<<Nkernel(i,j)<<std::endl;
  }

  }
  }

//B2, benchmarked
    for (index_t i =pphh_start; i <= pphh_end; i++){
           std::array<index_t, 4> &cf_bra = eom_confs.at(i);
           TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);
            Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
            Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
            size_t a=dbra1.p; size_t b=dbra1.q;
            Orbit &oa=modelspace->GetOrbit(a);
            Orbit &ob=modelspace->GetOrbit(b);
            if(oa.cvq != 1 && ob.cvq!=1)continue; // at least one need be valence
            double norm_fact1 = 1.;
            if(a==b) norm_fact1 = sqrt(2.);
    for (index_t j =pphh_start; j <= pphh_end; j++){
           std::array<index_t, 4> &cf_ket = eom_confs.at(j);
           if(cf_ket[2] != cf_bra[2])continue; // same channel
           if(cf_ket[1] != cf_bra[1])continue; // same holes
           TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(cf_ket[2]);
            Ket &dbra2 = tbc_ket.GetKet(cf_ket[0]);
            Ket &dket2 = tbc_ket.GetKet(cf_ket[1]);
            size_t c =dbra2.p; size_t d=dbra2.q;
            double norm_fact2 = 1.;
            if(c==d) norm_fact2 = sqrt(2.);

            Orbit &oc=modelspace->GetOrbit(c);
            Orbit &od=modelspace->GetOrbit(d);

           // std::cout<<" start "<< std::endl;
            double norm_fact=norm_fact1*norm_fact2;
            if(b==d){
            Nkernel(i,j)-=norm_fact*rdm.OneBody(c,a)*(2*tbc_bra.J+1.)/sqrt(oc.j2+1.);
       // std::cout<<" a cont "<< std::endl;
    }
            
            if(b!=a && ob.cvq==1 && a==d){
            Nkernel(i,j)-=norm_fact*rdm.OneBody(c,b)*(2*tbc_bra.J+1.)\
            *dbra1.Phase(tbc_bra.J)/sqrt(oc.j2+1.);
         //   std::cout<<" b cont"<< std::endl;
        }
            if( c!=d && od.cvq==1 && b==c){
            Nkernel(i,j)-=norm_fact*rdm.OneBody(d,a)*(2*tbc_bra.J+1.)\
            *dbra2.Phase(tbc_ket.J)/sqrt(od.j2+1.);
       // std::cout<<" c cont"<< std::endl;
    }
            if(b!=a && c!=d && od.cvq==1 && ob.cvq==1  && a==c){
            Nkernel(i,j)-=norm_fact*rdm.OneBody(d,b)*(2*tbc_bra.J+1.)\
            *dbra2.Phase(tbc_ket.J)*dbra1.Phase(tbc_bra.J)/sqrt(od.j2+1.);
       // std::cout<<" d cont"<< std::endl;
    }
    //if(abs(Nkernel(i,j))>0.00000001)std::cout<< i<<" " <<j<<" " << Nkernel(i,j) << std::endl;
     // if(abs(Nkernel(i,j))>0.001)std::cout << cf_bra[2]<<cf_bra[0]<<cf_bra[1]<<" "<<cf_ket[2]<<cf_ket[0]<<cf_ket[1]<<Nkernel(i,j)<<" "<<tbc_bra.J <<" "<<(2*tbc_bra.J+1.)/(oc.j2+1.)<<std::endl;
  }

  }

    // now off diagonal B3 benchmarked
 if(pphv_dim!=0 && ph_dim!=0){

      for (index_t i =pphv_start; i <= pphv_end; i++){
          std::array<index_t, 4> &cf_bra = eom_confs.at(i);
           TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

           Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
           Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
           size_t a=dbra1.p; size_t b=dbra1.q; // permute ec to make c the valence
           size_t c=dket1.p; size_t d=dket1.q;
           double norm_fact = 1;
           if(a==b)norm_fact=sqrt(2.);

           Orbit &oa=modelspace->GetOrbit(a);
           Orbit &oc=modelspace->GetOrbit(c);
           Orbit &ob=modelspace->GetOrbit(b);
           Orbit &od=modelspace->GetOrbit(d);
           double j1=tbc_bra.J;
           if(oa.cvq !=1 && ob.cvq !=1)continue;
           if(od.cvq !=1)continue;
           for (index_t j=ph_start; j<= ph_end; j++){
          std::array<index_t, 4> &cf_ket = eom_confs.at(j);
               size_t c1= cf_ket[0]; size_t b1=cf_ket[1];
               if(b1!=c)continue;

               if(c1== b ){
                   Nkernel(i,j)+=norm_fact*dbra1.Phase(j1)*rdm.OneBody(a,d)*(2*j1+1.)/sqrt(oa.j2+1.);
               }
               if(c1== a && a!=b){
                   Nkernel(i,j)+=norm_fact*rdm.OneBody(b,d)*(2*j1+1.)/sqrt(ob.j2+1.);
               }
          //     if(abs(Nkernel(i,j))>0.001)std::cout <<i<<" "<<j << " "<<Nkernel(i,j)<<std::endl;
           
           }
       }
 }
//
//
//
 if(pphv_dim!=0 && ph_dim!=0){

      for (index_t i =pphv_start; i <= pphv_end; i++){
          std::array<index_t, 4> &cf_bra = eom_confs.at(i);
           TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

           Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
           Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
           size_t a=dbra1.p; size_t b=dbra1.q; // permute ec to make c the valence
           size_t c=dket1.p; size_t d=dket1.q;
           Orbit &oa=modelspace->GetOrbit(a);
           Orbit &oc=modelspace->GetOrbit(c);
           Orbit &ob=modelspace->GetOrbit(b);
           Orbit &od=modelspace->GetOrbit(d);

           double j1=tbc_bra.J;
           if(oa.cvq !=1 || ob.cvq !=1 || od.cvq !=1)continue;
           for (index_t j=ph_start; j<= ph_end; j++){
          std::array<index_t, 4> &cf_ket = eom_confs.at(j);
               size_t c1= cf_ket[0]; size_t b1=cf_ket[1];
               if( b1==c){
                   double norm_fact =1.;
                   if(c1==d)norm_fact=sqrt(2.);
                   double val2= rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J,a,b,c1,d)*sqrt(j1*2.+1.);
                   Nkernel(i,j)-=val2*norm_fact;
               }
   //            if(abs(Nkernel(i,j))>0.000001) std::cout <<i<<" "<<j<<" " << Nkernel(i,j)<<std::endl;
           
           }
       }
 }

//// c2
if(ppvv_dim!=0 && qv_dim!=0){

      for (index_t i =ppvv_start; i <= ppvv_end; i++){
          std::array<index_t, 4> &cf_bra = eom_confs.at(i);
           TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(cf_bra[2]);

           Ket &dbra1 = tbc_bra.GetKet(cf_bra[0]);
           Ket &dket1 = tbc_bra.GetKet(cf_bra[1]);
           size_t a=dbra1.p; size_t b=dbra1.q; // permute ec to make c the valence
           size_t c=dket1.p; size_t d=dket1.q;
           Orbit &oa=modelspace->GetOrbit(a);
           Orbit &oc=modelspace->GetOrbit(c);
           Orbit &ob=modelspace->GetOrbit(b);
           Orbit &od=modelspace->GetOrbit(d);

           double j1=tbc_bra.J;
           if(oa.cvq !=1 || oc.cvq !=1 || od.cvq !=1 )continue;
           for (index_t j=qv_start; j<= qv_end; j++){
          std::array<index_t, 4> &cf_ket = eom_confs.at(j);
               size_t c1= cf_ket[0]; size_t b1=cf_ket[1];
               if( b==c1){
                   double norm_fact=1.;
                   if(a==b1)norm_fact=sqrt(2);
                   double val2= rdm.TwoBody.GetTBME_J_norm(tbc_bra.J, tbc_bra.J,a,b1,c,d)*sqrt(2*j1+1.);
                   Nkernel(i,j)+=val2*norm_fact;
               }
               if(abs(Nkernel(i,j))>0.000001) std::cout <<i<<" "<<j<<" " << Nkernel(i,j)<<std::endl;

           }
       }
 }
  //all done
  }



void EOM::PrintConfigs()
{
    for (std::array<index_t, 4> &cfs : eom_confs){
        std::cout<< cfs[0] <<" "<<cfs[1]<<" "<<cfs[2]<<" "<<cfs[3]<<std::endl;
    }
}

void EOM::ConstructHamiltonianMatrix()
{
}

void EOM::SolveEOM()
{
}

void EOM::Setup_rdm()
{
}
