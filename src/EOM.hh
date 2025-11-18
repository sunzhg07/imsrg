//////////////////////////////////////////////////////////////////////////////////
//    EOM.hh, part of  imsrg++
//    Copyright (C) 2021  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef EOM_hh
#define EOM_hh 1


#include "ModelSpace.hh"
#include "Operator.hh"
#include "Generator.hh"


class EOM
{
 public:
  // Fields
  ModelSpace * modelspace; ///< Pointer to the associated modelspace
  Operator H;
  Operator rdm;
  arma::sp_mat Nkernel; 
  arma::sp_mat Prj_kernel; 
  //arma::mat Nkernel;  
  arma::imat eom_configs;  
  arma::vec Energies;

  // configurations
  std::vector<std::array<index_t, 4>> eom_confs;
  index_t qv_start,qv_end,qv_dim;
  index_t ph_start,ph_end,ph_dim;
  index_t ppvv_start,ppvv_end,ppvv_dim;
  index_t pphv_start,pphv_end,pphv_dim;
  index_t pphh_start,pphh_end,pphh_dim;
  size_t channel;
  index_t eom_dims=0;


  // Methods
  EOM(); ///< Default constructor
  EOM(ModelSpace& ms, Operator& rdm);
  EOM(Operator& H);

  void ConstructConfigs();
  void PrintConfigs();
  void Setup_rdm();
  void ConstructNormMatrix();
  void ConstructProjectMatrix();
  void SolveEOM();
  double Core_Diagram(size_t a, size_t b,size_t c,size_t d,size_t e,size_t f,double j1, double j2 );
  arma::vec GetEnergies();
  void SqrtMat(arma::mat& Amat, size_t n);
  void ProjectOprator(Operator & Qin);

};



#endif
