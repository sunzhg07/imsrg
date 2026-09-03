///////////////////////////////////////////////////////////////////////////////////
//    FactorizedDoubleCommutator.hh, part of  imsrg++
//    Copyright (C) 2023 Bingcheng He and Ragnar Stroberg
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

#ifndef FactorizedDoubleCommutator_hh
#define FactorizedDoubleCommutator_hh 1

#include "Operator.hh"
#include <vector>

namespace Commutator {

namespace FactorizedDoubleCommutator {
extern bool use_goose_tank_1b;
extern bool use_goose_tank_2b;

extern bool use_1b_intermediates;
extern bool use_2b_intermediates;

extern bool use_goose_tank_only_1b;
extern bool use_goose_tank_only_2b;
extern bool use_TypeI_1b;    // 231 χ1b Diagram I  / f^I
extern bool use_TypeII_1b;   // 231 χ1b Diagram III / f^II
extern bool use_TypeIIIa_1b; // 231 χ2b II_a,II_c / f^III_a
extern bool use_TypeIII_1b;  // 231 χ2b II_b,II_d / f^III_b
extern bool use_TypeII_2b;
extern bool use_TypeIII_2b;
extern bool use_TypeGI_2b;    // 232 χ1b CHI_I  / Γ^I
extern bool use_TypeGII_2b;   // 232 χ1b CHI_II / Γ^II
extern bool use_TypeGIIIa_2b; // 232 IIa,IIc  χ^η ladder
extern bool use_TypeGIIIb_2b; // 232 IIb,IId  χ^η RC
extern bool use_TypeGIIIc_2b; // 232 IIe,IIf  χ^θ
extern bool use_TypeGIVa_2b;  // 232 IIIc,IIId χ^κ
extern bool use_TypeGIVb_2b;  // 232 χ^ι
extern bool use_TypeGIVc_2b;  // 232 IIIe,IIIf χ^λ

extern bool use__GT_TypeI_2b;
extern bool use__GT_TypeIV_2b;

void SetUse_GooseTank_1b(bool tf);
void SetUse_GooseTank_2b(bool tf);
void SetUse_1b_Intermediates(bool tf);
void SetUse_2b_Intermediates(bool tf);

void SetUse_GooseTank_only_1b(bool tf);
void SetUse_GooseTank_only_2b(bool tf);
void SetUse_TypeI_1b(bool tf);
void SetUse_TypeII_1b(bool tf);
void SetUse_TypeIIIa_1b(bool tf);
void SetUse_TypeIII_1b(bool tf);
void SetUse_TypeII_2b(bool tf);
void SetUse_TypeIII_2b(bool tf);
void SetUse_TypeGI_2b(bool tf);
void SetUse_TypeGII_2b(bool tf);
void SetUse_TypeGIIIa_2b(bool tf);
void SetUse_TypeGIIIb_2b(bool tf);
void SetUse_TypeGIIIc_2b(bool tf);
void SetUse_TypeGIVa_2b(bool tf);
void SetUse_TypeGIVb_2b(bool tf);
void SetUse_TypeGIVc_2b(bool tf);

void SetUse_GT_TypeI_2b(bool tf);
void SetUse_GT_TypeIV_2b(bool tf);

//    extern bool SlowVersion;
//    void UseSlowVersion(bool tf);
// factorize double commutator [Eta, [Eta, Gamma]]
void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z);

void comm223_231_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_231_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_232_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_132(const Operator &Eta, const Operator &Gamma, Operator &Z);
void comm223_132_cross(const Operator &Eta, const Operator &Gamma, Operator &Z);

} // namespace FactorizedDoubleCommutator
} // namespace Commutator

#endif
