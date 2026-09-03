///////////////////////////////////////////////////////////////////////////////////
//    TensorCommutators.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
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

#ifndef TensorCommutators_hh
#define TensorCommutators_hh 1

#include "Operator.hh"
#include <map>
#include <array>
#include <deque>
#include <armadillo>
#include <iostream>

namespace Commutator
{

    void DoTensorPandyaTransformation(const Operator &Z, std::map<std::array<index_t, 2>, arma::mat> &);
    void DoTensorPandyaTransformation_SingleChannel(const Operator &Z, arma::mat &X, int ch_bra_cc, int ch_ket_cc);
    void AddInverseTensorPandyaTransformation(Operator &Z, const std::map<std::array<index_t, 2>, arma::mat> &);
    std::deque<arma::mat> InitializePandya(Operator &Z, size_t nch, std::string orientation, int X_parity);

    void comm111st(const Operator &X, const Operator &Y, Operator &Z);
    void comm121st(const Operator &X, const Operator &Y, Operator &Z);
    void comm122st(const Operator &X, const Operator &Y, Operator &Z);
    void comm221st(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hh_221st(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_phst(const Operator &X, const Operator &Y, Operator &Z);

    // scalar-tensor with a 3b operator
    // 3n7
    void comm331st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm223st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    // Production tts (T×T→scalar leftover). GEMM: 111, 122, 220, 221, 222.
    // 121/132/231/223/232 wrap to Reference (ss has no extra BLAS; 232 leftover
    // 6js are not a drop-in of comm232ss_srs_optimized). Gold: Mscheme_Test_comm*tts.
    void comm223tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm231st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm231tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm110tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm220tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm111tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm121tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm122tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm221tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hhtts(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_phtts(const Operator &X, const Operator &Y, Operator &Z);
    void comm232st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm232tts(const Operator &X, const Operator &Y, Operator &Z);
    void comm133st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm132st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm132tts(const Operator &X, const Operator &Y, Operator &Z);
    // above 3n7
    void comm332_pphhst(const Operator &X, const Operator &Y, Operator &Z);       // PASS the unit test
    void comm332_ppph_hhhpst(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm233_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);      // PASS the unit test
    void comm233_phst(const Operator &X, const Operator &Y, Operator &Z);         // PASS the unit test
    void comm333_ppp_hhhst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test
    void comm333_pph_hhpst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test
}// namespace Commutator

#endif
