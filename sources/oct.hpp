/// \file
// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    QFib: Fast and Efficient Brain Tractogram Compression
//    C. Mercier*, S. Rousseau*, P. Gori, I. Bloch and T. Boubekeur
//    NeuroInformatics 2020
//    DOI: 10.1007/s12021-020-09452-0
//
// All rights reserved. Use of this source code is governed by a 
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------
#pragma once

#include "bitvalue.hpp"
#include <iostream>
#include <climits>

///
/// \brief octDecode Function to decode a point encoding using the octahedral quantization
/// \param projected Point to decode
/// \param vec Point decoded
///

template<typename T>
void octDecode(const Bitvalue<T> & projected, float vec[3]) {

    vec[0] = projected.getFloat(0);
    vec[1] = projected.getFloat(1);
    vec[2] = 1.0f - (fabs(vec[0]) + fabs(vec[1]));

    if (vec[2] < 0.0f) {
        float oldX = vec[0];
        vec[0] = ((1.0f) - fabs(vec[1])) * signNotZero(oldX);
        vec[1] = ((1.0f) - fabs(oldX))   * signNotZero(vec[1]);
    }
}

///
/// \brief octEncode Function to encode a point using the octahedral quantization
/// \param vec Point to encode
/// \param projected Point encoded
///

template<typename T>
void octEncode(const float vec[3], Bitvalue<T> & projected) {
    const float invL1Norm = (1.0f) / (fabs(vec[0]) + fabs(vec[1]) + fabs(vec[2]));

    if (vec[2] < 0.0f) {
        projected.setfValue((1.0f - float(fabs(vec[1] * invL1Norm))) * signNotZero(vec[0]),
                            (1.0f - float(fabs(vec[0] * invL1Norm))) * signNotZero(vec[1]));
    } else {
        projected.setfValue(vec[0] * invL1Norm, vec[1] * invL1Norm);
    }
}

