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

#include <cstdint>

constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;

inline float signNotZero(const float value) 
{
    return (value < 0.0f) ? -1.0f : 1.0f;
}

inline float clamp(const float value, const float low, const float high) 
{
    if (value <= low) 
    {
        return low;
    } else if (value >= high) 
    {
        return high;
    } else 
    {
        return value;
    }
}

//Number of bits inside the structure (that contains two elements of bitcount/2 bits
static uint8_t bitcount;

///
/// \brief The Bitvalue class Class to store two values of bitcount/2 bits in a bitcount bits element
/// This is especially useful to handle 2 int4 values that represent the coordinates of a point with our compression in 8 bits
///

template<typename T>
class Bitvalue 
{
private:
    T m_bits;
    explicit Bitvalue<T>(T a, T b) : m_bits((a << bitcount/2) | (b  & ((uint64_t(1) << (bitcount/2)) - 1))) {}
public:
    explicit Bitvalue() : m_bits(0) {}

    static Bitvalue fromBits(T a, T b) 
    {
        return Bitvalue<T>(a, b);
    }

    T operator[](const uint8_t elmt) const 
    {
        T returnValue;
        switch (elmt) {
        case 1:
            returnValue = T(m_bits << (bitcount/2)) >> (bitcount/2);
            break;
        case 0:
            returnValue = m_bits >> bitcount/2;
            break;
        default:
            std::cerr << "Access to a non existent element" << std::endl;
            exit(EXIT_FAILURE);
            break;
        }
        return returnValue;
    }

    void setfValue(const float f1, const float f2) 
    {
        T a =  (T)round(clamp(f1, -1.0f, 1.0f) * ((uint64_t(1) << (bitcount/2 - 1)) - 1));
        T b =  (T)round(clamp(f2, -1.0f, 1.0f) * ((uint64_t(1) << (bitcount/2 - 1)) - 1));
        m_bits = ((a << bitcount/2) | (b  & ((uint64_t(1) << (bitcount/2)) - 1)));
    }

    float getFloat(const uint8_t elmt) const 
    {
        return float(clamp(int((*this)[elmt]) * (1.0f / float((uint64_t(1) << (bitcount/2 - 1)) - 1)), -1.0f, 1.0f));
    }

    T bits(const uint8_t elmt) const 
    {
        return (*this)[elmt];
    }

    T getCombinedValue() const 
    {
        return m_bits;
    }

    void setCombinedValue(const T value) 
    {
        m_bits = value;
    }

};

