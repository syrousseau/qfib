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

// This file is part of the reference implementation for the paper
//    Fast Lossy Compression of 3D Unit Vector Sets
//    Sylvain Rousseau and Tamy Boubekeur
//    SIGGRAPH Asia 2017 Technical Briefs (SA '17)
//    DOI: https://doi.org/10.1145/3145749.3149436
//
// All rights reserved. Use of this source code is governed by a 
// MIT license


/// implementation of spherical fibonacci and inverse spherical fibonacci inspired by 
/// https://www.shadertoy.com/view/lllXz4 witch using the method described in the paper: 
/// Spherical Fibonacci Mapping
/// Benjamin Keinert, Matthias Innmann, Michael Snger and Marc Stamminger
/// Proc. of SIGGRAPH Asia 2015
/// The original source code was realized by Inigo Quilez and was released under the MIT license

#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include <cmath>

#include "bitvalue.hpp"

using namespace Eigen;

namespace fib
{
	#undef PI
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
	const double TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846;
	const double PHI = 1.6180339887498948482045868343656381177203091798057628621354486227;

    double myfract(const double x)
    {
        return x - std::floor(x);
    }

    Vector2d myfract(const Vector2d x)
    {
        return x - Vector2d(std::floor(x.x()), std::floor(x.y()));
    }

    Vector2d myfloor(const Vector2d v)
    {
        return Vector2d(floor(v.x()), floor(v.y()));
    }

    ///
    /// \brief SF Get the coordinates of point $id of the spherical fibonacci point set containing $nbPoints points
    /// \param id Point in spherical Fibonacci coordinates of which we want to obtain the Cartesian coordinates
    /// \param nbPoints Number of point in the Fibonacci point set
    /// \return Cartesian coordinates of the point
    ///

    template<typename T>
    Vector3f SF(const Bitvalue<T> id, const uint32_t nbPoints)
	{
        double m = 1.0 - 1.0 / nbPoints;
        double phi = TWO_PI * myfract(id.getCombinedValue() * PHI);
        double cosTheta = m - 2.0 * double(id.getCombinedValue()) / nbPoints;
		double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        return Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta).normalized();
	}




    // Note that the inverseSF method have numerical precision issue when using more than 2^22 points. If you try to develop a 32-bits version 
    // of our compression code, you may have some issue with the spherica fibonacci quantization. The octahedral quantization is more advised in that case.

    ///
    /// \brief inverseSF Get the id of the closest point to $p$ in the spherical fibonacci point set containing $n$ points
    /// due to numerical instabilities, keep n under 2^24
    /// \param vec Point in Cartesian coordinates
    /// \param n Number of points in the Fibonacci point set
    /// \return Spherical Fibonacci coordinates of the point
    ///

    template<typename T>
    T inverseSF(const Vector3f & vec, const unsigned n)
    {
        Vector3d p(vec.x(), vec.y(), vec.z());
        double m = 1.0 - 1.0 / n;
        double phi = std::min(std::atan2(p.y(), p.x()), PI), cosTheta = p.z();
        double k = std::max(2.0, floor(log(n * PI * sqrt(5.0) * (1.0 - cosTheta*cosTheta)) / log(PHI * PHI)));
        double Fk = pow(PHI, k) / sqrt(5.0);
        Vector2d  F = Vector2d(round(Fk), round(Fk * PHI));

        Vector2d ka = 2.0 * F / double(n);
        Vector2d kb = TWO_PI *(myfract((F + Vector2d(1.0, 1.0)) * PHI) - Vector2d(PHI - 1.0, PHI - 1.0));
        Matrix2d iB;
        iB << ka.y(), kb.y(), -ka.x(), -kb.x();
        iB /= (ka.y() * kb.x() - ka.x() * kb.y());

        Vector2d c = myfloor(iB * Vector2d(phi, cosTheta - m));
        double d = 8.0;
        T j = 0;
        for (unsigned s = 0; s < 4; ++s)
        {
            Vector2d uv = Vector2d(double(s - 2 * (s / 2)), double(s / 2));
            double i = F.dot(uv + c);

            double phi = TWO_PI * myfract(i * PHI);
            double cosTheta = m - 2.0 * i / n;
            double sinTheta = sqrt(1.0 - cosTheta*cosTheta);

            Vector3d q = Vector3d(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
            double squaredDistance = (q - p).dot(q - p);
            if (squaredDistance < d)
            {
                d = squaredDistance;
                j = T(i);
            }
        }
        return j;
    }


}
