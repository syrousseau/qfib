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


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "saveload.hpp"
#include <vector>

namespace utl
{
///
/// \brief meanMaxError Function to compute average and maximum error between two sets of fibers containing the exact same number of points and fibers
/// \param fibersA Vector containing the first set of fibers
/// \param fibersB Vector containing the second set of fibers
/// \param meanError Average error
/// \param maxError Maximum error
/// \param totalNbPts Variable containing the total number of points in the bundle, used for average error computation
///

void meanMaxError(const std::vector< std::vector<Vector3f> > & fibersA, const std::vector< std::vector<Vector3f> > & fibersB, double & meanError, float & maxError, uint64_t & totalNbPts)
{
	meanError = 0;
	maxError = 0.f;
	assert(fibersA.size() == fibersB.size());
	for(unsigned f = 0; f < fibersA.size(); ++f)
	{
		assert(fibersA[f].size() == fibersB[f].size());
		for(unsigned p = 0; p < fibersA[f].size(); ++p)
		{
			float err = (fibersA[f][p] - fibersB[f][p]).norm();
			maxError = err > maxError ? err : maxError;
			meanError  += err / double(totalNbPts);
		}
	}
}

///
/// \brief error Function to compute the error between two sets of fibers with the exact same number of points and fibers
/// This error is computed depending on the length of the fibers
/// \param fibersA Vector containing the first set of fibers
/// \param fibersB Vector containing the second set of fibers
/// \param minError Minimum error
/// \param meanError Average error
/// \param maxError Maximum error
/// \param minEndPointError Minimum error on every end point
/// \param meanEndPointError Average error on every end point
/// \param maxEndPointError Maximum error on every end point
/// \param totalNbPts Variable containing the total number of points in the bundle, used for average error computation
/// \param steps Number of length separation to perform
/// \param minValue Minimum length to consider
/// \param maxValue Maximum length to consider
///

void error(const std::vector<std::vector<Vector3f> > & fibersA, const std::vector< std::vector<Vector3f> > & fibersB,
		   vector<float>& minError, vector<double>& meanError, vector<float>& maxError,
		   float & minEndPointError, float & meanEndPointError, float & maxEndPointError, uint64_t & totalNbPts,
		   unsigned steps = 10, float minValue = 35, float maxValue = 260)
{
	minEndPointError = std::numeric_limits<float>::max();
	meanEndPointError = 0.f;
	maxEndPointError = 0.f;
	unsigned nbSteps = (maxValue - minValue)/steps;
	minError.resize(nbSteps, std::numeric_limits<float>::max());
	meanError.resize(nbSteps, 0);
	maxError.resize(nbSteps, 0);

	unsigned nbEndPoints = fibersA.size();
	assert(fibersA.size() == fibersB.size());
	for(unsigned f = 0; f < fibersA.size(); ++f)
	{
		unsigned tablePos = (float)((fibersA[f].size() - 1) * (fibersA[f][1]-fibersA[f][0]).norm() - minValue) / (float)(maxValue - minValue) * (nbSteps-1);
		if (tablePos>nbSteps-1) tablePos = nbSteps - 1;
		assert(fibersA[f].size() == fibersB[f].size());
		for(unsigned p = 0; p < fibersA[f].size()-1; ++p)
		{
			float err = (fibersA[f][p] - fibersB[f][p]).norm();
			meanError[tablePos]  += err / double(totalNbPts);
			maxError[tablePos] = err > maxError[tablePos] ? err : maxError[tablePos];
			minError[tablePos] = err < minError[tablePos] ? err : minError[tablePos];
		}
		float err = (fibersA[f][fibersA[f].size()-1] - fibersB[f][fibersA[f].size()-1]).norm();
		meanError[tablePos]  += err / double(totalNbPts);
		maxError[tablePos] = err > maxError[tablePos] ? err : maxError[tablePos];
		minError[tablePos] = err < minError[tablePos] ? err : minError[tablePos];
		minEndPointError = err < minEndPointError ? err : minEndPointError;
		meanEndPointError += err / double(nbEndPoints);
		maxEndPointError = err > maxEndPointError ? err : maxEndPointError;
	}
}
}
