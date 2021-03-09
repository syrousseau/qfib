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
#include <bitset>

#include "oct.hpp"
#include "sfibonacci.h"
#include "bitvalue.hpp"
#include "parallel.hpp"

using namespace Eigen;

//struct CompressedVector {Bitvalue<T>pv;};
template<typename T>
using CompressedVector = Bitvalue<T>;

template<typename T>
struct CompressedFiber
{
	Vector3f origin;
	Vector3f second;
	std::vector<CompressedVector<T>> data;
};

float epsMapping = 0.08f;

enum QuantizationMethod
{
	OCTAHEDRAL,
	SPHERICALFIBONACCI
};

QuantizationMethod quantizationMethod = QuantizationMethod::SPHERICALFIBONACCI;

namespace fc
{
template<typename T>
Vector3f unquantize(const CompressedVector<T> & cvec)
{
	switch(quantizationMethod)
	{
	float vec[3];
	case(QuantizationMethod::OCTAHEDRAL):
		octDecode(cvec, vec);
		return Vector3f(vec[0], vec[1], vec[2]).normalized();
	case(QuantizationMethod::SPHERICALFIBONACCI):
		static uint32_t nbFiboPoints = 1 << (bitcount);
		return fib::SF(cvec, nbFiboPoints);
	default :
		std::cerr << "ERROR: Quantization Method not handled !" << std::endl;
		exit(EXIT_FAILURE);
		break;

	}
}

template<typename T>
void quantize(const Vector3f & vec, CompressedVector<T> & cvec)
{
	switch(quantizationMethod)
	{
	case(QuantizationMethod::OCTAHEDRAL) :
		octEncode(vec.data(), cvec);
		break;
	case(QuantizationMethod::SPHERICALFIBONACCI):
		static uint32_t nbFiboPoints = 1 << (bitcount);
		cvec.setCombinedValue(fib::inverseSF<T>(vec, nbFiboPoints));
		break;
	default:
		std::cerr << "ERROR: Quantization Method not handled !" << std::endl;
		exit(EXIT_FAILURE);
		break;
	}
}

///
/// \brief uniformMapping Uniform mapping from uniquant (Fast lossy compression of 3d unit vectors)
/// \param v Vector to map
/// \param axis Axis to use for the spherical cap
/// \param ratio Ratio to use for the mapping
/// \return The mapped vector
///
Vector3f uniformMapping(const Vector3f & v, const Vector3f & axis, const float ratio)
{
	Vector3d dv = v.cast<double>();
	Vector3d daxis = axis.cast<double>();
	dv.normalize();
	daxis.normalize();

	double K = 1.0 / ratio;
	double c = dv.dot(daxis);

	// numerical instabilities
	if (c > 1.0)
		c = 1.0;
	if (c < -1.0)
		c = -1.0;
	Vector3d p1 = (dv - c * daxis).normalized();
	double delta = (1.0 - ((1.0 - c) / K));
	return  ((p1 * std::sqrt(1.0 - delta * delta) + delta * daxis).cast<float>()).normalized();
}

///
/// \brief inverseUniformMapping Inverse uniform mapping from uniquant (Fast lossy compression of 3d unit vectors)
/// \param v Vector to unmap
/// \param axis Axis to use for the spherical cap
/// \param ratio Ratio to use for the mapping
/// \return The unmapped vector
///

Vector3f inverseUniformMapping(const Vector3f & v, const Vector3f & axis, const float ratio)
{
	//return v;
	Vector3d v2 = v.cast<double>();
	Vector3d axis2 = axis.cast<double>();
	v2.normalize();
	axis2.normalize();

	double K = ratio;
	double c = v2.dot(axis2);

	/// numerical instabilities
	if (c > 1.0)
		c = 1.0;
	if (c < -1.0)
		c = -1.0;
	Vector3d p1 = (v2 - c * axis2).normalized();
	double delta = (1.0 - ((1.0 - c) / K));
	return  ((p1 * std::sqrt(1.0 - delta * delta) + delta * axis2).cast<float>()).normalized();
}

///
/// \brief computeRatio Compute the minimum dot product of the input fiber (only one)
/// \param fiber Fiber from which to compute the minimum dot product
/// \return Minimum dot product of the unique fiber
///

float computeRatio(const std::vector<Vector3f> & fibers)
{
	Vector3f v1, v2;
	float mindot = 2.f;
	v1 = (fibers[1] - fibers[0]).normalized();
	for(unsigned p = 2; p < fibers.size(); ++p)
	{
		v2 = (fibers[p] - fibers[p-1]).normalized();
		float dot = v2.dot(v1);
		mindot = (dot <  mindot) ? dot : mindot;
		v1 = v2;
	}

	return mindot;
}

///
/// \brief computeRatioVerif Compute the ratio of the input fibers and verify that stepsizes are almost constant
/// \param fibers Fibers from which to compute the ratio
/// \param verbose Parameter to display the maximum angle if true
/// \param totalNbPts Variable that will contain the total number of points in the bundle, used for average error computation
/// \return Ratio of the bundle of fibers
///

float computeRatioVerif(const std::vector<std::vector<Vector3f> > & fibers, bool verbose, uint64_t & totalNbPts)
{
	unsigned maxThreads = std::max(unsigned(1), std::thread::hardware_concurrency());
	std::vector<uint64_t> nbPts(maxThreads, 0);
	std::vector<float> mindot(maxThreads, 2.0f);
	std::vector<float> maxStep(maxThreads, 0.0f);
	std::vector<float> minStep(maxThreads, std::numeric_limits<float>::max());
#ifdef USE_OPENMP
#pragma omp parallel for
	for(int f = 0; f < fibers.size(); ++f)
	{
		nbPts[omp_get_thread_num()] += fibers[f].size();
		Vector3f v1, v2;
		float v2norm;
		v1 = (fibers[f][1] - fibers[f][0]).normalized();
		for(unsigned p = 2; p < fibers[f].size(); ++p)
		{
			v2 = (fibers[f][p] - fibers[f][p-1]);
			v2norm = v2.norm();
			maxStep[omp_get_thread_num()] = (v2norm > maxStep[omp_get_thread_num()]) ? v2norm : maxStep[omp_get_thread_num()];
			minStep[omp_get_thread_num()] = (v2norm < minStep[omp_get_thread_num()]) ? v2norm : minStep[omp_get_thread_num()];
			float dot = (v2 / v2norm).dot(v1);
			mindot[omp_get_thread_num()] = (dot <  mindot[omp_get_thread_num()]) ? dot : mindot[omp_get_thread_num()];
			v1 = v2 / v2norm;
		}
	}
#else
	parallel.For<unsigned>(0, fibers.size(), 1, [&](unsigned f)
	{
		nbPts[parallel.getID()] += fibers[f].size();
		Vector3f v1, v2;
		float v2norm;
		v1 = (fibers[f][1] - fibers[f][0]).normalized();
		for(unsigned p = 2; p < fibers[f].size(); ++p)
		{
			v2 = (fibers[f][p] - fibers[f][p-1]);
			v2norm = v2.norm();
			maxStep[parallel.getID()] = (v2norm > maxStep[parallel.getID()]) ? v2norm : maxStep[parallel.getID()];
			minStep[parallel.getID()] = (v2norm < minStep[parallel.getID()]) ? v2norm : minStep[parallel.getID()];
			float dot = (v2 / v2norm).dot(v1);
			mindot[parallel.getID()] = (dot <  mindot[parallel.getID()]) ? dot : mindot[parallel.getID()];
			v1 = v2 / v2norm;
		}
	});
#endif
	float mdot = mindot[0], maxS = maxStep[0], minS = minStep[0];
	totalNbPts=nbPts[0];
	for (unsigned int i=1; i<mindot.size(); i++)
	{
		totalNbPts += nbPts[i];
		mdot = (mindot[i] < mdot) ? mindot[i] : mdot;
		maxS = (maxStep[i] > maxS) ? maxStep[i] : maxS;
		minS = (minStep[i] < minS) ? minStep[i] : minS;
	}
	float diff = maxS - minS;

	if (verbose)
	{
		std::cout << "Max angle: " << acos(mdot) * 180 / fib::PI << "°" << std::endl;
		std::cout << "Max length between segments: " << maxS << std::endl;
		std::cout << "Min length between segments: " << minS << std::endl;
		std::cout << "Difference: " << diff << std::endl;
	}

	if (diff > 0.1f * (maxS + minS) / 2.0f)
	{
		std::cerr << std::endl << "-------------!ERROR!--------------" << std::endl;
		std::cerr << "Difference in segments length too big (" << diff << "), no compression possible" << std::endl;
		std::cerr << "Try to resample your data with a constant stepsize before compression" << std::endl;
		std::cerr << "-------------!ERROR!--------------" << std::endl << std::endl;
		exit(EXIT_FAILURE);
	}

	epsMapping = 0.f;
	for(unsigned i = 0; i < 20; ++i)
		epsMapping = 3 * sqrt(2) * pow(2, -bitcount/2.f) * ((1 - (mdot))/2+epsMapping);
	//Error maximization proportionaly to itself (epsMapping < 1)
	epsMapping = sqrt(epsMapping);
	return std::min(1.f, (1.f - mdot) / 2.f + epsMapping);
}

///
/// \brief computeRatioVerif Compute the ratio of the input fibers without verifying that stepsizes are almost constant
/// \param fibers Fibers from which to compute the ratio
/// \param verbose Parameter to display the maximum angle if true
/// \param totalNbPts Variable that will contain the total number of points in the bundle, used for average error computation
/// \return Ratio of the bundle of fibers
///

float computeRatio(const std::vector<std::vector<Vector3f> > & fibers, bool verbose, uint64_t & totalNbPts)
{
	unsigned maxThreads = std::max(unsigned(1), std::thread::hardware_concurrency());
	std::vector<uint64_t> nbPts(maxThreads, 0);
	std::vector<float> mindot(maxThreads, 2.0f);

#ifdef USE_OPENMP
#pragma omp parallel for
	for(int f = 0; f < fibers.size(); ++f)
	{
		nbPts[omp_get_thread_num()] += fibers[f].size();
		Vector3f v1, v2;
		float v2norm;
		v1 = (fibers[f][1] - fibers[f][0]).normalized();
		for(unsigned p = 2; p < fibers[f].size(); ++p)
		{
			v2 = (fibers[f][p] - fibers[f][p-1]).normalized();
			float dot = v2.dot(v1);
			mindot[omp_get_thread_num()] = (dot <  mindot[omp_get_thread_num()]) ? dot : mindot[omp_get_thread_num()];
			v1 = v2;
		}
	}
#else
	parallel.For<unsigned>(0, fibers.size(), 1, [&](unsigned f)
	{
		nbPts[parallel.getID()] += fibers[f].size();
		Vector3f v1, v2;
		float v2norm;
		v1 = (fibers[f][1] - fibers[f][0]).normalized();
		for(unsigned p = 2; p < fibers[f].size(); ++p)
		{
			v2 = (fibers[f][p] - fibers[f][p-1]).normalized();
			float dot = v2.dot(v1);
			mindot[parallel.getID()] = (dot <  mindot[parallel.getID()]) ? dot : mindot[parallel.getID()];
			v1 = v2;
		}
	});
#endif
	float mdot=mindot[0];
	totalNbPts=nbPts[0];
	for (unsigned int i=1; i<mindot.size(); i++)
	{
		totalNbPts += nbPts[i];
		mdot = (mindot[i] < mdot) ? mindot[i] : mdot;
	}

	if (verbose)
		std::cout << "Max angle: " << acos(mdot) * 180 / fib::PI << "°" << std::endl;

	epsMapping = 0.f;
	for(unsigned i = 0; i < 20; ++i)
		epsMapping = 3 * sqrt(2) * pow(2, -bitcount/2.f) * ((1 - (mdot))/2+epsMapping);
	//Error maximization proportionaly to itself (epsMapping < 1)
	epsMapping = sqrt(epsMapping);
	return std::min(1.f, (1.f - mdot) / 2.f + epsMapping);
}

///
/// \brief compressFiber Compression of the fiber considering its ratio
/// \param fiber Fiber to compress
/// \param ratio Ratio to use for the compression
/// \param cfiber Fiber compressed
/// \param template T is the type of the compressed data, either int8 or int16 depending on the precision asked
///

template<typename T>
void compressFiber(const std::vector<Vector3f> & fiber, const float ratio, CompressedFiber<T> & cfiber)
{
	assert(fiber.size() > 1);

	//encode the first two points
	cfiber.origin = fiber[0];
	cfiber.second = fiber[1];

	cfiber.data.resize(fiber.size() - 2);

	//compute the first unit vector and the stepsize
	Vector3f axis = (fiber[1] - fiber[0]).normalized();
	float stepsize = (fiber[1] - fiber[0]).norm();
	Vector3f currentpoint = fiber[1];

	//compress the rest of the fiber
	Vector3f v;
	constexpr float eps = 0.00000001;
	for(unsigned i = 2; i < fiber.size(); ++i)
	{
		v = (fiber[i] - currentpoint).normalized();
		if(((1.f - v.dot(axis)) / 2.f) > ratio) //Conservative mapping in case the epsilon is not enough
		{
			fc::quantize((-((1.0 - eps) * axis) + (eps * v)).normalized(), cfiber.data[i-2]);
		}
		else //Otherwise, general case
		{
			if(abs(v.dot(axis)) < 1.f)
				fc::quantize(fc::inverseUniformMapping(v, axis, ratio), cfiber.data[i-2]);
			else
			{
				fc::quantize(v, cfiber.data[i-2]);
			}
		}
		axis = fc::uniformMapping(fc::unquantize(cfiber.data[i-2]), axis, ratio);
		currentpoint += (stepsize * axis);
		//currentpoint = fiber[i]; //Uncomment if you want to disable error propagation !!
	}
}

///
/// \brief compressFiber Compression of the fiber considering its ratio, computing the obtained error at the same time
/// \param fiber Fiber to compress
/// \param ratio Ratio to use for the compression
/// \param cfiber Fiber compressed
/// \param maxError Maximum error (at point level) when compressing the fiber
/// \param meanError Average error (at point level) when compressing the fiber
/// \param template T is the type of the compressed data, either int8 or int16 depending on the precision asked
///

template<typename T>
void compressFiber(const std::vector<Vector3f> & fiber, const float ratio, CompressedFiber<T> & cfiber, float & maxError, double & meanError)
{
	assert(fiber.size() > 1);

	//encode the first two points
	cfiber.origin = fiber[0];
	cfiber.second = fiber[1];

	cfiber.data.resize(fiber.size() - 2);

	//compute the first unit vector and the stepsize
	Vector3f axis = (fiber[1] - fiber[0]).normalized();
	float stepsize = (fiber[1] - fiber[0]).norm();
	Vector3f currentpoint = fiber[1];

	//compress the rest of the fiber
	Vector3f v;
	constexpr float eps = 0.00000001;
	for(unsigned i = 2; i < fiber.size(); ++i)
	{
		v = (fiber[i] - currentpoint).normalized();
		if(((1.f - v.dot(axis)) / 2.f) > ratio) //Conservative mapping in case the epsilon is not enough
		{
			fc::quantize((-((1.0 - eps) * axis) + (eps * v)).normalized(), cfiber.data[i-2]);
		}
		else //Otherwise, general case
		{
			if(abs(v.dot(axis)) < 1.f)
				fc::quantize(fc::inverseUniformMapping(v, axis, ratio), cfiber.data[i-2]);
			else
				fc::quantize(v, cfiber.data[i-2]);
		}
		axis = fc::uniformMapping(fc::unquantize(cfiber.data[i-2]), axis, ratio);
		currentpoint += (stepsize * axis);
		float err = (fiber[i] - currentpoint).norm();
		maxError = err > maxError ? err : maxError;
		meanError += err;
	}
}

///
/// \brief decompressFiber Decompress the input fiber
/// \param cfiber Compressed fiber to decompress
/// \param ratio Ratio to use for decompression (same than for compression)
/// \param fiber Decompressed fiber
///

template<typename T>
void decompressFiber(const CompressedFiber<T> & cfiber, const float ratio, std::vector<Vector3f> & fiber)
{
	fiber.reserve(cfiber.data.size() + 2);
	fiber.push_back(cfiber.origin);
	fiber.push_back(cfiber.second);

	float stepsize = (cfiber.second - cfiber.origin).norm();
	Vector3f axis = (fiber[1] - fiber[0]).normalized();

	for(unsigned i = 0; i < cfiber.data.size(); ++i)
	{
		Vector3f v = fc::unquantize(cfiber.data[i]);
		if(abs(v.dot(axis)) < 1.f)
			v = fc::uniformMapping(v, axis, ratio);

		fiber.push_back(fiber[i+1] + stepsize * v);
		axis = v;
	}
}
}
