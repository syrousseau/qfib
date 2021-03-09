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

#include "saveload.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

namespace vc
{
int my_round(const unsigned type, const float value)
{
	switch (type) {
	case 0:
		return floor(value);
		break;
	case 1:
		return ceil(value);
		break;
	default:
		return round(value);
		break;
	}
}

///
/// \brief computeTrilinearInterpolatedValue Compute an interpolated value at the position, considering the grid
/// \param grid Grid to consider for the interpolated value
/// \param dim Dimension of the grid
/// \param gridPosition Position at which to compute the interpolated value
/// \return Return the computed value
///

float computeTrilinearInterpolatedValue(std::vector<float> const & grid, std::vector<unsigned> const & dim, Vector3f const & gridPosition)
{
	std::vector<float> value;
	for (unsigned pos2=0; pos2<2; pos2++)
		for (unsigned pos3=0; pos3<2; pos3++)
			for (unsigned pos1=0; pos1<2; pos1++)
				value.push_back(grid[(int)my_round(pos1,gridPosition(0))+dim[0]*(int)my_round(pos2,gridPosition(1))+dim[0]*dim[1]*(int)my_round(pos3,gridPosition(2))]);
	float xd = gridPosition(0) - floor(gridPosition(0));
	float yd = gridPosition(1) - floor(gridPosition(1));
	float zd = gridPosition(2) - floor(gridPosition(2));

	float c00 = value[0]*(1-xd)+value[1]*xd;
	float c01 = value[2]*(1-xd)+value[3]*xd;
	float c10 = value[4]*(1-xd)+value[5]*xd;
	float c11 = value[6]*(1-xd)+value[7]*xd;

	float c0 = c00*(1-yd)+c10*yd;
	float c1 = c01*(1-yd)+c11*yd;

	float val = c0*(1-zd)+c1*zd;

	return val;
}

//positive filter
int getSign(const float value)
{
	if (value>0)
		return 1;
	else
		return 0;
}

float k(const float x, const float dirX, const float dimV)
{
	return (getSign(dirX) + round(x/dimV)) * dimV;
}

float k2(const float x, const float dirX)
{
	return floor(x) + getSign(dirX);
}

///
/// \brief getValueFromSegment Compute the value along a segment
/// \param grid Grid to consider
/// \param dim Dimensions of the grid
/// \param p1 First point of the segment
/// \param p2 Second point of the segment
/// \return Return the value of the segment
///

float getValueFromSegment(std::vector<float> const & grid, std::vector<unsigned> const & dim, Vector3f const & p1, Vector3f const & p2)
{
	Vector3f direction = p2-p1;
	Vector3f directionN = direction.normalized();
	float totalDist = direction.norm();
	//Number of voxels changes along each axis
	unsigned voxelsX = static_cast<unsigned>(fabs(floor(p2.x())-floor(p1.x())));
	unsigned voxelsY = static_cast<unsigned>(fabs(floor(p2.y())-floor(p1.y())));
	unsigned voxelsZ = static_cast<unsigned>(fabs(floor(p2.z())-floor(p1.z())));
	unsigned totalVoxelsCrossed = 1 + voxelsX + voxelsY + voxelsZ;
	Vector3f pTemp = p1;
	float value = 0;
	float epsilon = 0.0001f;
	float gammaX, gammaY, gammaZ, gamma;

	for (unsigned i=0; i<totalVoxelsCrossed - 1; i++)
	{
		if (fabs(directionN.x())<epsilon)
			gammaX =  numeric_limits<float>::max();
		else
			gammaX = (k2(pTemp.x(), direction.x()) - pTemp.x()) / direction.x();
		if (fabs(directionN.y())<epsilon)
			gammaY =  numeric_limits<float>::max();
		else
			gammaY = (k2(pTemp.y(), direction.y()) - pTemp.y()) / direction.y();
		if (fabs(directionN.z())<epsilon)
			gammaZ =  numeric_limits<float>::max();
		else
			gammaZ = (k2(pTemp.z(), direction.z()) - pTemp.z()) / direction.z();
		gamma = min(min(gammaX, gammaY), gammaZ);
		pTemp += (gamma) * direction;
		value += (grid[(int)floor(pTemp.x())+dim[0]*(int)floor(pTemp.y())+(dim[0]*dim[1]*(int)floor(pTemp.z()))]);
	}
	value += (grid[(int)floor(p2.x())+dim[0]*(int)floor(p2.y())+(dim[0]*dim[1]*(int)floor(p2.z()))]);
	return (value / (float)(totalVoxelsCrossed));
}

///
/// \brief compareFA Function to compare the FA values and display the results
/// \param fibers1 First set of fibers
/// \param fibers2 Second set of fibers
/// \param faFilename FA file to use for FA computation
///

void compareFA(std::vector<std::vector<Vector3f>> const & fibers1, std::vector<std::vector<Vector3f>> const & fibers2, const string & faFilename)
{
	cout << faFilename << endl;
	std::vector<float> faValues, vox;
	std::vector<unsigned> dim;
	MatrixXf transform;
	fls::loadFAValues(faFilename, faValues, dim, transform, vox);

	transform.conservativeResize(4,4);
	transform.row(3) << 0,0,0,1;
	transform = transform.inverse();

	unsigned maxThreads = std::max(unsigned(1), std::thread::hardware_concurrency());
	std::vector<float> ecartMoyenTab(maxThreads, 0);
	std::vector<uint64_t> nbPointsTab(maxThreads, 0);
	std::vector<float> ecartMaxTab(maxThreads, 0);
	std::vector<float> averageFiber1Error(maxThreads, 0);
	std::vector<float> averageFiber2Error(maxThreads, 0);
#ifdef USE_OPENMP
#pragma omp parallel for
	for (unsigned i=0; i<fibers1.size(); i++)
	{
		if (fibers1[i].size() == fibers2[i].size())
		{
			nbPointsTab[omp_get_thread_num()] += fibers1[i].size();
			for (unsigned j=0; j<fibers1[i].size(); j++)
			{
				Vector3f p = fibers1[i][j];
				for (unsigned k=0; k<3; k++)
					p(k)/=vox[k];
				Vector3f gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
				float value1 = computeTrilinearInterpolatedValue(faValues, dim, gridPosition);
				Vector3f p2 = fibers2[i][j];
				for (unsigned k=0; k<3; k++)
					p2(k)/=vox[k];
				Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
				float value2=computeTrilinearInterpolatedValue(faValues, dim, gridPosition2);
				float ecart = fabs(value1-value2);
				ecartMaxTab[omp_get_thread_num()] = max(ecart, ecartMaxTab[omp_get_thread_num()]);
				ecartMoyenTab[omp_get_thread_num()]+=ecart;
			}
		}

		//Per fiber error
		float fiber1Value = 0.f;
		float fiber2Value = 0.f;
		Vector3f p = fibers1[i][0];
		for (unsigned k=0; k<3; k++)
		{
			p(k)/=vox[k];
		}
		Vector3f gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
		fiber1Value += faValues[(int)floor(gridPosition.x())+dim[0]*(int)floor(gridPosition.y())+(dim[0]*dim[1]*(int)floor(gridPosition.z()))];
		for (unsigned j=0; j<fibers1[i].size()-1; j++)
		{
			Vector3f p1 = fibers1[i][j];
			Vector3f p2 = fibers1[i][j+1];
			for (unsigned k=0; k<3; k++)
			{
				p1(k)/=vox[k];
				p2(k)/=vox[k];
			}
			Vector3f gridPosition1 = (p1*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			fiber1Value += getValueFromSegment(faValues, dim, gridPosition1, gridPosition2);
		}
		p = fibers2[i][0];
		for (unsigned k=0; k<3; k++)
		{
			p(k)/=vox[k];
		}
		gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
		fiber2Value += faValues[(int)floor(gridPosition.x())+dim[0]*(int)floor(gridPosition.y())+(dim[0]*dim[1]*(int)floor(gridPosition.z()))];
		for (unsigned j=0; j<fibers2[i].size()-1; j++)
		{
			Vector3f p1 = fibers2[i][j];
			Vector3f p2 = fibers2[i][j+1];
			for (unsigned k=0; k<3; k++)
			{
				p1(k)/=vox[k];
				p2(k)/=vox[k];
			}
			Vector3f gridPosition1 = (p1*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			fiber2Value += getValueFromSegment(faValues, dim, gridPosition1, gridPosition2);
		}
		fiber1Value /= float(fibers1[i].size());
		fiber2Value /= float(fibers2[i].size());
		averageFiber1Error[omp_get_thread_num()] += fiber1Value;
		averageFiber2Error[omp_get_thread_num()] += fiber2Value;
	}
#else
	parallel.For<unsigned>(0, fibers1.size(), 1, [&](unsigned i)
	{
		if (fibers1[i].size() == fibers2[i].size())
		{
			nbPointsTab[parallel.getID()] += fibers1[i].size();
			for (unsigned j=0; j<fibers1[i].size(); j++)
			{
				Vector3f p = fibers1[i][j];
				for (unsigned k=0; k<3; k++)
					p(k)/=vox[k];
				Vector3f gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
				float value1 = computeTrilinearInterpolatedValue(faValues, dim, gridPosition);
				Vector3f p2 = fibers2[i][j];
				for (unsigned k=0; k<3; k++)
					p2(k)/=vox[k];
				Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
				float value2=computeTrilinearInterpolatedValue(faValues, dim, gridPosition2);
				float ecart = fabs(value1-value2);
				ecartMaxTab[parallel.getID()] = max(ecart, ecartMaxTab[parallel.getID()]);
				ecartMoyenTab[parallel.getID()]+=ecart;
			}
		}

		//Per fiber error
		float fiber1Value = 0.f;
		float fiber2Value = 0.f;
		Vector3f p = fibers1[i][0];
		for (unsigned k=0; k<3; k++)
		{
			p(k)/=vox[k];
		}
		Vector3f gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
		fiber1Value += faValues[(int)floor(gridPosition.x())+dim[0]*(int)floor(gridPosition.y())+(dim[0]*dim[1]*(int)floor(gridPosition.z()))];
		for (unsigned j=0; j<fibers1[i].size()-1; j++)
		{
			Vector3f p1 = fibers1[i][j];
			Vector3f p2 = fibers1[i][j+1];
			for (unsigned k=0; k<3; k++)
			{
				p1(k)/=vox[k];
				p2(k)/=vox[k];
			}
			Vector3f gridPosition1 = (p1*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			fiber1Value += getValueFromSegment(faValues, dim, gridPosition1, gridPosition2);
		}
		p = fibers2[i][0];
		for (unsigned k=0; k<3; k++)
		{
			p(k)/=vox[k];
		}
		gridPosition = (p*transform.block(0,0,3,3)+transform.block(0,3,3,1));
		fiber2Value += faValues[(int)floor(gridPosition.x())+dim[0]*(int)floor(gridPosition.y())+(dim[0]*dim[1]*(int)floor(gridPosition.z()))];
		for (unsigned j=0; j<fibers2[i].size()-1; j++)
		{
			Vector3f p1 = fibers2[i][j];
			Vector3f p2 = fibers2[i][j+1];
			for (unsigned k=0; k<3; k++)
			{
				p1(k)/=vox[k];
				p2(k)/=vox[k];
			}
			Vector3f gridPosition1 = (p1*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			Vector3f gridPosition2 = (p2*transform.block(0,0,3,3)+transform.block(0,3,3,1));
			fiber2Value += getValueFromSegment(faValues, dim, gridPosition1, gridPosition2);
		}
		fiber1Value /= float(fibers1[i].size());
		fiber2Value /= float(fibers2[i].size());
		averageFiber1Error[parallel.getID()] += fiber1Value;
		averageFiber2Error[parallel.getID()] += fiber2Value;
	});
#endif

	float ecartMoyen = ecartMoyenTab[0], ecartMax = ecartMaxTab[0];
	float averageFiber1 = averageFiber1Error[0], averageFiber2 = averageFiber2Error[0];
	uint64_t nbPoints = nbPointsTab[0];
	for (unsigned i=1; i<maxThreads; i++)
	{
		nbPoints += nbPointsTab[i];
		ecartMoyen += ecartMoyenTab[i];
		ecartMax = max(ecartMax, ecartMaxTab[i]);
		averageFiber1 += averageFiber1Error[i];
		averageFiber2 += averageFiber2Error[i];
	}

	ecartMoyen /= nbPoints;
	averageFiber1 /= (float)maxThreads;
	averageFiber2 /= (float)maxThreads;
	cout << "fa per fiber mean error: " << fabs(averageFiber2 - averageFiber1) / averageFiber1 * 100 << "%" << endl;
	cout << "fa max error: " << ecartMax << endl;
	cout << "fa mean error: " << ecartMoyen << endl;
}
}
