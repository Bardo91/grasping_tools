/*
 * utilsGPIS.h
 *
 *  Created on: 5 Dec 2015
 *      Author: WolframMartens
 *
 *      static utilities for Gaussian process implicit surface (GPIS) computations
 */

#include "utilsGPIS.h"

using namespace arma;

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	uvec indsToDerInds(const uvec& inds, const int numDimsPlusOne) {
		uvec derInds(inds.size() * numDimsPlusOne);

		int indsIndex = -1;
		for (uvec::const_iterator indsIter = inds.begin(); indsIter != inds.end(); indsIter++) {
			indsIndex++;

			int dimension = -1;
			for (uvec::iterator derIndsIter = derInds.begin() + (indsIndex * numDimsPlusOne);
			derIndsIter != derInds.begin() + (indsIndex * numDimsPlusOne) + numDimsPlusOne;
				derIndsIter++) {
				dimension++;
				*derIndsIter = (*indsIter * numDimsPlusOne) + dimension;
			}
		}

		return derInds;
	}

	//-----------------------------------------------------------------------------------------------------------------
	uvec indsToDerInds(const int partIndex, const int numDimsPlusOne) {
		uvec derInds(numDimsPlusOne);
		int derIndsVal = partIndex * numDimsPlusOne;
		for (uvec::iterator iter = derInds.begin(); iter != derInds.end(); iter++) {
			*iter = derIndsVal;
			derIndsVal++;
		}

		return derInds;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void indsToDerInds(const int partIndex, const int numDimsPlusOne, std::vector<unsigned int> &derIndsOut) {
		int derIndsVal = partIndex * numDimsPlusOne;
		for (std::vector<unsigned int>::iterator iter = derIndsOut.begin(); iter != derIndsOut.end(); iter++) {
			*iter = derIndsVal;
			derIndsVal++;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	uvec indsToOnlyDerInds(const uvec& inds, const int numDims) {
		int numParts = inds.size();
		uvec onlyDerInds(numParts * (numDims));

		for (int n = 0; n < numParts; ++n) {
			for (int d = 0; d < numDims; ++d) {
				onlyDerInds(n * (numDims)+d) = inds(n) * (numDims + 1) + d + 1;
			}
		}

		return onlyDerInds;
	}

	//-----------------------------------------------------------------------------------------------------------------
	uvec indsToOnlyDerInds(const int partIndex, const int numDims){
		uvec onlyDerInds(numDims);
		for (int d = 0; d < numDims; ++d) {
			onlyDerInds(d) = partIndex * (numDims + 1) + d + 1;
		}
		return onlyDerInds;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void indsToOnlyDerInds(const int partIndex, const int numDims, std::vector<unsigned int> &onlyDerIndsOut) {
		for (int d = 0; d < numDims; ++d) {
			onlyDerIndsOut[d] = partIndex * (numDims + 1) + d + 1;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	// Compute covariance function between gradients and function values two points
	double covxixj(double _sigma, double _gamma, const vec &x1, const vec &x2, const int i, const int j) {
		if (i == 0 && j == 0)
			return covIsoSE(_sigma, _gamma, x1, x2);
		else if (i > 0 && j == 0)
			return -_gamma * (x1(i - 1) - x2(i - 1)) * covIsoSE(_sigma, _gamma, x1, x2);
		else if (i == 0 && j > 0)
			return _gamma * (x1(j - 1) - x2(j - 1)) * covIsoSE(_sigma, _gamma, x1, x2);
		else if (i == j)
			return _gamma * (1 - _gamma * (x1(i - 1) - x2(i - 1)) * (x1(j - 1) - x2(j - 1))) * covIsoSE(_sigma, _gamma, x1, x2);
		else
			return -_gamma*_gamma * (x1(i - 1) - x2(i - 1)) * (x1(j - 1) - x2(j - 1)) * covIsoSE(_sigma, _gamma, x1, x2);
	}

	//-----------------------------------------------------------------------------------------------------------------
	// Subfunction for isotropic squared exponential cov-function
	double covIsoSE(double _sigma, double _gamma, const vec &x1, const vec &x2){
		return _sigma * _sigma * exp(as_scalar(-2* _gamma * trans(x1 - x2) * (x1 - x2)));
	}

}	// namespace gpis