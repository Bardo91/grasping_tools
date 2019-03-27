/*
 * utilsGPIS.h
 *
 *  Created on: 5 Dec 2015
 *      Author: WolframMartens
 *
 *      static utilities for Gaussian process implicit surface (GPIS) computations
 */


#ifndef UTILSGPIS_H
#define UTILSGPIS_H

#include <iostream>

#define ARMA_NO_DEBUG
#include <armadillo>

namespace gpis {
	arma::uvec indsToDerInds(const arma::uvec& inds, const int numDimsPlusOne);

	arma::uvec indsToDerInds(const int partIndex, const int numDimsPlusOne);

	void indsToDerInds(const int partIndex, const int numDimsPlusOne, std::vector<unsigned int> &derIndsOut);

	arma::uvec indsToOnlyDerInds(const arma::uvec& inds, const int numDims);

	arma::uvec indsToOnlyDerInds(const int partIndex, const int numDims);

	void indsToOnlyDerInds(const int partIndex, const int numDims, std::vector<unsigned int> &onlyDerIndsOut);

	double covxixj(double _sigma, double _gamma, const arma::vec &x1, const arma::vec &x2, const int d1, const int d2);

	double covIsoSE(double _sigma, double _gamma, const arma::vec &x1, const arma::vec &x2);
}	// namespace gpis

#endif // UTILSGPIS_H
