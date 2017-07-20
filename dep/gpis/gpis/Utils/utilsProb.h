/*
 * utilsProb.h
 *
 *  Created on: 5 Dec 2015
 *      Author: WolframMartens
 *
 *      static utilities for probability computations
 */

#ifndef INCLUDE_DIRICHLETGPIS_UTILSPROB_H_
#define INCLUDE_DIRICHLETGPIS_UTILSPROB_H_

#define ARMA_NO_DEBUG
#include <armadillo>

namespace gpis {
	double computeGaussianPdf(const arma::vec &x, const arma::vec &mean, const arma::mat &covMat);

	double computeGaussianLogPdf(const arma::vec &x, const arma::vec &mean, const arma::mat &covMat);

	arma::vec productOf1DGaussians(const arma::vec &muVec, const double prec);

	int sampleUniVarCat(const arma::colvec &prob);

	void multiVarGaussSample(arma::vec &outVec);

	int CRPSample(const arma::vec &probVec, const arma::uvec &numAssocPartsVec, const double &alpha);

	void reinitRandomNumbers();

	double evalGammaDist(double _x, double _alpha, double _beta);

	double sampleGammaDist(double _alpha, double _beta);
}	// namespace gpis

#endif /* INCLUDE_DIRICHLETGPIS_UTILSPROB_H_ */
