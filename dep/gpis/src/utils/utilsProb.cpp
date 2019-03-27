/*
 * utilsProb.cpp
 *
 *  Created on: 5 Dec 2015
 *      Author: WolframMartens
 *
 *      static utilities for probability computations
 */

#include <gpis/utils/utilsProb.h>
#ifdef _WIN32
	#define _USE_MATH_DEFINES
#endif

#include <math.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <mutex>

using namespace arma;

namespace gpis {
	const double two_M_PI = 2.0 * M_PI;

	boost::mt19937 engineU;
	boost::uniform_real<> uDistr(0.0, 1.0);
	boost::variate_generator<boost::mt19937, boost::uniform_real<>> genU(engineU, uDistr);
	std::mutex uniformMutex;
	boost::mt19937 engineN;
	boost::normal_distribution<double> nDistr(0.0, 1.0);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<double>> genN(engineN, nDistr);
	std::mutex gaussianMutex;

	std::default_random_engine sGenerator;

	//unsigned int SEED_DEFAULT = 5489u;
	unsigned int SEED_DEFAULT = 5488u;

	//-----------------------------------------------------------------------------------------------------------------
	void reinitRandomNumbers() {
		genU.engine().seed(SEED_DEFAULT);
		genU.distribution().reset();
		genN.engine().seed(SEED_DEFAULT);
		genN.distribution().reset();
	}


	//-----------------------------------------------------------------------------------------------------------------
	double computeGaussianPdf(const arma::vec &x, const arma::vec &mean, const arma::mat &covMat){
		//	numDims = x.n_rows;
		double factor = 1 / sqrt(pow(two_M_PI, x.n_rows) * det(covMat));
		arma::vec diff;
		diff = x - mean;
		return factor * exp(as_scalar(-0.5 * trans(diff) * inv(covMat) * diff));
	}

	//-----------------------------------------------------------------------------------------------------------------
	double computeGaussianLogPdf(const arma::vec &x, const arma::vec &mean, const arma::mat &covMat) {
		arma::vec diff;
		diff = x - mean;
		double logPdf = - 0.5 * (as_scalar(trans(diff) * inv(covMat) * diff)
				+ log(det(covMat))
				+ x.n_rows * log(2 * M_PI));
		return logPdf;
	}

	//-----------------------------------------------------------------------------------------------------------------
	vec productOf1DGaussians(const arma::vec &muVec, const double prec) {
		// Compute the mean and precision resulting from the product of N
		// gaussians with identical precision and different means
		int N = muVec.n_elem;
		vec muAndPrec(2);
		muAndPrec(1) = N * prec;
		muAndPrec(0) = 1/muAndPrec(1) * sum(prec * muVec);
		return muAndPrec;
	}

	//-----------------------------------------------------------------------------------------------------------------
	int sampleUniVarCat(const arma::colvec &prob) {
		//// Remove randomness
		//unsigned index;
		//prob.max(index);
		//return index;
		// Regular sampling
		arma::colvec probNormed;
		probNormed = prob / sum(prob);
		//    cout << probNormed.t() << endl;
		uniformMutex.lock();
		double uVar = genU();
		uniformMutex.unlock();
		double sum = probNormed.at(0);
		int i = 0;
		while (sum < uVar) {
			++i;
			sum += probNormed.at(i);
		}
		return i;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void multiVarGaussSample(arma::vec &outVec) {
		//// Remove randomness
		//outVec.zeros();return;
		// Regular sampling // Why not randn()?
		outVec.imbue([&]() { return genN();});
		//outVec.randn(outVec.n_elem);
	}

	//-----------------------------------------------------------------------------------------------------------------
	int CRPSample(const vec &probVec, const uvec &numAssocPartsVec, const double &alpha) {
		return sampleUniVarCat(join_cols(numAssocPartsVec % probVec, alpha * vec(1).ones()));
	}

	//-----------------------------------------------------------------------------------------------------------------
	double evalGammaDist(double _x, double _alpha, double _beta) {
		if (_x <= 0 || _alpha <= 0 || _beta <= 0)
			return 0.0;
		else
			return pow(_beta, _alpha) * pow(_x, _alpha - 1)*exp(-_beta*_x) / tgamma(_alpha);
	}

	//-----------------------------------------------------------------------------------------------------------------
	double sampleGammaDist(double _alpha, double _beta) {
		std::gamma_distribution<double> distribution(_alpha, _beta);
		return distribution(sGenerator);
	}

}	//	namespace gpis
