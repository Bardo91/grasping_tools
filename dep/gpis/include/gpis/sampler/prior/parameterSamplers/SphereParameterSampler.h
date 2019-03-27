#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREPARAMETERSAMPLER_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREPARAMETERSAMPLER_H_

#include <gpis/sampler/prior/ParameterSampler.h>

namespace gpis {
	namespace parameterSamplers {
		class SphereParametersSampler {
		public:
			SphereParametersSampler(double radius) : radius_(radius) { }
			arma::vec operator()(const std::vector<Part> &allParts, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &_oldParameters) {
				// compute the mean and covariance matrix for the location distribution for this object
				unsigned numDims = allParts[0].location().size();
				unsigned numDimsPlusOne = numDims + 1;
				unsigned numAssocParts = assocParts.n_elem;

				arma::mat eyeRepMat = arma::zeros<arma::mat>(numDims *numAssocParts, numDims);
				arma::vec muTilde(numDims *numAssocParts);
				arma::uvec indicesCholDec(numDims *numAssocParts);
				for (int n = 0; n <numAssocParts; ++n) {
					for (int d = 0; d < numDims; ++d) {
						eyeRepMat.at((n * numDims) + d, d) = 1;
						muTilde((n * numDims) + d) = allParts[assocParts[n]].location()(d) -
							radius_* allParts[assocParts[n]].surfaceNormals(d);
						indicesCholDec[n*numDims + d] = n*numDimsPlusOne + d + 1;
					}
				}

				auto lChol = lCholDec(indicesCholDec, indicesCholDec);

				arma::mat P = arma::solve(trimatl(lChol), eyeRepMat);
				arma::mat Q = arma::solve(trimatl(lChol), muTilde);
				arma::vec B = arma::trans(P) * Q / radius_ / radius_;

				arma::mat locPrecMat = arma::trans(P) * P / radius_ / radius_;
				arma::mat locCovMat = arma::inv(locPrecMat);
				arma::vec locMean = locCovMat * B;

				// resample the location of the object
				arma::vec nVar(numDims);
				multiVarGaussSample(nVar);
				return locMean + chol(locCovMat, "lower") * nVar;
			}
		protected:
			double radius_;
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif
