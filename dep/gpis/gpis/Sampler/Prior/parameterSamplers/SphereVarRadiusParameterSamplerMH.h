#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREVARRADIUSPARAMETERSAMPLERMH_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREVARRADIUSPARAMETERSAMPLERMH_H_

#include <Sampler/Prior/ParameterSampler.h>
#include <Utils/utilsProb.h>

#include <Sampler/Prior/means/SphereVarRadiusMean.h>

#include <random>

namespace gpis {
	namespace parameterSamplers {
		class SphereVarRadiusParametersSamplerMH {
		public:
			SphereVarRadiusParametersSamplerMH(double _alpha, double _beta) : mAlpha(_alpha), mBeta(_beta) { }
			arma::vec operator()(const std::vector<Part> &partList, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &oldParameters) {
				if (oldParameters.n_elem == 0) { // First sample of parameters, init.
					
					double sampledRadius = sampleGammaDist(mAlpha, mBeta);

					unsigned randomIndexPart = double(rand()) / RAND_MAX * assocParts.n_elem;

					auto &part = partList[assocParts[randomIndexPart]];

					arma::vec params(4);
					params.subvec(0, 2) = part.location() - sampledRadius*part.surfaceNormals();;
					params[3] = sampledRadius;

					return params;
				}
				else {

					// Sample for the position using previous radius
					arma::vec newParams(oldParameters.n_elem);
					multiVarGaussSample(newParams);
					newParams *= oldParameters[3] * 0.2;
					newParams += oldParameters;
					newParams[3] = oldParameters[3];	// Set old radius

					arma::vec newOldParams = mhStep(partList, assocParts, newParams, oldParameters, lCholDec);

					// Sample for the radius
					newParams = newOldParams;
					arma::vec newRadius(1);
					multiVarGaussSample(newRadius);
					newRadius *= newOldParams[3]*0.25;
					newRadius += newOldParams[3];
					newParams[3] = newRadius[0];

					return mhStep(partList, assocParts, newParams, newOldParams, lCholDec);
				}
			}

		private:
			arma::vec mhStep(const std::vector<Part> &partList, const arma::uvec &assocParts, const arma::vec &newParams, const arma::vec &oldParams, const arma::mat &lCholDec) {
				// Compute de likelihoods
				means::SphereVarRadiusMean mean;
				const unsigned dimensions = partList[0].location().size();
				arma::vec deviationOld(assocParts.n_elem*(dimensions + 1));
				arma::vec deviationNew(assocParts.n_elem*(dimensions + 1));

				for (unsigned i = 0; i < assocParts.n_elem; i++) {
					arma::vec observedVals(4);
					observedVals.head(1) = 0;
					observedVals.tail(dimensions) = partList[assocParts(i)].surfaceNormals();

					deviationOld.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(partList[assocParts(i)].location(), oldParams, true);
					deviationNew.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(partList[assocParts(i)].location(), newParams, true);
				}

				auto covariance = trans(lCholDec)*lCholDec;	// 666 TODO: not needed but need to modify the interface.
				double oldLL = -arma::as_scalar(0.5*trans(deviationOld)*trans(solve(covariance, deviationOld)));
				double newLL = -arma::as_scalar(0.5*trans(deviationNew)*trans(solve(covariance, deviationNew)));

				// Weight probabilities unsing gamma as prior knowledge about the radius
				double gammaOld = 1;//evalGammaDist(oldParams[3], mAlpha, mBeta);
				double gammaNew = 1;//evalGammaDist(newParams[3], mAlpha, mBeta);

				// Compute acceptance ratio
				double alpha = (gammaNew * exp(newLL)) / (gammaOld * exp(oldLL));

				if (alpha > 1) {
					return newParams;
				}
				else {
                    double u = (rand() / (double)(RAND_MAX));
					if (u < alpha)
						return newParams;
					else
						return oldParams;
				}
			}

		private:
			// Parameters of the gamma distribution for the radius of the object.
			double mAlpha, mBeta;
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif	//	GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREVARRADIUSPARAMETERSAMPLERMH_H_
