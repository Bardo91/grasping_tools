#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREPARAMETERSAMPLERMH_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_SPHEREPARAMETERSAMPLERMH_H_

#include <Sampler/Prior/ParameterSampler.h>
#include <Utils/utilsProb.h>

#include <Sampler/Prior/means/SphereMean.h>

namespace gpis {
	namespace parameterSamplers {
		class SphereParametersSamplerMH {
		public:
			SphereParametersSamplerMH(double radius) : radius_(radius) { }
			arma::vec operator()(const std::vector<Part> &partList, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &oldParameters) {
				if (oldParameters.n_elem == 0) { // First sample of parameters, init.
					unsigned index = double(rand()) / RAND_MAX * assocParts.n_elem;

					auto &part = partList[assocParts[index]];
					return part.location() - radius_*part.surfaceNormals();
				}
				else {

					// Sample a new set of parameters
					arma::vec newParams(oldParameters.n_elem);
					multiVarGaussSample(newParams);
					newParams *= radius_/10;
					newParams += oldParameters;

					// Compute de likelihoods
					means::SphereMean mean(radius_);
					const unsigned dimensions = partList[0].location().size();
					arma::vec deviationOld(assocParts.n_elem*(dimensions + 1));
					arma::vec deviationNew(assocParts.n_elem*(dimensions + 1));

					for (unsigned i = 0; i < assocParts.n_elem; i++) {
						arma::vec observedVals(4);
						observedVals.head(1) = 0;
						observedVals.tail(dimensions) = partList[assocParts(i)].surfaceNormals();

						deviationOld.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(partList[assocParts(i)].location(), oldParameters, true);
						deviationNew.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(partList[assocParts(i)].location(), newParams, true);
					}

					auto covariance = trans(lCholDec)*lCholDec;	// 666 TODO: not needed but need to modify the interface.
					double oldLL = -arma::as_scalar(0.5*trans(deviationOld)*trans(solve(covariance, deviationOld)));
					double newLL = -arma::as_scalar(0.5*trans(deviationNew)*trans(solve(covariance, deviationNew)));


					// Compute acceptance ratio
					double alpha = exp(newLL)/ exp(oldLL);
					if (alpha > 1) {
						return newParams;
					}
					else {
                        double u = (rand() / (double)(RAND_MAX));
						if (u < alpha)
							return newParams;
						else
							return oldParameters;
					}
				}
			}
		protected:
			double radius_;
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif
