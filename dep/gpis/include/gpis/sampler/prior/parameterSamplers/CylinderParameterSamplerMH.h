#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_CYLINDERPARAMETERSAMPLERMH_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_CYLINDERPARAMETERSAMPLERMH_H_

#include <gpis/sampler/prior/ParameterSampler.h>
#include <gpis/utils/utilsProb.h>

#include <gpis/sampler/prior/means/CylinderMean.h>

namespace gpis {
	namespace parameterSamplers {
		class CylinderParametersSamplerMH {
		public:
			CylinderParametersSamplerMH(double _a, double _b) : mA(_a), mB(_b) {
				mMatrix.zeros(3, 3);
				mMatrix(0, 0) = 1 / (_a*_a);
				mMatrix(1, 1) = 1 / (_b*_b);
			}

			arma::vec operator()(const std::vector<Part> &allParts, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &oldParameters) {
				if (oldParameters.n_elem == 0) { // First sample of parameters, init.
					unsigned index = double(rand()) / RAND_MAX * assocParts.n_elem;	// Pick random part

					auto &part = allParts[assocParts[index]];

					arma::vec parameters(6);
					parameters.subvec(0, 2) = part.location() - mA*part.surfaceNormals();
					arma::vec randVec;
					randVec.randu(3);
					parameters.subvec(3, 5) = arma::cross(randVec, part.surfaceNormals());
					parameters.subvec(3, 5) /= norm(parameters.subvec(3, 5));
					return parameters;
				}
				else {

					// Sample a new set of parameters
					arma::vec newParams(oldParameters.n_elem);
					multiVarGaussSample(newParams);
					newParams.subvec(0, 2) *= mA;
					newParams.subvec(3, 5) *= mA;
					newParams += oldParameters;
					newParams.subvec(3, 5) /= norm(newParams.subvec(3, 5));


					// Compute de likelihoods
					means::CylinderMean mean(mA,mB);
					const unsigned dimensions = allParts[0].location().size();
					arma::vec deviationOld(assocParts.n_elem*(dimensions + 1));
					arma::vec deviationNew(assocParts.n_elem*(dimensions + 1));

					for (unsigned i = 0; i < assocParts.n_elem; i++) {
						arma::vec observedVals(4);
						observedVals.head(1) = 0;
						observedVals.tail(dimensions) = allParts[assocParts(i)].surfaceNormals();

						deviationOld.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(allParts[assocParts(i)].location(), oldParameters, true);
						deviationNew.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(allParts[assocParts(i)].location(), newParams, true);
					}

					double oldLL = -arma::as_scalar(0.5*trans(deviationOld)*solve(trimatu(trans(lCholDec)), solve(trimatl(lCholDec), deviationOld)));
					double newLL = -arma::as_scalar(0.5*trans(deviationNew)*solve(trimatu(trans(lCholDec)), solve(trimatl(lCholDec), deviationNew)));


					// Compute acceptance ratio
					double alpha = exp(newLL) / exp(oldLL);
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

		private:
			double mA, mB;
			arma::mat mMatrix;
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif	//	GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_CYLINDERPARAMETERSAMPLERMH_H_
