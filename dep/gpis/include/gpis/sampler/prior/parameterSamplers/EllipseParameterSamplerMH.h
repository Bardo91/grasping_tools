#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_ELLIPSEPARAMETERSAMPLER_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_ELLIPSEPARAMETERSAMPLER_H_

#include <gpis/sampler/prior/ParameterSampler.h>
#include <gpis/sampler/prior/means/EllipseMean.h>

namespace gpis {
	namespace parameterSamplers {
		class EllipseParametersSampler {
		public:
			EllipseParametersSampler(double _a, double _b, double _c) : mA(_a), mB(_b), mC(_c) {
				mMatrix.zeros(3, 3);
				mMatrix(0, 0) = 1 / (_a*_a);
				mMatrix(1, 1) = 1 / (_b*_b);
				mMatrix(2, 2) = 1 / (_c*_c);
			}

			arma::vec operator()(const std::vector<Part> &allParts, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &_oldParameters) {
				if (_oldParameters.n_elem == 0) { // First sample of parameters, init.
					unsigned index = double(rand()) / RAND_MAX * assocParts.n_elem;

					auto &part = allParts[assocParts[index]];

					arma::vec parameters(6);
					parameters.subvec(0, 2) = part.location() - mC*part.surfaceNormals();
					arma::vec randVec;
					randVec.randu(3);
					parameters.subvec(3, 5) = part.surfaceNormals();
					parameters.subvec(3, 5) /= norm(parameters.subvec(3, 5));
					return parameters;
				}
				else {
					assert(_oldParameters.n_elem == 6);

					// Sample a new set of parameters
					arma::vec newParams(_oldParameters.n_elem);
					multiVarGaussSample(newParams);
					newParams.subvec(0, 2) *= mC;
					newParams.subvec(3, 5) *= mC/2;
					newParams += _oldParameters;
					newParams.subvec(3, 5) /= norm(newParams.subvec(3, 5));


					// Compute de likelihoods
					means::EllipseMean mean(mA, mB, mC);
					const unsigned dimensions = allParts[0].location().size();
					arma::vec deviationOld(assocParts.n_elem*(dimensions + 1));
					arma::vec deviationNew(assocParts.n_elem*(dimensions + 1));

					for (unsigned i = 0; i < assocParts.n_elem; i++) {
						arma::vec observedVals(4);
						observedVals.head(1) = 0;
						observedVals.tail(dimensions) = allParts[assocParts(i)].surfaceNormals();

						deviationOld.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(allParts[assocParts(i)].location(), _oldParameters, true);
						deviationNew.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - mean(allParts[assocParts(i)].location(), newParams, true);
					}

					double oldLL = arma::as_scalar(-0.5*trans(deviationOld)*solve(trimatu(trans(lCholDec)), solve(trimatl(lCholDec), deviationOld)));
					double newLL = arma::as_scalar(-0.5*trans(deviationNew)*solve(trimatu(trans(lCholDec)), solve(trimatl(lCholDec), deviationNew)));
					
					// Compute acceptance ratio
					double alpha = exp(newLL) / exp(oldLL);
					if (alpha > 1) {
						return newParams;
					}
					else {
						double u = double(rand())/RAND_MAX;
						if (u < alpha) {
							return newParams;
						}
						else {
							return _oldParameters;
						}
					}
				}
			}

		private:
			double mA, mB, mC;
			arma::mat mMatrix;
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif	//	GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_ELLIPSEPARAMETERSAMPLER_H_
