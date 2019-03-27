

#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_PLANEPARAMETERSAMPLER_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLERS_PLANEPARAMETERSAMPLER_H_

#include <gpis/sampler/prior/ParameterSampler.h>

namespace gpis {
	namespace parameterSamplers {
		class PlaneParameterSampler {// 666 TODO implement.
		public:
			arma::vec operator()(const std::vector<Part> &allParts, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &_oldParameters) {
				// update the objects. at this point the parts are associated to some objects.
				// we look at this one object and find the parts associated with it.
				// based on the parts we sample the parameters of that object.

				unsigned numDims = allParts[0].location().size();
				unsigned numDimsPlusOne = numDims + 1;
				unsigned numAssocParts = assocParts.n_elem;

				arma::mat rHSMat = arma::zeros<arma::mat>(numDimsPlusOne *  numAssocParts, numDimsPlusOne);
				arma::vec muTilde = arma::zeros<arma::vec>(numDimsPlusOne * numAssocParts);
				for (int n = 0; n < numAssocParts; ++n) {
					rHSMat.row((n * numDimsPlusOne)) = join_vert(allParts[assocParts[n]].location(), arma::vec{ 1 }).t();
					for (int d = 0; d < numDims; ++d) {
						rHSMat.at((n * numDimsPlusOne) + d + 1, d) = 1;
						muTilde((n * numDimsPlusOne) + d + 1) = allParts[assocParts[n]].surfaceNormals(d);
					}
				}

				arma::mat P = solve(trimatl(lCholDec), rHSMat);
				arma::mat Q = solve(trimatl(lCholDec), muTilde); // TODO: recycle solution
				arma::vec B = trans(P) * Q;

				arma::mat planeParamPrecMat = trans(P) * P;
				arma::mat planeParamCovMat = inv(planeParamPrecMat);
				arma::vec planeParamMean = planeParamCovMat * B;

				// resample the parameters of the object
				arma::vec nVar(numDimsPlusOne);
				multiVarGaussSample(nVar);
				//nVar.zeros();
                arma::vec planeParams = planeParamMean + arma::chol(planeParamCovMat, "lower") * nVar;
				planeParams /= arma::as_scalar(arma::norm(planeParams));
				return planeParams;

			}
		};
	}	//	namespace parameterSamplers
}	// namespace gpis

#endif
