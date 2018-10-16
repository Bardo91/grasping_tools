
#include <gpis/sampler/prior/GPISPrior.h>
#include <gpis/utils/utilsProb.h>
#include <gpis/utils/utilsGPIS.h>
#include <cassert>

using namespace arma;
using namespace std;

namespace gpis {
	//--------------------------------------------------------------------------------------------------------------------
	GPISPrior::GPISPrior(	Kernel *kernel,
							Mean *mean,
							function<double(vec &)> &geomLikelihood,
							ParameterSampler &parameterSampler):
								kernel_(kernel),
								mean_(mean),
								geomLikelihood_(geomLikelihood),
								parameterSampler_(parameterSampler){

	}
	
	//--------------------------------------------------------------------------------------------------------------------
	bool GPISPrior::init(const vector<Part>& parts, bool useGradients) {
		if (parts.size() == 0)
			return false;

		mat points(parts[0].location().n_elem, parts.size());
		for (unsigned i = 0; i < parts.size(); i++) {
			points.col(i) = parts[i].location();
		}

		fullCovMatrix_ = (*kernel_)(points, points, useGradients);
		return true;
	}

	//--------------------------------------------------------------------------------------------------------------------
	vec GPISPrior::computeDataDeviation(const vector<Part> &parts,
											  const uvec indexList,
											  const vec &_priorParameters) {
		assert(parts.size() != 0);

		const unsigned dimensions = parts[0].location().size();
		vec deviation(indexList.n_elem*(dimensions+1));

		for(unsigned i = 0 ; i < indexList.n_elem; i++){
			vec observedVals(4);
			observedVals.head(1) = 0;
			observedVals.tail(parts[0].location().n_elem) = parts[indexList(i)].surfaceNormals();

			deviation.subvec(i*(dimensions + 1), i*(dimensions + 1) + dimensions) = observedVals - (*mean_)(parts[indexList(i)].location(), _priorParameters, true);
		}

		return deviation;
	}

	//--------------------------------------------------------------------------------------------------------------------
	double GPISPrior::geometricLogLikelihood(const vector<Part> &parts,
											 const uvec indexList,
											 const vec &pose) {
		double logLikelihood = 0.0;
		for(unsigned i = 0 ; i < indexList.n_elem; i++){
			vec loc = parts[i].location();
			logLikelihood += log(geomLikelihood_(loc));	// 666 TODO: prior parameters!
		}
		return logLikelihood;
	}

	//--------------------------------------------------------------------------------------------------------------------
	vec  GPISPrior::sampleParameters(const vector<Part> &allParts, const arma::uvec &assocParts,const arma::mat &lCholDec, const arma::vec &oldParameters) {
		return parameterSampler_(allParts, assocParts, lCholDec, oldParameters);
	}    
	
}   //namespace gpis
