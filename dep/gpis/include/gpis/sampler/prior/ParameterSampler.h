
#ifndef GPIS_SAMPLER_PRIOR_PARAMETERSAMPLER_H_
#define GPIS_SAMPLER_PRIOR_PARAMETERSAMPLER_H_

#include <functional>
#include <armadillo>
#include <gpis/sampler/part/Part.h>
#include <gpis/utils/utilsProb.h>

namespace gpis {
	typedef std::function<arma::vec(const std::vector<Part> &, const arma::uvec &, const arma::mat &, const arma::vec &)> ParameterSampler;
	
}

#endif //GPIS_SAMPLER_PRIOR_PARAMETERSAMPLER_H_
