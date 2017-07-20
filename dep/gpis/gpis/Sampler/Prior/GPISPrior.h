#ifndef GPISPRIOR_H
#define GPISPRIOR_H

#include <Sampler/Part/Part.h>
#include <Sampler/Object/MatrixData.h>
#include <Sampler/Prior/Kernel.h>
#include <Sampler/Prior/Mean.h>
#include <Sampler/Prior/ParameterSampler.h>

namespace gpis {
	/// Prior object shape and surface properties
	class GPISPrior {
	public:
		/// Main constructor.
		/// \param parts: list of parts in current scene.
		/// \param kernelParams: kernel parameters.
		GPISPrior(	Kernel *kernel,
					Mean *mean,
					std::function<double(arma::vec &)> &geomLikelihood, 
					ParameterSampler &parameterSampler);

		/// Initialize prior with given parts. Compute the full covariance matrix.
		/// \param parts: list of parts of the scene.
		bool init(const std::vector<Part> &parts, bool useGradients = false);

		/// Given an object pose, computes the deviation of each datapoint from the prior function
		/// \param parts:
		/// \param pose:
		arma::vec computeDataDeviation(const std::vector<Part> &parts,
									   const arma::uvec indexList,
									   const arma::vec &parameters);

		/// Compute the geometric loglikelihood with given data.
		/// \param parts:
		/// \param pose:
		double geometricLogLikelihood(const std::vector<Part> &parts, const arma::uvec indexList, const arma::vec &pose);

		/// Sample a new position depending on the prior behaviour.
		arma::vec sampleParameters(const std::vector<Part> &partList, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &oldParameters);
		
		/// Get the covariance matrix
		const arma::mat &fullCovariance() const { return fullCovMatrix_; }

		/// Print prior info
		void printInfo();

		/// Get kernel
        Kernel *kernel() { return kernel_; }

		/// Get mean
        Mean *mean() { return mean_; }

	protected:
		Kernel *kernel_;
		Mean *mean_;
		std::function<double(arma::vec &)> geomLikelihood_;
		ParameterSampler parameterSampler_;

		arma::mat fullCovMatrix_;
	};
}   // namespace gpis

#endif // GPISPRIOR_H
