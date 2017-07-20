//
//
//
//
//

#ifndef GPIS_SAMPLER_PRIOR_KERNELS_SQUAREDEXPONENTIALKERNEL_H_
#define GPIS_SAMPLER_PRIOR_KERNELS_SQUAREDEXPONENTIALKERNEL_H_

#include "../Kernel.h"

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace kernels {
		//-----------------------------------------------------------------------------------------------------------------
		/// Class for making "configurable" kernel function for squared-exponential covariance function (val-val)
		class SquaredExponentialKernel :public Kernel {
		public:
			/// Constructor.
			///	\param sigma:
			///	\param gamma:
			///	\param sigmaNoiseVals:
			/// \param sigmaNoiseGrads:
			SquaredExponentialKernel(double sigma, double gamma, double sigmaNoiseVals = 0, double sigmaNoiseGrads = 0)
				: sigma_(sigma), gamma_(gamma), sigmaNoiseVals_(sigmaNoiseVals), sigmaNoiseGrads_(sigmaNoiseGrads) {}


		public:
			double covarianceValVal(const arma::vec& _point1, const arma::vec& _point2, const bool isSameObs) const {
				arma::vec diff = _point1 - _point2;
				double noise = isSameObs*pow(sigmaNoiseVals_, 2.0);
				return pow(sigma_, 2.0) * exp(as_scalar(-0.5 * gamma_ * arma::trans(diff) * diff)) + noise;
			}

			arma::vec covarianceValDer(const arma::vec& _point1, const arma::vec& _point2) const {
				return gamma_ * (_point1 - _point2) * covarianceValVal(_point1, _point2, false);
			}

			arma::mat covarianceDerDer(const arma::vec& _point1, const arma::vec& _point2, const bool isSameObs) const {
				auto I = arma::eye<arma::mat>(_point1.n_elem, _point1.n_elem);
				auto diff = _point1 - _point2;
				auto noise = isSameObs*pow(sigmaNoiseGrads_, 2.0);
				return	gamma_ * (I - gamma_ * (diff)* arma::trans(diff)) * covarianceValVal(_point1, _point2, false) + I*noise;
			}
		private:
			double sigma_ = 0;
			double gamma_ = 0;
			double sigmaNoiseVals_ = 0;
			double sigmaNoiseGrads_ = 0;
		};
	}
}

#endif