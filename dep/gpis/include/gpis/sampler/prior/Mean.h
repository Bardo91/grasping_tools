
#ifndef GPIS_SAMPLER_PRIOR_MEAN_H_
#define GPIS_SAMPLER_PRIOR_MEAN_H_

#include <functional>
#include <armadillo>

namespace gpis{
	class Mean{
	public:
		/// Compute the mean of a set of points..
		///	\param points:
		/// \param computeGradients:
		arma::mat operator()(const arma::mat &_points, const arma::vec &_parameters, bool computeGradients = false) const;
	
		/// Abstract method. Given a point and the parameters, it computes the mean value funtion.
		virtual double meanFun(const arma::vec& _point, const arma::vec &_parameters) const = 0;

		/// Abstract method. Given a point and the parameters, it computes the derivatives of the mean value funtion.
		virtual arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const = 0;
	};
}	// namespace gpis

#endif	//	GPIS_SAMPLER_PRIOR_MEAN_H_
