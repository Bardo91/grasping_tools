
#ifndef GPIS_SAMPLER_PRIOR_KERNEL_H_
#define GPIS_SAMPLER_PRIOR_KERNEL_H_

#include <functional>
#include <armadillo>

namespace gpis{
	class Kernel{
	public:
		arma::mat operator()(const arma::mat &points1, bool computeGradients = false) const;
		arma::mat operator()(const arma::mat &points1, const arma::mat &points2, bool computeGradients = false) const;

	public:
		virtual double covarianceValVal(const arma::vec& _point1, const arma::vec& _point2, const bool isSameObs) const = 0;
		virtual arma::vec covarianceValDer(const arma::vec& _point1, const arma::vec& _point2) const = 0;
		virtual arma::mat covarianceDerDer(const arma::vec& _point1, const arma::vec& _point2, const bool isSameObs) const = 0;
	};
}

#endif
