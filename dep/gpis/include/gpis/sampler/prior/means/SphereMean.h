
#ifndef GPIS_SAMPLER_PRIOR_MEANS_SPHEREMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_SPHEREMEAN_H_

#include <gpis/sampler/prior/Mean.h>

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		//-----------------------------------------------------------------------------------------------------------------
		/// Implementation of sphere mean prior type. This prior is constructed with fixed radius and the evaluation parameters
		/// are the centroid's coordinates.
		class SphereMean : public Mean {
		public:
			/// Contructor
			/// \ param _radius: radius of the prior
			SphereMean(double _radius) : radius_(_radius) {}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				auto diff = _point - _parameters;
				return (arma::as_scalar(arma::trans(diff)*diff) - radius_*radius_) / 2 / radius_;
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// derivatives of the mean value funtion.
			arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const {
				auto diff = _point - _parameters;
				return diff / radius_;
			}

			double radius_;
		};
	}	//	namespace means
}	//	namespace gpis


#endif
