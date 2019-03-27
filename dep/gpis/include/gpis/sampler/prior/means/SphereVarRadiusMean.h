
#ifndef GPIS_SAMPLER_PRIOR_MEANS_SPHERERADIUSMHMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_SPHERERADIUSMHMEAN_H_


#include <gpis/sampler/prior/Mean.h>

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		//-----------------------------------------------------------------------------------------------------------------
		/// Implementation of sphere mean prior type without fixed radius, parameters are p = {x, y, z, r}
		class SphereVarRadiusMean : public Mean {
		public:
			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				auto position = _parameters.subvec(0, 2);
				auto &radius = _parameters[3];
				auto diff = _point - position;
				return (arma::as_scalar(arma::trans(diff)*diff) - radius*radius) / 2 / radius;
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// derivatives of the mean value funtion.
			arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const {
				auto position = _parameters.subvec(0, 2);
				auto &radius = _parameters[3];
				auto diff = _point - position;
				return diff / radius;
			}
		}; //class SphereVarRadiusMean
	}	//	namespace means
}	//	namespace gpis


#endif	//	GPIS_SAMPLER_PRIOR_MEANS_SPHERERADIUSMHMEAN_H_
