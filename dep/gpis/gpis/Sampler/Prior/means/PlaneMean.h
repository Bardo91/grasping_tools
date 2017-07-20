

#ifndef GPIS_SAMPLER_PRIOR_MEANS_PLANEMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_PLANEMEAN_H_

#include <Sampler/Prior/Mean.h>
namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		/// Implementation of sphere mean prior type. This prior doesn't need any extra information. The evaluation parameters are
		/// the [a,b,c,d] parameters for the plane equation "ax+by+cz+d=0".
		class PlaneMean : public Mean {

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				return _parameters(0)*_point(0) + _parameters(1)*_point(1) + _parameters(2)*_point(2) + _parameters(3);
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// derivatives of the mean value funtion.
			arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const {
				return{ _parameters(0), _parameters(1), _parameters(2) };
			}
		};
	}	//	namespace means
}	//	namespace gpis


#endif	//	GPIS_SAMPLER_PRIOR_MEANS_PLANEMEAN_H_