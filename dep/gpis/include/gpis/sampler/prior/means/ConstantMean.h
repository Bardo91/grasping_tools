
#ifndef GPIS_SAMPLER_PRIOR_MEANS_CONSTANTMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_CONSTANTMEAN_H_


#include <gpis/sampler/prior/Mean.h>

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		//-----------------------------------------------------------------------------------------------------------------
		/// Implementation of sphere mean prior type without fixed radius, parameters are p = {x, y, z, r}
		class ConstantMean : public Mean {
		public:
			/// Contructor
			/// \ param _radius: radius of the prior
			ConstantMean(double _mean) : mean_(_mean) {}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				return mean_;
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// derivatives of the mean value funtion.
			arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const {
				return {0,0,0};
			}

		private:
			double mean_;
		}; //class SphereVarRadiusMean
	}	//	namespace means
}	//	namespace gpis


#endif	//	GPIS_SAMPLER_PRIOR_MEANS_CONSTANTMEAN_H_
