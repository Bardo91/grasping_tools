

#ifndef GPIS_SAMPLER_PRIOR_MEANS_CYLINDERMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_CYLINDERMEAN_H_

#include <Sampler/Prior/Mean.h>
#include <Utils/utilsMath.h>

#include <cassert>

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		//-----------------------------------------------------------------------------------------------------------------
		/// This class represent the mean function of an ellipse with fixed shape parameters. 
		/// For computing its mean, the user need to provide 6 parameters [x,y,z,avx, avy, avz].
		/// The first 3 elements correspond to a translation from the origin to some point on
		/// the semiaxis. The last 3 elements correspond to a vector in the direction of the 
		/// axis of the cylinder.
		class CylinderMean: public Mean {
		public:
			/// Contructor infinite cylinder with semiaxis a and b, being a the largest semiaxis
			/// \ param _a: lenght of X semiaxis
			/// \ param _b: lenght of Y semiaxis
			CylinderMean(double _a, double _b) : mA(_a), mB(_b) {
				mMatrix.zeros(3, 3);
				mMatrix(0, 0) = 1 / (_a*_a);
				mMatrix(1, 1) = 1 / (_b*_b);
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				assert(_parameters.size() == 6);
				// Compute rotation to align cylinder matrix.
				auto R = rotUnitaryVectors(_parameters.subvec(3, 5), {0,0,1});
				auto diff = _point - _parameters.subvec(0, 2);
				return arma::as_scalar((arma::trans(diff)*trans(R)*mMatrix*R*diff)*mA/2);
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// derivatives of the mean value funtion.
			arma::vec meanGradsFun(const arma::vec& _point, const arma::vec &_parameters) const {
				assert(_parameters.size() == 6);
				// Compute rotation to align cylinder matrix.
				auto R = rotUnitaryVectors(_parameters.subvec(3, 5), { 0,0,1 });
				auto diff = _point - _parameters.subvec(0, 2);
				return trans(R)*mMatrix*R*diff*mA;
			}

			double mA, mB;
			arma::mat mMatrix;
		};
	}	//	namespace means
}	//	namespace gpis


#endif	//	GPIS_SAMPLER_PRIOR_MEANS_CYLINDERMEAN_H_