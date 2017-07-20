

#ifndef GPIS_SAMPLER_PRIOR_MEANS_ELLIPSEMEAN_H_
#define GPIS_SAMPLER_PRIOR_MEANS_ELLIPSEMEAN_H_

#include <Sampler/Prior/Mean.h>
#include <Utils/utilsMath.h>

#include <cassert>

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	namespace means {
		//-----------------------------------------------------------------------------------------------------------------
		/// This class represent the mean function of an ellipse with fixed shape parameters. 
		/// For computing its mean, the user need to provide 6 parameters [x,y,z,rx,ry,rz]
		/// The first 3 elements correspond to a translation from the origin to the center of the
		/// ellipse. The last 3 elements correspond to a vector in the direction of the 
		/// major axis of the ellipse.
		class EllipseMean : public Mean{
		public:
			/// Contructor
			/// \param _a: X semi-axis
			/// \param _b: Y semi-axis
			/// \param _c: Z semi-axis
			EllipseMean(double _a, double _b, double _c) : mA(_a), mB(_b), mC(_c) {
				mMatrix.zeros(3, 3);
				mMatrix(0, 0) = 1 / (_a*_a);
				mMatrix(1, 1) = 1 / (_b*_b);
				mMatrix(2, 2) = 1 / (_c*_c);
			}

			/// Definition of abract method in Mean. Given a point and the parameters, it computes the 
			/// mean value funtion.
			double meanFun(const arma::vec& _point, const arma::vec &_parameters)  const {
				assert(_parameters.size() == 6);
				// Compute rotation to align cylinder matrix.
				auto R = rotUnitaryVectors(_parameters.subvec(3, 5), { 0,0,1 });
				auto diff = _point - _parameters.subvec(0, 2);
				return arma::as_scalar((arma::trans(diff)*trans(R)*mMatrix*R*diff -1)*mA / 2);
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

		private:
			double mA, mB, mC;
			arma::mat mMatrix;
		};

	}	//	namespace means
}	//	namespace gpis


#endif	//	GPIS_SAMPLER_PRIOR_MEANS_ELLIPSEMEAN_H_