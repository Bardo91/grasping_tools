
#ifndef GPIS_SAMPLER_PRIOR_PRIORFACTORY_H_
#define GPIS_SAMPLER_PRIOR_PRIORFACTORY_H_

#include "GPISPrior.h"

namespace gpis{
	class PriorFactory{
	public:	// Public interface
		enum class eShapeTypes {	Sphere,			/// Spherical prior. Mean function -> m(x) = ((x-centroid)'*(x-centroid) - R*R)/2/R. Parameters: Radio (R).
									SphereMH,		/// Spherical prior. Mean function -> m(x) = ((x-centroid)'*(x-centroid) - R*R)/2/R. Parameters: Radio (R).
													/// Use a Metropolis hasting alhorithm for the sampler
									SphereRadiusMH,	/// Spherical prior with unknow radius but prior knowledge using gamma distribution. 
													/// Mean function -> m(x) = ((x-centroid)'*(x-centroid) - R*R)/2/R. Parameters: alpha y beta parameters of gamma distribution.
													/// Use a Metropolis hasting alhorithm for sampling all the parameters
									Plane,			/// Planar prior. Mean function -> m(x) = ax + by +cz + d. Parameters: no parameters 
									EllipseMH,		/// Ellipsoidal prior. Mean function -> m(x) = (x-c)RAR(x-c). Parameters: X, Y and Z semiaxis. 
                                    CylinderMH,		/// Cylindrical prior. Mean function -> m(x) = (x-c)RAR(x-c). Parameters: X and Y semiaxis.
                                    Constant        /// Constant mean. Mean function -> m(x) = c. Parameters: c mean value.
								};

		enum class eKernelTypes {	Linear,				/// Linear kernel. k(x,y) = x'*y. Parameters: empty.
									SquaredExponential	/// Squared Exponential kernel. k(x,y) = sigma^2*exp(-gamma*(x-y)'*(x-y)). Parameters: Sigma, gamma.
								};


		/// Create a new prior with the given information.
		///	\param shapeType: enum with desired shape for the prior function. It is used to configure the mean and the geometrical function.
		///	\param shapeParams:	vector containing the variables for the shape definition. Initiallizer list can be used.
		///	\param kernelType: enum with the desired kernel.
		/// \param KernelParams: vector containing the variables for the kernel definition. Initiallizer list can be used.
		///
		/// Example:
		/// \code
		///		auto prior = create(eShapeTypes::Sphere, {0.3}, SquaredExponential, {0.1,0.1,0.0}
		///	\endcode
		static GPISPrior* create(eShapeTypes shapeType, const arma::vec &shapeParams, eKernelTypes kernelType, const arma::vec &kernelParams);

	private:	// Private interface
	};
	
}

#endif
