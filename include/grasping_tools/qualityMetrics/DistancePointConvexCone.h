///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_QUALITYMETRICS_DISTANCEPOINTCONVEXCONE_H_
#define GPISGRASPING_QUALITYMETRICS_DISTANCEPOINTCONVEXCONE_H_

#include <armadillo>

namespace grasping_tools {
	/// Based on Distance Between a Point and a Convex Cone in n - Dimensional  Space : Computation and Applications, 
	/// by Yu Zheng and Chee Meng Chew.
	class DistancePointConvexCone {
	public:		// Public interface
		/// Set the compact set of points used to define the convex cone from the Origin.
		/// \params _A: Matrix containing points on each column. Size: n x (num points).
		void basePoints(const arma::mat &_A);

		/// Set point to compute the distance.
		/// \params _b: n-dimensional column vector.
		void point(const arma::colvec &_b);

		/// Set tolerance epsilon.
		/// \param _epsilon:
		void epsilon(double _epsilon);

		/// Get result minimun distance.
		double minDistance() const;

		/// Get closest point in the cone
		arma::colvec closestPoint() const;
		
		/// Get reduced set
		arma::mat baseSet() const;

		/// Perform the algorithm described in ther article to compute the minimun distance
		/// between the given point and the given convex cone.
		bool compute();

		/// Set maximum number of iterations. Default 100.
		void maxIters(unsigned _iters);

	private:	// Private methods
		bool checkInputData();

		double supportFunction(const arma::colvec &_r);
		arma::colvec mappingFunction(const arma::colvec &_r);

		// Subalgorithm described in the article
		void updateSubConvexCone();

	private:	// Members
		arma::mat mA;
		arma::colvec mB;

		double mEpsilon;
		unsigned mMaxIters = 100;
		double mMinDistance = std::numeric_limits<double>::max();
		arma::colvec mClosestPoint;

		// Temporal variables used on the algorithm shared between the functions
		arma::mat mAi, mAmi;
		arma::colvec mV, mR;
	};
}

#endif