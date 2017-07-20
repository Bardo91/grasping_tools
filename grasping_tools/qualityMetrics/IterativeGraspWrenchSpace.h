///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_ITERATIVEGRASPWRENCHSPACE_H_
#define GPISGRASPING_ITERATIVEGRASPWRENCHSPACE_H_

#include <armadillo>
#include <vector>

#include "ContactPoint.h"
#include "Grasp.h"

namespace gpisGrasping {
	/// Based on the article "An efficient Algorithm for a Grasp Quality Measure" by Yu Zheng.
	class IterativeGraspWrenchSpace {
	public:		// Public interface
		enum class eForceLimitType { PerFinger, TotalSource };
	
		/// Contructor. Algorithm based on An Efficient Algorithm for a Grasp Quality Measure Yu Zheng
		IterativeGraspWrenchSpace(const Grasp &_grasp, eForceLimitType _type = eForceLimitType::TotalSource);

		void epsilon(double _epsilon);

		/// Compute the largest-minimum resisted wrench quality metric.
		double lmrw();
	private:	// Private methods
		bool computeInitialSet();
		bool computeFacetsData();
		bool updateFacets();

		double supportFunction(const arma::colvec &_u) const;	// hw
		arma::vec mappingFunction(const arma::colvec &_u) const;	// sw

		// Create an orthogonal vector to a set of vectors
		arma::colvec orthogonalVector(arma::mat _space) const;

	private:	// Members
		Grasp mGrasp;
		eForceLimitType mType;
		double mEpsilon;

		unsigned mWrenchDims;	//	Dimensions of the wrench space
		arma::mat mWrenches;

		arma::mat				mPolytope;	// So called Vk in the article
		arma::umat				mFacets;		// Matrix that contains on each column the indices of the vertes of the polytope that form a facet.

		std::vector<arma::colvec>	mNks;
		std::vector<double>			mDks;
		arma::colvec	mNk;
		double			mDk;

		double mLastLMRW;
	};
}	//	namespace gpisGrasping

#endif //	GPISGRASPING_GRASPWRENCHSPACE_H_
