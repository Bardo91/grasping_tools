///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_GRASP_H_
#define GPISGRASPING_GRASP_H_

#include <armadillo>
#include <vector>

#include <grasping_tools/ContactPoint.h>
#include <grasping_tools/Object.h>

namespace grasping_tools {
	class Grasp {
	public:
		Grasp();
		Grasp(const arma::colvec3& _objectCenter, const std::vector<ContactPoint> & _contactPoints);

		void contactPoints(std::vector<ContactPoint> &_contactPoints);

		void objectCentroid(const arma::colvec3& _objectCentroid);

		/// Get list of contact points
		std::vector<ContactPoint> contactPoints() const;

		/// Compute grasp matrix using current contact points.
		bool computeGraspMatrix();

		/// Get grasp matrix
		arma::mat graspMatrix() const;

		/// Get rank of grasp matrix
		unsigned rank();

		/// Compute quality metrics that take into account the algebraic properties of the Grasp Matrix G.
		bool computeQualityMetricsofG();
		double minSingularValue();
		double areaEllipsoid();
		double graspIsotropyIndex();

		/// Compute the 6-dimensional cone resulting of the union of the friction cones of the contact points.
		/// \param _pointsPerContact: number of points used to approximate the each of the contact friction crones
		arma::mat aprxWrenchCone(unsigned _pointsPerContact);

		/// Determine if the current grasp has force closure.
		bool hasForceClosure();

		/// Compute approximate largest-minimum resisted wrench (lmrw).
		double apprLmrw();

		/// Compute exact largest-minimum resisted wrench (lmrw) from convex hull.
		double lmrw();

		/// Compute probability of force closure
		double probabilityForceClosure(unsigned _samples = 100);

		/// Overrided operator=.
		void operator=(const Grasp &_grasp);

	private:
		void checkForceClosure();

	private:
		arma::colvec3				mObjectCenter;
		std::vector<ContactPoint>	mContactPoints;
		arma::mat mGraspMatrix;

		unsigned mRank;

		double mMinSingularValueG;
		double mMaxSingularValueG;
		double mAreaG;
		double mGraspIsotropyIndexG;

		const double cEpsilonDistance = 1e-5;
		arma::mat mConvexConeGrasp;
		bool mHasForceClosure;
		double mLmrw;

		bool mComputedGraspMatrix		= false;
		bool mComputedQualityMetrics	= false;
		bool mComputedForceClosure		= false;
		bool mComputedConvexCone		= false;
		bool mComputedApprxLmrw			= false;
		bool mComputedLmrw			= false;
	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_GRASP_H_
