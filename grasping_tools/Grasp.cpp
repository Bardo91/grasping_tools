///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Grasp.h"
#include "mathTools.h"
#include "qualityMetrics/DistancePointConvexCone.h"
#include "qualityMetrics/IterativeGraspWrenchSpace.h"

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/RboxPoints.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullQh.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullLinkedList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullSet.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <cassert>

namespace gpisGrasping {
	//-----------------------------------------------------------------------------------------------------------------
	Grasp::Grasp() {
		mObjectCenter = arma::zeros(3, 1);
	}
	//-----------------------------------------------------------------------------------------------------------------
	Grasp::Grasp(const arma::colvec3 & _objectCenter, const std::vector<ContactPoint>& _contactPoints){
		mObjectCenter = _objectCenter;
		mContactPoints = _contactPoints;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void Grasp::contactPoints(std::vector<ContactPoint>& _contactPoints) {
		mContactPoints = _contactPoints;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void Grasp::objectCentroid(const arma::colvec3 & _objectCentroid) {
		mObjectCenter = _objectCentroid;
	}


	//-----------------------------------------------------------------------------------------------------------------
	std::vector<ContactPoint> Grasp::contactPoints() const{
		return mContactPoints;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool Grasp::computeGraspMatrix(){
		// If currently computed return true;
		if (mComputedGraspMatrix)
			return true;

		// Need at least 1 contact point
		if (mContactPoints.size() == 0) {
			return false;
		}
		else {
			int colsSelectionMatrix = mContactPoints[1].selectionMatrix().n_cols;
			mGraspMatrix = arma::mat(6, colsSelectionMatrix*mContactPoints.size());
			for (unsigned i = 0; i < mContactPoints.size(); i++) {
				auto &cp = mContactPoints[i];
				arma::vec3 r = cp.position() - mObjectCenter;
				
				arma::mat S = { {0, -r(2), r(1)},
								{r(2), 0, -r(0)},
								{-r(1), r(0), 0} };
				arma::mat R = cp.frame();

				arma::mat Ai = arma::zeros(6, 6);
				Ai.submat(0, 0, 2, 2) = R;
				Ai.submat(3, 0, 5, 2) = -S*R;
				Ai.submat(3, 3, 5, 5) = R;

				arma::mat Gi = Ai*cp.selectionMatrix(); 

				mGraspMatrix.submat(0, i*colsSelectionMatrix, 5, i*colsSelectionMatrix + colsSelectionMatrix - 1) = Gi;
			}

			// 666 TODO: check efficient ways
			mRank = arma::rank(mGraspMatrix);

			mComputedGraspMatrix = true;
			mComputedQualityMetrics = false;
			return true;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	arma::mat Grasp::graspMatrix() const{
		return mGraspMatrix;
	}

	//-----------------------------------------------------------------------------------------------------------------
	unsigned Grasp::rank() {
		assert(mComputedGraspMatrix);
		return mRank;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool Grasp::computeQualityMetricsofG(){
		// If currently computed return true;
		if (mComputedQualityMetrics)
			return true;

		if (mContactPoints.size() == 0)
			return false;

		if (!mComputedGraspMatrix) {
			computeGraspMatrix();
			mComputedGraspMatrix = true;
		}

		// Check rank G, if rank(G) < 6, then grasp is not stable.
		if (mRank < 6) {
			return false;
		}
		else {
			// Compute SVD of G.
			arma::vec singularValues = arma::svd(mGraspMatrix);

			mMinSingularValueG = arma::min(singularValues);
			mMaxSingularValueG = arma::max(singularValues);
			mAreaG = arma::prod(singularValues);
			mGraspIsotropyIndexG = mMinSingularValueG/ mMaxSingularValueG;

			mComputedQualityMetrics = true;
			return true;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::minSingularValue(){
		return mMinSingularValueG;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::areaEllipsoid(){
		return mAreaG;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::graspIsotropyIndex(){
		return mGraspIsotropyIndexG;
	}

	//-----------------------------------------------------------------------------------------------------------------
	arma::mat Grasp::aprxWrenchCone(unsigned _pointsPerContact) {
		if (!mComputedConvexCone) {
			for (auto &cp: mContactPoints) {
				auto wrenches = cp.aprxWrenchCone(_pointsPerContact);
				mConvexConeGrasp.insert_cols(mConvexConeGrasp.n_cols, wrenches);
			}

			mComputedConvexCone = true;
		}
		return mConvexConeGrasp;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool Grasp::hasForceClosure() {
		if (mContactPoints.size() == 0)
			return false;

		if (!mComputedForceClosure) {
			checkForceClosure();
		}

		return mHasForceClosure;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::apprLmrw() {
		if (mComputedApprxLmrw) {
			return mLmrw;
		}
		else {
			if (mContactPoints.size() == 0)
				return false;

			IterativeGraspWrenchSpace gws(*this, IterativeGraspWrenchSpace::eForceLimitType::TotalSource);
			gws.epsilon(1e-3);
			mLmrw = gws.lmrw();

			if (mLmrw == 0.0) {
				mHasForceClosure = false;
			}
			else {
				mHasForceClosure = true;
			}

			mComputedForceClosure = true;
			mComputedApprxLmrw = true;

			return mLmrw;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::lmrw(){
		mComputedLmrw = false;
		if (mComputedLmrw) {
			return mLmrw;
		}
		else {
			if (mContactPoints.size() == 0)
				return false;

			// Compute wrenches
			auto wrenchesCones = aprxWrenchCone(8);

			// Compute convex hull

			orgQhull::Qhull convexHull;
			orgQhull::PointCoordinates *points = new orgQhull::PointCoordinates(6);
			std::vector<double> concatenationPoints;
			for (unsigned i = 0; i < wrenchesCones.n_cols; i++) {
				
				concatenationPoints.push_back(wrenchesCones.col(i)[0]);
				concatenationPoints.push_back(wrenchesCones.col(i)[1]);
				concatenationPoints.push_back(wrenchesCones.col(i)[2]);
				concatenationPoints.push_back(wrenchesCones.col(i)[3]);
				concatenationPoints.push_back(wrenchesCones.col(i)[4]);
				concatenationPoints.push_back(wrenchesCones.col(i)[5]);
				
			}
			points->append(concatenationPoints);

			try {
				convexHull.runQhull(*points, "Qt");
			}
			catch (orgQhull::QhullError e) {
				std::cout << "Error computing convex hull. Error msg: \n" /*<< e.what() << std::endl*/ << std::endl;
				mLmrw = 0;
				mHasForceClosure = false;
				return mLmrw;
			}

			orgQhull::QhullFacetList facets = convexHull.facetList();
			double minDistance = 99999999;
			int minIdx = 0;
			int counterIdx = 0;
			orgQhull::PointCoordinates originCoordinates(6);
			std::vector<double> originstd = { 0.0,0.0,0.0,0.0,0.0,0.0, };
			originCoordinates.append(originstd);

			orgQhull::QhullPoint origin6d(6, &originstd[0]);
			//std::cout << "Num of facets of CH is " << facets.size() << std::endl;
			for (orgQhull::QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it) {
				orgQhull::QhullHyperplane hyperPlane = (*it).hyperplane();
				double distance = hyperPlane.distance(origin6d); 
				if (distance > 0) {
					//std::cout << "ERROR! Origin is not contained in the QH" << std::endl;
					minDistance = 0.0;
					break;
				}else{
					distance *= -1;
				}

				// From Qhull, negative distances means that point is inside and positive that is outside.
				if (distance < minDistance) {
					minDistance = distance;
					minIdx = counterIdx;
				}
				counterIdx++;
			}

			mLmrw = minDistance;

			if (mLmrw <= 0.0) {
				mHasForceClosure = false;
			}
			else {
				mHasForceClosure = true;
			}

			mComputedForceClosure = true;
			mComputedLmrw = true;

			return mLmrw;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	double Grasp::probabilityForceClosure(unsigned _samples) {
		double pqf = 0.0;

		for (unsigned i = 0; i < _samples; i++) {
			auto cps = this->contactPoints();
			// Sample cps
			std::vector<ContactPoint> newCps;
			for (auto cp : cps) {
				newCps.push_back(cp.sample());
			}

			Grasp grasp;
			grasp.contactPoints(newCps);
			if (grasp.lmrw() > 0) {
				pqf++;
			}
		}

		return pqf/_samples;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void Grasp::operator=(const Grasp & _grasp) {
		mContactPoints		= _grasp.mContactPoints;
		mGraspMatrix		= _grasp.mGraspMatrix;
		
		mRank				= _grasp.mRank;
		
		mMinSingularValueG		= _grasp.mMinSingularValueG;
		mMaxSingularValueG		= _grasp.mMaxSingularValueG;
		mAreaG					= _grasp.mAreaG;
		mGraspIsotropyIndexG	= _grasp.mGraspIsotropyIndexG;
		
		mConvexConeGrasp		= _grasp.mConvexConeGrasp;
		mHasForceClosure		= _grasp.mHasForceClosure;
		mLmrw					= _grasp.mLmrw;
		
		mComputedGraspMatrix	= _grasp.mComputedGraspMatrix;
		mComputedQualityMetrics	= _grasp.mComputedQualityMetrics;
		mComputedForceClosure	= _grasp.mComputedForceClosure;
		mComputedConvexCone		= _grasp.mComputedConvexCone;
		mComputedApprxLmrw		= _grasp.mComputedApprxLmrw;
		mComputedLmrw			= _grasp.mComputedLmrw;
	}

	//-----------------------------------------------------------------------------------------------------------------
	void Grasp::checkForceClosure() {
		if (!mComputedConvexCone)
			aprxWrenchCone(8);


		if (!mComputedGraspMatrix) {
			computeGraspMatrix();
			mComputedGraspMatrix = true;
		}

		arma::colvec ws = arma::zeros(6, 1);
		for (auto cp : mContactPoints) {
			arma::colvec wsi = arma::zeros(6, 1);
			wsi.rows(0, 2) = cp.normal();
			arma::vec p = cp.position();
			arma::vec n = cp.normal();
			wsi.rows(3, 5) = cross(p,n);
			ws += wsi;
		}

		// d(ws, coneW)
		DistancePointConvexCone distCone;
		distCone.basePoints(mConvexConeGrasp);
		distCone.epsilon(1e-6);
		distCone.point(-ws);

		if (distCone.compute())
			mComputedForceClosure = true;

		mHasForceClosure = distCone.minDistance() < cEpsilonDistance;

	}
}
