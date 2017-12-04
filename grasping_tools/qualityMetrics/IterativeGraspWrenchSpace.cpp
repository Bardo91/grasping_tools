///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cassert>

#include "IterativeGraspWrenchSpace.h"
#include "mathTools.h"
#include "DistancePointConvexCone.h"

using namespace arma;
using namespace std;

namespace grasping_tools {
	//-----------------------------------------------------------------------------------------------------------------
	IterativeGraspWrenchSpace::IterativeGraspWrenchSpace(const Grasp & _grasp, eForceLimitType _type) :
		mGrasp(_grasp),
		mType(_type) {
		mWrenchDims = _grasp.graspMatrix().n_rows;
	}
	
	//-----------------------------------------------------------------------------------------------------------------
	void IterativeGraspWrenchSpace::epsilon(double _epsilon) {
		mEpsilon = _epsilon;
	}
	
	//-----------------------------------------------------------------------------------------------------------------
	double IterativeGraspWrenchSpace::lmrw(){
		// Check if G is sujertive.
		if (mGrasp.rank() < mWrenchDims) {
			return -1.0;
		}

		// Initialization
		if (!computeInitialSet())
			return 0.0;

		// Compute initial number of facets
		uvec indices(mPolytope.n_cols, 1);
		for (unsigned i = 0; i < indices.n_rows; i++) {
			indices[i] = i;
		}

		auto combinations = generateCombinationsWithoutRepetitions(indices, 6);
		mFacets.resize(mPolytope.n_rows, combinations.size());
		for (unsigned i = 0; i < combinations.size(); i++) {	// 666 Possible efficiency hot spot
			mFacets.col(i) = combinations[i];
		}

		// Compute facets data.
		if (!computeFacetsData())
			return 0.0;


		// Get index of min dk
		unsigned minj = std::min_element(mDks.begin(), mDks.end()) - mDks.begin();

		mNk = mNks[minj];
		mDk = mDks[minj];

		// Iterate
		while (supportFunction(mNk) - mDk > mEpsilon) {
			mPolytope.insert_cols(mPolytope.n_cols, mappingFunction(mNk));
			// Compute new facets and delete old ones.
			updateFacets();
			// Compute facets data.
			if (!computeFacetsData())
				return 0.0;


			// Get index of min dk
			unsigned minj = std::min_element(mDks.begin(), mDks.end()) - mDks.begin();

			mNk = mNks[minj];
			mDk = mDks[minj];

			std::cout << supportFunction(mNk) << ", " << mDk << std::endl;
		}
		
		return mDk;
	}

	//-----------------------------------------------------------------------------------------------------------------
	//	Private interface
	//-----------------------------------------------------------------------------------------------------------------
	arma::colvec IterativeGraspWrenchSpace::orthogonalVector(arma::mat _space) const {
		mat B = null(_space.t());	// 666 TODO: suboptimal, check if needed more than one, to use directly the hole null space.
		return B.col(0);
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool IterativeGraspWrenchSpace::computeInitialSet() {
		// Compute wrenches
		mWrenches = mGrasp.aprxWrenchCone(8); 

		colvec wc = zeros(mWrenchDims);
		for (auto cp : mGrasp.contactPoints()) {
			arma::colvec wci = arma::zeros(6, 1);
			wci.rows(0, 2) = cp.normal();
			arma::vec p = cp.position();
			arma::vec n = cp.normal();
			wci.rows(3, 5) = cross(p, n);
			wc += wci;
		}
		wc /= mGrasp.contactPoints().size();


		colvec a0;
		if (sum(wc) == 0) {
			colvec u = randn(mWrenchDims);
			u /= norm(u);
			a0 = mappingFunction(u);
		}
		else {
			a0 = wc;
		}

		// Compute distance between a0 and the wrench cone Wco
		DistancePointConvexCone dist;
		dist.basePoints(mGrasp.aprxWrenchCone(8));
		dist.point(-a0);
		dist.epsilon(1e-5);

		if (!dist.compute())
			return false;

		if (dist.minDistance() > 1e-10)	// Approx 0
			return false;	// If not return 0.

						// Compute initial set V
		mPolytope = dist.baseSet();

		// Check that the initial set is at least linearly independent.
		assert(arma::rank(mPolytope) == mPolytope.n_cols);

		while (mPolytope.n_cols < 6) {
			// Compute ortogonal u'
			auto u = orthogonalVector(mPolytope);
			// Increase V with the supportFunction(u');
			mPolytope.insert_cols(mPolytope.n_cols, mappingFunction(u));
			assert(arma::rank(mPolytope) == mPolytope.n_cols);
		}
		assert(arma::rank(mPolytope) == mPolytope.n_rows);

		mPolytope.insert_cols(0, a0);
		mat sub = mPolytope.cols(0, 5);
		assert(arma::rank(sub) == mPolytope.n_rows);


		return true;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool IterativeGraspWrenchSpace::computeFacetsData() {
		// For each facet compute its data
		mNks.resize(mFacets.n_cols);
		mDks.resize(mFacets.n_cols);
		for (unsigned j = 0; j < mFacets.n_cols; j++) {
			mat subV = mPolytope.cols(mFacets.col(j));
			
			assert(subV.n_cols == 6);
			assert(arma::rank(subV) == 6);
			arma::colvec a0 = subV.col(0);
			subV.shed_col(0);
			for (unsigned l = 0; l < subV.n_cols; l++) {
				subV.col(l) -= a0;
			}

			assert(arma::rank(subV) == 5);

			auto u = orthogonalVector(subV);
			u = u / norm(u);
			
			double dotProd = dot(u, a0);
			if (dotProd < 0) {
				mNks[j] = -u;
				mDks[j] = -dotProd;
			}
			else {
				mNks[j] = u;
				mDks[j] = dotProd;
			}
		}
		return true;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool IterativeGraspWrenchSpace::updateFacets() {
		// Get indices of facets (Jk) that are contained in the polytope Vk+1 that belonged to Vk
		arma::ucolvec Jk(mFacets.n_cols);
		arma::colvec swnk = mappingFunction(mNk);
		unsigned realSize = 0;
		for (unsigned j = 0; j < mFacets.n_cols; j++) {
			if (arma::dot(mNks[j], swnk) > mDks[j]) {
				Jk[realSize] = j;
				realSize++;
			}
		}
		Jk.resize(realSize);

		// For each facet in Jk, check if they share ridges with other facet in Jk. If it is only in one
		// facet in Jk, then the ridge and sw(nk) form a new facet in Vk+1. If not, is not a facet in Vk+1.
		arma::umat newFacets;
		for (unsigned j = 0; j < Jk.n_rows; j++) {
			ucolvec vertices = mFacets.col(j);
			assert(vertices.n_elem == 6);
			for(unsigned l = 0; l < vertices.n_rows; l++) {
				// Check in all the other possible facets if someone contains the same ridge.
				ucolvec ridge = vertices;
				ridge.shed_row(l);	// 666 TODO: possible optimization, just skip the index.
				
				bool isRidgeShared = false;
				for (unsigned j2 = 0; j2 < Jk.n_rows; j2++) {
					if (j == j2)	// dont compare with the same set.
						continue;
					
					ucolvec vertices2 = mFacets.col(j2);
					unsigned numSharedIndices = 0;
					for (auto ind1 : ridge) {
						for(auto ind2: vertices2){
							if (ind1 == ind2)
								numSharedIndices++;
						}
					}
					
					if (numSharedIndices == 5) {
						isRidgeShared = true;
						break;
					}
				}

				if (!isRidgeShared) {
					ridge.insert_rows(ridge.n_rows, arma::uvec{mPolytope.n_cols-1});
					newFacets.insert_cols(newFacets.n_cols, ridge);
				}
			}
		}

		// Remove facets on Jk and add the new facets
		for (int j = Jk.n_rows - 1; j >= 0; j --) {
			mFacets.shed_col(Jk[j]);
		}

		mFacets.insert_cols(mFacets.n_cols, newFacets);

		return true;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double IterativeGraspWrenchSpace::supportFunction(const vec & _u) const {
		double maxVal = 0;
		for (unsigned i = 0; i < mWrenches.n_cols; i++) {
			double currentVal = dot(mWrenches.col(i), _u);
			if (currentVal > maxVal) {
				maxVal = currentVal;
			}
		}
		return maxVal;
	}

	//-----------------------------------------------------------------------------------------------------------------
	vec IterativeGraspWrenchSpace::mappingFunction(const vec & _u) const {
		double maxVal = 0;
		unsigned maxIndex = 0;
		for (unsigned i = 0; i < mWrenches.n_cols; i++) {
			double currentVal = dot(mWrenches.col(i), _u);
			if (currentVal > maxVal) {
				maxIndex = i;
				maxVal = currentVal;
			}
		}
		return mWrenches.col(maxIndex);
	}

}	//	namespace grasping_tools

