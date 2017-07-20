///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "DistancePointConvexCone.h"
#include "mathTools.h"

using namespace arma;

namespace gpisGrasping {
	//--------------------------------------------------------------------------------------------------------------------
	void DistancePointConvexCone::basePoints(const mat & _A) {
		mA = _A;
	}
	
	//--------------------------------------------------------------------------------------------------------------------
	void DistancePointConvexCone::point(const colvec & _b) {
		mB = _b;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void DistancePointConvexCone::epsilon(double _epsilon) {
		mEpsilon = _epsilon;
	}

	//--------------------------------------------------------------------------------------------------------------------
	double DistancePointConvexCone::minDistance() const {
		return mMinDistance;
	}

	//--------------------------------------------------------------------------------------------------------------------
	colvec DistancePointConvexCone::closestPoint() const{
		return mClosestPoint;
	}

	arma::mat DistancePointConvexCone::baseSet() const {
		return mAi;	// or mAmi??
	}

	//--------------------------------------------------------------------------------------------------------------------
	bool DistancePointConvexCone::compute() {
		// Init variables
		unsigned dim = mB.size();
		mV = zeros(dim);
		mR = mB;
		mAi = arma::mat();
		mAmi = arma::mat();

		// Iterate
		unsigned iters;
		while (supportFunction(mR) > mEpsilon) {
			if (iters++ > mMaxIters)
				return false;

			auto sA = mappingFunction(mR);
			mAi = mAmi;
			mAi.insert_cols(mAi.n_cols, sA);
			// Subalgorithm.
			updateSubConvexCone();

			mR = mB - mV;
		}


		// Set final values;
		mClosestPoint = mV;
		mMinDistance = norm(mR);
		return true;
	}

	//--------------------------------------------------------------------------------------------------------------------
	bool DistancePointConvexCone::checkInputData() {
		if (mA.n_cols < 1 || mB.size() < 1 || mA.n_cols != mB.size())
			return false;
		else
			return true;
	}

	//--------------------------------------------------------------------------------------------------------------------
	double DistancePointConvexCone::supportFunction(const arma::colvec & _r) {
		double maxVal = 0;

		for (unsigned i = 0; i < mA.n_cols; i++) {
			double val = dot(_r, mA.col(i));
			if (val > maxVal)
				maxVal = val;
		}
		
		return maxVal;
	}
	
	//--------------------------------------------------------------------------------------------------------------------
	arma::colvec DistancePointConvexCone::mappingFunction(const arma::colvec & _r) {
		double maxVal = 0;
		unsigned maxIndex = 0;
		for (unsigned i = 0; i < mA.n_cols; i++) {
			double val = dot(_r, mA.col(i));
			if (val > maxVal) {
				maxVal = val;
				maxIndex = i;
			}
		}

		return mA.col(maxIndex);
	}

	//--------------------------------------------------------------------------------------------------------------------
	void DistancePointConvexCone::updateSubConvexCone() {
		for (;;) {
			// ----- STEP 1 -----
			// Compute the associated coefficient vector of the projection of b on the linear space spanned by Ai.
			arma::colvec c = solve(mAi.t()*mAi, mAi.t()*mB);
			
			// Create submatrices of A that correspong to negative, cero and positive indices of c;
			arma::uvec posIndices, ceroIndices, negIndices;
			for (unsigned i = 0; i < c.size(); i++) {
				if (c(i) < 0) {
					arma::uvec index = { i };
					negIndices.insert_rows(negIndices.n_rows, index);
				}
				else if (c(i) == 0) {
					arma::uvec index = { i };
					ceroIndices.insert_rows(ceroIndices.n_rows, index);
				}
				else {
					arma::uvec index = { i };
					posIndices.insert_rows(posIndices.n_rows, index);
				}
			}

			if (negIndices.size() == 0) { // If all the coefficients are positives or cero
				// Compute projection of b on the linear space spanned by Ai.
				arma::colvec p = mAi*c;
				mV = p;
				mAmi = mAi.cols(posIndices);
				return;
			}
			// ----- STEP 2 -----
			else if (negIndices.size() == 1) { // It the subspace of the negative coefficients contains only one point.
				// Remove that point from the matrix and repeat the process.
				mAi.shed_col(negIndices(0));
				continue;
			}
			else {
				// ----- STEP 3 -----
				// Remove sA as possible option. Is always the last cols, so go from c = 1 ... N-1 (or c = 0 ... N - 2).
				for (unsigned i = 0; i < negIndices.size(); i++) {
					if (negIndices(i) == mAi.n_cols - 1) // sA needs to be in the system, so ignore it.
						continue;

					arma::mat Aip = mAi;
					Aip.shed_col(negIndices(i));
					// ----- STEP 4 -----
					// If subset satisfied conditions of theorem 4:
					//		1) All the components of cp are strictly positives.
					//		2) a'*rp(b) < 0 for all a in cone(mAi) (the article also says: "mAi\Aip" ??).
					arma::colvec cp = solve(Aip.t()*Aip, Aip.t()*mB);
					arma::colvec pp = Aip*cp;
					arma::colvec rp = mB - pp;

					// Checking first condition.
					for (unsigned j = 0; j < cp.size(); j++) {
						if (cp[j] <= 0) // if any element of cp is less or equal to cero then cp is not strictly positive.
							continue;
					}

					// Checking second condition
					for (unsigned j = 0; j < mAi.n_cols; j++) {
						if (arma::dot(mAi.col(j), rp) > 0) // If second condition is not satisfied for some of the elements, then continue with other posible cones.
							continue;
					}

					// Finally, if all conditions are satisfied.
					mV = pp;
					mAmi = Aip;
					return;
				}
			}
		}

	}
}