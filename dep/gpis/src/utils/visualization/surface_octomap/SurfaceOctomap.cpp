///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <gpis/utils/visualization/surface_octomap/SurfaceOctomap.h>
#include <gpis/utils/visualization/surface_octomap/expandCell.h>
#include <gpis/utils/visualization/surface_octomap/validateCell.h>

using namespace arma;

namespace gpis {
	//--------------------------------------------------------------------------------------------------------------------
	SurfaceOctomap::SurfaceOctomap(Mean * _mean, Kernel * _kernel): mMean(_mean), mKernel(_kernel) {
	}

	//--------------------------------------------------------------------------------------------------------------------
	mat SurfaceOctomap::data() const {
		return mData;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceOctomap::data(const mat & _data) {
		mComputed = false;
		mData = _data;
	}

	//--------------------------------------------------------------------------------------------------------------------
	bool SurfaceOctomap::compute(const arma::vec &_params, unsigned _preExpandingIterations, unsigned _iterations) {
		if (mData.n_cols == 0)
			return false;

		// Compute covariance of datapoints.
		mat K = (*mKernel)(mData.rows(0, 2), true);

		// Compute mean of data points.
		mat mu = (*mMean)(mData.rows(0, 2), _params, true);
		mu.reshape(mu.n_cols*mu.n_rows, 1);

		// Compute deviation
		vec f(mData.n_cols*4);

		for (unsigned i = 0; i < mData.n_cols; i++) {
			vec observedVals(4);
			observedVals.head(1) = 0;
			observedVals.tail(3) = mData.col(i).rows(3,5);
			f.subvec(i*4, i*4 + 3) = observedVals;
		}
		vec deviation = f - mu;

		mat qMat = solve(K, deviation);

		GpisCell::evalFun evalFunction = [&](const vec3 &_point)->double {
			return ((*mMean)(_point, _params, true) + (*mKernel)(_point, mData.rows(0, 2), true) * qMat).at(0,0);
		};

		// Init rootcell
		
		mRoot = new GpisCell(	{ mMinX , mMaxX },
								{ mMinY , mMaxY },
								{ mMinZ , mMaxZ },
								evalFunction);

		for (unsigned i = 0; i < _preExpandingIterations; i++) {
			expandCell(*mRoot, evalFunction);
		}
		validateCell(*mRoot, *mRoot);
		
		for (unsigned i = 0; i < _iterations; i++) {
			expandCell(*mRoot, evalFunction);
			validateCell(*mRoot, *mRoot);
		}

		mComputed = true;
		return mComputed;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceOctomap::surfacePoints(mat &_cloud) const {
		assert(false); // Unimplemented.
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceOctomap::limits(double _minX, double _maxX, double _minY, double _maxY, double _minZ, double _maxZ) {
		mMinX = _minX;
		mMaxX = _maxX;
		mMinY = _minY;
		mMaxY = _maxY;
		mMinZ = _minZ;
		mMaxZ = _maxZ;
	}

	//--------------------------------------------------------------------------------------------------------------------
}
