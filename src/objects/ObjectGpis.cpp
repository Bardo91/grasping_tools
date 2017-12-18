///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/objects/ObjectGpis.h>
#include <gpis/Utils/SurfaceGpis.h>

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectGpis::ObjectGpis(const arma::mat & _dataPoints, gpis::Mean * _gpMean, gpis::Kernel * _gpKernel, const arma::vec &_gpMeanParams):
	mDataPoints(_dataPoints.rows(0,2)),
	mDataNormals(_dataPoints.rows(3, 5)),
	mGpMean(_gpMean), 
	mGpKernel(_gpKernel),
	mGpMeanParams(_gpMeanParams){

}

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectGpis::ObjectGpis(const arma::mat &_dataPoints, const arma::mat &_observations, gpis::Mean *_gpMean, gpis::Kernel *_gpKernel, const arma::vec &_gpMeanParams, double _sigmaData):
	mDataPoints(_dataPoints.rows(0, 2)),
	mObservationPoints(_observations),
	mSigmaAlignment(_sigmaData),
	mDataNormals(_dataPoints.rows(3, 5)),
	mGpMean(_gpMean),
	mGpKernel(_gpKernel),
	mGpMeanParams(_gpMeanParams) {

}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 grasping_tools::ObjectGpis::center(){
	if (mCenter.n_elem == 0) {
		mCenter.resize(3);
		mCenter[0] = arma::accu(mDataPoints.row(0)) / mDataPoints.n_cols;
		mCenter[1] = arma::accu(mDataPoints.row(1)) / mDataPoints.n_cols;
		mCenter[2] = arma::accu(mDataPoints.row(2)) / mDataPoints.n_cols;
	}

	return mCenter;
}

//---------------------------------------------------------------------------------------------------------------------
double grasping_tools::ObjectGpis::evaluate(const arma::vec3 & _point){
	if (!mPrecomputedGpData) {
		precomputeGpData();
		mPrecomputedGpData = true;
	}

	return mEvalFunctionOnlyVal(_point);
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectGpis::evaluate(const arma::vec3 & _point, arma::colvec4 & _values) {
	if (!mPrecomputedGpData) {
		precomputeGpData();
		mPrecomputedGpData = true;
	}

	_values = mEvalFunction(_point);
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectGpis::evaluate(const arma::vec3 & _point, arma::colvec4 & _values, arma::mat44 & _covariances) {
	if (!mPrecomputedGpData) {
		precomputeGpData();
		mPrecomputedGpData = true;
	}

	_values = mEvalFunction(_point);
	_covariances = mEvalCovarianceFunction(_point);

}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectGpis::precomputeGpData() {
	mCovarianceData = (*mGpKernel)(mDataPoints, true);

	// Check if model is aligned, i.e., if has observations and stdDev and add the varying noise.
	if (mObservationPoints.n_cols != 0) {
		for (unsigned i = 0; i < mDataPoints.n_cols; i++) {
			double dist = closestDistanceData(arma::vec(mDataPoints.col(i)));
			mCovarianceData(i, i) += mSigmaAlignment*exp(dist);
		}
	}
	//

	mSpanVectorVal.resize(mDataPoints.n_cols);

	for (unsigned i = 0; i < mDataPoints.n_cols; i++) {
		mSpanVectorVal[i] = i * 4;
	}

	// Compute F  and R for the data
	arma::mat observedVals = mDataNormals;
	observedVals.insert_rows(0, 1, true);
	arma::mat fPlusData = observedVals - (*mGpMean)(mDataPoints, mGpMeanParams, true);
	fPlusData.reshape(fPlusData.n_elem, 1);
	mRVector = arma::solve(mCovarianceData, fPlusData);

	
	mEvalFunctionOnlyVal = [&](const arma::colvec3 &_point)->double {
		auto crossCovariance = (*mGpKernel)(_point, mDataPoints, true);
		return ((*mGpMean)(_point, mGpMeanParams, true) + crossCovariance * mRVector)[0];
	};

	mEvalFunction = [&](const arma::colvec3 &_point)->arma::vec {
		auto crossCovariance = (*mGpKernel)(_point, mDataPoints, true);
		return (*mGpMean)(_point, mGpMeanParams, true) + crossCovariance * mRVector;
	};

	// Check if model is aligned, i.e., if has observations and stdDev and add the varying noise.
	if (mObservationPoints.n_cols != 0) {
		mEvalCovarianceFunction = [&](const arma::colvec3 &_point)->arma::mat {
			auto crossCovariance = (*mGpKernel)(_point, mDataPoints, true);
			arma::mat noiseSigma = arma::eye(4, 4)*mSigmaAlignment*exp(closestDistanceData(arma::vec({ _point[0], _point[1],_point[2] }))); // 666 TODO include normals!!!
			arma::mat pointCovariance = (*mGpKernel)(_point, true) + noiseSigma;
			return pointCovariance - crossCovariance*arma::solve(mCovarianceData, crossCovariance.t());
		};
	}
	else {
		mEvalCovarianceFunction = [&](const arma::colvec3 &_point)->arma::mat {
			auto crossCovariance = (*mGpKernel)(_point, mDataPoints, true);
			auto pointCovariance = (*mGpKernel)(_point, true);
			return pointCovariance - crossCovariance*arma::solve(mCovarianceData, crossCovariance.t());
		};
	}
	//
	
}


double grasping_tools::ObjectGpis::closestDistanceData (const arma::vec &_point) {	// 666 TODO do it with kdtree!
	double minDist = 99999;
	int minIdx = 0;
	for (unsigned i = 0; i < mObservationPoints.n_cols; i++) {
		double dist;
		if ((dist = sqrt(arma::dot(_point, mObservationPoints.col(i)))) < minDist) {
			minDist = dist;
			minIdx = 0;
		}
	}
	return minDist;
}

//---------------------------------------------------------------------------------------------------------------------
arma::mat grasping_tools::ObjectGpis::data() const {
	return mDataPoints;
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectGpis::mesh(pcl::PointCloud<pcl::PointXYZ> &_vertices, std::vector<pcl::Vertices> &_faces){

}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectGpis::mesh(pcl::PolygonMesh &_mesh){

}
