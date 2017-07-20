///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_OBJECTGPIS_H_
#define GPISGRASPING_OBJECTGPIS_H_

#include "../Object.h"

#include <Sampler/Prior/Kernel.h>
#include <Sampler/Prior/Mean.h>

namespace gpisGrasping {
	class ObjectGpis: private Object{
	public:
		/// Constructor Object GPIS 
		ObjectGpis(const arma::mat &_dataPoints, gpis::Mean *_gpMean, gpis::Kernel *_gpKernel, const arma::vec &_gpMeanParams);

		/// Constructor Object GPIS 
		ObjectGpis(const arma::mat &_dataPoints, const arma::mat &_observations, gpis::Mean *_gpMean, gpis::Kernel *_gpKernel, const arma::vec &_gpMeanParams, double &_sigmaData);

		/// Get inferred center of objects
		arma::colvec3 center(); 

		/// Get data associated to gpis object
		arma::mat data() const;

		/// Evaluate point. 
		/// \return its inferred value with the constructed data
		double evaluate(const arma::colvec3 &_point);

		/// Evaluate point. 
		/// \return its inferred value with the constructed data including gradients
		void evaluate(const arma::colvec3 &_point, arma::colvec4 &_values);

		/// Evaluate point. 
		/// \return its inferred value with the constructed data including gradients and also compute covariances.
		void evaluate(const arma::colvec3 &_point, arma::colvec4 &_values, arma::mat44 &_covariances);

	private:
		void precomputeGpData();
		double closestDistanceData(const arma::vec &_point);
	private:
		arma::mat		mDataPoints;
		arma::mat		mDataNormals;
		arma::mat		mObservationPoints;
		double			mSigmaAlignment;
		arma::vec		mCenter;

		gpis::Mean		*mGpMean;
		gpis::Kernel	*mGpKernel;
		arma::vec		mGpMeanParams;

		bool			mPrecomputedGpData = false;
		arma::mat		mCovarianceData;
		arma::uvec		mSpanVectorVal;

		arma::mat		mRVector;


		std::function<double (const arma::colvec3 & )>			mEvalFunctionOnlyVal;
		std::function<arma::colvec4(const arma::colvec3 &)>		mEvalFunction;
		std::function<arma::mat44(const arma::colvec3 &)>		mEvalCovarianceFunction;

	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_OBJECT_H_