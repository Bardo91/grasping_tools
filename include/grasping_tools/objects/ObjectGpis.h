///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_OBJECTGPIS_H_
#define GPISGRASPING_OBJECTGPIS_H_

#include <grasping_tools/Object.h>

#include <gpis/sampler/prior/Kernel.h>
#include <gpis/sampler/prior/Mean.h>

#include <pcl/PolygonMesh.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace grasping_tools {
	class ObjectGpis: public Object{
	public:
		/// Constructor Object GPIS 
		ObjectGpis(const arma::mat &_dataPoints, gpis::Mean *_gpMean, gpis::Kernel *_gpKernel, const arma::vec &_gpMeanParams);

		/// Constructor Object GPIS 
        ObjectGpis(const arma::mat &_dataPoints, const arma::mat &_observations, gpis::Mean *_gpMean, gpis::Kernel *_gpKernel, const arma::vec &_gpMeanParams, double _sigmaData);

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

        /// Get mesh info
        void mesh(pcl::PointCloud<pcl::PointXYZ> &_vertices, std::vector<pcl::Vertices> &_faces);

        /// Get mesh info
        void mesh(pcl::PolygonMesh &_mesh);

	protected:
		virtual void moveObject();

	private:
		void precomputeGpData();
		double closestDistanceData(const arma::vec &_point);
	private:
		arma::mat		mDataPoints;
		arma::mat		mDataNormals;
		arma::mat		mObservationPoints;
        double			mSigmaAlignment = 0.1;
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
