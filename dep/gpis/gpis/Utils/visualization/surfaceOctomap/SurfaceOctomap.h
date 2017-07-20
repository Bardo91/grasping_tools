///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_SURFACEOCTOMAP_H_
#define GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_SURFACEOCTOMAP_H_

#include "GpisCell.h"
#include "getLeafsCell.h"

#include <Sampler/Prior/Mean.h>
#include <Sampler/Prior/Kernel.h>

#include <armadillo>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace gpis {
	class SurfaceOctomap {
	public:
		/// Constructor.
		/// \param _mean:
		/// \param _kernel:
		SurfaceOctomap(Mean *_mean, Kernel *_kernel);

		/// Get current observations associated to the GPIS
		arma::mat data() const;

		/// Set associated data to the GPIS.
		void data(const arma::mat &_data);

		/// Compute surface.
		/// \param _preExpandingIterations:
		/// \param _iterations:
		bool compute(const arma::vec &_params, unsigned _preExpandingIterations, unsigned _iterations);

		/// Get last computed surface.
		void surfacePoints(arma::mat &_cloud) const;

		/// Get last computed surface.
		template<typename PointType_>
		void surfacePoints(pcl::PointCloud<PointType_> &_cloud) const;

		/// Set limits for the computation
		void limits(double _minX, double _maxX, double _minY, double _maxY, double _minZ, double _maxZ);

	private:
		GpisCell *mRoot;

		Mean *mMean;
		Kernel *mKernel;

		arma::mat mData;

		bool mComputed = false;
		double mMinX, mMaxX, mMinY, mMaxY, mMinZ, mMaxZ;

	};

	// Inline definitions
	template<typename PointType_>
	void SurfaceOctomap::surfacePoints(pcl::PointCloud<PointType_> &_cloud) const {
		getLeafsCell(*mRoot, _cloud);
	}
}


#endif	//	GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_SURFACEOCTOMAP_H_
