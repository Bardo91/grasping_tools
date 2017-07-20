/*
 * Part.h
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      represents a part of a scene. a set of parts makes an object.
 */

#ifndef PART_H_
#define PART_H_

#define ARMA_NO_DEBUG

#include "CholUpdateData.h"

#include <pcl/visualization/pcl_visualizer.h>
#include <armadillo>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace gpis {
	class Part {
	public:
		////////////////////////////////
		// constructors and destructors
        Part(   unsigned index,
                const arma::vec &location,
                const arma::vec &surfaceNormals,
                const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &voxels);

		///////////////////
        // public methods
        inline const arma::vec& location() const { return location_; }
        inline const double& location(const int d) const { return location_(d); }
        inline const arma::vec& surfaceNormals() const { return surfaceNormals_; }

        inline pcl::PointXYZRGB &voxelPoint(const int index) const { return (*voxels_)[index]; }
        inline pcl::PointCloud<pcl::PointXYZRGB> voxelPointCloud() const { return *voxels_; }
        inline int numVoxels() const { return voxels_->size(); }
        inline double surfaceNormals(const int index) const { return surfaceNormals_[index]; }
        inline void addPoint(const pcl::PointXYZRGB &point) { voxels_->push_back(point); }
        inline unsigned index() const {return index_;}

	private:
		///////////////////
		// private members
        arma::vec location_;
        arma::vec surfaceNormals_;
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr voxels_;
        unsigned index_;
	public:
		///////////////////
		// public members
		CholUpdateData cholUpdateData_;
	};
}	//	namespace gpis
#endif /* PART_H_ */
