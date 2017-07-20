/*
 * OverSegmentation.h
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      reads in a list of parts from a file. this represents the initial oversegmented
 *      view of the world
 */

#ifndef OVERSEGMENTATION_H_
#define OVERSEGMENTATION_H_

#include <vector>

#define ARMA_NO_DEBUG
#include <armadillo>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace gpis {

	struct PartInfo
	{
		arma::mat location;
		arma::mat surfaceNormal;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr voxels = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	};

	class OverSegmentation {
	public:
		/////////////////////////////////
		// constructors and destructors
		OverSegmentation(const std::string &filename);
		OverSegmentation(const pcl::PointCloud<pcl::PointNormal> &_cloud);
		OverSegmentation();

		///////////////
		// constants
		static const int NUM_DIMS;
		const int numDims_;

		///////////////////
		// public methods
		int numParts() const { return numParts_; };
		const std::vector<PartInfo>& allParts() const { return allParts_; };
		const PartInfo& allParts(const int partIndex) const { return allParts_.at(partIndex); };

		void printInfo() const;
		bool isInitialised() const;
		arma::colvec getSceneLimits(const unsigned int colIndex) const;

		void fillWithCloud(const pcl::PointCloud<pcl::PointNormal> &_cloud);
		void fillWithCloud(const pcl::PointCloud<pcl::PointNormal> &_cloud, const std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> &_pointsPerCluster);

	private:
		///////////////////
		// private members
		int numParts_;
		std::vector<PartInfo> allParts_;
		bool initialised_;
		arma::mat sceneLimits_;

		///////////////////
		// private methods
		bool readFromFile(const std::string &filename);
	};
}	// namespace gpis

#endif /* OVERSEGMENTATION_H_ */
