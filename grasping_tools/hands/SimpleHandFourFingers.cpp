///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SimpleHandFourFingers.h"

namespace grasping_tools {
	SimpleHandFourFingers::SimpleHandFourFingers(double _radius):ps1(new pcl::PointCloud<pcl::PointXYZRGB>) ,ps2(new pcl::PointCloud<pcl::PointXYZRGB>), line(new pcl::PointCloud<pcl::PointXYZRGB>) {
		mRadius = _radius;
	}
}
