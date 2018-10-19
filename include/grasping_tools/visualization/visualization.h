///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_VISUALIZATION_H_
#define GPISGRASPING_VISUALIZATION_H_

#include "../ContactPoint.h"
#include <pcl/visualization/pcl_visualizer.h>

namespace grasping_tools {
	/// Simple method to plot contact points as cones in pcl visualizers.
	void plotContactPoint(const ContactPoint &_cp, pcl::visualization::PCLVisualizer &_viewer, double _scaleFactor = 1.0, std::string _tag = "cone", unsigned _viewport = 0);
}

#endif	//	GPISGRASPING_VISUALIZATION_H_