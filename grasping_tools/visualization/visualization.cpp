///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "visualization.h"
#include <pcl/ModelCoefficients.h>

namespace grasping_tools {
	void plotContactPoint(const ContactPoint & _cp, pcl::visualization::PCLVisualizer & _viewer, double _scaleFactor, std::string _tag, unsigned _viewport) {
		pcl::ModelCoefficients cone;
		cone.values.resize(7);
		cone.values[0] = _cp.position()[0];
		cone.values[1] = _cp.position()[1];
		cone.values[2] = _cp.position()[2];
		cone.values[3] = _cp.normal()[0]*_scaleFactor;
		cone.values[4] = _cp.normal()[1]*_scaleFactor;
		cone.values[5] = _cp.normal()[2]*_scaleFactor;
		cone.values[6] = 15;

		if (_viewer.contains(_tag)) {
			_viewer.removeShape(_tag);
		}
		
		_viewer.addCone(cone, _tag, _viewport);
		_viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, _tag, _viewport);

	}
	
}
