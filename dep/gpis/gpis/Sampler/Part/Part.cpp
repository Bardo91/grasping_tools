/*
 * Part.cpp
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      represents a part of a scene. a set of parts makes an object.
 */

#include "Part.h"


namespace gpis {
	///////////////////////////////////////////////////////////////////////////////////
	// constructors and destructors
	//
    Part::Part(	unsigned index,
				const arma::vec &location,
				const arma::vec &surfaceNormals,
                const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &voxels)
        :	location_(location),
            surfaceNormals_(surfaceNormals),
            voxels_(voxels),
            index_(index),
            cholUpdateData_(location_.n_elem)
	{
        //empty
	}
}	//	namespace gpis
