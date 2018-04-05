///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_HANDS_GRIPPERHANDTWOFINGERS_H_
#define GPISGRASPING_HANDS_GRIPPERHANDTWOFINGERS_H_

#include <grasping_tools/Hand.h>
#include <pcl/PolygonMesh.h>
#include <pcl/point_types.h>
namespace grasping_tools {
	class GripperHandTwoFingers : public Hand {
	public:
		/// Constructor of hand. 
		/// \param _aperture: defines maximum aperture of the gripper. 
		GripperHandTwoFingers(double _aperture);

		/// Templatized method for generating grasps
		/// \param _object: templatized object to be grasped.
		/// \return: generated grasp.
		template<typename ObjectType_>
		Grasp generate(ObjectType_ &_object);

		/// Templatized method for generating a list of all possible grasps
		/// \param _object: templatized object to be grasped.
		/// \return: generated grasp.
		template<typename ObjectType_>
		std::vector<Grasp> generateGrasps(ObjectType_ &_object, double _downsampleFactor);

		/// Given a grasp, generate a 3d meshed model for visualization purposes
		/// \param _grasp: grasp to be shown
		/// \param _rotation: define the orietnation of the grasp.
		pcl::PolygonMesh generateHandModel(const Grasp &_grasp, double _rotation = 0);

	private:
		double mAperture;
		std::string mModel;
		pcl::PolygonMesh mMeshBase;
		pcl::PolygonMesh mLeftGripper;
		pcl::PolygonMesh mRightGripper;
        bool mLoadedModel = false;
    };
}

#include "GripperHandTwoFingers.inl"

#endif 
