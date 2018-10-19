///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_HANDS_ARM3DOF_H_
#define GPISGRASPING_HANDS_ARM3DOF_H_

#include<grasping_tools/Hand.h>
#include <pcl/PolygonMesh.h>

namespace grasping_tools {
	/// Class that inherit from Hand parent class that is used to generate grasps taking into account the restrictions
	/// of a 3DoF arm.
	///
	/// Class defined in `#include <grasping_tools/hands/Arm3DOF.h>`
	class Arm3DOF : public Hand {
	public:
		/// Constructor of hand. 
		/// \param _aperture: defines maximum aperture of the gripper. 
		Arm3DOF(double _aperture);

		/// Templatized method for generating grasps
		/// \param _object: templatized object to be grasped.
		/// \param _pointIK: return final point from IK.
		/// \return: generated grasp.
		template<typename ObjectType_>
		Grasp generate(ObjectType_ &_object, arma::colvec3 &_pointIK);

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

	public:
		arma::mat puntosCirculo;
		std::vector<bool> validPointsIk;
    };
}

#include "Arm3DOF.inl"

#endif 
