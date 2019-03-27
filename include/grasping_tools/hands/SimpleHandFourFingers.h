///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_HANDS_SIMPLEHANDTWOFINGERS_H_
#define GPISGRASPING_HANDS_SIMPLEHANDTWOFINGERS_H_

#include <grasping_tools/Hand.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
namespace grasping_tools {
	/// Class that inherit from Hand parent class that is able to generate grasps for objects using a four fingers
	/// gripper model in cross morphology.
	///
	/// Class defined in `#include <grasping_tools/hands/GripperHandTwoFingers.h>`
	class SimpleHandFourFingers: public Hand {
	public:
		/// Constructor
		/// \param _radius: radius used for computing contact points
		SimpleHandFourFingers(double _radius);

		/// Generate a random grasp for the given object.
		template<typename ObjectType_>
		Grasp generate(ObjectType_ &_object);

		/// Generate a grasp given an initial approximation point
		template<typename ObjectType_>
		Grasp generate(ObjectType_ &_object, const arma::colvec3 &_initialPoint);


		pcl::PointCloud<pcl::PointXYZRGB>::Ptr ps1, ps2, line;
		pcl::PointNormal cp1;

	private:
		template<typename ObjectType_>
		Grasp generateGrasp(const arma::colvec3 &_initialPoint, ObjectType_ &_object);

	private:
		double mRadius;
	};


}	//	gpisGrasping

#include "SimpleHandFourFingers.inl"

#endif	//	GPISGRASPING_HANDS_SIMPLEHANDTWOFINGERS_H_