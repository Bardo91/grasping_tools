///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_HANDS_SIMPLEHANDTWOFINGERS_H_
#define GPISGRASPING_HANDS_SIMPLEHANDTWOFINGERS_H_

#include"../Hand.h"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
namespace grasping_tools {
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