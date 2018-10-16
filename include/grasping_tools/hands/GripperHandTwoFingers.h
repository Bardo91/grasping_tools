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
		std::vector<Grasp> generateGrasps(ObjectType_ &_object, double _resolution);

		/// Given a grasp, generate a 3d meshed model for visualization purposes
		/// \param _grasp: grasp to be shown
		/// \param _rotation: define the orietnation of the grasp.
		///
		///	The method performs ray. 
		///	Let be a ray tracing the first object. It generates 4 intersections called a, b, c and d.
		/// The number of intersections should be always even.
		/// From the figure it is intuitive that 3 "compressing" grasps are possible (a,b), (c,d) and (a,d).
		/// If intersections are arranged in distance from the source. 
		///
		///					_________		___________
		///				 a	|		| b	  c	|		   | d
		///	>---------	 <--|		|--> <--|		   |-->
		///					|		|		|		   |
		///					|		|_______|		   |
		///					|				 		   |
		///					|__________________________|
		///
		/// In this other example we have (a,b), (c,d), (e,f), (a,d), (c,f), (a,f)
		///					_________		___________				  _______________
		///				 a	|		| b	  c	|		   | d			e |				| f
		///	>---------	 <--|		|--> <--|		   |-->		   <--|				|-->
		///					|		|		|		   |			  |				|
		///					|		|_______|		   |______________|				|
		///					|				 		   								|
		///					|_______________________________________________________|
		///
		/// In general terms, 	N_intersections = N_folds*2
		///						N_grasps = SUM(i=1 to N_folds) { i }
		///
		///
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
