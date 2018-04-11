///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_OBJECT_H_
#define GPISGRASPING_OBJECT_H_

#include <armadillo>
#include <Eigen/Eigen>

namespace grasping_tools {
	class Object {
	public:
		/// Return the centroid of the object.
		virtual arma::colvec3 center() = 0;

		/// Get the min and max point that encloses boundaries of the object.
		virtual void minMax(arma::colvec3 &_min, arma::colvec&_max) {std::cout << "Not implemented" << std::endl;};

		/// Set a pose to an object.
		void move(arma::mat44 &_pose){
			mPose = _pose*mPose;
			for(auto i = 0; i < 4; i++){
				for(auto j=0; j <4; j++){
					mPoseEigen(i,j) = mPose(i,j);
				}
			}
			moveObject();
		}

		/// Get object pose
		arma::mat44 move(){
			return mPose;
		}

	protected:
		virtual void moveObject() = 0;

	protected:
		arma::mat44 	mPose = arma::eye(4,4);
		Eigen::Matrix4f mPoseEigen;		// Computed once 
	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_OBJECT_H_