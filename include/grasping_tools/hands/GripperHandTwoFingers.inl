///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/objects/ObjectGpis.h>
#include <grasping_tools/objects/ObjectMesh.h>
#include <grasping_tools/Grasp.h>
#include <grasping_tools/mathTools.h>
#include <deque>

namespace grasping_tools {
	template<>
	inline Grasp GripperHandTwoFingers::generate<ObjectGpis>(ObjectGpis &_object) {
		// Get random point on a sphere with radius 1.5 of max distance from object's origin.
		arma::colvec3 center = _object.center();
		arma::mat data = _object.data();

        double minX = arma::min(data.row(0)), maxX = arma::max(data.row(0));
        double minY = arma::min(data.row(1)), maxY = arma::max(data.row(1));
        double minZ = arma::min(data.row(2)), maxZ = arma::max(data.row(2));

        double radius = std::max(std::max(fabs(maxX - minX) * 3,
                            fabs(maxY - minY) * 3),
                            fabs(maxZ - minZ) * 3);


		int randIdx = double(rand())/RAND_MAX*(data.n_cols -1);
		arma::colvec cpPos1 = data.col(randIdx).head(3);;
		// Get opposite point
		arma::colvec p0 = cpPos1 - data.col(randIdx).tail(3)*0.01;
		arma::colvec p1 = cpPos1 - data.col(randIdx).tail(3)*0.5;
		do {
			arma::colvec3 newPoint = (p1 + p0) / 2;
			double val = _object.evaluate(newPoint);
			if (val > 0) {
				p1 = newPoint;
			}
			else {
				p0 = newPoint;
			}
		} while (norm(p1 - p0) > 1e-3);
		arma::colvec cpPos2 = p0;

		if(arma::norm(cpPos1.head(3)-cpPos2.head(3)) > mAperture){
			return Grasp();
		}

		std::vector<ContactPoint> cps;
		arma::colvec4 values;
		arma::mat44 covariances;
		_object.evaluate(cpPos1, values, covariances);
		// Normalize normal from GPIS
		values.tail(3) /= norm(values.tail(3));
        cps.push_back(ContactPoint(cpPos1, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 1, 1));

		_object.evaluate(cpPos2, values, covariances);
		// Normalize normal from GPIS
		values.tail(3) /= norm(values.tail(3));
        cps.push_back(ContactPoint(cpPos2, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 1, 1));

		Grasp grasp;
		grasp.contactPoints(cps);

		return grasp;

	}

	//-----------------------------------------------------------------------------------------------------------------
	template<>
	inline Grasp GripperHandTwoFingers::generate<ObjectMesh>(ObjectMesh &_object) {
        // Generate grasp
		Grasp grasp;
        std::vector<ContactPoint> cps;

        auto candidatePoints = _object.centroidFaces();

        // Get a random candidate point // 666 TODO improve candidate selection
        arma::imat id = arma::randi(1, 1, arma::distr_param(0, candidatePoints.n_cols - 1));
        arma::colvec6 candidatePoint = candidatePoints.col(id(0,0));

        arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2)+ 2*candidatePoint.rows(3, 5), candidatePoint.rows(0, 2) - candidatePoint.rows(3, 5));
		//std::cout << intersections << std::endl;
		if(intersections.n_cols >= 2){
			if(arma::norm(intersections.col(0).head(3)-intersections.col(1).head(3)) > mAperture){
				return Grasp();
			}

			for(unsigned int i = 0; i < intersections.n_cols; i++){
				ContactPoint cp(intersections.col(i).head(3), arma::eye(3, 3), intersections.col(i).tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1);
				cps.push_back(cp);
			}

			grasp.contactPoints(cps);

			return grasp;
		}else{
			return Grasp();
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	template<>
	inline  std::vector<Grasp> GripperHandTwoFingers::generateGrasps(ObjectMesh &_object, double _downsampleFactor){
		// Generate grasp
		std::vector<Grasp> grasps;

        auto candidatePoints = _object.centroidFaces();

		unsigned step = candidatePoints.n_cols*_downsampleFactor;
		step = step<1? 1: step;

		for(unsigned i = 0; i < candidatePoints.n_cols; i += step) {
			arma::colvec6 candidatePoint = candidatePoints.col(i);
			arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2)+ 2*candidatePoint.rows(3, 5), candidatePoint.rows(0, 2) - candidatePoint.rows(3, 5));
			//std::cout << intersections << std::endl;
			if(intersections.n_cols >= 2){
				if(arma::norm(intersections.col(0).head(3)-intersections.col(1).head(3)) < mAperture){
        			std::vector<ContactPoint> cps;
					for(unsigned int i = 0; i < intersections.n_cols; i++){
						ContactPoint cp(intersections.col(i).head(3), arma::eye(3, 3), intersections.col(i).tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1);
						cps.push_back(cp);
					}
					Grasp grasp;
					grasp.contactPoints(cps);
					grasps.push_back(grasp);
				}
			}
		}

		return grasps;
	}

}
