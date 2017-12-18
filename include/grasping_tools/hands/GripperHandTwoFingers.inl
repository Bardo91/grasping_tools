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

		arma::colvec3 initPoint = generateRandomPointSphere(center, radius);
		// Get point on surface - line thought center
		arma::colvec p0 = _object.center();
		arma::colvec p1 = initPoint;
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
		arma::colvec cpPos1 = p0;


		// Get opposite point
		p0 = _object.center();
        p1 = _object.center() - 2*(cpPos1-_object.center());
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

		std::vector<ContactPoint> cps;
		arma::colvec4 values;
		arma::mat44 covariances;
		_object.evaluate(cpPos1, values, covariances);
		//std::cout << "normal 1: " << values.tail(3).t() << ". Norm: " << norm(values.tail(3)) << std::endl;
		//std::cout << "pos1: " << cpPos1.t();
		// Normalize normal from GPIS
		values.tail(3) /= norm(values.tail(3));
        cps.push_back(ContactPoint(cpPos1, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 1, 1));

		_object.evaluate(cpPos2, values, covariances);
		//std::cout << "normal 2: " << values.tail(3).t() << ". Norm: " << norm(values.tail(3)) << std::endl;
		//std::cout << "pos2: " << cpPos2.t();
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
		// Is it possible to check if point is inside of the mesh?
		// Section analisys? slice with random planes and use that as simplified grasp version
		// Somekind of shape analisys? 
		// Pregeneration of contact points and choose grasp with them
		// Machine learning for grasping?
/*
		auto candidatePoints = _object.centroidFaces();

		// Get a random candidate point // 666 TODO improve candidate selection
		arma::imat id = arma::randi(1, 1, arma::distr_param(0, candidatePoints.n_cols - 1));

        arma::colvec6 candidatePoint = candidatePoints.col(id(0,0));

        arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2), candidatePoint.rows(0, 2) + candidatePoint.rows(3, 5));

        typedef std::pair<int, double> PairRayData;
        std::vector<PairRayData> arrangedPoints;
        for(unsigned int i = 0; i < intersections.n_cols; i++){
            arma::colvec3 v = candidatePoint.rows(0,2) - intersections.col(i);
            if(arma::dot(v, candidatePoint.rows(3,5)) <0){
                arrangedPoints.push_back(PairRayData(i, -arma::norm(v)));
            }else{
                arrangedPoints.push_back(PairRayData(i, arma::norm(v)));
            }
        }

        std::sort(arrangedPoints.begin(), arrangedPoints.end(),[](const PairRayData &_a, const PairRayData &_b) -> bool{
            return _a.second < _b.second;
        });

		// Analisys of cut points (finger size, aperture etc...)
        // 666 TODO
*/
        // Generate grasp
		Grasp grasp;
        std::vector<ContactPoint> cps;

        auto candidatePoints = _object.centroidFaces();


        // Get a random candidate point // 666 TODO improve candidate selection
        arma::imat id = arma::randi(1, 1, arma::distr_param(0, candidatePoints.n_cols - 1));
        arma::colvec6 candidatePoint = candidatePoints.col(id(0,0));

        arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2)+ 2*candidatePoint.rows(3, 5), candidatePoint.rows(0, 2) - candidatePoint.rows(3, 5));
        //typedef std::pair<int, double> PairRayData;
        //std::vector<PairRayData> arrangedPoints;
        for(unsigned int i = 0; i < intersections.n_cols; i++){
            //arma::colvec3 v = candidatePoint.rows(0,2) - intersections.col(i);
            //if(arma::dot(v, candidatePoint.rows(3,5)) <0){
            //    arrangedPoints.push_back(PairRayData(i, -arma::norm(v)));
            //}else{
            //    arrangedPoints.push_back(PairRayData(i, arma::norm(v)));
            //}
            ContactPoint cp(intersections.col(i).head(3), arma::eye(3, 3), intersections.col(i).tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1);
            cps.push_back(cp);
        }

        //ContactPoint cp1(cpPos2, arma::eye(3, 3), values.tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 0.4, 0.4);
        //ContactPoint cp2(cpPos2, arma::eye(3, 3), values.tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 0.4, 0.4);
        //cps.push_back(cp1);
        //cps.push_back(cp2);
        grasp.contactPoints(cps);

		return grasp;
	}

}
