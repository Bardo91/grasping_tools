///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <objects/ObjectGpis.h>
#include <objects/ObjectMesh.h>
#include <Grasp.h>
#include <mathTools.h>

namespace gpisGrasping {
	template<>
	inline Grasp GripperHandTwoFingers::generate<ObjectGpis>(ObjectGpis &_object) {
		// Get random point on a sphere with radius 1.5 of max distance from object's origin.
		arma::colvec3 center = _object.center();
		arma::mat data = _object.data();

		double minX = min(data.row(0)), maxX = max(data.row(0));
		double minY = min(data.row(1)), maxY = max(data.row(1));
		double minZ = min(data.row(2)), maxZ = max(data.row(2));
		arma::vec possibleRadius = { abs(maxX - minX) * 1.5,
			abs(maxY - minY) * 1.5,
			abs(maxZ - minZ) * 1.5 };

		double radius = max(possibleRadius);

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
		p1 = _object.center() - (initPoint-_object.center());
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
		cps.push_back(ContactPoint(cpPos1, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 0.4, 0.4));

		_object.evaluate(cpPos2, values, covariances);
		//std::cout << "normal 2: " << values.tail(3).t() << ". Norm: " << norm(values.tail(3)) << std::endl;
		//std::cout << "pos2: " << cpPos2.t();
		// Normalize normal from GPIS
		values.tail(3) /= norm(values.tail(3));
		cps.push_back(ContactPoint(cpPos2, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 0.4, 0.4));

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

		auto candidatePoints = _object.centroidFaces();

		// Get a random candidate point // 666 TODO improve candidate selection
		arma::imat id = arma::randi(1, 1, arma::distr_param(0, candidatePoints.n_cols - 1));

		arma::colvec6 candatePoint = candidatePoints.col(id(0,0));

		arma::mat intersections = _object.intersectRay(candidatePoints.rows(0, 2), candidatePoints.rows(0, 2) + candidatePoints.rows(3, 5));

		// Analisys of cut points (finger size, aperture etc...)

		Grasp grasp;
		return grasp;
	}

}