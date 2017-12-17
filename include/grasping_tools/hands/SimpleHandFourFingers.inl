///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/objects/ObjectGpis.h>
#include <grasping_tools/mathTools.h>
#include <grasping_tools/Grasp.h>
#include <cassert>

namespace grasping_tools {
	template<>
	inline Grasp SimpleHandFourFingers::generate<ObjectGpis>(ObjectGpis &_object) {
		// Get random point on a sphere with radius 1.5 of max distance from object's origin.
		arma::colvec3 center = _object.center();
		arma::mat data = _object.data();

		double minX = min(data.row(0)), maxX = max(data.row(0));
		double minY = min(data.row(1)), maxY = max(data.row(1));
		double minZ = min(data.row(2)), maxZ = max(data.row(2));
		arma::vec possibleRadius = {	abs(maxX - minX) * 1.5,
										abs(maxY - minY) * 1.5,
										abs(maxZ - minZ) * 1.5};

		double radius = max(possibleRadius);

		arma::colvec3 point = generateRandomPointSphere(center, radius);
		
		return generateGrasp<ObjectGpis>(point, _object);
	}

	//--------------------------------------------------------------------------------------------------------------------
	template<>
	inline Grasp SimpleHandFourFingers::generate<ObjectGpis>(ObjectGpis &_object, const arma::colvec3 &_initialPoint) {
		return generateGrasp<ObjectGpis>(_initialPoint, _object);
	}

	//--------------------------------------------------------------------------------------------------------------------
	template<typename ObjectType_>
	inline Grasp SimpleHandFourFingers::generateGrasp(const arma::colvec3 &_initialPoint, ObjectType_ &_object) {
		line->clear();
		// Get gradient pointing to the object.
		const double epsilon = 1e-3;

		auto data = _object.data();

		unsigned minIndex;
		double minDist = 9999999999999999.9;
		for (unsigned i = 0; i < data.n_cols; i++) {
			double dist = arma::norm(_initialPoint - data.col(i));
			if (dist < minDist) {
				minDist = dist;
				minIndex = i;
			}
		}


		arma::colvec p0 = data.col(minIndex);
		arma::colvec p1 = _initialPoint;
		do {
			arma::colvec3 newPoint = (p1 + p0) / 2;
			double val = _object.evaluate(newPoint);
			pcl::PointXYZRGB newP;
			if (val > 0) {
				p1 = newPoint;
				newP.g = 255;
			}
			else {
				p0 = newPoint;
				newP.r = 255;
			}
			newP.x = newPoint[0];
			newP.y = newPoint[1];
			newP.z = newPoint[2];
			line->push_back(newP);
		} while (norm(p1 - p0) > epsilon);


		// Get fingers' trajectories.
		unsigned nPoints = 100;
		arma::mat points1 = pointsInCircle(mRadius, p1, 0.0, 3.1416 / 4, nPoints);
		arma::mat points2 = pointsInCircle(mRadius, p1, 0.0, -3.1416 / 4, nPoints);

		// Get normal of p1.
		arma::colvec4 contactGp;
		_object.evaluate(p1, contactGp);
		arma::mat R = rotUnitaryVectors({ 0,1,0 }, contactGp.rows(1, 3));

		cp1.x = p1[0];
		cp1.y = p1[1];
		cp1.z = p1[2];
		cp1.normal_x = contactGp[1];
		cp1.normal_y = contactGp[2];
		cp1.normal_z = contactGp[3];

		// 666 Do it more efficiently;
		arma::vec eval1(nPoints), eval2(nPoints);
		for (unsigned i = 0; i < points1.n_cols; i++) {
			points1.col(i) = p1 + R*(points1.col(i) - p1);
			eval1[i] = _object.evaluate(points1.col(i));

			points2.col(i) = p1 + R*(points2.col(i) - p1);
			eval2[i] = _object.evaluate(points2.col(i));
		}

		ps1->clear();
		ps2->clear();
		for (unsigned i = 0; i < nPoints; i++) {
			pcl::PointXYZRGB p1, p2;
			p1.x = points1.col(i)[0];
			p1.y = points1.col(i)[1];
			p1.z = points1.col(i)[2];
			if (eval1[i] < 0) {
				p1.r = 255;
			}
			else {
				p1.g = 255;
			}

			p2.x = points2.col(i)[0];
			p2.y = points2.col(i)[1];
			p2.z = points2.col(i)[2];
			if (eval2[i] < 0) {
				p2.r = 255;
			}
			else {
				p2.g = 255;
			}

			ps1->push_back(p1);
			ps2->push_back(p2);

		}

		std::vector<arma::colvec3> surfacePoints1, surfacePoints2;
		// Get surface points.
		for (unsigned i = 0; i < nPoints - 1; i++) {
			if (eval1[i] < 0 != eval1[i + 1] < 0) {
				surfacePoints1.push_back(points1.col(i));
			}

			if (eval2[i] < 0 != eval2[i + 1] < 0) {
				surfacePoints2.push_back(points2.col(i));
			}
		}

		// Choose contact points
		std::vector<ContactPoint> cps;
		if (surfacePoints1.size() == 2) {
			for (auto point : surfacePoints1) {
				arma::colvec4 values;
				arma::mat44 covariances;
				_object.evaluate(point, values, covariances);

				ContactPoint cp(point, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 0.2, 0.2);
				cps.push_back(cp);
			}
		}
		else {
			assert(false);	// 666 TODO:
		}

		if (surfacePoints2.size() == 2) {
			for (auto point : surfacePoints2) {
				arma::colvec4 values;
				arma::mat44 covariances;
				_object.evaluate(point, values, covariances);

				ContactPoint cp(point, arma::eye(3, 3)*covariances(0, 0), values.tail(3), covariances.submat(1, 1, 3, 3), eContactTypes::SFC, 1, 0.2, 0.2);
				cps.push_back(cp);
			}
		}
		else {
			assert(false);	// 666 TODO:
		}
		// Generate grasp
		return 	Grasp(_object.center(), cps);
	}

}	//	gpisGrasping
