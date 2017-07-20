///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "mathTools.h"
#include <Eigen/Eigen>

using namespace arma;

namespace gpisGrasping {
	//--------------------------------------------------------------------------------------------------------------------
	mat samplePointsEllipsoidEquidistributed(double _a, double _b, double _c, unsigned _nPoints) {
		mat points(3, _nPoints+1);
		// Generate points on a Sphere of radius 1 and scale them using the semiaxis.
		double R = 1;
		
		double a = 4 * M_PI * R * R / _nPoints;
		double d = sqrt(a);

		unsigned mt = unsigned(round(M_PI / d));
		double dt = M_PI / mt;
		double dp = a / dt;

		unsigned counter = 0;
		for (unsigned m = 0; m < mt; m++) {
			double t = M_PI*(m + 0.5) / mt;
			double mp = round(2 * M_PI*sin(t) / dp);
			for (unsigned n = 0; n < mp; n++) {
				double p = 2 * M_PI*n / mp;
				colvec3 newPoint = {_a*sin(t)*cos(p), _b*sin(t)*sin(p), _c*cos(t)};
				points.col(counter) = newPoint;
				counter++;
			}
		}

		points.resize(points.n_rows, _nPoints);

		return points;
	}

	//--------------------------------------------------------------------------------------------------------------------
	mat pointsInCircle(double _radius, colvec3 _center, double _azimut, double _zenit, unsigned _nPoints) {
		


		colvec3 n = { cos(_azimut)*sin(_zenit), sin(_zenit)*sin(_azimut), cos(_azimut) };
		colvec3 u = { -sin(_azimut), cos(_azimut), 0 };

		colvec3 nxu = { cos(_zenit)*cos(_azimut), cos(_zenit)*sin(_azimut), -sin(_zenit) };

		auto pCircle = [&](double t)->arma::colvec3 {
			return _radius*cos(t)*u + _radius*sin(t)*nxu + _center;
		};

		mat points(3, _nPoints);

		double incT = 2 * M_PI / _nPoints;
		for (unsigned i = 0; i < _nPoints; i++) {
			points.col(i) = pCircle(i*incT);
		}

		return points;
	}

	//--------------------------------------------------------------------------------------------------------------------
	mat rotUnitaryVectors(const vec & _a, const vec & _b) {
		auto a = _a / norm(_a);
		auto b = _b / norm(_b);

		if (norm(a - b) < 0.001) {
			return eye<mat>(3, 3);
		}


		vec v = cross(a, b);
		double s = norm(v);
		double c = dot(a, b);
		mat vx = { { 0, -v[2], v[1] },{ v[2], 0, -v[0] },{ -v[1], v[0], 0 } };

		mat R = eye<mat>(3, 3) + vx + vx*vx* (1 - c) / (s*s);

		return R;

	}

	//-----------------------------------------------------------------------------------------------------------------
	arma::colvec3 generateRandomPointSphere(const arma::colvec3 & _sphereCenter, const double _sphereRadius) {

		double theta = as_scalar(arma::randu(1, 1))*3.14159265359 * 2 - 3.14159265359;
		double phi = acos(2 * as_scalar(arma::randu(1, 1)) - 1);

		arma::colvec3 point, newPoint;
		point[0] = _sphereRadius * cos(theta)*sin(phi) + _sphereCenter[0];
		point[1] = _sphereRadius * sin(theta)*sin(phi) + _sphereCenter[1];
		point[2] = _sphereRadius * cos(phi) + _sphereCenter[2];

		return point;
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool intersectRayPlane(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _n, arma::colvec3 &_intersection) {
		double den = arma::dot(_n, (_p1 - _p0));
		if (den == 0) {
			return false;
		}
		else {
			double r = arma::dot(_n, (_v0 - _p0)) / den;
			_intersection = _p0 + r*(_p1 - _p0);
			return true;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool intersectRayTriangle(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _v1, arma::colvec3 _v2, arma::colvec3 &_intersection) {
		arma::colvec3 u = (_v1 - _v0);
		arma::colvec3 v = (_v2 - _v0);
		if (intersectRayPlane(_p0, _p1, _v1, arma::cross(u, v), _intersection)) {
			arma::colvec3 w = (_intersection - _v0);
			double s =	(arma::dot(u, v)*arma::dot(w, v) - arma::dot(v, v)*arma::dot(w, u)) /
						(arma::dot(u, v)*arma::dot(u, v) - arma::dot(u, u)*arma::dot(v, v));
			double t =	(arma::dot(u, v)*arma::dot(w, u) - arma::dot(u, u)*arma::dot(w, v)) /
						(arma::dot(u, v)*arma::dot(u, v) - arma::dot(u, u)*arma::dot(v, v));

			if (s >= 0 && t >= 0 && s + t <= 1) {
				return true;
			}else{
				return false;
			}
		}else {
			return false;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool intersectSegmentPlane(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _n, arma::colvec3 &_intersection) {
		double den = arma::dot(_n, (_p1 - _p0));
		if (den == 0) {
			return false;
		}
		else {
			double r = arma::dot(_n, (_v0 - _p0)) / den;
			if (r >= 0 && r <= 1) {
				_intersection = _p0 + r*(_p1 - _p0);
				return true;
			}else{
				return false;
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	bool intersectSegmentTriangle(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _v1, arma::colvec3 _v2, arma::colvec3 &_intersection){
		arma::colvec3 u = (_v1 - _v0);
		arma::colvec3 v = (_v2 - _v0);
		if (intersectSegmentPlane(_p0, _p1, _v1, arma::cross(u, v), _intersection)) {
			arma::colvec3 w = (_intersection - _v0);
			double s =	(arma::dot(u, v)*arma::dot(w, v) - arma::dot(v, v)*arma::dot(w, u)) /
						(arma::dot(u, v)*arma::dot(u, v) - arma::dot(u, u)*arma::dot(v, v));
			double t =	(arma::dot(u, v)*arma::dot(w, u) - arma::dot(u, u)*arma::dot(w, v)) /
						(arma::dot(u, v)*arma::dot(u, v) - arma::dot(u, u)*arma::dot(v, v));

			if (s >= 0 && t >= 0 && s + t <= 1) {
				return true;
			} else {
				return false;
			}
		}
		else {
			return false;
		}
	}
}
