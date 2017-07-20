#include "utilsMath.h"
#include <cassert>

using namespace arma;

mat gpis::rotUnitaryVectors(const vec & _a, const vec & _b) {
	auto a = _a / norm(_a);
	auto b = _b / norm(_b);

	if (norm(a - b) < 0.001) {
		return eye<mat>(3, 3);
	}


	vec v = cross(a, b);
	double s = norm(v);
	double c = dot(a, b);
	mat vx = { {0, -v[2], v[1]}, {v[2], 0, -v[0]}, {-v[1], v[0], 0} };
	
	mat R = eye<mat>(3,3) + vx + vx*vx* (1 - c) / (s*s);

	return R;

}
