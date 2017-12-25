///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/qualityMetrics/DistancePointConvexCone.h>

#include <cassert>
#include <Eigen/Eigen>
#include <iostream>

#ifndef M_PI
	#define M_PI 3.14159265359
#endif


void test3dCone();
void test6dConeSimple();
void test6dCone();

const double epsilon = 0.01;
const double maxError = 0.1;
int main(void) {

	test3dCone();
	test6dConeSimple();
	test6dCone();

}

//---------------------------------------------------------------------------------------------------------------------
void test3dCone() {
	// Create support points for the cone
	arma::mat coneA;
	const unsigned nPoints = 100;
	const double radius = 1;
	const double incAngle = 2 * M_PI / (nPoints + 1);
	for (unsigned i = 0; i < nPoints ; i++) {
		double x = radius*cos(incAngle*i);
		double y = radius*sin(incAngle*i);
		
		arma::colvec p = { x, y + 1, 1 };
		coneA.insert_cols(coneA.n_cols, p);
	}

	grasping_tools::DistancePointConvexCone distCone;
	distCone.basePoints(coneA);
	distCone.epsilon(epsilon);

	// Test point [0,0,1] --> distance should be 0
	arma::colvec b = { 0, 0, 1 };
	distCone.point(b);
	distCone.compute();
	assert(abs(distCone.minDistance() - 0) <= maxError);
	// Test point [0,2,1] --> distance should be 0
	b = {0, 2, 1};
	distCone.point(b);
	distCone.compute();
	assert(abs(distCone.minDistance() - 0) <= maxError);
	// Test point [0,-1,1] --> distance should be 1
	b = {0, -1, 1};
	distCone.point(b);
	distCone.compute();
	assert(abs(distCone.minDistance() - 1) <= maxError);
	// Test point [0,-2,1] --> distance should be 2
	b = {0, -2, 1};
	distCone.point(b);
	distCone.compute();
	assert(abs(distCone.minDistance() - 2) <= maxError);
	// Test point [1,1,1] --> distance should be 0
	b = {1, 1, 1};
	distCone.point(b);
	distCone.compute();
	assert(abs(distCone.minDistance() - 0) <= maxError);
}

//---------------------------------------------------------------------------------------------------------------------
void test6dConeSimple() {

}

//---------------------------------------------------------------------------------------------------------------------
void test6dCone() {

}
