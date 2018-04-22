///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/Grasp.h>

#include <cassert>
#include <Eigen/Eigen>
#include <iostream>

#ifndef M_PI
	#define M_PI 3.14159265359
#endif



const double epsilon = 0.01;
const double maxError = 0.1;
int main(void) {

	std::vector<grasping_tools::ContactPoint> cps;
	cps.push_back(grasping_tools::ContactPoint({1,0,0}, arma::eye(3, 3),  {1,0,0}, arma::eye(3, 3), grasping_tools::eContactTypes::SFC, 1, 1, 0.5));
	cps.push_back(grasping_tools::ContactPoint({-1,0,0}, arma::eye(3, 3), {-1,0,0}, arma::eye(3, 3), grasping_tools::eContactTypes::SFC, 1, 1, 0.5));
	//cps.push_back(grasping_tools::ContactPoint({0,1,0}, arma::eye(3, 3),  {0,1,0}, arma::eye(3, 3), grasping_tools::eContactTypes::SFC, 1, 1, 1));
	//cps.push_back(grasping_tools::ContactPoint({0,-1,0}, arma::eye(3, 3), {0,-1,0}, arma::eye(3, 3), grasping_tools::eContactTypes::SFC, 1, 1, 1));

	grasping_tools::Grasp grasp;
	grasp.contactPoints(cps);

	std::cout << "lmrw: " << grasp.lmrw()  << std::endl;

	auto tasks = grasp.taskWrenches();

	for(auto &task: tasks){
		std::cout << task << ",";
	}
	std::cout << std::endl;

}
