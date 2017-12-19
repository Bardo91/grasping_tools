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
#include <deque>
#include <vector>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/RboxPoints.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullQh.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullLinkedList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullSet.h>
#include <libqhullcpp/QhullVertexSet.h>

#include <math.h>

namespace grasping_tools {
	template<>
	inline Grasp Arm3DOF::generate<ObjectGpis>(ObjectGpis &_object) {
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
	inline Grasp Arm3DOF::generate<ObjectMesh>(ObjectMesh &_object) {
        // Generate grasp
		Grasp grasp;
        std::vector<ContactPoint> cps;

        float radio, humero, base;
        radio = 0.35;
        humero = 0.15;
        base = 0.08;
/*
        // Workspace Arm
		std::vector<double> WSA;
		int cont = 0;
        for(int theta1 = -90; theta1 < 91; theta1 = theta1 + 10){
        
            for(int theta2 = 0; theta2 < 91; theta2 = theta2 + 10){
        
                for(int theta3 = 0; theta3 < 91; theta3 = theta3 + 10){
                    std::vector<double> TA0;
                    TA0.at(0) = radio*(cos(theta1*M_PI/180)*cos(theta2*M_PI/180 + M_PI/2)*cos(theta3*M_PI/180) - cos(theta1*M_PI/180)*sin(theta2*M_PI/180 + M_PI/2)*sin(theta3*M_PI/180) ) + humero*(cos(theta2*M_PI/180 + M_PI/2)*cos(theta3*M_PI/180));
                    TA0.at(1) = radio*(sin(theta1*M_PI/180)*cos(theta2*M_PI/180 + M_PI/2)*cos(theta3*M_PI/180)) + humero*sin(theta1*M_PI/180)*cos(theta2*M_PI/180 + M_PI/2);
                    TA0.at(2) = radio*(sin(theta2*M_PI/180 + M_PI/2)*cos(theta3*M_PI/180) + cos(theta2*M_PI/180 + M_PI/2)*sin(theta3*M_PI/180)) + humero*(sin(theta2*M_PI/180 + M_PI/2)) + base;
        
                    WSA.push_back(TA0.at(0));
                    WSA.push_back(TA0.at(1));
                    WSA.push_back(TA0.at(2));
					
					cont++;

                }
            }
        }
    */

        auto candidatePoints = _object.centroidFaces();

        // Get a random candidate point // 666 TODO improve candidate selection
        arma::imat id = arma::randi(1, 1, arma::distr_param(0, candidatePoints.n_cols - 1));
        arma::colvec6 candidatePoint = candidatePoints.col(id(0,0));

        arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2)+ 2*candidatePoint.rows(3, 5), candidatePoint.rows(0, 2) - candidatePoint.rows(3, 5));

        //std::cout << " Matriz intersection: " << intersections << std::endl;

        for(unsigned int i = 0; i < intersections.n_cols; i++){


            // Calculate inverse kinematic for centroid point
            std::vector<int> IKArm(3);
            float valx = intersections.at(0,i);
            float valy = intersections.at(1,i);
            float valz = intersections.at(2,i);


            IKArm.at(0) = atan2(valy, valx);
 
            float y = sqrt(valx*valx + valy*valy);
            float x = valz - base;
        
            float D = (x*x+y*y-humero*humero-radio*radio)/(2*humero*radio);
            float auxD = 1 - D*D;
            if(auxD < 0){
                auxD = 0;
                std::cout << "1-D^2 < 0, unrecheable point, fitting to 0\n";
            }
        
            float theta2a = atan2(sqrt(auxD), D);
        
            float k1a = humero + radio*cos(theta2a);
            float k2a = radio*sin(theta2a);
            float theta1a = atan2(y,x)-atan2(k2a,k1a);
        
            IKArm.at(1) = theta1a;
            IKArm.at(2) = theta2a;

            // If these point is a valid angle for Workspace Arm, add it to ContactPoint
            if((IKArm.at(0)>=(-90*M_PI/180) && IKArm.at(0)<=(90*M_PI/180)) && (IKArm.at(1)>=(-90*M_PI/180) && IKArm.at(1)<=(90*M_PI/180)) && (IKArm.at(2)>=(-90*M_PI/180) && IKArm.at(2)<=(90*M_PI/180)) ){
                std::cout << "Valid point" << std::endl;
                ContactPoint cp(intersections.col(i).head(3), arma::eye(3, 3), intersections.col(i).tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1);
                cps.push_back(cp);
            }
            else{
                std::cout << "Invalid point" << std::endl;
            }
            
        }

        grasp.contactPoints(cps);

		return grasp;
	}

}
