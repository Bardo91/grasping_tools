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
	inline Grasp Arm3DOF::generate<ObjectGpis>(ObjectGpis &_object, arma::colvec3 &_pointIK) {
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
	inline Grasp Arm3DOF::generate<ObjectMesh>(ObjectMesh &_object, arma::colvec3 &_pointIK) {
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
        

        arma::mat intersections = _object.intersectRay(candidatePoint.rows(0, 2)+ 2*candidatePoint.rows(3, 5), candidatePoint.rows(0, 2) + candidatePoint.rows(3, 5));

        std::cout << " Matriz intersections: " << intersections << std::endl;

        if( intersections.n_cols == 2 ){
        // Cogemos punto medio
        arma::colvec3 puntoMedio = { (intersections(0,0) + intersections(0,1))/2, (intersections(1,0) + intersections(1,1))/2, (intersections(2,0) + intersections(2,1))/2 };

        std::vector<double> auxRadio(3);
        auxRadio.at(0) = intersections(0,1) - intersections(0,0);
        auxRadio.at(1) = intersections(1,1) - intersections(1,0);
        auxRadio.at(2) = intersections(2,1) - intersections(2,0);

        // Calculo Apertura garra
        double apert_garra = arma::norm( auxRadio);
        if(apert_garra > mAperture){
            std::cout << "Apertura mas grande de la permitida" << std::endl;
            return grasp;
        }
        else{
            std::cout << "Apertura permitida" << std::endl;
        }

        // CHECK THAT RESULT LIES OUTSIDE OF THE MESH
        QHULL_LIB_CHECK;
        pcl::PointCloud<pcl::PointXYZ> verticesObject;
        _object.vertices(verticesObject);
        // Compute convex hull
        orgQhull::PointCoordinates points(3, "Object convex hull");
        std::vector<double> concatenationPoints;
        for (auto &p:verticesObject) {
            concatenationPoints.push_back(p.x);
            concatenationPoints.push_back(p.y);
            concatenationPoints.push_back(p.z);
        }
        points.append(concatenationPoints);
        try {
            orgQhull::Qhull convexHull("", 3, wrenchesCones.n_cols, &concatenationPoints[0], "Qt");
            orgQhull::QhullFacetList facets = convexHull.facetList();
            bool isOutside = true;
            for (orgQhull::QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it) {
                orgQhull::QhullFacet f = *it;
                if (!f.isGood()) continue;
                
            }
        }catch (orgQhull::QhullError e) {
            std::cout << "Error computing convex hull. Error msg: \n" << e.what()  << std::endl;
            return Grasp;
        }

        // Obtencion vectores N y U para Circulo 3D
        std::vector<double> vec_n(3);
        vec_n.at(0) = intersections(0,1) - auxRadio.at(0);
        vec_n.at(1) = intersections(1,1) - auxRadio.at(1);
        vec_n.at(2) = intersections(2,1) - auxRadio.at(2);

        std::vector<double> proy_nu(3);
        arma::mat vec_r = arma::eye(3, 3);
        double min_producto = 9999999;

        for(unsigned int k = 0; k < 3; k++){
            double producto_e = vec_r.at(k,0)*vec_n.at(0) + vec_r.at(k,1)*vec_n.at(1) + vec_r.at(k,2)*vec_n.at(2);

            if( producto_e < min_producto ){
                min_producto = producto_e;
                proy_nu.at(0) = producto_e*vec_n.at(0);
                proy_nu.at(1) = producto_e*vec_n.at(1);
                proy_nu.at(2) = producto_e*vec_n.at(2);
            }
        }

        std::vector<double> vec_u(3);
        vec_u.at(0) = vec_n.at(0) -  proy_nu.at(0);
        vec_u.at(1) = vec_n.at(1) -  proy_nu.at(1);
        vec_u.at(2) = vec_n.at(2) -  proy_nu.at(2);

        double mod_u = sqrt( vec_n.at(0)*vec_n.at(0) + vec_n.at(1)*vec_n.at(1) + vec_n.at(2)*vec_n.at(2) );
        vec_u.at(0) = vec_u.at(0)/mod_u;
        vec_u.at(1) = vec_u.at(1)/mod_u;
        vec_u.at(2) = vec_u.at(2)/mod_u;

        // Creacion Circulo 3D y obtencion de sus puntos
        int npCircle = 40;
        double radioPuntos = 0.05;
        arma::mat puntosCirculo(3, npCircle);
        arma::colvec3 vec_n_arma = {vec_n.at(0), vec_n.at(1), vec_n.at(2)};
        arma::colvec3 vec_u_arma = {vec_u.at(0), vec_u.at(1), vec_u.at(2)};
        puntosCirculo = pointsInCircleNU(radioPuntos, puntoMedio, vec_n_arma, vec_u_arma, npCircle);

        std::vector<double> validPoints;

        for(unsigned int i = 0; i < npCircle; i++){

            // Calculate inverse kinematic for points
            std::vector<int> IKArm(3);
            float valx = puntosCirculo.at(0,i);
            float valy = puntosCirculo.at(1,i);
            float valz = puntosCirculo.at(2,i);

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
                //std::cout << "Valid point" << std::endl;
                validPoints.push_back(puntosCirculo.at(0,i));
                validPoints.push_back(puntosCirculo.at(1,i));
                validPoints.push_back(puntosCirculo.at(2,i));
            }
            else{
                //std::cout << "Invalid point" << std::endl;
            }
            
        }

        // De los que obtenemos,
        // nos quedamos con el que menor distancia tenga a la posicion home del brazo
        // Home brazo -> (0.35, 0, 0.23)
        double min_dist = 9999999;
        arma::colvec3 pointFinal;

        for(unsigned int j = 0; j < validPoints.size() ; j=j+3){

            std::vector<double> auxDist(3);
            auxDist.at(0) = (validPoints.at(j) - 0.35)*(validPoints.at(j) - 0.35);
            auxDist.at(1) = (validPoints.at(j+1) - 0.0)*(validPoints.at(j+1) - 0.0);
            auxDist.at(2) = (validPoints.at(j+2) - 0.23)*(validPoints.at(j+2) - 0.23);
            double dist = sqrt(auxDist.at(0) + auxDist.at(1) + auxDist.at(2));

            if(dist < min_dist){
                pointFinal = { validPoints.at(j), validPoints.at(j+1), validPoints.at(j+2)};
                min_dist = dist;
            }
        }

        // return point final IK
        _pointIK = pointFinal;
        std::cout << "Punto final: " << pointFinal << std::endl;

       for(unsigned int i = 0; i < intersections.n_cols; i++){
            ContactPoint cp(intersections.col(i).head(3), arma::eye(3, 3), intersections.col(i).tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1);
            cps.push_back(cp);
        }
        grasp.contactPoints(cps);

        }

		return grasp;
	}

}
