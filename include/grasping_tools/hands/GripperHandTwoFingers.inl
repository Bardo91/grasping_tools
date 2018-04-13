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

#include <hecatonquiros/model_solvers/ModelSolverOpenRave.h>

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
	inline  std::vector<Grasp> GripperHandTwoFingers::generateGrasps<ObjectMesh>(ObjectMesh &_object, double _resolution){

		// Get boundary of object
		arma::colvec3 min, max;
		_object.minMax(min, max);
		arma::colvec3 sizes = 	{(max[0] - min[0]), (max[1]-min[1]), (max[2] - min[2])};
		//arma::colvec3 center = 	{max[0] - min[0], max[1]-min[1], max[2] - min[2]}; center /= 2;

		// Generate N points equally distributed on the faces of a centered cube.
		arma::mat initPointsPlusNormals;
		
		int nX = (sizes[0])/_resolution;
		int nY = (sizes[1])/_resolution;
		int nZ = (sizes[2])/_resolution;

		// (-)X normal faces
		for(unsigned i = 0; i < nY; i++){
			for(unsigned j = 0; j < nZ; j++){
				arma::colvec6 pointPlusNormal = {-sizes[0]/2,
												-nY*_resolution/2 + i*_resolution,
												-nZ*_resolution/2 + j*_resolution,
												1,
												0,
												0};
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
				pointPlusNormal[0]*=-1;
				pointPlusNormal[3]*=-1;
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
			}	
		}

		// (-)Y normal faces
		for(unsigned i = 0; i < nX; i++){
			for(unsigned j = 0; j < nZ; j++){
				arma::colvec6 pointPlusNormal = {-nX*_resolution/2 + i*_resolution,
												-sizes[1]/2,
												-nZ*_resolution/2 + j*_resolution,
												0,
												1,
												0};
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
				pointPlusNormal[1]*=-1;
				pointPlusNormal[4]*=-1;
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
			}	
		}

		// (-)Z normal faces
		for(unsigned i = 0; i < nX; i++){
			for(unsigned j = 0; j < nY; j++){
				arma::colvec6 pointPlusNormal = {-nX*_resolution/2 + i*_resolution,
												-nY*_resolution/2 + j*_resolution,
												-sizes[2]/2,
												0,
												0,
												1};
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
				pointPlusNormal[2]*=-1;
				pointPlusNormal[5]*=-1;
				initPointsPlusNormals.insert_cols(initPointsPlusNormals.n_cols, pointPlusNormal);
			}	
		}

		// Rotate points with pose
		arma::mat44 objectPose = _object.move();
		initPointsPlusNormals.rows(0,2) = objectPose.cols(0,2).rows(0,2)*initPointsPlusNormals.rows(0,2);
		for(unsigned i = 0; i < initPointsPlusNormals.n_cols;i++){
			initPointsPlusNormals.col(i).rows(0,2) +=  objectPose.col(3).rows(0,2);
		}
		initPointsPlusNormals.rows(3,5) = objectPose.cols(0,2).rows(0,2)*initPointsPlusNormals.rows(3,5);

		//#ifdef HAS_OPENRAVE
		//	//auto boxHandle = hecatonquiros::ModelSolverOpenRave::getEnvironment()->drawbox(	OpenRAVE::RaveVector<float>(objectPose.col(3)[0], objectPose.col(3)[1], objectPose.col(3)[2]),
		//	//																	OpenRAVE::RaveVector<float>(sizes[0]/2, sizes[1]/2, sizes[2]/2));
		//	std::vector<OpenRAVE::GraphHandlePtr> raysHandle;
		//#endif


		// Generate candidate grasps
		std::vector<Grasp> grasps;

		// Trace rays
		for(unsigned i = 0; i < initPointsPlusNormals.n_cols; i++) {
			arma::colvec6 initPointPlusNormal = initPointsPlusNormals.col(i);

			//#ifdef HAS_OPENRAVE
			//	OpenRAVE::RaveVector<float> p1 = {initPointPlusNormal[0], initPointPlusNormal[1], initPointPlusNormal[2]}, p2;
			//	OpenRAVE::RaveVector<float> dir = {initPointPlusNormal[3], initPointPlusNormal[4], initPointPlusNormal[5]};
			//	p2 = p1 + dir*0.03;
			//	raysHandle.push_back(hecatonquiros::ModelSolverOpenRave::getEnvironment()->drawarrow(p1, p2,0.001,OpenRAVE::RaveVector< float >(1, 0, 0, 1)));
			//#endif
			arma::mat intersections = _object.intersectRay(initPointPlusNormal.rows(0,2), initPointPlusNormal.rows(0,2)+initPointPlusNormal.rows(3,5));
			
			if(intersections.n_cols > 1 && intersections.n_cols%2 == 0){
				//	Example of ray tracing. Let be a ray tracing an object, it generates 4 intersections called a, b, c and d.
				//  The number of intersections should be always even.
				//  From the figure it is intuitive that 3 "compressing" grasps are possible (a,b), (c,d) and (a,d).
				//  If intersections are arranged in distance from the source. 
				//
				//					_________		___________
				//				 a	|		| b	  c	|		   | d
				//	>---------	 <--|		|--> <--|		   |-->
				//					|		|		|		   |
				//					|		|_______|		   |
				//					|				 		   |
				//					|__________________________|
				//
				// In this other example we have (a,b), (c,d), (e,f), (a,d), (c,f), (a,f)
				//					_________		___________				  _______________
				//				 a	|		| b	  c	|		   | d			e |				| f
				//	>---------	 <--|		|--> <--|		   |-->		   <--|				|-->
				//					|		|		|		   |			  |				|
				//					|		|_______|		   |______________|				|
				//					|				 		   								|
				//					|_______________________________________________________|
				//
				// In general terms, 	N_intersections = N_folds*2
				//						N_grasps = SUM(i=1 to N_folds) { i }

				// Compute distances to source point.
				std::vector<std::pair<int, float>> arrangedIds(intersections.n_cols);
				
				for(unsigned idx = 0; idx < intersections.n_cols; idx++){
					arrangedIds[idx] = {idx, arma::norm(initPointPlusNormal.head(3) - intersections.col(idx).head(3))};
				}

				// Sort intersections by distance.
				std::sort(arrangedIds.begin(), arrangedIds.end(), 
						[](std::pair<int, float> &_a, std::pair<int, float> &_b){ return _a.second < _b.second; });

				int nFolds = intersections.n_cols/2;
				for(unsigned jumper = 1; jumper <= nFolds; jumper+=2){
					for(unsigned initId = 0; initId +jumper < intersections.n_cols; initId += 2){
						auto p1 = intersections.col(initId);
						auto p2 = intersections.col(initId+jumper);
						
						// Check aperture of gripper
						if(arma::norm(p1.head(3)-p2.head(3)) < mAperture){	// 666 TODO take into account finger size
							std::vector<ContactPoint> cps;
							cps.push_back(ContactPoint(p1.head(3), arma::eye(3, 3), p1.tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1));
							cps.push_back(ContactPoint(p2.head(3), arma::eye(3, 3), p2.tail(3), arma::eye(3, 3), eContactTypes::SFC, 1, 1, 1));
							
							Grasp grasp;
							grasp.contactPoints(cps);
							grasps.push_back(grasp);
						}
					}
				}
			}
		}
		return grasps;
	}

}
