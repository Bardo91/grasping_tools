/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/objects/ObjectPointCloud.h>
#include <grasping_tools/mathTools.h>

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_cloud.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/octree/octree.h>
#include <pcl/features/normal_3d.h>

#include <thread>
#include <chrono>
#include <algorithm>

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectPointCloud::ObjectPointCloud(std::string _filename) {
    // Load information from file.
	std::cout << "Loading file: " << _filename << std::endl;
	pcl::PCLPointCloud2 cloud;
	if (_filename.find(".pcd") != std::string::npos) {
        if (pcl::io::loadPCDFile(_filename, cloud) != 0) {
			std::cerr << "Error loading point cloud of object" << std::endl;
			return;
		}
	}else if(_filename.find(".ply") != std::string::npos){
		if (pcl::io::loadPLYFile(_filename, cloud) != 0) {
			std::cerr << "Error loading point cloud of object" << std::endl; 
			return;
		}
	}else if(_filename.find(".obj") != std::string::npos){
		if (pcl::io::loadOBJFile(_filename, cloud) != 0) {
			std::cerr << "Error loading point cloud of object" << std::endl; 
			return;
		}
	}
	else{
		std::cerr << "Unsupported filetype, please provide another filetype." << std::endl;
		assert(false);
		return;
	}

	pcl::fromPCLPointCloud2(cloud, mVertices);
	bool hasCurvature = false;
	for(auto &field: cloud.fields){
		if(field.name == "curvature")
			hasCurvature = true;
	}
	if(!hasCurvature){
		std::cout << "Computing curvature" << std::endl;
		pcl::NormalEstimation<pcl::PointNormal, pcl::Normal> ne;
		ne.setInputCloud (mVertices.makeShared());
		pcl::search::KdTree<pcl::PointNormal>::Ptr tree (new pcl::search::KdTree<pcl::PointNormal> ());
		ne.setSearchMethod (tree);
		pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
		// Use all neighbors in a sphere of radius 3cm
		ne.setRadiusSearch (0.015);	// 666 Make adaptative value
		// Compute the features
		ne.compute (*cloud_normals);
		for(unsigned i = 0; i< mVertices.size();i++){
			mVertices[i].normal_x = -cloud_normals->points[i].normal_x;
			mVertices[i].normal_y = -cloud_normals->points[i].normal_y;
			mVertices[i].normal_z = -cloud_normals->points[i].normal_z;
		}
	}

	// Compute centroid.
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

	mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
	
	pcl::getMinMax3D(mVertices, mMinPoint, mMaxPoint);
}

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectPointCloud::ObjectPointCloud(pcl::PointCloud<pcl::PointNormal>& _vertices) {
	mVertices = _vertices;

	// Compute centroid.
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

    mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 grasping_tools::ObjectPointCloud::center() {
	return mCentroid;
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectPointCloud::minMax(arma::colvec3 &_min, arma::colvec&_max){
	_min = {mMinPoint.x, mMinPoint.y, mMinPoint.z};
	_max = {mMaxPoint.x, mMaxPoint.y, mMaxPoint.z};
}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 grasping_tools::ObjectPointCloud::intersect(arma::colvec3 _initPoint, arma::colvec3 _dir) {
	pcl::PointNormal initPoint;
	initPoint.x = _initPoint[0];
	initPoint.y = _initPoint[1];
	initPoint.z = _initPoint[2];

	auto vertex = mVertices[closestVertexId(initPoint)];
	
	arma::colvec3 closestPoint = {	vertex.x,vertex.y, vertex.z};

	// Return point
	return closestPoint;
}

//---------------------------------------------------------------------------------------------------------------------
arma::mat grasping_tools::ObjectPointCloud::intersectRay(arma::colvec3 _p1, arma::colvec3 _p2) {
	pcl::octree::OctreePointCloudSearch<pcl::PointNormal> rayTracer(0.01);
    rayTracer.setInputCloud(mVertices.makeShared());
    rayTracer.addPointsFromInputCloud();

	Eigen::Vector3f cameraOrigin =  {_p1[0], _p1[1], _p1[2]};
	Eigen::Vector3f rayEnd = {_p2[0], _p2[1], _p2[2]};

	Eigen::Vector3f rayDir = rayEnd - cameraOrigin;

	std::vector<int> indices;
	rayTracer.getIntersectedVoxelIndices(cameraOrigin, rayDir, indices);
	
	arma::mat intersections;

	for(unsigned i = 0; i < indices.size(); i++){
		auto p = mVertices[i];
		arma::colvec6 point;
		point[0] = p.x;
		point[1] = p.y;
		point[2] = p.z;
		point[3] = p.normal_x;
		point[4] = p.normal_y;
		point[5] = p.normal_z;
		point.tail(3) /= arma::norm(point.tail(3));
		intersections.insert_cols(intersections.n_cols, point);
	}

	return intersections;
}

//---------------------------------------------------------------------------------------------------------------------
pcl::PointNormal grasping_tools::ObjectPointCloud::closestVertex(const pcl::PointNormal &_p) {
	return mVertices[closestVertexId(_p)];
}

//---------------------------------------------------------------------------------------------------------------------
int grasping_tools::ObjectPointCloud::closestVertexId(const pcl::PointNormal &_p) {
	if (!mIsKdtreeInit) {
		std::cout << mVertices.size() << std::endl;
		if(mVertices.size() == 0){
			return -1;
		}
		mKdTree.setInputCloud(mVertices.makeShared());
		mIsKdtreeInit = true;
	}
	int K = 1;
	std::vector<int> pointIdxNKNSearch(K);
	std::vector<float> pointNKNSquaredDistance(K);

	mKdTree.nearestKSearch(_p, K, pointIdxNKNSearch, pointNKNSquaredDistance);
	return pointIdxNKNSearch[0];
}

//---------------------------------------------------------------------------------------------------------------------
pcl::PointNormal grasping_tools::ObjectPointCloud::vertex(int _id) {
	return mVertices[_id];
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectPointCloud::moveObject(){
	pcl::transformPointCloud(mVertices, mVertices, mPoseEigen);
	
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

	mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
	
	pcl::getMinMax3D(mVertices, mMinPoint, mMaxPoint);
}