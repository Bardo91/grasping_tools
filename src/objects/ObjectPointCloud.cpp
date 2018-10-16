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

#include <thread>
#include <chrono>

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectPointCloud::ObjectPointCloud(std::string _filename) {
    // Load information from file.
	std::cout << "Loading file: " << _filename << std::endl;
	pcl::PCLPointCloud2 cloud;
	if (_filename.find(".pcd") != std::string::npos) {
        if (pcl::io::loadPCDFile(_filename, cloud) == 0) {
			std::cerr << "Error loading point cloud of object" << std::endl;
			return;
		}else{
			pcl::fromPCLPointCloud2(cloud, mVertices);
			#warning Cloud is assumed to have normals !!! Need to fix
			// if(cloud.fields.contains("curvature")){
			// 	// do nothing
			// }else{
			// 	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
  			// 	ne.setInputCloud (cloud);
			// 	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
			// 	ne.setSearchMethod (tree);
			// 	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
			// 	// Use all neighbors in a sphere of radius 3cm
			// 	ne.setRadiusSearch (0.03);
			// 	// Compute the features
			// 	ne.compute (*cloud_normals);
			// 	for(unsigend i = 0; i< cloud.size();i++){
			// 		cloud.normal_x = cloud_normals->points[i].normal_x;
			// 		cloud.normal_y = cloud_normals->points[i].normal_y;
			// 		cloud.normal_z = cloud_normals->points[i].normal_z;
			// 	}
			// }

		}
	}else if(_filename.find(".ply") != std::string::npos){
		if (pcl::io::loadPLYFile(_filename, mVertices) == 0) {
			std::cerr << "Error loading point cloud of object" << std::endl;
			return;
		}else{
			#warning Cloud is assumed to have normals !!! Need to fix
			// if(cloud.fields.contains("curvature")){
			// 	// do nothing
			// }else{
			// 	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
			// 	ne.setInputCloud (cloud);
			// 	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
			// 	ne.setSearchMethod (tree);
			// 	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
			// 	// Use all neighbors in a sphere of radius 3cm
			// 	ne.setRadiusSearch (0.03);
			// 	// Compute the features
			// 	ne.compute (*cloud_normals);
			// 	for(unsigend i = 0; i< cloud.size();i++){
			// 		cloud.normal_x = cloud_normals->points[i].normal_x;
			// 		cloud.normal_y = cloud_normals->points[i].normal_y;
			// 		cloud.normal_z = cloud_normals->points[i].normal_z;
			// 	}
			// }
		}
	}
	else{
		std::cerr << "Unsupported filetype, please provide another filetype." << std::endl;
		assert(false);
		return;
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
    pcl::PointNormal initPoint;
	initPoint.x = _p1[0];
	initPoint.y = _p1[1];
	initPoint.z = _p1[2];

	auto vertex = mVertices[closestVertexId(initPoint)];
	
	arma::mat intersections = {	vertex.x,vertex.y, vertex.z,
								vertex.normal_x, vertex.normal_y, vertex.normal_z};
	
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