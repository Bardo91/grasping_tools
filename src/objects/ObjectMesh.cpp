/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/objects/ObjectMesh.h>
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
grasping_tools::ObjectMesh::ObjectMesh(std::string _filename) {
    // Load information from file.
	std::cout << "Loading file: " << _filename << std::endl;
	if (_filename.find(".ply") != std::string::npos) {
        if(pcl::io::loadPLYFile(_filename, mMesh) != 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else if (_filename.find(".obj") != std::string::npos) {
        if (pcl::io::loadOBJFile(_filename, mMesh) != 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else if (_filename.find(".stl") != std::string::npos) {
        if (pcl::io::loadPolygonFileSTL(_filename, mMesh) == 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else{
		std::cerr << "Unsupported filetype, please provide another filetype." << std::endl;
		assert(false);
		return;
	}

	// Transform cloud.
    pcl::fromPCLPointCloud2(mMesh.cloud, mVertices);
    mFaces = mMesh.polygons;

	// Compute centroid.
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

	mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
	
	pcl::getMinMax3D(mVertices, mMinPoint, mMaxPoint);

    //viewer.registerKeyboardCallback(&grasping_tools::ObjectMesh::callbackKeyboard3dViewer,  *this, (void*)&viewer);
}

//---------------------------------------------------------------------------------------------------------------------
grasping_tools::ObjectMesh::ObjectMesh(pcl::PointCloud<pcl::PointNormal>& _vertices, std::vector<pcl::Vertices>& _faces) {
	mVertices = _vertices;
	mFaces = _faces;

	// Compute centroid.
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

    mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 grasping_tools::ObjectMesh::center() {
	return mCentroid;
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectMesh::minMax(arma::colvec3 &_min, arma::colvec&_max){
	_min = {mMinPoint.x, mMinPoint.y, mMinPoint.z};
	_max = {mMaxPoint.x, mMaxPoint.y, mMaxPoint.z};
}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 grasping_tools::ObjectMesh::intersect(arma::colvec3 _initPoint, arma::colvec3 _dir) {
	// Look for closest point.
    pcl::PointNormal pInit; pInit.x = _initPoint[0]; pInit.y = _initPoint[1]; pInit.z = _initPoint[2];
    auto vertexId = closestVertexId(pInit);

	assert(vertexId != -1);

	std::vector<pcl::Vertices> facesSharedVertex;
	// Choose closest faces.
	for (auto &face : mFaces) {
		if (std::find(face.vertices.begin(), face.vertices.end(), vertexId) != face.vertices.end()) {
			facesSharedVertex.push_back(face);
		}
	}

	// Intersect planes.
	double minDist = std::numeric_limits<double>::max();
	arma::colvec3 closestPoint;
	for (auto &face : facesSharedVertex) {
		auto v1 = mVertices[face.vertices[0]];
		auto v2 = mVertices[face.vertices[1]];
		auto v3 = mVertices[face.vertices[2]];

		arma::colvec3 intersection;
		if (intersectRayTriangle(_initPoint, _initPoint + _dir, { v1.x,v1.y,v1.z }, { v2.x,v2.y,v2.z }, { v3.x,v3.y,v3.z }, intersection)) {
			double dist = arma::norm(_initPoint - intersection);
			if (dist < minDist) {
				minDist = dist;
				closestPoint = intersection;
			}
		}
		
	}

	// Return point
	return closestPoint;
}

//---------------------------------------------------------------------------------------------------------------------
arma::mat grasping_tools::ObjectMesh::intersectRay(arma::colvec3 _p1, arma::colvec3 _p2) {
	arma::mat intersections;
    for (auto &face : mFaces) { // 666 TODO pending parallelization
		auto v1 = mVertices[face.vertices[0]];
		auto v2 = mVertices[face.vertices[1]];
        auto v3 = mVertices[face.vertices[2]];

		arma::colvec3 intersection;
        arma::colvec ver1 = {v1.x,v1.y,v1.z}, ver2 = { v2.x,v2.y,v2.z }, ver3 = { v3.x,v3.y,v3.z };
        if (intersectRayTriangle(_p1, _p2, ver1, ver2, ver3, intersection)) {
            arma::colvec6 point;
            point.head(3) = intersection;
            point.tail(3) = -arma::cross(ver2-ver3,ver1-ver3);
            point.tail(3) /= arma::norm(point.tail(3));
            intersections.insert_cols(intersections.n_cols, point);
		}


	}
	return intersections;
}

//---------------------------------------------------------------------------------------------------------------------
arma::mat grasping_tools::ObjectMesh::centroidFaces() {
	if (mCentroidFaces.n_cols == 0) {
		for (auto &face : mFaces) {
            auto pclP1 = mVertices[face.vertices[0]];
            auto pclP2 = mVertices[face.vertices[1]];
            auto pclP3 = mVertices[face.vertices[2]];

            arma::colvec6 point;
            arma::colvec3 p1 = {pclP1.x, pclP1.y, pclP1.z};
            arma::colvec3 p2 = {pclP2.x, pclP2.y, pclP2.z};
            arma::colvec3 p3 = {pclP3.x, pclP3.y, pclP3.z};

            point.subvec(0,2) = (p1+p2+p3)/3;
            point.subvec(3,5) = -arma::cross(p2-p3,p1-p3);

            mCentroidFaces.insert_cols(mCentroidFaces.n_cols, point);
		}
		assert(mCentroidFaces.n_cols == mFaces.size());
	}

	return mCentroidFaces;
}

//---------------------------------------------------------------------------------------------------------------------
pcl::PointNormal grasping_tools::ObjectMesh::closestVertex(const pcl::PointNormal &_p) {
	return mVertices[closestVertexId(_p)];
}

//---------------------------------------------------------------------------------------------------------------------
int grasping_tools::ObjectMesh::closestVertexId(const pcl::PointNormal &_p) {
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

pcl::PointNormal grasping_tools::ObjectMesh::vertex(int _id) {
	return mVertices[_id];
}

//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectMesh::mesh(pcl::PointCloud<pcl::PointNormal>& _vertices, std::vector<pcl::Vertices>& _faces) {
	_vertices = mVertices;
	_faces = mFaces;
}


//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectMesh::mesh(pcl::PolygonMesh &_mesh){
    _mesh = mMesh;
}


//---------------------------------------------------------------------------------------------------------------------
void grasping_tools::ObjectMesh::moveObject(){
	pcl::PointCloud<pcl::PointNormal> points;
	pcl::fromPCLPointCloud2(mMesh.cloud, points);
	pcl::transformPointCloud(points, points, mPoseEigen);
	pcl::toPCLPointCloud2(points, mMesh.cloud);
	mVertices = points;

	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

	mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
	
	pcl::getMinMax3D(points, mMinPoint, mMaxPoint);
}