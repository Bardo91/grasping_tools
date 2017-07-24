/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <objects/ObjectMesh.h>
#include <mathTools.h>

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/common/centroid.h>

//---------------------------------------------------------------------------------------------------------------------
gpisGrasping::ObjectMesh::ObjectMesh(std::string _filename) {
	// Load information from file.
	pcl::PolygonMesh mesh;
	if (_filename.find(".ply") != std::string::npos) {
		if(pcl::io::loadPLYFile(_filename, mesh) != 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else if (_filename.find(".obj") != std::string::npos) {
		if (pcl::io::loadOBJFile(_filename, mesh) != 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else if (_filename.find(".stl") != std::string::npos) {
		if (pcl::io::loadPolygonFileSTL(_filename, mesh) == 0) {
			std::cerr << "Error loading mesh of object" << std::endl;
			return;
		}
	}
	else{
		assert(false);
		std::cerr << "Unsupported filetype, please provide another filetype." << std::endl;
		return;
	}

	// Transform cloud.
	pcl::fromPCLPointCloud2(mesh.cloud, mVertices);
	mFaces = mesh.polygons;

	// Compute centroid.
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(mVertices, centroid);

	mCentroid[0] = centroid[0];
	mCentroid[1] = centroid[1];
	mCentroid[2] = centroid[2];
}

//---------------------------------------------------------------------------------------------------------------------
gpisGrasping::ObjectMesh::ObjectMesh(pcl::PointCloud<pcl::PointXYZ>& _vertices, std::vector<pcl::Vertices>& _faces) {
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
arma::colvec3 gpisGrasping::ObjectMesh::center() {
	return mCentroid;
}

//---------------------------------------------------------------------------------------------------------------------
arma::colvec3 gpisGrasping::ObjectMesh::intersect(arma::colvec3 _initPoint, arma::colvec3 _dir) {
	// Look for closest point.
	auto vertexId = closestVertexId(pcl::PointXYZ(_initPoint[0], _initPoint[1], _initPoint[2]));

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
arma::mat gpisGrasping::ObjectMesh::intersectRay(arma::colvec3 _p1, arma::colvec3 _p2) {
	arma::mat intersections;
	for (auto &face : mFaces) {
		auto v1 = mVertices[face.vertices[0]];
		auto v2 = mVertices[face.vertices[1]];
		auto v3 = mVertices[face.vertices[2]];

		arma::colvec3 intersection;
		if (intersectRayTriangle(_p1, _p2, {v1.x,v1.y,v1.z}, { v2.x,v2.y,v2.z }, { v3.x,v3.y,v3.z }, intersection)) {
            intersections.insert_cols(intersections.n_cols, intersection);
		}

	}
	return intersections;
}

//---------------------------------------------------------------------------------------------------------------------
arma::mat gpisGrasping::ObjectMesh::centroidFaces() {
	if (mCentroidFaces.n_cols == 0) {
		for (auto &face : mFaces) {
			auto p1 = mVertices[face.vertices[0]];
			auto p2 = mVertices[face.vertices[1]];
			auto p3 = mVertices[face.vertices[2]];

			arma::colvec6 point;
			point[0] = ( p1.x + p2.x + p3.x )/3;
			point[1] = ( p1.y + p2.y + p3.y )/3;
			point[2] = ( p1.z + p2.z + p3.z )/3;

			arma::colvec3 v1 = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
			arma::colvec3 v2 = { p3.x - p1.x, p3.y - p1.y, p3.z - p1.z };
			point.subvec(3, 5) = arma::cross(v1, v2);

            mCentroidFaces.insert_cols(mCentroidFaces.n_cols, point);
		}
		assert(mCentroidFaces.n_cols == mFaces.size());
	}

	return mCentroidFaces;
}

//---------------------------------------------------------------------------------------------------------------------
pcl::PointXYZ gpisGrasping::ObjectMesh::closestVertex(const pcl::PointXYZ &_p) {
	return mVertices[closestVertexId(_p)];
}

//---------------------------------------------------------------------------------------------------------------------
int gpisGrasping::ObjectMesh::closestVertexId(const pcl::PointXYZ &_p) {
	if (mIsKdtreeInit) {
		mKdTree.setInputCloud(mVertices.makeShared());
		mIsKdtreeInit = true;
	}
	int K = 1;
	std::vector<int> pointIdxNKNSearch(K);
	std::vector<float> pointNKNSquaredDistance(K);

	mKdTree.nearestKSearch(_p, K, pointIdxNKNSearch, pointNKNSquaredDistance);
	return pointIdxNKNSearch[0];
}

pcl::PointXYZ gpisGrasping::ObjectMesh::vertex(int _id) {
	return mVertices[_id];
}

//---------------------------------------------------------------------------------------------------------------------
void gpisGrasping::ObjectMesh::mesh(pcl::PointCloud<pcl::PointXYZ>& _vertices, std::vector<pcl::Vertices>& _faces) {
	_vertices = mVertices;
	_faces = mFaces;
}
