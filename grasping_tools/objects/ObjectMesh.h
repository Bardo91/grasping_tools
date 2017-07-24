/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_OBJECTS_OBJECTMESH_H_
#define GPISGRASPING_OBJECTS_OBJECTMESH_H_

#include "../Object.h"

#include <pcl/PolygonMesh.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <vector>

namespace gpisGrasping {
	class ObjectMesh: public Object {
	public:
		/// Construct object from a file.
		/// \param _filename: path to file.
		ObjectMesh(std::string _filename);

		/// Construct object from vertices and faces indices.
		ObjectMesh(pcl::PointCloud<pcl::PointXYZ> &_vertices, std::vector<pcl::Vertices> &_faces);

		/// Return centroid 3D center of object.
		arma::colvec3 center();

		/// Intersect a directed ray with the mesh.
		/// \param _initpoint: point that defines the ray position.
		/// \param _dir: vector that defines the direction of the ray.
		arma::colvec3 intersect(arma::colvec3 _initPoint, arma::colvec3 _dir);

		/// Intersect a directed ray with the mesh.
		/// \param _initpoint: point that defines the ray position.
		/// \param _dir: vector that defines the direction of the ray.
		arma::mat intersectRay(arma::colvec3 _p1, arma::colvec3 _p2);

		/// Get pointcloud containing the centroid of the faces for any purpose.
		/// \return: centroid of faces.
		arma::mat centroidFaces();

		/// Get closest point on the mesh surface from a given point.
		/// \param _p: given point.
		/// \return: closest point to given point.
		pcl::PointXYZ closestVertex(const pcl::PointXYZ &_p);

		/// Get closest point on the mesh surface from a given point.
		/// \param _p: given point.
		/// \return: id of closest point to given point.
		int closestVertexId(const pcl::PointXYZ &_p);

		/// Get vertex point by id
		/// \param id
		/// \return point
		pcl::PointXYZ vertex(int _id);

		/// Get mesh info
		void mesh(pcl::PointCloud<pcl::PointXYZ> &_vertices, std::vector<pcl::Vertices> &_faces);

	private:
		pcl::PointCloud<pcl::PointXYZ> mVertices;
		std::vector<pcl::Vertices> mFaces;
		arma::mat mCentroidFaces;
		arma::colvec3 mCentroid;
		bool mIsKdtreeInit = false;
		pcl::KdTreeFLANN<pcl::PointXYZ> mKdTree;

	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_OBJECT_H_