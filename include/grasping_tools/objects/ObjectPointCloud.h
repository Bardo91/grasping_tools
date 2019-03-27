/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_OBJECTS_OBJECTPOINTCLOUD_H_
#define GPISGRASPING_OBJECTS_OBJECTPOINTCLOUD_H_

#include <grasping_tools/Object.h>

#include <pcl/PolygonMesh.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <vector>
//#include <pcl/visualization/pcl_visualizer.h>


namespace grasping_tools {
	/// Class that inherit from Object parent class. This kind of objects uses a point cloud as model for object's surface
	///
	/// Class defined in `#include <grasping_tools/objects/ObjecPointCloud.h>`
	class ObjectPointCloud: public Object {
	public:
		/// Construct object from a file.
		/// \param _filename: path to file.
		ObjectPointCloud(std::string _filename);

		/// Construct object from vertices and faces indices.
        ObjectPointCloud(pcl::PointCloud<pcl::PointNormal> &_vertices);

		/// Return centroid 3D center of object.
		arma::colvec3 center();

		/// Get the min and max point that encloses boundaries of the object.
		virtual void minMax(arma::colvec3 &_min, arma::colvec&_max);

		/// Intersect a directed ray with the mesh.
		/// \param _initpoint: point that defines the ray position.
		/// \param _dir: vector that defines the direction of the ray.
		arma::colvec3 intersect(arma::colvec3 _initPoint, arma::colvec3 _dir);

		/// Intersect a directed ray with the mesh.
		/// \param _initpoint: point that defines the ray position.
		/// \param _dir: vector that defines the direction of the ray.
		arma::mat intersectRay(arma::colvec3 _p1, arma::colvec3 _p2);

		/// Get closest point on the mesh surface from a given point.
		/// \param _p: given point.
		/// \return: closest point to given point.
        pcl::PointNormal closestVertex(const pcl::PointNormal &_p);

		/// Get closest point on the mesh surface from a given point.
		/// \param _p: given point.
		/// \return: id of closest point to given point.
        int closestVertexId(const pcl::PointNormal &_p);

		/// Get vertex point by id
		/// \param id
		/// \return point
        pcl::PointNormal vertex(int _id);

		/// Get mesh info
		template<typename PointType_>
        void vertices(pcl::PointCloud<PointType_> &_vertices){
			for(auto &p:mVertices){
				PointType_ point;
				point.x = p.x;
				point.y = p.y;
				point.z = p.z;
				_vertices.push_back(point);
			}
		}

	protected:
		virtual void moveObject();

	private:
        pcl::PointCloud<pcl::PointNormal> mVertices;
		arma::colvec3 mCentroid;
		bool mIsKdtreeInit = false;
        pcl::KdTreeFLANN<pcl::PointNormal> mKdTree;

		pcl::PointNormal mMinPoint, mMaxPoint;
	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_OBJECTS_OBJECTPOINTCLOUD_H_
