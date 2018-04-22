///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_QUALITYMETRICS_MATHTOOLS_H_
#define GPISGRASPING_QUALITYMETRICS_MATHTOOLS_H_

#include <armadillo>
#include <vector>

#include <pcl/conversions.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>

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

#include <cassert>

#ifndef M_PI
	#define M_PI 3.14159265359
#endif

namespace grasping_tools {
	//--------------------------------------------------------------------------------------------------------------------
	/// Generate all possible combinations between the elements in the vector with varying size of combinations,
	/// without repetitions and without order. Example:
	///	\param: vector with numbers to combine.
	/// \code
	///		arma::vec v = {1, 2, 3};
	///		auto res = generateCombinationVaryingSize(v);
	///		// res variable contains:
	///		//	{ {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3} }.
	///	\endcode
	template<typename _vectorType>
	std::vector<_vectorType> generateCombinationsVaryingSize(const _vectorType &_vec);

	//--------------------------------------------------------------------------------------------------------------------
	/// Generate all possible combinations between the elements in the vector without repetitions and without order. Example:
	///	\param: vector with numbers to combine.
	/// \code
	///		arma::vec v = {1, 2, 3};
	///		auto res = generateCombinationVaryingSize(v);
	///		// res variable contains:
	///		//	{ {1, 2, 3},  {2, 3, 4}, {1, 2, 4}, {1, 3, 4} }.
	///	\endcode
	template<typename _vectorType>
	std::vector<_vectorType> generateCombinationsWithoutRepetitions(const _vectorType &_vec);


	//--------------------------------------------------------------------------------------------------------------------
	/// Generate points uniformly on the surface on an unrotated ellipsoid with semi-axis _a, _b, _c.
	///	\param _a: semiaxis A.
	///	\param _b: semiaxis B.
	///	\param _c: semiaxis C.
	arma::mat samplePointsEllipsoidEquidistributed(double _a, double _b, double _c, unsigned _nPoints);

	//--------------------------------------------------------------------------------------------------------------------
	/// Generate points on a circle on 3D
	/// \param _radius:
	/// \param _center:
	/// \param _azimut:
	/// \param _zenit:
	/// \param _nPoint:
	/// \return matrix that holds all the points on each row.
	arma::mat pointsInCircle(double _radius, arma::colvec3 _center, double _azimut, double _zenit, unsigned _nPoints);
		
	//--------------------------------------------------------------------------------------------------------------------
	/// Generate points on a circle on 3D 666 TODO RENAME; SAME MANE DIFFERENT INPUTS!
	arma::mat pointsInCircleNU(double _radius, arma::colvec3 _center, arma::colvec3 n, arma::colvec3 u, unsigned _nPoints);

	//--------------------------------------------------------------------------------------------------------------------
	/// Generate points on a circle on 3D
	arma::mat pointsInCircle(double _radius, arma::colvec3 _center, arma::colvec3 _n, unsigned _nPoints);

	/// Compute the rotation matrix to tansform an unitary vector A into B.
	/// \param _a: vector to rotate.
	/// \param _b: target vector.
	arma::mat rotUnitaryVectors(const arma::vec &_a, const arma::vec &_b);


	/// Generate a random point on the surface of a sphere
	/// \param _sphereCenter
	/// \param _sphereRadius
	arma::colvec3 generateRandomPointSphere(const arma::colvec3 & _sphereCenter, const double _sphereRadius);

	/// Intersect a ray with a plane
	/// \param _p0: first point of ray
	/// \param _p1: second point of ray
	/// \param _v0: point on the plane
	/// \param _n: normal of plane
	/// \param _intersection: intersectedpoint
	/// \return true if intersects false if not
	bool intersectRayPlane(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _n, arma::colvec3 &_intersection);

	/// Intersect a ray with a triangle
	/// \param _p0: first point of ray
	/// \param _p1: second point of ray
	/// \param _v0: point on the plane
	/// \param _n: normal of plane
	/// \param _intersection: intersectedpoint
	/// \return true if intersects false if not
	bool intersectRayTriangle(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _v1, arma::colvec3 _v2, arma::colvec3 &_intersection);

	/// Intersect a segment with a plane
	/// \param _p0: first point of sgment
	/// \param _p1: second point of segment
	/// \param _v0: point on the triangle
	/// \param _v1: point on the triangle
	/// \param _v2: point on the triangle
	/// \param _intersection: intersectedpoint
	/// \return true if intersects false if not
	bool intersectSegmentPlane(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _n, arma::colvec3 &_intersection);

	/// Intersect a segment with a trianle
	/// \param _p0: first point of segment
	/// \param _p1: second point of segment
	/// \param _v0: point on the triangle
	/// \param _v1: point on the triangle
	/// \param _v2: point on the triangle
	/// \param _intersection: intersectedpoint
	/// \return true if intersects false if not
	bool intersectSegmentTriangle(arma::colvec3 _p0, arma::colvec3 _p1, arma::colvec3 _v0, arma::colvec3 _v1, arma::colvec3 _v2, arma::colvec3 &_intersection);

	/// Returns the convex hull of given mesh
	bool convexHull(pcl::PolygonMesh &_mesh, pcl::PolygonMesh &_out);

	/// Intersect an n-dimensional facet with a n-dimensional ray
	/// \param _facetVertices: a matrix containing the vertices of the n-dimensional facet by columns
	/// \param _ray: n-dimensional ray to intersect with the facet (Assumed to be in the origin)
	/// \param _cutPoint: cut point if intersected
	/// \return true if the ray intersects
	bool facetIntersection(arma::mat &_facetVertices, arma::colvec &_ray, arma::colvec &_cutPoint);

	/// Intersect an n-dimensional facet with a n-dimensional ray
	/// \param _facetVertices: a matrix containing the vertices of the n-dimensional facet by columns
	/// \param _rayOrigin: origin of n-dimensional ray to intersect with the facet.
	/// \param _rayEnd: end of n-dimensional ray to intersect with the facet.
	/// \param _cutPoint: cut point if intersected
	/// \return true if the ray intersects
	bool facetIntersection(arma::mat &_facetVertices, arma::colvec &_rayOrigin, arma::colvec &_rayEnd, arma::colvec &_cutPoint);

	/// Returns the convex hull of given point cloud
	template<typename PointType_>
	bool convexHull(pcl::PointCloud<PointType_> &_cloud, pcl::PolygonMesh &_mesh);
}

#include "mathTools.inl"

#endif