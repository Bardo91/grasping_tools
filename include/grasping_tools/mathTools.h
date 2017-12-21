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
		
	arma::mat pointsInCircleNU(double _radius, arma::colvec3 _center, arma::colvec3 n, arma::colvec3 u, unsigned _nPoints);

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
}

#include "mathTools.inl"

#endif