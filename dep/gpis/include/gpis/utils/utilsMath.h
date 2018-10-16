//
//
//
//
//

#include <armadillo>

namespace gpis {
	/// Compute the rotation matrix to tansform an unitary vector A into B.
	/// \param _a: vector to rotate.
	/// \param _b: target vector.
	arma::mat rotUnitaryVectors(const arma::vec &_a, const arma::vec &_b);


}