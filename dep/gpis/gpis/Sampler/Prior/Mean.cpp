#include "Mean.h"
#include <cassert>

using namespace std;
using namespace arma;

namespace gpis{
	mat Mean::operator()(const mat & points, const arma::vec &_parameters, bool computeGradients) const {
		if (!computeGradients) { // If don't need gradients.
			mat meanMat(1, points.n_cols);
			for (unsigned i = 0; i < points.n_cols; i++) {
				meanMat.col(i).head(1) = meanFun(points.col(i), _parameters);
			}
			return meanMat;
		}
		else {	// If gradients are needed.
			const unsigned dimensions = points.n_rows;
			mat meanMat(dimensions + 1, points.n_cols);
			for (unsigned i = 0; i < points.n_cols; i++) {
				meanMat.col(i).head(1) = meanFun(points.col(i), _parameters);
				meanMat.col(i).tail(dimensions) = meanGradsFun(points.col(i), _parameters);
			}
			return meanMat;
		}
	}

}	//	namespace gpis
