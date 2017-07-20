#include "Kernel.h"

#include<cassert>

using namespace arma;
using namespace std;

namespace gpis {
	//-----------------------------------------------------------------------------------------------------------------
	mat Kernel::operator()(const mat & points1, bool computeGradients) const {
		// Check that both vectors of points have the same dimension.
		if (computeGradients) { // If gradients are needed.
			const unsigned dimensions = points1.n_rows;
			const unsigned dimPlusOne = dimensions + 1;
			mat covarianceMatrix(points1.n_cols*dimPlusOne, points1.n_cols*dimPlusOne);

			for (unsigned i = 0; i < points1.n_cols; i++) {
				for (unsigned j = i; j < points1.n_cols; j++) {
					auto iPoint = points1.col(i);
					auto jPoint = points1.col(j);
					covarianceMatrix(i*dimPlusOne, j*dimPlusOne) = covarianceValVal(iPoint, jPoint, (i == j));
					auto valder = covarianceValDer(iPoint, jPoint);
					unsigned igs = i*dimPlusOne + 1; //gradient i-index start
					unsigned ige = igs + dimensions - 1; //gradient i-index end
					unsigned jgs = j*dimPlusOne + 1; //gradient j-index start
					unsigned jge = jgs + dimensions - 1; //gradient i-index end

					covarianceMatrix.submat(i*dimPlusOne, jgs, i*dimPlusOne, jge) = trans(valder);
					covarianceMatrix.submat(igs, j*dimPlusOne, ige, j*dimPlusOne) = -valder;
					auto derder = covarianceDerDer(iPoint, jPoint, (i == j));
					covarianceMatrix.submat(igs, jgs, ige, jge) = derder;
					// fill up lower triangular matrix exploiting symmetry
					if (i != j) {
						covarianceMatrix(j*dimPlusOne, i*dimPlusOne) = covarianceMatrix(i*dimPlusOne, j*dimPlusOne);
						covarianceMatrix.submat(jgs, i*dimPlusOne, jge, i*dimPlusOne) = valder;
						covarianceMatrix.submat(j*dimPlusOne, igs, j*dimPlusOne, ige) = trans(-valder);
						covarianceMatrix.submat(jgs, igs, jge, ige) = trans(derder);
					}
				}
			}
			return covarianceMatrix;
		}
		else {	// If don't need gradients.
			mat covarianceMatrix(points1.n_cols, points1.n_cols);
			for (unsigned i = 0; i < points1.n_cols; i++) {
				for (unsigned j = 0; j < points1.n_cols; j++) {
					covarianceMatrix(i, j) = covarianceValVal(points1.col(i), points1.col(j), (i == j));
				}
			}
			return covarianceMatrix;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	mat Kernel::operator()(const mat & points1, const mat & points2, bool computeGradients) const {
		// Check that both vectors of points have the same dimension.
		assert(points1.n_rows == points2.n_rows);
		if (computeGradients) { // If gradients are needed.
			const unsigned dimensions = points1.n_rows;
			const unsigned dimPlusOne = dimensions + 1;
			mat covarianceMatrix(points1.n_cols*dimPlusOne, points2.n_cols*dimPlusOne);

			for (unsigned i = 0; i < points1.n_cols; i++) {
				for (unsigned j = 0; j < points2.n_cols; j++) {
					auto iPoint = points1.col(i);
					auto jPoint = points2.col(j);
					covarianceMatrix(i*dimPlusOne, j*dimPlusOne) = covarianceValVal(iPoint, jPoint, (i == j));
					auto valder = covarianceValDer(iPoint, jPoint);
					unsigned igs = i*dimPlusOne + 1; //gradient i-index start
					unsigned ige = igs + dimensions - 1; //gradient i-index end
					unsigned jgs = j*dimPlusOne + 1; //gradient j-index start
					unsigned jge = jgs + dimensions - 1; //gradient i-index end

					covarianceMatrix.submat(i*dimPlusOne, jgs, i*dimPlusOne, jge) = trans(valder);
					covarianceMatrix.submat(igs, j*dimPlusOne, ige, j*dimPlusOne) = -valder;
					auto derder = covarianceDerDer(iPoint, jPoint, (i == j));
					covarianceMatrix.submat(igs, jgs, ige, jge) = derder;\
				}
			}
			return covarianceMatrix;
		}
		else {	// If don't need gradients.
			mat covarianceMatrix(points1.n_cols, points2.n_cols);
			for (unsigned i = 0; i < points1.n_cols; i++) {
				for (unsigned j = 0; j < points2.n_cols; j++) {
					covarianceMatrix(i, j) = covarianceValVal(points1.col(i), points2.col(j), (i == j));
				}
			}
			return covarianceMatrix;
		}
	}

}	//	namespace gpis
