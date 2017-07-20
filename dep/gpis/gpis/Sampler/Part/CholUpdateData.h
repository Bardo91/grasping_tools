/*
 * CholUpdateData.h
 *
 *  Created on: 26/01/2016
 *      Author: wolfram
 *
 *      stores Choleski update data (for a single part)
 */

#ifndef CHOLUPDATEDATA_H
#define CHOLUPDATEDATA_H

#define ARMA_NO_DEBUG
#include <armadillo>

namespace gpis {

	class CholUpdateData
	{
	public:
		//////////////////////////////////
		// constructors and destructors
		CholUpdateData(const int numDims);

		// These need to be accessed and updated from composition.cpp, so need to be public?
		int currentSize_;
		arma::cube updLChol_;
		arma::mat cosVec_;
		arma::mat sinVec_;
		arma::mat xVec_;
		int firstDifferentTrIndex_;
		int firstDifferentTrIndexOld_;

	private:
	};
}	// namespace gpis

#endif // CHOLUPDATEDATA_H
