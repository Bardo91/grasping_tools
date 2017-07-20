/*
 * CholUpdateData.cpp
 *
 *  Created on: 26/01/2016
 *      Author: wolfram
 *
 *      stores Choleski update data (for a single part)
 */

#include "CholUpdateData.h"

using namespace arma;
namespace gpis {
	CholUpdateData::CholUpdateData(const int numDims)
		: currentSize_(0),
		updLChol_(cube(currentSize_, currentSize_, numDims + 1)),
		cosVec_(mat(currentSize_, numDims + 1)),
		sinVec_(mat(currentSize_, numDims + 1)),
		xVec_(mat(currentSize_, numDims + 1)),
		firstDifferentTrIndex_(0),
		firstDifferentTrIndexOld_(0)
	{
		// empty
	}
}