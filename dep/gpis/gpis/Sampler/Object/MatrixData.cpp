#include "MatrixData.h"
#include <algorithm>

#include <Utils/utilsProb.h>

using namespace std;
using namespace arma;

namespace gpis {
	MatrixData::MatrixData(const int numDimsPlusOne, const int numParts, int numAssocParts, arma::mat covMat)
		: numDimsPlusOne_(numDimsPlusOne),
		numDims_(numDimsPlusOne - 1),
		numParts_(numParts),
		numAssocParts_(numAssocParts),
		size_(numAssocParts_ * numDimsPlusOne_),
		fullSize_(size_),
		auxPMat_(size_, numParts_ * numDimsPlusOne_),
		auxPMatRem_(size_, numParts_ * numDimsPlusOne_),
		auxQVec_(size_),
		auxQVecRem_(size_),
		lCholDec_(size_, size_),
		firstDifferentIndexPMat_(numParts_),
		firstDifferentIndexQVec_(0),
		firstDifferentTrIndexQVec_(0)
	{
		auxPMat_.zeros();
		auxPMatRem_.zeros();
		auxQVec_.zeros();
		auxQVecRem_.zeros();
		lCholDec_ = symmatl(chol(covMat, "lower"));
		firstDifferentIndexPMat_.zeros();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	double MatrixData::gpisLogLikelihood(const int &partIndex, const mat &covMatAssocAndCandidate, const vec &fPlusVecAssocAndCandidate) {
		// Update the auxiliary matrix
		updateAuxPMat(partIndex, covMatAssocAndCandidate.head_rows(size_));
		// Update the auxiliary vector
		updateAuxQVec(fPlusVecAssocAndCandidate.head_rows(size_));

		span assocPartsDerInds = span(0, size_ - 1);
		span partDerInds = span(partIndex * numDimsPlusOne_, (partIndex + 1) * numDimsPlusOne_ - 1);

		// compute the association probability
		mat fStarCovMatAndMean = trans(auxPMat_(assocPartsDerInds, partDerInds)) *
			join_horiz(auxPMat_(assocPartsDerInds, partDerInds), auxQVec_(assocPartsDerInds));

		// return fStarCovMatAndMean;	
		return computeGaussianLogPdf(	fPlusVecAssocAndCandidate.tail_rows(numDimsPlusOne_),
										fStarCovMatAndMean.tail_cols(1),
										covMatAssocAndCandidate.tail_rows(numDimsPlusOne_) - fStarCovMatAndMean.head_cols(numDimsPlusOne_));;

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	double MatrixData::gpisLogLikelihood(const int &partIndex,
										 const mat &covMatAssocAndCandidate,
										 const vec &fPlusVecAssocAndCandidate,
										 const int locationInAssocParts,
										 CholUpdateData &cholUpdateData) {

		span partDerInds = span(partIndex * numDimsPlusOne_, (partIndex + 1) * numDimsPlusOne_ - 1);

		choleskyRemovePartialUpdate(locationInAssocParts, cholUpdateData);
		// Update the auxiliary matrix
		updateAuxPMatRem(partIndex, covMatAssocAndCandidate.head_rows(size_), cholUpdateData, locationInAssocParts);
		// Update the auxiliary vector
		updateAuxQVecRem(fPlusVecAssocAndCandidate.head_rows(size_), cholUpdateData, locationInAssocParts);

		span persistentDerInds = span(0, locationInAssocParts * numDimsPlusOne_ - 1);
		span nonPersistentDerInds = span(locationInAssocParts * numDimsPlusOne_, size_ - numDimsPlusOne_ - 1);

		// auxPMat is now put together from auxPMat and auxPMatRem
		mat fStarCovMatAndMean;
		fStarCovMatAndMean = trans(auxPMat_(persistentDerInds, partDerInds)) *
				join_horiz(auxPMat_(persistentDerInds, partDerInds),
						   auxQVec_.subvec(persistentDerInds)) +
				trans(auxPMatRem_(nonPersistentDerInds, partDerInds)) *
				join_horiz(auxPMatRem_(nonPersistentDerInds, partDerInds),
						   auxQVecRem_(nonPersistentDerInds));

		return computeGaussianLogPdf(fPlusVecAssocAndCandidate.tail_rows(numDimsPlusOne_),
									 fStarCovMatAndMean.tail_cols(1),
									 covMatAssocAndCandidate.tail_rows(numDimsPlusOne_) - fStarCovMatAndMean.head_cols(numDimsPlusOne_));
	}



	/////////////////////////////////////////////////////////////////////////////////////////////////////
	void MatrixData::removePartPermanently(const int partIndex, const int locationInAssocParts, const CholUpdateData &cholUpdateData) {
		numAssocParts_--;
		size_ = numAssocParts_ * numDimsPlusOne_;

        span headDerInds ;
        if (locationInAssocParts > 0){
            headDerInds = span(0, locationInAssocParts * numDimsPlusOne_ - 1);
        }
        span tailDerInds = span(locationInAssocParts * numDimsPlusOne_, size_ - 1);
		span tailDerIndsShifted = span((locationInAssocParts + 1) * numDimsPlusOne_, size_ + numDimsPlusOne_ - 1);
        int numTail = size_ - locationInAssocParts * numDimsPlusOne_;
		span cholUpdateInds = span(0, numTail - 1);
		span partDerInds = span(partIndex * numDimsPlusOne_, (partIndex + 1) * numDimsPlusOne_ - 1);
		
		if (numTail > 0) {
			auxPMat_(tailDerInds, partDerInds) =
				auxPMatRem_(tailDerInds, partDerInds);
			auxQVec_(tailDerInds) =
				auxQVecRem_(tailDerInds);

            if (locationInAssocParts > 0){
                lCholDec_(tailDerInds, headDerInds) = lCholDec_(tailDerIndsShifted, headDerInds);
            }
			lCholDec_(tailDerInds, tailDerInds) = cholUpdateData.updLChol_.slice(0)(cholUpdateInds, cholUpdateInds);
		}
		firstDifferentIndexQVec_ = numAssocParts_;
		firstDifferentIndexPMat_(partIndex) = numAssocParts_;

		for (int n = 0; (n < numParts_); n++) {
			if (n != partIndex) {
				firstDifferentIndexPMat_(n) = std::min(locationInAssocParts, static_cast<int>(firstDifferentIndexPMat_(n)));
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	void MatrixData::removePartPermanently(const int partIndex, const int locationInAssocParts){
        choleskyRemoveFullUpdate(locationInAssocParts);
        numAssocParts_--;
        size_ = numAssocParts_ * numDimsPlusOne_;

		firstDifferentIndexQVec_ = locationInAssocParts;
		for (int n = 0; (n < numParts_); n++) {
			if (n != partIndex) {
				firstDifferentIndexPMat_(n) =  std::min(locationInAssocParts, static_cast<int>(firstDifferentIndexPMat_(n)));
			}
		}
	}

    void MatrixData::addAssociatedPart(const mat &covMatCandidate, unsigned partIndex) {
		numAssocParts_++;
		size_ = numAssocParts_ * numDimsPlusOne_;
		if (size_ > fullSize_) {
			fullSize_ = size_;
			auxPMat_.resize(fullSize_, numParts_ * numDimsPlusOne_);
			auxPMatRem_.resize(fullSize_, numParts_ * numDimsPlusOne_);
			auxQVec_.resize(fullSize_);
			auxQVecRem_.resize(fullSize_);
			lCholDec_.resize(fullSize_, fullSize_);
		}
		// recompute the Cholesky decomposition for covariance matrix of the object (now including the new part)
        choleskyAppend(covMatCandidate, partIndex);

	}

    void MatrixData::addAssociatedPart(const mat &covMatAssocAndCandidate) {
        numAssocParts_++;
        size_ = numAssocParts_ * numDimsPlusOne_;
        if (size_ > fullSize_) {
            fullSize_ = size_;
            auxPMat_.resize(fullSize_, numParts_ * numDimsPlusOne_);
            auxPMatRem_.resize(fullSize_, numParts_ * numDimsPlusOne_);
            auxQVec_.resize(fullSize_);
            auxQVecRem_.resize(fullSize_);
            lCholDec_.resize(fullSize_, fullSize_);
        }
        // recompute the Cholesky decomposition for covariance matrix of the object (now including the new part)
        choleskyAppend(covMatAssocAndCandidate);

    }


    void MatrixData::choleskyAppend(const arma::mat &covMatCandidate, unsigned partIndex) {
		span oldIndices = span(0, size_ - numDimsPlusOne_ - 1);
		span appIndices = span(size_ - numDimsPlusOne_, size_ - 1);
        lCholDec_(appIndices, oldIndices) = auxPMat_.submat(0,
                                                            partIndex * numDimsPlusOne_,
                                                            size_ - numDimsPlusOne_ - 1,
                                                            (partIndex + 1) * numDimsPlusOne_ - 1).t();
        lCholDec_(appIndices, appIndices) = chol(covMatCandidate - lCholDec_(appIndices, oldIndices) *
                                                 lCholDec_(appIndices, oldIndices).t(), "lower");
	}

    void MatrixData::choleskyAppend(const arma::mat &covMatAssocAndCandidate) {
        span oldIndices = span(0, size_ - numDimsPlusOne_ - 1);
        span appIndices = span(size_ - numDimsPlusOne_, size_ - 1);
        lCholDec_(appIndices, oldIndices) = solve(trimatl(lCholDec_(oldIndices, oldIndices)), covMatAssocAndCandidate(oldIndices, span())).t();
        lCholDec_(appIndices, appIndices) = chol(covMatAssocAndCandidate(appIndices, span()) - lCholDec_(appIndices, oldIndices) * lCholDec_(appIndices, oldIndices).t(), "lower");
    }



	void MatrixData::choleskyRemovePartialUpdate(const int locationOfPartToRemove, CholUpdateData &cholUpdateData) {
        if (locationOfPartToRemove == numAssocParts_ - 1) {
			// The part to be removed is the last of all parts associated to this object.
			// This means that the cholesky decomposition does not change (except removing the last line)
			return;
		}

		int nBRC = numAssocParts_ - locationOfPartToRemove - 1;
		int nBRCDer = nBRC * numDimsPlusOne_;

		if (cholUpdateData.firstDifferentTrIndex_ == nBRC) {
			// We are done and can use the old Cholesky update
			return;
		}


		// Check if the cholUpdateData-matrices need to be extended
		if (nBRCDer > cholUpdateData.currentSize_) {
			// They do, so extend to new number of trailing derInds
			cholUpdateData.updLChol_.resize(nBRCDer, nBRCDer, numDimsPlusOne_);
			cholUpdateData.xVec_.resize(nBRCDer, numDimsPlusOne_);
			cholUpdateData.sinVec_.resize(nBRCDer, numDimsPlusOne_);
			cholUpdateData.cosVec_.resize(nBRCDer, numDimsPlusOne_);
			cholUpdateData.currentSize_ = nBRCDer; //currentSize_ is now always >= nBRCDer
		}


		// some vector spans we will need more often
		span remainingBRCDerInds = span(cholUpdateData.firstDifferentTrIndex_*numDimsPlusOne_, nBRCDer - 1);
		span allBRCDerInds = span(0, nBRCDer - 1);
		span remainingDerInds = span((locationOfPartToRemove + cholUpdateData.firstDifferentTrIndex_ + 1) * numDimsPlusOne_, size_ - 1);
		span allDerInds = span((locationOfPartToRemove + 1) * numDimsPlusOne_, size_ - 1);



		// The trailing derIndices have changed since the last computation and some of the
		// Cholesky update needs to be recomputed

		// The remaining indices are filled up with the original CholDec and x-Vector: Start with the last one (will be the first to be evaluated)

		cholUpdateData.updLChol_.slice(numDims_)(remainingBRCDerInds, allBRCDerInds) =
			lCholDec_(remainingDerInds, allDerInds);
		cholUpdateData.xVec_(remainingBRCDerInds, numDims_) =
			lCholDec_(remainingDerInds, locationOfPartToRemove * numDimsPlusOne_ + numDims_);

		// cache to avoid repeated arma element access
		double lCholDec_kk = 0.0;
		double lCholDec_kj = 0.0;

		double t = 0.0;
		double norm = 0.0;
		for (int d = numDims_; d >= 0; --d) // needs to be done for value and all gradient components. Start with last index
		{
			if (d < numDims_)
			{
				// The remaining indices are filled up with the previous CholDec and x-Vector
				cholUpdateData.updLChol_.slice(d)(remainingBRCDerInds, allBRCDerInds) =
					cholUpdateData.updLChol_.slice(d + 1)(remainingBRCDerInds, allBRCDerInds);
				cholUpdateData.xVec_(remainingBRCDerInds, d) = lCholDec_(remainingDerInds, locationOfPartToRemove * numDimsPlusOne_ + d);
			}

			for (int k = (cholUpdateData.firstDifferentTrIndex_ * numDimsPlusOne_); k < nBRCDer; ++k) // start with first non-persistent index
			{
				// (Part of) the actual cholesky update
				for (int j = 0; j < k; ++j)
				{
					lCholDec_kj = cholUpdateData.updLChol_(k, j, d);
					t = cholUpdateData.cosVec_(j, d) * lCholDec_kj
						+ cholUpdateData.sinVec_(j, d) * cholUpdateData.xVec_(k, d);
					cholUpdateData.xVec_(k, d) =
						-cholUpdateData.sinVec_(j, d) * lCholDec_kj
						+ cholUpdateData.cosVec_(j, d) * cholUpdateData.xVec_(k, d);
					cholUpdateData.updLChol_(k, j, d) = t;
					//						cholUpdateData.updLChol_(j,k,d) = t;
				}
				lCholDec_kk = cholUpdateData.updLChol_(k, k, d);
				norm = sqrt(pow(lCholDec_kk, 2) + pow(cholUpdateData.xVec_(k, d), 2));
				cholUpdateData.cosVec_(k, d) = lCholDec_kk / norm;
				cholUpdateData.sinVec_(k, d) = cholUpdateData.xVec_(k, d) / norm;
				cholUpdateData.updLChol_(k, k, d) = norm;
			}
		}

		cholUpdateData.firstDifferentTrIndexOld_ = cholUpdateData.firstDifferentTrIndex_;
		cholUpdateData.firstDifferentTrIndex_ = nBRC;
	}


	// 666 CHECK METHOD
	void MatrixData::choleskyRemoveFullUpdate(const int locationOfPartToRemove) { 

		if (locationOfPartToRemove == numAssocParts_ - 1) {
			// The part to be removed is the last of all parts associated to this object.
			// This means that the cholesky decomposition does not change
			// The cholesky update data also doesn't change. Leave everything as it is and return.
		}
		else {

			int nBRC = (numAssocParts_ - 1) - locationOfPartToRemove;
			int nBRCDer = nBRC * numDimsPlusOne_;

			double t = 0.0;
			double norm = 0.0;
            vec x(nBRCDer);
            vec c(nBRCDer);
            vec s(nBRCDer);

            unsigned firstBRCDerInd = (locationOfPartToRemove + 1) * numDimsPlusOne_;
            for (int d = numDims_; d >= 0; --d) { // needs to be done for value and all gradient components. Start with last index
				unsigned indToRemove = locationOfPartToRemove * numDimsPlusOne_ + d;
                x = lCholDec_(span(firstBRCDerInd, firstBRCDerInd + nBRCDer - 1), indToRemove);
				for (int k = 0; k < nBRCDer; ++k) {
					// The cholesky update
					for (int j = 0; j < k; ++j) {
                        t = c(j) * lCholDec_(firstBRCDerInd + k, firstBRCDerInd + j) + s(j) * x(k);
                        x(k) = -s(j) * lCholDec_(firstBRCDerInd + k, firstBRCDerInd + j) + c(j) * x(k);
                        lCholDec_(firstBRCDerInd + k, firstBRCDerInd + j) = t;
					}

                    norm = sqrt(pow(lCholDec_(firstBRCDerInd + k, firstBRCDerInd + k), 2) + pow(x(k), 2));
                    c(k) = lCholDec_(firstBRCDerInd + k, firstBRCDerInd + k) / norm;
                    s(k) = x(k) / norm;
                    lCholDec_(firstBRCDerInd + k, firstBRCDerInd + k) = norm;
				}
			}
            span tailDerInds = span(firstBRCDerInd - numDimsPlusOne_, size_ - numDimsPlusOne_ - 1);
            span tailDerIndsShifted = span(firstBRCDerInd, size_ - 1);
            span headDerInds = span(0, firstBRCDerInd - numDimsPlusOne_);
            lCholDec_(tailDerInds, headDerInds) = lCholDec_(tailDerIndsShifted, headDerInds);
            lCholDec_(tailDerInds, tailDerInds) = lCholDec_(tailDerIndsShifted, tailDerIndsShifted);
		}

	}



	void MatrixData::updateAuxPMat(const int partIndex, const mat &covMatCols) {
		// This auxiliary matrix stores inv(fullLCholDec(derIndices,derIndices)) * fullCovMat(derIndices,partDerIndices) for this object
		int firstPartDerInd = partIndex * numDimsPlusOne_;

		int firstDifferentDerIndex = firstDifferentIndexPMat_(partIndex) * numDimsPlusOne_;
		// iterators are MUCH faster than element access for vectors
		double sum;
		int numLHS = size_;
		for (int iRHS = 0; iRHS < numDimsPlusOne_; iRHS++) {
			for (int iLHS = firstDifferentDerIndex; iLHS < numLHS; iLHS++) {
				sum = 0;
				for (int jLHS = 0; jLHS < iLHS; jLHS++) {
					sum += lCholDec_(iLHS, jLHS) * auxPMat_(jLHS, firstPartDerInd + iRHS);
				}
				auxPMat_(iLHS, firstPartDerInd + iRHS) = (covMatCols(iLHS, iRHS) - sum) / lCholDec_(iLHS, iLHS);
			}
		}
		firstDifferentIndexPMat_(partIndex) = numAssocParts_;

	}

	void MatrixData::updateAuxPMatRem(const int partIndex, const mat &covMatCols, const CholUpdateData &cholUpdateData, const int locationInAssocParts) {

		int firstPartDerInd = partIndex * numDimsPlusOne_;

		int firstDifferentDerIndex = firstDifferentIndexPMat_(partIndex) * numDimsPlusOne_;
		int locationInAssocPartsDerIndex = locationInAssocParts * numDimsPlusOne_;
		int firstDifferentTrDerIndex = cholUpdateData.firstDifferentTrIndexOld_ * numDimsPlusOne_;
		double sum;

		int numLHS = size_ - numDimsPlusOne_;


		for (int iRHS = 0; iRHS < numDimsPlusOne_; iRHS++) {
			for (int iLHS = firstDifferentDerIndex; iLHS < locationInAssocPartsDerIndex; iLHS++) {
				sum = 0;
				for (int jLHS = 0; jLHS < iLHS; jLHS++) {
					sum += lCholDec_(iLHS, jLHS) * auxPMat_(jLHS, firstPartDerInd + iRHS);
				}
				auxPMat_(iLHS, firstPartDerInd + iRHS) =
					(covMatCols(iLHS, iRHS) - sum) /
					lCholDec_(iLHS, iLHS);
			}
			for (int iLHS = locationInAssocPartsDerIndex + firstDifferentTrDerIndex; iLHS < numLHS; iLHS++) {
				sum = 0;
				for (int jLHS = 0; jLHS < locationInAssocPartsDerIndex; jLHS++) {
					sum += lCholDec_(iLHS + numDimsPlusOne_, jLHS) * auxPMat_(jLHS, firstPartDerInd + iRHS);
				}
				for (int jLHS = locationInAssocPartsDerIndex; jLHS < iLHS; jLHS++) {
					sum += cholUpdateData.updLChol_.slice(0)
						(iLHS - locationInAssocPartsDerIndex,
							jLHS - locationInAssocPartsDerIndex) *
						auxPMatRem_(jLHS, firstPartDerInd + iRHS);
				}
				auxPMatRem_(iLHS, firstPartDerInd + iRHS) =
					(covMatCols(iLHS + numDimsPlusOne_, iRHS) - sum) /
					cholUpdateData.updLChol_.slice(0)
					(iLHS - locationInAssocPartsDerIndex,
						iLHS - locationInAssocPartsDerIndex);
			}
		}
		firstDifferentIndexPMat_(partIndex) = std::max(locationInAssocParts, int(firstDifferentIndexPMat_(partIndex)));
	}

	void MatrixData::updateAuxQVec(const vec &auxFPlusVecAssoc) {
		int firstDifferentDerIndex = firstDifferentIndexQVec_ * numDimsPlusOne_;

		double sum;

		int numLHS = size_;

		for (int i = firstDifferentDerIndex; i < numLHS; i++) {
			sum = 0;
			for (int j = 0; j < i; j++) {
				sum += lCholDec_(i, j) * auxQVec_(j);
			}
			auxQVec_(i) = (auxFPlusVecAssoc(i) - sum) / lCholDec_(i, i);
		}
		firstDifferentIndexQVec_ = numAssocParts_;

	}

	void MatrixData::updateAuxQVecRem(const vec &auxFPlusVecAssoc,
		const CholUpdateData &cholUpdateData,
		const int locationInAssocParts) {

		int firstDifferentDerIndex = firstDifferentIndexQVec_ * numDimsPlusOne_;
		int locationInAssocPartsDerIndex = locationInAssocParts * numDimsPlusOne_;
		int firstDifferentTrDerIndex = firstDifferentTrIndexQVec_ * numDimsPlusOne_;

		double sum;

		int num = size_ - numDimsPlusOne_;

		for (int i = firstDifferentDerIndex; i < locationInAssocPartsDerIndex; i++) {
			sum = 0;
			for (int j = 0; j < i; j++) {
				sum += lCholDec_(i, j) * auxQVec_(j);
			}
			auxQVec_(i) =
				(auxFPlusVecAssoc(i) - sum) /
				lCholDec_(i, i);
		}
		for (int i = locationInAssocPartsDerIndex + firstDifferentTrDerIndex; i < num; i++) {
			sum = 0;
			for (int j = 0; j < locationInAssocPartsDerIndex; j++) {
				sum += lCholDec_(i + numDimsPlusOne_, j) * auxQVec_(j);
			}
			for (int j = locationInAssocPartsDerIndex; j < i; j++) {
				sum += cholUpdateData.updLChol_.slice(0)
					(i - locationInAssocPartsDerIndex,
						j - locationInAssocPartsDerIndex)
					* auxQVecRem_(j);
			}

			auxQVecRem_(i) =
				(auxFPlusVecAssoc(i + numDimsPlusOne_) - sum) /
				cholUpdateData.updLChol_.slice(0) (i - locationInAssocPartsDerIndex, i - locationInAssocPartsDerIndex);
		}
		firstDifferentIndexQVec_ = std::max(locationInAssocParts, int(firstDifferentIndexQVec_));
		firstDifferentTrIndexQVec_ = cholUpdateData.firstDifferentTrIndexOld_;
	}

}	//	namespace gpis
