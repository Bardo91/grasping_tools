#ifndef MATRIXDATA_H
#define MATRIXDATA_H


#define ARMA_NO_DEBUG
#include <armadillo>
#include <vector>

#include <gpis/sampler/part/CholUpdateData.h>

namespace gpis {
	/// MatrixData does all the GP related matrix operations
	class MatrixData {
	public:
		/// Constructor
		MatrixData(const int numDimsPlusOne, const int numParts, int numAssocParts, arma::mat covMat);

		// public methods
		/// Get the current cholesky decomposition matrix.
		const arma::mat lCholDec() const { return lCholDec_(arma::span(0, size_ - 1), arma::span(0, size_ - 1)); };

		/// Compute association probability for queried part with the given data, part does not currently belong to object
		double gpisLogLikelihood(const int &partIndex,
									const arma::mat &covMatPartAndAssoc,
									const arma::vec &fPlusVecPartAndAssoc);

		/// Compute association probability for queried part with the given data, part currently belongs to object
		double gpisLogLikelihood(const int &partIndex,
									const arma::mat &covMatPartAndAssoc,
									const arma::vec &fPlusVecPartAndAssoc,
									const int locationInAssocParts,
									CholUpdateData &cholUpdateData);

		/// Remove an associated part from the object when prior is the current priorType_ of the object
		void removePartPermanently(const int partIndex,
			const int locationInAssocParts,
			const CholUpdateData &cholUpdateData);

		/// Remove an associated part from the object when prior is not the current priorType_ of the object
		void removePartPermanently(const int partIndex,
			const int locationInAssocParts);

        /// Associate a new part to the object.
        void addAssociatedPart(const arma::mat &covMatCandidate, unsigned partIndex);

        /// Associate a new part to the object.
        void addAssociatedPart(const arma::mat &covMatAssocAndCandidate);

		/// Reset the indeces.
		void resetFirstDifferentIndexQVec() { firstDifferentIndexQVec_ = 0; firstDifferentTrIndexQVec_ = 0; }

//        /// Update fPlusVector
//        void fPlusVec(const arma::vec &fPlusVec) {fPlusVec_ = fPlusVec;}

	private:
		int numDimsPlusOne_; //TODO: Should be const, but problem with copy/move constructor
		int numDims_; //TODO: Should be const, but problem with copy/move constructor
		int numParts_; //TODO: Should be const, but problem with copy/move constructor
		int numAssocParts_;
		int size_;
		int fullSize_;
		arma::mat auxPMat_;
		arma::mat auxPMatRem_;
		arma::vec auxQVec_;
		arma::vec auxQVecRem_;
		arma::mat lCholDec_;
		arma::uvec firstDifferentIndexPMat_;
		unsigned int firstDifferentIndexQVec_;
		unsigned int firstDifferentTrIndexQVec_;

        void choleskyAppend(const arma::mat &covMatCandidate, unsigned partIndex);
        void choleskyAppend(const arma::mat &covMatAssocAndCandidate);
		void choleskyRemoveFullUpdate(const int locationOfPartToRemove);
		void choleskyRemovePartialUpdate(const int locationOfPartToRemove, CholUpdateData &cholUpdateData);
		void updateAuxPMat(const int partIndex, const arma::mat &covMatCols);
		void updateAuxPMatRem(const int partIndex, const arma::mat &covMatCols, const CholUpdateData &cholUpdateData, const int locationInAssocParts);
		void updateAuxQVec(const arma::vec &auxFPlusVecAssoc);
		void updateAuxQVecRem(const arma::vec &auxFPlusVecAssoc, const CholUpdateData &cholUpdateData, const int locationInAssocParts);

	};
}	//	namespace gpis
#endif // MATRIXDATA_H
