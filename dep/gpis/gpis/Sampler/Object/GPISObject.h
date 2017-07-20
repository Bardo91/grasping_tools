/*
 * Object.h
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      represents our belief about a geometric object in the world (like an apple)
 */

#ifndef OBJECT_H_
#define OBJECT_H_

#include "../Part/Part.h"
#include "../Prior/GPISPrior.h"

#define ARMA_NO_DEBUG
#include <armadillo>
#include <vector>

namespace gpis {
	class GPISObject {
	public:
		////////////////////////////////
		/// Constructor
		GPISObject(const std::vector<GPISPrior*> &gpisPriors, const std::vector<Part> &allParts, const arma::uvec &associatedParts, const int numDims);

		//////////////////
		// public methods
		/// Returns vector of associated parts' indices
		const arma::uvec assocParts() const { return assocParts_.head(numAssocParts_); }

		/// Returns the index of a single associated part
		unsigned int assocParts(const int index) const { return assocParts_[index]; }

		/// Return pose
		const arma::vec &parameters() const { return priorParameters_[priorType_]; }

		/// Compute association probability for one part
		double computeAssocLogProb(std::vector<Part> &parts, const unsigned partIndex, bool isSameAsPrevious);

		/// Find part's location in association vector
		int findPartInAssocParts(const int partIndex);

		/// Add new part to object
		void addAssociatedPart(const int partIndex);

		/// Remove associated part from object
		void permanentlyRemoveAssociatedPart(const int partIndex, const int locationInAssocParts, const CholUpdateData &cholUpdateData);

		/// Update the poses for all the priors
		void updatePriorParameters(const std::vector<Part> &parts);

		/// Update the prior type
		void updatePriorType(std::vector<Part> &parts);

		/// Get prior type
		int type() { return priorType_; }

	private:
		int numDims_;
		int numDimsPlusOne_;
		arma::uvec assocParts_;
		int numAssocParts_;
		int maxNumAssocParts_;
		std::vector<GPISPrior*> gpisPriors_;
		std::vector<arma::vec>  priorParameters_;
		std::vector<MatrixData> priorMatrixData_;
		int priorType_;
	};
}	//	namespace gpis


#endif /* OBJECT_H_ */
