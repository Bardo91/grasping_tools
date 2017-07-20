/*
 * Object.cpp
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      represents our belief about a geometric object in the world (like an apple)
 */

#include "GPISObject.h"
#include <Utils/utilsMisc.h>
#include <Utils/utilsGPIS.h>
#include <Utils/utilsProb.h>

#include <iostream>

using namespace std;
using namespace arma;

namespace gpis {
	///////////////////////////////////////////////////////////////////////////////////
	// constructors and destructors
	GPISObject::GPISObject(const std::vector<GPISPrior*> &gpisPriors,
						   const std::vector<Part> &allParts,
						   const arma::uvec &associatedParts,
						   const int numDims)
										: numDims_(numDims),
										numDimsPlusOne_(numDims + 1),
										numAssocParts_(0),
										maxNumAssocParts_(0),
										gpisPriors_(gpisPriors),
										priorType_(0) {
		for (int i = 0; i < gpisPriors_.size(); i++) {

			auto covarianceMatrix = gpisPriors_[i]->fullCovariance();
			auto partDerInds = indsToDerInds(associatedParts, numDimsPlusOne_);
			priorMatrixData_.push_back(MatrixData(numDimsPlusOne_, allParts.size(), 0, covarianceMatrix.submat(partDerInds, partDerInds)));
			priorParameters_.push_back(arma::vec());
		}

		for (unsigned i = 0; i < associatedParts.n_elem; i++) {
			addAssociatedPart(associatedParts[i]);
		}

		for (int i = 0; i < gpisPriors_.size(); i++) {
			priorParameters_[i] = gpisPriors_[i]->sampleParameters(allParts, associatedParts, priorMatrixData_[i].lCholDec(), priorParameters_[i]);
		}
	}    

	///////////////////////////////////////////////////////////////////////////////////
	// public methods
	//


	//-----------------------------------------------------------------------------------------------------------------
	double GPISObject::computeAssocLogProb(std::vector<Part> &parts, const unsigned partIndex, bool isSameAsPrevious) {
		uvec partIndexUvec = {partIndex};
		uvec assocAndCandidateInds = join_vert(assocParts(), partIndexUvec);
		uvec assocAndCandidateDerInds = indsToDerInds(assocAndCandidateInds, numDimsPlusOne_);
		uvec partDerInds = indsToDerInds(partIndex, numDimsPlusOne_);


		double logPGeom = gpisPriors_[priorType_]->geometricLogLikelihood(parts, { partIndex }, priorParameters_[priorType_]);


		if (!isSameAsPrevious) {
			if (exp(logPGeom) < 0.000001) {
				return -std::numeric_limits<double>::max();
			}

			mat covMatAssocAndCandidate = gpisPriors_[priorType_]->fullCovariance()(assocAndCandidateDerInds, partDerInds);
			vec fPlusVecAssocAndCandidate = gpisPriors_[priorType_]->computeDataDeviation(parts, assocAndCandidateInds, priorParameters_[priorType_]);

			return logPGeom + priorMatrixData_[priorType_].gpisLogLikelihood(	partIndex,
																				covMatAssocAndCandidate,
																				fPlusVecAssocAndCandidate);
		}
		else {
			mat covMatAssocAndCandidate = gpisPriors_[priorType_]->fullCovariance()(assocAndCandidateDerInds, partDerInds);
			vec fPlusVecAssocAndCandidate = gpisPriors_[priorType_]->computeDataDeviation(parts, assocAndCandidateInds, priorParameters_[priorType_]);

			// Find location of part in object's part list
			int locationInAssocParts = findPartInAssocParts(partIndex);
			if (locationInAssocParts == -1)
				return -std::numeric_limits<double>::max();
			else
				return logPGeom + priorMatrixData_[priorType_].gpisLogLikelihood(	partIndex,
																				covMatAssocAndCandidate,
																				fPlusVecAssocAndCandidate,
																				locationInAssocParts,
																				parts[partIndex].cholUpdateData_);
		}
	}


	//-----------------------------------------------------------------------------------------------------------------
	int GPISObject::findPartInAssocParts(const int partIndex){
		uvec locationInAssocParts = find(assocParts() == partIndex);
		if (locationInAssocParts.is_empty()) {
			std::ostringstream errorString;
			errorString << "Object::findPartinAssocParts called for part that isn't associated. Part index: "
				<< partIndex;
			error(errorString);
			return -1;
		}
		return as_scalar(locationInAssocParts);
	}

	//-----------------------------------------------------------------------------------------------------------------
	void GPISObject::addAssociatedPart(const int partIndex) {
		numAssocParts_++;
		if (numAssocParts_ > maxNumAssocParts_) {
			assocParts_.resize(numAssocParts_);
			maxNumAssocParts_ = numAssocParts_;
		}
		assocParts_(numAssocParts_ - 1) = partIndex;

		// Update the matrixData for each prior
		uvec partIndexUvec(1);
		partIndexUvec.at(0) = partIndex;

//		uvec assocAndCandidateDerInds = indsToDerInds(assocParts(), numDimsPlusOne_);
		uvec candidateDerInds = indsToDerInds(partIndex, numDimsPlusOne_);

		for (unsigned i = 0; i < gpisPriors_.size(); i++){
			const mat &fullCov = gpisPriors_[i]->fullCovariance();
//            mat data = fullCov(assocAndCandidateDerInds, partDerInds);
			mat data = fullCov(candidateDerInds, candidateDerInds);
			priorMatrixData_[i].addAssociatedPart(data, partIndex);
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	void GPISObject::permanentlyRemoveAssociatedPart(const int partIndex, const int locationInAssocParts, const CholUpdateData &cholUpdateData) {
		numAssocParts_--;
		for (int n = locationInAssocParts; n < numAssocParts_; n++) {
			assocParts_(n) = assocParts_(n + 1);
		}
		// Update the matrixData for each prior
		for (unsigned i = 0; i < gpisPriors_.size(); i++){
			if (i == priorType_){
				priorMatrixData_[i].removePartPermanently(partIndex, locationInAssocParts, cholUpdateData);
			}
			else {
				priorMatrixData_[i].removePartPermanently(partIndex, locationInAssocParts);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	void GPISObject::updatePriorParameters(const std::vector<Part> &parts) {
		// Update the pose for each prior
		for (unsigned i = 0; i < gpisPriors_.size(); i++) {
			priorParameters_[i] = gpisPriors_[i]->sampleParameters(parts, assocParts(), priorMatrixData_[i].lCholDec(), priorParameters_[i]);
			priorMatrixData_[i].resetFirstDifferentIndexQVec();
		}
	}	
	
	//-----------------------------------------------------------------------------------------------------------------
	void GPISObject::updatePriorType(std::vector<Part> &parts) {
		if (numAssocParts_ == 0)
			return;

		auto vecIndsAssocParts = indsToDerInds(assocParts(), 4);

		// Compute Log-likelihood of all priors
		arma::vec logLikelihoods(gpisPriors_.size());
		for(unsigned i = 0; i < gpisPriors_.size(); i++) {
			vec fPlusVecAssoc = gpisPriors_[i]->computeDataDeviation(parts, assocParts(), priorParameters_[i]);
			mat covariance = gpisPriors_[i]->fullCovariance()(vecIndsAssocParts, vecIndsAssocParts);
			mat lCholDec = priorMatrixData_[i].lCholDec();
			mat t1 = -0.5*trans(fPlusVecAssoc)*solve(trimatu(trans(lCholDec)), solve(trimatl(lCholDec), fPlusVecAssoc));

			//arma::vec logDiag = 2*log(diagvec(lCholDec));
			//mat t2 = -0.5*sum(logDiag);
			//logLikelihoods[i] = as_scalar(t1 + t2);

			logLikelihoods[i] = as_scalar(t1);
		}
		// Normalize and compute probabilities
		arma::vec probs = exp(logLikelihoods - max(logLikelihoods));

		// Sample from probabilities
		unsigned newPriorType = sampleUniVarCat(probs);

		if (priorType_ != newPriorType) {
			priorType_ = newPriorType;
			for (unsigned i = 0; i < numAssocParts_; i++) {
				parts[assocParts_(i)].cholUpdateData_.firstDifferentTrIndex_ = 0;
			}
		}


//		priorType_ = 0;
	}



}	//	namespace gpis 
