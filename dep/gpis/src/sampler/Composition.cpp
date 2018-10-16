/*
 * Composition.cpp
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      stores and manages a composition of parts and objects -- basically the
 *      belief about a scene
 *
 */

#include <gpis/sampler/Composition.h>
#include <gpis/utils/utilsMisc.h>
#include <gpis/utils/utilsProb.h>
#include <gpis/utils/utilsGPIS.h>
#include <gpis/sampler/part/CholUpdateData.h>

#ifdef PARALLEL_FOR_ENABLED
	#include <omp.h>
#endif
#include <chrono>
#include <numeric>
#include <string>
#include <algorithm>


using namespace std;
using namespace arma;

namespace gpis {
	///////////////////////////////////////////////////////////////////////////////////
	// constructors and destructors
	//
	Composition::Composition(const OverSegmentation &overSeg, const std::vector<GPISPrior*> &gpisPriors)
		: numDims_(overSeg.numDims_),
		numDimsPlusOne_(numDims_ + 1),
		numObjects_(0),
		maxNumObjects_(numObjects_),
		numParts_(overSeg.numParts()),
		gpisPriors_(gpisPriors) {

		assoc_.zeros(numParts_);
		numAssocPartsVec_.zeros(0);

		// loop through parts - all this does is interpret the information from the oversegmentation
		for (int n = 0; n < numParts_; ++n) {
			Part part(n, overSeg.allParts(n).location, overSeg.allParts(n).surfaceNormal, overSeg.allParts(n).voxels);
			parts_.push_back(part);
		}

		// Initialize all the priors 
		for (auto prior : gpisPriors_) {
			prior->init(parts_, true);
		}

		// now create an object for each part
		for (int k = 0; k < numParts_; ++k) {
			addObject(k);
		}
	}


	///////////////////////////////////////////////////////////////////////////////////
	// public methods
	//
	void Composition::updateParts(const double alpha) {
		for (int partIndex = 0; partIndex < numParts_; ++partIndex) {
			updatePart(partIndex, alpha);
		}
	}

	void Composition::updateObjects() {
		//#ifdef PARALLEL_FOR_ENABLED
		//	#pragma omp parallel for
		//	for (int objectIndex = 0; objectIndex < numObjects_; objectIndex++) {
		//		updateObject(objectIndex);
		//	}
		//#else
			for (int objectIndex = 0; objectIndex < numObjects_; objectIndex++) {
				updateObject(objectIndex);
			}
		//#endif
	}


	///////////////////////////////////////////////////////////////////////////////////
	// private methods
	//


	void Composition::updatePart(const int _partIndex, const double _alpha) {
		// Figure out which object the part belongs to
		int oldObjectIndex = assoc_(_partIndex);

		// Temporarily remove from that object
		numAssocPartsVec_(oldObjectIndex)--;

		// compute probabilities for all objects
		vec assocProbs = computeAllAssocProbs(_partIndex, oldObjectIndex);

		// Sample new association
		const int newObjectIndex = CRPSample(assocProbs,
											 numAssocPartsVec_.head(numObjects_),
											 _alpha);
		numAssocPartsVec_(newObjectIndex)++;

		if (oldObjectIndex == newObjectIndex) { // part is reassociated to old object. Nothing changes.
			return;
		}
		// Part gets associated to different object (new or existing)

		//Add part to new object (do this first, in case the object vector gets messed up due to old object's death)
		addPartToObject(_partIndex, newObjectIndex);
		// Update association
		assoc_(_partIndex) = newObjectIndex;
		// Clean up CholUpdateData for this part
		parts_[_partIndex].cholUpdateData_.firstDifferentTrIndex_ = 0;

		if ((numAssocPartsVec_(oldObjectIndex) > 0)){ //old object lives on
			// Deal with the old object
			// Find part in assocParts for that object
			int locationInAssocParts = gpisObjects_[oldObjectIndex].findPartInAssocParts(_partIndex);
 
			gpisObjects_[oldObjectIndex].permanentlyRemoveAssociatedPart(_partIndex, locationInAssocParts, parts_[_partIndex].cholUpdateData_);
			// Clean up CholUpdateData of all parts in old object
			// Adapt first different trailing index for all parts associated to that object
			for (int n = 0; n < numAssocPartsVec_(oldObjectIndex); n++) {
				// If part n comes before the removed part, the distance to that part indicates the maximum new first different index
				// If part n comes after the removed part, all is lost and needs to be recomputed (first different index = 0)
				parts_[gpisObjects_[oldObjectIndex].assocParts(n)].cholUpdateData_.firstDifferentTrIndex_ =
					std::max(std::min(locationInAssocParts - n - 1, parts_[gpisObjects_[oldObjectIndex].assocParts(n)].cholUpdateData_.firstDifferentTrIndex_), 0);
			}
		}
		else
		{
			//The old object has to die.
			removeObject(oldObjectIndex);
		}

	}



	void Composition::addPartToObject(const int _partIndex, const int _objectIndex){
		if (_objectIndex < numObjects_){ // Add to existing object
			// add part to  object
			gpisObjects_[_objectIndex].addAssociatedPart(_partIndex);
		}
		else {// This wants to add to a new object. We need to create one.
			addObject(_partIndex);
		}

		// done
	}


	///////////////////////////////////////////////////////////////////////////////////
	// private methods
	//

	void Composition::addObject(const int _partIndex) {
		GPISObject obj(gpisPriors_, parts_, {(unsigned) _partIndex }, numDims_);
		gpisObjects_.push_back(obj);
		assoc_(_partIndex) = numObjects_;
		numObjects_++;
		if (numObjects_ > maxNumObjects_) {
			numAssocPartsVec_.resize(numObjects_);
			maxNumObjects_ = numObjects_;
		}
		numAssocPartsVec_(numObjects_ - 1) = 1;
		assoc_(_partIndex) = numObjects_ - 1;
		// and update the associated data
		updateObject(numObjects_ - 1);
	}


	void Composition::removeObject(const int _objectIndex) {
		// blow away the object
		gpisObjects_.erase(gpisObjects_.begin() + _objectIndex);
		numObjects_--;

		for (int k = _objectIndex; k < numObjects_; k++) {
			numAssocPartsVec_(k) = numAssocPartsVec_(k + 1);
		}
		// now fix up the indices in the assoc vector. All the object indices
		// greater than the one we just blew away will be decremented. So we just fix them
		// up in the assoc vector.
		uvec indicesToDecrement = find(assoc_ > _objectIndex);
		assoc_.elem(indicesToDecrement) -= 1;
	}



	vec Composition::computeAllAssocProbs(const int _partIndex, int _prevObject) {
		vec assocLogProbs(numObjects_);
		// loop through all the objects
		for (int objectIndex = 0; objectIndex < numObjects_; objectIndex++) {
			if (numAssocPartsVec_(objectIndex) > 0) {
				assocLogProbs(objectIndex) = gpisObjects_[objectIndex].computeAssocLogProb( parts_, _partIndex, (_prevObject == objectIndex));
			}
			else {
				assocLogProbs(objectIndex) = -std::numeric_limits<double>::max(); // set probability to zero (won't be used anyway)
			}
		}
		assocLogProbs -= assocLogProbs.max();
		return exp(assocLogProbs);

	}

	void Composition::updateObject(const int objectIndex) {
		gpisObjects_[objectIndex].updatePriorParameters(parts_);
		gpisObjects_[objectIndex].updatePriorType(parts_);
	}


}	//	namespace gpis
