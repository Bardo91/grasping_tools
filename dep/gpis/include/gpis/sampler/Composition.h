/*
 * Composition.h
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      stores and manages a composition of parts and objects -- basically the
 *      belief about a scene
 */

#ifndef GPIS_SAMPLER_COMPOSITION_H_
#define GPIS_SAMPLER_COMPOSITION_H_

#include <gpis/sampler/object/GPISObject.h>
#include <gpis/sampler/part/Part.h>
#include <gpis/utils/OverSegmentation.h>

namespace gpis {
	class Composition {
	public:
		/// Constructor
		Composition(const OverSegmentation &overSeg, const std::vector<GPISPrior*> &gpisPriors);

		/// Get number of dimensions of the problems 2D or 3D usually.
		int numDims() const { return numDims_; }

		/// Get number of objects in the scene.
		int	numObjects() const { return numObjects_; }

		/// Get a const reference to the objects in the scene.
		const std::vector<GPISObject> & objects() const { return gpisObjects_ ; }

		/// Get number of parts in the current composition.
		int	numParts() const { return numParts_; }

		/// Return list of parts
		std::vector<Part> & parts() { return parts_; }
		/// Return list of parts
		Part & part(unsigned _n) { return parts_[_n]; }

		/// Return a vector containing the associations
		const arma::uvec & associationVector() { return assoc_;  }

		/// Update parts' association.
		/// \param alpha:
		void updateParts(const double alpha);

		/// Update list of objects.
		void updateObjects();

	private:	// private members
		const unsigned numDims_;
		const unsigned numDimsPlusOne_;
		int numObjects_;
		int maxNumObjects_;
		const unsigned numParts_;
		arma::uvec numAssocPartsVec_;
		arma::uvec assoc_;
		std::vector<GPISObject> gpisObjects_;
		std::vector<Part> parts_;
		std::vector<GPISPrior*> gpisPriors_;

	private:	// private methods
		void updatePart(const int _partIndex, const double _alpha);
		void addPartToObject(const int _partIndex, const int _objectIndex);
		void addObject(const int _partIndex);
		void removeObject(const int _objectIndex);
		arma::vec computeAllAssocProbs(const int _partIndex, int _prevObject);
		void updateObject(const int objectIndex);
	};
}	//	namespace gpis
#endif /* COMPOSITION_H_ */
