/*
 * SamplerOutput.h
 *
 *  Created on: 18/01/2016
 *      Author: rfitch
 */

#ifndef SAMPLEROUTPUT_H_
#define SAMPLEROUTPUT_H_

#define ARMA_NO_DEBUG
#include <armadillo>
#include<string>

namespace gpis {
	class SamplerOutput {
	public:
		/////////////////////////////////
		// constructors and destructors
		SamplerOutput(const int numObjects, const arma::mat &objectLocations);
		SamplerOutput(const SamplerOutput &other);

		//////////////////
		// public methods
		const int numObjects() const { return numObjects_; };
		const arma::mat& objectLocations() const { return objectLocations_; };

		const std::string toString() const;

	private:
		///////////////////
		// private members
		const int numObjects_;
		const arma::mat objectLocations_;
	};
}	// namespace gpis

#endif /* SAMPLEROUTPUT_H_ */
