/*
 * SamplerOutput.cpp
 *
 *  Created on: 18/01/2016
 *      Author: rfitch
 */

#include <gpis/sampler/SamplerOutput.h>

#include<sstream>

using namespace std;

namespace gpis {
	//////////////////////////////////////////////////////////////////////////////////
	// constructors and initialisation
	//
	SamplerOutput::SamplerOutput(const int numObjects, const arma::mat &objectLocations)
		: numObjects_(numObjects), objectLocations_(objectLocations)
	{
		// empty
	}

	SamplerOutput::SamplerOutput(const SamplerOutput &other)
		: numObjects_(other.numObjects_), objectLocations_(other.objectLocations_)
	{
		// empty
	}


	//////////////////////////////////////////////////////////////////////////////////
	// public methods
	//
	const std::string SamplerOutput::toString() const
	{
		std::stringstream ss;
		ss << "Num objects: " << numObjects_ << endl;
		ss << "Locations: " << objectLocations_;
		return ss.str();
	}
}	// namespace gpis
