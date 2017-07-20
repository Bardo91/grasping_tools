/*
 * Sampler.h
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *      the guts of the GPIS sampler algorithm
 */

#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "../Utils/OverSegmentation.h"
#include "Composition.h"
#include "SamplerOutput.h"

#include <string>
#include <stdio.h>
#include <cstdlib>

#ifdef __linux__
#include <unistd.h>
#elif _WIN32

#endif

#include <pcl/visualization/pcl_visualizer.h>

#define ARMA_NO_DEBUG
#include <armadillo>

namespace gpis {
	struct PriorData {
		std::string			mShapeName;
		std::vector<double> mShapeParams;
		std::string			mKernelName;
		std::vector<double> mKernelParams;
	};

	class Sampler {
	public:
		/////////////////////////////////
		// constructors and destructors
		Sampler(const OverSegmentation &overSeg,
				const double alpha,
				const std::vector<PriorData> &priorParameterList,
				boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = nullptr);

		//////////////////
		// public methods
		int generateSamples(unsigned numSamples);
		SamplerOutput getOutput() const;

		Composition &composition();
		const std::vector<GPISPrior*>& priors() const { return gpisPriors_; }
	
	private:
		///////////////////
		// private methods
		void updateSample();
		void writeToFile();
		void printInfo() const;
		void plotSample();
		void plotEntropyMap();
	private:
		//////////////////////
		// private constants
		static const double VIEW_ANGLE_DEFAULT;
		static const int PLOT_DIMS_1_DEFAULT;
		static const int PLOT_DIMS_2_DEFAULT;

		///////////////////
		// private members
		std::vector<GPISPrior*> gpisPriors_;
		Composition *composition_;
		double alpha_;
		boost::shared_ptr<pcl::visualization::PCLVisualizer> m3dViewer;

		arma::mat assocMap_;
		arma::mat entropyMap_;
		arma::vec entropyParts_;

	};
}	// namespace gpis

#endif /* SAMPLER_H_ */
