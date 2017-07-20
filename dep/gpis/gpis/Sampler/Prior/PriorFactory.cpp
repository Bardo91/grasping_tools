//
//
//
//
//


#include "PriorFactory.h"

#include <Sampler/Prior/parameterSamplers/PlaneParameterSampler.h>
#include <Sampler/Prior/parameterSamplers/SphereParameterSampler.h>
#include <Sampler/Prior/parameterSamplers/SphereParameterSamplerMH.h>
#include <Sampler/Prior/parameterSamplers/CylinderParameterSamplerMH.h>
#include <Sampler/Prior/parameterSamplers/SphereVarRadiusParameterSamplerMH.h>
#include <Sampler/Prior/parameterSamplers/EllipseParameterSamplerMH.h>

// Mean types
#include <Sampler/Prior/means/SphereMean.h>
#include <Sampler/Prior/means/SphereVarRadiusMean.h>
#include <Sampler/Prior/means/PlaneMean.h>
#include <Sampler/Prior/means/CylinderMean.h>
#include <Sampler/Prior/means/EllipseMean.h>
#include <Sampler/Prior/means/ConstantMean.h>

// Kernel types
#include <Sampler/Prior/kernels/SquaredExponentialKernel.h>

using namespace arma;
using namespace std;

namespace gpis{

	GPISPrior* PriorFactory::create(eShapeTypes shapeType, const vec &shapeParams, eKernelTypes kernelType, const vec &kernelParams){
		// Configure mean function and geometric likelihood
		Mean *mean;
		function<double(vec &)> geomLikelihood;
		ParameterSampler parameterSampler;

		// Configure mean and parameter sampler
		switch (shapeType) {
		case eShapeTypes::Sphere:
			assert(shapeParams.size() == 1);
			mean = new means::SphereMean(shapeParams[0]);
			geomLikelihood = [](const vec& _point){return 1.0;}; // for testing with previous results (fruit scene)
			parameterSampler = parameterSamplers::SphereParametersSampler(shapeParams[0]);
			break;
		case eShapeTypes::SphereMH:
			assert(shapeParams.size() == 1);
			mean = new means::SphereMean(shapeParams[0]);
			geomLikelihood = [](const vec& _point) {return 1.0; }; // for testing with previous results (fruit scene)
			parameterSampler = parameterSamplers::SphereParametersSamplerMH(shapeParams[0]);
			break;
		case eShapeTypes::SphereRadiusMH:
			assert(shapeParams.size() == 2);
			mean = new means::SphereVarRadiusMean();
			geomLikelihood = [](const vec& _point) {return 1.0; }; // for testing with previous results (fruit scene)
			parameterSampler = parameterSamplers::SphereVarRadiusParametersSamplerMH(shapeParams[0], shapeParams[1]);
			break;
		case eShapeTypes::Plane:
			assert(shapeParams.size() == 0);
			mean = new means::PlaneMean();
			geomLikelihood = [](const vec& _point) {return 1.0; };
			parameterSampler = parameterSamplers::PlaneParameterSampler();
			break;
		case eShapeTypes::EllipseMH:
			assert(shapeParams.size() == 3);
			mean = new means::EllipseMean(shapeParams[0], shapeParams[1], shapeParams[2]);
			geomLikelihood = [](const vec& _point) {return 1.0;};
			parameterSampler = parameterSamplers::EllipseParametersSampler(shapeParams[0], shapeParams[1], shapeParams[2]);
			break;
		case eShapeTypes::CylinderMH:
			assert(shapeParams.size() == 2);
			mean = new means::CylinderMean(shapeParams[0], shapeParams[1]);
			geomLikelihood = [](const vec& _point) {return 1.0; }; // for testing with previous results (fruit scene)
			parameterSampler = parameterSamplers::CylinderParametersSamplerMH(shapeParams[0], shapeParams[1]);
			break;
        case eShapeTypes::Constant:
            assert(shapeParams.size() == 1);
            mean = new means::ConstantMean(shapeParams[0]);
            geomLikelihood = [](const vec& _point) {return 1.0; }; // for testing with previous results (fruit scene)
            parameterSampler = [](const std::vector<Part> &allParts, const arma::uvec &assocParts, const arma::mat &lCholDec, const arma::vec &_oldParameters){return _oldParameters;};
            break;

		default:
			return nullptr;
			break;
		}

		// Configure kernel
		Kernel *kernel;
		switch (kernelType) {
		case eKernelTypes::SquaredExponential:
			assert(kernelParams.size() == 4);
			kernel = new kernels::SquaredExponentialKernel(kernelParams[0], kernelParams[1], kernelParams[2], kernelParams[3]);
			break;
		case eKernelTypes::Linear:
		default:
			return nullptr;
			break;
		}

		// Configure geometric likelihood
		return new GPISPrior(kernel, mean, geomLikelihood, parameterSampler);
	}
}
