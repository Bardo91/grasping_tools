/*
 * Sampler.cpp
 *
 *  Created on: 26/10/2015
 *  Updated on: 9/1/2016
 *      Author: wolfram
 *
 *      the guts of the GPIS sampler algorithm
 */

#include "Sampler.h"
#include <Utils/utilsGPIS.h>
#include <Utils/utilsMisc.h>
#include <Sampler/Prior/PriorFactory.h>

#include <chrono>
#include <thread>

#include <algorithm>
#include <math.h>

using namespace std;

namespace gpis {
	//--------------------------------------------------------------------------------------------------------------------
	// constants
	//--------------------------------------------------------------------------------------------------------------------
	const double Sampler::VIEW_ANGLE_DEFAULT = 270.0;
	const int Sampler::PLOT_DIMS_1_DEFAULT = OverSegmentation::NUM_DIMS;
	const int Sampler::PLOT_DIMS_2_DEFAULT = 2;


	//--------------------------------------------------------------------------------------------------------------------
	Sampler::Sampler(	const OverSegmentation &overSeg,
						const double alpha,
						const std::vector<PriorData> &priorParameterList,
						boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer) : alpha_(alpha) {

		// Generate priores
		for (auto data : priorParameterList) {
			PriorFactory::eShapeTypes shapeType;
			if (data.mShapeName == "sphere") {
				shapeType = PriorFactory::eShapeTypes::Sphere;
			}
			else if (data.mShapeName == "sphereMH") {
				shapeType = PriorFactory::eShapeTypes::SphereMH;
			}
			else if (data.mShapeName == "sphereRadiusMH") {
				shapeType = PriorFactory::eShapeTypes::SphereRadiusMH;
			}
			else if (data.mShapeName == "plane") {
				shapeType = PriorFactory::eShapeTypes::Plane;
			}
			else if (data.mShapeName == "cylinderMH") {
				shapeType = PriorFactory::eShapeTypes::CylinderMH;
			}
			else if (data.mShapeName == "ellipseMH") {
				shapeType = PriorFactory::eShapeTypes::EllipseMH;
			}
			else {
				assert(false);
			}

			PriorFactory::eKernelTypes kernelType;
			if (data.mKernelName == "squaredExp") {
				kernelType = PriorFactory::eKernelTypes::SquaredExponential;
			}
			else {
				assert(false);
			}

			auto prior = PriorFactory::create(shapeType, data.mShapeParams, kernelType, data.mKernelParams);
			gpisPriors_.push_back(prior);
		}

		// Create composition
		composition_ = new Composition(overSeg, gpisPriors_);

		// Resize association map.
		assocMap_.zeros(composition_->numParts(), composition_->numParts());

		if (viewer) {
			m3dViewer = viewer;
			m3dViewer->addCoordinateSystem(0.1, "XYZ_map");
		}
	}

	//--------------------------------------------------------------------------------------------------------------------
	int Sampler::generateSamples(unsigned numSamples) {
		double start = getUnixTimeMillis();
		for (unsigned i = 0; i < numSamples; i++) {
			std::cout << "---------------------------------------" << std::endl;
			cout << "Sample: " << i << endl;
			updateSample();
		}
		double stop = getUnixTimeMillis();

		entropyParts_.zeros(assocMap_.n_rows);
		entropyMap_.zeros(assocMap_.n_rows, assocMap_.n_cols);
		for (unsigned i = 0; i < assocMap_.n_rows; i++) {
			for (unsigned j = 0; j < assocMap_.n_cols; j++) {
				auto prob = assocMap_(i,j) / numSamples;
				if (prob == 0) {
					entropyMap_(i, j) = 0;
				}
				else {
					entropyMap_(i, j) = -prob*log(prob) - (1 - prob)*log(1 - prob);
				}
			}
			entropyParts_[i] = arma::max(entropyMap_.row(i));
		}
		//if (m3dViewer) { // Testing entropy
		//	plotEntropyMap();
		//}

		cout << "Total time: " << (stop - start) << " ms" << endl;

		std::cout << "---------------------------------------" << std::endl;
		std::cout << "Locations: " << std::endl;
		std::cout << getOutput().objectLocations() << std::endl;
		std::cout << "---------------------------------------" << std::endl;

		return 1;
	}

	//--------------------------------------------------------------------------------------------------------------------
	SamplerOutput Sampler::getOutput() const {
		arma::mat objectLocations(composition_->objects()[0].parameters().size(), composition_->numObjects());
		int objCounter = -1;
		std::vector<GPISObject>::const_iterator iter;
		for (iter = composition_->objects().begin(); iter != composition_->objects().end(); ++iter) {
			objCounter++;
			objectLocations.col(objCounter) = iter->parameters();
		}

		SamplerOutput out(composition_->numObjects(), objectLocations);

		return out;
	}

	//--------------------------------------------------------------------------------------------------------------------
	Composition & Sampler::composition() {
		return *composition_;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void Sampler::updateSample() {
		double sampleStartTime = getUnixTimeMillis();
		double startTime;
		double stopTime;
		// Update parts
		startTime = getUnixTimeMillis();
		composition_->updateParts(alpha_);
		stopTime = getUnixTimeMillis();
		cout << "Parts time: " << (stopTime - startTime) << " ms" << endl;

		// update objects
		startTime = getUnixTimeMillis();
		composition_->updateObjects();
		stopTime = getUnixTimeMillis();
		cout << "Objects time: " << (stopTime - startTime) << " ms" << endl;

		// Store result in association map.
		auto & objects = composition_->objects();
		for (auto object : objects) {
			auto assocParts = object.assocParts();
			assocMap_(assocParts, assocParts) += arma::ones<arma::mat>(assocParts.size(), assocParts.size());
		}

		if (m3dViewer) {
			plotSample();
		}

	}

	//--------------------------------------------------------------------------------------------------------------------
	void Sampler::plotSample() {
		m3dViewer->removeAllPointClouds();

		unsigned numObjs = composition_->numObjects();
		std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, Eigen::aligned_allocator<pcl::PointCloud<pcl::PointXYZRGB>::Ptr>> voxelsPerObj;
		std::vector<std::array<int, 3>> colorList;
		for (int k = 0; k < numObjs; ++k) {
			colorList.push_back({ rand() % 255, rand() % 255, rand() % 255 });
			voxelsPerObj.push_back(pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>));
		}

		for (int k = 0; k < numObjs; ++k) {
			// Fill point
			pcl::PointXYZRGB p(colorList[k][0], colorList[k][1], colorList[k][2]);
			auto location = composition_->objects()[k].parameters();
			p.x = location[0];
			p.y = location[1];
			p.z = location[2];

			// Insert point into cloud
			pcl::PointCloud<pcl::PointXYZRGB> centroidCloud;
			centroidCloud.push_back(p);

			// Draw centroid
			m3dViewer->addPointCloud(centroidCloud.makeShared(), "centroid_" + std::to_string(k));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "centroid_" + std::to_string(k));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "centroid_" + std::to_string(k));
		}

		for (int n = 0; n < composition_->numParts(); ++n) {
			unsigned objIndex = composition_->associationVector()[n];
			if (composition_->parts()[n].numVoxels() == 0) {
				pcl::PointXYZRGB point;
				point.x = composition_->part(n).location()[0];
				point.y = composition_->part(n).location()[1];
				point.z = composition_->part(n).location()[2];
				point.r = colorList[objIndex][0];
				point.g = colorList[objIndex][1];
				point.b = colorList[objIndex][2];
				composition_->part(n).addPoint(point);
			}
			else {
				for (int i = 0; i < composition_->parts()[n].numVoxels(); i++) {
					// Fill point
					composition_->part(n).voxelPoint(i).r = colorList[objIndex][0];
					composition_->part(n).voxelPoint(i).g = colorList[objIndex][1];
					composition_->part(n).voxelPoint(i).b = colorList[objIndex][2];

				}
			}
			m3dViewer->addPointCloud(composition_->parts()[n].voxelPointCloud().makeShared(), "parts_" + std::to_string(n));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "parts_" + std::to_string(n));
		}

		m3dViewer->spinOnce(30);
		std::this_thread::sleep_for(std::chrono::microseconds(30));

	}
	//--------------------------------------------------------------------------------------------------------------------
	void heatmap(float value, int &red, int &green, int &blue) {
		int aR = 0;   int aG = 0; int aB = 255;  // RGB for our 1st color (blue in this case).
		int bR = 255; int bG = 0; int bB = 0;    // RGB for our 2nd color (red in this case).

		red = (float)(bR - aR) * value + aR;      // Evaluated as -255*value + 255.
		green = (float)(bG - aG) * value + aG;      // Evaluates as 0.
		blue = (float)(bB - aB) * value + aB;      // Evaluates as 255*value + 0.
	}

	void Sampler::plotEntropyMap() {
		m3dViewer->removeAllPointClouds();
		std::vector<std::array<int, 3>> colorList;
		unsigned numObjs = composition_->numObjects();
		std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, Eigen::aligned_allocator<pcl::PointCloud<pcl::PointXYZRGB>::Ptr>> voxelsPerObj;
		for (int k = 0; k < numObjs; ++k) {
			int r, g, b;
			heatmap(entropyParts_[k], r, g, b);
			colorList.push_back({ r,g,b });
			voxelsPerObj.push_back(pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>));
		}

		for (int k = 0; k < numObjs; ++k) {
			// Fill point
			pcl::PointXYZRGB p(colorList[k][0], colorList[k][1], colorList[k][2]);
			auto location = composition_->objects()[k].parameters();
			p.x = location[0];
			p.y = location[1];
			p.z = location[2];

			// Insert point into cloud
			pcl::PointCloud<pcl::PointXYZRGB> centroidCloud;
			centroidCloud.push_back(p);

			// Draw centroid
			m3dViewer->addPointCloud(centroidCloud.makeShared(), "centroid_" + std::to_string(k));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "centroid_" + std::to_string(k));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "centroid_" + std::to_string(k));
		}

		for (int n = 0; n < composition_->numParts(); ++n) {
			unsigned objIndex = composition_->associationVector()[n];
			if (composition_->parts()[n].numVoxels() == 0) {
				pcl::PointXYZRGB point;
				point.x = composition_->part(n).location()[0];
				point.y = composition_->part(n).location()[1];
				point.z = composition_->part(n).location()[2];
				point.r = colorList[objIndex][0];
				point.g = colorList[objIndex][1];
				point.b = colorList[objIndex][2];
				composition_->part(n).addPoint(point);
			}
			else {
				for (int i = 0; i < composition_->parts()[n].numVoxels(); i++) {
					// Fill point
					composition_->part(n).voxelPoint(i).r = colorList[objIndex][0];
					composition_->part(n).voxelPoint(i).g = colorList[objIndex][1];
					composition_->part(n).voxelPoint(i).b = colorList[objIndex][2];

				}
			}
			m3dViewer->addPointCloud(composition_->parts()[n].voxelPointCloud().makeShared(), "parts_" + std::to_string(n));
			m3dViewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "parts_" + std::to_string(n));
		}

		m3dViewer->spinOnce(30);
		std::this_thread::sleep_for(std::chrono::microseconds(30));

	}
}	// namespace gpis
