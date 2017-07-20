/*
 * OverSegmentation.cpp
 *
 *  Created on: 26/10/2015
 *      Author: wolfram
 *
 *		reads in a list of parts from a file. this represents the initial oversegmented
 *      view of the world
 *
 */

#include "OverSegmentation.h"

using namespace std;

namespace gpis {
	//////////////////////////////////////////////////////////////////////////////////
	// Constants
	//
	const int OverSegmentation::NUM_DIMS = 3;


	//////////////////////////////////////////////////////////////////////////////////
	// Constructors and initialization
	//
	OverSegmentation::OverSegmentation(const std::string &filename)
		: numDims_(NUM_DIMS), numParts_(0), allParts_(0), initialised_(false),
		sceneLimits_({ {10000,-10000,0},
			  {10000,-10000,0},
			  {10000,-10000,0} })
	{
		// read in the data from a file
		initialised_ = readFromFile(filename);
	}

	//---------------------------------------------------------------------------------------------------------------------
	OverSegmentation::OverSegmentation(const pcl::PointCloud<pcl::PointNormal>& _cloud) :
		numDims_(3),
		numParts_(0),
		allParts_(0),
		initialised_(false),
		sceneLimits_({ {10000,-10000,0},
						{10000,-10000,0},
						{10000,-10000,0} }) {
		// Parse cloud
		fillWithCloud(_cloud);
	}

	//---------------------------------------------------------------------------------------------------------------------
	OverSegmentation::OverSegmentation() : numDims_(3),
		numParts_(0),
		allParts_(0),
		initialised_(false),
		sceneLimits_({ {10000,-10000,0},
						{10000,-10000,0},
						{10000,-10000,0} }) {
	}


	//////////////////////////////////////////////////////////////////////////////////
	// Public methods
	//
	bool OverSegmentation::isInitialised() const {
		return initialised_;
	}

	//---------------------------------------------------------------------------------------------------------------------
	arma::colvec OverSegmentation::getSceneLimits(const unsigned int colIndex) const {
		return sceneLimits_.col(colIndex);
	}

	//---------------------------------------------------------------------------------------------------------------------
	void OverSegmentation::fillWithCloud(const pcl::PointCloud<pcl::PointNormal>& _cloud) {
		PartInfo singlePart;
		singlePart.location.resize(numDims_);
		singlePart.surfaceNormal.resize(numDims_);
		for (auto point : _cloud) {
			if (isnan(point.normal_x) || isnan(point.normal_y) || isnan(point.normal_z))
			{
				continue;
			}
			singlePart.location(0) = point.x;
			singlePart.location(1) = point.y;
			singlePart.location(2) = point.z;
			singlePart.surfaceNormal(0) = point.normal_x;
			singlePart.surfaceNormal(1) = point.normal_y;
			singlePart.surfaceNormal(2) = point.normal_z;

			for (int d = 0; d < numDims_; d++) {
				if (singlePart.location(d) < sceneLimits_(d, 0)) {
					sceneLimits_(d, 0) = singlePart.location(d);
				}
				if (singlePart.location(d) > sceneLimits_(d, 1)) {
					sceneLimits_(d, 1) = singlePart.location(d);
				}

			}

			sceneLimits_.col(2) = sceneLimits_.col(1) - sceneLimits_.col(0);
			allParts_.push_back(singlePart);
		}
		numParts_ = allParts_.size();
		initialised_ = true;
	}

	void OverSegmentation::fillWithCloud(const pcl::PointCloud<pcl::PointNormal>& _cloud, const std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr>& _pointsPerCluster) {
		PartInfo singlePart;
		singlePart.location.resize(numDims_);
		singlePart.surfaceNormal.resize(numDims_);
		for (unsigned i = 0; i < _cloud.size(); i++) {
			auto point = _cloud[i];
			if (isnan(point.normal_x) || isnan(point.normal_y) || isnan(point.normal_z))
			{
				continue;
			}
			singlePart.location(0) = point.x;
			singlePart.location(1) = point.y;
			singlePart.location(2) = point.z;
			singlePart.surfaceNormal(0) = point.normal_x;
			singlePart.surfaceNormal(1) = point.normal_y;
			singlePart.surfaceNormal(2) = point.normal_z;

			for (int d = 0; d < numDims_; d++) {
				if (singlePart.location(d) < sceneLimits_(d, 0)) {
					sceneLimits_(d, 0) = singlePart.location(d);
				}
				if (singlePart.location(d) > sceneLimits_(d, 1)) {
					sceneLimits_(d, 1) = singlePart.location(d);
				}

			}

			singlePart.voxels = _pointsPerCluster[i];
			sceneLimits_.col(2) = sceneLimits_.col(1) - sceneLimits_.col(0);
			allParts_.push_back(singlePart);
		}
		numParts_ = allParts_.size();
		initialised_ = true;

	}

	//---------------------------------------------------------------------------------------------------------------------
	void OverSegmentation::printInfo() const
	{
		cout << "Oversegmentation data:" << endl;
		cout << "Number of parts: " << numParts_ << endl;
		cout << "Number of dimensions: " << numDims_ << endl;

		int partIndex = -1;
		std::vector<PartInfo>::const_iterator iter;
		for (iter = allParts_.begin(); iter != allParts_.end(); ++iter)
		{
			partIndex++;
			cout << "Part " << partIndex << ": " << iter->location << iter->surfaceNormal << endl;
		}
	}


	//////////////////////////////////////////////////////////////////////////////////
	// Private methods
	//
	bool OverSegmentation::readFromFile(const std::string &filename)
	{
		ifstream inputFile(filename.c_str());

		if (inputFile.is_open())
		{
			string line;
			while (getline(inputFile, line)) {
				PartInfo singlePart;
				singlePart.location.resize(numDims_);
				singlePart.surfaceNormal.resize(numDims_);
				stringstream linestream(line);
				char numstr[20]; //CHAR? STRING? SIZE?

				// read the first theree columns -> location
				for (int d = 0; d < numDims_; d++)
				{
					linestream >> numstr;
					singlePart.location(d) = atof(numstr);
					if (singlePart.location(d) < sceneLimits_(d, 0))
					{
						sceneLimits_(d, 0) = singlePart.location(d);
					}
					if (singlePart.location(d) > sceneLimits_(d, 1))
					{
						sceneLimits_(d, 1) = singlePart.location(d);
					}
				}
				sceneLimits_.col(2) = sceneLimits_.col(1) - sceneLimits_.col(0);

				// read the second three columns -> surface normals.
				for (int d = 0; d < numDims_; d++) {
					linestream >> numstr;
					singlePart.surfaceNormal(d) = atof(numstr);
				}

				// Ignore points with invalid normals.
				if (singlePart.surfaceNormal.has_nan()) {
					continue;
				}

				// read until the end of the line in triplets -> voxels.
				while (!linestream.eof()) {
					pcl::PointXYZRGB point;
					for (int d = 0; d < numDims_; d++) {
						linestream >> numstr;
						if (d == 0)
							point.x = atof(numstr);
						else if (d == 1)
							point.y = atof(numstr);
						else if (d == 2)
							point.z = atof(numstr);
					}
					singlePart.voxels->push_back(point);
				}
				allParts_.push_back(singlePart);

			}
			inputFile.close();
			numParts_ = allParts_.size();
		}
		else
		{
			cout << "File " << filename.c_str() << " not found." << endl;
			return false;
		}

		// must be ok
		return true;
	}
}