///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "SurfaceGpis.h"
#include <pcl/visualization/pcl_visualizer.h>
#include <chrono>

using namespace arma;
using namespace pcl;
using namespace std;

namespace gpis {
	//--------------------------------------------------------------------------------------------------------------------
	//	Public interface
	//--------------------------------------------------------------------------------------------------------------------
	SurfaceGpis::SurfaceGpis(Kernel *_kernel, Mean *_mean) {
		mMean = _mean;
		mKernel = _kernel;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::addData(const mat & _data) {
		mData.insert_cols(mData.n_cols, _data.rows(0, 2));
		mDataGrads.insert_cols(mDataGrads.n_cols, _data.rows(3, 5));

		for (unsigned i = 0; i < mData.n_cols; i++) {
			mCloudData.push_back(PointXYZ(mData.col(i)[0], mData.col(i)[1], mData.col(i)[2]));
		}
	}

	void SurfaceGpis::removeData() {
		mData = arma::mat();
		mDataGrads = arma::mat();

		if (mPolygons.size() != 0) {
			mPolygons.clear();
			mCloud.clear();
			mCloudData.clear();
		}
	}

	bool SurfaceGpis::checkTimeLimit(double _maxSeconds){
		auto t1 = std::chrono::high_resolution_clock::now();
		double timeSpent = double(std::chrono::duration_cast<std::chrono::milliseconds>(t1-mT0).count())/1000.0;
		if(timeSpent > _maxSeconds)
			return true;
		else
			return false;
	}

	//--------------------------------------------------------------------------------------------------------------------
	bool SurfaceGpis::compute(const arma::vec &_priorParameters, mat _initPoints, double _minDist, unsigned _maxIters, double _maxSeconds, pcl::visualization::PCLVisualizer *_viewer) {
		mT0 = std::chrono::high_resolution_clock::now();
		if (_viewer) { // Plot data points
			pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorData(mCloudData.makeShared(), 255, 0, 0);
			_viewer->addPointCloud(mCloudData.makeShared(), colorData, "DataPoints");
			_viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "DataPoints");
			_viewer->spinOnce(30);
		}

		mPriorParameters = _priorParameters;

		// Compute F  and R for the data
		arma::mat observedVals = mDataGrads;
		observedVals.insert_rows(0,1,true);
		arma::mat fPlusData = observedVals - (*mMean)(mData, _priorParameters, true);
		fPlusData.reshape(fPlusData.n_elem, 1);
		arma::mat covariance = (*mKernel)(mData, mData, true);

		arma::vec rVector = solve(covariance, fPlusData);

		auto f = [&](const vec &_x) -> vec {
			return (*mMean)(_x, _priorParameters, true) + (*mKernel)(_x, mData, true)*rVector;
		};

		mat gradX;
		mat x;
		gradX.insert_rows(0,3);
		x.insert_rows(0,3);

		// Initial 3 points
		if (_initPoints.n_cols == 0) {
			if (mData.n_cols == 0) {
				return 0;
			}
			_initPoints.insert_cols(0, mData.col(0));
		}

		// Initialize frontiers

		// Create a face
		Vertices triangle;
		Frontier initialFrontier;
		vector<Frontier> frontiers;
		for (unsigned i = 0; i < _initPoints.n_cols; i++) {
			vector<unsigned> indices = {3*i, 3*i + 1, 3*i + 2};
			initFrontier(_initPoints.col(i), f, _minDist, indices, x, gradX, initialFrontier);
			frontiers.push_back(initialFrontier);
			triangle.vertices = indices;
			mPolygons.push_back(triangle);
			if (_viewer) {// Plot initial faces
				redrawCloud(_viewer, gradX);
				auto vertices = mPolygons.back().vertices;
				plotFace(vertices, _viewer);
			}
		}

		// Iterate while exist frontiers
		unsigned  iters = 0;
		while (frontiers.size() != 0 && iters < _maxIters) {
			if(checkTimeLimit(_maxSeconds))
				return false;

			iters++;
			// Iterate over frontiers
			for (auto iter = frontiers.begin(); iter != frontiers.end();) {
				if(checkTimeLimit(_maxSeconds))
					return false;

				// Get active edge
				unsigned index1 = iter->mIndices.back();
				unsigned index2 = iter->mIndices.front();

				if (_viewer) {	// plot active edge
					drawLine(_viewer, mCloud[index1], mCloud[index2], "activeEdge");
				}

				if (iter->mAngles.back() < M_PI / 3) {
					// Angle at ending edge node small enough to close gap immediately
					triangle.vertices = { index1, *(iter->mIndices.end() - 2), index2 };
					mPolygons.push_back(triangle);

					iter->mAngles.pop_back();
					iter->mIndices.pop_back();
					updateEdgeAngles(*iter, x, gradX, { iter->mIndices.size() - 1, 0});
					
					if (_viewer) {
						redrawCloud(_viewer, gradX);
						auto vertices = mPolygons.back().vertices;
						plotFace(vertices, _viewer);
					}
				}
				else if (iter->mAngles.front() < M_PI / 3) {
					// Angle at beginning edge node small enough to close gap immediately
					triangle.vertices = { index1, *(iter->mIndices.begin() + 1), index2 };
					mPolygons.push_back(triangle);

					iter->mAngles.erase(iter->mAngles.begin());
					iter->mIndices.erase(iter->mIndices.begin());
					updateEdgeAngles(*iter, x, gradX, { 0, iter->mIndices.size() - 1 });
					if (_viewer) {
						redrawCloud(_viewer, gradX);
						auto vertices = mPolygons.back().vertices;
						plotFace(vertices, _viewer);
					}
				}
				else {
					// Compute candidate point
					auto candidate = thirdPoint(x.col(index1), x.col(index2), gradX.col(index1), gradX.col(index2), _minDist, -1, 0.5);
					int nearIndex = checkDistance(candidate, x, _minDist, iter->mIndices);
					if (_viewer) {
						drawPoint(_viewer, PointXYZ(candidate[0], candidate[1], candidate[2]), "candidate");
						/*cv::namedWindow("lala");
						cv::waitKey();*/
					}
					if (nearIndex == -1) { // Check for intersection with other frontiers
						bool intersectWithOther = false;
						vector<Frontier>::iterator iter2;
						for (iter2 = frontiers.begin(); iter2 != frontiers.end(); iter2++) {
							if(checkTimeLimit(_maxSeconds))
								return false;
							if (iter != iter2) {
								nearIndex = checkDistance(candidate, x, _minDist, iter2->mIndices);
								if (nearIndex != -1) {
									intersectWithOther = true;
									break;
								}
							}
						}

						if (intersectWithOther) {
							// Two frontiers are merged
							triangle.vertices = { index1, iter2->mIndices[nearIndex], index2 };
							mPolygons.push_back(triangle);

							unsigned oldEnd = iter->mIndices.size() - 1;

							iter->mIndices.insert(iter->mIndices.end(), iter2->mIndices.begin() + nearIndex, iter2->mIndices.end());
							iter->mIndices.insert(iter->mIndices.end(), iter2->mIndices.begin(), iter2->mIndices.begin() + nearIndex + 1);
							iter->mAngles.insert(iter->mAngles.end(), iter2->mAngles.begin() + nearIndex, iter2->mAngles.end());
							iter->mAngles.insert(iter->mAngles.end(), iter2->mAngles.begin(), iter2->mAngles.begin() + nearIndex + 1);


                            updateEdgeAngles(*iter, x, gradX, { 0, oldEnd, oldEnd + 1, iter->mIndices.size() - 1 });
							iter = frontiers.erase(iter2);
							if (_viewer) {
								redrawCloud(_viewer, gradX);
								auto vertices = mPolygons.back().vertices;
								plotFace(vertices, _viewer);
							}
						}
						else {
							// We add a new point to the grid
							candidate = thirdPoint(x.col(index1), x.col(index2), gradX.col(index1), gradX.col(index2), _minDist, -1, sqrt(3) / 2);

							unsigned newIndex = x.n_cols;
							vec xNew, gradNew;
							newtonStep(candidate, f, xNew, gradNew);

							// Extend data 
							x.insert_cols(x.n_cols, xNew);
							gradX.insert_cols(gradX.n_cols, gradNew);
							triangle.vertices = { index1, newIndex, index2 };
							mPolygons.push_back(triangle);
							mCloud.push_back(PointXYZ(x.tail_cols(1)[0], x.tail_cols(1)[1], x.tail_cols(1)[2]));

							// Update frontier
							iter->mIndices.push_back(newIndex);
							iter->mAngles.push_back(M_PI);
                            updateEdgeAngles(*iter, x, gradX, { 0, iter->mIndices.size() - 2 });
							if (_viewer) {
								redrawCloud(_viewer, gradX);
								auto vertices = mPolygons.back().vertices;
								plotFace(vertices, _viewer);
							}
						}
					}
					else if (nearIndex == iter->mIndices.size() - 2) {
						// Candidate point close enough to previous edge point
						triangle.vertices = { index1, *(iter->mIndices.end() - 2), index2 };
						mPolygons.push_back(triangle);

						iter->mAngles.pop_back();
						iter->mIndices.pop_back();
                        updateEdgeAngles(*iter, x, gradX, { iter->mIndices.size() - 1, 0 });
						if (_viewer) {
							redrawCloud(_viewer, gradX);
							auto vertices = mPolygons.back().vertices;
							plotFace(vertices, _viewer);
						}
					}
					else if (nearIndex == 1) {
						triangle.vertices = { index1, *(iter->mIndices.begin() + 1), index2 };
						mPolygons.push_back(triangle);

						iter->mAngles.erase(iter->mAngles.begin());
						iter->mIndices.erase(iter->mIndices.begin());
                        updateEdgeAngles(*iter, x, gradX, { 0, iter->mIndices.size() - 1 });
						if (_viewer) {
							redrawCloud(_viewer, gradX);
							auto vertices = mPolygons.back().vertices;
							plotFace(vertices, _viewer);
						}
					}
					else if (nearIndex != - 1) {
						//Candidate point close enough to other point on the
						//surface: Split surface
						triangle.vertices = { index1, *(iter->mIndices.begin() + nearIndex), index2 };
						mPolygons.push_back(triangle);

						Frontier newFrontier;
						newFrontier.mIndices.insert(newFrontier.mIndices.begin(), iter->mIndices.begin() + nearIndex, iter->mIndices.end());
						newFrontier.mAngles.insert(newFrontier.mAngles.begin(), iter->mAngles.begin() + nearIndex, iter->mAngles.end());
                        updateEdgeAngles(newFrontier, x, gradX, { 0, newFrontier.mIndices.size() - 1 });

						iter->mIndices.resize(nearIndex + 1);
						iter->mAngles.resize(nearIndex + 1);
                        updateEdgeAngles(*iter, x, gradX, { 0, iter->mIndices.size() - 1 });

						iter = frontiers.insert(frontiers.end(), newFrontier);
						if (_viewer) {
							redrawCloud(_viewer, gradX);
							auto vertices = mPolygons.back().vertices;
							plotFace(vertices, _viewer);
						}
					}
				}

				if (iter != frontiers.end()) {
					iter++;
				}
			}

			for (auto iter = frontiers.begin(); iter != frontiers.end();) {
				if(checkTimeLimit(_maxSeconds))
					return false;
				if (iter->mIndices.size() < 3) {
					iter = frontiers.erase(iter);
				}
				else if (iter->mIndices.size() == 3) {
					triangle.vertices[0] = iter->mIndices[0];
					triangle.vertices[1] = iter->mIndices[1];
					triangle.vertices[2] = iter->mIndices[2];
					mPolygons.push_back(triangle);
					iter = frontiers.erase(iter);
				}
				else {
					++iter;
				}
			}
		}

		if (_viewer) {	
			_viewer->removeShape("activeEdge");
			_viewer->removePointCloud("candidate");
		}
		return true;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::data(mat & _data) const {
		_data = mData;
		_data.insert_rows(_data.n_rows, mDataGrads);
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::data(pcl::PointCloud<pcl::PointXYZ>& _data) const {
		_data = mCloudData;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::mesh(PointCloud<PointXYZ> &_cloud, vector<Vertices> &_polygons) const {
		_cloud = mCloud;
		_polygons = mPolygons;
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::mesh(pcl::PolygonMesh &_mesh) const{
		_mesh.polygons = mPolygons;
        pcl::toPCLPointCloud2(mCloud, _mesh.cloud);

	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::colorMap(double _val, double &_r, double &_g, double &_b) const {
		if(std::isnan(_val)){
			_val = 1;
		}else{
			_val = _val<0?0:_val;
			_val = _val>1?1:_val;
		}
		double size = (cParulaMap.size()) / 3;
		int idx = int(size*_val);

		_r = cParulaMap[idx * 3];
		_g = cParulaMap[idx * 3+1];
		_b = cParulaMap[idx * 3+2];
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::mesh(pcl::PointCloud<pcl::PointXYZRGB> &_cloud, std::vector<pcl::Vertices> &_polygons) const {
		arma::mat covariance = (*mKernel)(mData, mData, true);

		auto evalCovariance = [&](const arma::colvec3 &_point)->arma::mat {
			auto crossCovariance = (*mKernel)(_point, mData, true);
			auto pointCovariance = (*mKernel)(_point, true);
			return pointCovariance - crossCovariance*arma::solve(covariance, crossCovariance.t());
		};
		
		std::vector<double> covariances;
		double maxCov=0, minCov = 99999;
		for (auto p : mCloud) {
			// Compute variance.
			auto pCov = evalCovariance({ p.x,p.y,p.z });
			covariances.push_back(pCov(0, 0));
			if (pCov(0, 0) > maxCov) maxCov = pCov(0, 0);
			if (pCov(0, 0) < minCov) minCov = pCov(0, 0);
		}
		for (unsigned i = 0; i < mCloud.size(); i++) {
			double r, g, b;
			colorMap((covariances[i] - minCov) / (maxCov-minCov), r, g, b);
			pcl::PointXYZRGB cPoint(r*255,g*255,b*255);
			cPoint.x = mCloud[i].x;
			cPoint.y = mCloud[i].y;
			cPoint.z = mCloud[i].z;

			_cloud.push_back(cPoint);

		}
		_polygons = mPolygons;
	}


	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::points(mat &_points) const {

	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::points(pcl::PointCloud<pcl::PointXYZ> &_points) const {
		_points = mCloud;
	}


	//--------------------------------------------------------------------------------------------------------------------
	// Private interface
	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::initFrontier(const arma::vec & _point, std::function<arma::vec(const arma::vec&)> _fplus, double _dist, const std::vector<unsigned>& _indices, arma::mat & _x, arma::mat & _grad, Frontier & _frontier) {
		unsigned oldSize = _x.n_cols;
		_x.insert_cols(oldSize, 3);
		_grad.insert_cols(oldSize, 3);

		vec x, grad;
		newtonStep(_point, _fplus, x, grad);
		_x.col(oldSize) = x;
		_grad.col(oldSize) = grad;
		mCloud.push_back(PointXYZ(x[0], x[1], x[2]));
		vec perpVector = cross(_grad.col(oldSize), randu<vec>(grad.n_rows));
		x += _dist * perpVector / norm(perpVector);
		newtonStep(x, _fplus, x, grad);
		_x.col(oldSize + 1) = x;
		_grad.col(oldSize + 1) = grad;
		mCloud.push_back(PointXYZ(x[0], x[1], x[2]));

		x = thirdPoint(_x.col(oldSize), _x.col(oldSize + 1), _grad.col(oldSize), _grad.col(oldSize + 1), _dist, 1, sqrt(3) / 2);
		newtonStep(x, _fplus, x, grad);
		_x.col(oldSize + 2) = x;
		_grad.col(oldSize + 2) = grad;
		mCloud.push_back(PointXYZ(x[0], x[1], x[2]));

		_frontier.mIndices = _indices;
		_frontier.mAngles = {M_PI, M_PI, M_PI};
	}
	
	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::newtonStep(const arma::vec &_x0, std::function<arma::vec(const arma::vec&)> _fPlus, arma::vec &_x, arma::vec &_grad) {
		arma::vec fPlus = _fPlus(_x0);
		_grad = fPlus.tail(3);
		_x = _x0 - _grad*as_scalar(fPlus.head(1) / norm(_grad));

	}

	//--------------------------------------------------------------------------------------------------------------------
	vec gpis::SurfaceGpis::thirdPoint(const vec &_p1, const vec &_p2, const vec &_grad1, const vec &_grad2, double _minDist, int _dir, double _height) {
		vec perp = cross(_p2 - _p1, (_grad1 + _grad2) / 2);
		perp = perp / norm(perp) * _dir;
		return (_p1 + _p2) / 2 + perp * _height * _minDist;
	}

	//--------------------------------------------------------------------------------------------------------------------
	unsigned SurfaceGpis::checkDistance(const vec & _point, const mat &_x, double _minDist, const std::vector<unsigned> &_indices) {	
		double minDist = 99999;
		 unsigned minIndex = -1;
		 for (unsigned i = 1; i < _indices.size() - 1; i++) {
			 if (_indices[i] == _indices[0] || _indices[i] == _indices.back())
				 continue;

			 double dist = norm(_x.col(_indices[i]) - _point);
			 if (dist < minDist) {
				 minDist = dist;
				 minIndex = i;
			 }
		 }

		 if (minDist < _minDist) {
			 return minIndex;
		 }
		 else {
			 return -1;
		 }
	 }
	
	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::updateEdgeAngles(Frontier & _frontier, const arma::mat & _x, const arma::mat & _grads, const std::vector<unsigned>& _updateIndices) {
		std::vector<double> updatedAngles;
		for (unsigned i = 0; i < _updateIndices.size(); i++) {
			unsigned ind, ind1, ind2;
			if (_updateIndices[i] == 0) {
				ind = 0;
				ind1 = _frontier.mIndices.size() - 1;
				ind2 = 1;
			}
			else if (_updateIndices[i] == _frontier.mIndices.size() - 1) {
				ind = _frontier.mIndices.size() - 1;
				ind1 = _frontier.mIndices.size() - 2;
				ind2 = 0;
			}
			else {
				ind = _updateIndices[i];
				ind1 = _updateIndices[i] - 1;
				ind2 = _updateIndices[i] + 1;
			}

			_frontier.mAngles[ind] = edgeAngle(	_x.col(_frontier.mIndices[ind]), 
												_x.col(_frontier.mIndices[ind1]), 
												_x.col(_frontier.mIndices[ind2]), 
												_grads.col(_frontier.mIndices[ind]));
		}
	}

	//--------------------------------------------------------------------------------------------------------------------
	double SurfaceGpis::edgeAngle(const arma::vec &_x, const arma::vec &_x1, const arma::vec &_x2, const arma::vec &_gradX){
		auto dVec1 = _x1 - _x;
		auto dVec2 = _x2 - _x;

		if (as_scalar(trans(cross(dVec1, dVec2))*_gradX) < 0) {
			return acos(as_scalar(trans(dVec1)* dVec2/(norm(dVec1) * norm(dVec2))));
		}
		else {
			return M_PI;
		}

	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::plotFace(const std::vector<unsigned> &_vertices, pcl::visualization::PCLVisualizer *_viewer) {
		for (unsigned i = 0; i < _vertices.size() - 1; i++) {
			_viewer->addLine<pcl::PointXYZ>(mCloud[_vertices[i]], mCloud[_vertices[i + 1]], "line_"+to_string(_vertices[i])+"_"+ to_string(_vertices[i+1]));
		}
		_viewer->addLine<pcl::PointXYZ>(mCloud[_vertices.back()], mCloud[_vertices.front()], "line_" + to_string(_vertices.back()) + "_" + to_string(_vertices.front()));
		_viewer->spinOnce(30);
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::redrawCloud(pcl::visualization::PCLVisualizer * _viewer, arma::mat _grads) {
		pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorData(mCloud.makeShared(), 0, 0, 255);
		if (_viewer->contains("cloud"))
			_viewer->removePointCloud("cloud");
		
		_viewer->addPointCloud(mCloud.makeShared(), colorData, "cloud");
		_viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud");


		for (unsigned i = 0; i < mCloud.size(); i++) {
			auto p1 = mCloud[i];
			auto grad = _grads.col(i);
			auto p2 = p1;
			p2.x += 0.3*grad[0];
			p2.y += 0.3*grad[1];
			p2.z += 0.3*grad[2];
			drawLine(_viewer, p1, p2, "grad_" + to_string(i));
		}

		_viewer->spinOnce(1);
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::drawPoint(pcl::visualization::PCLVisualizer * _viewer, const pcl::PointXYZ & _point, const std::string & _name){
		PointCloud<PointXYZ> cloud;
		cloud.push_back(_point);
		pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorData(cloud.makeShared(), 0, 255, 0);
		if (_viewer->contains(_name))
			_viewer->removePointCloud(_name);

		_viewer->addPointCloud(cloud.makeShared(), colorData, _name);

		_viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, _name);
		_viewer->spinOnce(1);
	}

	//--------------------------------------------------------------------------------------------------------------------
	void SurfaceGpis::drawLine(pcl::visualization::PCLVisualizer * _viewer, const pcl::PointXYZ & _p1, const pcl::PointXYZ & _p2, const std::string & _name, unsigned _r, unsigned _g, unsigned _b) {
		if (_viewer->contains(_name))
			_viewer->removeShape(_name);

		_viewer->addLine<pcl::PointXYZ>(_p1, _p2, _r, _g, _b,_name);
		_viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, _name);
		_viewer->spinOnce(1);
	}
}
