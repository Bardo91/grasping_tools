///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_GPISCELL_H_
#define GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_GPISCELL_H_

#include <armadillo>
#include <functional>
#include <vector>

namespace gpis {
	struct GpisCell {
		enum class eCellState {disabled, expanded, unexpanded};
		
		typedef std::function<double(const arma::vec3 &_point)> evalFun;
		
		GpisCell(arma::vec2 _xLimits, arma::vec2 _yLimits, arma::vec2 _zLimits, evalFun &_evalFun) {
			xLimits = _xLimits;
			yLimits = _yLimits;
			zLimits = _zLimits;
			centroid[0] = (_xLimits[0] + _xLimits[1]) / 2;
			centroid[1] = (_yLimits[0] + _yLimits[1]) / 2;
			centroid[2] = (_zLimits[0] + _zLimits[1]) / 2;
			gpValue = _evalFun(centroid);
			state = eCellState::unexpanded;
			childs = nullptr;
		}

		GpisCell() {
			state = eCellState::disabled;
		}

		//-------------------------------------------------------------------------------------------------------------
		arma::vec2	xLimits;
		arma::vec2	yLimits;
		arma::vec2	zLimits;

		arma::vec3	centroid;

		double		gpValue;

		eCellState	state;

		GpisCell*	childs;
	};

}	//	namespace gpis

#endif	//	GPIS_UTILS_VISUALIZATION_SURFACEOCTOMAP_GPISCELL_H_
