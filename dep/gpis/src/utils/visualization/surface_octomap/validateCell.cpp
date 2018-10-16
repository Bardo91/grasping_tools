///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <gpis/utils/visualization/surface_octomap/validateCell.h>
#include <gpis/utils/visualization/surface_octomap/closestCell.h>

namespace gpis {
	void validateCell(GpisCell & _cell, GpisCell &_root) {
		if (_cell.state == GpisCell::eCellState::expanded) {
			// If cell has childs, go deeper
			for (unsigned i = 0; i < 8; i++) {
				validateCell(_cell.childs[i], _root);
			}
		}else{
			// Compute bound points.
			double incX = abs(_cell.xLimits[0] - _cell.xLimits[1]);
			double incY = abs(_cell.yLimits[0] - _cell.yLimits[1]);
			double incZ = abs(_cell.zLimits[0] - _cell.zLimits[1]);

			arma::vec3 points[6];
			for (auto &vec : points) {
				vec = _cell.centroid;
			}
			points[0][0]	+= incX;
			points[1][0]	-= incX;
			points[2][1]	+= incY;
			points[3][1]	-= incY;
			points[4][2]	+= incZ;
			points[5][2]	-= incZ;

			// Check value of bound points;
			for (auto &vec : points) {
				// Check if point is in bounds
				if (vec[0] < _root.xLimits[0] || vec[0] > _root.xLimits[1] ||
					vec[1] < _root.yLimits[0] || vec[1] > _root.yLimits[1] ||
					vec[2] < _root.zLimits[0] || vec[2] > _root.zLimits[1] ) {
					continue;
				}

				if (closestCell(vec, _root).gpValue * _cell.gpValue < 0) {
					return;
				}
			}

			// If reached this point, disable cell
			_cell.state = GpisCell::eCellState::disabled;

		}
	}
}
