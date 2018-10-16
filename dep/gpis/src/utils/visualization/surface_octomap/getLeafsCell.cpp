///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <gpis/utils/visualization/surface_octomap/getLeafsCell.h>

using namespace arma;

namespace gpis {
	void getLeafsCell(const GpisCell &_cell, arma::mat &_cloud){
		switch (_cell.state) {
		case GpisCell::eCellState::unexpanded: {	// Leaf
			vec3 point;
			point[0] = _cell.centroid[0];
			point[1] = _cell.centroid[1];
			point[2] = _cell.centroid[2];
			_cloud.insert_cols(_cloud.n_cols, point);
			break;
		}
		case GpisCell::eCellState::expanded:	// Has childs
			for (unsigned i = 0; i < 8; i++) {
				getLeafsCell(_cell.childs[i], _cloud);
			}
			break;
		case GpisCell::eCellState::disabled:	// do nothing
			break;
		}

	}
}
