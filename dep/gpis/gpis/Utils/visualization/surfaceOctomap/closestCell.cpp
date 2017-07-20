///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "closestCell.h"

namespace gpis {
    GpisCell & closestCell(const arma::vec3 & _point, GpisCell &_root) {
		if (_root.state == GpisCell::eCellState::expanded) {
			arma::vec3 vector = _point - _root.centroid;
			unsigned childIndex = 0;

			if (vector[0] > 0)
				childIndex += 1;

			if (vector[1] > 0)
				childIndex += 2;

			if (vector[2] > 0)
				childIndex += 4;

			return closestCell(_point, _root.childs[childIndex]);
		}
		else {
			return _root;
		}
	}
}
