///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace gpis {
	template<typename PointType_>
	inline void getLeafsCell(const GpisCell &_cell, pcl::PointCloud<PointType_> &_cloud) {
		switch (_cell.state) {
		case GpisCell::eCellState::unexpanded: {	// Leaf
			PointType_ point;
			point.x = _cell.centroid[0];
			point.y = _cell.centroid[1];
			point.z = _cell.centroid[2];
			if (_cell.gpValue > 0)
				point.r = 255;
			else
				point.g = 255;

			_cloud.push_back(point);
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