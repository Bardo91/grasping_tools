///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <gpis/utils/visualization/surface_octomap/expandCell.h>
#include <cassert>

namespace gpis {
	void expandCell(GpisCell & cell, GpisCell::evalFun &_evalFun){
		switch (cell.state) {
		case GpisCell::eCellState::unexpanded:		// Expand cell
			cell.childs = new GpisCell[8];
			// Compute new cells;
			cell.childs[0] = GpisCell(	{ cell.xLimits[0], cell.centroid[0] },
										{ cell.yLimits[0], cell.centroid[1] },
										{ cell.zLimits[0], cell.centroid[2] },
										_evalFun);

			cell.childs[1] = GpisCell(	{ cell.centroid[0], cell.xLimits[1] },
										{ cell.yLimits[0], cell.centroid[1] },
										{ cell.zLimits[0], cell.centroid[2] },
										_evalFun);

			cell.childs[2] = GpisCell(	{ cell.xLimits[0], cell.centroid[0] },
										{ cell.centroid[1], cell.yLimits[1] },
										{ cell.zLimits[0], cell.centroid[2] },
										_evalFun);

			cell.childs[3] = GpisCell(	{ cell.centroid[0], cell.xLimits[1] },
										{ cell.centroid[1], cell.yLimits[1] },
										{ cell.zLimits[0], cell.centroid[2] },
										_evalFun);

			cell.childs[4] = GpisCell(	{ cell.xLimits[0], cell.centroid[0] },
										{ cell.yLimits[0], cell.centroid[1] },
										{ cell.centroid[2], cell.zLimits[1] },
										_evalFun);

			cell.childs[5] = GpisCell(	{ cell.centroid[0], cell.xLimits[1] },
										{ cell.yLimits[0], cell.centroid[1] },
										{ cell.centroid[2], cell.zLimits[1] },
										_evalFun);

			cell.childs[6] = GpisCell(	{ cell.xLimits[0], cell.centroid[0] },
										{ cell.centroid[1], cell.yLimits[1] },
										{ cell.centroid[2], cell.zLimits[1] },
										_evalFun);

			cell.childs[7] = GpisCell(	{ cell.centroid[0], cell.xLimits[1] },
										{ cell.centroid[1], cell.yLimits[1] },
										{ cell.centroid[2], cell.zLimits[1] },
										_evalFun);

			cell.state = GpisCell::eCellState::expanded;

			break;
		case GpisCell::eCellState::expanded:	// Go deeper
			for (unsigned i = 0; i < 8; i++) {
				expandCell(cell.childs[i],_evalFun);
			}
			break;
		case GpisCell::eCellState::disabled:	// Do nothing
			break;
		}
	}
}
