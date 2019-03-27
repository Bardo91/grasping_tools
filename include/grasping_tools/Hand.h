///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_HAND_H_
#define GPISGRASPING_HAND_H_

#include <grasping_tools/Grasp.h>
#include <grasping_tools/Object.h>

namespace grasping_tools {
	/// Parent class from which all hands inherit
	///
	/// Class defined in `#include <grasping_tools/Hand.h>`
	class Hand {
	public:
		template<typename ObjectType_>
		Grasp generate(ObjectType_ &_object);

		template<typename ObjectType_>
		std::vector<Grasp> generateGrasps(ObjectType_ &_object, double _resolution);

	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_HAND_H_