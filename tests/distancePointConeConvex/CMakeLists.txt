if(UNIX)
	# Testing class DistancePointConvexCone, that computes the distance from a n-dimensional point to a convex n-dimensional cone.
	set(TEST_NAME testCones)

	add_executable(${TEST_NAME} testCones.cpp)

	target_include_directories(${TEST_NAME} PUBLIC "../../dep/gpis")
	target_link_libraries(${TEST_NAME} LINK_PUBLIC gpis)

	target_include_directories(${TEST_NAME} PUBLIC "../../grasping_tools")
	target_link_libraries(${TEST_NAME} LINK_PUBLIC grasping_tools)

	add_test(testCones1 testCones)
endif(UNIX)