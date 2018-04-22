if(UNIX)
	# Test.
	set(TEST_NAME test_task_wrench)

	add_executable(${TEST_NAME} test_task.cpp)

	target_include_directories(${TEST_NAME} PUBLIC "../../dep/gpis")
	target_link_libraries(${TEST_NAME} LINK_PUBLIC gpis)

	target_include_directories(${TEST_NAME} PUBLIC "../../grasping_tools")
	target_link_libraries(${TEST_NAME} LINK_PUBLIC grasping_tools)

	add_test(test_task1 test_task_wrench)
endif(UNIX)