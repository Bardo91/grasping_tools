macro(build_example example_name example_files)
    add_executable(${example_name} ${${example_files}})
    target_include_directories(${example_name} PUBLIC ../dep/gpis/include)
    target_link_libraries(${example_name} LINK_PUBLIC gpis)
    target_include_directories(${example_name} PUBLIC ../include)
    target_link_libraries(${example_name} LINK_PUBLIC grasping_tools)


    target_include_directories(${example_name} PUBLIC "dep/qhull/src")
    if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        target_link_libraries(${example_name} LINK_PUBLIC qhullcpp_d)
        target_link_libraries(${example_name} LINK_PUBLIC qhullstatic_rd)
    else()
        target_link_libraries(${example_name} LINK_PUBLIC qhullcpp)
        target_link_libraries(${example_name} LINK_PUBLIC qhullstatic_r)
    endif()

    target_link_libraries(${example_name} LINK_PUBLIC ${qhull_STATIC})


endmacro()

set(MESHGRASPING_EXAMPLE_FILES mesh_grasping/main.cpp)
build_example(example_mesh_grasping MESHGRASPING_EXAMPLE_FILES)

set(POINTCLOUDGRASPING_EXAMPLE_FILES point_cloud_grasping/main.cpp)
build_example(example_point_cloud_grasping POINTCLOUDGRASPING_EXAMPLE_FILES)

