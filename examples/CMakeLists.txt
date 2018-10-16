macro(build_example example_name example_files)
    add_executable(${example_name} ${${example_files}})
    target_include_directories(${example_name} PUBLIC ../dep/gpis/include)
    target_link_libraries(${example_name} LINK_PUBLIC gpis)
    target_include_directories(${example_name} PUBLIC ../include)
    target_link_libraries(${example_name} LINK_PUBLIC grasping_tools)
endmacro()

set(MESHGRASPING_EXAMPLE_FILES mesh_grasping/main.cpp)
build_example(example_mesh_grasping MESHGRASPING_EXAMPLE_FILES)

