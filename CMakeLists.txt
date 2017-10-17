###############################################################
######                                                   ######
######  	GPIS library                                 ######
######			Author: Pablo R.S. (aka. Bardo91)        ######
######                                                   ######
###############################################################

cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
project(grasping_tools)

#################################################
######### Custom project configuration ##########
#################################################

#################################################
######## Automatic project configuration ########
#################################################

# Clasify sources according to folder structure.
# This creates a nice layout for the project in visual studio
# This macro is courtesy of http://www.cmake.org/pipermail/cmake/2013-November/056336.html
macro(GroupSources curdir)
	file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/${curdir}
		${PROJECT_SOURCE_DIR}/${curdir}/*)
	foreach(child ${children})
		if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/${curdir}/${child})
			GroupSources(${curdir}/${child})
		else()
			string(REPLACE "/" "\\" groupname ${curdir})
			source_group(${groupname} FILES ${PROJECT_SOURCE_DIR}/${curdir}/${child})
		endif()
	endforeach()
endmacro()

GroupSources(grasping_tools)

#########################################
###### Loading 3rd party libraries.######
#########################################
if(WIN32)
	set(ARMADILLO_INCLUDE_DIRS $ENV{ARMADILLO_INCLUDE_DIRS})
	set(ARMADILLO_LIBRARIES $ENV{ARMADILLO_LIBRARIES})
	set(Boost_USE_STATIC_LIBS ON)
else()
    find_package(Armadillo REQUIRED HINTS "/usr/local /usr/")
endif()
find_package(Boost REQUIRED COMPONENTS system)

## Load PCL
find_package(PCL 1.8 QUIET REQUIRED)
if(${PCL_FOUND} AND ${VTK_FOUND} AND ${Boost_FOUND} AND ${FLANN_FOUND} STREQUAL "TRUE")
	MESSAGE( STATUS "Adding properly PCL and its dependencies")
	add_definitions(-DENABLE_PCL)
	add_definitions(${PCL_DEFINITIONS})
else()
	MESSAGE( FATAL_ERROR "Couldnt Add PCL and/or it's dependencies")
endif()

#########################################
######			Compile extra Libraries 		   ######
#########################################

add_subdirectory(dep/qhull)
add_subdirectory(dep/gpis)

#########################################
######			Library 		   ######
#########################################

# Collect all sources
file(GLOB_RECURSE GRASPTOOLS_SOURCE_FILES "grasping_tools/*.cpp" "grasping_tools/*.h" "grasping_tools/*.inl")

set(PROJECT_NAME grasping_tools)

add_library(${PROJECT_NAME} STATIC ${GRASPTOOLS_SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/grasping_tools)

target_include_directories(${PROJECT_NAME} PUBLIC "dep/gpis/gpis")
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC gpis)

target_include_directories(${PROJECT_NAME} PUBLIC "dep/qhull/src")
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC qhull)

target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${ARMADILLO_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${Boost_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PCL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${PCL_LIBRARIES})

#########################################
######			   TESTS           ######
#########################################
enable_testing()

add_subdirectory(tests)


set_property(TARGET grasping_tools PROPERTY FOLDER "dep")
