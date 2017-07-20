###############################################################
######                                                   ######
######  	GPIS library                                 ######
######			Author: Pablo R.S. (aka. Bardo91)        ######
######                                                   ######
###############################################################

cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
project(gpis)

#################################################
######### Custom project configuration ##########
#################################################

# Modify following flags to enable/disable features
set(ENABLE_PARALLEL_FOR OFF)

#################################################
######## Automatic project configuration ########
#################################################
# Collect all sources
file(GLOB_RECURSE GPIS_SOURCE_FILES "gpis/*.cpp" "gpis/*.h" "gpis/*.inl")

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

GroupSources(gpis)

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

## Load OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions("-DPARALLEL_FOR_ENABLED")
    message(STATUS "Found OpenMP, enabling parallelization")
else()
    message(STATUS "Not Found OpenMP, disabling parallelization")
endif()

#########################################
######			Library 		   ######
#########################################
set(PROJECT_NAME gpis)

add_library(${PROJECT_NAME} ${GPIS_SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/gpis)

target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${ARMADILLO_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${Boost_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PCL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${PCL_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${OpenCV_INCLUDE_DIRS})  ## REMOVE
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${OpenCV_LIBS})
