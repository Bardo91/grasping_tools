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
find_package(PCL QUIET REQUIRED)
if(${PCL_FOUND})
	MESSAGE( STATUS "Adding properly PCL and its dependencies")
	# -march=native -msse4.2 -mfpmath=sse -DDISABLE_PCAP-DDISABLE_PNG-DDISABLE_LIBUSB_1_0
	list(APPEND GRASPING_TOOLS_DEFINITIONS "ENABLE_PCL")
	list(APPEND GRASPING_TOOLS_DEFINITIONS "DISABLE_PCAP-DDISABLE_PNG")
	list(APPEND GRASPING_TOOLS_DEFINITIONS "DISABLE_LIBUSB_1_0")
	set(GRASPING_TOOLS_COMPILE_OPTIONS "${GRASPING_TOOLS_COMPILE_OPTIONS} -march=native -msse4.2 -mfpmath=sse -L/usr/lib/x86_64-linux-gnu")
	if(${PCL_VERSION} VERSION_GREATER 1.8)
		list(APPEND GRASPING_TOOLS_DEFINITIONS "HAS_PCL_1_8")
	else()
		list(APPEND GRASPING_TOOLS_DEFINITIONS "HAS_PCL_1_7")
	endif()
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


# Collect all sources
file(GLOB_RECURSE GPIS_SOURCE_FILES "src/*.cpp")
file(GLOB_RECURSE GPIS_HEADER_FILES "include/*.h" "include/*.inl")

add_library(${PROJECT_NAME} ${GPIS_SOURCE_FILES} ${GPIS_HEADER_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "include") 

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${ARMADILLO_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${Boost_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PCL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${PCL_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${OpenCV_INCLUDE_DIRS})  ## REMOVE
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${OpenCV_LIBS})


#########################################
######			INSTALL 		   ######
#########################################
include(CMakePackageConfigHelpers)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif() 

install(TARGETS gpis
        EXPORT gpisTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}
        )

write_basic_package_version_file(
    gpisConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )


# This installs the include folder
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.h")
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.inl")

export(TARGETS gpis NAMESPACE gpis:: FILE gpisTargets.cmake)

# This generates gpisTargets.cmake
install(EXPORT gpisTargets
    FILE gpisTargets.cmake
    NAMESPACE gpis::
    DESTINATION lib/cmake/gpis
    )

configure_file(gpisConfig.cmake.in gpisConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/gpisConfig.cmake"
        DESTINATION lib/cmake/gpis
        )



