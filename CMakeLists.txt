###############################################################
######                                                   ######
######  	GPIS library                                 ######
######			Author: Pablo R.S. (aka. Bardo91)        ######
######                                                   ######
###############################################################

cmake_minimum_required (VERSION 3.8 FATAL_ERROR)
project(grasping_tools VERSION 0.1)

#################################################
######### Custom project configuration ##########
#################################################

set(USE_CATKIN OFF)



if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE )
endif()


list(APPEND GRASPING_TOOLS_DEFINITIONS "")
set(GRASPING_TOOLS_COMPILE_OPTIONS "-Wall")
#########################################
###### Loading 3rd party libraries.######
#########################################
if(WIN32)
    find_package(Armadillo REQUIRED HINTS "C:/Program Files (x86)/armadillo/share/Armadillo/CMake")
else()
    find_package(Armadillo REQUIRED HINTS "/usr/local /usr/")
endif()

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

#Load catkin
if(${USE_CATKIN})
	find_package(catkin)
	if(${catkin_FOUND})
		catkin_package(
			CATKIN_DEPENDS roscpp
			INCLUDE_DIRS include
		)
		list(APPEND GRASPING_TOOLS_DEFINITIONS "RGBDTOOLS_USE_ROS")
	endif(${catkin_FOUND})
endif(${USE_CATKIN})


#########################################
######			Compile extra Libraries 		   ######
#########################################

add_subdirectory(dep/qhull)
add_subdirectory(dep/gpis)

#########################################
######			Library 		   ######
#########################################

# Collect all sources
file(GLOB_RECURSE GRASPTOOLS_SOURCE_FILES "src/*.cpp")
file(GLOB_RECURSE GRASPTOOLS_HEADER_FILES "include/*.h" "include/*.inl")

add_library(${PROJECT_NAME} STATIC ${GRASPTOOLS_SOURCE_FILES} ${GRASPTOOLS_HEADER_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include) 

target_compile_definitions(${PROJECT_NAME} PUBLIC ${GRASPING_TOOLS_DEFINITIONS})
string(REPLACE " " ";" GRASPING_TOOLS_COMPILE_OPTIONS_REPLACED_FLAGS ${GRASPING_TOOLS_COMPILE_OPTIONS})
target_compile_options(${PROJECT_NAME} PUBLIC ${GRASPING_TOOLS_COMPILE_OPTIONS_REPLACED_FLAGS})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/dep/gpis/include) 
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC gpis)

target_include_directories(${PROJECT_NAME} PUBLIC "dep/qhull/src")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    target_link_libraries(${PROJECT_NAME} LINK_PUBLIC qhullcpp_d)
    target_link_libraries(${PROJECT_NAME} LINK_PUBLIC qhullstatic_rd)
else()
    target_link_libraries(${PROJECT_NAME} LINK_PUBLIC qhullcpp)
    target_link_libraries(${PROJECT_NAME} LINK_PUBLIC qhullstatic_r)
endif()

target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${qhull_STATIC})

target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${ARMADILLO_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${Boost_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PCL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${PCL_LIBRARIES})

if(${catkin_FOUND})
	target_include_directories(${PROJECT_NAME} PUBLIC ${catkin_INCLUDE_DIRS})
	target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${catkin_LIBRARIES} )
endif(${catkin_FOUND})

#########################################
######			   TESTS           ######
#########################################
#enable_testing()

#add_subdirectory(tests)

#########################################
######			   EXAMPLES        ######
#########################################

add_subdirectory(examples)

#########################################
######			Documentation	   ######
#########################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(	doc_grasping_tools
						${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
						WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
						COMMENT "Generating API documentation with Doxygen" VERBATIM )
endif(DOXYGEN_FOUND)


#########################################
######			INSTALL 		   ######
#########################################
include(CMakePackageConfigHelpers)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif() 

install(TARGETS grasping_tools
        EXPORT grasping_toolsTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}
        )

write_basic_package_version_file(
    grasping_toolsConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )


# This installs the include folder
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.h")
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.inl")

export(TARGETS grasping_tools NAMESPACE grasping_tools:: FILE grasping_toolsTargets.cmake)

# This generates grasping_toolsTargets.cmake
install(EXPORT grasping_toolsTargets
    FILE grasping_toolsTargets.cmake
    NAMESPACE grasping_tools::
    DESTINATION lib/cmake/grasping_tools
    )

configure_file(grasping_toolsConfig.cmake.in grasping_toolsConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/grasping_toolsConfig.cmake"
        DESTINATION lib/cmake/grasping_tools
        )



# uninstall target 
if(NOT TARGET uninstall) 
  configure_file( 
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in" 
      "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" 
      IMMEDIATE @ONLY) 
 
  add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake) 
endif() 
 
