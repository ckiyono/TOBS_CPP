cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TOBS_CPP)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

SET(TOBS_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TOBS_LIBRARY_DIR "${TOBS_PARENT_DIR}/libs")

set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


########################
########################
# Boost Libraries
#set(BOOST_ROOT "${VP_LIBRARY_DIR}/boost_1_60_0")
#set(BOOST_ROOT "E:/Downloads/libs/boost_1_60_0")

#set (Boost_USE_STATIC_LIBS ON)
#set (Boost_USE_MULTITHREADED ON)
#unset (Boost_INCLUDE_DIR CACHE)
#unset (Boost_LIBRARY_DIRS CACHE)

#INCLUDE(InstallRequiredSystemLibraries)

########################
########################
# Eigen3 Libraries
#set(EIGEN3_INCLUDE_DIR "${VP_LIBRARY_DIR}/Eigen3")

########################
########################

add_subdirectory(TOBS)
add_subdirectory(test)

