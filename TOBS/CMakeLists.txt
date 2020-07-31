cmake_minimum_required(VERSION 2.8)
PROJECT(TOBS_libs)

file(GLOB CPPSRC
    "source/*.cpp")

file(GLOB INCLSRC
    "include/*.h")

########################
########################
# Cplex Libraries

find_package(Cplex)

#######################################################################################
 
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CPLEX_INCLUDE_DIRS})

add_library(TOBS_libs ${CPPSRC} ${HPPSRC} ${INCLSRC}) 

TARGET_LINK_LIBRARIES(TOBS_libs
	PRIVATE ${CPLEX_LIBRARIES}
	)

TARGET_INCLUDE_DIRECTORIES(TOBS_libs
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
	)

 