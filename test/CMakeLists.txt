cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TOBS_TEST)

file(GLOB CPPSRC
    "*.cpp")

########################
########################
# Cplex Libraries

find_package(Cplex)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CPLEX_INCLUDE_DIRS})
add_executable(TOBS_TEST ${CPPSRC}) 

TARGET_LINK_LIBRARIES(TOBS_TEST 
	PUBLIC TOBS_libs
	PRIVATE ${CPLEX_LIBRARIES}
	)

