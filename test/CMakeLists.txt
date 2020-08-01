cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TOBS_TEST)

if(MSVC)
	# With C89 we're forced to use non-secure functions
	ADD_DEFINITIONS (/bigobj /D _CRT_SECURE_NO_WARNINGS)

endif(MSVC)

file(GLOB CPPSRC
    "*.cpp")

########################
########################
# Cplex Libraries

find_package(Cplex REQUIRED)

########################
########################
# Gtest Libraries

find_package(GTest REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CPLEX_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS})
add_executable(TOBS_TEST ${CPPSRC}) 

TARGET_LINK_LIBRARIES(TOBS_TEST 
	PUBLIC TOBS_libs
	PRIVATE ${CPLEX_LIBRARIES} ${GTEST_LIBRARIES}
	)

