cmake_minimum_required(VERSION 3.25)

project(RGNFilter)

file(GLOB headers *.h)
aux_source_directory(. sources)

# Eigen
if (WIN32)
include_directories(${CMAKE_SOURCE_DIR}/3rd/eigen-3.3.4/)  
elseif(UNIX)
	find_package(Eigen3 REQUIRED)
	if(EIGEN3_FOUND)
		message("-- Eigen version: ${EIGEN3_VERSION}")
		include_directories(${EIGEN3_INCLUDE_DIR})
	endif(EIGEN3_FOUND)
endif(WIN32)

set (TARGET_NAME "RGNFilter")
add_executable(${TARGET_NAME} ${headers} ${sources})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
