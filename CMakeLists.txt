cmake_minimum_required(VERSION 2.8.8)
project(LIBGWASPP)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set (LIBGWASPP_VERSION_MAJOR 0)
set (LIBGWASPP_VERSION_MINOR 1)
set (LIBGWASPP_VERSION_REVISION 0)

set (LIBGWASPP_COMPILE_LEVEL 0)
set (LIBGWASPP_DEBUG_LEVEL 1)

set (LIBGWASPP_COMPRESSION_LEVEL 3)

set (PROCESSOR_WORD_SIZE 64)

set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set (R_LIBRARIES_DIR "/usr/lib/")
set (R_RMATH Rmath)

configure_file (
	"${PROJECT_SOURCE_DIR}/libgwasppConfig.h.in"
	"${PROJECT_SOURCE_DIR}/src/libgwaspp/libgwasppConfig.h"
	)

find_package( ZLIB REQUIRED )
find_package( Boost COMPONENTS program_options REQUIRED )
find_package( GSL REQUIRED )

message(STATUS "Boost Include Path: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost Library Path: ${Boost_LIBRARIES}")

set(CMAKE_CXX_FLAGS "-g -fPIC -Wall -msse3 -O3")
add_subdirectory(src)
