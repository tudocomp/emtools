############################################################################
#  CMakeLists.txt
#
#  Part of a simple STXXL example. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

# require cmake 2.6.4 (but please use 2.8.x)
cmake_minimum_required(VERSION 2.6.4)

# TODO: too many system specific paths...
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} /home/niki/opt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /scratch/usr/lib/cmake/stxxl)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /scratch/usr/lib)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /scratch/usr/include)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I/scratch/usr/include" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L/home/USERNAME/opt/lib" )

SET(STXXL_DIR "/scratch/usr/lib/cmake/stxxl" CACHE PATH "STXXL build path")
SET(SDSL_LIBRARY "/home/niki/opt/lib/libsdsl.a" CACHE FILEPATH "SDSL library")

# we first give our project a name
project(myproject)
set(CXX_STANDARD c++11)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")


# search for stxxl-config.cmake which contains the library's configuration
find_package(STXXL REQUIRED)

# print some info (this can be removed)
message(STATUS "STXXL_DIR: ${STXXL_DIR}")
message(STATUS "STXXL_CXX_FLAGS: ${STXXL_CXX_FLAGS}")
message(STATUS "STXXL_INCLUDE_DIRS: ${STXXL_INCLUDE_DIRS}")
message(STATUS "STXXL_LIBRARIES: ${STXXL_LIBRARIES}")
message(STATUS "SDSL_LIBRARY: ${SDSL_LIBRARY}")

# apply CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -W -Wall -O3 -DNDEBUG -march=native -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -W -Wall -O0 -ggdb -fopenmp")

# add STXXL include directory
include_directories(${STXXL_INCLUDE_DIRS})

# create and executable and linke with STXXL
add_executable(isaandbwt isaandbwt.cpp)
target_link_libraries(isaandbwt ${STXXL_LIBRARIES})

add_executable(readplcp readplcp.cpp)
target_link_libraries(readplcp ${STXXL_LIBRARIES} ${SDSL_LIBRARY})

add_executable(standardize standardize.cpp)
target_link_libraries(standardize ${STXXL_LIBRARIES})

add_executable(unstandardize unstandardize.cpp)
target_link_libraries(unstandardize ${STXXL_LIBRARIES})
