# Base CMake file for STXXL

project(stxxl-lcpscan)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CheckIncludeFileCXX)

# disallow in-source builds

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

###############################################################################
# compilation options 

option(BUILD_TESTS "Build all stxxl test programs" OFF)

option(USE_BOOST "Use Boost libraries for threads,config,filesystem,random and date_time" OFF)

option(USE_GNU_PARALLEL "Use GNU parallel STL extensions" OFF)

###############################################################################
# enable use of "make test"

enable_testing()

###############################################################################
# enable warnings

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
endif()

###############################################################################
# find thread library

find_package(Threads)

###############################################################################
# determine large file support

include(gmxTestLargeFiles)
GMX_TEST_LARGE_FILES (GMX_LARGEFILES)

if (GMX_LARGEFILES)
  add_definitions(-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES)
else (NOT GMX_LARGEFILES)
  message(FATAL_ERROR "Large file support was not detectable.")
endif (GMX_LARGEFILES)

###############################################################################
# optional Boost libraries

if(USE_BOOST)

  find_package( Boost 1.36.0 COMPONENTS thread date_time iostreams system filesystem )

  if(Boost_FOUND)
    # globally add boost include directories
    include_directories(${Boost_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTXXL_BOOST_CONFIG -DSTXXL_BOOST_FILESYSTEM -DSTXXL_BOOST_RANDOM -DSTXXL_BOOST_THREADS -DSTXXL_BOOST_TIMESTAMP")
  else()
    message(FATAL_ERROR "Boost libraries not found. Try compilation without them.")
  endif()
  
endif(USE_BOOST)

###############################################################################
# optional GNU parallel STL mode

if(USE_GNU_PARALLEL)

  include(FindOpenMP)

  if(NOT OPENMP_FOUND)
    message(FATAL_ERROR "OpenMP not found. Try compilation without GNU parallel mode.")
  else()
    check_include_file_cxx(parallel/algorithm HAVE_PARALLEL_ALGORITHM_H) 
   
    if (NOT HAVE_PARALLEL_ALGORITHM_H)
      message(FATAL_ERROR "GNU parallel mode header not found. Try compilation without parallel mode.")
    else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

      set(STXXL_PARALLEL_MODE_EXPLICIT "1")

    endif()
  endif()

endif(USE_GNU_PARALLEL)

###############################################################################
# optional CppUnit for more tests
# -- from http://www.cmake.org/pipermail/cmake/2006-October/011446.html

find_path(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
  /usr/local/include
  /usr/include
  )

# With Win32, important to have both
if(WIN32)
  find_library(CPPUNIT_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  find_library(CPPUNIT_DEBUG_LIBRARY cppunitd
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
else(WIN32)
  # On unix system, debug and release have the same name
  find_library(CPPUNIT_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  find_library(CPPUNIT_DEBUG_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
endif(WIN32)

if(CPPUNIT_INCLUDE_DIR)
  if(CPPUNIT_LIBRARY)
    set(CPPUNIT_FOUND "YES")
    set(CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY} ${CMAKE_DL_LIBS})
    set(CPPUNIT_DEBUG_LIBRARIES ${CPPUNIT_DEBUG_LIBRARY} ${CMAKE_DL_LIBS})
  endif(CPPUNIT_LIBRARY)
endif(CPPUNIT_INCLUDE_DIR)

if (NOT CPPUNIT_FOUND)
  message("CppUnit libraries not found. Some tests will not be compiled.")
else (CPPUNIT_FOUND)
  message("CppUnit libraries found: ${CPPUNIT_INCLUDE_DIR}.")
endif (NOT CPPUNIT_FOUND)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static -funroll-loops")

###############################################################################
# configure environment for building

# create config.h with define switches
configure_file(${PROJECT_SOURCE_DIR}/include/stxxl/bits/config.h.in ${PROJECT_BINARY_DIR}/include/stxxl/bits/config.h)

# globally adds top-level include directories
include_directories( "${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include" )

# for targets using stxxl library
set(STXXL_LIBS ${STXXL_LIBS} stxxl-lcpscan ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# build libstxxl in /lib
add_subdirectory(lib)

# macro for building stxxl programs

macro(stxxl_build PROGNAME)

  add_executable(${PROGNAME} ${PROGNAME}.cpp)
  target_link_libraries(${PROGNAME} ${STXXL_LIBS})

endmacro (stxxl_build)

# macro for building stxxl tests

macro(stxxl_test TESTNAME)

  if(BUILD_TESTS)
    stxxl_build(${TESTNAME})
    add_test(${TESTNAME} ${TESTNAME})
  endif(BUILD_TESTS)

endmacro (stxxl_test)

# macro for building extra stxxl tests

macro(stxxl_extra TESTNAME)

  if(BUILD_TESTS)
    stxxl_build(${TESTNAME})
  endif(BUILD_TESTS)

endmacro (stxxl_extra)

# macro for setting additional defines for targets

macro (add_define PROGNAME)

  if(TARGET ${PROGNAME})
    set_property(TARGET ${PROGNAME} APPEND PROPERTY COMPILE_DEFINITIONS ${ARGN})
  else()
    #message("Ignoring add_define(${PROGNAME} ${ARGN})")
  endif()

endmacro (add_define TESTNAME)

# build extra binaries and tests
add_subdirectory(common)
add_subdirectory(algo)
add_subdirectory(containers)
add_subdirectory(io)
add_subdirectory(mng)
add_subdirectory(stream)
add_subdirectory(utils)
