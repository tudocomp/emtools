cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
project (emtools)

include(ExternalProject)

# C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -march=native")

# Warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# cmakemodules
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# Find OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Tools
add_subdirectory(pSAscan-0.1.0)                # SA
add_subdirectory(EM-SuccinctIrreducible-0.1.0) # PLCP
add_subdirectory(LCPscan-0.2.0)                # LCP

# Find STXXL
set(STXXL_DIR "" CACHE PATH "STXXL build path")
find_package(STXXL REQUIRED)

# stxxltools
include_directories(${STXXL_INCLUDE_DIRS})
add_subdirectory(stxxltools)

# gen.sh Script
set(TDC_BIN_DIR  "" CACHE PATH "tudocomp binary directory")
if(TDC_BIN_DIR STREQUAL "")
    message("WARNING: TDC_BIN_DIR is empty, please set it to the path where tudocomp's binaries 'tdc' and 'plcp' are located")
else()
    set(TDC_BIN  ${TDC_BIN_DIR}/tdc)
    set(PLCP_BIN ${TDC_BIN_DIR}/plcp)

    configure_file(
        ${CMAKE_SOURCE_DIR}/gen.sh.in
        ${CMAKE_BINARY_DIR}/gen.sh
        @ONLY)
endif()

