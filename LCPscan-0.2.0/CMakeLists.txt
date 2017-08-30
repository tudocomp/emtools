cmake_minimum_required(VERSION 2.8)

project(LCPscan)

# Disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# Find OpenMP
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP not found.")
endif()

# Enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

# Include debug symbols
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g2")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g2")

# Add C++0x support
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++0x")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")

# Descend into subdirectories
add_subdirectory(src)

include_directories("${PROJECT_SOURCE_DIR}/src/stxxl/include" "${PROJECT_BINARY_DIR}/src/stxxl/include")

add_executable(construct_lcp ${PROJECT_SOURCE_DIR}/src/main.cc
                             ${PROJECT_SOURCE_DIR}/src/lcpscan_src/sort.cc
                             ${PROJECT_SOURCE_DIR}/src/lcpscan_src/utils.cc
                             ${PROJECT_SOURCE_DIR}/src/lcpscan_src/stxxl_utils.cc
)

target_link_libraries(construct_lcp stxxl)
