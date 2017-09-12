cmake_minimum_required(VERSION 2.8)

# Descend into subdirectories
add_subdirectory(src)

add_executable(construct_lcp ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/lcpscan_src/sort.cc
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/lcpscan_src/utils.cc
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/lcpscan_src/stxxl_utils.cc
)

target_include_directories(construct_lcp PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/stxxl/include"
    "${CMAKE_CURRENT_BINARY_DIR}/src/stxxl/include")
target_link_libraries(construct_lcp stxxl-lcpscan)

