add_executable(isaandbwt isaandbwt.cpp)
target_link_libraries(isaandbwt stxxl)

add_executable(readplcp readplcp.cpp)
target_link_libraries(readplcp stxxl)

add_executable(standardize standardize.cpp)
add_executable(unstandardize unstandardize.cpp)

set(STANDARDIZE_BIN ${CMAKE_CURRENT_BINARY_DIR}/standardize PARENT_SCOPE)
set(UNSTANDARDIZE_BIN ${CMAKE_CURRENT_BINARY_DIR}/unstandardize PARENT_SCOPE)
set(ISAANDBWT_BIN ${CMAKE_CURRENT_BINARY_DIR}/isaandbwt PARENT_SCOPE)
set(READPLCP_BIN ${CMAKE_CURRENT_BINARY_DIR}/readplcp PARENT_SCOPE)

