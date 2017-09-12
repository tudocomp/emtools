add_executable(isaandbwt isaandbwt.cpp)
target_link_libraries(isaandbwt stxxl)

add_executable(readplcp readplcp.cpp)
target_link_libraries(readplcp stxxl)

add_executable(standardize standardize.cpp)
add_executable(unstandardize unstandardize.cpp)

