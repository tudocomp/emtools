ExternalProject_Add(
    EM-SuccinctIrreducible
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src
    PREFIX "EM-SuccinctIrreducible-0.1.0"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "make"
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
)

set(CONSTRUCT_PLCP_SEQ_BIN
    ${CMAKE_CURRENT_SOURCE_DIR}/src/construct_plcp_sequential
    PARENT_SCOPE)

