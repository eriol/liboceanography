# Top-level CMakeLists.txt for liboceanography

cmake_minimum_required(VERSION 2.8)

project(liboceanography C)

set(oceanography_version 1.0.0)

# Prevent build in source directory
if ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
    message (FATAL_ERROR "Building in the source directory is not supported.")
endif ()

include_directories(src)

add_definitions(-Wall -W -Wextra -Wstrict-prototypes -pedantic -ansi -std=c89)

add_subdirectory(src)
add_subdirectory(tests)

enable_testing()

# Packaging stuff
set(CPACK_SOURCE_PACKAGE_FILE_NAME liboceanography-${oceanography_version})
set(CPACK_SOURCE_IGNORE_FILES
    "~$"
    "^\\\\.*$"
    "^${PROJECT_SOURCE_DIR}/build/"
    "^${PROJECT_SOURCE_DIR}/docs/_build/"
)
set(CPACK_SOURCE_GENERATOR TGZ)
include(CPack)
