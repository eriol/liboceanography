cmake_minimum_required(VERSION 2.8)

project(liboceanography C)

include_directories(src)

add_definitions(-Wall -W -Wextra -Wstrict-prototypes -pedantic -ansi -std=c89)

add_subdirectory(src)
add_subdirectory(tests)

enable_testing()
