cmake_minimum_required(VERSION 2.8)
include_directories (. ../src/)
add_definitions(-Wall -Wextra -pedantic -ansi -std=c89)
link_libraries(check oceanography)

add_executable(test_oceanography ${COMMON} test_oceanography.c)

enable_testing()
add_test(test_oceanography test_oceanography)
