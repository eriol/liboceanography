cmake_minimum_required(VERSION 2.8)
include_directories (. ../src/)
add_definitions(-Wall -Wextra -pedantic -ansi -std=c89)
link_libraries(check oceanography)

add_executable(check_oceanography ${COMMON} check_oceanography.c)

enable_testing()
add_test(check_oceanography check_oceanography)
