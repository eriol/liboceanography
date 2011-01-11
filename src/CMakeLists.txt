cmake_minimum_required(VERSION 2.8)
include_directories(.)
add_definitions(-Wall -Wextra -pedantic -ansi -std=c89)
link_libraries(m)

add_library(oceanography SHARED oceanography.c)
