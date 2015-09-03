link_libraries(check oceanography pthread rt subunit)

add_executable(test_oceanography test_oceanography.c)

enable_testing()
add_test(test_oceanography test_oceanography)
