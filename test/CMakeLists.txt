cmake_minimum_required(VERSION 3.10)

add_executable(test_biotsavart test_biotsavart.f90)
add_library(test_util SHARED test_util.f90)
target_link_libraries(test_biotsavart PRIVATE biotsavart test_util)
