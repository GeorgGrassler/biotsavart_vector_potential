cmake_minimum_required(VERSION 3.10)
project(BiotSavart Fortran)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_Fortran_COMPILER gfortran)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

include_directories(${SRC_DIR})
add_library(biotsavart SHARED ${SRC_DIR}/biotsavart.f90)

# Add the test program
add_executable(test_biotsavart ${SRC_DIR}/test_biotsavart.f90)
add_library(test_util SHARED ${SRC_DIR}/test_util.f90)
target_link_libraries(test_biotsavart PRIVATE biotsavart test_util)
target_include_directories(test_biotsavart PRIVATE ${SRC_DIR})