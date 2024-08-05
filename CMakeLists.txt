cmake_minimum_required(VERSION 3.10)
project(BiotSavart Fortran)

set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} \
                            -Wall \
                            -Wextra \
                            -Werror \
                            -fbacktrace \
                            -fbounds-check \
                            -fcheck=all \
                            -g \
                            -O0"
                        )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

add_subdirectory(${CMAKE_SOURCE_DIR}/test)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(${CMAKE_SOURCE_DIR}/src)
