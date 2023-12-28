cmake_minimum_required(VERSION 3.12)
project(lab6 VERSION 1.0.0)

add_library(solver math/solver.cpp math/solver.hpp)
add_library(matrix math/matrix.cpp math/matrix.hpp)
add_library(lab6 labs/lab6.cpp labs/lab6.hpp)

add_executable(prog labs/main.cpp)

target_link_libraries(prog PRIVATE lab6 matrix solver)