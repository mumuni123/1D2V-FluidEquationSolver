cmake_minimum_required(VERSION 3.5)
project(FluidMaxwell1D CXX)

# Compatible with GCC 4.8.5 and CMake 2.8.12.x
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall")

add_executable(fluid_maxwell_1d
    src/main.cpp
    src/fluid_maxwell_1d.cpp
    src/maxwell_solver.cpp
    src/fluid_solver.cpp
    src/result_output.cpp
)
