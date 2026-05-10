cmake_minimum_required(VERSION 3.5)
project(FluidMaxwell1D CXX)

# Compatible with GCC 4.8.5 and CMake 2.8.12.x
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall")

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
else()
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

add_executable(fluid_maxwell_1d
    src/main.cpp
    src/fluid_maxwell_1d.cpp
    src/laser_driver.cpp
    src/maxwell_solver.cpp
    src/fluid_solver.cpp
    src/result_output.cpp
)
