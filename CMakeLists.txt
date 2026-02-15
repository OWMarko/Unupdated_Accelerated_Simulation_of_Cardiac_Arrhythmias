cmake_minimum_required(VERSION 3.15)
project(CardiacSim-FEM VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -O3 -march=native)

_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include)

set(SOURCES
    src/solvers/AlievPanfilov.cpp
    src/solvers/ExplicitEuler.cpp
    src/solvers/RK4.cpp
    src/main.cpp
)

add_executable(CardiacSim ${SOURCES})


find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(CardiacSim PUBLIC OpenMP::OpenMP_CXX)
endif()