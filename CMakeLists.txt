cmake_minimum_required(VERSION 3.14)
project(CardiacSimFEM VERSION 1.0 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(MSVC)
    add_compile_options(/O2 /fp:fast)
else()
    add_compile_options(-O3 -march=native -ffast-math -Wall -Wextra)
endif()


find_package(OpenMP REQUIRED)


include_directories(${CMAKE_SOURCE_DIR}/external/eigen)


include_directories(include)


file(GLOB_RECURSE SOURCES "src/*.cpp")


add_executable(CardiacSim ${SOURCES})


target_link_libraries(CardiacSim PUBLIC OpenMP::OpenMP_CXX)


# add_subdirectory(tests)
