cmake_minimum_required(VERSION 3.14)
project(NBodySimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to release build for performance
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Treat warnings as errors, enable extensive warnings to show C++ rigor
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Werror")

# Release flags for performance optimization (targeting Apple Silicon/ARM)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcpu=apple-m1 -DNDEBUG")

# OpenMP Setup for Apple Clang via Homebrew
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES omp)
    set(OpenMP_omp_LIBRARY /opt/homebrew/opt/libomp/lib/libomp.dylib)
    include_directories(/opt/homebrew/opt/libomp/include)
endif()

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(include)

# Core Simulation Library
add_library(nbody_core src/nbody.cpp src/quad_tree.cpp)
target_link_libraries(nbody_core PRIVATE OpenMP::OpenMP_CXX)

# Main Executable
add_executable(nbody_sim src/main.cpp)
target_link_libraries(nbody_sim PRIVATE nbody_core)

# GoogleTest setup via FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.12.1
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(nbody_tests tests/test_nbody.cpp)
target_link_libraries(nbody_tests PRIVATE nbody_core GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(nbody_tests)
