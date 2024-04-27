# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (sdcsim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake >= 3.24:
cmake_policy(SET CMP0135 NEW)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

# added include
include_directories(include)

if(DEFINED ACC_TARGET)

  message(STATUS "ACC_TARGET: ${ACC_TARGET}")

  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -acc=${ACC_TARGET} -Minfo=accel -Xcompiler -fopenmp -mp=multicore")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -acc=${ACC_TARGET} -Minfo=accel -Xcompiler -fopenmp -mp=multicore")

  if(ACC_TARGET MATCHES "gpu")
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -gpu=cc61,debug -cudalib=curand")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -gpu=cc61,fastmath -cudalib=curand")
  endif()

  # Add an executable
  add_executable(sdcsim
    src/main.cpp
    src/domain.cpp
    src/strand.cpp
    src/molecule.cpp
    src/molecule_gpu.cpp
    src/output.cpp
    src/assembly.cpp
    src/register.cpp
    src/svg.cpp
    src/nucleotides.cpp
    src/error.cpp
  )

  target_link_libraries(sdcsim argparse)
endif()

if(NOT DEFINED ACC_TARGET)
  message(STATUS "CPU only compiler. Use NVHPC compiler If you want compile GPU version. ${CMAKE_CXX_COMPILER_ID}")

  # google test
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )

  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  add_library(project_lib
    src/domain.cpp
    src/strand.cpp
    src/molecule.cpp
    src/molecule_gpu.cpp
    src/output.cpp
    src/assembly.cpp
    src/register.cpp
    src/svg.cpp
    src/nucleotides.cpp
    src/error.cpp
  )

  # Add an executable
  add_executable(sdcsim
    src/main.cpp
  )

  target_link_libraries(sdcsim project_lib)
  target_link_libraries(sdcsim argparse)

  ##############
  # Unit Tests
  ##############
  enable_testing()

  add_executable(unitTests
    test/unit/domain.cpp
    test/unit/strand.cpp
    test/unit/molecule.cpp
    test/unit/assembly.cpp
  )

  # Standard linking to gtest stuff.
  target_link_libraries(unitTests gtest gtest_main)

  # Extra linking for the project.
  target_link_libraries(unitTests project_lib)
endif()