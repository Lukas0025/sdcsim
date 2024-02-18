# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (sdcsim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# google test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# added include
include_directories(include)

add_library(project_lib
    src/domain.cpp
    src/strand.cpp
    src/molecule.cpp
    src/output.cpp
    src/assembly.cpp
    src/register.cpp
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