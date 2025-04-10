cmake_minimum_required(VERSION 3.18)

project(ATSP)

include_directories(include)

add_executable(atsp
    src/adding_edges.cpp)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
