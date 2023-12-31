cmake_minimum_required(VERSION 3.10)

project(pe827)

set(CMAKE_CXX_STANDARD 11)

add_executable(pe827 main.cpp)

target_link_libraries(pe827 PRIVATE gmp)
