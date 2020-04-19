cmake_minimum_required(VERSION 3.15)
project(haptools)

set(CMAKE_CXX_STANDARD 14)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}/htslib-1.10.2 -lz -lm -lbz2 -llzma -lcurl" CACHE STRING "Set C++ Compiler Flags" FORCE)

add_executable(haptools main.cpp)

target_link_libraries(haptools ${CMAKE_SOURCE_DIR}/htslib-1.10.2/libhts.dylib)