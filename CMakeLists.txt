cmake_minimum_required(VERSION 3.15)
project(haptools)

set(CMAKE_CXX_STANDARD 14)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}/htslib-1.10.2 " CACHE STRING "Set C++ Compiler Flags" FORCE)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0") #for test
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0") #for test
#
#set(GCC_OMIT_FRAME_POINTER "-fno-omit-frame-pointer") #for test
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_OMIT_FRAME_POINTER}" ) #for test

add_executable(haptools haptk.cpp)

set(SOURCE_FILES convert.cpp)

target_link_libraries(haptools ${CMAKE_SOURCE_DIR}/htslib-1.10.2/libhts.a)

target_link_libraries(haptools z)

target_link_libraries(haptools m)

target_link_libraries(haptools bz2)

target_link_libraries(haptools lzma)

target_link_libraries(haptools curl)
