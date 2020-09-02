cmake_minimum_required(VERSION 2.6)
project(mhaptools)

set(CMAKE_CXX_STANDARD 14)

INCLUDE_DIRECTORIES(
        include
        ./htslib-1.10.2/htslib/
)

AUX_SOURCE_DIRECTORY(
        ${CMAKE_SOURCE_DIR} DIR_SRCS
)

add_executable(mhaptools haptk.cpp convert.cpp hap.cpp merge.cpp beta.cpp summary.cpp)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0") #for test
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0") #for test
#
#set(GCC_OMIT_FRAME_POINTER "-fno-omit-frame-pointer") #for test
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_OMIT_FRAME_POINTER}" ) #for test

#set(SOURCE_FILES hap.cpp)

target_link_libraries(mhaptools ${CMAKE_SOURCE_DIR}/htslib-1.10.2/libhts.a)

target_link_libraries(mhaptools z)

target_link_libraries(mhaptools m)

target_link_libraries(mhaptools bz2)

target_link_libraries(mhaptools lzma)

target_link_libraries(mhaptools curl)
