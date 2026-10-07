cmake_minimum_required(VERSION 2.6)

project(mhaptools C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")

INCLUDE_DIRECTORIES(
        ./include
        ./htslib-1.10.2/htslib
        ./htslib-1.10.2
)

AUX_SOURCE_DIRECTORY(
        ${CMAKE_SOURCE_DIR} DIR_SRCS
)

add_executable(mhaptools haptk.cpp convert.c bsread.c mhap.cpp merge.cpp beta.cpp summary.cpp utils.cpp)

# libhts.so or libhts.dylib built in htslib-1.10.2 (see README); the shared
# library brings the codecs it was configured with (bzip2, lzma, curl)
find_library(HTS_LIB hts PATHS ${CMAKE_SOURCE_DIR}/htslib-1.10.2 NO_DEFAULT_PATH)

target_link_libraries(mhaptools ${HTS_LIB})

target_link_libraries(mhaptools z)

target_link_libraries(mhaptools m)

target_link_libraries(mhaptools pthread)
