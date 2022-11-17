cmake_minimum_required(VERSION 3.0)
project(NarrowbandSSB)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

find_package(PkgConfig)
pkg_check_modules(VOLK REQUIRED volk)
pkg_check_modules(FFTW3 REQUIRED fftw3f)
pkg_check_modules(SNDFILE REQUIRED sndfile)

link_directories(${VOLK_LIBRARY_DIRS})
include_directories(${VOLK_INCLUDE_DIRS})
link_libraries(${VOLK_LIBRARIES})

link_directories(${FFTW3_LIBRARY_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
link_libraries(${FFTW3_LIBRARIES})

link_directories(${SNDFILE_LIBRARY_DIRS})
include_directories(${SNDFILE_INCLUDE_DIRS})
link_libraries(${SNDFILE_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/dsp)

add_executable(modulate modulate.cpp)
add_executable(demodulate demodulate.cpp)
