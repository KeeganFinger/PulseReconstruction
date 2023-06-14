cmake_minimum_required(VERSION 3.14)
project(PulseReconstruction)

set(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED COMPONENTS CXX HL C)

add_library(multiphotonIonization STATIC physics_functions.cpp sub_matrix_elements.cpp h5_functions.cpp h5_functions.h)

link_directories(${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

add_executable(PulseReconstruction main.cpp autocorrelation.cpp)

target_link_libraries(PulseReconstruction multiphotonIonization ${HDF5_LIBRARIES})
