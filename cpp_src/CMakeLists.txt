cmake_minimum_required(VERSION 3.25)
project(PulseReconstruction)

set(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED COMPONENTS CXX)

if(MINGW)
    # See if GCC's -gdwarf-4 is supported
    # See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101377 for why this is needed on Windows
#    if(HAVE_GDWARF_4)
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-4")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-4")
#    endif()
endif()

add_library(multiphotonIonization STATIC physics_functions.cpp sub_matrix_elements.cpp)

add_executable(PulseReconstruction main.cpp autocorrelation.cpp)

target_link_libraries(PulseReconstruction multiphotonIonization hdf5::hdf5_cpp)