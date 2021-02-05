cmake_minimum_required(VERSION 3.11)
project(mhd_transfer)

include_directories("/usr/include/eigen3", "include")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lCGAL -lz -std=c++14 -fopenmp -march=native -DNDEBUG -O3 -fext-numeric-literals")
set(SOURCE_FILES src/main.cpp src/utils.cpp src/Ray.cpp src/Image.cpp
        src/ImagePlane.cpp src/Geometry.cpp include/MyExceptions.h
        src/Intersection.cpp src/BField.cpp
        src/SimulationInterpolater.cpp src/NField.cpp src/VField.cpp src/Jet.cpp
        src/Observation.cpp src/System.cpp)
add_executable(mhd_transfer ${SOURCE_FILES})


# Require MPI for this project:
find_package(MPI REQUIRED)


set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
message("${MPI_INCLUDE_PATH}")
target_link_libraries(mhd_transfer ${MPI_LIBRARIES})