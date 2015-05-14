cmake_minimum_required(VERSION 2.8)
project(sph)

# Not ideal but works
set(CMAKE_C_FLAGS "-std=c99")

file(GLOB SRC_FILES src/*.c src/*.cpp)

add_executable(sph ${SRC_FILES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph ${MPI_LIBRARIES})

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

find_library(M_LIB m)
target_link_libraries(sph ${M_LIB})

install(TARGETS sph DESTINATION bin)

