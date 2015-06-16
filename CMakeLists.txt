cmake_minimum_required(VERSION 3.1)
project(sph)

option(OPENACC "Enable OpenACC support" ON)

set(CMAKE_C_STANDARD 99)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DDEBUG -DTHRUST_DEBUG)
endif()

file(GLOB SRC_FILES src/*.c src/*.cc)

# A bit of a hack to enable easy Thrust/OpenACC compile on Titan
# Compile thrust/cuda file seperatly
find_package(CUDA)
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_HOST_COMPILER /usr/bin/gcc)
file(GLOB NVCC_FILES src/*.cu)
cuda_compile(NVCC_OBJECTS ${NVCC_FILES})

add_executable(sph ${SRC_FILES} ${NVCC_OBJECTS})

find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph ${MPI_LIBRARIES})

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

find_library(M_LIB m)
target_link_libraries(sph ${M_LIB})

if (OPENACC)
  set(CMAKE_C_FLAGS "-acc -Minfo")
  set(CMAKE_CXX_FLAGS "-acc -Minfo")
endif()

install(TARGETS sph DESTINATION bin)
