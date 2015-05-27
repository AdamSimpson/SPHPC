cmake_minimum_required(VERSION 3.1)
project(sph)

set(CMAKE_C_STANDARD 99)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DDEBUG=1)
endif()

file(GLOB SRC_FILES src/*.c src/*.cc)

add_executable(sph ${SRC_FILES})

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

install(TARGETS sph DESTINATION bin)
