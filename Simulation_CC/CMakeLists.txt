cmake_minimum_required(VERSION 3.1)
project(sph)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../CMake/)

set(CMAKE_CXX_STANDARD 14)
set(CXX_STANDARD_REQUIRED TRUE)

file(GLOB SRC_FILES src/*.cpp)
include_directories(".")

add_executable(sph ${SRC_FILES})

# For when compiling without NVCC
add_definitions(-DHEMI_CUDA_DISABLE -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP)
MESSAGE("DISABLING CUDA")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../Common)

#set(MPI_CXX_COMPILER mpic++)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph ${MPI_CXX_LIBRARIES})

set(Boost_NO_BOOST_CMAKE TRUE) # Ignore system boost and use module system boost
find_package(Boost COMPONENTS mpi system filesystem serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(sph ${Boost_LIBRARIES})

set(HEMI_INCLUDE_DIR "external_libraries/hemi")
message("Setting HEMI_INCLUDE_DIR manually")
find_package(Hemi REQUIRED)
include_directories(${HEMI_INCLUDE_DIRS})

find_package(Thrust REQUIRED)
include_directories(${THRUST_INCLUDE_DIRS})

set(ENV{ADIOS_DIR} "/usr")
message("Setting ADIOS_DIR manually")
find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})
target_link_libraries(sph ${ADIOS_LIBRARIES})

# CLION header hack: https://youtrack.jetbrains.com/issue/CPP-270
FILE (GLOB_RECURSE clion_all_headers ${CMAKE_SOURCE_DIR}/*.h ${THRUST_INCLUDE_DIR}/thrust/*.h)
message(${clion_all_headers})
ADD_CUSTOM_TARGET(all_clion SOURCES ${clion_all_headers})

install(TARGETS sph DESTINATION bin)
