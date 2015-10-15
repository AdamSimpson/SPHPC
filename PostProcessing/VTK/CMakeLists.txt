cmake_minimum_required(VERSION 3.1)
project(sph_post)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/)

set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED TRUE)

file(GLOB SRC_FILES src/*.cpp)

add_executable(sph_post ${SRC_FILES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../Common)

#set(MPI_C_COMPILER mpicc)
#set(MPI_CXX_COMPILER mpic++)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph_post ${MPI_C_LIBRARIES})
target_link_libraries(sph_post ${MPI_CXX_LIBRARIES})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(sph_post ${VTK_LIBRARIES})

set(Boost_NO_BOOST_CMAKE TRUE) # Ignore system boost and use module system boost
find_package(Boost COMPONENTS mpi system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(sph_post ${Boost_LIBRARIES})

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})
target_link_libraries(sph_post ${ADIOS_LIBRARIES})

install(TARGETS sph_post DESTINATION bin)
