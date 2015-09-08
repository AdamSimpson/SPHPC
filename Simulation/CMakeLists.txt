cmake_minimum_required(VERSION 3.1)
project(sph_sim)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake/)

option(USE_OPENACC "Enable OpenACC support" ON)
option(USE_ADIOS "Enable Adios file writing" OFF)

set(CMAKE_C_STANDARD 99)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DDEBUG -DTHRUST_DEBUG)
endif()

file(GLOB SRC_FILES src/*.c src/*.cc)
include_directories(../Common)
 
# Remove adios files from glob if not set
if (${USE_ADIOS})
  file(GLOB MPI_IO_REMOVE src/fileio_MPI_IO.c)
  list(REMOVE_ITEM SRC_FILES ${MPI_IO_REMOVE})
  add_definitions( -DUSE_ADIOS )
else()
  file(GLOB ADIOS_REMOVE src/fileio_adios.c)
  list(REMOVE_ITEM SRC_FILES ${ADIOS_REMOVE})
endif()

# A bit of a hack to enable easy Thrust/OpenACC compile on Titan
# Compile thrust/cuda file separately
find_package(CUDA)
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_HOST_COMPILER /usr/bin/gcc)
#set(CUDA_NVCC_FLAGS "-O3 -arch sm_35 --compiler-options -Wall")
file(GLOB NVCC_FILES src/*.cu)
# Use local thrust library as version currently on Titan has bug in lower/upper bounds
CUDA_INCLUDE_DIRECTORIES( $ENV{THRUST_INCLUDE} ${CUDA_INCLUDE_DIRS})

cuda_compile(NVCC_OBJECTS ${NVCC_FILES})

add_executable(sph_sim ${SRC_FILES} ${NVCC_OBJECTS})

find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph_sim ${MPI_LIBRARIES})

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

if (${USE_ADIOS})
find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})
target_link_libraries(sph_sim ${ADIOS_LIBRARIES})
endif()

# Machine specific settings
SITE_NAME(myhostname)
if (${myhostname} MATCHES "(.*)titan(.*)" OR ${myhostname} MATCHES "(.*)chester(.*)") # Titan(or Chester)
  find_library(M_LIB m)
  target_link_libraries(sph_sim ${M_LIB})

  if (USE_OPENACC)
    set(ACC_FLAGS "-acc -Minfo=acc -O3 -ta=tesla:cc35,cuda7.0 -Mcuda")
  endif()

  set(CMAKE_C_FLAGS "${ACC_FLAGS}")
  set(CMAKE_CXX_FLAGS "${ACC_FLAGS}")

else() # Adam's workstation
  target_link_libraries(sph_sim ${CUDA_LIBRARIES})

  if (USE_OPENACC)
    set(CMAKE_C_FLAGS "-acc -Minfo=acc -O3 -ta=tesla:cuda7.5,cc50 -Mcuda")
    set(CMAKE_CXX_FLAGS "-acc -Minfo=acc -O3 -ta=tesla:cuda7.5,cc50 -Mcuda")
  endif()
endif()

install(TARGETS sph_sim DESTINATION bin)
