cmake_minimum_required(VERSION 2.8)
project(sph_post)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../CMake/)

set(CMAKE_C_STANDARD 99)

file(GLOB SRC_FILES src/*.c src/*.cc)

add_executable(sph_post ${SRC_FILES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../Common)

find_package(OpenVDB REQUIRED)
include_directories(${OPENVDB_INCLUDE_DIRS})
target_link_libraries(sph_post ${OPENVDB_LIBRARIES})

find_package(Boost 1.58.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
Message("${Boost_INCLUDE_DIRS}")

install(TARGETS sph_post DESTINATION bin)
