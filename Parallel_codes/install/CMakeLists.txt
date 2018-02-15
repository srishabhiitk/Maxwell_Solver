#Cmake setup
cmake_minimum_required (VERSION 2.8.3)
project (maxwell)

set (CMAKE_CXX_COMPILER mpic++)
get_filename_component (INSTALL_DIR ${CMAKE_BINARY_DIR} PATH)
get_filename_component (PARENT_DIR ${INSTALL_DIR} PATH)

set (HEADER_DIR ${PARENT_DIR}/lib)
set (SOURCE_DIR ${PARENT_DIR}/src)

include_directories(${HEADER_DIR})

find_package(OpenMP)
message(status ${OpenMP_C_FLAGS})

SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -fopenmp")

#Adding sources
set (SOURCES ${SOURCE_DIR}/test.cc
             ${SOURCE_DIR}/reader.cc
             ${SOURCE_DIR}/grid.cc
             ${SOURCE_DIR}/parallel.cc
             ${SOURCE_DIR}/mpidata.cc
             ${SOURCE_DIR}/field.cc
             ${SOURCE_DIR}/writer.cc
             ${SOURCE_DIR}/vfield.cc)
#Compile
set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})
add_executable(maxwell ${SOURCES})
#set_target_properties(maxwell PROPERTIES COMPILE_FLAGS "-Wall -g -DDEBUG_RUN -DBZ_DEBUG")
set_target_properties(maxwell PROPERTIES COMPILE_FLAGS "-O3 -std=c++11 -O2 -Wall -pedantic -fopenmp")
set( CMAKE_CXX_FLAGS "-Wall -fopenmp -O3" )
#set_target_properties(maxwell PROPERTIES COMPILE_FLAGS "-Wall -O3")
#target_link_libraries(maxwell yaml-cpp /usr/local/lib/libblitz.a)
target_link_libraries(maxwell yaml-cpp hdf5)
