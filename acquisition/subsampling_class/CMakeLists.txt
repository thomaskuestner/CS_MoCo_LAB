#---------------------------
#The code for 'subsampling' can be compiled easily.
#This CMakeLists.txt is provided just for easier on Ubuntu 14.04.
#Thanks to siatjiasen

project(subsampling)

cmake_minimum_required(VERSION 2.8)

set(include_file Approx.h LinkedList.h Point.h SomeFunctions.h VariableDensity.h VDSamplingUpper.h SubsampleMain.h)
set(src_file Approx.cpp LinkedList.cpp Point.cpp SomeFunctions.cpp VariableDensity.cpp VDSamplingUpper.cpp SubsampleMain.cpp)
set(CMAKE_INSTALL_PATH ${CMAKE_HOME_DIRECTORY})

add_executable(subsampling ${src_file} ${include_file})

install(TARGETS subsampling RUNTIME DESTINATION ${CMAKE_INSTALL_PATH})
