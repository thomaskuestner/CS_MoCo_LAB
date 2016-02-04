cmake_minimum_required(VERSION 2.6)

project(CS_LAB)

if (WIN32)
ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif (WIN32)

###############################################################
#Bootstrap search for libraries 
# (We need to find cmake modules in Gadgetron)
###############################################################
find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
$ENV{GADGETRON_HOME}/cmake
C:/Program Files/GADGETRON/gadgetron)

if (NOT GADGETRON_CMAKE_MODULES)
  MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})
###############################################################

find_package(Ismrmrd REQUIRED)
find_package(Gadgetron REQUIRED)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
find_package(FFTW3 REQUIRED)

set(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})

INCLUDE_DIRECTORIES(${ACE_INCLUDE_DIR} 
     ${Boost_INCLUDE_DIR}
     ${GADGETRON_INCLUDE_DIR}
     ${ISMRMRD_INCLUDE_DIR}
     ${ISMRMRD_SCHEMA_DIR}
     ${ISMRMRD_XSD_INCLUDE_DIR}
	 ${FFTW3}
     )

LINK_DIRECTORIES(${GADGETRON_LIB_DIR})

ADD_LIBRARY(CS_LAB SHARED CS_LAB.cpp)

TARGET_LINK_LIBRARIES(CS_LAB 
                      cpucore 
                      optimized ${ACE_LIBRARIES} 
                      debug ${ACE_DEBUG_LIBRARY})

INSTALL (FILES CS_LAB.h
         CS_LAB_export.h
         DESTINATION include)

INSTALL(TARGETS CS_LAB DESTINATION lib)

INSTALL(FILES CS_LAB.xml DESTINATION config)