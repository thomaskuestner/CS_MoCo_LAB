cmake_minimum_required(VERSION 2.6)
project(CS_LAB)

#Set the build type to Release if not specified
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# build options for 64 bits system
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("64bit system is found")
  set( HAS_64_BIT On CACHE BOOL "64bit build")
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("32bit system is found")
  set( HAS_64_BIT Off CACHE BOOL "64bit build")
endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# build options for OpenMP support
find_package(OpenMP)
option(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if(USE_OPENMP) 
        message("OpenMP multithreading enabled")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSE_OMP)
    else (USE_OPENMP)
        message("OpenMP multithreading is supported, but disabled")
    endif(USE_OPENMP) 
else (OPENMP_FOUND)
  message("OpenMP multithreading not supported")
endif (OPENMP_FOUND)

if (WIN32)
    add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS)
#    ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_VARIADIC_MAX=10) #to fix compiler limitations in Visual Studio Express
    if ( HAS_64_BIT )
        add_definitions(-DWIN64 -D_WIN64)
    endif ( HAS_64_BIT )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
    set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL:NO")
    set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} /INCREMENTAL:NO")
    set (CMAKE_EXE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set (CMAKE_SHARED_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set (CMAKE_STATIC_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set (CMAKE_MODULE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    # The two flags below is to fix Windows problems in relation to multiple defined operators new/delete and some constructors that are defined in our headers
    #SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE") 
    #SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE") 
else (WIN32)
    if (UNIX)
        if (APPLE)
        else (APPLE)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libgfortran")
        endif (APPLE)
    endif (UNIX)
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/CS_LAB)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# whether to suppress compilation warnings
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    if (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    elseif (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
    endif (WIN32)
endif (BUILD_SUPPRESS_WARNINGS)

# whether to compile toolboxes as static library
option(BUILD_TOOLBOX_STATIC "Build static library for toolboxes" Off)

if ( BUILD_TOOLBOX_STATIC )
    message("Building static toolbox libray")
    add_definitions(-DBUILD_TOOLBOX_STATIC)
    set(LIBTYPE STATIC)
    if ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
    endif ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
else ( BUILD_TOOLBOX_STATIC )
    message("Building dynamic toolbox libray")
    set(LIBTYPE SHARED)
endif ( BUILD_TOOLBOX_STATIC )

###############################################################
# BOOST
###############################################################
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# necessary for Windows and RHEL <=6 systems
set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif(WIN32)

if(WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif(WIN32)

# We actually only use system and thread explicitly, but they require linking in date_time and chrono
if (WIN32)
  find_package(Boost COMPONENTS system thread date_time chrono program_options filesystem REQUIRED)
else(WIN32)
  find_package(Boost COMPONENTS system thread program_options filesystem REQUIRED)
endif(WIN32)

###############################################################
# ACE
###############################################################
find_package(ACE)
if(NOT ACE_FOUND)
  message("ACE not found.")
endif(NOT ACE_FOUND)

###############################################################
# CUDA
###############################################################
find_package(CUDA_advanced)
if (CUDA_FOUND)
  add_definitions(-DUSE_CUDA)
  message("CUDA found, GPU components are compiled if CULA is also found")
  set( GADGETRON_CUDA_FOUND_BOOL 1 )
  include_directories( ${CUDA_INCLUDE_DIRS} )
  #set(CUDA_VERBOSE_BUILD ON)

  # Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
  # Support compute model 3.0 from Cuda 4.2 and up
  # Support compute model 3.5 from Cuda 5 and up
  message("Compiling with ${CUDA_NVCC_FLAGS}")

else (CUDA_FOUND)
  message("CUDA not found. CUDA components will not be compiled.")
  set( GADGETRON_CUDA_FOUND_BOOL 0 )
endif (CUDA_FOUND)

###############################################################
# Armadillo
###############################################################
find_package(Armadillo)
# check whether ILP64 MKL should is used
if(ARMADILLO_FOUND)
    add_definitions(-DUSE_ARMADILLO)
    set(ARMADILLO_BLAS_LONG_LONG FALSE)
    if(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp")
        # Read and parse armadillo config.hpp to find out whether BLAS uses long long
        file(STRINGS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp" _armadillo_blas_long_long REGEX "// #define ARMA_BLAS_LONG_LONG")
        if ( NOT _armadillo_blas_long_long )
            set(ARMADILLO_BLAS_LONG_LONG TRUE)
            message("Armadillo is found to use long long for BLAS calls")
        else ( NOT _armadillo_blas_long_long )
            message("Armadillo found. Note that ARMADILLO_BLAS_LONG_LONG _must_ be defined in ${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp to link against the MKL ILP64 interface.")
        endif ( NOT _armadillo_blas_long_long )
        unset(_armadillo_blas_long_long)
    endif()
else()
  message("Armadillo not found. This will disable many toolboxes and gadgets.")
endif()

###############################################################
# HDF
###############################################################
find_package(HDF5 1.8 COMPONENTS C CXX HL)
if (HDF5_FOUND)
message("HDF5 Found")
else()
message("HDF5 not found")
endif()

###############################################################
# FFTW
###############################################################
find_package(FFTW3 COMPONENTS single double)
if (FFTW3_FOUND)
message("FFTW3 Found")
else()
message("FFTW3 not found")
endif()

###############################################################
# ISMRMRD
###############################################################
find_package(Ismrmrd)
if (ISMRMRD_FOUND)
message("ISMRMRD Found")
else()
message("ISMRMRD not found")
endif()

###############################################################
# XSD/XersesC
###############################################################
find_package(XSD)
find_package(XercesC)

###############################################################
# MATLAB
###############################################################
find_package(GMatlab)
if (MATLAB_FOUND)
    add_definitions(-DUSE_MATLAB)
    include_directories(${MATLAB_INCLUDE_DIR})
    target_link_libraries(CS_LAB ${MATLAB_MEX_LIBRARY})
endif (MATLAB_FOUND)

if (ACE_FOUND AND ISMRMRD_FOUND AND XSD_FOUND AND XERCESC_FOUND AND FFTW3_FOUND AND HDF5_FOUND)
else()
  message("Required dependencies for gadget compilation not found (ACE, ISMRMRD, XSD, XERCESC, FFTW3, HDF5).")
endif()

add_subdirectory(cmake)

IF (WIN32)
    add_definitions(-D__BUILD_CS_LAB__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(Gadgetron REQUIRED)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
find_package(FFTW3 REQUIRED)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${ARMADILLO_INCLUDE_DIRS}
    ${ISMRMRD_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src/FOCUSS
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${GADGETRON_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${ISMRMRD_SCHEMA_DIR}
    ${ISMRMRD_XSD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${XSD_INCLUDE_DIR}
    ${XERCESC_INCLUDE_DIR}
)

set( cs_lab_header_files 	src/FOCUSS/CS_Control.h 
                            src/FOCUSS/CS_FOCUSS.h 
                            src/FOCUSS/GlobalVar_FOCUSS.h 
                            src/FOCUSS/hoNDDCT.h 
                            src/FOCUSS/SomeFunctions.h 
                            src/FOCUSS/SomeFunctions.hxx 
                            src/FOCUSS/Transform.h 
							src/FOCUSS/TransformWrapper.h
							CS_LAB.h
							 )

set( cs_lab_src_files		src/FOCUSS/CS_Control.cpp 
                            src/FOCUSS/CS_FOCUSS.cpp 
                            src/FOCUSS/CS_FOCUSS_2D.cpp 
                            src/FOCUSS/CS_FOCUSS_2Dt.cpp 
							src/FOCUSS/CS_FOCUSS_3D.cpp
							CS_LAB.cpp
							src/FOCUSS/GlobalVar_FOCUSS.cpp
							src/FOCUSS/hoNDDCT.cpp
							src/FOCUSS/Transform.cpp
							)

set( cs_lab_config_files CS_LAB.xml)

add_library(CS_LAB SHARED 
    CS_LAB_export.h 
    ${cs_lab_header_files} 
    ${cs_lab_src_files}
    ${cs_lab_config_files}
    ${ISMRMRD_XSD_SOURCE}
)

if(WIN32)
	set(ADDITIONAL_LIBRARIES  
						${GADGETRON_HOME}/lib/cpucore.lib
						${GADGETRON_HOME}/lib/cpufft.lib
						${GADGETRON_HOME}/lib/cpucore_math.lib
						${Boost_LIBRARY_DIR}/boost_system-vc100-mt-1_55.lib
						${Boost_LIBRARY_DIR}/boost_thread-vc100-mt-1_55.lib
						${Boost_LIBRARY_DIR}/boost_date_time-vc100-mt-1_55.lib
						${Boost_LIBRARY_DIR}/boost_chrono-vc100-mt-1_55.lib
						${Boost_LIBRARY_DIR}/boost_program_options-vc100-mt-1_55.lib
						${Boost_LIBRARY_DIR}/boost_filesystem-vc100-mt-1_55.lib
		)
else()

endif()

target_link_libraries(	CS_LAB					
						${ISMRMRD_LIBRARIES} 
						${FFTW3_LIBRARIES} 
						${BOOST_LIBRARIES}
						${XERCESC_LIBRARIES} 
						${ACE_LIBRARY}
						${ADDITIONAL_LIBRARIES}
)



install (FILES 
         CS_LAB_export.h
		 ${cs_lab_header_files}
         DESTINATION include)

install(TARGETS CS_LAB DESTINATION lib)

install(FILES CS_LAB.xml DESTINATION config)

link_directories(${GADGETRON_LIB_DIR})

#add_library(CS_LAB SHARED CS_LAB.cpp)
