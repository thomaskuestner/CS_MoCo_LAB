# add gadgetron flags
set(__GADGETRON__ 3)
set(__GADGETRON_MINOR__ 8)
set(__GADGETRON_PATCHLEVEL__ 2)

add_definitions(-D__GADGETRON_VERSION_HIGHER_3_6__=0)
IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
  add_definitions(-D__GADGETRON_VERSION_HIGHER_3_6__=1)
ENDIF()

cmake_minimum_required(VERSION 2.6)
project(CS_LAB)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

#Set the build type to Release if not specified
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# build options for 64 bits system
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("64bit system is found")
  set( HAS_64_BIT On CACHE BOOL "64bit build")
ELSE( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("32bit system is found")
  set( HAS_64_BIT Off CACHE BOOL "64bit build")
ENDIF( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# build options for OpenMP support
find_package(OpenMP)
option(USE_OPENMP "Use OpenMP" On)
IF (OPENMP_FOUND)
    IF(USE_OPENMP)
        message("OpenMP multithreading enabled")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSE_OMP)
    ELSE (USE_OPENMP)
        message("OpenMP multithreading is supported, but disabled")
    ENDIF(USE_OPENMP)
ELSE (OPENMP_FOUND)
  message("OpenMP multithreading not supported")
ENDIF (OPENMP_FOUND)

IF (WIN32)
    add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS)
#    ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_VARIADIC_MAX=10) #to fix compiler limitations in Visual Studio Express
    IF ( HAS_64_BIT )
        add_definitions(-DWIN64 -D_WIN64)
    ENDIF ( HAS_64_BIT )
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
ELSE (WIN32)
    IF (UNIX)
        IF (APPLE)
        ELSE (APPLE)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libgfortran")
        ENDIF (APPLE)
    ENDIF (UNIX)
ENDIF (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/CS_LAB)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# whether to suppress compilation warnings
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
IF (BUILD_SUPPRESS_WARNINGS)
    IF (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    ELSEIF (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
    ENDIF (WIN32)
ENDIF (BUILD_SUPPRESS_WARNINGS)

# whether to compile toolboxes as static library
option(BUILD_TOOLBOX_STATIC "Build static library for toolboxes" Off)

IF ( BUILD_TOOLBOX_STATIC )
    message("Building static toolbox libray")
    add_definitions(-DBUILD_TOOLBOX_STATIC)
    set(LIBTYPE STATIC)
    IF ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
    ENDIF ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
ELSE ( BUILD_TOOLBOX_STATIC )
    message("Building dynamic toolbox libray")
    set(LIBTYPE SHARED)
ENDIF ( BUILD_TOOLBOX_STATIC )

###############################################################
# BOOST
###############################################################
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# necessary for Windows and RHEL <=6 systems
set(Boost_NO_BOOST_CMAKE ON)

IF(WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
ENDIF(WIN32)

IF(WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
ENDIF(WIN32)

# We actually only use system and thread explicitly, but they require linking in date_time and chrono
IF (WIN32)
  find_package(Boost COMPONENTS system thread date_time chrono program_options filesystem REQUIRED)
ELSE(WIN32)
  find_package(Boost COMPONENTS system thread program_options filesystem REQUIRED)
ENDIF(WIN32)

###############################################################
# ACE
###############################################################
find_package(ACE)
IF(NOT ACE_FOUND)
  message("ACE not found.")
ENDIF(NOT ACE_FOUND)

###############################################################
# CUDA
###############################################################
find_package(CUDA_advanced)
IF (CUDA_FOUND)
  add_definitions(-DUSE_CUDA)
  message("CUDA found, GPU components are compiled if CULA is also found")
  set( GADGETRON_CUDA_FOUND_BOOL 1 )
  include_directories( ${CUDA_INCLUDE_DIRS} )
  #set(CUDA_VERBOSE_BUILD ON)

  # Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
  # Support compute model 3.0 from Cuda 4.2 and up
  # Support compute model 3.5 from Cuda 5 and up
  message("Compiling with ${CUDA_NVCC_FLAGS}")

ELSE (CUDA_FOUND)
  message("CUDA not found. CUDA components will not be compiled.")
  set( GADGETRON_CUDA_FOUND_BOOL 0 )
ENDIF (CUDA_FOUND)

###############################################################
# Armadillo
###############################################################
find_package(Armadillo)
# check whether ILP64 MKL should is used
IF(ARMADILLO_FOUND)
    add_definitions(-DUSE_ARMADILLO)
    set(ARMADILLO_BLAS_LONG_LONG FALSE)
    IF(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp")
        # Read and parse armadillo config.hpp to find out whether BLAS uses long long
        file(STRINGS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp" _armadillo_blas_long_long REGEX "// #define ARMA_BLAS_LONG_LONG")
        IF ( NOT _armadillo_blas_long_long )
            set(ARMADILLO_BLAS_LONG_LONG TRUE)
            message("Armadillo is found to use long long for BLAS calls")
        ELSE ( NOT _armadillo_blas_long_long )
            message("Armadillo found. Note that ARMADILLO_BLAS_LONG_LONG _must_ be defined in ${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp to link against the MKL ILP64 interface.")
        ENDIF ( NOT _armadillo_blas_long_long )
        unset(_armadillo_blas_long_long)
    ENDIF()
ELSE()
  message("Armadillo not found. This will disable many toolboxes and gadgets.")
ENDIF()

###############################################################
# HDF
###############################################################
find_package(HDF5 1.8 COMPONENTS C CXX HL)
IF (HDF5_FOUND)
message("HDF5 Found")
ELSE()
message("HDF5 not found")
ENDIF()

###############################################################
# FFTW
###############################################################
find_package(FFTW3 COMPONENTS single double)
IF (FFTW3_FOUND)
message("FFTW3 Found")
ELSE()
message("FFTW3 not found")
ENDIF()

###############################################################
# ISMRMRD
###############################################################
find_package(Ismrmrd)
IF (ISMRMRD_FOUND)
message("ISMRMRD Found")
ELSE()
message("ISMRMRD not found")
ENDIF()

###############################################################
# XSD/XersesC
###############################################################
IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	
ELSE()
	IF (WIN32)
		find_package(XSD)
		find_package(XercesC)
	ENDIF()	
ENDIF()


###############################################################
# MATLAB
###############################################################
find_package(GMatlab)
IF (MATLAB_FOUND)
    add_definitions(-DUSE_MATLAB)
ENDIF (MATLAB_FOUND)

IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	IF (ACE_FOUND AND ISMRMRD_FOUND AND FFTW3_FOUND AND HDF5_FOUND)
	ELSE()
	  message("Required dependencies for gadget compilation not found (ACE, ISMRMRD, FFTW3, HDF5).")
	ENDIF()
ELSE()
	IF (ACE_FOUND AND ISMRMRD_FOUND AND XSD_FOUND AND XERCESC_FOUND AND FFTW3_FOUND AND HDF5_FOUND)
	ELSE()
	  message("Required dependencies for gadget compilation not found (ACE, ISMRMRD, XSD, XERCESC, FFTW3, HDF5).")
	ENDIF()
ENDIF()

add_subdirectory(cmake)

IF (WIN32)
    add_definitions(-D__BUILD_CS_LAB__)
ENDIF (WIN32)

IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	find_package(Gadgetron REQUIRED)
ELSE()
	find_package(GadgetronOld REQUIRED)
ENDIF()

find_package(Ismrmrd REQUIRED)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
find_package(FFTW3 REQUIRED)

find_package(ITK REQUIRED)
INCLUDE( ${ITK_USE_FILE} )

###############################################################
# Elastix
###############################################################
# include elastix files and library
SET(Elastix_DIR /opt/code/elastix/build)	# THIS LINE IS MACHINE DEPENDENT. YOU HAVE TO ADAPT THE DIRECTORY!
SET( ELASTIX_USE_FILE ${ELASTIX_BINARY_DIR}/UseElastix.cmake )
find_package(Elastix REQUIRED)

MESSAGE( STATUS "Including Elastix settings.")
INCLUDE( ${ELASTIX_USE_FILE} )

IF (WIN32)
	link_libraries(param)
	link_libraries(elastix)
	link_libraries(transformix)
ENDIF()

IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	include_directories(
		${CMAKE_SOURCE_DIR}/include
		${CMAKE_SOURCE_DIR}/include/new
		${HDF5_INCLUDE_DIR}
		${HDF5_INCLUDE_DIR}/cpp
		${ARMADILLO_INCLUDE_DIRS}
		${CMAKE_SOURCE_DIR}/src/FOCUSS
		${CMAKE_SOURCE_DIR}/src/GADGET_TOOLS
		${ACE_INCLUDE_DIR}
		${Boost_INCLUDE_DIR}
		${GADGETRON_INCLUDE_DIR}
		${ISMRMRD_INCLUDE_DIR}
		${ISMRMRD_INCLUDE_DIR}/ismrmrd
		${ISMRMRD_SCHEMA_DIR}
		${ISMRMRD_XSD_INCLUDE_DIR}
		${FFTW3_INCLUDE_DIR}
	)
ELSE()
	if (WIN32)
		include_directories(
			${CMAKE_SOURCE_DIR}/include
			${CMAKE_SOURCE_DIR}/include/old
			${HDF5_INCLUDE_DIR}
			${HDF5_INCLUDE_DIR}/cpp
			${ARMADILLO_INCLUDE_DIRS}
			${CMAKE_SOURCE_DIR}/src/FOCUSS
			${CMAKE_SOURCE_DIR}/src/GADGET_TOOLS
			${ACE_INCLUDE_DIR}
			${Boost_INCLUDE_DIR}
			${GADGETRON_INCLUDE_DIR}
			${ISMRMRD_INCLUDE_DIR}
			${ISMRMRD_INCLUDE_DIR}/ismrmrd
			${ISMRMRD_SCHEMA_DIR}
			${ISMRMRD_XSD_INCLUDE_DIR}
			${FFTW3_INCLUDE_DIR}
			${XSD_INCLUDE_DIR}
			${XERCESC_INCLUDE_DIR}
		)
	ENDIF()
ENDIF()

IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))

	set(cs_lab_header_files CS_LAB.h
		src/FOCUSS/CS_Control.h
		src/FOCUSS/CS_FOCUSS.h
		src/FOCUSS/Transform.h
		src/FOCUSS/TransformWrapper.h
		src/FOCUSS/GlobalVar.h
		src/FOCUSS/SomeFunctions.h
		src/FOCUSS/SomeFunctions.hxx
		src/GADGET_TOOLS/AcquisitionFinishGadget.h
		src/GADGET_TOOLS/CS_AccumulatorGadget.h
		src/GADGET_TOOLS/CS_MRIImageToAcquisitionGadget.h
		src/GADGET_TOOLS/SlicerGadget.h
		src/GADGET_TOOLS/CS_CombineGadget.h
		src/GADGET_TOOLS/hoNDFFT_CS.h
		src/GADGET_TOOLS/hoNDKLT_CS.h
		src/GADGET_TOOLS/hoNDDCT.h
# 		src/REGISTRATION/RegistrationGadget.h
		src/LAP_GADGET/LAPRegistrationGadget.h
		src/LAP_GADGET/lap3d.h
		src/LAP_GADGET/gaussianfilterbasis.h
		src/LAP_GADGET/mynamespaces.h
		src/LAP_GADGET/shiftengine3d.h
		src/GADGET_TOOLS/ImageAccumulatorGadget.h
		src/GADGET_TOOLS/ImageSlicerGadget.h
		src/RETRO/CS_Retro_AccumulatorGadget.h
		src/RETRO/CS_Retro_NavigatorGadget.h
		src/RETRO/CS_Retro_PopulationGadget.h
		include/new/hoNDArray_math_util.h
		include/new/GadgetronCommon.h
		include/new/hoNDArray_operators.h
		include/new/hoNDArray_blas.h
		include/new/hoMatrix_util.h
	)

	set(cs_lab_src_files CS_LAB.cpp
		src/FOCUSS/CS_Control.cpp
		src/FOCUSS/CS_FOCUSS.cpp
		src/FOCUSS/CS_FOCUSS_2D.cpp
		src/FOCUSS/CS_FOCUSS_2Dt.cpp
		src/FOCUSS/CS_FOCUSS_3D.cpp
		src/FOCUSS/CS_FOCUSS_4D.cpp
		src/FOCUSS/Transform.cpp
		src/FOCUSS/GlobalVar.cpp
		src/GADGET_TOOLS/AcquisitionFinishGadget.cpp
		src/GADGET_TOOLS/CS_AccumulatorGadget.cpp
		src/GADGET_TOOLS/CS_MRIImageToAcquisitionGadget.cpp
		src/GADGET_TOOLS/SlicerGadget.cpp
		src/GADGET_TOOLS/CS_CombineGadget.cpp
		src/GADGET_TOOLS/hoNDFFT_CS.cpp
		src/GADGET_TOOLS/hoNDKLT_CS.cpp
		src/GADGET_TOOLS/hoNDDCT.cpp
# 		src/REGISTRATION/RegistrationGadget.cpp
		src/LAP_GADGET/LAPRegistrationGadget.cpp
		src/LAP_GADGET/lap3d.cpp
		src/LAP_GADGET/gaussianfilterbasis.cpp
		src/LAP_GADGET/mynamespaces.cpp
		src/LAP_GADGET/shiftengine3d.cpp
		src/GADGET_TOOLS/ImageAccumulatorGadget.cpp
		src/GADGET_TOOLS/ImageSlicerGadget.cpp
		src/RETRO/CS_Retro_AccumulatorGadget.cpp
		src/RETRO/CS_Retro_NavigatorGadget.cpp
		src/RETRO/CS_Retro_PopulationGadget.cpp
		include/new/hoNDArray_math_util.cpp
		include/new/hoNDArray_operators.cpp
	)

ELSE()
	set(cs_lab_header_files CS_LAB.h
		src/FOCUSS/CS_Control.h
        src/FOCUSS/CS_FOCUSS.h
		src/FOCUSS/Transform.h
		src/FOCUSS/TransformWrapper.h
		src/FOCUSS/GlobalVar.h
		src/FOCUSS/SomeFunctions.h
        src/FOCUSS/SomeFunctions.hxx
		src/GADGET_TOOLS/AcquisitionFinishGadget.h
		src/GADGET_TOOLS/CS_AccumulatorGadget.h
		src/GADGET_TOOLS/CS_MRIImageToAcquisitionGadget.h
		src/GADGET_TOOLS/SlicerGadget.h
		src/GADGET_TOOLS/CS_CombineGadget.h
		src/GADGET_TOOLS/hoNDFFT_CS.h
		src/GADGET_TOOLS/hoNDKLT_CS.h
        src/GADGET_TOOLS/hoNDDCT.h
        src/REGISTRATION/RegistrationGadget.h
		src/LAP_GADGET/LAPRegistrationGadget.h
		src/LAP_GADGET/lap3d.h
		src/LAP_GADGET/gaussianfilterbasis.h
		src/LAP_GADGET/mynamespaces.h
		src/LAP_GADGET/shiftengine3d.h
		src/GADGET_TOOLS/ImageAccumulatorGadget.h
		src/GADGET_TOOLS/ImageSlicerGadget.h
		src/RETRO/CS_Retro_AccumulatorGadget.h
		src/RETRO/CS_Retro_NavigatorGadget.h
		src/RETRO/CS_Retro_PopulationGadget.h
		include/old/hoNDArray_linalg.h
	)

	set(cs_lab_src_files CS_LAB.cpp
		src/FOCUSS/CS_Control.cpp
		src/FOCUSS/CS_FOCUSS.cpp
		src/FOCUSS/CS_FOCUSS_2D.cpp
		src/FOCUSS/CS_FOCUSS_2Dt.cpp
		src/FOCUSS/CS_FOCUSS_3D.cpp
		src/FOCUSS/CS_FOCUSS_4D.cpp
		src/FOCUSS/Transform.cpp
		src/FOCUSS/GlobalVar.cpp
		src/GADGET_TOOLS/AcquisitionFinishGadget.cpp
		src/GADGET_TOOLS/CS_AccumulatorGadget.cpp
		src/GADGET_TOOLS/CS_MRIImageToAcquisitionGadget.cpp
		src/GADGET_TOOLS/SlicerGadget.cpp
		src/GADGET_TOOLS/CS_CombineGadget.cpp
		src/GADGET_TOOLS/hoNDFFT_CS.cpp
		src/GADGET_TOOLS/hoNDKLT_CS.cpp
		src/GADGET_TOOLS/hoNDDCT.cpp
		src/REGISTRATION/RegistrationGadget.cpp
		src/LAP_GADGET/LAPRegistrationGadget.cpp
		src/LAP_GADGET/lap3d.cpp
		src/LAP_GADGET/gaussianfilterbasis.cpp
		src/LAP_GADGET/mynamespaces.cpp
		src/LAP_GADGET/shiftengine3d.cpp
		src/GADGET_TOOLS/ImageAccumulatorGadget.cpp
		src/GADGET_TOOLS/ImageSlicerGadget.cpp
		src/RETRO/CS_Retro_AccumulatorGadget.cpp
		src/RETRO/CS_Retro_NavigatorGadget.cpp
		src/RETRO/CS_Retro_PopulationGadget.cpp
		include/old/hoNDArray_linalg.cpp
	)
ENDIF()

set(cs_lab_config_files config/CS_LAB_offline.xml config/CS_LAB_online.xml config/CS_Retro_offline.xml config/CS_Retro_online_1.xml config/CS_Retro_online_2.xml)

IF((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	add_library(CS_LAB SHARED
		CS_LAB_export.h
		${cs_lab_header_files}
		${cs_lab_src_files}
		${cs_lab_config_files}
	)
ELSE()
	IF (WIN32)
			add_library(CS_LAB SHARED
				CS_LAB_export.h
				${cs_lab_header_files}
				${cs_lab_src_files}
				${cs_lab_config_files}
				${ISMRMRD_XSD_SOURCE}
			)
	ENDIF()
ENDIF()

IF (WIN32)
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
ELSE()
	set(ADDITIONAL_LIBRARIES
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpucore.so
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpufft.so
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpucore_math.so
	)
ENDIF()

IF ((${__GADGETRON__} GREATER 2) AND (${__GADGETRON_MINOR__} GREATER 5))
	target_link_libraries(	CS_LAB
							${ISMRMRD_LIBRARIES}
							${FFTW3_LIBRARIES}
							${BOOST_LIBRARIES}
							${ACE_LIBRARY}
							${ADDITIONAL_LIBRARIES}
							#${ITK_LIBRARIES}
							#${ELASTIX_LIBRARY}
							#${ELASTIX_ALL_COMPONENT_LIBS}
	)
ELSE ()
	target_link_libraries(	CS_LAB
							${ISMRMRD_LIBRARIES}
							${FFTW3_LIBRARIES}
							${BOOST_LIBRARIES}
							${XERCESC_LIBRARIES}
							${ACE_LIBRARY}
							${ADDITIONAL_LIBRARIES}
							${ITK_LIBRARIES}
							${ELASTIX_LIBRARY}
							${ELASTIX_ALL_COMPONENT_LIBS}
	)
ENDIF()
#set_property(TARGET CS_LAB PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install (FILES CS_LAB_export.h ${cs_lab_header_files} DESTINATION ${GADGETRON_HOME}/include)
install(TARGETS CS_LAB DESTINATION ${GADGETRON_HOME}/lib)

IF (WIN32)
	install(FILES ${cs_lab_config_files} DESTINATION ${GADGETRON_HOME}/config)
ELSE()
	install(FILES ${cs_lab_config_files} DESTINATION ${GADGETRON_HOME}/share/gadgetron/config)
ENDIF()

link_directories(${GADGETRON_LIB_DIR})


