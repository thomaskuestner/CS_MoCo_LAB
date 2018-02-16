cmake_minimum_required(VERSION 2.6)
project(CS_LAB)

###############################################################
# Choose some options
###############################################################
# build static library? (Off=dynamic library)
option(BUILD_TOOLBOX_STATIC "Build static library for toolboxes" Off)
# Use OpenMP?
option(USE_OPENMP "Use OpenMP" On)
# Supress compilation warnings?
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# build options for 64 bits system
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
	message("64bit system is found")
	set(HAS_64_BIT On CACHE BOOL "64bit build")
else()
	message("32bit system is found")
	set(HAS_64_BIT Off CACHE BOOL "64bit build")
endif()

# build options for OpenMP support
find_package(OpenMP)
if (OPENMP_FOUND)
	if (USE_OPENMP)
		message("OpenMP multithreading enabled")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		add_definitions(-DUSE_OMP)
	else()
		message("OpenMP multithreading is supported, but disabled")
	endif()
else()
	message("OpenMP multithreading not supported")
endif()

if (WIN32)
	add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS)
# 	ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-D_VARIADIC_MAX=10) #to fix compiler limitations in Visual Studio Express
	if (HAS_64_BIT)
		add_definitions(-DWIN64 -D_WIN64)
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
	set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL:NO")
	set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} /INCREMENTAL:NO")
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
	set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
	set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
	set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
	# The two flags below is to fix Windows problems in relation to multiple defined operators new/delete and some constructors that are defined in our headers
# 	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE")
# 	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE")
else()
	if (UNIX)
		if (NOT APPLE)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
			set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libgfortran")
		endif()
	endif()
endif()

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/CS_LAB)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# whether to suppress compilation warnings
if (BUILD_SUPPRESS_WARNINGS)
	if (WIN32)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
	endif()
endif()

# whether to compile toolboxes as static library
if (BUILD_TOOLBOX_STATIC)
	message("Building static toolbox libray")
	add_definitions(-DBUILD_TOOLBOX_STATIC)
	set(LIBTYPE STATIC)
	if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
	endif()
else()
	message("Building dynamic toolbox libray")
	set(LIBTYPE SHARED)
endif()

###############################################################
# BOOST
###############################################################
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# necessary for Windows and RHEL <=6 systems
set(Boost_NO_BOOST_CMAKE ON)

if (WIN32)
	add_definitions(-DBOOST_ALL_NO_LIB)
	add_definitions(-DBOOST_ALL_DYN_LINK)
endif()

# We actually only use system and thread explicitly, but they require linking in date_time and chrono
if (WIN32)
	find_package(Boost COMPONENTS system thread date_time chrono program_options filesystem REQUIRED)
else()
	find_package(Boost COMPONENTS system thread program_options filesystem REQUIRED)
endif()

###############################################################
# ACE
###############################################################
find_package(ACE REQUIRED)

###############################################################
# CUDA
###############################################################
find_package(CUDA_advanced)
if (CUDA_FOUND)
	add_definitions(-DUSE_CUDA)
	set(GADGETRON_CUDA_FOUND_BOOL 1)
	include_directories(${CUDA_INCLUDE_DIRS})
	#set(CUDA_VERBOSE_BUILD ON)

	# Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
	# Support compute model 3.0 from Cuda 4.2 and up
	# Support compute model 3.5 from Cuda 5 and up
	message("Compiling with ${CUDA_NVCC_FLAGS}")
else()
	message("CUDA not found. CUDA components will not be compiled.")
	set(GADGETRON_CUDA_FOUND_BOOL 0)
endif()

###############################################################
# Armadillo
###############################################################
find_package(Armadillo)
# check whether ILP64 MKL should is used
if (ARMADILLO_FOUND)
	add_definitions(-DUSE_ARMADILLO)
	set(ARMADILLO_BLAS_LONG_LONG FALSE)
	if (EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp")
		# Read and parse armadillo config.hpp to find out whether BLAS uses long long
		file(STRINGS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp" _armadillo_blas_long_long REGEX "// #define ARMA_BLAS_LONG_LONG")
		if (NOT _armadillo_blas_long_long)
			set(ARMADILLO_BLAS_LONG_LONG TRUE)
			message("Armadillo is found to use long long for BLAS calls")
		else()
			message("Armadillo found. Note that ARMADILLO_BLAS_LONG_LONG _must_ be defined in ${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp to link against the MKL ILP64 interface.")
		endif()
		unset(_armadillo_blas_long_long)
	endif()
else()
	message("Armadillo not found. This will disable many toolboxes and gadgets.")
endif()

###############################################################
# HDF
###############################################################
find_package(HDF5 1.8 COMPONENTS C CXX HL)

###############################################################
# FFTW
###############################################################
find_package(FFTW3 COMPONENTS single double REQUIRED)

###############################################################
# ISMRMRD
###############################################################
find_package(Ismrmrd REQUIRED)

###############################################################
# GADGETRON
###############################################################
find_package(Gadgetron REQUIRED)
if (Gadgetron_FOUND)
	if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
		add_definitions(-D__GADGETRON_VERSION_HIGHER_3_6__)
	endif()
endif()

###############################################################
# XSD/XersesC
###############################################################
if (NOT ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5))))
	if (WIN32)
		find_package(XSD)
		find_package(XercesC)
	endif()
endif()

###############################################################
# MATLAB
###############################################################
find_package(GMatlab)
if (MATLAB_FOUND)
	add_definitions(-DUSE_MATLAB)
endif()

if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
	if (NOT (ACE_FOUND AND ISMRMRD_FOUND AND FFTW3_FOUND AND HDF5_FOUND))
		message("Required dependencies for gadget compilation not found (ACE, ISMRMRD, FFTW3, HDF5).")
	endif()
else()
	if (NOT (ACE_FOUND AND ISMRMRD_FOUND AND XSD_FOUND AND XERCESC_FOUND AND FFTW3_FOUND AND HDF5_FOUND))
		message("Required dependencies for gadget compilation not found (ACE, ISMRMRD, XSD, XERCESC, FFTW3, HDF5).")
	endif()
endif()

add_subdirectory(cmake)

if (WIN32)
	add_definitions(-D__BUILD_CS_LAB__)
endif()

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

###############################################################
# Elastix
###############################################################
# include elastix files and library
find_package(Elastix REQUIRED)

if (EXISTS ${ELASTIX_USE_FILE})
	include(${ELASTIX_USE_FILE})

	# These libraries are necessary for ElastixRegistrationGadget to find the elastixlib.h!
	# Please make sure to have elastix compiled as shared library if linking fails.
	# See CMakeLists.txt in elastix or its documentation for more information.
	# If these lines are commented out the resulting libCS_LAB.so could not be found by Gadgetron!
	link_libraries(param)
	link_libraries(elastix)
	link_libraries(transformix)
endif()

# include directories
include_directories(
	${CMAKE_SOURCE_DIR}/include
	${HDF5_INCLUDE_DIR}
	${HDF5_INCLUDE_DIR}/cpp
	${ARMADILLO_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/src
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

# include some specific directories
if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
	include_directories(
		${CMAKE_SOURCE_DIR}/include/new
	)
else()
	include_directories(
		${CMAKE_SOURCE_DIR}/include/old
		${XSD_INCLUDE_DIR}
		${XERCESC_INCLUDE_DIR}
	)
endif()

# set header files
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
	src/ELASTIX_GADGET/ElastixRegistrationGadget.h
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
)

# set source files
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
	src/ELASTIX_GADGET/ElastixRegistrationGadget.cpp
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
)

# set some specific header and source files
if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
	set(cs_lab_header_files ${cs_lab_header_files}
		include/new/hoNDArray_math_util.h
		include/new/GadgetronCommon.h
		include/new/hoNDArray_operators.h
		include/new/hoNDArray_blas.h
		include/new/hoMatrix_util.h
	)

	set(cs_lab_src_files ${cs_lab_src_files}
		include/new/hoNDArray_math_util.cpp
		include/new/hoNDArray_operators.cpp
	)
else()
	set(cs_lab_header_files ${cs_lab_header_files}
		include/old/hoNDArray_linalg.h
	)

	set(cs_lab_src_files ${cs_lab_src_files}
		include/old/hoNDArray_linalg.cpp
	)
endif()

set(cs_lab_config_files config/CS_LAB_offline.xml config/CS_LAB_online.xml config/CS_Retro_offline.xml config/CS_Retro_online_1.xml config/CS_Retro_online_2.xml)

if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
	add_library(CS_LAB SHARED
		CS_LAB_export.h
		${cs_lab_header_files}
		${cs_lab_src_files}
		${cs_lab_config_files}
	)
else()
	add_library(CS_LAB SHARED
		CS_LAB_export.h
		${cs_lab_header_files}
		${cs_lab_src_files}
		${cs_lab_config_files}
		${ISMRMRD_XSD_SOURCE}
	)
endif()

if (WIN32)
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
	set(ADDITIONAL_LIBRARIES
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpucore.so
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpufft.so
		${GADGETRON_HOME}/lib/libgadgetron_toolbox_cpucore_math.so
	)
endif()

if ((${GADGETRON_VERSION_MAJOR} GREATER 3) OR ((${GADGETRON_VERSION_MAJOR} EQUAL 3) AND (${GADGETRON_VERSION_MINOR} GREATER 5)))
	target_link_libraries( CS_LAB
		${ISMRMRD_LIBRARIES}
		${FFTW3_LIBRARIES}
		${BOOST_LIBRARIES}
		${ACE_LIBRARY}
		${ADDITIONAL_LIBRARIES}
		#${ITK_LIBRARIES}
		#${ELASTIX_LIBRARY}
		#${ELASTIX_ALL_COMPONENT_LIBS}
	)
else ()
	target_link_libraries(CS_LAB
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
endif()

#set_property(TARGET CS_LAB PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install (FILES CS_LAB_export.h ${cs_lab_header_files} DESTINATION ${GADGETRON_HOME}/include)
install(TARGETS CS_LAB DESTINATION ${GADGETRON_HOME}/lib)

if (WIN32)
	install(FILES ${cs_lab_config_files} DESTINATION ${GADGETRON_HOME}/config)
else()
	install(FILES ${cs_lab_config_files} DESTINATION ${GADGETRON_HOME}/share/gadgetron/config)
endif()

link_directories(${GADGETRON_LIB_DIR})


