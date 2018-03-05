
project(LAPCpp)
cmake_minimum_required(VERSION 2.8)
#if(UNIX)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
#endif()
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})



aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
INCLUDE_DIRECTORIES(/usr/local/include/ITK-4.10)


TARGET_LINK_LIBRARIES(LAPCpp lapack openblas armadillo ${ITK_LIBRARIES})

#set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP)
if (OPENMP_FOUND AND CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -Wall -std=gnu++0x")
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
elseif(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3 -Wall -std=gnu++0x")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()



