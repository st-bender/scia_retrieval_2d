cmake_minimum_required(VERSION 2.8)

project(scia_retrieval CXX C)

file(GLOB scia_SRC *.cpp *.c)

find_library(SUNPERF sunperf /opt/SUNWspro/lib/amd64)
if (SUNPERF)
	message("sunperf found at ${SUNPERF}.")
	set(linalg_LIBS openblas sunperf)
	set(sp_PATH /opt/SUNWspro/lib/amd64)
else ()
	message("sunperf not found, using lapack and blas.")
	set(linalg_LIBS lapack blas)
endif ()

find_library(NOEM NOEM
	${CMAKE_SOURCE_DIR}/lib_solaris
	${CMAKE_SOURCE_DIR}/lib_linux
	/home/msk-data/Scia/software/lib_linux)
if (NOEM)
	message("NOEM model found at ${NOEM}.")
	get_filename_component(noem_DIR ${NOEM} PATH)
	add_definitions(-DHAVE_NOEM)
	set(noem_LIB NOEM)
else ()
	message("NOEM model not found.")
endif ()

find_library(NCDFF netcdff)
if (NCDFF)
	set(ncdf_LIB netcdff)
else ()
	set(ncdf_LIB netcdf)
endif ()

#find_path(loc_LIBS NAMES libNOEM.a libgzstream.a
#	PATHS
#	${CMAKE_SOURCE_DIR}/lib_solaris
#	${CMAKE_SOURCE_DIR}/lib_linux_clang
#	/home/msk-data/Scia/software/lib_solaris
#	/home/msk-data/Scia/software/lib_linux
#	DOC "local libs")

link_directories(/usr/lib64
	/usr/local/lib/amd64
	/home/Processor/Work_SB/openblas_64.i/lib
	${noem_DIR} ${sp_PATH})

include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES} ${EIGEN3_INCLUDES}")
check_include_file_cxx("Eigen/Dense" HAVE_EIGEN3)
if (HAVE_EIGEN3)
	include_directories(${EIGEN3_INCLUDES})
	add_definitions(-DHAVE_EIGEN3)
endif ()

include(FindHDF5)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})

include_directories(.)
add_definitions(-O2 -march=native -g -pedantic -Wall -Wextra -Wshadow)
add_definitions(-std=gnu++0x)
add_definitions(-ffast-math -fno-builtin)
add_definitions(-ftree-vectorize -mveclibabi=svml)
#add_definitions(-DDEBUG_NO_MATRICES)
#add_definitions(-DDEBUG_RETRIEVAL_MATRICES)
add_definitions(-fopenmp -D_GLIBCXX_PARALLEL)
#add_definitions(-Wconversion -Wformat)
#add_definitions(-Wfloat-equal)
add_definitions(`${CMAKE_CURRENT_SOURCE_DIR}/git_version_flags.sh`)
add_executable(scia_retrieval ${scia_SRC})

target_link_libraries(scia_retrieval ${linalg_LIBS} ${noem_LIB} gfortran netcdf ${ncdf_LIB} z gomp)
