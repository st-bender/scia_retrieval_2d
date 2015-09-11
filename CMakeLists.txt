cmake_minimum_required(VERSION 2.6)

project(scia_retrieval CXX C)

file(GLOB scia_SRC *.cpp *.c)
list(REMOVE_ITEM scia_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Polarization_correction.cpp")

find_library(SUNPERF sunperf /opt/SUNWspro/lib/amd64)
if (SUNPERF)
	message("sunperf found at ${SUNPERF}.")
	set(linalg_LIBS openblas sunperf)
	set(sp_PATH /opt/SUNWspro/lib/amd64)
	find_path(loc_LIBS NAMES libNOEM.a libgzstream.a
		PATHS ${CMAKE_SOURCE_DIR}/lib_solaris
		/home/msk-data/Scia/software/lib_solaris
		DOC "local libs")
else ()
	message("sunperf not found, using lapack and blas.")
	set(linalg_LIBS lapack blas)
	find_path(loc_LIBS NAMES libNOEM.a libgzstream.a
		PATHS ${CMAKE_SOURCE_DIR}/lib_linux
		/home/msk-data/Scia/software/lib_linux
		DOC "local libs")
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

#link_directories(${CMAKE_SOURCE_DIR}/lib ${sp_PATH})
link_directories(/usr/lb64 /usr/local/lib/amd64 ${loc_LIBS} ${sp_PATH})

add_definitions(-O2 -m64 -march=native -g -pedantic -Wall -Wextra -Wshadow)
add_definitions(-std=gnu++0x)
add_definitions(-ffast-math -fno-builtin)
add_definitions(-ftree-vectorize -mveclibabi=svml)
#add_definitions(-DDEBUG_NO_MATRICES)
#add_definitions(-DDEBUG_RETRIEVAL_MATRICES)
#add_definitions(-fopenmp -D_GLIBCXX_PARALLEL)
#add_definitions(-Wconversion -Wformat)
#add_definitions(-Wfloat-equal)
add_definitions(`${CMAKE_CURRENT_SOURCE_DIR}/git_version_flags.sh`)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")
add_executable(scia_retrieval ${scia_SRC})

target_link_libraries(scia_retrieval ${linalg_LIBS} NOEM gfortran ${ncdf_LIB} gzstream z gomp)
