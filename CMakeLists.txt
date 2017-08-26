cmake_minimum_required(VERSION 2.6)
project(CFL3D)
set(CFL3D_VERSION_MAJOR 6)
set(CFL3D_VERSION_MINOR 7)
enable_language(Fortran C)
SET (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include CACHE PATH "Single Directory for all fortran modules.")

## output some information
message("CMAKE_SYSTEM: ${CMAKE_SYSTEM}")
message("MSYS: ${MSYS}")
## options
option(USE_SINGLE "use single precision? " OFF)
option(USE_CGNS "use CGNS? " ON)
option(USE_NO_REDIRECT "use no-redirect? " OFF)
option(USE_FASTIO "use fast-io of mpi? It is effective when mpi is found. " ON)
option(BUILD_MPI "build MPI version? " ON)
option(BUILD_CMPLX "build complex version? " ON)

## default Build Type is Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

## no redirect
if(USE_NO_REDIRECT)
	add_definitions(-DNOREDIRECT)
endif(USE_NO_REDIRECT)


## use options to setup flags and find libraries
if(USE_SINGLE)
	# do nothing
else(USE_SINGLE)
	add_definitions(-DDBLE_PRECSN)

	if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
		message(STATUS "Using the Intel Fortran compiler")
		add_compile_options(-r8)
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
		message(STATUS "Using the GNU Fortran compiler")
		add_compile_options(-fdefault-real-8)
	else(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
		message(WARNING "compiler do not match Intel/GNU! use Intel's setting instead.")
		add_compile_options(-r8)
	endif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

endif(USE_SINGLE)


# no FindCGNS in CMake release
# using Paraview's FindCGNS.cmake
# https://github.com/Kitware/ParaView/raw/master/CMake/FindCGNS.cmake
# modified to find `cgnslib_f.h` and `libcgns.a`
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(USE_CGNS)
	find_package(CGNS 2.4 REQUIRED)
	if(CGNS_FOUND)
		message("Found CGNS for CFL3D")
		message("CGNS include dir: ${CGNS_INCLUDE_DIR}")
		message("CGNS library: ${CGNS_LIBRARY}")
		add_library(cgns STATIC IMPORTED)
		set_property(TARGET cgns PROPERTY IMPORTED_LOCATION ${CGNS_LIBRARY})
		set_property(TARGET cgns PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES Fortran)
		set_property(TARGET cgns PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${CGNS_INCLUDE_DIR}>)
		set_property(TARGET cgns PROPERTY INTERFACE_COMPILE_DEFINITIONS CGNS)
		# If using libcgns.so in system lib, it is required to find and use libhdf5.so
		#find_package(HDF5 REQUIRED)
		#if(HDF5_FOUND)
		#	message("HDF5_lib:${HDF5_hdf5_LIBRARY_RELEASE}")
		#	set_property(TARGET cgns PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_hdf5_LIBRARY_RELEASE})
		#else(HDF5_FOUND)
		#	message(WARNING "HDF5 not found, there might be some error when linking to libcgns.so")
		#endif(HDF5_FOUND)
	else(CGNS_FOUND)
		message(WARNING "system CGNS not found for CFL3D")
		set(USE_CGNS OFF CACHE BOOL "use CGNS? (not available)" FORCE) # set the CGNS not available
	endif(CGNS_FOUND)
else(USE_CGNS)
endif(USE_CGNS)

if(BUILD_MPI)
	find_package(MPI 1.0 REQUIRED)
	if(MPI_Fortran_FOUND)
		message("Found MPI for Fortran")
		## add mpi library
		add_library(mpi STATIC IMPORTED)
		set_property(TARGET mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
		message("Choose the first element of libs: ${MPI_Fortran_LIBRARIES}")
		set(other_elems ${MPI_Fortran_LIBRARIES})
		list(GET other_elems 0 first_elem)
		list(REMOVE_AT other_elems 0)
		set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${first_elem})
		set_property(TARGET mpi PROPERTY INTERFACE_LINK_LIBRARIES ${other_elems})
		set_property(TARGET mpi APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "DIST_MPI")
		set_property(TARGET mpi APPEND PROPERTY LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
		set_property(TARGET mpi APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${MPI_Fortran_COMPILE_FLAGS}")
		## FAST-io
		if(USE_FASTIO)
			set_property(TARGET mpi APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "FASTIO")
		endif(USE_FASTIO)
		## if using MinGW in Windows 
		if(MINGW)
			message("Using MinGW compilers, in Windows, MPICH2 1.4.1p1 is supported.")
			get_filename_component(MPI_Fortran_LIB_DIR ${first_elem} DIRECTORY)
			find_library(MPICH2_WIN libfmpich2g.a ${MPI_Fortran_LIB_DIR})
			if(MPICH2_WIN)
			message("libfmpich2g.a of mpich2 for MinGW fortran 77 is found and used: ${MPICH2_WIN}")
			set_property(TARGET mpi PROPERTY IMPORTED_LOCATION "${MPICH2_WIN}")
			else(MPICH2_WIN)
			message("libfmpich2g.a of mpich2 for MinGW fortran 77 is not found, trying to use ${first_elem}")
			endif(MPICH2_WIN)
		endif(MINGW)

	else(MPI_Fortran_FOUND)
		## FAST-io
		if(USE_FASTIO)
			message(WARNING "no mpi, no fast-io.")
			set(USE_FASTIO OFF CACHE BOOL "build MPI? (not available)" FORCE) # set the fastio not available
		endif(USE_FASTIO)
		set(BUILD_MPI OFF CACHE BOOL "build MPI? (not available)" FORCE) # set the MPI not available
	endif(MPI_Fortran_FOUND)
else(BUILD_MPI)
endif(BUILD_MPI)

# in LINUX system
if(UNIX AND NOT APPLE)
	message("in Linux system")
	add_definitions(-DLINUX)
endif(UNIX AND NOT APPLE)


if(MINGW)
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
else(MINGW)
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z muldefs")
endif(MINGW)

## project compiler dependent settings
if(${CMake_Fortran_COMPILER_ID} MATCHES Intel)
	add_definitions(-DINTEL)
	add_compile_options(-xHost )
	add_compile_options(-traceback)
	add_compile_options(-fpe0) # 
	#add_compile_options(-w) #bad behavior, old codes usually use common block incorrectly. However, the warning is harmless. 
	add_compile_options(-ip) # interprocedural optimization
	add_compile_options(-fno-alias) #gcc by default uses `-fargument-noalias-global` 
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
	add_compile_options(-march=native)
	add_compile_options(-Wno-align-commons) #suppress some warnings about the padding of commons
	add_compile_options(-fbacktrace)
else(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
	message(WARNING "Unsupported compiler iD: " ${CMAKE_Fortran_COMPILER_ID})
	message(WARNING "Supported compiler iD: Intel/GNU")
	message(WARNING "Try to use GNU options")
	add_compile_options(-march=native)
	add_compile_options(-Wno-align-commons) #suppress some warnings about the padding of commons
	add_compile_options(-fbacktrace)
endif(${CMake_Fortran_COMPILER_ID} MATCHES Intel)

# it seems always use P3D_SINGLE in original Install script
add_definitions(-DP3D_SINGLE)

message("---------------------CFL3D programs summary--------------------------")
message("Build type? " ${CMAKE_BUILD_TYPE})
message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message("CMAKE Fortran compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
message("---------------------------------------------------------------------")

# set output library
# use CMAKE_XXX_OUTPUT_DIRECTORY to set the default dirs when create targets
# set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib/fortran)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(source) 
#add_subdirectory(testcases) 

