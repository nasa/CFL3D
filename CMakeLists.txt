cmake_minimum_required(VERSION 3.0.2)
cmake_policy(SET CMP0048 NEW) #used to set Version Number
project(CFL3D VERSION 6.7)
enable_language(Fortran C)
SET (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include CACHE PATH "Single Directory for all fortran modules.")

## options
option(USE_SINGLE "use single precision? " OFF)
option(USE_CGNS "use CGNS? " OFF)
option(USE_NO_REDIRECT "use no-redirect? " OFF)
option(USE_FASTIO "use fast-io of mpi? It is effective when mpi is found. " ON)
option(BUILD_MPI "build MPI version? " OFF)
option(BUILD_CMPLX "build complex version? " OFF)

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
# use Paraview's FindCGNS.cmake
# https://github.com/Kitware/ParaView/raw/master/CMake/FindCGNS.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
if(USE_CGNS)
	find_package(CGNS 2.4 REQUIRED)
	if(CGNS_FOUND)
		message("Found CGNS")
		add_definitions(-DCGNS)
	else(CGNS_FOUND)
		message(WARNING "CGNS not found, try to use https://github.com/Kitware/ParaView/raw/master/CMake/FindCGNS.cmake and place it in /cmake/Modules/")
		set(USE_CGNS OFF CACHE BOOL "use CGNS? (not available)" FORCE) # set the CGNS not available
	endif(CGNS_FOUND)
else(USE_CGNS)
endif(USE_CGNS)

if(USE_NO_REDIRECT)
	add_definitions(-DNOREDIRECT)
endif(USE_NO_REDIRECT)

if(BUILD_MPI)
	find_package(MPI 1.0 REQUIRED)
	if(MPI_FOUND)
		message("Found MPI")
		## FAST-io
		if(USE_FASTIO)
			add_definitions(-DFASTIO)
		endif(USE_FASTIO)
	else(MPI_FOUND)
		## FAST-io
		if(USE_FASTIO)
			message(WARNING "no mpi, no fast-io.")
			set(USE_FASTIO OFF CACHE BOOL "build MPI? (not available)" FORCE) # set the fastio not available
		endif(USE_FASTIO)

		set(BUILD_MPI OFF CACHE BOOL "build MPI? (not available)" FORCE) # set the MPI not available
	endif(MPI_FOUND)
else(BUILD_MPI)
endif(BUILD_MPI)

# always use LINUX
add_definitions(-DLINUX)

## compiler dependent settings
if(${CMake_Fortran_COMPILER_ID} MATCHES Intel)
	add_definitions(-DINTEL)
	add_compile_options(-xHost )
	add_compile_options(-traceback)
	add_compile_options(-fpe0)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
	add_compile_options(-march=native -Wno-align-commons)
	add_compile_options(-fbacktrace)
else(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
	message(WARNING "Unsupported compiler iD: " ${CMAKE_Fortran_COMPILER_ID})
	message(WARNING "Supported compiler iD: Intel/GNU")
	message(WARNING "Try to use GNU options")
	add_compile_options(-march=native -Wno-align-commons)
	add_compile_options(-fbacktrace)
endif(${CMake_Fortran_COMPILER_ID} MATCHES Intel)

add_definitions(-DP3D_SINGLE)

message("---------------------CFL3D programs summary--------------------------")
message("Build type? " ${CMAKE_BUILD_TYPE})
message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message("CMAKE Fortran compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
message("---------------------------------------------------------------------")

# set output library
# use CMAKE_XXX_OUTPUT_DIRECTORY to set the default dirs when create targets
# set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib/fortran)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(source) 

