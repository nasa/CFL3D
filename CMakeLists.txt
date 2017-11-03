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
option(USE_FASTIO "use fast-io of mpi? It is effective when mpi is found. " OFF)
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
		##add_definitions(-DCGNS)
	else(CGNS_FOUND)
		message(WARNING "system CGNS not found for CFL3D")
		set(USE_CGNS OFF CACHE BOOL "use CGNS? (not available)" FORCE) # set the CGNS not available
	endif(CGNS_FOUND)
else(USE_CGNS)
	message("NO CGNS for CFL3D")
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

		## if using MinGW in Windows 
		if(MINGW)
			message(WARNING "fast-io is not available in MinGW.")
			set(USE_FASTIO OFF CACHE BOOL "fast-io is not available in MinGW, turn it off" FORCE) # set the fastio not available

			message("Using MinGW compilers, in Windows, MPICH2 1.4.1p1 is supported.")
			get_filename_component(MPI_Fortran_LIB_DIR ${first_elem} DIRECTORY)
			find_library(MPICH2_WIN libfmpich2g.a ${MPI_Fortran_LIB_DIR})
			if(MPICH2_WIN)
				message("libfmpich2g.a of mpich2 for MinGW fortran 77 is found and used: ${MPICH2_WIN}")
				set_property(TARGET mpi PROPERTY IMPORTED_LOCATION "${MPICH2_WIN}")
				set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} "-localroot")
			else(MPICH2_WIN)
				message("libfmpich2g.a of mpich2 for MinGW fortran 77 is not found, trying to use ${first_elem}")
			endif(MPICH2_WIN)
		else(MINGW)
			## FAST-io
			if(USE_FASTIO)
				set_property(TARGET mpi APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "FASTIO")
			endif(USE_FASTIO)
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
if(UNIX)
	message("in Linux or Linux-like system")
	add_definitions(-DLINUX)
endif(UNIX)

if(MINGW)
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
else(MINGW)
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z muldefs")
endif(MINGW)

## use options to setup flags and find libraries
if(USE_SINGLE)
	# do nothing
else(USE_SINGLE)
	add_definitions(-DDBLE_PRECSN)
	if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
		set(CFL3D_FLAGS ${CFL3D_FLAGS} -fdefault-real-8)
	else(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
		set(CFL3D_FLAGS ${CFL3D_FLAGS} -r8)
	endif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
endif(USE_SINGLE)


## project compiler dependent settings
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	add_definitions(-DINTEL)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -xHost)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -traceback)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -fpe0)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -ip)	
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -fno-alias)
	# no `-w` option, do not turn warnings off.
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -march=native)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -Wno-align-commons)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -fbacktrace)
else()
	message(WARNING "Untested compiler iD: " ${CMAKE_Fortran_COMPILER_ID})
	message(WARNING "Supported compiler iD: Intel/GNU")
	message(WARNING "Try to use GNU options")
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -march=native)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -Wno-align-commons)
	set(CFL3D_FLAGS ${CFL3D_FLAGS} -fbacktrace)
endif()

message("Final compile flags for fortran is: ${CFL3D_FLAGS}")


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

#test
enable_testing() 
add_test(
	NAME 2D_SEQ_FLATPLATE
	COMMAND sh -c "$<TARGET_FILE:splitter><split.inp_1blk && $<TARGET_FILE:cfl3d_seq><grdflat5.inp"
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/Flatplate
)

add_test(
	NAME 2D_RONNIE_BACKSTEP
	COMMAND sh -c "$<TARGET_FILE:splitter> < split.inp_1blk && $<TARGET_FILE:ronnie> < ron_step_grdgen.inp && $<TARGET_FILE:cfl3d_seq> < step_grdgen.inp"
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/Backstep
)
add_test(
	NAME 2D_MAGGIE_NACA4412
	COMMAND sh -c "$<TARGET_FILE:splitter> < split.inp_xmera && $<TARGET_FILE:maggie> < maggie.inp && $<TARGET_FILE:cfl3d_seq> < cfl3d.inp_xmera"
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/NACA_4412
)

if(BUILD_MPI)
	add_test(
		NAME 2.5D_MPI_AXIBUMP
		COMMAND sh -c "$<TARGET_FILE:splitter> < split.inp_3blk  && ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfl3d_mpi> ${MPIEXEC_POSTFLAGS} < bumpperiodic.inp_3blk"
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/Axibump
	)
endif(BUILD_MPI)

if(USE_CGNS)
	add_test(
		NAME 3D_SEQ_DELTA_PLOT3D_TO_CGNS
		COMMAND sh -c "$<TARGET_FILE:plot3dg_to_cgns> < plot3dg_to_cgns.inp"
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/Delta_cgns
	)
	add_test(
		NAME 3D_SEQ_DELTA_CGNS_RUN
		COMMAND sh -c "$<TARGET_FILE:cfl3d_seq> < delta_cgns.inp"
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testcase/Delta_cgns
	)
endif(USE_CGNS)
