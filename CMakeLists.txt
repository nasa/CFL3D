project(cfl3d_software Fortran)
cmake_minimum_required(VERSION 3.9)

option(USE_MPI "compile parallel version of cfl3d" "ON")
option(IPO_OPT "build cfl3d program with Interprocedural Optimization option" "ON")
option(PLOT3D_SINGLE "write plot 3d With SINGLE precision real type" "ON")

##########################处理IPO##########################################
if(IPO_OPT)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "set build type to release" FORCE)
  if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      string(REGEX REPLACE "^lib" "xilib" CMAKE_Fortran_CREATE_STATIC_LIBRARY "${CMAKE_Fortran_CREATE_STATIC_LIBRARY}")
    endif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  endif(CMAKE_SYSTEM_NAME MATCHES "Windows")

  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fast")
    set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -ipo")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto")
    #set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -flto")
  endif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
endif(IPO_OPT)

#########################平台相关宏定义###################################
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
  add_definitions(-DWINDOWS)
endif(CMAKE_SYSTEM_NAME MATCHES "Windows")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions(-DGNU)
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

########################静态链接Fortran库#############################
if(${STATIC_FORTRAN_LIB})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -libs:static")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static-libgfortran")
  endif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
endif(${STATIC_FORTRAN_LIB})

###################单精度和双精度参数###########################
set(SINGLE_TYPE_FLAG "")
set(DOUBLE_TYPE_FLAG "-DDBLE_PRECSN")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(SINGLE_TYPE_FLAG "${SINGLE_TYPE_FLAG} -real-size:32")
    set(DOUBLE_TYPE_FLAG "${DOUBLE_TYPE_FLAG} -real-size:64")
  else(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(SINGLE_TYPE_FLAG "${SINGLE_TYPE_FLAG} -real-size 32")
    set(DOUBLE_TYPE_FLAG "${DOUBLE_TYPE_FLAG} -real-size 64")
  endif(CMAKE_SYSTEM_NAME MATCHES "Windows")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(SINGLE_TYPE_FLAG "${SINGLE_TYPE_FLAG}")
  #set(SINGLE_TYPE_FLAG "${SINGLE_TYPE_FLAG} -fdefault-real-4")
  set(DOUBLE_TYPE_FLAG "${DOUBLE_TYPE_FLAG} -fdefault-real-8")
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

#####################   MPI  ###########################
if(${USE_MPI})
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  if("${MPI_Fortran_LIBRARIES}" MATCHES "^.*/msmpi.lib$") #MSPPI
  	string(REGEX REPLACE "/msmpi.lib$" "/msmpifec.lib" TMP1 ${MPI_Fortran_LIBRARIES})
  	string(REGEX REPLACE "/msmpi.lib$" "/msmpifmc.lib" TMP2 ${MPI_Fortran_LIBRARIES})
  	set(MPI_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES};${TMP1};${TMP2}")
  endif()
else(${USE_MPI})
  message(WARNING "USE_MPI not set. Parallel version CFL3D will NOT Compile")
endif(${USE_MPI})
#write plot3d
if(PLOT3D_SINGLE)
  add_definitions(-DP3D_SINGLE)
endif(PLOT3D_SINGLE)

# set CFL3D program, module and library name(seq|mpi,)
# add excuate target,echo excuate shoud start with S,D,C,Z
# represent single,double,single complex,double complex
function(add_exe)
  cmake_parse_arguments(MY "" "EXEC;CFLAG" "SRC" ${ARGN})
  add_executable(${MY_EXEC} ${MY_SRC})
  # get type from first character of execute name(judge type)
  string(SUBSTRING ${MY_EXEC} 0 1 TYPE)
  # set compile flag
  set_target_properties(${MY_EXEC} PROPERTIES COMPILE_FLAGS "${MY_CFLAG}")
  # set fortran modules path
  if((${TYPE} STREQUAL S) OR (${TYPE} STREQUAL D))
    set_target_properties(${MY_EXEC} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${TYPE}mod)
  endif()
  # link libcommon library
  target_link_libraries(${MY_EXEC} ${TYPE}libcommon)
  # link mpi library
  string(REGEX MATCH "_mpi$" MA ${MY_EXEC})
  if(MA)
    target_link_libraries(${MY_EXEC} ${MPI_Fortran_LIBRARIES})
  endif(MA)
  # set install information
  string(SUBSTRING ${MY_EXEC} 1 -1 NAME)
  if((${TYPE} STREQUAL S) OR (${TYPE} STREQUAL C))
    set(DIR sbin)
  else()
    set(DIR dbin)
  endif()
  install(TARGETS ${MY_EXEC} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${DIR} RENAME ${NAME})
endfunction(add_exe)
add_subdirectory(source)

