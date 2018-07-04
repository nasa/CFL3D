set(SD ${PROJECT_SOURCE_DIR}/source/cfl3d/libs)
set(DD ${PROJECT_SOURCE_DIR}/source/cfl3d/dist)
set(TOOLSLIB_SRC
        ${DD}/umalloc_c.F ${DD}/umalloc_r.F ${SD}/readkey.F    ${SD}/global0.F
        ${SD}/cgnstools.F ${DD}/termn8.F    ${SD}/parser.F     ${SD}/my_flush.F
        ${SD}/outbuf.F)

add_library(Slibtools STATIC ${TOOLSLIB_SRC})
set_target_properties(Slibtools PROPERTIES COMPILE_FLAGS "${SINGLE_TYPE_FLAG}")

add_library(Dlibtools STATIC ${TOOLSLIB_SRC})
set_target_properties(Dlibtools PROPERTIES COMPILE_FLAGS "${DOULBE_TYPE_FLAG}")

add_library(Clibtools STATIC ${TOOLSLIB_SRC})
set_target_properties(Clibtools PROPERTIES COMPILE_FLAGS "-DCMPLX ${SINGLE_TYPE_FLAG}")

add_library(Zlibtools STATIC ${TOOLSLIB_SRC})
set_target_properties(Zlibtools PROPERTIES COMPILE_FLAGS "-DCMPLX ${DOULBE_TYPE_FLAG}")
set(TOOL_EXEC grid_perturb       Get_FD            moovmaker
              plot3dg_to_cgns    XINTOUT_to_ovrlp  p3d_to_INGRID     INGRID_to_p3d
              cfl3d_to_pegbc     everyother_xyz    initialize_field  p3d_to_cfl3drst
              v6_restart_mod     v6inpdoubhalf     cgns_to_cfl3dinput)

add_custom_target(tools)
# set CFL3D program, module and library name(seq|mpi,)
# add excuate target,echo excuate shoud start with S,D,C,Z
# represent single,double,single complex,double complex
function(add_tools)
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
  # link libtools library
  target_link_libraries(${MY_EXEC} ${TYPE}libtools)
  # set install information
  string(SUBSTRING ${MY_EXEC} 1 -1 NAME)
  if((${TYPE} STREQUAL S) OR (${TYPE} STREQUAL C))
    set(DIR sbin)
  else()
    set(DIR dbin)
  endif()
  install(TARGETS ${MY_EXEC} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${DIR} RENAME ${NAME})
endfunction(add_tools)

# 串行版本
foreach(TE ${TOOL_EXEC})
  add_tools(EXEC "S${TE}" SRC ${TE}.F CFLAG "${SINGLE_TYPE_FLAG}")
  add_tools(EXEC "D${TE}" SRC ${TE}.F CFLAG "${DOUBLE_TYPE_FLAG}")
  add_dependencies(tools "S${TE}" "D${TE}")
endforeach(TE ${TOOL_EXEC})
add_tools(EXEC "Cgrid_perturb_cmplx" SRC grid_perturb_cmplx.F CFLAG "-DCMPLX ${SINGLE_TYPE_FLAG}")
add_tools(EXEC "Zgrid_perturb_cmplx" SRC grid_perturb_cmplx.F CFLAG "-DCMPLX ${DOUBLE_TYPE_FLAG}")
add_dependencies(tools Cgrid_perturb_cmplx Zgrid_perturb_cmplx)
