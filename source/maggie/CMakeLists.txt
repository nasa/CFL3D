set(SD ${PROJECT_SOURCE_DIR}/source/cfl3d/libs)
set(DD ${PROJECT_SOURCE_DIR}/source/cfl3d/dist)
set(MAG_SRC maggie.F cputim.F wkstn.F ${SD}/my_flush.F)
# 串行版本
add_exe(EXEC Smaggine SRC ${MAG_SRC} CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dmaggine SRC ${MAG_SRC} CFLAG "${DOUBLE_TYPE_FLAG}")
set_target_properties(Smaggine PROPERTIES INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/header")
set_target_properties(Dmaggine PROPERTIES INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/header")