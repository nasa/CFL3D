set(SD ${PROJECT_SOURCE_DIR}/source/cfl3d/libs)
set(DD ${PROJECT_SOURCE_DIR}/source/cfl3d/dist)
set(RON_SRC ${SD}/arc.F    ${SD}/diagnos.F    ${SD}/extrae.F    ${SD}/newfit.F    ${SD}/rp3d.F
         ${SD}/avgint.F     ${SD}/direct.F   ${SD}/global2.F   ${DD}/patcher.F         setup.F
       ${SD}/collapse.F      ${SD}/dsmin.F    ${SD}/invert.F   ${SD}/project.F    ${SD}/shear.F
         ${SD}/expand.F       ${SD}/lead.F     ${SD}/rechk.F     ${SD}/topol.F         usrint.F
          ${SD}/collx.F      ${SD}/extra.F    ${SD}/loadgr.F          ronnie.F     ${SD}/topol2.F
          ${SD}/trace.F         ${SD}/xe.F       ${SD}/xe2.F    ${DD}/cputim.F     ${DD}/writ_buf.F
               termn8.F     ${SD}/transp.F    ${SD}/rotatp.F   ${DD}/umalloc.F
      ${SD}/pre_patch.F            sizer.F  ${SD}/ccomplex.F            main.F)
# 串行版本
add_exe(EXEC Sronnie SRC ${RON_SRC} CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dronnie SRC ${RON_SRC} CFLAG "${DOUBLE_TYPE_FLAG}")

set(PRERON_SRC        ${SD}/global2.F       termn8.F    usrint.F    ${DD}/umalloc.F
         sizer.F     ${SD}/ccomplex.F     main_pre.F)
add_exe(EXEC Spreronnie SRC ${PRERON_SRC} CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dpreronnie SRC ${PRERON_SRC} CFLAG "${DOUBLE_TYPE_FLAG}")
