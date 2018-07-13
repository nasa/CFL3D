set(SD ${PROJECT_SOURCE_DIR}/source/cfl3d/libs)
set(DD ${PROJECT_SOURCE_DIR}/source/cfl3d/dist)
set(PRECFL_SRC
      ${DD}/sizer.F   get_cmd_args.F    ${SD}/rpatch0.F
    ${SD}/getibk0.F   ${SD}/global.F    ${SD}/global2.F
   ${DD}/pointers.F    ${DD}/compg2n.F  ${DD}/umalloc.F   ${SD}/cntsurf.F
       ${SD}/lead.F  ${SD}/pre_patch.F ${SD}/pre_blockbc.F
  ${SD}/pre_embed.F ${SD}/pre_period.F   ${DD}/usrint.F    ${DD}/termn8.F
    ${SD}/getdhdr.F    ${SD}/global0.F  ${SD}/readkey.F    ${SD}/parser.F
   ${SD}/ccomplex.F  ${SD}/cgnstools.F   ${SD}/setseg.F            main.F)
# 串行版本
add_exe(EXEC Sprecfl3d SRC ${PRECFL_SRC} CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dprecfl3d SRC ${PRECFL_SRC} CFLAG "${DOUBLE_TYPE_FLAG}")
