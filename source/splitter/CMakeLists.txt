set(SD ${PROJECT_SOURCE_DIR}/source/cfl3d/libs)
set(DD ${PROJECT_SOURCE_DIR}/source/cfl3d/dist)
set(SPL_SRC 
         cfl3d_to_tlns3d.F  mapblk.F      movseg.F     spltsg.F
         getdhdr.F          mapin.F       outgrd.F     termn8.F
         ipadd.F            mapout.F      renmbr.F     tlns3d_to_cfl3d.F
         ipfix.F            mapsplt.F     shortinp.F   splitter.F
         iprvsg.F           mkintr.F      sizer.F      ${DD}/umalloc.F
         mkseg.F            splt2.F       ingrd.F
         ${SD}/readkey.F    ${SD}/parser.F  main.F)
set(SPLCMX_SRC
         cfl3d_to_tlns3d.F  mapblk.F      movseg.F        spltsg.F
         getdhdr.F          mapin.F       outgrd_cmplx.F  termn8.F
         ipadd.F            mapout.F      renmbr.F        tlns3d_to_cfl3d.F
         ipfix.F            mapsplt.F     shortinp.F      splitter_cmplx.F
         iprvsg.F           mkintr.F      sizer.F         ${DD}/umalloc.F
         mkseg.F            splt2.F       ingrd_cmplx.F
         ${DD}/umalloc_c.F  ${SD}/readkey.F  ${SD}/parser.F main_cmplx.F)
# 串行版本
add_exe(EXEC Ssplitter      SRC ${SPL_SRC}    CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dsplitter      SRC ${SPL_SRC}    CFLAG "${DOUBLE_TYPE_FLAG}")
add_exe(EXEC Csplittercmplx SRC ${SPLCMX_SRC} CFLAG "-DCMPLX ${SINGLE_TYPE_FLAG}")
add_exe(EXEC Zsplittercmplx SRC ${SPLCMX_SRC} CFLAG "-DCMPLX ${DOUBLE_TYPE_FLAG}")
