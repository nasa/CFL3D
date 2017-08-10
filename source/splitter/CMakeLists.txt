cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(
 SPLIT_FSRCS
 cfl3d_to_tlns3d.F  mapblk.F      movseg.F     spltsg.F 
 getdhdr.F          mapin.F       outgrd.F     termn8.F 
 ipadd.F            mapout.F      renmbr.F     tlns3d_to_cfl3d.F 
 ipfix.F            mapsplt.F     shortinp.F   splitter.F 
 iprvsg.F           mkintr.F      sizer.F       
 mkseg.F            splt2.F       ingrd.F
 main.F
)

add_executable(splitter ${SPLIT_FSRCS})
target_link_libraries(splitter PUBLIC common)
target_link_libraries(splitter PUBLIC dist)
if(BUILD_CMPLX)
	add_executable(splittercmplx ${SPLIT_FSRCS})
	target_link_libraries(splittercmplx PUBLIC common_cmplx)
	target_link_libraries(splittercmplx PUBLIC dist_cmplx)
endif(BUILD_CMPLX)

