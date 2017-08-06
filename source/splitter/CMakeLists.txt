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

#set (SPLIT_FSRCS ${SPLIT_FSRCS}
#../cfl3d/dist/umalloc.F 
#../cfl3d/libs/readkey.F          ../cfl3d/libs/parser.F      ../cfl3d/libs/outbuf.F     ../cfl3d/libs/my_flush.F 
#)

add_executable(splitter ${SPLIT_FSRCS})
target_link_libraries(splitter PUBLIC common)
target_link_libraries(splitter PUBLIC dist)

##if(BUILD_MPI)
##	add_executable(cfl3d_mpi ${CFL3D_MAIN_FSRCS})
##	target_link_libraries(cfl3d_mpi PUBLIC dist)
##	target_link_libraries(cfl3d_mpi PUBLIC common)
##endif(BUILD_MPI)
