cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(
DIST_FSRCS
bc_blkint.F    findmin_new.F  plot3d.F       rrest.F 
bc_patch.F     forceout.F     plot3t.F       rrestg.F 
calyplus.F     pointers.F     setup.F        writ_buf.F 
mgblk.F        qinter.F       prntcp.F       newalpha.F 
cputim.F       patcher.F      qout.F         termn8.F 
dynptch.F      plot3c.F       resp.F         usrint.F 
wrest.F        wrestg.F       pre_bc.F       bc_embed.F 
updateg.F      compg2n.F      resetg.F       bc_period.F 
yplusout.F     sizer.F        cfl3d.F        trnsfr_vals.F 
updatedg.F     ae_corr.F      mgbl.F         setslave.F 
umalloc.F      reass.F        qoutavg.F      plot3davg.F 
qout_coarse.F  qout_2d.F      plot3d_2d.F    plot3d_coarse.F
)

add_library(dist ${DIST_FSRCS})
target_link_libraries(dist PUBLIC common)
if(USE_CGNS)
	target_link_libraries(dist PUBLIC cgns)
endif(USE_CGNS)

set(
CFL3D_MAIN_FSRCS
../libs/ccomplex.F development.F main.F
)
add_executable(cfl3d_seq ${CFL3D_MAIN_FSRCS})
target_link_libraries(cfl3d_seq PUBLIC dist)

if(BUILD_MPI)
	add_library(dist_mpi ${DIST_FSRCS})
	target_link_libraries(dist_mpi PUBLIC common)
	if(USE_CGNS)
		target_link_libraries(dist_mpi PUBLIC cgns)
	endif(USE_CGNS)
	target_link_libraries(dist_mpi PUBLIC mpi)
	add_executable(cfl3d_mpi ${CFL3D_MAIN_FSRCS})
	target_link_libraries(cfl3d_mpi PUBLIC dist_mpi)
endif(BUILD_MPI)

if(BUILD_CMPLX)
	add_library(dist_cmplx ${DIST_FSRCS})
	target_link_libraries(dist_cmplx PUBLIC common_cmplx)
	if(USE_CGNS)
		target_link_libraries(dist_cmplx PUBLIC cgns)
	endif(USE_CGNS)
	add_executable(cfl3dcmplx_seq ${CFL3D_MAIN_FSRCS})
	target_link_libraries(cfl3dcmplx_seq PUBLIC dist_cmplx)

	if(BUILD_MPI)
		add_library(dist_cmplx_mpi ${DIST_FSRCS})
		target_link_libraries(dist_cmplx_mpi PUBLIC common_cmplx)
		if(USE_CGNS)
			target_link_libraries(dist_cmplx_mpi PUBLIC cgns)
		endif(USE_CGNS)
		target_link_libraries(dist_cmplx_mpi PUBLIC mpi)
		add_executable(cfl3dcmplx_mpi ${CFL3D_MAIN_FSRCS})
		target_link_libraries(cfl3dcmplx_mpi PUBLIC dist_cmplx_mpi)
	endif(BUILD_MPI)
endif(BUILD_CMPLX)
