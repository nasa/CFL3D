cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(
	PRECFL_FSRCS
	get_cmd_args.F
	main.F
	usrint.F
)

set(PRECFL_SHARED_FSRCS
	${CFLSRC_D}/cputim.F 
	${CFLSRC_D}/writ_buf.F 
	${CFLSRC_D}/umalloc.F 
	${CFLSRC_D}/patcher.F 
	${CFLSRC_S}/lead.F 
	${CFLSRC_S}/pre_patch.F 
	${CFLSRC_S}/collx.F 
	${CFLSRC_S}/global2.F 
	${CFLSRC_S}/arc.F 
	${CFLSRC_S}/avgint.F 
	${CFLSRC_S}/collapse.F 
	${CFLSRC_S}/diagnos.F 
	${CFLSRC_S}/direct.F 
	${CFLSRC_S}/dsmin.F 
	${CFLSRC_S}/expand.F 
	${CFLSRC_S}/extra.F 
	${CFLSRC_S}/extrae.F 
	${CFLSRC_S}/invert.F 
	${CFLSRC_S}/loadgr.F 
	${CFLSRC_S}/newfit.F 
	${CFLSRC_S}/project.F 
	${CFLSRC_S}/rechk.F 
	${CFLSRC_S}/rp3d.F 
	${CFLSRC_S}/shear.F 
	${CFLSRC_S}/topol.F 
	${CFLSRC_S}/topol2.F 
	${CFLSRC_S}/trace.F 
	${CFLSRC_S}/xe.F 
	${CFLSRC_S}/xe2.F 
	${CFLSRC_S}/outbuf.F 
	${CFLSRC_S}/transp.F 
	${CFLSRC_S}/rotatp.F 
	${CFLSRC_S}/ccomplex.F 
	${CFLSRC_S}/my_flush.F 
)

add_executable(precfl3d ${PRECFL_FSRCS})
target_link_libraries(precfl3d PUBLIC common)
target_link_libraries(precfl3d PUBLIC dist)


