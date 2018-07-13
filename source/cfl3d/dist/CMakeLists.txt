set(DIST_SRC
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
main.F ../libs/ccomplex.F)
# serial version

add_exe(EXEC Scfl3d_seq SRC ${DIST_SRC} CFLAG "${SINGLE_TYPE_FLAG}")
add_exe(EXEC Dcfl3d_seq SRC ${DIST_SRC} CFLAG "${DOUBLE_TYPE_FLAG}")
add_exe(EXEC Ccfl3dcmplx_seq SRC ${DIST_SRC} CFLAG "-DCMPLX ${SINGLE_TYPE_FLAG}")
add_exe(EXEC Zcfl3dcmplx_seq SRC ${DIST_SRC} CFLAG "-DCMPLX ${DOUBLE_TYPE_FLAG}")
if(IPO_OPT)
  foreach(TRGT Scfl3d_seq Dcfl3d_seq Ccfl3dcmplx_seq Zcfl3dcmplx_seq)
  set_target_properties(${TRGT} PROPERTIES INTERPROCEDURAL_OPTIMIZATION True)
  endforeach(TRGT Scfl3d_seq Dcfl3d_seq Ccfl3dcmplx_seq Zcfl3dcmplx_seq)
endif(IPO_OPT)

# 并行版本
if(${USE_MPI})
  add_exe(EXEC      Scfl3d_mpi SRC ${DIST_SRC} ${MPI_Fortran_LIBRARIES} CFLAG "-DDIST_MPI ${SINGLE_TYPE_FLAG}")
  add_exe(EXEC      Dcfl3d_mpi SRC ${DIST_SRC} ${MPI_Fortran_LIBRARIES} CFLAG "-DDIST_MPI ${DOUBLE_TYPE_FLAG}")
  add_exe(EXEC Ccfl3dcmplx_mpi SRC ${DIST_SRC} ${MPI_Fortran_LIBRARIES} CFLAG "-DDIST_MPI -DCMPLX ${SINGLE_TYPE_FLAG}")
  add_exe(EXEC Zcfl3dcmplx_mpi SRC ${DIST_SRC} ${MPI_Fortran_LIBRARIES} CFLAG "-DDIST_MPI -DCMPLX ${DOUBLE_TYPE_FLAG}")
    if(IPO_OPT)
      foreach(TRGT Scfl3d_mpi Dcfl3d_mpi Ccfl3dcmplx_mpi Zcfl3dcmplx_mpi)
        set_target_properties(${TRGT} PROPERTIES INTERPROCEDURAL_OPTIMIZATION True)
      endforeach(TRGT Scfl3d_mpi Dcfl3d_mpi Ccfl3dcmplx_mpi Zcfl3dcmplx_mpi)
    endif(IPO_OPT)
endif(${USE_MPI})
