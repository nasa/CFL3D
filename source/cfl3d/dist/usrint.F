      subroutine usrint
c
c     $Id$
c
c***********************************************************************
c      Purpose: provide termination procedure in the event of a system
c      signal
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#endif
c
      character*80 errfile
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
      errfile = 'cfl3d.error'
      open(unit=99,file=errfile,form='formatted',status='unknown')
c
      ierrflg = -999
      write(99,99) ierrflg
   99 format(' error code:',/,i4)
c
      write(99,1)
    1 format(/,' abnormal termination due to receipt of',
     .         ' system signal',/,
     .         ' (kill, floating pt. exception,',
     .         ' segmentation fault, etc.)',/)
c
      call my_flush(3)
      call my_flush(4)
      call my_flush(11)
      call my_flush(12)
      call my_flush(13)
      call my_flush(14)
      call my_flush(15)
      call my_flush(17)
      call my_flush(20)
      call my_flush(23)
      call my_flush(24)
      call my_flush(25)
      call my_flush(99)
c
#if defined DIST_MPI
         write (6,2) myid
    2    format('node',i4,' is terminating the program ',
     .   'due to receipt of a system signal')
         call MPI_ABORT(MPI_COMM_WORLD, myid, mpierror)
         call MPI_Finalize (ierr)
#endif
c
      stop
      end
