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
#ifdef CRAY_TIME
#else
#   ifdef IBM
      call flush_(3)
      call flush_(4)
      call flush_(11)
      call flush_(12)
      call flush_(13)
      call flush_(14)
      call flush_(15)
      call flush_(17)
      call flush_(20)
      call flush_(23)
      call flush_(24)
      call flush_(25)
      call flush_(99)
#   else
      call flush(3)
      call flush(4)
      call flush(11)
      call flush(12)
      call flush(13)
      call flush(14)
      call flush(15)
      call flush(17)
      call flush(20)
      call flush(23)
      call flush(24)
      call flush(25)
      call flush(99)
#   endif
#endif
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
