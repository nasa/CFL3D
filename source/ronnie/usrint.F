      subroutine usrint
c
c     $Id$
c
c***********************************************************************
c      Purpose: provide termination procedure in the event of a system
c      signal
c***********************************************************************
#if defined DIST_MPI
c
      include 'mpif.h'
#endif
c
      character*80 errfile
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
      errfile = 'ronnie.error'
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
#   ifdef IBM
      call flush_(66)
      call flush_(25)
      call flush_(99)
#   else
      call flush(66)
      call flush(25)
      call flush(99)
#   endif
c
      close(66)
      close(25)
      close(99)
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
