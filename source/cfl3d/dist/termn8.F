      subroutine termn8(myid_stop,ierrflg,ibufdim,nbuf,bou,nou)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Flush output buffers (distributed version) so that last
c      information written to output files will be printed in case of
c      a program termination; also print error flag to signal
c      status at program termination, and close all files (except those
c      which are closed elsewhere).
c
c      error codes (ierrflg):
c
c       = 0  normal termination
c       < 0  abnormal termination
c
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
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
      write(99,99) ierrflg
   99 format(' error code:',/,i4)
c
      if (ierrflg.eq.0) then
            write(99,1)
    1       format(/,' execution terminated normally')
      else 
         if (ierrflg.ne.-999) then
            if (ierrflg.eq.-99) then
               if (nou(1).gt.0) then
                  write(99,22)
   22             format(/,' abnormal termination during array',
     .            ' sizing (precfl3d)',/,
     .            ' (error message follows)',/)
               else
                  write(99,33)
   33             format(/,' abnormal termination during array',
     .            ' sizing (precfl3d)',/,
     .            ' (see precfl3d.out for error message)',/)
               end if
            else
               iflag = 0
               do nn=1,nbuf
                  if (nou(nn).gt.0) iflag = 1
               end do
               if (iflag.gt.0) then
                  write(99,2)
    2             format(/,' abnormal termination due to cfl3d',
     .            ' error check',/,
     .            ' (error message follows)',/)
               else
                  write(99,3)
    3             format(/,' abnormal termination due to cfl3d',
     .            ' error check',/,
     .            ' (see main output file for error message)',/,
     .            ' (also check dynamic patch output file if',
     .            ' applicable)',/)
               end if
            end if
         else
            write(99,4)
    4       format(/,' abnormal termination due to receipt of',
     .      ' system signal',/,
     .      ' (kill, floating pt. exception,',
     .      ' segmentation fault, etc.)',/)
         end if
c
c        dump all internal buffers to make sure any relevant
c        info as to the cause of the abnormal termination 
c        is available to the user
c
         do nn=1,nbuf
            if (nou(nn).gt.0) then
               write(99,*)
               if (nn.eq.1) then
                  write(99,*)'dump of unit 11 (main output) buffer:'
               else if (nn.eq.2) then
                  write(99,*)'dump of unit 09 (fort.9, dyn. patch) ',
     .                       'buffer:'
               else if (nn.eq.3) then
                  write(99,*)'dump of unit 14 (baldwin-lomax) buffer:'
               else if (nn.eq.4) then
                  write(99,*)'dump of unit 25 (cfl3d.dynamic_patch) ',
     .                       'buffer:'
               else
                  write(99,*)'need a message for buffer',nn
               end if
               write(99,*)
               do kou = 1,nou(nn)
                  call outbuf(bou(kou,nn),99)
               end do
            end if
         end do
c
      end if
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
      call flush_(66)
      call flush_(96)
      call flush_(97)
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
      call flush(66)
      call flush(96)
      call flush(97)
      call flush(99)
#   endif
c
      if (myid .eq. myid_stop) then
         close(1)
         close(3)
         close(4)
         close(11)
         close(12)
         close(13)
         close(15)
         close(17)
         close(20)
         close(23)
         close(24)
         close(25)
         close(66)
         close(96)
         close(97)
         close(99)
      end if
#endif
c
#if defined DIST_MPI
c
      if (ierrflg.eq.0) then
         call MPI_Finalize (ierr)
      else if (ierrflg.eq.-99) then
         write (6,9) myid_stop
    9    format('node',i4,' is terminating the program ',
     .   'due to a precfl3d error check')
         write (6,*) 'see file precfl3d.error'
         call MPI_ABORT(MPI_COMM_WORLD, myid_stop, mpierror)
         call MPI_Finalize (ierr)
      else if (ierrflg.gt.-99) then
         write (6,10) myid_stop
   10    format('node',i4,' is terminating the program ',
     .   'due to a cfl3d error check')
         write (6,*) 'see file cfl3d.error'
         call MPI_ABORT(MPI_COMM_WORLD, myid_stop, mpierror)
         call MPI_Finalize (ierr)
      else if (ierrflg.eq.-999) then
         write (6,11) myid_stop
   11    format('node',i4,' is terminating the program ',
     .   'due to receipt of a system signal')
         write (6,*) 'see files precfl3d.error and cfl3d.error'
         call MPI_ABORT(MPI_COMM_WORLD, myid_stop, mpierror)
         call MPI_Finalize (ierr)
      end if
c
#else
      if (ierrflg.eq.0) then
         stop
      else if (ierrflg.eq.-99) then
         write (6,*) 'program termination ',
     .   'due to a precfl3d error check'
         write (6,*) 'see file precfl3d.error'
          stop
      else if (ierrflg.gt.-99) then
         write (6,*) 'program termination ',
     .   'due to a cfl3d error check'
         write (6,*) 'see file cfl3d.error'
         stop
      else if (ierrflg.eq.-999) then
         write (6,*) 'program termination ',
     .   'due to receipt of a system signal'
         write (6,*) 'see files precfl3d.error and cfl3d.error'
         stop
      end if
#endif
c
      stop
      end
