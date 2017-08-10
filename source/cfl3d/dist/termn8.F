c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
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
      call my_flush(66)
      call my_flush(96)
      call my_flush(97)
      call my_flush(98)
      call my_flush(99)
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
         close(98)
         close(99)
      end if
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
