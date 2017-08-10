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
