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
c     Purpose: output error flag to signal program status at termination
c
c      error codes (ierrflg):
c
c        0  normal termination
c
c       -1  input file error
c       -2  no cfl3d input file specified
c       -3  parameter mbloc too small
c       -4  parameter npmax too small
c       -5  inconsistant grid and grid sensitivity files
c       -6  internal file conversion error
c       -7  parameter mseqn too small
c
c      -99  miscellaneous error
c***********************************************************************
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      write(99,99) ierrflg
   99 format('error code:',/,i3)
c
c*******************************
c     normal termination section
c*******************************
c
      if (ierrflg.eq.0) then
         write(99,1)
    1    format('execution terminated normally')
         close(8,status='delete')
         close(20,status='delete')
         close(10)
         close(7)
         close(15)
         close(9)
         close(30)
         close(40)
         close(50)
         close(60)
         close(99)
         stop
      end if
c
c*********************************
c     abnormal termination section
c*********************************
c
      if (ierrflg.lt.0) then
         write(99,100)
  100    format('abnormal termination:')
      end if
      if (ierrflg.eq.-1) then
         write(99,101) 
  101    format('input file error')
      else if (ierrflg.eq.-2) then
         write(99,102)
  102    format('no cfl3d input file specified')
      else if (ierrflg.eq.-3) then
         write(99,103)
  103    format('parameter mbloc too small')
      else if (ierrflg.eq.-4) then
         write(99,104)
  104    format('parameter npmax too small')
      else if (ierrflg.eq.-5) then
         write(99,105)
  105    format('inconsistant grid and grid sensitivity files')
      else if (ierrflg.eq.-6) then
         write(99,106)
  106    format('internal file conversion error')
      else if (ierrflg.eq.-7) then
         write(99,107)
  107    format('parameter mseqn too small')
      else
         write(99,999)
  999    format('miscellaneous error')
      end if 
c
      do nn=1,nbuf
         if (nou(nn).gt.0) then
            write(99,*)
            do kou = 1,nou(nn)
               call outbuf(bou(kou,nn),99)
            end do
         end if
      end do
c
      write (6, *)  'terminating the program'
c
      if (ierrflg.eq.-7) then
c        conversion error: close, but don't remove temporary files
         close(8)
         close(20)
      else
c        close and delete temporary files
         close(8,status='delete')
         close(20,status='delete')
      end if
      close(10)
      close(7)
      close(15)
      close(9)
      close(30)
      close(40)
      close(50)
      close(60)
      close(99)
c
      stop
      end
