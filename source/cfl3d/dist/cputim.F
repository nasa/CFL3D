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
      subroutine  cputim(icall,nodes,string,myhost,myid,mycomm,iunit)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Determine timings for cpu usage...original version of 
c      this routine coutrtesy v. n. vatsa...mods by r. t. biedron
c
c      icall = 0 for initial call
c            < 0 for final call; total run time is output
c            + 1 for a call preceding an individual subroutine or
c                section of code to be timed, but after the initial
c                call; nothing is output
c            + 2 for a call following the individual subroutine or
c                section of code to be timed; time since last call
c                is output
c      nodes  = number of nodes, if distributed (may be 0 for icall=0/1)
c      string = a character string up to 50 characters, may be used to
c               identify particular item being timed
c***********************************************************************
c
#if defined ADP_OFF
#   ifdef CMPLX
#     ifdef DBLE_PRECSN
      implicit complex*8(a-h,o-z)
#     else
      implicit complex(a-h,o-z)
#     endif
#   else
#     ifdef DBLE_PRECSN
      implicit real*8 (a-h,o-z)
#     endif
#   endif
#else
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#endif
c
#if defined DIST_MPI
#     include "mpif.h"
#endif
c
#if defined CRAY
      real*4 t0,t1
#endif
#if defined DIST_MPI
      real*4 tim_nd(3,3)
#endif
      real*4 tim(3,3),tm(3)
      real*4 walltime,totaltime
c
      character*50 string
c
#if defined DIST_MPI
#   ifdef DBLE_PRECSN
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   else
#     define MY_MPI_REAL MPI_REAL
#   endif
      dimension istat(MPI_STATUS_SIZE)
#endif
      dimension ia(3)
c
c     common blocks time1/time2 preserve tim,tm,tim_nd between 
c     calls to this subroutine
c
      common/time1/ tim,tm
#if defined DIST_MPI
      common/time2/ tim_nd
#endif
c
c     etime returns elapsed time as:
c       tm(1) = user time
c       tm(2) = system time
c     itime returns wall clock time as:
c       ia(1) = hour (0-23)
c       ia(2) = minute (0-59)
c       ia(3) = second (0-59)
c
c     timing array modifed to track user and system time
c     tim(1,1)  = total user time
c     tim(2,1)  = total system time
c     tim(3,1)  = total wall clock time
c     tim(1,2)  = user time since last call to cputim()
c     tim(2,2)  = system time since last call to cputim()
c     tim(3,2)  = wall clock time since last call to cputim()
c     tim(1,3),tim(2,3) and tim(3,3) used for intermediate results
c
c     initialize 
c
      if (icall .eq. 0) then
         do j=1,3
            do i=1,3
               tim(i,j) = 0.
            enddo
         enddo
      end if

c
#if defined CRAY
c
c     timings for cray - note: user time is cpu time from the
c     intrinsic "second" function, system time set to zero(!), 
c     and wall clock time set to user time.
c
      t0 = second()
c
c     first time here
c
      if (real(tim(1,3)).eq.0.) then
         tim(1,3) = t0
         tim(2,3) = 0.
         tim(3,3) = t0
      endif
c
c     get totals
c
      do n=1,3,2
         tim(n,2) = t0 - tim(n,3)
         tim(n,1) = tim(n,1) + tim(n,2)
         tim(n,3) = t0
      enddo
c     zero out system time
      do n=1,3
         tim(2,n) = 0.
      enddo
c
#else
c
c     timing for workstations
c
#if defined IBM
      totaltime = etime_(tm)
      call itime_(ia)
#else
      totaltime = etime(tm)
      call itime(ia)
#endif
c
      tm(3) = ia(3) +60*(ia(2) +60*ia(1))
c
c     first time here
c
      if (real(tim(1,3)).eq.0.) then
         do n=1,3
            tim(n,3) = tm(n)
         enddo
      endif
c
c     check if wall clock passed 24:00 (midnight)
c     since the last call to this routine
c
      if (real(tm(3)).lt.real(tim(3,3))) tm(3) = tm(3) + 24*3600
c
c     get totals
c
      do n=1,3
         tim(n,2) = tm(n) - tim(n,3)
         tim(n,1) = tim(n,1) + tim(n,2)
         tim(n,3) = tm(n)
      enddo
#endif
c
c.....output the collective timings for this run
c
      if (icall .lt. 0) then
c
#if defined DIST_MPI
         if (myid.eq.myhost) then
c
#endif
c
c           best attempt at correct wall time for > 24  hour runs:
c           the problem seems to be that tim(3,2), supposedly the wall
c           time, gets reset after 24 hours. have to creative in the case
c           the total time - system+user - is not *quite* 24 hours and
c           the true wall time is just over 24 hours, so that ndays=0
c           but tim(3,2) is small. in that case use the total time as
c           the wall time.
c
            walltime  = tim(3,1)
            totaltime = tim(1,1) + tim(2,1)
            ndays = int(real(totaltime))/86400
            if (ndays .gt. 0) then
               walltime = walltime + ndays*86400
            else
               walltime = max(walltime,totaltime)
            end if
            write(iunit,'(/,a50,/,/,
     .      "   node      user    system     total    wall clock")')
     .      string
            write(iunit,'(" ",i6,3f10.2,f12.2)')
     .      0,real(tim(1,1)),real(tim(2,1)),
     .      real(tim(1,1))+real(tim(2,1)),real(walltime)
            totu = tim(1,1)
            tots = tim(2,1)
            tott = tim(1,1)+tim(2,1)
#if defined DIST_MPI
            do inode=1,nodes
               mytag = 1
               call MPI_Recv (tim_nd,9,MPI_REAL,inode,mytag,mycomm,
     .                        istat,ierr)
c
c              best attempt at correct wall time for > 24  hour runs
c
               walltime  = tim_nd(3,1)
               totaltime = tim_nd(1,1) + tim_nd(2,1)
               ndays = int(real(totaltime))/86400
               if (ndays .gt. 0) then
                  walltime = walltime + ndays*86400
               else
                  walltime = max(walltime,totaltime)
               end if
               write(iunit,'(" ",i6,3f10.2,f12.2)')
     .         inode,real(tim_nd(1,1)),real(tim_nd(2,1)),
     .         real(tim_nd(1,1))+real(tim_nd(2,1)),
     .         walltime
               totu = totu + tim_nd(1,1)
               tots = tots + tim_nd(2,1)
               tott = tott + tim_nd(1,1)+tim_nd(2,1)
            enddo
#endif
            write(iunit,'(" ------------------------------------")')
            write(iunit,'(" total:",3f10.2)')
     .      real(totu),real(tots),real(tott)
            nhrs  = walltime / 3600
            nsecs = walltime - nhrs * 3600
            nmins = nsecs / 60
            nsecs = nsecs - nmins * 60
            write(iunit,'(/," total run (wall) time = ",i4," hours ",
     .      i4," minutes ",i4," seconds")') nhrs,nmins,nsecs
#if defined DIST_MPI
c
         else
c
            mytag = 1
            call MPI_Send (tim,9,MPI_REAL,myhost,mytag,mycomm,ierr)
c
         endif
#endif
      end if
c
c.....print time since last call to this routine
c
      if (icall .eq. 2) then
c
#if defined DIST_MPI
         if (myid.eq.myhost) then
c
#endif
c
c           best attempt at correct wall time for > 24  hour runs
c
            walltime  = tim(3,2)
            totaltime = tim(1,1) + tim(2,1)
            ndays = int(real(totaltime))/86400
            if (ndays .gt. 0) then
               walltime = walltime + ndays*86400
            else
               walltime = max(walltime,totaltime)
            end if
            write(iunit,'(/,a50,/,/,
     .      "   node      user    system     total    wall clock")')
     .      string
            write(iunit,'(" ",i6,3f10.2,f12.2)')
     .      0,real(tim(1,2)),real(tim(2,2)),
     .      real(tim(1,2))+real(tim(2,2)),real(walltime)
            totu = tim(1,2)
            tots = tim(2,2)
            tott = tim(1,2)+tim(2,2)
#if defined DIST_MPI
c
            do inode=1,nodes
               mytag = 1
               call MPI_Recv (tim_nd,9,MPI_REAL,inode,mytag,mycomm,
     .                        istat,ierr)
c
c              best attempt at correct wall time for > 24  hour runs
c
               walltime  = tim_nd(3,2)
               totaltime = tim(1,1) + tim(2,1)
               ndays = int(real(totaltime))/86400
               if (ndays .gt. 0) then
                  walltime = walltime + ndays*86400
               else
                  walltime = max(walltime,totaltime)
               end if
               write(iunit,'(" ",i6,3f10.2,f12.2)')
     .         inode,real(tim_nd(1,2)),real(tim_nd(2,2)),
     .         real(tim_nd(1,2))+real(tim_nd(2,2)),real(walltime)
               totu = totu + tim_nd(1,2)
               tots = tots + tim_nd(2,2)
               tott = tott + tim_nd(1,2)+tim_nd(2,2)
            enddo
#endif
           write(iunit,'(" ------------------------------------")')
           write(iunit,'(" total:",3f10.2)')
     .     real(totu),real(tots),real(tott)
#if defined DIST_MPI
c
         else
c
            mytag = 1
            call MPI_Send (tim,9,MPI_REAL,myhost,mytag,mycomm,ierr)
c
         end if
c
#endif
      end if
c
      return
      end
