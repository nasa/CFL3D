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
#if defined ADPOFF
$autodbl off
#endif
      subroutine  cputim(icall)
c
c     $Id$
c
c=======================================================================
c
#if defined ADPOFF
      implicit real*8 (a-h,o-z)
#endif
c
c     common blocks preserve tim,tm,tsum between calls to this subroutine
      real*4 tim(3,3),tm(3)
      common/time1/ tim,tm
      dimension ia(3)
c
#if defined CRAY_TIME
      real*4 t0,t1
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
#if defined CRAY_TIME
c
c     timings for cray - note: user time is cpu time from the
c     intrinsic "second" function, system time set to zero(!), 
c     and wall clock time set to user time.
c
      t0 = second()
c
c     first time here
c
      if (tim(1,3).eq.0.) then
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
#if defined ETFLAG
      call etime_(tm)
      call itime_(ia)
#else
      call etime(tm)
      call itime(ia)
#endif
c
      tm(3) = ia(3) +60*(ia(2) +60*ia(1))
c
c     first time here
c
      if (tim(1,3).eq.0.) then
         do n=1,3
            tim(n,3) = tm(n)
         enddo
      endif
c
c     check if wall clock passed 24:00 (midnight)
c     since the last call to this routine
c
      if (tm(3).lt.tim(3,3)) tm(3) = tm(3) + 24*3600
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
      if (icall .lt. 0) then
c
c        output the collective timings for this run
c
         write (6,'(/,/,
     .   "                  time in seconds",/,/,
     .   "   node      user    system     total    wall clock")')
         write (6,'(" ",i6,3f10.2,f12.2)')
     .   0,tim(1,1),tim(2,1),tim(1,1)+tim(2,1),tim(3,1)
         totu = tim(1,1)
         tots = tim(2,1)
         tott = tim(1,1)+tim(2,1)
         write (6,'(" ------------------------------------")')
         write (6,'(" total:",3f10.2)')
     .   totu,tots,tott
         nhrs  = tim(3,1) / 3600
         nsecs = tim(3,1) - nhrs * 3600
         nmins = nsecs / 60
         nsecs = nsecs - nmins * 60
         write (6,'(/," total run (wall) time = ",i4," hours ",
     .   i4," minutes ",i4," seconds")') nhrs,nmins,nsecs
      end if
      return
      end
#if defined ADPOFF
$autodbl dblpad
#endif

