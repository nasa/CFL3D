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
      subroutine rsmooth(eps,idim,jdim,kdim,icall,dq,d,nou,bou,nbuf,
     .                   ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Implicit residual smoothing (constant coefficient).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension dq(jdim,kdim,icall,5),d(jdim,idim+kdim),eps(3)
c
      common /sklton/ isklton
c
      smoopi = eps(1)
      smoopj = eps(2)
      smoopk = eps(3)
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''smoothing residuals - eps='',3e12.4)') 
     .   real(smoopi),real(smoopj),real(smoopk)
      end if
c
c--------- smoothing in i direction 
c
      if (idim.gt.2 .and. abs(real(smoopi)).gt.0.) then
         do 940 k=1,kdim-1
cdir$ ivdep
         do 10 j=1,jdim-1
         a           = smoopi
         t           = 1./(1.+a+a)
         d(j,1)      = t*a
         dq(j,k,1,1) = t*dq(j,k,1,1) 
         dq(j,k,1,2) = t*dq(j,k,1,2) 
         dq(j,k,1,3) = t*dq(j,k,1,3) 
         dq(j,k,1,4) = t*dq(j,k,1,4) 
         dq(j,k,1,5) = t*dq(j,k,1,5) 
   10    continue
         do 20 i=2,idim-1
cdir$ ivdep
         do 20 j=1,jdim-1 
         t           = 1./(1.+a+a -a*d(j,i-1))
         d(j,i)      = t*a
         dq(j,k,i,1) = t*(dq(j,k,i,1)  +a*dq(j,k,i-1,1)) 
         dq(j,k,i,2) = t*(dq(j,k,i,2)  +a*dq(j,k,i-1,2)) 
         dq(j,k,i,3) = t*(dq(j,k,i,3)  +a*dq(j,k,i-1,3)) 
         dq(j,k,i,4) = t*(dq(j,k,i,4)  +a*dq(j,k,i-1,4)) 
         dq(j,k,i,5) = t*(dq(j,k,i,5)  +a*dq(j,k,i-1,5)) 
   20    continue
         do 40 i=idim-2,1,-1
cdir$ ivdep
         do 40 j=1,jdim-1
         dq(j,k,i,1) = dq(j,k,i,1)  +d(j,i)*dq(j,k,i+1,1)
         dq(j,k,i,2) = dq(j,k,i,2)  +d(j,i)*dq(j,k,i+1,2)
         dq(j,k,i,3) = dq(j,k,i,3)  +d(j,i)*dq(j,k,i+1,3)
         dq(j,k,i,4) = dq(j,k,i,4)  +d(j,i)*dq(j,k,i+1,4)
         dq(j,k,i,5) = dq(j,k,i,5)  +d(j,i)*dq(j,k,i+1,5)
   40    continue
  940    continue
      end if
c
c--------- smoothing in j direction 
c
      if (jdim.gt.2 .and. abs(real(smoopj)).gt.0.) then
         do 970 i=1,idim-1
cdir$ ivdep
         do 50 k=1,kdim-1
         a           = smoopj
         t           = 1./(1.+a+a)
         d(1,k)      = t*a
         dq(1,k,i,1) = t*dq(1,k,i,1) 
         dq(1,k,i,2) = t*dq(1,k,i,2) 
         dq(1,k,i,3) = t*dq(1,k,i,3) 
         dq(1,k,i,4) = t*dq(1,k,i,4) 
         dq(1,k,i,5) = t*dq(1,k,i,5) 
   50    continue
         do 60 j=2,jdim-1
cdir$ ivdep
         do 60 k=1,kdim-1
         t           = 1./(1.+a+a-a*d(j-1,k))
         d(j,k)      = t*a
         dq(j,k,i,1) = t*(dq(j,k,i,1)  +a*dq(j-1,k,i,1)) 
         dq(j,k,i,2) = t*(dq(j,k,i,2)  +a*dq(j-1,k,i,2)) 
         dq(j,k,i,3) = t*(dq(j,k,i,3)  +a*dq(j-1,k,i,3)) 
         dq(j,k,i,4) = t*(dq(j,k,i,4)  +a*dq(j-1,k,i,4)) 
         dq(j,k,i,5) = t*(dq(j,k,i,5)  +a*dq(j-1,k,i,5)) 
   60    continue
         do 70 j=jdim-2,1,-1
cdir$ ivdep
         do 70 k=1,kdim-1
         dq(j,k,i,1) = dq(j,k,i,1)  +d(j,k)*dq(j+1,k,i,1)
         dq(j,k,i,2) = dq(j,k,i,2)  +d(j,k)*dq(j+1,k,i,2)
         dq(j,k,i,3) = dq(j,k,i,3)  +d(j,k)*dq(j+1,k,i,3)
         dq(j,k,i,4) = dq(j,k,i,4)  +d(j,k)*dq(j+1,k,i,4)
         dq(j,k,i,5) = dq(j,k,i,5)  +d(j,k)*dq(j+1,k,i,5)
   70    continue
  970    continue
      end if
c
c--------- smoothing in k direction 
c
      if (kdim.gt.2 .and. abs(real(smoopk)).gt.0.) then
         do 1100 i=1,idim-1
         do 80  j=1,jdim-1
         a           = smoopk
         t           = 1./(1.+a+a)
         d(j,1)      = t*a
         dq(j,1,i,1) = t*dq(j,1,i,1) 
         dq(j,1,i,2) = t*dq(j,1,i,2) 
         dq(j,1,i,3) = t*dq(j,1,i,3) 
         dq(j,1,i,4) = t*dq(j,1,i,4) 
         dq(j,1,i,5) = t*dq(j,1,i,5) 
   80    continue
         do 90 k=2,kdim-1
cdir$ ivdep
         do 90 j=1,jdim-1
         t           = 1./(1.+a+a-a*d(j,k-1))
         d(j,k)      = t*a
         dq(j,k,i,1) = t*(dq(j,k,i,1)  +a*dq(j,k-1,i,1)) 
         dq(j,k,i,2) = t*(dq(j,k,i,2)  +a*dq(j,k-1,i,2)) 
         dq(j,k,i,3) = t*(dq(j,k,i,3)  +a*dq(j,k-1,i,3)) 
         dq(j,k,i,4) = t*(dq(j,k,i,4)  +a*dq(j,k-1,i,4)) 
         dq(j,k,i,5) = t*(dq(j,k,i,5)  +a*dq(j,k-1,i,5)) 
   90    continue
         do 100 k=kdim-2,1,-1
cdir$ ivdep
         do 100 j=1,jdim-1
         dq(j,k,i,1) = dq(j,k,i,1)  +d(j,k)*dq(j,k+1,i,1)
         dq(j,k,i,2) = dq(j,k,i,2)  +d(j,k)*dq(j,k+1,i,2)
         dq(j,k,i,3) = dq(j,k,i,3)  +d(j,k)*dq(j,k+1,i,3)
         dq(j,k,i,4) = dq(j,k,i,4)  +d(j,k)*dq(j,k+1,i,4)
         dq(j,k,i,5) = dq(j,k,i,5)  +d(j,k)*dq(j,k+1,i,5)
  100    continue
 1100    continue
      end if
c
      return
      end
