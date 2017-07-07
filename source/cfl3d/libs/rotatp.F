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
      subroutine rotatp(mdim,ndim,jmax1,kmax1,msub1,l,x1,y1,z1,
     .                  dthetx,dthety,dthetz,xorig,yorig,zorig,
     .                  mbl,nn,intmx,int)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Rotate "from" block to provide complete coverage for
c     interpolation for cases in which the complete physical domain is
c     not modeled.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(mdim,ndim,msub1),y1(mdim,ndim,msub1),
     .          z1(mdim,ndim,msub1),dthetx(intmx,msub1),
     .          dthety(intmx,msub1),dthetz(intmx,msub1),
     .          xorig(nn),yorig(nn),zorig(nn)
c
      xorg = xorig(mbl)
      yorg = yorig(mbl)
      zorg = zorig(mbl)
      dthx = dthetx(int,l)
      dthy = dthety(int,l)
      dthz = dthetz(int,l)
c
      if (abs(real(dthx)) .gt. 0.) then
c 
c        rotation about a line parallel to the x-axis
c
         ca = cos(dthx)
         sa = sin(dthx)
         do 10 j=1,jmax1
         do 10 k=1,kmax1
         ytemp = (y1(j,k,l) - yorg)*ca - (z1(j,k,l) - zorg)*sa + yorg
         ztemp = (y1(j,k,l) - yorg)*sa + (z1(j,k,l) - zorg)*ca + zorg
         y1(j,k,l) = ytemp
         z1(j,k,l) = ztemp
 10      continue
      end if
c
      if (abs(real(dthy)) .gt. 0.) then
c
c        rotation about a line parallel to the y-axis
c
         ca = cos(dthy)
         sa = sin(dthy)
         do 20 j=1,jmax1
         do 20 k=1,kmax1
         xtemp =  (x1(j,k,l) - xorg)*ca + (z1(j,k,l) - zorg)*sa + xorg
         ztemp = -(x1(j,k,l) - xorg)*sa + (z1(j,k,l) - zorg)*ca + zorg
         x1(j,k,l) = xtemp
         z1(j,k,l) = ztemp
 20      continue
      end if
c
      if (abs(real(dthz)) .gt. 0.) then
c
c        rotation about a line parallel to the z-axis
c
         ca = cos(dthz)
         sa = sin(dthx)
         do 30 j=1,jmax1
         do 30 k=1,kmax1
         xtemp = (x1(j,k,l) - xorg)*ca - (y1(j,k,l) - yorg)*sa + xorg
         ytemp = (x1(j,k,l) - xorg)*sa + (y1(j,k,l) - yorg)*ca + yorg
         x1(j,k,l) = xtemp
         y1(j,k,l) = ytemp
 30      continue
      end if
c
      return
      end
