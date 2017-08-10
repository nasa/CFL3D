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
      subroutine grdmove(nbl,jdim,kdim,idim,x,y,z,xorig,yorig,
     .                  zorig,xorg,yorg,zorg,thetax,thetay,
     .                  thetaz)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Move the grid from one position to another. On input,
c     (x,y,z) are the grid coordinates in the old position, with origin
c     (xorig,yorig,zorig). Given the new origin (xorg,yorg,zorg) and
c     rotational displacements (thetax,thetay,thetaz) the grid is first
c     rotated and then translated to the new position; on return,
c     (x,y,z) correspond to the new position. No translation occurs if
c     (xorg,yorg,zorg) = (xorig,yorig,zorig); no rotation occurs if
c     (thetax,thetay,thetaz) = (0,0,0)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
c     rotation about x-axis
      sa = sin(thetax)
      ca = cos(thetax)
      do 20 j=1,jdim
      do 20 k=1,kdim
      do 20 i=1,idim
      yy = (y(j,k,i)-yorig)*ca-(z(j,k,i)-zorig)*sa+yorig
      zz = (y(j,k,i)-yorig)*sa+(z(j,k,i)-zorig)*ca+zorig
      y(j,k,i) = yy
      z(j,k,i) = zz
   20 continue
c     
c     rotation about y-axis
      sa = sin(thetay)
      ca = cos(thetay)
      do 30 j=1,jdim
      do 30 k=1,kdim
      do 30 i=1,idim
      xx =  (x(j,k,i)-xorig)*ca+(z(j,k,i)-zorig)*sa+xorig
      zz = -(x(j,k,i)-xorig)*sa+(z(j,k,i)-zorig)*ca+zorig
      x(j,k,i) = xx
      z(j,k,i) = zz
   30 continue
c     
c     rotation about z-axis
      sa = sin(thetaz)
      ca = cos(thetaz)
      do 40 j=1,jdim
      do 40 k=1,kdim
      do 40 i=1,idim
      xx = (x(j,k,i)-xorig)*ca-(y(j,k,i)-yorig)*sa+xorig
      yy = (x(j,k,i)-xorig)*sa+(y(j,k,i)-yorig)*ca+yorig
      x(j,k,i) = xx
      y(j,k,i) = yy
   40 continue
c
c     translation
      do 10 j=1,jdim
      do 10 k=1,kdim
      do 10 i=1,idim
      x(j,k,i) = x(j,k,i) + xorg - xorig
      y(j,k,i) = y(j,k,i) + yorg - yorig
      z(j,k,i) = z(j,k,i) + zorg - zorig
  10  continue
c
      return
      end
