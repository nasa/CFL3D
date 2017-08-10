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
      subroutine loadgr(w,mgwk,lx,ly,lz,jindex,x,y,z,mdim,ndim,
     .                  idimg,jdimg,kdimg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Load the proper grid from 1-d storage array to 2-d
c     work array.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension w(mgwk),x(mdim,ndim),y(mdim,ndim),z(mdim,ndim)
c      
c     patch surface is an i=constant surface
c
      if (jindex/10.eq.1) then
         if (jindex.eq.11) then
            i = 1
         else
            i = idimg
         end if
         do 10 j=1,jdimg
         do 10 k=1,kdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(j,k) = w(lx+l)
         y(j,k) = w(ly+l)
         z(j,k) = w(lz+l)
   10    continue      
      end if
c      
c     patch surface is a j=constant surface
c
      if (jindex/10.eq.2) then
         if (jindex.eq.21) then
            j = 1
         else
            j = jdimg
         end if
         do 20 i=1,idimg
         do 20 k=1,kdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(k,i) = w(lx+l)
         y(k,i) = w(ly+l)
         z(k,i) = w(lz+l)
   20    continue      
      end if      
c      
c     patch surface is a k=constant surface
c
      if (jindex/10.eq.3) then
         if (jindex.eq.31) then
            k = 1
         else
            k = kdimg
         end if
         do 30 i=1,idimg
         do 30 j=1,jdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(j,i) = w(lx+l)
         y(j,i) = w(ly+l)
         z(j,i) = w(lz+l)
   30    continue      
      end if
      return
      end
