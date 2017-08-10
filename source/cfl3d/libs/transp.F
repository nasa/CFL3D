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
      subroutine transp(mdim,ndim,jmax1,kmax1,msub1,l,x1,y1,z1,
     .                  dx,dy,dz,intmx,int)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Translate "from" block to provide complete coverage for
c     interpolation for cases in which the complete physical domain is 
c     not modeled.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(mdim,ndim,msub1),y1(mdim,ndim,msub1),
     .          z1(mdim,ndim,msub1),dx(intmx,msub1),dy(intmx,msub1),
     .          dz(intmx,msub1)
c
      do 10 j=1,jmax1
      do 10 k=1,kmax1
      x1(j,k,l) = x1(j,k,l) + dx(int,l)
      y1(j,k,l) = y1(j,k,l) + dy(int,l)
      z1(j,k,l) = z1(j,k,l) + dz(int,l)
 10   continue
      return
      end
