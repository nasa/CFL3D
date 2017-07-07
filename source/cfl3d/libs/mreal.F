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
      subroutine mreal(xyz,mdim1,ndim1,mdim2,ndim2,x,y,z)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Unload grid data from transfer array xyz to x,y,z
c      work arrays
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xyz(mdim1,ndim1,3)
      dimension x(mdim2,ndim2),y(mdim2,ndim2),z(mdim2,ndim2)
c
c     initialize the individual x,y,z arrays, since they are typically
c     larger that the corresponding sections in xyz (i.e. mdim2*ndim2 
c     vs. mdim1*ndim1)
c
      do m=1,mdim2
         do n=1,ndim2
            x(m,n) = 0.
            y(m,n) = 0.
            z(m,n) = 0.
         end do
      end do
c
      do m=1,mdim1
         do n=1,ndim1
            x(m,n) = xyz(m,n,1)
            y(m,n) = xyz(m,n,2)
            z(m,n) = xyz(m,n,3)
         end do
      end do
c
      return
      end
