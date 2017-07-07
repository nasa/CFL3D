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
      subroutine xmukin(n,temp,visc,tinf)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Computes Sutherland's formula.  Note that this routine
c     is only called by ctime.  (Mu is generally computed in-line
c     throughout the code using this formula.)
c     Also note that an older version of this routine used a linear
c     law for low temperatures, but it was not consistently used
c     throughout the whole code, and so has been taken out.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension temp(n),visc(n)
c
      common /fluid2/ pr,prt,cbar
c
c      molecular viscosity
c
      c2b    = cbar/tinf
      c2bp   = c2b+1.e0
c
      do 1000 j=1,n
         t5      =  temp(j)
         t6      =  sqrt(t5)
         visc(j) =  c2bp*t5*t6/(c2b+t5)
 1000 continue
      return
      end
