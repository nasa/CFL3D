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
      subroutine collxtb(xtt,xttt,mdim,ndim,mm2,nn2,nbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict xtb and atb arrays containing grid boundary 
c     velocity and acceleration to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xtt(mdim,ndim,3,2),xttt(mm2,nn2,3,2)
c
c      restrict xtt  to coarser mesh
c
c      mdim,ndim  finer mesh
c      mm2,nn2    coarser mesh
c
      ninc = 2
      if (ndim.eq.2) ninc = 1
      minc = 2
      if (mdim.eq.2) minc = 1
      do 10 ii=1,2
      do 10 ll=1,3
      nn   = 0
      do 10 n=1,ndim,ninc
      nn   = nn+1
      mm   = 0
      do 10 m=1,mdim,minc
      mm   = mm+1
      xttt(mm,nn,ll,ii) = xtt(m,n,ll,ii)
   10 continue
      return
      end
