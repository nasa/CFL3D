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
      subroutine colldat(bcdata,mdim,ndim,bcdatac,mdimc,ndimc)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict auxilary boundary condition data arrays to 
c     coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension bcdata(mdim,ndim,2,12),bcdatac(mdimc,ndimc,2,12)
c
      do 10 l=1,12
      mm   = 0
      minc = 2
      mp   = 1
      if (mdim.eq.1) then
         minc = 1
         mp = 0
      end if
      do 10 m=1,mdim,minc
      mm   = mm+1
      mm   = min(mm,mdimc)
      nn   = 0
      ninc = 2
      np   = 1
      if (ndim.eq.1) then
         ninc = 1
         np = 0
      end if
      do 10 n=1,ndim,ninc
      nn   = nn+1
      nn   = min(nn,ndimc)
      nnn  = n+np
      nnn  = min(nnn,ndim)
      mmm  = m+mp
      mmm  = min(mmm,mdim)
      bcdatac(mm,nn,1,l) =
     .   .25*(bcdata(m,n,1,l) + bcdata(mmm,n,1,l) +
     .        bcdata(m,nnn,1,l) + bcdata(mmm,nnn,1,l))
      bcdatac(mm,nn,2,l) =
     .   .25*(bcdata(m,n,2,l) + bcdata(mmm,n,2,l) +
     .        bcdata(m,nnn,2,l) + bcdata(mmm,nnn,2,l))
   10 continue
      return 
      end
