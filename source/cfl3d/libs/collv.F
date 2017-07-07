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
      subroutine collv(vol,volc,nj,nk,ni,jj2,kk2,ii2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict volumes to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension vol(nj,nk,ni-1),volc(jj2,kk2,ii2-1)
c
      nj1 = nj-1
      nk1 = nk-1
      ni1 = ni-1
      jjl = jj2-1
      kkl = kk2-1
      iil = ii2-1
      ii  = 0
      if (ni.gt.2) then
         do 10 i=1,ni1,2
         ii  = ii+1
         kk  = 0
         do 10 k=1,nk1,2
         kk  = kk+1
         jj  = 0
         do 10 j=1,nj1,2
         jj  = jj+1
         volc(jj,kk,ii) = vol(j,k,i)+vol(j+1,k,i)
     .                   +vol(j,k+1,i)+vol(j+1,k+1,i)
     .                   +vol(j,k,i+1)+vol(j+1,k,i+1)
     .                   +vol(j,k+1,i+1)+vol(j+1,k+1,i+1)
   10    continue
      else
         i  = 1
         ii = 1
         kk = 0
         do 710 k=1,nk1,2
         kk = kk+1
         jj = 0
         do 710 j=1,nj1,2
         jj = jj+1
         volc(jj,kk,ii) = vol(j,k  ,i)+vol(j+1,k  ,i)
     .                   +vol(j,k+1,i)+vol(j+1,k+1,i)
  710    continue
      end if
      do 13 k=1,kkl
      do 13 i=1,iil
      volc(jj2,k,i) = volc(jjl,k,i)
   13 continue
      do 14 i=1,iil
      do 14 j=1,jj2
      volc(j,kk2,i) = volc(j,kkl,i)
   14 continue
      return
      end
