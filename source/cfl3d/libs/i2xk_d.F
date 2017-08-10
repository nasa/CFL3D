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
      subroutine i2xk_d(jc,kc,ic,qc,jf,kf,if,qkf,js,ks,is,
     .                  je,ke,ie,nblc,ldim,nbl,bckf,nface)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Interpolate primative variables from coarser 
c     meshes onto twice finer meshes, for k=constant surfaces.
c     This version of i2x requires only a 3-plane subset of the 
c     full qc array to be stored.
c***********************************************************************
c
c      interpolate from coarser mesh onto twice finer mesh
c      planes of constant k-index
c
c      jc,kc,ic    : dimension of coarser mesh (kc is unused)
c      qc          : 3-plane subset of q-array coarser mesh
c      jf,kf,if    : dimension of finer mesh
c      qkf         : q-array for interpolated points of finer mesh
c      js,ks,is    : starting indices of coarser mesh grid points
c                    defining boundary of finer mesh (ks is unused)
c      je,ke,ie    : ending indices of coarser mesh grid points
c                    defining boundary of finer mesh (ke is unused)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,3,ic,ldim)
      dimension qkf(jf,if-1,ldim,4)
      dimension bckf(jf,if-1,2)
      dimension q(3)
c
      jem = je-1
      iem = ie-1
c
      f1  = .75*.75
      f2  = .75*.25
      f4  = .25*.25
c
c     k = constant planes
c
      if (nface.eq.5) then
c
c     interpolate left boundary
c
         do 200 l=1,ldim
         do 200 i=is,iem
         do 200 j=js,jem
         do 200 jl=1,2
         jj = (j-js)*2+jl
         j2 = max(j-1+(jl-1)*2,1)
         j2 = min(jc-1,j2)
         do 200 il=1,2
         ii = (i-is)*2+il
         i2 = max(i-1+(il-1)*2,1)
         i2 = min(ic-1,i2)
         do 201 k=1,3
         q(k) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k,i2,l)+qc(j2,k,i,l))
     .        + f4*qc(j2,k,i2,l)
  201    continue
         qkf(jj,ii,l,1) = .25*q(1)+.75*q(2)
         qkf(jj,ii,l,2) = .75*q(2)+.25*q(3)
         bckf(jj,ii,1) = 0.0
  200    continue
      end if
c
      if (nface.eq.6) then
c
c     interpolate right boundary
c
         do 210 l=1,ldim
         do 210 i=is,iem
         do 210 j=js,jem
         do 210 jl=1,2
         jj = (j-js)*2+jl
         j2 = max(j-1+(jl-1)*2,1)
         j2 = min(jc-1,j2)
         do 210 il=1,2
         ii = (i-is)*2+il
         i2 = max(i-1+(il-1)*2,1)
         i2 = min(ic-1,i2)
         do 211 k=1,3
         q(k) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k,i2,l)+qc(j2,k,i,l))
     .        + f4*qc(j2,k,i2,l)
  211    continue
         qkf(jj,ii,l,3) = .25*q(1)+.75*q(2)
         qkf(jj,ii,l,4) = .75*q(2)+.25*q(3)
         bckf(jj,ii,2) = 0.0
  210    continue
      end if
c
c     **for safety**
c
      do 30 m=1,4
      do 30 l=1,ldim
      do 10 i=1,if-1
      qkf(jf,i,l,m) = qkf(jf-1,i,l,m)
   10 continue
   30 continue
c
      return
      end
