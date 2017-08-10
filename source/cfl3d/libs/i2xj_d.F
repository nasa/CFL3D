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
      subroutine i2xj_d(jc,kc,ic,qc,jf,kf,if,qjf,js,ks,is,
     .                  je,ke,ie,nblc,ldim,nbl,bcjf,nface)
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
c      planes of constant j-index
c
c      jc,kc,ic    : dimension of coarser mesh (jc is unused)
c      qc          : 3-plane subset of q-array coarser mesh
c      jf,kf,if    : dimension of finer mesh
c      qjf         : q-array for interpolated points of finer mesh
c      js,ks,is    : starting indices of coarser mesh grid points
c                    defining boundary of finer mesh (js is unused)
c      je,ke,ie    : ending indices of coarser mesh grid points
c                    defining boundary of finer mesh (je is unused)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(3,kc,ic,ldim)
      dimension qjf(kf,if-1,ldim,4)
      dimension bcjf(kf,if-1,2)
      dimension q(3)
c
      kem = ke-1
      iem = ie-1
c
      f1  = .75*.75
      f2  = .75*.25
      f4  = .25*.25
c
c     j = constant planes
c
      if (nface.eq.3) then
c
c     interpolate left boundary
c
         do 100 l=1,ldim
         do 100 i=is,iem
         do 100 k=ks,kem
         do 100 kl=1,2
         kk = (k-ks)*2+kl
         k2 = max(k-1+(kl-1)*2,1)
         k2 = min(kc-1,k2)
         do 100 il=1,2
         ii = (i-is)*2+il
         i2 = max(i-1+(il-1)*2,1)
         i2 = min(ic-1,i2)
         do 101 j=1,3
         q(j) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k,i2,l)+qc(j,k2,i,l))
     .        + f4*qc(j,k2,i2,l)
  101    continue
         qjf(kk,ii,l,1) = .25*q(1)+.75*q(2)
         qjf(kk,ii,l,2) = .75*q(2)+.25*q(3)
         bcjf(kk,ii,1) = 0.0
  100    continue
      end if      
c
      if (nface.eq.4) then
c
c     interpolate right boundary
c
         do 110 l=1,ldim
         do 110 i=is,iem
         do 110 k=ks,kem
         do 110 kl=1,2
         kk = (k-ks)*2+kl
         k2 = max(k-1+(kl-1)*2,1)
         k2 = min(kc-1,k2)
         do 110 il=1,2
         ii = (i-is)*2+il
         i2 = max(i-1+(il-1)*2,1)
         i2 = min(ic-1,i2)
         do 111 j=1,3
         q(j) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k,i2,l)+qc(j,k2,i,l))
     .        + f4*qc(j,k2,i2,l)
  111    continue
         qjf(kk,ii,l,3) = .25*q(1)+.75*q(2)
         qjf(kk,ii,l,4) = .75*q(2)+.25*q(3)
         bcjf(kk,ii,2) = 0.0
  110    continue
      end if
c
c     **for safety**
c
      do 30 m=1,4
      do 30 l=1,ldim
      do 10 i=1,if-1
      qjf(kf,i,l,m) = qjf(kf-1,i,l,m)
   10 continue
   30 continue
c
      return
      end
