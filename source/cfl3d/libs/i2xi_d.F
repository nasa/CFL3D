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
      subroutine i2xi_d(jc,kc,ic,qc,jf,kf,if,qif,js,ks,is,
     .                  je,ke,ie,nblc,ldim,nbl,bcif,nface)
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
c      planes of constant i-index
c
c      jc,kc,ic    : dimension of coarser mesh (ic is unused)
c      qc          : 3-plane subset of q-array coarser mesh
c      jf,kf,if    : dimension of finer mesh
c      qif         : q-array for interpolated points of finer mesh
c      js,ks,is    : starting indices of coarser mesh grid points
c                    defining boundary of finer mesh (is is unused)
c      je,ke,ie    : ending indices of coarser mesh grid points
c                    defining boundary of finer mesh (ie is unused)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,kc,3,ldim)
      dimension qif(jf,kf,ldim,4)
      dimension bcif(jf,kf,2)
      dimension q(3)
c
      kem = ke-1
      jem = je-1
c
      f1  = .75*.75
      f2  = .75*.25
      f4  = .25*.25
c
c     i = constant plane
c
      if (nface.eq.1) then
c
c     interpolate left boundary
c
         do 300 l=1,ldim
         do 300 k=ks,kem
         do 300 j=js,jem
         do 300 jl=1,2
         jj = (j-js)*2+jl
         j2 = max(j-1+(jl-1)*2,1)
         j2 = min(jc-1,j2)
         do 300 kl=1,2
         kk = (k-ks)*2+kl
         k2 = max(k-1+(kl-1)*2,1)
         k2 = min(kc-1,k2)
         do 301 i=1,3
         q(i) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k2,i,l)+qc(j2,k,i,l))
     .        + f4*qc(j2,k2,i,l)
  301    continue
         qif(jj,kk,l,1) = .25*q(1)+.75*q(2)
         qif(jj,kk,l,2) = .75*q(2)+.25*q(3)
         bcif(jj,kk,1) = 0.0
  300    continue
      end if
c
      if (nface.eq.2) then
c
c     interpolate right boundary
c
         do 310 l=1,ldim
         do 310 k=ks,kem
         do 310 j=js,jem
         do 310 jl=1,2
         jj = (j-js)*2+jl
         j2 = max(j-1+(jl-1)*2,1)
         j2 = min(jc-1,j2)
         do 310 kl=1,2
         kk = (k-ks)*2+kl
         k2 = max(k-1+(kl-1)*2,1)
         k2 = min(kc-1,k2)
         do 311 i=1,3
         q(i) = f1*qc(j,k,i,l)
     .        + f2*(qc(j,k2,i,l)+qc(j2,k,i,l))
     .        + f4*qc(j2,k2,i,l)
  311    continue
         qif(jj,kk,l,3) = .25*q(1)+.75*q(2)
         qif(jj,kk,l,4) = .75*q(2)+.25*q(3)
         bcif(jj,kk,2) = 0.0
  310    continue
      end if
c
c     **for safety**
c
      do 30 m=1,4
      do 30 l=1,ldim
      do 20 k=1,kf-1
      qif(jf,k,l,m) = qif(jf-1,k,l,m)
   20 continue
      do 30 j=1,jf-1
      qif(j,kf,l,m) = qif(j,kf-1,l,m)
   30 continue
c
      return
      end
