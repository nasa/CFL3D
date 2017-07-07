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
      subroutine tau2x(jdim,kdim,idim,res,qr,js,ks,is,je,ke,ie,kode)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Put the restricted residual from a finer embedded
c     mesh into a coarser mesh.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qr(je-js+1,ke-ks+1,ie-is,5)
      dimension res(jdim,kdim,idim-1,5)
c
c      put restricted  r from finer embedded mesh into coarser mesh
c
c      jdim,kdim,idim  coarser mesh indices
c      js,ks,is        coarser mesh starting indices
c      je,ke,ie        coarser mesh ending indices
c
      if (kode.ge.2) then
         do 40 n=1,5
         kk = 0
         do 40 k=ks,ke-1
         kk = kk+1
         ii = 0
         do 40 i=is,ie-1
         ii = ii+1
         jj = 0
         do 40 j=js,je-1
         jj = jj+1
         qr(jj,kk,ii,n) =  qr(jj,kk,ii,n) - res(j,k,i,n)  
   40    continue
      end if
c
      do 60 n=1,5
      kk = 0
      do 60 k=ks,ke-1
      kk = kk+1
      ii = 0
      do 60 i=is,ie-1
      ii = ii+1
      jj = 0
      do 60 j=js,je-1
      jj = jj+1
      res(j,k,i,n) = res(j,k,i,n) + qr(jj,kk,ii,n)
   60 continue
      return
      end
