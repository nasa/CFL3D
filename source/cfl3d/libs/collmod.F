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
      subroutine collmod(xmdj,xmdk,xmdi,xmdjc,xmdkc,xmdic,jdim,
     .                   kdim,idim,jj2,kk2,ii2,nm,nmds,iaes,maxaes)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict modal surface definitions to coarser meshes
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xmdj(kdim,idim,6,nmds,maxaes),
     .          xmdk(jdim,idim,6,nmds,maxaes),
     .          xmdi(jdim,kdim,6,nmds,maxaes)
      dimension xmdjc(kk2,ii2,6,nmds,maxaes),
     .          xmdkc(jj2,ii2,6,nmds,maxaes),
     .          xmdic(jj2,kk2,6,nmds,maxaes)
c
c     j-surfaces
c
      ii   = 0
      iinc = 2
      if (idim.eq.2) iinc = 1
      do i=1,idim,iinc
         ii = ii+1
         kk = 0
         do k=1,kdim,2
            kk = kk+1
            do ll=1,6
               xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes)
               xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes)
               xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes)
            end do
         end do
      end do
c
c     k-surfaces
c
      ii   = 0
      iinc = 2
      if (idim.eq.2) iinc = 1
      do i=1,idim,iinc
         ii = ii+1
         jj = 0
         do j=1,jdim,2
            jj = jj+1
            do ll=1,6
               xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes)
               xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes)
               xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes)
            end do
         end do
      end do
c
c     i-surfaces
c
      kk   = 0
      do k=1,kdim,2
         kk = kk+1
         jj = 0
         do j=1,jdim,2
            jj = jj+1
            do ll=1,6
               xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes)
               xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes)
               xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes)
            end do
         end do
      end do
c
      return
      end
