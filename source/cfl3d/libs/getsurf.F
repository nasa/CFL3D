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
      subroutine getsurf(x,y,z,deltj,deltk,delti,xnm1,ynm1,znm1,
     .                   icsi,icsf,jcsi,jcsf,kcsi,kcsf,wk,nwork,nbl,
     .                   idim,jdim,kdim,nsurf,nsurfb,nsegdfrm,maxbl,
     .                   idfrmseg,maxsegdg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set up a list (array) containing, for all deforming solid
c     surface points, the x,y,z values at the surface point, as well as
c     delta displacement at those points. the list structure is as 
c     folllows: each successive 9 entries in the list give,
c     in order, the x, y, z, deltx, delty, deltz, xnm1, ynm1,
c     znm1 values for the surface point (where xnm1, etc. are
c     are the surface points at time n-1); these 9 data are repeated 
c     for each solid surface point that undergoes deformation. Thus, 
c     if there are a total of nsurf solid surface points that
c     undergo deformation, the list will be of dimension 9*nsurf
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xnm1(jdim,kdim,idim),ynm1(jdim,kdim,idim),
     .          znm1(jdim,kdim,idim)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension wk(nwork),nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg)
c
      ll = nsurf*9
c
      do is=1,nsegdfrm(nbl)
c
c        put i-surface points in list
c
         if (icsi(nbl,is) .eq. icsf(nbl,is)) then
c
            i  = icsi(nbl,is)
            ii = 1
            if (icsi(nbl,is).eq.idim) ii=2
            do j=jcsi(nbl,is),jcsf(nbl,is)
               do k=kcsi(nbl,is),kcsf(nbl,is)
                  wk(ll+1) = x(j,k,i)
                  wk(ll+2) = y(j,k,i)
                  wk(ll+3) = z(j,k,i)
                  wk(ll+4) = delti(j,k,1,ii)
                  wk(ll+5) = delti(j,k,2,ii)
                  wk(ll+6) = delti(j,k,3,ii)
                  wk(ll+7) = xnm1(j,k,i)
                  wk(ll+8) = ynm1(j,k,i)
                  wk(ll+9) = znm1(j,k,i)
                  ll       = ll + 9
                  nsurfb   = nsurfb + 1
                end do
            end do
c
         end if
c
c        put j-surface points in list
c
         if (jcsi(nbl,is) .eq. jcsf(nbl,is)) then
c
            j  = jcsi(nbl,is)
            jj = 1
            if (jcsi(nbl,is).eq.jdim) jj=2
            do i=icsi(nbl,is),icsf(nbl,is)
               do k=kcsi(nbl,is),kcsf(nbl,is)
                  wk(ll+1) = x(j,k,i)
                  wk(ll+2) = y(j,k,i)
                  wk(ll+3) = z(j,k,i)
                  wk(ll+4) = deltj(k,i,1,jj)
                  wk(ll+5) = deltj(k,i,2,jj)
                  wk(ll+6) = deltj(k,i,3,jj)
                  wk(ll+7) = xnm1(j,k,i)
                  wk(ll+8) = ynm1(j,k,i)
                  wk(ll+9) = znm1(j,k,i)
                  ll       = ll + 9
                  nsurfb   = nsurfb + 1
                end do
            end do
c
         end if
c
c        put k-surface points in list
c
         if (kcsi(nbl,is) .eq. kcsf(nbl,is)) then
c
            k  = kcsi(nbl,is)
            kk = 1
            if (kcsi(nbl,is).eq.kdim) kk=2
            do j=jcsi(nbl,is),jcsf(nbl,is)
               do i=icsi(nbl,is),icsf(nbl,is)
                  wk(ll+1) = x(j,k,i)
                  wk(ll+2) = y(j,k,i)
                  wk(ll+3) = z(j,k,i)
                  wk(ll+4) = deltk(j,i,1,kk)
                  wk(ll+5) = deltk(j,i,2,kk)
                  wk(ll+6) = deltk(j,i,3,kk)
                  wk(ll+7) = xnm1(j,k,i)
                  wk(ll+8) = ynm1(j,k,i)
                  wk(ll+9) = znm1(j,k,i)
                  ll       = ll + 9
                  nsurfb   = nsurfb + 1
                end do
            end do
c
         end if
c
      end do
c
      return
      end
