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
      subroutine bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,jcsf,
     .                   kcsi,kcsf,icsi,icsf,jdim,kdim,idim,maxbl,
     .                   maxsegdg,nsegdfrm)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Install boundary deltas from forced oscillation or
c     aeroelastic motion into the 3d displacement array
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension dx(jdim,kdim,idim),dy(jdim,kdim,idim),dz(jdim,kdim,idim)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl)
c
      do iseg=1,nsegdfrm(nbl)
c
         if (icsi(nbl,iseg).eq.icsf(nbl,iseg)) then
            i = icsi(nbl,iseg)
            if (icsi(nbl,iseg).eq.1) then
               ii = 1
            else
               ii = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)
               do k=kcsi(nbl,iseg),kcsf(nbl,iseg)
                  dx(j,k,i)    = delti(j,k,1,ii)
                  dy(j,k,i)    = delti(j,k,2,ii)
                  dz(j,k,i)    = delti(j,k,3,ii)
               end do
            end do
         end if
c
         if (jcsi(nbl,iseg).eq.jcsf(nbl,iseg)) then
            j = jcsi(nbl,iseg)
            if (jcsi(nbl,iseg).eq.1) then
               jj = 1
            else
               jj = 2
            end if
            do k=kcsi(nbl,iseg),kcsf(nbl,iseg)
               do i=icsi(nbl,iseg),icsf(nbl,iseg)
                  dx(j,k,i)    = deltj(k,i,1,jj)
                  dy(j,k,i)    = deltj(k,i,2,jj)
                  dz(j,k,i)    = deltj(k,i,3,jj)
               end do
            end do
         end if
         if (kcsi(nbl,iseg).eq.kcsf(nbl,iseg)) then
            k = kcsi(nbl,iseg)
            if (kcsi(nbl,iseg).eq.1) then
               kk = 1
            else
               kk = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)
               do i=icsi(nbl,iseg),icsf(nbl,iseg)
                  dx(j,k,i)    = deltk(j,i,1,kk)
                  dy(j,k,i)    = deltk(j,i,2,kk)
                  dz(j,k,i)    = deltk(j,i,3,kk)
               end do
            end do
         end if
c
      end do
c
      return
      end
