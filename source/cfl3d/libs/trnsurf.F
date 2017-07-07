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
      subroutine trnsurf(jdim,kdim,idim,x,y,z,deltj,deltk,delti,nbl,
     .                   idef,xorg,yorg,zorg,utran,vtran,wtran,
     .                   rfreqt,ici,icf,jci,jcf,kci,kcf,time,nou,bou,
     .                   nbuf,ibufdim,myid,wkj,wkk,wki,xorg0,yorg0,
     .                   zorg0)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determines increment to delta displacement due to 
c     surface translation.
c
c     idefrm...modulation for mesh deformation
c              = 0 no deformation
c              = 1 sinusoidal variation of surface translation
c              = 2 sinusoidal variation of surface rotation
c              = 999 block undergoes deformation, but not by surface
c                rotation or translation (not handled by this routine)
c
c     surface translation set in the range (ici,icf), (jci,jcf), (kci,kcf)
c     one pair of the indicies must be identical; this set of constant
c     indicies determines which surface in the grid is deformed 
c
c     deltj/k/i...arrays for storage of delta displacements due to
c                 surface translation; upon entering this routine, 
c                 these are assumed to be zero.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf),wkj(kdim,idim,2),wkk(jdim,idim,2),
     .          wki(jdim,kdim,2)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     .          deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
c
      common /sklton/ isklton
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     determine which block face is to be translated
c
      if (ici .eq. icf) then
         isurf = 1
         ii    = ici
         ll    = 1
         if (ici.eq.idim) ll = 2
      else if (jci .eq. jcf) then
         isurf = 2
         jj    = jci
         ll    = 1
         if (jci.eq.jdim) ll = 2
      else if (kci .eq. kcf) then
         isurf = 3
         kk    = kci
         ll    = 1
         if (kci.eq.kdim) ll = 2
      end if
c
c     ft modulates the translation
c     dfdt is the time derivative of ft
c     d2fdt2 is the second time derivative of ft
c
      if (idef .eq. 0)  then
         return
      else if (idef .eq. 1)  then
         ft     = sin(rfreqt*time)
         dfdt   = rfreqt*cos(rfreqt*time)
         d2fdt2 = -(rfreqt)**2*sin(rfreqt*time)
      else if (idef .eq. 999)  then
         return
      end if
c
      xold = xorg
      yold = yorg
      zold = zorg
c
      xnew = utran*ft + xorg0
      ynew = vtran*ft + yorg0
      znew = wtran*ft + zorg0
c
      dx   = xnew - xold
      dy   = ynew - yold
      dz   = znew - zold
c
      if (isurf .eq. 1) then
         do j=jci,jcf
            do k=kci,kcf
               delti(j,k,1,ll) = wki(j,k,ll)*dx + delti(j,k,1,ll)
               delti(j,k,2,ll) = wki(j,k,ll)*dy + delti(j,k,2,ll)
               delti(j,k,3,ll) = wki(j,k,ll)*dz + delti(j,k,3,ll)
               wki(j,k,ll)     = 0.
            end do
         end do
      else if (isurf .eq. 2) then
         do k=kci,kcf
            do i=ici,icf
               deltj(k,i,1,ll) = wkj(k,i,ll)*dx + deltj(k,i,1,ll)
               deltj(k,i,2,ll) = wkj(k,i,ll)*dy + deltj(k,i,2,ll)
               deltj(k,i,3,ll) = wkj(k,i,ll)*dz + deltj(k,i,3,ll)
               wkj(k,i,ll)     = 0.
            end do
         end do
      else if (isurf .eq. 3) then
         do j=jci,jcf
            do i=ici,icf
               deltk(j,i,1,ll) = wkk(j,i,ll)*dx + deltk(j,i,1,ll)
               deltk(j,i,2,ll) = wkk(j,i,ll)*dy + deltk(j,i,2,ll)
               deltk(j,i,3,ll) = wkk(j,i,ll)*dz + deltk(j,i,3,ll)
               wkk(j,i,ll)     = 0.
            end do
         end do
      end if
c
      xorg = xorg+dx
      yorg = yorg+dy
      zorg = zorg+dz
c
      return
      end
