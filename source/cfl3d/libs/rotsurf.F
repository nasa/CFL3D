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
      subroutine rotsurf(jdim,kdim,idim,x,y,z,deltj,deltk,delti,nbl,
     .                   idef,xorg,yorg,zorg,omegx,omegy,omegz,thetax,
     .                   thetay,thetaz,rfreqr,ici,icf,jci,jcf,kci,kcf,
     .                   time,nou,bou,nbuf,ibufdim,myid,wkj,wkk,wki)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determines increment to delta displacement due to 
c     surface rotation. 
c
c     idefrm...modulation for mesh deformation
c              = 0 no deformation
c              = 1 sinusoidal variation of surface translation
c              = 2 sinusoidal variation of surface rotation
c              = 999 block undergoes deformation, but not by surface
c                rotation or translation (not handled by this routine)
c
c     surface rotation set in the range (ici,icf), (jci,jcf), (kci,kcf)
c     one pair of the indicies must be identical; this set of constant
c     indicies determines which surface in the grid is deformed 
c
c     deltj/k/i...arrays for storage of delta displacements due to
c                 surface rotation; upon entering this routine, 
c                 these will contain any deltas from surface translation
c
c     rotations are taken with positive angular displacement following
c     the right hand rule
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
c     determine which block face is to be rotated
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
c     ft modulates the rotation
c     dfdt is the time derivative of ft
c     d2fdt2 is the second time derivative of ft
c
      if (idef .eq. 0)  then
         return
      else if (idef .eq. 2)  then
         ft     = sin(rfreqr*time)
         dfdt   = rfreqr*cos(rfreqr*time)
         d2fdt2 = -(rfreqr)**2*sin(rfreqr*time)
      else if (idef .eq. 999)  then
         return
      end if
c
      if (abs(real(omegx)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the x-axis
c***********************************************************************
c
         if (time .ne. 0.) then
            theold = thetax
         else
            theold = 0.e0
         end if
c
c***************************************************
c        calculate rotated y and z surface points
c        delt(1)=0 (unaltered) delt(2)=dy delt(3)=dz 
c***************************************************
c
         theta    = omegx*ft
         dthedt   = omegx*dfdt
         d2thedt2 = omegx*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
c
         if (isurf .eq. 1) then
            do j=jci,jcf
               do k=kci,kcf
                  tempy           = y(j,k,ii)
                  tempz           = z(j,k,ii)
                  delti(j,k,2,ll) = wki(j,k,ll)*((tempy-yorg)*ca
     .                            - (tempz-zorg)*sa+yorg - tempy)
     .                            + delti(j,k,2,ll)
                  delti(j,k,3,ll) = wki(j,k,ll)*((tempy-yorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + delti(j,k,3,ll)
                  wki(j,k,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 2) then
            do k=kci,kcf
               do i=ici,icf
                  tempy           = y(jj,k,i)
                  tempz           = z(jj,k,i)
                  deltj(k,i,2,ll) = wkj(k,i,ll)*((tempy-yorg)*ca
     .                            - (tempz-zorg)*sa+yorg - tempy)
     .                            + deltj(k,i,2,ll)
                  deltj(k,i,3,ll) = wkj(k,i,ll)*((tempy-yorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + deltj(k,i,3,ll)
                  wkj(k,i,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 3) then
            do j=jci,jcf
               do i=ici,icf
                  tempy           = y(j,kk,i)
                  tempz           = z(j,kk,i)
                  deltk(j,i,2,ll) = wkk(j,i,ll)*((tempy-yorg)*ca
     .                            - (tempz-zorg)*sa+yorg - tempy)
     .                            + deltk(j,i,2,ll)
                  deltk(j,i,3,ll) = wkk(j,i,ll)*((tempy-yorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + deltk(j,i,3,ll)
                  wkk(j,i,ll)     = 0.
               end do
            end do
         end if
c
         thetax = theta
c
      else if (abs(real(omegy)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the y-axis
c***********************************************************************
c
         if (time .ne. 0.) then
            theold = thetay
         else
            theold = 0.e0
         end if
c
c***************************************************
c        calculate rotated x and z surface points
c        delt(1)=dx delt(2)=0 (unaltered) delt(3)=dz               
c***************************************************
c
         theta    = omegy*ft
         dthedt   = omegy*dfdt
         d2thedt2 = omegy*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
         if (isurf .eq. 1) then
            do j=jci,jcf
               do k=kci,kcf
                  tempx           = x(j,k,ii)
                  tempz           = z(j,k,ii)
                  delti(j,k,1,ll) = wki(j,k,ll)*((tempx-xorg)*ca
     .                            + (tempz-zorg)*sa+xorg - tempx)
     .                            + delti(j,k,1,ll)
                  delti(j,k,3,ll) = wki(j,k,ll)*(-(tempx-xorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + delti(j,k,3,ll)
                  wki(j,k,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 2) then
            do k=kci,kcf
               do i=ici,icf
                  tempx           = x(jj,k,i)
                  tempz           = z(jj,k,i)
                  deltj(k,i,1,ll) = wkj(k,i,ll)*((tempx-xorg)*ca
     .                            + (tempz-zorg)*sa+xorg - tempx)
     .                            + deltj(k,i,1,ll)
                  deltj(k,i,3,ll) = wkj(k,i,ll)*(-(tempx-xorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + deltj(k,i,3,ll)
                  wkj(k,i,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 3) then
            do j=jci,jcf
               do i=ici,icf
                  tempx           = x(j,kk,i)
                  tempz           = z(j,kk,i)
                  deltk(j,i,1,ll) = wkk(j,i,ll)*((tempx-xorg)*ca
     .                            + (tempz-zorg)*sa+xorg - tempx)
     .                            + deltk(j,i,1,ll)
                  deltk(j,i,3,ll) = wkk(j,i,ll)*(-(tempx-xorg)*sa
     .                            + (tempz-zorg)*ca+zorg - tempz)
     .                            + deltk(j,i,3,ll)
                  wkk(j,i,ll)     = 0.
               end do
            end do
         end if
c
         thetay = theta
c
      else if (abs(real(omegz)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the z-axis
c***********************************************************************
c
         if (time .ne. 0.) then
            theold = thetaz
         else
            theold = 0.e0
         end if
c
c***************************************************
c        calculate rotated x and y surface points
c        delt(1)=dx delt(2)=dy delt(3)=0 (unaltered)
c***************************************************
c
         theta    = omegz*ft
         dthedt   = omegz*dfdt
         d2thedt2 = omegz*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
         if (isurf .eq. 1) then
            do j=jci,jcf
               do k=kci,kcf
                  tempx           = x(j,k,ii)
                  tempy           = y(j,k,ii)
                  delti(j,k,1,ll) = wki(j,k,ll)*((tempx-xorg)*ca
     .                            - (tempy-yorg)*sa+xorg - tempx)
     .                            + delti(j,k,1,ll)
                  delti(j,k,2,ll) = wki(j,k,ll)*((tempx-xorg)*sa
     .                            + (tempy-yorg)*ca+yorg - tempy)
     .                            + delti(j,k,2,ll)
                  wki(j,k,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 2) then
            do k=kci,kcf
               do i=ici,icf
                  tempx           = x(jj,k,i)
                  tempy           = y(jj,k,i)
                  deltj(k,i,1,ll) = wkj(k,i,ll)*((tempx-xorg)*ca
     .                            - (tempy-yorg)*sa+xorg - tempx)
     .                            + deltj(k,i,1,ll)
                  deltj(k,i,2,ll) = wkj(k,i,ll)*((tempx-xorg)*sa
     .                            + (tempy-yorg)*ca+yorg - tempy)
     .                            + deltj(k,i,2,ll)
                  wkj(k,i,ll)     = 0.
               end do
            end do
         else if (isurf .eq. 3) then
            do j=jci,jcf
               do i=ici,icf
                  tempx           = x(j,kk,i)
                  tempy           = y(j,kk,i)
                  deltk(j,i,1,ll) = wkk(j,i,ll)*((tempx-xorg)*ca
     .                            - (tempy-yorg)*sa+xorg - tempx)
     .                            + deltk(j,i,1,ll)
                  deltk(j,i,2,ll) = wkk(j,i,ll)*((tempx-xorg)*sa
     .                            + (tempy-yorg)*ca+yorg - tempy)
     .                            + deltk(j,i,2,ll)
                  wkk(j,i,ll)     = 0.
               end do
            end do
         end if
c
         thetaz = theta
c
      else 
c
         if (isklton .gt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),101)
         end if
 101     format(40h WARNING: this block has zero rotational,
     .          21h surface displacement)
      end if
c
      return
      end
