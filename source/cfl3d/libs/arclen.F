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
      subroutine arclen(idim,jdim,kdim,arci,arcj,arck,x,y,z,nbl,
     .                  nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose: compute normalized arc lengths for use in TFI blending
c     functions.
c
c     If the total arc length of a grid line is < tol, (presumably
c     because it wraps around a singular axis), then replace the arc
c     length in physical space by the length in computational space
c     (i.e. the grid index). This will effective replace the arc-length
c     based blending function with a linear blending function based on
c     grid index. 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension arci(jdim,kdim,idim),arcj(jdim,kdim,idim),
     .          arck(jdim,kdim,idim)
c
      common /zero/ iexp
      common /sklton/ isklton
c
c     tolerance for switch to linear blending function
c     (10.**(-iexp) is machine zero)
c
      tol = max(1.e-07,10.**(-iexp+1))
c
c     arc length measured from i=1 surface
c
      do j=1,jdim
         do k=1,kdim
            arci(j,k,1) = 0.0
         end do
      end do
      nlinear = 0
      do j=1,jdim
         do k=1,kdim
            do i=2,idim
               ds = sqrt((x(j,k,i)-x(j,k,i-1))**2 + 
     .                   (y(j,k,i)-y(j,k,i-1))**2 + 
     .                   (z(j,k,i)-z(j,k,i-1))**2)
               arci(j,k,i) = arci(j,k,i-1) + ds
            end do
            if (real(arci(j,k,idim)) .lt. real(tol)) then
               nlinear = nlinear + 1
               do i=2,idim
                  arci(j,k,i) = float(i-1)
               end do
            end if
         end do
      end do
      if (nlinear.gt.0 .and. isklton.eq.1) then
c           nou(1) = min(nou(1)+1,ibufdim)
c           write(bou(nou(1),1),'('' linear blending functions''
c    .      '' in i-direction, block '',i3)') nlinear,nbl
      end if
c
c     arc length measured from j=1 surface
c
      do i=1,idim
         do k=1,kdim
            arcj(1,k,i) = 0.0
         end do
      end do
      nlinear = 0
      do i=1,idim
         do k=1,kdim
            do j=2,jdim
               ds = sqrt((x(j,k,i)-x(j-1,k,i))**2 + 
     .                   (y(j,k,i)-y(j-1,k,i))**2 + 
     .                   (z(j,k,i)-z(j-1,k,i))**2)
               arcj(j,k,i) = arcj(j-1,k,i) + ds
            end do
            if (real(arcj(jdim,k,i)) .lt. real(tol)) then
               nlinear = nlinear + 1
               do j=2,jdim
                  arcj(j,k,i) = float(j-1)
               end do
            end if
         end do
      end do
      if (nlinear.gt.0 .and. isklton.eq.1) then
c           nou(1) = min(nou(1)+1,ibufdim)
c           write(bou(nou(1),1),'('' linear blending functions''
c    .      '' in j-direction, block '',i3)') nlinear,nbl
      end if
c
c     arc length measured from k=1 surface
c
      do j=1,jdim
         do i=1,idim
            arck(j,1,i) = 0.0
         end do
      end do
      nlinear = 0
      do j=1,jdim
         do i=1,idim
            do k=2,kdim
               ds = sqrt((x(j,k,i)-x(j,k-1,i))**2 + 
     .                   (y(j,k,i)-y(j,k-1,i))**2 + 
     .                   (z(j,k,i)-z(j,k-1,i))**2)
               arck(j,k,i) = arck(j,k-1,i) + ds
            end do
            if (real(arck(j,kdim,i)) .lt. real(tol)) then
               nlinear = nlinear + 1
               do k=2,kdim
                  arck(j,k,i) = float(k-1)
               end do
            end if
         end do
      end do
      if (nlinear.gt.0 .and. isklton.eq.1) then
c           nou(1) = min(nou(1)+1,ibufdim)
c           write(bou(nou(1),1),'('' linear blending functions''
c    .      '' in k-direction, block '',i3)') nlinear,nbl
      end if
c
c     normalize arc lengths
c 
      do i=1,idim
         do j=1,jdim
            do k=1,kdim
               arci(j,k,i) = arci(j,k,i)/arci(j,k,idim)
               arcj(j,k,i) = arcj(j,k,i)/arcj(jdim,k,i)
               arck(j,k,i) = arck(j,k,i)/arck(j,kdim,i)
            end do
         end do
      end do
c
      return
      end
