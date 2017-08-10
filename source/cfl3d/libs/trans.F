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
      subroutine trans(jdim,kdim,idim,t,tti,ttj,ttk,x,y,z,itran,rfreqt,
     .                 utran,vtran,wtran,xorg,yorg,zorg,xold,yold,
     .                 zold,xorg0,yorg0,zorg0,iupdat,time2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determines increment to grid position due to translation
c
c     itran....modulation for translational motion
c              =  0 no translation
c              =  1 constant velocity
c              =  2 sinusoidal variation of displacement
c              =  3 smooth increase in displacement, asypmtotically
c                   reaching a fixed displacement
c              = 99 modulation driven external to this routine: on
c                   input xorg,yorg,zorg already contain new 
c                   translational displacements, and utran,vtran,wtran
c                   already contain new rotational rates.
c
c     iupdat..flag to update grid position
c             = 0 don't update position
c             > 0 update position
c
c     t   is a temp array for storage of grid point velocities 
c     tti is a temp array for storage of i-boundary point accelerations
c     ttj is a temp array for storage of j-boundary point accelerations
c     ttk is a temp array for storage of k-boundary point accelerations
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim),
     .          t(jdim*kdim*idim,3),tti(jdim*kdim,3,2),
     .          ttj(kdim*idim,3,2),ttk(jdim*idim,3,2)
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
c     ft modulates the displacement
c     dfdt is the time derivative of ft
c     d2fdt2 is the second time derivative of ft
c
      if (itran .eq. 0)  then
         return
      else if (itran .eq. 1)  then
         ft     = time2
         dfdt   = 1.0
         d2fdt2 = 0.
      else if (itran .eq. 2)  then
         ft     = sin(rfreqt*time2)
         dfdt   = rfreqt*cos(rfreqt*time2)
         d2fdt2 = -(rfreqt)**2*sin(rfreqt*time2)
      else if (itran .eq. 3)  then
         expt   = exp(-rfreqt*time2)
         ft     = 1.-expt
         dfdt   = rfreqt*expt
         d2fdt2 = -(rfreqt)**2*expt
      else if (itran .eq. 99)  then
         ft     = 0.
         dfdt   = 1.
         d2fdt2 = 0.
      end if
c
c     xnew = utran*ft + xorg0
c     ynew = vtran*ft + yorg0
c     znew = wtran*ft + zorg0

      if (itran .ne. 99) then
         xnew = utran*ft + xorg0
         ynew = vtran*ft + yorg0
         znew = wtran*ft + zorg0
      else
         xnew = xorg
         ynew = yorg
         znew = zorg
      end if
c
      dx     = xnew - xold
      dy     = ynew - yold
      dz     = znew - zold
      dxdt   = utran*dfdt
      dydt   = vtran*dfdt
      dzdt   = wtran*dfdt
      d2xdt2 = utran*d2fdt2
      d2ydt2 = vtran*d2fdt2
      d2zdt2 = wtran*d2fdt2
c
c**************************************************
c     calculate increment to speed of grid points
c     due to translation and add to current values
c     t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
c**************************************************
c
      n = jdim*kdim
      do 11 i=1,idim
      js = jdim*kdim*(i-1)+1
      do 1001 izz=1,n
      t(izz+js-1,1) = t(izz+js-1,1) + dxdt
      t(izz+js-1,2) = t(izz+js-1,2) + dydt
      t(izz+js-1,3) = t(izz+js-1,3) + dzdt
 1001 continue
   11 continue
c
c**************************************************
c     calculate increment to acceleration of grid 
c     points on the block boundaries due to trans-
c     lation and add to current values
c     tti(1)=d2x/dt2 tti(2)=d2y/dt2 tti(3)=d2z/dt2
c**************************************************
c
c     i0/idim boundaries
c
      n = jdim*kdim
      i = 1
      do 12 ii=1,2
      do 1002 izz=1,n
c     neglect acceleration for now
c     tti(izz,1,ii) = tti(izz,1,ii) + d2xdt2
c     tti(izz,2,ii) = tti(izz,2,ii) + d2ydt2
c     tti(izz,3,ii) = tti(izz,3,ii) + d2zdt2  
      tti(izz,1,ii) = tti(izz,1,ii) + 0.0
      tti(izz,2,ii) = tti(izz,2,ii) + 0.0
      tti(izz,3,ii) = tti(izz,3,ii) + 0.0
 1002 continue
      i = i + idim - 1
   12 continue   
c
c     j0/jdim boundaries
c
      n = idim*kdim
      j = 1
      do 13 jj=1,2
      do 1003 izz=1,n
c     neglect acceleration for now
c     ttj(izz,1,jj) = ttj(izz,1,jj) + d2xdt2
c     ttj(izz,2,jj) = ttj(izz,2,jj) + d2ydt2
c     ttj(izz,3,jj) = ttj(izz,3,jj) + d2zdt2
      ttj(izz,1,jj) = ttj(izz,1,jj) + 0.
      ttj(izz,2,jj) = ttj(izz,2,jj) + 0.
      ttj(izz,3,jj) = ttj(izz,3,jj) + 0.
 1003 continue
      j = j + jdim - 1
   13 continue
c
c     k0/kdim boundaries
c
      n = jdim*idim
      k = 1
      do 14 kk=1,2
      do 1004 izz=1,n
c     neglect acceleration for now
c     ttk(izz,1,kk) = ttk(izz,1,kk) + d2xdt2
c     ttk(izz,2,kk) = ttk(izz,2,kk) + d2ydt2
c     ttk(izz,3,kk) = ttk(izz,3,kk) + d2zdt2
      ttk(izz,1,kk) = ttk(izz,1,kk) + 0.
      ttk(izz,2,kk) = ttk(izz,2,kk) + 0.
      ttk(izz,3,kk) = ttk(izz,3,kk) + 0.
 1004 continue
      k = k + kdim - 1
   14 continue
c
c*************************************************
c     update grid to new position; update rotation
c     center to new position
c*************************************************
c
      if (iupdat .gt. 0) then
         n=jdim*kdim*idim
         do 50 izz=1,n
         x(izz) = x(izz)+dx
         y(izz) = y(izz)+dy
         z(izz) = z(izz)+dz
   50    continue
c
         xorg = xold+dx
         yorg = yold+dy
         zorg = zold+dz
c
      end if
c
      return
      end
