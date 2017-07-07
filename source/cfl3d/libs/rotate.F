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
      subroutine rotate(jdim,kdim,idim,t,tti,ttj,ttk,t1,x,y,z,nbl,
     .                  irot,rfreqr,omegx,omegy,omegz,xorg,yorg,zorg,
     .                  thetax,thetay,thetaz,thxold,thyold,thzold,
     .                  iupdat,time2,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determines increment to grid position due to rotation
c
c     irot....modulation for rotational motion
c             =  0 no rotation
c             =  1 constant rotation speed
c             =  2 sinusoidal variation of rotational displacement
c             =  3 smooth increase in rotational displacement, 
c                  asypmtotically reaching a fixed rotational displacement
c             = 99 modulation driven external to this routine: on 
c                  input thetax,thetay,thetaz already contain new 
c                  angular displacements, and omegx,omegy,omegz 
c                  already contain new rotational rates.
c
c     iupdat..flag to update grid position
c             = 0 don't update position
c             > 0 update position
c
c     t1  is a temp array for storage of rotated grid points
c     t   is a temp array for storage of grid point velocities
c     tti is a temp array for storage of i-boundary point accelerations
c     ttj is a temp array for storage of j-boundary point accelerations
c     ttk is a temp array for storage of k-boundary point accelerations
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
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     .          t(jdim*kdim*idim,3),tti(jdim*kdim,3,2),
     .          ttj(kdim*idim,3,2),ttk(jdim*idim,3,2),
     .          t1(jdim*kdim*idim,3) 
c
      common /sklton/ isklton
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     ft modulates the rotation
c     dfdt is the time derivative of ft
c     d2fdt2 is the second time derivative of ft
c
      if (irot .eq. 0)  then
         return
      else if (irot .eq. 1)  then
         ft     = time2
         dfdt   = 1.0
         d2fdt2 = 0.
      else if (irot .eq. 2)  then
         ft     = sin(rfreqr*time2)
         dfdt   = rfreqr*cos(rfreqr*time2)
         d2fdt2 = -(rfreqr)**2*sin(rfreqr*time2)
      else if (irot .eq. 3)  then
         expt   = exp(-rfreqr*time2)
         ft     = 1.-expt
         dfdt   = rfreqr*expt
         d2fdt2 = -(rfreqr)**2*expt
      else if (irot .eq. 99)  then
         ft     = 0.
         dfdt   = 1.0
         d2fdt2 = 0.
      end if
c
      if (abs(real(omegx)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the x-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thxold
         else 
            theold = 0.e0
         end if
c
c**************************************************
c        calculate rotated y and z at grid points
c        t1(1)=x (unaltered) t1(2)=y t1(3)=z 
c**************************************************
c
         if (irot .ne. 99) then
            theta = omegx*ft
         else
            theta = thetax
         end if
         dthedt   = omegx*dfdt
         d2thedt2 = omegx*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
         n  = jdim*kdim
         do 10 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 1000 izz=1,n
         t1(izz+js-1,1) =   x(izz,1,i)
         t1(izz+js-1,2) =  (y(izz,1,i)-yorg)*ca
     .                    -(z(izz,1,i)-zorg)*sa+yorg
         t1(izz+js-1,3) =  (y(izz,1,i)-yorg)*sa
     .                    +(z(izz,1,i)-zorg)*ca+zorg
 1000    continue
   10    continue
c
c*************************************************
c     update grid to new position
c*************************************************
c
         if (iupdat .gt. 0) then
            thetax = theta
            do 12 i=1,idim
            js = jdim*kdim*(i-1)+1
            do 2012 izz=1,n
            x(izz,1,i) =  t1(izz+js-1,1)
            y(izz,1,i) =  t1(izz+js-1,2)
            z(izz,1,i) =  t1(izz+js-1,3)
 2012       continue
   12       continue
c
         end if
c
c***************************************************
c        calculate increment to speed of grid points
c        due to rotation and add to current values
c        t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
c***************************************************
c
         do 11 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 1001 izz=1,n
         t(izz+js-1,1)  = t(izz+js-1,1) + 0.e0
         t(izz+js-1,2)  = t(izz+js-1,2) - (z(izz,1,i)-zorg)*dthedt
         t(izz+js-1,3)  = t(izz+js-1,3) + (y(izz,1,i)-yorg)*dthedt
 1001    continue
   11    continue
c
c***************************************************
c        calculate increment to acceleration of grid 
c        points on the block boundaries due to rotation
c        and add to current values
c        tti(1)=d2x/dt2 tti(2)=d2y/dt2 tti(3)=d2z/dt2
c***************************************************
c
c        i0/idim boundaries
c
         n = jdim*kdim
         i = 1
         do 13 ii=1,2
         do 1002 izz=1,n
c        neglect acceleration for now
         tti(izz,1,ii) = tti(izz,1,ii) + 0.e0
         tti(izz,2,ii) = tti(izz,2,ii) + 0.e0
         tti(izz,3,ii) = tti(izz,3,ii) + 0.e0
 1002    continue
         i = i + idim - 1
   13    continue   
c
c        j0/jdim boundaries
c
         n = idim*kdim
         j = 1
         do 14 jj=1,2
         do 1003 izz=1,n
c        neglect acceleration for now
         ttj(izz,1,jj) = ttj(izz,1,jj) + 0.e0
         ttj(izz,2,jj) = ttj(izz,2,jj) + 0.e0
         ttj(izz,3,jj) = ttj(izz,3,jj) + 0.e0
 1003    continue
         j = j + jdim - 1
   14    continue
c
c        k0/kdim boundaries
c
         n = jdim*idim
         k = 1
         do 15 kk=1,2
         do 1004 izz=1,n
c        neglect acceleration for now
         ttk(izz,1,kk) = ttk(izz,1,kk) + 0.e0
         ttk(izz,2,kk) = ttk(izz,2,kk) + 0.e0
         ttk(izz,3,kk) = ttk(izz,3,kk) + 0.e0
 1004    continue
         k = k + kdim - 1
   15    continue
c
      else if (abs(real(omegy)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the y-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thyold

         else
            theold = 0.e0
         end if
c
c**************************************************
c        calculate rotated x and z at grid points
c        t1(1)=x t1(2)=y (unaltered) t1(3)=z               
c**************************************************
c
         if (irot .ne. 99) then
            theta = omegy*ft
         else
            theta = thetay
         end if
         dthedt   = omegy*dfdt
         d2thedt2 = omegy*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
         n  = jdim*kdim
         do 20 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 2000 izz=1,n
         t1(izz+js-1,1) =  (x(izz,1,i)-xorg)*ca
     .                    +(z(izz,1,i)-zorg)*sa+xorg
         t1(izz+js-1,2) =   y(izz,1,i)
         t1(izz+js-1,3) = -(x(izz,1,i)-xorg)*sa
     .                    +(z(izz,1,i)-zorg)*ca+zorg
 2000    continue
   20    continue
c
c*************************************************
c     update grid to new position
c*************************************************
c
         if (iupdat .gt. 0) then
            thetay = theta
            do 22 i=1,idim
            js = jdim*kdim*(i-1)+1
            do 2022 izz=1,n
            x(izz,1,i) =  t1(izz+js-1,1)
            y(izz,1,i) =  t1(izz+js-1,2)
            z(izz,1,i) =  t1(izz+js-1,3)
 2022       continue
   22       continue
c
         end if
c
c***************************************************
c        calculate increment to speed of grid points
c        due to rotation and add to current values
c        t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
c***************************************************
c
         do 21 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 2001 izz=1,n
         t(izz+js-1,1)  = t(izz+js-1,1) + (z(izz,1,i)-zorg)*dthedt
         t(izz+js-1,2)  = t(izz+js-1,2) + 0.e0
         t(izz+js-1,3)  = t(izz+js-1,3) - (x(izz,1,i)-xorg)*dthedt
 2001    continue
   21    continue
c
c***************************************************
c        calculate increment to acceleration of grid 
c        points on the block boundaries due to rotation
c        and add to current values
c        tti(1)=d2x/dt2 tti(2)=d2y/dt2 tti(3)=d2z/dt2
c***************************************************
c
c        i0/idim boundaries
c
         n = jdim*kdim
         i = 1
         do 23 ii=1,2
         do 2002 izz=1,n
c        neglect acceleration for now
         tti(izz,1,ii) = tti(izz,1,ii) + 0.e0
         tti(izz,2,ii) = tti(izz,2,ii) + 0.e0
         tti(izz,3,ii) = tti(izz,3,ii) + 0.e0
 2002    continue
         i = i + idim - 1
   23    continue   
c
c        j0/jdim boundaries
c
         n = idim*kdim
         j = 1
         do 24 jj=1,2
         do 2003 izz=1,n
c        neglect acceleration for now
         ttj(izz,1,jj) = ttj(izz,1,jj) + 0.e0
         ttj(izz,2,jj) = ttj(izz,2,jj) + 0.e0
         ttj(izz,3,jj) = ttj(izz,3,jj) + 0.e0
 2003    continue
         j = j + jdim - 1
   24    continue
c
c        k0/kdim boundaries
c
         n = jdim*idim
         k = 1
         do 25 kk=1,2
         do 2004 izz=1,n
c        neglect acceleration for now
         ttk(izz,1,kk) = ttk(izz,1,kk) + 0.e0
         ttk(izz,2,kk) = ttk(izz,2,kk) + 0.e0
         ttk(izz,3,kk) = ttk(izz,3,kk) + 0.e0
 2004    continue
         k = k + kdim - 1
   25    continue
c
      else if (abs(real(omegz)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the z-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thzold

         else
            theold = 0.e0
         end if
c
c**************************************************
c        calculate rotated x and y at grid points
c        t1(1)=x t1(2)=y t1(3)=z (unaltered)
c**************************************************
c
         if (irot .ne. 99) then
            theta = omegz*ft
         else
            theta = thetaz
         end if
         dthedt   = omegz*dfdt
         d2thedt2 = omegz*d2fdt2
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
         n  = jdim*kdim
         do 30 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 3000 izz=1,n
         t1(izz+js-1,1) =  (x(izz,1,i)-xorg)*ca
     .                    -(y(izz,1,i)-yorg)*sa+xorg
         t1(izz+js-1,2) =  (x(izz,1,i)-xorg)*sa
     .                    +(y(izz,1,i)-yorg)*ca+yorg
         t1(izz+js-1,3) =   z(izz,1,i)
 3000    continue
   30    continue
c
c*************************************************
c     update grid to new position
c*************************************************
c
         if (iupdat .gt. 0) then
            thetaz = theta
            do 32 i=1,idim
            js = jdim*kdim*(i-1)+1
            do 2032 izz=1,n
            x(izz,1,i) =  t1(izz+js-1,1)
            y(izz,1,i) =  t1(izz+js-1,2)
            z(izz,1,i) =  t1(izz+js-1,3)
 2032       continue
   32       continue
c
         end if
c
c***************************************************
c        calculate increment to speed of grid points
c        due to rotation and add to current values
c        t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
c***************************************************
c
         do 31 i=1,idim
         js = jdim*kdim*(i-1)+1
         do 3001 izz=1,n
         t(izz+js-1,1)  = t(izz+js-1,1) - (y(izz,1,i)-yorg)*dthedt
         t(izz+js-1,2)  = t(izz+js-1,2) + (x(izz,1,i)-xorg)*dthedt
         t(izz+js-1,3)  = t(izz+js-1,3) + 0.e0
 3001    continue
   31    continue
c
c***************************************************
c        calculate increment to acceleration of grid 
c        points on the block boundaries due to rotation
c        and add to current values
c        tti(1)=d2x/dt2 tti(2)=d2y/dt2 tti(3)=d2z/dt2
c***************************************************
c
c        i0/idim boundaries
c
         n = jdim*kdim
         i = 1
         do 33 ii=1,2
         do 3002 izz=1,n
c        neglect acceleration for now
         tti(izz,1,ii) = tti(izz,1,ii) + 0.e0
         tti(izz,2,ii) = tti(izz,2,ii) + 0.e0
         tti(izz,3,ii) = tti(izz,3,ii) + 0.e0
 3002    continue
         i = i + idim - 1
   33    continue   
c
c        j0/jdim boundaries
c
         n = idim*kdim
         j = 1
         do 34 jj=1,2
         do 3003 izz=1,n
c        neglect acceleration for now
         ttj(izz,1,jj) = ttj(izz,1,jj) + 0.e0
         ttj(izz,2,jj) = ttj(izz,2,jj) + 0.e0
         ttj(izz,3,jj) = ttj(izz,3,jj) + 0.e0
 3003    continue
         j = j + jdim - 1
   34    continue
c
c        k0/kdim boundaries
c
         n = jdim*idim
         k = 1
         do 35 kk=1,2
         do 3004 izz=1,n
c        neglect acceleration for now
         ttk(izz,1,kk) = ttk(izz,1,kk) + 0.e0
         ttk(izz,2,kk) = ttk(izz,2,kk) + 0.e0
         ttk(izz,3,kk) = ttk(izz,3,kk) + 0.e0
 3004    continue
         k = k + kdim - 1
   35    continue
c
c
      else 
c
         if (isklton .gt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),101)
         end if
 101     format(40h WARNING: this block has zero rotational,
     .          13h displacement)
      end if
c
      return
      end
