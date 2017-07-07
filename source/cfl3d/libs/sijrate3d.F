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
      subroutine sijrate3d(idim,jdim,kdim,q,ux,vol,si,sj,sk,vx)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute 3-D material derivatives DSij/Dt,
c     ignoring the time terms: i.e., uj*d(Sij)/dxj.
c     Unlike, sijrate2d, which finds DSij/Dt directly from velocity
c     vectors (taking 2nd derivatives), this routine works with
c     the ux() velocity derivatives already obtained at cell centers,
c     taking central differences of them.  As a result of this,
c     at boundaries, one-sided differencing is employed 
c     (i.e., unlike sijrate2d, this routine
c     does NOT account for BC information, so the derivatives are
c     lower order at all block boundaries).  If the index in any 
c     direction is 2, such as would occur for 2-D, then the derivative 
c     in that index direction is set to zero.
c       vx(1)=DS11/Dt, vx(2)=DS12/Dt, vx(3)=DS13/Dt,
c       vx(4)=DS22/Dt, vx(5)=DS23/Dt, vx(6)=DS33/Dt
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),
     +  vol(jdim,kdim,idim-1),si(jdim,kdim,idim,5),
     +  sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     +  vx(0:jdim,0:kdim,idim-1,6)
      dimension ux(jdim-1,kdim-1,idim-1,9)
c
c     initialize
      do m=1,6
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              vx(j,k,i,m)=0.
            enddo
          enddo
        enddo
      enddo
c     Note:
c       s11 = ux(j,k,i,1)
c       s22 = ux(j,k,i,5)
c       s33 = ux(j,k,i,9)
c       s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
c       s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
c       s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
c     j-direction:
        if (jdim .gt. 2) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                xc=0.5*(sj(j  ,k  ,i  ,1)*sj(j  ,k  ,i  ,4)+
     +                  sj(j+1,k  ,i  ,1)*sj(j+1,k  ,i  ,4))/
     +                  vol(j,k,i)
                yc=0.5*(sj(j  ,k  ,i  ,2)*sj(j  ,k  ,i  ,4)+
     +                  sj(j+1,k  ,i  ,2)*sj(j+1,k  ,i  ,4))/
     +                  vol(j,k,i)
                zc=0.5*(sj(j  ,k  ,i  ,3)*sj(j  ,k  ,i  ,4)+
     +                  sj(j+1,k  ,i  ,3)*sj(j+1,k  ,i  ,4))/
     +                  vol(j,k,i)
                tc=0.5*(sj(j  ,k  ,i  ,5)*sj(j  ,k  ,i  ,4)+
     +                  sj(j+1,k  ,i  ,5)*sj(j+1,k  ,i  ,4))/
     +                  vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                if (j .ge. 2 .and. j .le. jdim-2) then
                  fac=2.
                  s11p = ux(j+1,k  ,i  ,1)
                  s22p = ux(j+1,k  ,i  ,5)
                  s33p = ux(j+1,k  ,i  ,9)
                  s12p = 0.5*(ux(j+1,k  ,i  ,2) + ux(j+1,k  ,i  ,4))
                  s13p = 0.5*(ux(j+1,k  ,i  ,3) + ux(j+1,k  ,i  ,7))
                  s23p = 0.5*(ux(j+1,k  ,i  ,6) + ux(j+1,k  ,i  ,8))
                  s11m = ux(j-1,k  ,i  ,1)
                  s22m = ux(j-1,k  ,i  ,5)
                  s33m = ux(j-1,k  ,i  ,9)
                  s12m = 0.5*(ux(j-1,k  ,i  ,2) + ux(j-1,k  ,i  ,4))
                  s13m = 0.5*(ux(j-1,k  ,i  ,3) + ux(j-1,k  ,i  ,7))
                  s23m = 0.5*(ux(j-1,k  ,i  ,6) + ux(j-1,k  ,i  ,8))
                else if (j .eq. 1) then
                  fac=1.
                  s11p = ux(j+1,k  ,i  ,1)
                  s22p = ux(j+1,k  ,i  ,5)
                  s33p = ux(j+1,k  ,i  ,9)
                  s12p = 0.5*(ux(j+1,k  ,i  ,2) + ux(j+1,k  ,i  ,4))
                  s13p = 0.5*(ux(j+1,k  ,i  ,3) + ux(j+1,k  ,i  ,7))
                  s23p = 0.5*(ux(j+1,k  ,i  ,6) + ux(j+1,k  ,i  ,8))
                  s11m = ux(j  ,k  ,i  ,1)
                  s22m = ux(j  ,k  ,i  ,5)
                  s33m = ux(j  ,k  ,i  ,9)
                  s12m = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13m = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23m = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                else if (j .eq. jdim-1) then
                  fac=1.
                  s11p = ux(j  ,k  ,i  ,1)
                  s22p = ux(j  ,k  ,i  ,5)
                  s33p = ux(j  ,k  ,i  ,9)
                  s12p = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13p = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23p = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                  s11m = ux(j-1,k  ,i  ,1)
                  s22m = ux(j-1,k  ,i  ,5)
                  s33m = ux(j-1,k  ,i  ,9)
                  s12m = 0.5*(ux(j-1,k  ,i  ,2) + ux(j-1,k  ,i  ,4))
                  s13m = 0.5*(ux(j-1,k  ,i  ,3) + ux(j-1,k  ,i  ,7))
                  s23m = 0.5*(ux(j-1,k  ,i  ,6) + ux(j-1,k  ,i  ,8))
                end if
                vx(j,k,i,1)=vx(j,k,i,1)+uu*(s11p-s11m)/fac
                vx(j,k,i,2)=vx(j,k,i,2)+uu*(s12p-s12m)/fac
                vx(j,k,i,3)=vx(j,k,i,3)+uu*(s13p-s13m)/fac
                vx(j,k,i,4)=vx(j,k,i,4)+uu*(s22p-s22m)/fac
                vx(j,k,i,5)=vx(j,k,i,5)+uu*(s23p-s23m)/fac
                vx(j,k,i,6)=vx(j,k,i,6)+uu*(s33p-s33m)/fac
              enddo
            enddo
          enddo
        end if
c     k-direction:
        if (kdim .gt. 2) then
          do i=1,idim-1
            do j=1,jdim-1
              do k=1,kdim-1
                xc=0.5*(sk(j  ,k  ,i  ,1)*sk(j  ,k  ,i  ,4)+
     +                  sk(j  ,k+1,i  ,1)*sk(j  ,k+1,i  ,4))/
     +                  vol(j,k,i)
                yc=0.5*(sk(j  ,k  ,i  ,2)*sk(j  ,k  ,i  ,4)+
     +                  sk(j  ,k+1,i  ,2)*sk(j  ,k+1,i  ,4))/
     +                  vol(j,k,i)
                zc=0.5*(sk(j  ,k  ,i  ,3)*sk(j  ,k  ,i  ,4)+
     +                  sk(j  ,k+1,i  ,3)*sk(j  ,k+1,i  ,4))/
     +                  vol(j,k,i)
                tc=0.5*(sk(j  ,k  ,i  ,5)*sk(j  ,k  ,i  ,4)+
     +                  sk(j  ,k+1,i  ,5)*sk(j  ,k+1,i  ,4))/
     +                  vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                if (k .ge. 2 .and. k .le. kdim-2) then
                  fac=2.
                  s11p = ux(j  ,k+1,i  ,1)
                  s22p = ux(j  ,k+1,i  ,5)
                  s33p = ux(j  ,k+1,i  ,9)
                  s12p = 0.5*(ux(j  ,k+1,i  ,2) + ux(j  ,k+1,i  ,4))
                  s13p = 0.5*(ux(j  ,k+1,i  ,3) + ux(j  ,k+1,i  ,7))
                  s23p = 0.5*(ux(j  ,k+1,i  ,6) + ux(j  ,k+1,i  ,8))
                  s11m = ux(j  ,k-1,i  ,1)
                  s22m = ux(j  ,k-1,i  ,5)
                  s33m = ux(j  ,k-1,i  ,9)
                  s12m = 0.5*(ux(j  ,k-1,i  ,2) + ux(j  ,k-1,i  ,4))
                  s13m = 0.5*(ux(j  ,k-1,i  ,3) + ux(j  ,k-1,i  ,7))
                  s23m = 0.5*(ux(j  ,k-1,i  ,6) + ux(j  ,k-1,i  ,8))
                else if (k .eq. 1) then
                  fac=1.
                  s11p = ux(j  ,k+1,i  ,1)
                  s22p = ux(j  ,k+1,i  ,5)
                  s33p = ux(j  ,k+1,i  ,9)
                  s12p = 0.5*(ux(j  ,k+1,i  ,2) + ux(j  ,k+1,i  ,4))
                  s13p = 0.5*(ux(j  ,k+1,i  ,3) + ux(j  ,k+1,i  ,7))
                  s23p = 0.5*(ux(j  ,k+1,i  ,6) + ux(j  ,k+1,i  ,8))
                  s11m = ux(j  ,k  ,i  ,1)
                  s22m = ux(j  ,k  ,i  ,5)
                  s33m = ux(j  ,k  ,i  ,9)
                  s12m = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13m = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23m = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                else if (k .eq. kdim-1) then
                  fac=1.
                  s11p = ux(j  ,k  ,i  ,1)
                  s22p = ux(j  ,k  ,i  ,5)
                  s33p = ux(j  ,k  ,i  ,9)
                  s12p = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13p = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23p = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                  s11m = ux(j  ,k-1,i  ,1)
                  s22m = ux(j  ,k-1,i  ,5)
                  s33m = ux(j  ,k-1,i  ,9)
                  s12m = 0.5*(ux(j  ,k-1,i  ,2) + ux(j  ,k-1,i  ,4))
                  s13m = 0.5*(ux(j  ,k-1,i  ,3) + ux(j  ,k-1,i  ,7))
                  s23m = 0.5*(ux(j  ,k-1,i  ,6) + ux(j  ,k-1,i  ,8))
                end if
                vx(j,k,i,1)=vx(j,k,i,1)+uu*(s11p-s11m)/fac
                vx(j,k,i,2)=vx(j,k,i,2)+uu*(s12p-s12m)/fac
                vx(j,k,i,3)=vx(j,k,i,3)+uu*(s13p-s13m)/fac
                vx(j,k,i,4)=vx(j,k,i,4)+uu*(s22p-s22m)/fac
                vx(j,k,i,5)=vx(j,k,i,5)+uu*(s23p-s23m)/fac
                vx(j,k,i,6)=vx(j,k,i,6)+uu*(s33p-s33m)/fac
              enddo
            enddo
          enddo
        end if
c     i-direction:
        if (idim .gt. 2) then
          do k=1,kdim-1
            do j=1,jdim-1
              do i=1,idim-1
                xc=0.5*(si(j  ,k  ,i  ,1)*si(j  ,k  ,i  ,4)+
     +                  si(j  ,k  ,i+1,1)*si(j  ,k  ,i+1,4))/
     +                  vol(j,k,i)
                yc=0.5*(si(j  ,k  ,i  ,2)*si(j  ,k  ,i  ,4)+
     +                  si(j  ,k  ,i+1,2)*si(j  ,k  ,i+1,4))/
     +                  vol(j,k,i)
                zc=0.5*(si(j  ,k  ,i  ,3)*si(j  ,k  ,i  ,4)+
     +                  si(j  ,k  ,i+1,3)*si(j  ,k  ,i+1,4))/
     +                  vol(j,k,i)
                tc=0.5*(si(j  ,k  ,i  ,5)*si(j  ,k  ,i  ,4)+
     +                  si(j  ,k  ,i+1,5)*si(j  ,k  ,i+1,4))/
     +                  vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                if (i .ge. 2 .and. i .le. idim-2) then
                  fac=2.
                  s11p = ux(j  ,k  ,i+1,1)
                  s22p = ux(j  ,k  ,i+1,5)
                  s33p = ux(j  ,k  ,i+1,9)
                  s12p = 0.5*(ux(j  ,k  ,i+1,2) + ux(j  ,k  ,i+1,4))
                  s13p = 0.5*(ux(j  ,k  ,i+1,3) + ux(j  ,k  ,i+1,7))
                  s23p = 0.5*(ux(j  ,k  ,i+1,6) + ux(j  ,k  ,i+1,8))
                  s11m = ux(j  ,k  ,i-1,1)
                  s22m = ux(j  ,k  ,i-1,5)
                  s33m = ux(j  ,k  ,i-1,9)
                  s12m = 0.5*(ux(j  ,k  ,i-1,2) + ux(j  ,k  ,i-1,4))
                  s13m = 0.5*(ux(j  ,k  ,i-1,3) + ux(j  ,k  ,i-1,7))
                  s23m = 0.5*(ux(j  ,k  ,i-1,6) + ux(j  ,k  ,i-1,8))
                else if (i .eq. 1) then
                  fac=1.
                  s11p = ux(j  ,k  ,i+1,1)
                  s22p = ux(j  ,k  ,i+1,5)
                  s33p = ux(j  ,k  ,i+1,9)
                  s12p = 0.5*(ux(j  ,k  ,i+1,2) + ux(j  ,k  ,i+1,4))
                  s13p = 0.5*(ux(j  ,k  ,i+1,3) + ux(j  ,k  ,i+1,7))
                  s23p = 0.5*(ux(j  ,k  ,i+1,6) + ux(j  ,k  ,i+1,8))
                  s11m = ux(j  ,k  ,i  ,1)
                  s22m = ux(j  ,k  ,i  ,5)
                  s33m = ux(j  ,k  ,i  ,9)
                  s12m = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13m = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23m = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                else if (i .eq. idim-1) then
                  fac=1.
                  s11p = ux(j  ,k  ,i  ,1)
                  s22p = ux(j  ,k  ,i  ,5)
                  s33p = ux(j  ,k  ,i  ,9)
                  s12p = 0.5*(ux(j  ,k  ,i  ,2) + ux(j  ,k  ,i  ,4))
                  s13p = 0.5*(ux(j  ,k  ,i  ,3) + ux(j  ,k  ,i  ,7))
                  s23p = 0.5*(ux(j  ,k  ,i  ,6) + ux(j  ,k  ,i  ,8))
                  s11m = ux(j  ,k  ,i-1,1)
                  s22m = ux(j  ,k  ,i-1,5)
                  s33m = ux(j  ,k  ,i-1,9)
                  s12m = 0.5*(ux(j  ,k  ,i-1,2) + ux(j  ,k  ,i-1,4))
                  s13m = 0.5*(ux(j  ,k  ,i-1,3) + ux(j  ,k  ,i-1,7))
                  s23m = 0.5*(ux(j  ,k  ,i-1,6) + ux(j  ,k  ,i-1,8))
                end if
                vx(j,k,i,1)=vx(j,k,i,1)+uu*(s11p-s11m)/fac
                vx(j,k,i,2)=vx(j,k,i,2)+uu*(s12p-s12m)/fac
                vx(j,k,i,3)=vx(j,k,i,3)+uu*(s13p-s13m)/fac
                vx(j,k,i,4)=vx(j,k,i,4)+uu*(s22p-s22m)/fac
                vx(j,k,i,5)=vx(j,k,i,5)+uu*(s23p-s23m)/fac
                vx(j,k,i,6)=vx(j,k,i,6)+uu*(s33p-s33m)/fac
              enddo
            enddo
          enddo
        end if
c
      return
      end
