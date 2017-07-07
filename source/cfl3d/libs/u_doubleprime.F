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
      subroutine u_doubleprime(idim,jdim,kdim,q,ux,vol,si,sj,sk,vx)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute 3-D term U'', defined as:
c     sqrt((d2u/dx2+d2u/dy2+d2u/dz2)**2+
c          (d2v/dx2+d2v/dy2+d2v/dz2)**2+
c          (d2w/dx2+d2w/dy2+d2w/dz2)**2)
c     Put in vx(1)
c     vx(2),vx(3),vx(4) used for temporary storage
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),
     +  vol(jdim,kdim,idim-1),si(jdim,kdim,idim,5),
     +  sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     +  vx(0:jdim,0:kdim,idim-1,4)
      dimension ux(jdim-1,kdim-1,idim-1,9)
c
c     initialize
      do m=1,4
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              vx(j,k,i,m)=0.
            enddo
          enddo
        enddo
      enddo
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
c               tc=0.5*(sj(j  ,k  ,i  ,5)*sj(j  ,k  ,i  ,4)+
c    +                  sj(j+1,k  ,i  ,5)*sj(j+1,k  ,i  ,4))/
c    +                  vol(j,k,i)
                if (j .ge. 2 .and. j .le. jdim-2) then
                  fac=2.
                  dudxp = ux(j+1,k  ,i  ,1)
                  dudyp = ux(j+1,k  ,i  ,2)
                  dudzp = ux(j+1,k  ,i  ,3)
                  dvdxp = ux(j+1,k  ,i  ,4)
                  dvdyp = ux(j+1,k  ,i  ,5)
                  dvdzp = ux(j+1,k  ,i  ,6)
                  dwdxp = ux(j+1,k  ,i  ,7)
                  dwdyp = ux(j+1,k  ,i  ,8)
                  dwdzp = ux(j+1,k  ,i  ,9)
                  dudxm = ux(j-1,k  ,i  ,1)
                  dudym = ux(j-1,k  ,i  ,2)
                  dudzm = ux(j-1,k  ,i  ,3)
                  dvdxm = ux(j-1,k  ,i  ,4)
                  dvdym = ux(j-1,k  ,i  ,5)
                  dvdzm = ux(j-1,k  ,i  ,6)
                  dwdxm = ux(j-1,k  ,i  ,7)
                  dwdym = ux(j-1,k  ,i  ,8)
                  dwdzm = ux(j-1,k  ,i  ,9)
                else if (j .eq. 1) then
                  fac=1.
                  dudxp = ux(j+1,k  ,i  ,1)
                  dudyp = ux(j+1,k  ,i  ,2)
                  dudzp = ux(j+1,k  ,i  ,3)
                  dvdxp = ux(j+1,k  ,i  ,4)
                  dvdyp = ux(j+1,k  ,i  ,5)
                  dvdzp = ux(j+1,k  ,i  ,6)
                  dwdxp = ux(j+1,k  ,i  ,7)
                  dwdyp = ux(j+1,k  ,i  ,8)
                  dwdzp = ux(j+1,k  ,i  ,9)
                  dudxm = ux(j  ,k  ,i  ,1)
                  dudym = ux(j  ,k  ,i  ,2)
                  dudzm = ux(j  ,k  ,i  ,3)
                  dvdxm = ux(j  ,k  ,i  ,4)
                  dvdym = ux(j  ,k  ,i  ,5)
                  dvdzm = ux(j  ,k  ,i  ,6)
                  dwdxm = ux(j  ,k  ,i  ,7)
                  dwdym = ux(j  ,k  ,i  ,8)
                  dwdzm = ux(j  ,k  ,i  ,9)
                else if (j .eq. jdim-1) then
                  fac=1.
                  dudxp = ux(j  ,k  ,i  ,1)
                  dudyp = ux(j  ,k  ,i  ,2)
                  dudzp = ux(j  ,k  ,i  ,3)
                  dvdxp = ux(j  ,k  ,i  ,4)
                  dvdyp = ux(j  ,k  ,i  ,5)
                  dvdzp = ux(j  ,k  ,i  ,6)
                  dwdxp = ux(j  ,k  ,i  ,7)
                  dwdyp = ux(j  ,k  ,i  ,8)
                  dwdzp = ux(j  ,k  ,i  ,9)
                  dudxm = ux(j-1,k  ,i  ,1)
                  dudym = ux(j-1,k  ,i  ,2)
                  dudzm = ux(j-1,k  ,i  ,3)
                  dvdxm = ux(j-1,k  ,i  ,4)
                  dvdym = ux(j-1,k  ,i  ,5)
                  dvdzm = ux(j-1,k  ,i  ,6)
                  dwdxm = ux(j-1,k  ,i  ,7)
                  dwdym = ux(j-1,k  ,i  ,8)
                  dwdzm = ux(j-1,k  ,i  ,9)
                end if
                vx(j,k,i,2)=vx(j,k,i,2)+(xc*(dudxp-dudxm)
     +                                 + yc*(dudyp-dudym)
     +                                 + zc*(dudzp-dudzm))/fac
                vx(j,k,i,3)=vx(j,k,i,3)+(xc*(dvdxp-dvdxm)
     +                                 + yc*(dvdyp-dvdym)
     +                                 + zc*(dvdzp-dvdzm))/fac
                vx(j,k,i,4)=vx(j,k,i,4)+(xc*(dwdxp-dwdxm)
     +                                 + yc*(dwdyp-dwdym)
     +                                 + zc*(dwdzp-dwdzm))/fac
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
c               tc=0.5*(sk(j  ,k  ,i  ,5)*sk(j  ,k  ,i  ,4)+
c    +                  sk(j  ,k+1,i  ,5)*sk(j  ,k+1,i  ,4))/
c    +                  vol(j,k,i)
                if (k .ge. 2 .and. k .le. kdim-2) then
                  fac=2.
                  dudxp = ux(j  ,k+1,i  ,1)
                  dudyp = ux(j  ,k+1,i  ,2)
                  dudzp = ux(j  ,k+1,i  ,3)
                  dvdxp = ux(j  ,k+1,i  ,4)
                  dvdyp = ux(j  ,k+1,i  ,5)
                  dvdzp = ux(j  ,k+1,i  ,6)
                  dwdxp = ux(j  ,k+1,i  ,7)
                  dwdyp = ux(j  ,k+1,i  ,8)
                  dwdzp = ux(j  ,k+1,i  ,9)
                  dudxm = ux(j  ,k-1,i  ,1)
                  dudym = ux(j  ,k-1,i  ,2)
                  dudzm = ux(j  ,k-1,i  ,3)
                  dvdxm = ux(j  ,k-1,i  ,4)
                  dvdym = ux(j  ,k-1,i  ,5)
                  dvdzm = ux(j  ,k-1,i  ,6)
                  dwdxm = ux(j  ,k-1,i  ,7)
                  dwdym = ux(j  ,k-1,i  ,8)
                  dwdzm = ux(j  ,k-1,i  ,9)
                else if (k .eq. 1) then
                  fac=1.
                  dudxp = ux(j  ,k+1,i  ,1)
                  dudyp = ux(j  ,k+1,i  ,2)
                  dudzp = ux(j  ,k+1,i  ,3)
                  dvdxp = ux(j  ,k+1,i  ,4)
                  dvdyp = ux(j  ,k+1,i  ,5)
                  dvdzp = ux(j  ,k+1,i  ,6)
                  dwdxp = ux(j  ,k+1,i  ,7)
                  dwdyp = ux(j  ,k+1,i  ,8)
                  dwdzp = ux(j  ,k+1,i  ,9)
                  dudxm = ux(j  ,k  ,i  ,1)
                  dudym = ux(j  ,k  ,i  ,2)
                  dudzm = ux(j  ,k  ,i  ,3)
                  dvdxm = ux(j  ,k  ,i  ,4)
                  dvdym = ux(j  ,k  ,i  ,5)
                  dvdzm = ux(j  ,k  ,i  ,6)
                  dwdxm = ux(j  ,k  ,i  ,7)
                  dwdym = ux(j  ,k  ,i  ,8)
                  dwdzm = ux(j  ,k  ,i  ,9)
                else if (k .eq. kdim-1) then
                  fac=1.
                  dudxp = ux(j  ,k  ,i  ,1)
                  dudyp = ux(j  ,k  ,i  ,2)
                  dudzp = ux(j  ,k  ,i  ,3)
                  dvdxp = ux(j  ,k  ,i  ,4)
                  dvdyp = ux(j  ,k  ,i  ,5)
                  dvdzp = ux(j  ,k  ,i  ,6)
                  dwdxp = ux(j  ,k  ,i  ,7)
                  dwdyp = ux(j  ,k  ,i  ,8)
                  dwdzp = ux(j  ,k  ,i  ,9)
                  dudxm = ux(j  ,k-1,i  ,1)
                  dudym = ux(j  ,k-1,i  ,2)
                  dudzm = ux(j  ,k-1,i  ,3)
                  dvdxm = ux(j  ,k-1,i  ,4)
                  dvdym = ux(j  ,k-1,i  ,5)
                  dvdzm = ux(j  ,k-1,i  ,6)
                  dwdxm = ux(j  ,k-1,i  ,7)
                  dwdym = ux(j  ,k-1,i  ,8)
                  dwdzm = ux(j  ,k-1,i  ,9)
                end if
                vx(j,k,i,2)=vx(j,k,i,2)+(xc*(dudxp-dudxm)
     +                                 + yc*(dudyp-dudym)
     +                                 + zc*(dudzp-dudzm))/fac
                vx(j,k,i,3)=vx(j,k,i,3)+(xc*(dvdxp-dvdxm)
     +                                 + yc*(dvdyp-dvdym)
     +                                 + zc*(dvdzp-dvdzm))/fac
                vx(j,k,i,4)=vx(j,k,i,4)+(xc*(dwdxp-dwdxm)
     +                                 + yc*(dwdyp-dwdym)
     +                                 + zc*(dwdzp-dwdzm))/fac
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
c               tc=0.5*(si(j  ,k  ,i  ,5)*si(j  ,k  ,i  ,4)+
c    +                  si(j  ,k  ,i+1,5)*si(j  ,k  ,i+1,4))/
c    +                  vol(j,k,i)
                if (i .ge. 2 .and. i .le. idim-2) then
                  fac=2.
                  dudxp = ux(j  ,k  ,i+1,1)
                  dudyp = ux(j  ,k  ,i+1,2)
                  dudzp = ux(j  ,k  ,i+1,3)
                  dvdxp = ux(j  ,k  ,i+1,4)
                  dvdyp = ux(j  ,k  ,i+1,5)
                  dvdzp = ux(j  ,k  ,i+1,6)
                  dwdxp = ux(j  ,k  ,i+1,7)
                  dwdyp = ux(j  ,k  ,i+1,8)
                  dwdzp = ux(j  ,k  ,i+1,9)
                  dudxm = ux(j  ,k  ,i-1,1)
                  dudym = ux(j  ,k  ,i-1,2)
                  dudzm = ux(j  ,k  ,i-1,3)
                  dvdxm = ux(j  ,k  ,i-1,4)
                  dvdym = ux(j  ,k  ,i-1,5)
                  dvdzm = ux(j  ,k  ,i-1,6)
                  dwdxm = ux(j  ,k  ,i-1,7)
                  dwdym = ux(j  ,k  ,i-1,8)
                  dwdzm = ux(j  ,k  ,i-1,9)
                else if (i .eq. 1) then
                  fac=1.
                  dudxp = ux(j  ,k  ,i+1,1)
                  dudyp = ux(j  ,k  ,i+1,2)
                  dudzp = ux(j  ,k  ,i+1,3)
                  dvdxp = ux(j  ,k  ,i+1,4)
                  dvdyp = ux(j  ,k  ,i+1,5)
                  dvdzp = ux(j  ,k  ,i+1,6)
                  dwdxp = ux(j  ,k  ,i+1,7)
                  dwdyp = ux(j  ,k  ,i+1,8)
                  dwdzp = ux(j  ,k  ,i+1,9)
                  dudxm = ux(j  ,k  ,i  ,1)
                  dudym = ux(j  ,k  ,i  ,2)
                  dudzm = ux(j  ,k  ,i  ,3)
                  dvdxm = ux(j  ,k  ,i  ,4)
                  dvdym = ux(j  ,k  ,i  ,5)
                  dvdzm = ux(j  ,k  ,i  ,6)
                  dwdxm = ux(j  ,k  ,i  ,7)
                  dwdym = ux(j  ,k  ,i  ,8)
                  dwdzm = ux(j  ,k  ,i  ,9)
                else if (i .eq. idim-1) then
                  fac=1.
                  dudxp = ux(j  ,k  ,i  ,1)
                  dudyp = ux(j  ,k  ,i  ,2)
                  dudzp = ux(j  ,k  ,i  ,3)
                  dvdxp = ux(j  ,k  ,i  ,4)
                  dvdyp = ux(j  ,k  ,i  ,5)
                  dvdzp = ux(j  ,k  ,i  ,6)
                  dwdxp = ux(j  ,k  ,i  ,7)
                  dwdyp = ux(j  ,k  ,i  ,8)
                  dwdzp = ux(j  ,k  ,i  ,9)
                  dudxm = ux(j  ,k  ,i-1,1)
                  dudym = ux(j  ,k  ,i-1,2)
                  dudzm = ux(j  ,k  ,i-1,3)
                  dvdxm = ux(j  ,k  ,i-1,4)
                  dvdym = ux(j  ,k  ,i-1,5)
                  dvdzm = ux(j  ,k  ,i-1,6)
                  dwdxm = ux(j  ,k  ,i-1,7)
                  dwdym = ux(j  ,k  ,i-1,8)
                  dwdzm = ux(j  ,k  ,i-1,9)
                end if
                vx(j,k,i,2)=vx(j,k,i,2)+(xc*(dudxp-dudxm)
     +                                 + yc*(dudyp-dudym)
     +                                 + zc*(dudzp-dudzm))/fac
                vx(j,k,i,3)=vx(j,k,i,3)+(xc*(dvdxp-dvdxm)
     +                                 + yc*(dvdyp-dvdym)
     +                                 + zc*(dvdzp-dvdzm))/fac
                vx(j,k,i,4)=vx(j,k,i,4)+(xc*(dwdxp-dwdxm)
     +                                 + yc*(dwdyp-dwdym)
     +                                 + zc*(dwdzp-dwdzm))/fac
              enddo
            enddo
          enddo
        end if
c
        do k=1,kdim-1
          do j=1,jdim-1
            do i=1,idim-1
              vx(j,k,i,1)=sqrt(vx(j,k,i,2)**2 + vx(j,k,i,3)**2 +
     +                         vx(j,k,i,4)**2)
            enddo
          enddo
        enddo
c
      return
      end
