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
      subroutine sijrate2d(idim,jdim,kdim,q,qj0,qk0,
     . bcj,bck,vol,sj,sk,vx)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute 2-D material derivatives DS11/Dt and DS13/Dt,
c     ignoring the time terms and any contributions from the i-index
c     direction (x-z plane only is assumed).
c     Output of this routine is vx(3)=DS11/Dt and vx(4)=DS13/Dt.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),qj0(kdim,idim-1,5,4),
     +  qk0(jdim,idim-1,5,4),bcj(kdim,idim-1,2),bck(jdim,idim-1,2),
     +  vol(jdim,kdim,idim-1),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),vx(0:jdim,0:kdim,idim-1,4)
c
c   Put q (u and w) velocities in vx(1) and (2), and fill edges with
c   appropriate BC values:
        do i=1,idim-1
        do j=1,jdim-1
        do k=1,kdim-1
          vx(j,k,i,1)=q(j,k,i,2)
          vx(j,k,i,2)=q(j,k,i,4)
        enddo
        enddo
c     Get ghost cell values
        do j=1,jdim-1
          vx(j,0,i,1)=qk0(j,i,2,1)*(1.-bck(j,i,1))+
     +      (2.*qk0(j,i,2,1)-q(j,1,i,2))*bck(j,i,1)
          vx(j,0,i,2)=qk0(j,i,4,1)*(1.-bck(j,i,1))+
     +      (2.*qk0(j,i,4,1)-q(j,1,i,4))*bck(j,i,1)
          vx(j,kdim,i,1)=qk0(j,i,2,3)*(1.-bck(j,i,2))+
     +      (2.*qk0(j,i,2,3)-q(j,kdim-1,i,2))*bck(j,i,2)
          vx(j,kdim,i,2)=qk0(j,i,4,3)*(1.-bck(j,i,2))+
     +      (2.*qk0(j,i,4,3)-q(j,kdim-1,i,4))*bck(j,i,2)
        enddo
        do k=1,kdim-1
          vx(0,k,i,1)=qj0(k,i,2,1)*(1.-bcj(k,i,1))+
     +      (2.*qj0(k,i,2,1)-q(1,k,i,2))*bcj(k,i,1)
          vx(0,k,i,2)=qj0(k,i,4,1)*(1.-bcj(k,i,1))+
     +      (2.*qj0(k,i,4,1)-q(1,k,i,4))*bcj(k,i,1)
          vx(jdim,k,i,1)=qj0(k,i,2,3)*(1.-bcj(k,i,2))+
     +      (2.*qj0(k,i,2,3)-q(jdim-1,k,i,2))*bcj(k,i,2)
          vx(jdim,k,i,2)=qj0(k,i,4,3)*(1.-bcj(k,i,2))+
     +      (2.*qj0(k,i,4,3)-q(jdim-1,k,i,4))*bcj(k,i,2)
        enddo
        enddo
c   compute material derivative of principal strain directions (ignore time term)
c   the following assumes that the boundary terms are at ghost cells
c     do for each component (2-D for now only):
c     j-direction:
        do i=1,idim-1
          do j=1,jdim-1
            if (j .eq. 1) then
              jm=j
            else
              jm=j-1
            end if
            if (j .eq. jdim-1) then
              jp=j
            else
              jp=j+1
            end if
            do k=1,kdim-1
              voll=vol(jm,k,i)
              volu=vol(jp,k,i)
              skm1u=sk(jm,k+1,i,1)
              skm3u=sk(jm,k+1,i,3)
              skm4u=sk(jm,k+1,i,4)
              skm1x=sk(jm,k,i,1)
              skm3x=sk(jm,k,i,3)
              skm4x=sk(jm,k,i,4)
              skp1u=sk(jp,k+1,i,1)
              skp3u=sk(jp,k+1,i,3)
              skp4u=sk(jp,k+1,i,4)
              skp1x=sk(jp,k,i,1)
              skp3x=sk(jp,k,i,3)
              skp4x=sk(jp,k,i,4)
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volu))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volu))
              xm=sj(j  ,k,i,1)*sj(j  ,k,i,4)/(0.5*(vol(j,k,i)+voll))
              zm=sj(j  ,k,i,3)*sj(j  ,k,i,4)/(0.5*(vol(j,k,i)+voll))
              uu=xc*q(j,k,i,2)+zc*q(j,k,i,4)+tc
              d2udn2=xp*(vx(j+1,k,i,1)-vx(j  ,k,i,1)) -
     +               xm*(vx(j  ,k,i,1)-vx(j-1,k,i,1))
              d2wdn2=zp*(vx(j+1,k,i,2)-vx(j  ,k,i,2)) -
     +               zm*(vx(j  ,k,i,2)-vx(j-1,k,i,2))
              d2udn3=zp*(vx(j+1,k,i,1)-vx(j  ,k,i,1)) -
     +               zm*(vx(j  ,k,i,1)-vx(j-1,k,i,1))
              d2wdn3=xp*(vx(j+1,k,i,2)-vx(j  ,k,i,2)) -
     +               xm*(vx(j  ,k,i,2)-vx(j-1,k,i,2))
              vx(j,k,i,3)=0.5*uu*(d2udn2-d2wdn2)
              vx(j,k,i,4)=0.5*uu*(d2udn3+d2wdn3)
c  x-deriv terms:
              xp=0.25*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                skp1u*skp4u+sk(j,k,i,1)*sk(j,k,i,4)+
     +                skp1x*skp4x)/(0.5*(vol(j,k,i)+volu))
              zp=0.25*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                skp3u*skp4u+sk(j,k,i,3)*sk(j,k,i,4)+
     +                skp3x*skp4x)/(0.5*(vol(j,k,i)+volu))
              xm=0.25*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                skm1u*skm4u+sk(j,k,i,1)*sk(j,k,i,4)+
     +                skm1x*skm4x)/(0.5*(vol(j,k,i)+voll))
              zm=0.25*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                skm3u*skm4u+sk(j,k,i,3)*sk(j,k,i,4)+
     +                skm3x*skm4x)/(0.5*(vol(j,k,i)+voll))
              dxunx=0.25*(xp*(vx(jp ,k+1,i,1)-vx(jp ,k-1,i,1)-
     +                        vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1))
     +                   -xm*(vx(jm ,k+1,i,1)-vx(jm ,k-1,i,1)-
     +                        vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1)))
              dxwnz=0.25*(zp*(vx(jp ,k+1,i,2)-vx(jp ,k-1,i,2)-
     +                        vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2))
     +                   -zm*(vx(jm ,k+1,i,2)-vx(jm ,k-1,i,2)-
     +                        vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2)))
              dxunz=0.25*(zp*(vx(jp ,k+1,i,1)-vx(jp ,k-1,i,1)-
     +                        vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1))
     +                   -zm*(vx(jm ,k+1,i,1)-vx(jm ,k-1,i,1)-
     +                        vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1)))
              dxwnx=0.25*(xp*(vx(jp ,k+1,i,2)-vx(jp ,k-1,i,2)-
     +                        vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2))
     +                   -xm*(vx(jm ,k+1,i,2)-vx(jm ,k-1,i,2)-
     +                        vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2)))
              vx(j,k,i,3)=vx(j,k,i,3)+0.5*uu*(dxunx-dxwnz)
              vx(j,k,i,4)=vx(j,k,i,4)+0.5*uu*(dxunz+dxwnx)
            enddo
          enddo
        enddo
c     k-direction:
        do i=1,idim-1
          do k=1,kdim-1
            if (k .eq. 1) then
              km=k
            else
              km=k-1
            end if
            if (k .eq. kdim-1) then
              kp=k
            else
              kp=k+1
            end if
            do j=1,jdim-1
              voll=vol(j,km,i)
              volu=vol(j,kp,i)
              sjm1u=sj(j+1,km,i,1)
              sjm3u=sj(j+1,km,i,3)
              sjm4u=sj(j+1,km,i,4)
              sjm1x=sj(j,km,i,1)
              sjm3x=sj(j,km,i,3)
              sjm4x=sj(j,km,i,4)
              sjp1u=sj(j+1,kp,i,1)
              sjp3u=sj(j+1,kp,i,3)
              sjp4u=sj(j+1,kp,i,4)
              sjp1x=sj(j,kp,i,1)
              sjp3x=sj(j,kp,i,3)
              sjp4x=sj(j,kp,i,4)
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volu))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volu))
              xm=sk(j,k  ,i,1)*sk(j,k  ,i,4)/(0.5*(vol(j,k,i)+voll))
              zm=sk(j,k  ,i,3)*sk(j,k  ,i,4)/(0.5*(vol(j,k,i)+voll))
              uu=xc*q(j,k,i,2)+zc*q(j,k,i,4)+tc
              d2udn2=xp*(vx(j,k+1,i,1)-vx(j,k  ,i,1)) -
     +               xm*(vx(j,k  ,i,1)-vx(j,k-1,i,1))
              d2wdn2=zp*(vx(j,k+1,i,2)-vx(j,k  ,i,2)) -
     +               zm*(vx(j,k  ,i,2)-vx(j,k-1,i,2))
              d2udn3=zp*(vx(j,k+1,i,1)-vx(j,k  ,i,1)) -
     +               zm*(vx(j,k  ,i,1)-vx(j,k-1,i,1))
              d2wdn3=xp*(vx(j,k+1,i,2)-vx(j,k  ,i,2)) -
     +               xm*(vx(j,k  ,i,2)-vx(j,k-1,i,2))
              vx(j,k,i,3)=vx(j,k,i,3)+0.5*uu*(d2udn2-d2wdn2)
              vx(j,k,i,4)=vx(j,k,i,4)+0.5*uu*(d2udn3+d2wdn3)
c  x-deriv terms:
              xp=0.25*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sjp1u*sjp4u+sj(j,k,i,1)*sj(j,k,i,4)+
     +                sjp1x*sjp4x)/(0.5*(vol(j,k,i)+volu))
              zp=0.25*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sjp3u*sjp4u+sj(j,k,i,3)*sj(j,k,i,4)+
     +                sjp3x*sjp4x)/(0.5*(vol(j,k,i)+volu))
              xm=0.25*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sjm1u*sjm4u+sj(j,k,i,1)*sj(j,k,i,4)+
     +                sjm1x*sjm4x)/(0.5*(vol(j,k,i)+voll))
              zm=0.25*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sjm3u*sjm4u+sj(j,k,i,3)*sj(j,k,i,4)+
     +                sjm3x*sjm4x)/(0.5*(vol(j,k,i)+voll))
              dxunx=0.25*(xp*(vx(j+1,kp ,i,1)-vx(j-1,kp ,i,1)-
     +                        vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1))
     +                   -xm*(vx(j+1,km ,i,1)-vx(j-1,km ,i,1)-
     +                        vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1)))
              dxwnz=0.25*(zp*(vx(j+1,kp ,i,2)-vx(j-1,kp ,i,2)-
     +                        vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2))
     +                   -zm*(vx(j+1,km ,i,2)-vx(j-1,km ,i,2)-
     +                        vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2)))
              dxunz=0.25*(zp*(vx(j+1,kp ,i,1)-vx(j-1,kp ,i,1)-
     +                        vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1))
     +                   -zm*(vx(j+1,km ,i,1)-vx(j-1,km ,i,1)-
     +                        vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1)))
              dxwnx=0.25*(xp*(vx(j+1,kp ,i,2)-vx(j-1,kp ,i,2)-
     +                        vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2))
     +                   -xm*(vx(j+1,km ,i,2)-vx(j-1,km ,i,2)-
     +                        vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2)))
              vx(j,k,i,3)=vx(j,k,i,3)+0.5*uu*(dxunx-dxwnz)
              vx(j,k,i,4)=vx(j,k,i,4)+0.5*uu*(dxunz+dxwnx)
            enddo
          enddo
        enddo
c
      return
      end
