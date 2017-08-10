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
      subroutine getdelt(maxbl,maxsegdg,idim,jdim,kdim,delti,deltj,
     .                   deltk,x,y,z,nbl,icsi,icsf,jcsi,jcsf,kcsi,
     .                   kcsf,iseg,nou,bou,nbuf,ibufdim,wkj,wkk,wki)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate the displacements between the new surface
c     (stored in deltj/deltk,delti when the routine is entered) and
c     the current one, and store the displacements in the 
c     deltj/deltk,delti arrays.
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
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
      is = icsi(nbl,iseg)
      ie = icsf(nbl,iseg)
      js = jcsi(nbl,iseg)
      je = jcsf(nbl,iseg)
      ks = kcsi(nbl,iseg)
      ke = kcsf(nbl,iseg)
c
      if (is .eq. ie) then
          mm = 1
          if (is .eq. idim) mm = 2
          do k=ks,ke
             do j=js,je
                delti(j,k,1,mm) = delti(j,k,1,mm) 
     .                          - x(j,k,is)*wki(j,k,mm)
                delti(j,k,2,mm) = delti(j,k,2,mm) 
     .                          - y(j,k,is)*wki(j,k,mm)
                delti(j,k,3,mm) = delti(j,k,3,mm) 
     .                          - z(j,k,is)*wki(j,k,mm)
                wki(j,k,mm)     = 0.
             end do
          end do
      else if (js .eq. je) then
          mm = 1
          if (js .eq. jdim) mm = 2
          do i=is,ie
             do k=ks,ke
                deltj(k,i,1,mm) = deltj(k,i,1,mm) 
     .                          - x(js,k,i)*wkj(k,i,mm)
                deltj(k,i,2,mm) = deltj(k,i,2,mm) 
     .                          - y(js,k,i)*wkj(k,i,mm)
                deltj(k,i,3,mm) = deltj(k,i,3,mm) 
     .                          - z(js,k,i)*wkj(k,i,mm)
                wkj(k,i,mm)     = 0.
             end do
          end do
      else if (ks .eq. ke) then
          mm = 1
          if (ks .eq. kdim) mm = 2
          do i=is,ie
             do j=js,je
                deltk(j,i,1,mm) = deltk(j,i,1,mm) 
     .                          - x(j,ks,i)*wkk(j,i,mm)
                deltk(j,i,2,mm) = deltk(j,i,2,mm) 
     .                          - y(j,ks,i)*wkk(j,i,mm)
                deltk(j,i,3,mm) = deltk(j,i,3,mm) 
     .                          - z(j,ks,i)*wkk(j,i,mm)
                wkk(j,i,mm)     = 0.
             end do
          end do
      else
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' error in getdelt...one surface'',
     .    '' dimension must be 1'')')
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      return
      end
