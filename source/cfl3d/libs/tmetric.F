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
      subroutine tmetric(jdim,kdim,idim,sj,sk,si,x,y,z,t,t1,t2,t3,nbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate time-metric terms for a grid in motion
c
c     si   - spatial and temporal metrics for i-faces
c     sj   - spatial and temporal metrics for j-faces
c     sk   - spatial and temporal metrics for k-faces
c
c     t   is a temp array containing grid velocity at grid points
c     t1 is a temp array for storage of cell face-average grid velocity
c
c     note: in a jdim*kdim*idim grid there are really only
c                   idim*(jdim-1)*(kdim-1) i-faces
c                   jdim*(idim-1)*(kdim-1) j-faces
c                   kdim*(idim-1)*(jdim-1) k-faces
c     however the metric arrays are dimensioned larger than this -
c     for efficient vectorization, loops throughout the code run
c     over these fictitious faces. thus, for safety, the time metrics
c     for these fictitious faces are set to zero.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension sj(jdim*kdim,idim-1,5),sk(jdim*kdim,idim-1,5),
     .          si(jdim*kdim,idim,5)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     .          t(jdim*kdim*idim,3),t1(jdim*kdim*idim,3),
     .          t2(jdim*kdim,6),t3(jdim*kdim,idim,5)
c
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 50 i=1,idim1
c
c**********************************************
c     unsteady sj terms
c**********************************************
c
c     face-average grid speeds
      js    = jdim*kdim*(i-1)+1
      jk1   = js+jdim
      jki1  = jdim*kdim*i+1
      jk1i1 = jki1+jdim
      n     = jdim*kdim1
      do 1002 izz=1,n
      t1(izz,1)    = 0.25*(t(izz+js-1,1)  +t(izz+jk1-1,1)
     .                    +t(izz+jki1-1,1)+t(izz+jk1i1-1,1))
      t1(izz,2)    = 0.25*(t(izz+js-1,2)  +t(izz+jk1-1,2)
     .                    +t(izz+jki1-1,2)+t(izz+jk1i1-1,2))
      t1(izz,3)    = 0.25*(t(izz+js-1,3)  +t(izz+jk1-1,3)
     .                    +t(izz+jki1-1,3)+t(izz+jk1i1-1,3))
c
      sj(izz,i,5) = - t1(izz,1)*sj(izz,i,1) - t1(izz,2)*sj(izz,i,2)
     .              - t1(izz,3)*sj(izz,i,3)
 1002 continue
c
c     zero out sj(j,k=kdim,i,5) 
      do 1003 izz=1,jdim
      sj(izz+n,i,5) = 0.e0
 1003 continue
c
c**********************************************
c     unsteady sk terms 
c**********************************************
c
c     face-average grid speeds
      j1k   = js+1
      j1ki1 = jki1+1
      n     = jdim*kdim-1
      do 1004 izz=1,n
      t1(izz,1) = 0.25*(t(izz+js-1,1)  +t(izz+j1k-1,1)
     .                 +t(izz+jki1-1,1)+t(izz+j1ki1-1,1))
      t1(izz,2) = 0.25*(t(izz+js-1,2)  +t(izz+j1k-1,2)
     .                 +t(izz+jki1-1,2)+t(izz+j1ki1-1,2))
      t1(izz,3) = 0.25*(t(izz+js-1,3)  +t(izz+j1k-1,3)
     .                 +t(izz+jki1-1,3)+t(izz+j1ki1-1,3))
 1004 continue
c
      n       = jdim*kdim
      do 1005 izz=1,n
      sk(izz,i,5) = - t1(izz,1)*sk(izz,i,1) - t1(izz,2)*sk(izz,i,2)
     .              - t1(izz,3)*sk(izz,i,3)
 1005 continue
c
c     zero out sk(j=jdim,k,i,5)
      do 70 k=1,kdim
      jk      = jdim*k
      sk(jk,i,5) = 0.0
   70 continue
c
   50 continue
c
      do 40 i=1,idim
c
c**********************************************
c     unsteady si terms
c**********************************************
c
c     face-average grid speeds
      js   = jdim*kdim*(i-1)+1
      j1k  = js+1
      jk1  = js+jdim
      j1k1 = jk1+1
      n    = jdim*(kdim-1)-1
      do 1006 izz=1,n
      t1(izz,1) = 0.25*(t(izz+js-1,1) +t(izz+j1k-1,1)
     .                 +t(izz+jk1-1,1)+t(izz+j1k1-1,1))
      t1(izz,2) = 0.25*(t(izz+js-1,2) +t(izz+j1k-1,2)
     .                 +t(izz+jk1-1,2)+t(izz+j1k1-1,2))
      t1(izz,3) = 0.25*(t(izz+js-1,3) +t(izz+j1k-1,3)
     .                 +t(izz+jk1-1,3)+t(izz+j1k1-1,3))
 1006 continue
c
      n       = jdim*kdim1
      do 1007 izz=1,n
      si(izz,i,5) = - t1(izz,1)*si(izz,i,1) - t1(izz,2)*si(izz,i,2)
     .              - t1(izz,3)*si(izz,i,3)
 1007 continue
c
c     zero out si(j,k=kdim,i,5) 
      do 1008 izz=1,jdim
      si(izz+n,i,5) = 0.e0
 1008 continue
c
c     zero out si(j=jdim,k,i,5) 
      do 1009 k=1,kdim1
      jk = jdim*k
      si(jk,i,5) = 0.e0
 1009 continue
c
   40 continue
c
      return
      end
