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
      subroutine qface(jdim,kdim,idim,dum,dumj0,dumk0,dumi0,bcj,bck,bci,
     .                 blank,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine variables at cell-face centers on the edges
c     of the grid and install in the j0/k0/i0 arrays for use in output 
c     routines. The input data in the dum,dumj0,dumk0, and dumi0
c     arrays are assumed to have the same form as q,qj0,qk0, and qi0
c     arrays:
c
c     dum(j=1,jdim;k=1,kdim;i=1,idim;l=1,ldim)
c
c     dumj0(k=1,kdim;i=1,idim-1;l=1,ldim;mm=1,4)
c
c     dumk0(j=1,jdim;i=1,idim-1;l=1,ldim;mm=1,4)
c
c     dumi0(j=1,jdim;k=1,kdim;  l=1,ldim;mm=1,4)
c
c     ldim is the dimension covering the number of solution variables:
c
c     ldim=5 for q/qj0/qk0/qi0
c     ldim=1 for vist3d/vj0/vk0/vi0
c     ldim=nummem for tursav/tj0/tk0/ti0
c
c     The bcj/bck/bci arrays indicate whether the data stored in the
c     j0/k0/i0 boundary data arrays in the mm=1 and mm=3 arrays are
c     already face-center based (array entry = 1) or at ghost cell
c     centers (array entry = 1).
c
c     NOTE: The resultant edge values are stored in the mm=2 or
c     mm=4 locations for min and max faces, respectively.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension dum(jdim,kdim,idim,ldim), dumj0(kdim,idim-1,ldim,4),
     .          dumk0(jdim,idim-1,ldim,4),dumi0(jdim,kdim,ldim,4)
      dimension blank(jdim,kdim,idim)
c
      common /twod/ i2d
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 9000 l=1,ldim
c
c     j - boundaries
c
      do 6400 m=1,2
c
      imin = 1
      imax = idim-1
      kmin = 1
      kmax = kdim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         do 6120 ii=imin,imax
         do 6120 kk=kmin,kmax
         aa = 1. +  bcj(kk,ii,m)
         bb = 1. -  bcj(kk,ii,m) 
         cc = 0.
         if (blank(1,kk,ii).eq.0) then
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumj0(kk,ii,l,2) = 0.5*(aa*dumj0(kk,ii,l,1)
     .                    +      bb*dum(1,kk,ii,l)
     .                    +      cc*dum(2,kk,ii,l))
 6120    continue
c
      else
c
c        right boundary
c
         do 6330 ii=imin,imax
         do 6330 kk=kmin,kmax
         aa = 1. +  bcj(kk,ii,m)
         bb = 1. -  bcj(kk,ii,m) 
         cc = 0.
         if (blank(jdim-1,kk,ii).eq.0) then 
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumj0(kk,ii,l,4) = 0.5*(aa*dumj0(kk,ii,l,3)
     .                    +      bb*dum(jdim-1,kk,ii,l)
     .                    +      cc*dum(jdim-2,kk,ii,l))
 6330    continue
c
      end if
 6400 continue
c
c     k - boundaries
c
      do 7400 m=1,2
c
      imin = 1
      imax = idim-1
      jmin = 1
      jmax = jdim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         do 7120 ii=imin,imax
         do 7120 jj=jmin,jmax
         aa = 1. +  bck(jj,ii,m)
         bb = 1. -  bck(jj,ii,m)
         cc = 0.
         if (blank(jj,1,ii).eq.0) then 
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumk0(jj,ii,l,2) = 0.5*(aa*dumk0(jj,ii,l,1)
     .                    +      bb*dum(jj,1,ii,l)
     .                    +      cc*dum(jj,2,ii,l))
 7120    continue
c
      else
c
c        right boundary
c
         do 7330 ii=imin,imax
         do 7330 jj=jmin,jmax
         aa = 1. +  bck(jj,ii,m)
         bb = 1. -  bck(jj,ii,m)
         cc = 0.
         if (blank(jj,kdim-1,ii).eq.0) then 
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumk0(jj,ii,l,4) = 0.5*(aa*dumk0(jj,ii,l,3)
     .                    +      bb*dum(jj,kdim-1,ii,l)
     .                    +      cc*dum(jj,kdim-2,ii,l))
 7330    continue
      end if
 7400 continue
c
c     i - boundaries
c
      do 8400 m=1,2
      kmin = 1
      kmax = kdim-1
      jmin = 1
      jmax = jdim-1
      i2   = 2
      if (i2d.eq.1) i2 = 1
      id2  = idim-2
      if (i2d.eq.1 .or. idim.eq.2) id2 = idim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         do 8120 kk=kmin,kmax
         do 8120 jj=jmin,jmax
         aa = 1. +  bci(jj,kk,m)
         bb = 1. -  bci(jj,kk,m)
         cc = 0.
         if (blank(jj,kk,1).eq.0) then
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumi0(jj,kk,l,2) = 0.5*(aa*dumi0(jj,kk,l,1)
     .                    +      bb*dum(jj,kk,1,l)
     .                    +      cc*dum(jj,kk,i2,l))
 8120    continue
c
      else
c
c      right boundary
c
         do 8330 kk=kmin,kmax
         do 8330 jj=jmin,jmax
         aa = 1. +  bci(jj,kk,m)
         bb = 1. -  bci(jj,kk,m)
         cc = 0.
         if (blank(jj,kk,idim-1).eq.0) then
            aa = 0.
            bb = 3.
            cc = -1.
         end if
         dumi0(jj,kk,l,4)  = 0.5*(aa*dumi0(jj,kk,l,3) 
     .                     +      bb*dum(jj,kk,idim-1,l)
     .                     +      cc*dum(jj,kk,id2,l))
 8330    continue
      end if
 8400 continue
c
 9000 continue   
      return
      end 
