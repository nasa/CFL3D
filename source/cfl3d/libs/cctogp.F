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
      subroutine cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,dum,
     .                  dumi0,dumj0,dumk0,jdw,kdw,idw,dumgp,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Convert from cell-center data to grid-point data and 
c               load into a single precision array for plot3d output
c               Cell-center data is assumed to have the format of
c               q/qj0/qk0/qi0 arrays or vist3d/vj0/vk0/vi0 arrays
c
c     NOTE: corners and edges are not done very well!
c
c     dum...input array containing cell-center interior data
c
c     dumi0/dumj0/dumk0...input arrays of cell-center boundary data
c
c     dumgp...output array with grid point data (single precision)
c***********************************************************************
c
#if defined ADP_OFF
#   ifdef CMPLX
#     ifdef DBLE_PRECSN
      implicit complex*8(a-h,o-z)
#     else
      implicit complex(a-h,o-z)
#     endif
#   else
#     ifdef DBLE_PRECSN
      implicit real*8 (a-h,o-z)
#     endif
#   endif
#else
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#endif
c
      dimension dumgp(jdw,kdw,idw,ldim)
      dimension dum(jdim,kdim,idim,ldim), dumi0(jdim,kdim,ldim,4),
     .          dumj0(kdim,idim-1,ldim,4),dumk0(jdim,idim-1,ldim,4)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     determine values at grid points
c
      iw = 0
      do 1000 i=i1,i2,i3
      iw = iw+1
      kw = 0
      do 1001 k=k1,k2,k3
      kw = kw+1
      jw = 0
      do 1002 j=j1,j2,j3
      jw = jw+1
      do 1003 l=1,ldim
      id  = i
      id1 = id-1
      if (id1.le.0)   id1 = 1
      if (id.gt.idim1) id = idim1
      jd = j
      kd = k
c
      if (k.eq.1 .or. k.eq.kdim) then
c     k=1/k=kdim faces
         if (k.eq.kdim) kd = kdim1
         if (j.eq.1 .or. j.eq.jdim) then
c           edge points and corner points
            if (j.eq.jdim) jd = jdim1
            dumgp(jw,kw,iw,l) = 0.5*(dum(jd,kd,id,l) + dum(jd,kd,id1,l))
         else
c           interior points on  k=1/k=kdim faces
            m = 2
            if (k.eq.kdim) m = 4
            dumgp(jw,kw,iw,l) = .25*(dumk0(jd,id,l,m)
     .                        +      dumk0(jd-1,id,l,m)
     .                        +      dumk0(jd,id1,l,m) 
     .                        +      dumk0(jd-1,id1,l,m))
         end if
      else if (j.eq.1 .or. j.eq.jdim) then
c        interior points on j=1/j=jdim faces
         if (j.eq.jdim) jd = jdim1
         m = 2
         if (j.eq.jdim) m = 4
         dumgp(jw,kw,iw,l) = .25*(dumj0(kd,id,l,m)
     .                     +      dumj0(kd-1,id,l,m)
     .                     +      dumj0(kd,id1,l,m)  
     .                     +      dumj0(kd-1,id1,l,m))
      else
c     interior points
         if (i.eq.1 .or. i.eq.idim) go to 1003
         dumgp(jw,kw,iw,l) = .125*(dum(j,k,i,l)
     .                     +       dum(j,k,i-1,l)
     .                     +       dum(j-1,k,i,l)
     .                     +       dum(j-1,k,i-1,l)
     .                     +       dum(j,k-1,i,l)
     .                     +       dum(j,k-1,i-1,l)
     .                     +       dum(j-1,k-1,i,l)
     .                     +       dum(j-1,k-1,i-1,l))
      end if
 1003 continue
 1002 continue
 1001 continue
 1000 continue
c
c     interior points on i=1/i=idim faces
      iw = 0
      do 1008 i=i1,i2,i3
      iw = iw+1
      if (i.ne.1 .and. i.ne.idim) go to 1008
      m  = 2
      if (i.eq.idim) m = 4
      kw = 0
      do 1007 k=k1,k2,k3
      kw = kw+1
      jw = 0
      do 1006 j=j1,j2,j3
      jw = jw+1
      if (j.eq.1 .or. j.eq.jdim .or. k.eq.1 .or. k.eq.kdim) go to 1006
c     interior points
      do 1005 l=1,ldim
      dumgp(jw,kw,iw,l) = .25*(dumi0(j,k,l,m)   + dumi0(j,k-1,l,m)
     .                  +      dumi0(j-1,k,l,m) + dumi0(j-1,k-1,l,m))
 1005 continue
 1006 continue
 1007 continue
 1008 continue
c
      return
      end
