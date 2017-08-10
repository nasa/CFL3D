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
      subroutine xtbatb(jdim,kdim,idim,xtbj,xtbk,xtbi,atbj,atbk,atbi,
     .                  t,tti,ttj,ttk)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Store grid speeds and accelerations on the boundaries
c     for use in setting no-slip and wall pressure boundary conditions
c     for dynamic meshes
c
c     xtbi - velocity components of i0/idim boundaries (face centers)
c     xtbj - velocity components of j0/jdim boundaries (face centers)
c     xtbk - velocity components of k0/kdim boundaries (face centers)
c     atbi - acceleration components of i0/idim boundaries (face centers)
c     atbj - acceleration components of j0/jdim boundaries (face centers)
c     atbk - acceleration components of k0/kdim boundaries (face centers)
c
c     t   is a temp array containing grid velocity at grid points
c     tti is a temp array containing acceleration at i-boundary points
c     ttj is a temp array containing acceleration at j-boundary points
c     ttk is a temp array containing acceleration at k-boundary points
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension t(jdim*kdim*idim,3),tti(jdim*kdim,3,2),
     .          ttj(kdim*idim,3,2),ttk(jdim*idim,3,2)
c
      kdim1 = kdim-1
      idim1 = idim-1
      jdkd  = jdim*kdim
      jdid  = jdim*idim
      idkd  = idim*kdim
c
      do 10 i=1,idim1
c
c**********************************************
c     j=1/jdim boundary velocities
c     use face-averaged values
c**********************************************
c
      izz1 = (i-1)*jdkd
      do 100 nn = 1,2
      j = 1
      if (nn.eq.2) j = jdim
      do 110 k=1,kdim1
      izz     = izz1 + (k-1)*jdim + j
      izzk1   = izz  + jdim
      izzi1   = izz  + jdkd
      izzi1k1 = izz  + jdkd+jdim
      xtbj(k,i,1,nn) = 0.25*( t(izz,1)   + t(izzk1,1)  +
     .                        t(izzi1,1) + t(izzi1k1,1) )
      xtbj(k,i,2,nn) = 0.25*( t(izz,2)   + t(izzk1,2)  +
     .                        t(izzi1,2) + t(izzi1k1,2) )
      xtbj(k,i,3,nn) = 0.25*( t(izz,3)   + t(izzk1,3)  +
     .                        t(izzi1,3) + t(izzi1k1,3) )
  110 continue
c
c     zero out extra layer of boundary faces for saftey
      xtbj(kdim,i,1,nn) = 0.
      xtbj(kdim,i,2,nn) = 0.
      xtbj(kdim,i,3,nn) = 0.
c
  100 continue
c
c**********************************************
c     j=1/jdim boundary accelerations
c     use face-averaged values
c**********************************************
c
      izz0 = (i-1)*kdim
      do 120 nn = 1,2
      do 130 k=1,kdim-1
      izz     = izz0 + k
      izzk1   = izz  + 1
      izzi1   = izz  + kdim
      izzi1k1 = izz  + kdim + 1 
      atbj(k,i,1,nn) = 0.25*( ttj(izz,1,nn)   + ttj(izzk1,1,nn)  +
     .                        ttj(izzi1,1,nn) + ttj(izzi1k1,1,nn) )
      atbj(k,i,2,nn) = 0.25*( ttj(izz,2,nn)   + ttj(izzk1,2,nn)  +
     .                        ttj(izzi1,2,nn) + ttj(izzi1k1,2,nn) )
      atbj(k,i,3,nn) = 0.25*( ttj(izz,3,nn)   + ttj(izzk1,3,nn)  +
     .                        ttj(izzi1,3,nn) + ttj(izzi1k1,3,nn) )
  130 continue
c
c     zero out extra layer of boundary faces for saftey
      atbj(kdim,i,1,nn) = 0.
      atbj(kdim,i,2,nn) = 0.
      atbj(kdim,i,3,nn) = 0.
c
  120 continue
c
c**********************************************
c     k=1/kdim boundary velocities
c     use face-averaged values
c**********************************************
c
      izz0 = (i-1)*jdkd
      do 140 nn = 1,2
      k = 1
      if (nn.eq.2) k = kdim
      do 150 j=1,jdim-1
      izz     = izz0 + (k-1)*jdim + j
      izzi1   = izz  + jdkd
      izzj1   = izz  + 1
      izzi1j1 = izz  + jdkd + 1
      xtbk(j,i,1,nn) = 0.25*( t(izz,1)   + t(izzj1,1)  +
     .                        t(izzi1,1) + t(izzi1j1,1) )
      xtbk(j,i,2,nn) = 0.25*( t(izz,2)   + t(izzj1,2)  +
     .                        t(izzi1,2) + t(izzi1j1,2) )
      xtbk(j,i,3,nn) = 0.25*( t(izz,3)   + t(izzj1,3)  +
     .                        t(izzi1,3) + t(izzi1j1,3) )
  150 continue
c
c     zero out extra layer of boundary faces for saftey
      xtbk(jdim,i,1,nn) = 0.
      xtbk(jdim,i,2,nn) = 0.
      xtbk(jdim,i,3,nn) = 0.
c
  140 continue
c
c**********************************************
c     k=1/kdim boundary accelerations
c     use face-averaged values
c**********************************************
c
      do 160 nn = 1,2
      do 170 j=1,jdim-1
      izz     = (i-1)*jdim+j
      izzj1   = izz + 1
      izzi1   = izz + jdim
      izzi1j1 = izz + jdim + 1
      atbk(j,i,1,nn) = 0.25*( ttk(izz,1,nn)   + ttk(izzj1,1,nn)  +
     .                        ttk(izzi1,1,nn) + ttk(izzi1j1,1,nn) )
      atbk(j,i,2,nn) = 0.25*( ttk(izz,2,nn)   + ttk(izzj1,2,nn)  +
     .                        ttk(izzi1,2,nn) + ttk(izzi1j1,2,nn) )
      atbk(j,i,3,nn) = 0.25*( ttk(izz,3,nn)   + ttk(izzj1,3,nn)  +
     .                        ttk(izzi1,3,nn) + ttk(izzi1j1,3,nn) )
  170 continue
c
c     zero out extra layer of boundary faces for saftey
      atbk(jdim,i,1,nn) = 0.
      atbk(jdim,i,2,nn) = 0.
      atbk(jdim,i,3,nn) = 0.
c
  160 continue
c
   10 continue
c
c**********************************************
c     i=1/idim boundary velocities
c     use face-averaged values
c**********************************************
c
      do 180 nn = 1,2
      i = 1
      if (nn.eq.2) i = idim
      izz0 = (i-1)*jdkd
      do 190 j=1,jdim-1
      do 190 k=1,kdim-1
      izz = izz0 + (k-1)*jdim + j
      izzj1 = izz + 1
      izzk1 = izz + jdim
      izzj1k1 = izz + 1 + jdim
      xtbi(j,k,1,nn) = 0.25*( t(izz,1)   + t(izzj1,1)  +
     .                        t(izzk1,1) + t(izzj1k1,1) )
      xtbi(j,k,2,nn) = 0.25*( t(izz,2)   + t(izzj1,2)  +
     .                        t(izzk1,2) + t(izzj1k1,2) )
      xtbi(j,k,3,nn) = 0.25*( t(izz,3)   + t(izzj1,3)  +
     .                        t(izzk1,3) + t(izzj1k1,3) )
  190 continue
c
c     zero out extra layer of boundary faces for saftey
      do 200 j=1,jdim
      xtbi(j,kdim,1,nn) = 0.
      xtbi(j,kdim,2,nn) = 0.
      xtbi(j,kdim,3,nn) = 0.
  200 continue
c
c     zero out extra layer of boundary faces for saftey
      do 210 k=1,kdim
      xtbi(jdim,k,1,nn) = 0.
      xtbi(jdim,k,2,nn) = 0.
      xtbi(jdim,k,3,nn) = 0.
  210 continue
c
  180 continue
c
c**********************************************
c     i=1/idim boundary accelerations
c     use face-averaged values
c**********************************************
c
      do 220 nn=1,2     
      do 230 j=1,jdim-1
      do 230 k=1,kdim-1
      izz     = (k-1)*jdim+j
      izzj1   = izz + 1
      izzk1   = izz + jdim
      izzj1k1 = izz + jdim + 1
      atbi(j,k,1,nn) = 0.25*( tti(izz,1,nn)   + tti(izzj1,1,nn)  +
     .                        tti(izzk1,1,nn) + tti(izzj1k1,1,nn) )
      atbi(j,k,2,nn) = 0.25*( tti(izz,2,nn)   + tti(izzj1,2,nn)  +
     .                        tti(izzk1,2,nn) + tti(izzj1k1,2,nn) )
      atbi(j,k,3,nn) = 0.25*( tti(izz,3,nn)   + tti(izzj1,3,nn)  +
     .                        tti(izzk1,3,nn) + tti(izzj1k1,3,nn) )
  230 continue
c
c     zero out extra layer of boundary faces for saftey
      do 240 k=1,kdim
      atbi(jdim,k,1,nn) = 0.
      atbi(jdim,k,2,nn) = 0.
      atbi(jdim,k,3,nn) = 0.
  240 continue
c
c     zero out extra layer of boundary faces for saftey
      do 250 j=1,jdim
      atbi(j,kdim,1,nn) = 0.
      atbi(j,kdim,2,nn) = 0.
      atbi(j,kdim,3,nn) = 0.
  250 continue
c
  220 continue
c
      return
      end
