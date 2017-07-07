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
      subroutine ffluxl(k,npl,xkap,idf,jdim,kdim,idim,res,q,qi0,si,
     .                  dfp,dfm,t,nvtq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute the left-hand flux contributions due to the
c     inviscid terms for the I-direction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension si(jdim*kdim,idim,5)
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),t(nvtq,39)
      dimension dfp(npl*(jdim-1)*idim,5,5),dfm(npl*(jdim-1)*idim,5,5)
c
      common /fvfds/ rkap0(3),ifds(3)
c
c     f-flux-p-m jacobians
c
      jdim1 = jdim-1
      idim1 = idim-1
      jv    = npl*jdim1
      js    = jv+1
c
      do 8008 kpl=1,npl
      kk = k+kpl-1
      do 8008 l=1,5
      do 8008 j=1,jdim
      q(j,kk,idim,l) = qi0(j,kk,l,3)
 8008 continue
c
      do 8000 i=1,idim
      jiv  = (i-1)*jdim1*npl + 1
      do 8000 kpl=1,npl
      jiv1 = jiv + (kpl-1)*jdim1
      kk   = k+kpl-1
      ji   = (kk-1)*jdim + 1
      do 9000 l=1,4
cdir$ ivdep
      do 1000 izz=1,jdim1
      t(izz+jiv1-1,35+l) = si(izz+ji-1,i,l)
 1000 continue
 9000 continue
cdir$ ivdep
      do 1001 izz=1,jdim1
      t(izz+jiv1-1,6) = si(izz+ji-1,i,5)
 1001 continue
 8000 continue
c
c     create transposed q
c
      n   = npl*jdim1*idim
      do 100 kpl=1,npl
      kk  = k+kpl-1
      jiv = (kpl-1)*jdim1 + 1
      do 100 l=1,5
cdir$ ivdep
      do 1002 izz=1,jdim1
      t(izz+jiv-1,15+l) = qi0(izz,kk,l,1)
 1002 continue
      do 100 i=1,idim
      jiv1 = (i-1)*jdim1*npl + jiv
cdir$ ivdep
      do 1003 izz=1,jdim1
      t(izz+jiv1-1,25+l) = q(izz,kk,i,l)
 1003 continue
  100 continue
c
      n0 = n-jv
c
c     first order
c
      do 40 l=1,5
cdir$ ivdep
      do 1004 izz=1,n0
      t(izz+js-1,20+l) = t(izz,25+l)
 1004 continue
cdir$ ivdep
      do 1005 izz=1,jv
      t(izz,20+l) = t(izz,15+l)
 1005 continue
   40 continue
c  put unsteady metrics into t(20)
cdir$ ivdep
      do 1006 izz=1,n
      t(izz,20) = t(izz,6)
 1006 continue
c
      if (ifds(1).eq.0) then
         call dfluxpm(t(1,1), t(1,2), t(1,36),t(1,37),t(1,38),t(1,39),
     .                t(1,20),t(1,21),dfp,n,n,nvtq,+1)
      else
         call dfhat(t(1,36),t(1,37),t(1,38),t(1,39),t(1,20),t(1,21),dfp,
     .              n,nvtq,+1)
      end if
c
      if (ifds(1).eq.0) then
         call dfluxpm(t(1,1), t(1,2), t(1,36),t(1,37),t(1,38),t(1,39),
     .                t(1,20),t(1,26),dfm,n,n,nvtq,-1)
      else
         call dfhat(t(1,36),t(1,37),t(1,38),t(1,39),t(1,20),t(1,26),dfm,
     .              n,nvtq,-1)
      end if
      return
      end
