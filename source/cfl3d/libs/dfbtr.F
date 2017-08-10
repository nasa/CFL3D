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
      subroutine dfbtr(nvmax,n,nmax,il,iu,a,b,c,f)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform the back substitution for a scalar
c     tridiagonal system of equations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax),f(nvmax,nmax)
c
      il1 = il+1
      is  = il
c
c      f=binv*f
c
cdir$ ivdep
      do 1000 izz=1,n
      f(izz,is) = b(izz,is)*f(izz,is)
 1000 continue
c
c      forward sweep
c
      do 9100 is=il1,iu
      ir = is-1
      it = is+1
c      first row reduction
cdir$ ivdep
      do 1001 izz=1,n
c
c      f=binv*f
c
      f(izz,is) = b(izz,is)*(f(izz,is)-a(izz,is)*f(izz,ir)) 
 1001 continue
 9100 continue
c      back substitution
      do 9180 ii=il1,iu
      is = il+iu-ii
      it = is+1
cdir$ ivdep
      do 1002 izz=1,n
      f(izz,is) = f(izz,is)-c(izz,is)*f(izz,it)
 1002 continue
 9180 continue
      return
      end
