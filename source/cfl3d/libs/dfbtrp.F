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
      subroutine dfbtrp(nvmax,n,nmax,il,iu,a,b,c,f,g,h)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform the back substitution for a periodic scalar
c     tridiagonal system of equations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax)
      dimension g(nvmax,nmax),h(nvmax,nmax),f(nvmax,nmax)
c
c      periodic solver  5-3-85
c
c      inversion of block tridiagonal...a,b,c are scalars
c      f is forcing function and solution is output in f
c      solution is by upper triangularization with unity diagonal
c      block inversions use nonpivoted lu decomposition
c      il and iu are starting and finishing indices
c      b,c,and e are overloaded
c
      il1 = il+1
      is  = il
c
c      f=binv*f
c
cdir$ ivdep
      do 1000 izz=1,n
      f(izz,is) = b(izz,is)*(f(izz,is))
 1000 continue
c
c      forward sweep
c
      iux = iu
      iux = iu-1
      do 9100 is=il1,iux
      ir  = is-1
      it  = is+1
c      first row reduction
cdir$ ivdep
      do 1001 izz=1,n
      f(izz,is) = b(izz,is)*(f(izz,is)-a(izz,is)*f(izz,ir)) 
c
c      f=binv*f
c
 1001 continue
 9100 continue
      is = iu
      ir = is-1
      it = is+1
c      first row reduction
cdir$ ivdep
      do 1002 izz=1,n
      f(izz,is) = f(izz,is)-a(izz,is)*f(izz,ir)
 1002 continue
 9121 continue
      iu2 = iu-2
      do 9101 ix=il,iu2
cdir$ ivdep
      do 1003 izz=1,n
      f(izz,is) = f(izz,is)-h(izz,ix)*f(izz,ix)
 1003 continue
 9101 continue
c
c      f=binv*f
c
cdir$ ivdep
      do 1004 izz=1,n
      f(izz,is) = b(izz,is)*(f(izz,is))
 1004 continue
c
c      back substitution
c
      iux = iu
      iux = il1
      do 9179 ii=il1,iux
      is  = il+iu-ii
      it  = is+1
cdir$ ivdep
      do 1005 izz=1,n
      f(izz,is) = f(izz,is)-c(izz,is)*f(izz,it)
 1005 continue
 9179 continue
      il11 = il1+1
      do 9180 ii=il11,iu
      is   = il+iu-ii
      it   = is+1
cdir$ ivdep
      do 1006 izz=1,n
      f(izz,is) = f(izz,is)-c(izz,is)*f(izz,it)-g(izz,is)*f(izz,iu)
 1006 continue
 9180 continue
      return
      end
