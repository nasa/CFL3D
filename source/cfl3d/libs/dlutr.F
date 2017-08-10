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
      subroutine dlutr(nvmax,n,nmax,il,iu,a,b,c)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform the scalar tridiagonal (LU) decomposition.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax)
c
c      inversion of scalar tridiagonal...a,b,c are scalars
c      f is forcing function and solution is output in f
c      solution is by upper triangularization with unity diagonal
c      block inversions use nonpivoted lu decomposition
c      il and iu are starting and finishing indices
c      b and c are overloaded
c
      il1 = il+1
      i   = il
c
c      l-u decomposition
c
cdir$ ivdep
      do 1000 izz=1,n
      b(izz,i) = 1.e0/b(izz,i)
 1000 continue
c
      if (i.eq.iu) go to 1030
c
c      c=ainv*c
c
cdir$ ivdep
      do 1001 izz=1,n
      c(izz,i) = b(izz,i)*c(izz,i)
 1001 continue
 1030 continue
c      forward sweep
      do 100 i=il1,iu
      ir = i-1
      it = i+1
c      first row reduction
cdir$ ivdep
      do 1002 izz=1,n
      b(izz,i) = b(izz,i)-a(izz,i)*c(izz,ir)
c
c      l-u decomposition
c
      b(izz,i) = 1.e0/b(izz,i)
 1002 continue
      if (i.eq.iu) go to 100
c
c      c=ainv*c
c
cdir$ ivdep
      do 1003 izz=1,n
      c(izz,i) = b(izz,i)*c(izz,i)
 1003 continue
  100 continue
      return
      end
