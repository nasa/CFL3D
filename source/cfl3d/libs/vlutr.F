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
      subroutine vlutr(nvmax,n,nmax,il,iu,a,b,c,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform the LU decomposition for a block 5X5 
c     tridiagonal system of equations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension a(nvmax,nmax,5,5),b(nvmax,nmax,5,5),c(nvmax,nmax,5,5)
c
c     inversion of block tridiagonal...a,b,c are 5x5 blocks
c     f is forcing function and solution is output in f
c     solution is by upper triangularization with unity diagonal
c     block inversions use nonpivoted lu decomposition
c     il and iu are starting and finishing indices
c     b,c,and e are overloaded
c
c     do 7777 k=1,5
c     do 7777 l=1,5
c     a(1,il,k,l;n) = 0.e0
c     c(1,iu,k,l;n) = 0.e0
c7777 continue
      il1 = il+1
      i   = il
c
c      l-u decomposition
c
c     nou(1) = min(nou(1)+1,ibufdim)
c     write(bou(nou(1),1),1783) nvmax,n,nmax,il,iu,(b(jv,i,1,1),jv=1,n)
c1783 format(1x,5i5/(1x,10e12.5))
cdir$ ivdep
      do 1000 izz=1,n
      b(izz,i,1,1) = 1.e0/b(izz,i,1,1)
      b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2)
      b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3)
      b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4)
      b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5)
c     b(1,i,2,1;n) = b(1,i,2,1;n)
      b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2))
      b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)
     .              -b(izz,i,2,1)* b(izz,i,1,3))
      b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)
     .              -b(izz,i,2,1)* b(izz,i,1,4))
      b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)
     .              -b(izz,i,2,1)* b(izz,i,1,5))
c     b(1,i,3,1;n) = b(1,i,3,1;n)
      b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2)
      b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
     .                                 -b(izz,i,3,2)*b(izz,i,2,3))
      b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
     .              -b(izz,i,3,1)* b(izz,i,1,4)
     .              -b(izz,i,3,2)* b(izz,i,2,4))
      b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
     .              -b(izz,i,3,1)*b(izz,i,1,5)
     .              -b(izz,i,3,2)*b(izz,i,2,5))
c     b(1,i,4,1;n) = b(1,i,4,1;n)
 1000 continue
cdir$ ivdep
      do 3000 izz=1,n
      b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2)
      b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
     .                           -b(izz,i,4,2)*b(izz,i,2,3)
      b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
     .                    -b(izz,i,4,2)*b(izz,i,2,4)
     .                    -b(izz,i,4,3)*b(izz,i,3,4))
      b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
     .                            -b(izz,i,4,1)*b(izz,i,1,5)
     .                            -b(izz,i,4,2)* b(izz,i,2,5)
     .                            -b(izz,i,4,3)*b(izz,i,3,5))
c     b(1,i,5,1;n) = b(1,i,5,1;n)
      b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2)
      b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
     .                           -b(izz,i,5,2)*b(izz,i,2,3)
      b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
     .                           -b(izz,i,5,2)*b(izz,i,2,4)
     .                           -b(izz,i,5,3)*b(izz,i,3,4)
      b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
     .                                 -b(izz,i,5,2)*b(izz,i,2,5)
     .                                 -b(izz,i,5,3)*b(izz,i,3,5)
     .                                 -b(izz,i,5,4)*b(izz,i,4,5))
 3000 continue
      if (i.eq.iu) go to 1030
c
c      c=ainv*c
c
      do 1010 m=1,5
cdir$ ivdep
      do 1001 izz=1,n
      c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m))
      c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)
     .              -b(izz,i,2,1)* c(izz,i,1,m))
      c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
     .              -b(izz,i,3,1)* c(izz,i,1,m)
     .              -b(izz,i,3,2)* c(izz,i,2,m))
      c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
     .              -b(izz,i,4,1)* c(izz,i,1,m)
     .              -b(izz,i,4,2)* c(izz,i,2,m)
     .              -b(izz,i,4,3)* c(izz,i,3,m))
      c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
     .              -b(izz,i,5,1)* c(izz,i,1,m)
     .              -b(izz,i,5,2)* c(izz,i,2,m)
     .              -b(izz,i,5,3)* c(izz,i,3,m)
     .              -b(izz,i,5,4)* c(izz,i,4,m))
      c(izz,i,4,m) = c(izz,i,4,m)-b(izz,i,4,5)*c(izz,i,5,m)
      c(izz,i,3,m) = c(izz,i,3,m)-b(izz,i,3,5)*c(izz,i,5,m)
     .                           -b(izz,i,3,4)*c(izz,i,4,m)
      c(izz,i,2,m) = c(izz,i,2,m)-b(izz,i,2,5)*c(izz,i,5,m)
     .                           -b(izz,i,2,4)*c(izz,i,4,m)
     .                           -b(izz,i,2,3)*c(izz,i,3,m)
      c(izz,i,1,m) = c(izz,i,1,m)-b(izz,i,1,5)*c(izz,i,5,m)
     .                           -b(izz,i,1,4)*c(izz,i,4,m)
     .                           -b(izz,i,1,3)*c(izz,i,3,m)
     .                           -b(izz,i,1,2)*c(izz,i,2,m)
 1001 continue
 1010 continue
 1030 continue
c
c      forward sweep
      do 100 i=il1,iu
      ir = i-1
      it = i+1
c      first row reduction
      do 20 m=1,5
      do 20 l=1,5
cdir$ ivdep
      do 1002 izz=1,n
      b(izz,i,m,l) = b(izz,i,m,l)-a(izz,i,m,1)*c(izz,ir,1,l)
     .                           -a(izz,i,m,2)*c(izz,ir,2,l)
     .                           -a(izz,i,m,3)*c(izz,ir,3,l)
     .                           -a(izz,i,m,4)*c(izz,ir,4,l)
     .                           -a(izz,i,m,5)*c(izz,ir,5,l)
 1002 continue
   20 continue
c
c      l-u decomposition
c
cdir$ ivdep
      do 1003 izz=1,n
      b(izz,i,1,1) = 1.e0/b(izz,i,1,1)
      b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2)
      b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3)
      b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4)
      b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5)
c     b(1,i,2,1;n) = b(1,i,2,1;n)
      b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2))
      b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)
     .              -b(izz,i,2,1)* b(izz,i,1,3))
      b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)
     .              -b(izz,i,2,1)* b(izz,i,1,4))
      b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)
     .              -b(izz,i,2,1)* b(izz,i,1,5))
c     b(1,i,3,1;n) = b(1,i,3,1;n)
      b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2)
      b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
     .                                 -b(izz,i,3,2)*b(izz,i,2,3))
      b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
     .              -b(izz,i,3,1)* b(izz,i,1,4)
     .              -b(izz,i,3,2)* b(izz,i,2,4))
      b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
     .              -b(izz,i,3,1)* b(izz,i,1,5)
     .              -b(izz,i,3,2)* b(izz,i,2,5))
c     b(1,i,4,1;n) = b(1,i,4,1;n)
 1003 continue
cdir$ ivdep
      do 3003 izz=1,n
      b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2)
      b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
     .                           -b(izz,i,4,2)*b(izz,i,2,3)
      b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
     .                                 -b(izz,i,4,2)*b(izz,i,2,4)
     .                                 -b(izz,i,4,3)*b(izz,i,3,4))
      b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
     .              -b(izz,i,4,1)* b(izz,i,1,5)
     .              -b(izz,i,4,2)* b(izz,i,2,5)
     .              -b(izz,i,4,3)* b(izz,i,3,5))
c     b(1,i,5,1;n) = b(1,i,5,1;n)
      b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2)
      b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
     .                           -b(izz,i,5,2)*b(izz,i,2,3)
      b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
     .                           -b(izz,i,5,2)*b(izz,i,2,4)
     .                           -b(izz,i,5,3)*b(izz,i,3,4)
      b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
     .                                 -b(izz,i,5,2)*b(izz,i,2,5)
     .                                 -b(izz,i,5,3)*b(izz,i,3,5)
     .                                 -b(izz,i,5,4)*b(izz,i,4,5))
 3003 continue
      if (i.eq.iu) go to 1130
c
c      c=ainv*c
c
      do 1110 m=1,5
cdir$ ivdep
      do 1004 izz=1,n
      c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m))
      c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)
     .              -b(izz,i,2,1)* c(izz,i,1,m))
      c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
     .              -b(izz,i,3,1)* c(izz,i,1,m)
     .              -b(izz,i,3,2)* c(izz,i,2,m))
      c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
     .              -b(izz,i,4,1)* c(izz,i,1,m)
     .              -b(izz,i,4,2)* c(izz,i,2,m)
     .              -b(izz,i,4,3)* c(izz,i,3,m))
      c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
     .              -b(izz,i,5,1)* c(izz,i,1,m)
     .              -b(izz,i,5,2)* c(izz,i,2,m)
     .              -b(izz,i,5,3)* c(izz,i,3,m)
     .              -b(izz,i,5,4)*c(izz,i,4,m))
      c(izz,i,4,m) = c(izz,i,4,m)-b(izz,i,4,5)*c(izz,i,5,m)
      c(izz,i,3,m) = c(izz,i,3,m)-b(izz,i,3,5)*c(izz,i,5,m)
     .                           -b(izz,i,3,4)*c(izz,i,4,m)
      c(izz,i,2,m) = c(izz,i,2,m)-b(izz,i,2,5)*c(izz,i,5,m)
     .                           -b(izz,i,2,4)*c(izz,i,4,m)
     .                           -b(izz,i,2,3)*c(izz,i,3,m)
      c(izz,i,1,m) = c(izz,i,1,m)-b(izz,i,1,5)*c(izz,i,5,m)
     .                           -b(izz,i,1,4)*c(izz,i,4,m)
     .                           -b(izz,i,1,3)*c(izz,i,3,m)
     .                           -b(izz,i,1,2)*c(izz,i,2,m)
 1004 continue
 1110 continue
 1130 continue
  100 continue
      return
      end
