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
      subroutine bsub(id1,id2,a,b,c,f,i1,i2,il,iu)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Performs the back substitution for a block 5x5 tridi-
c     agonal matrix equation solution.  The vectorization is over points
c     i1-i2 and the tridiagonal matrix equation spans points il-iu.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(id1,id2,5,5),b(id1,id2,5,5),c(id1,id2,5,5)
      dimension f(id1,id2,5)
c
c     solve matrix equation
c
      il1 = il+1
      is  = il
c
c      f=binv*f
c
cdir$ ivdep
      do 100 i=i1,i2
      t1        = b(i,is,1,1)*(f(i,is,1))
      t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1)
      t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2)
      t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
     .                                  -b(i,is,4,3)*t3)
      f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
     .                                  -b(i,is,5,3)*t3-b(i,is,5,4)*t4)
      f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5)
      f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4)
      f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
     .              -b(i,is,2,3)*f(i,is,3)
      f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
     .              -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2)
  100 continue
c
c      forward sweep
c 
      do 200 is=il1,iu
      ir = is-1
      it = is+1
c      first row reduction
      do 190 m=1,5
cdir$ ivdep
      do 102 i=i1,i2
      f(i,is,m) = f(i,is,m)-a(i,is,m,1)*f(i,ir,1)
     .                     -a(i,is,m,2)*f(i,ir,2)
     .                     -a(i,is,m,3)*f(i,ir,3)
     .                     -a(i,is,m,4)*f(i,ir,4)
     .                     -a(i,is,m,5)*f(i,ir,5)
  102 continue
  190 continue
c
c      f=binv*f
c
cdir$ ivdep
      do 300 i=i1,i2
      t1        = b(i,is,1,1)*(f(i,is,1))
      t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1)
      t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2)
      t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
     .                                  -b(i,is,4,3)*t3)
      f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
     .                                -b(i,is,5,3)*t3-b(i,is,5,4)*t4)
      f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5)
      f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4)
      f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
     .              -b(i,is,2,3)*f(i,is,3)
      f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
     .              -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2)
  300 continue
  200 continue
c
c      back substitution
c
      do 400 iqq=il1,iu
      is = il+iu-iqq
      it = is+1
      do 415 m=1,5
cdir$ ivdep
      do 405 i=i1,i2
      f(i,is,m) =  f(i,is,m)-c(i,is,m,1)*f(i,it,1)
     .                      -c(i,is,m,2)*f(i,it,2)
     .                      -c(i,is,m,3)*f(i,it,3)
     .                      -c(i,is,m,4)*f(i,it,4)
     .                      -c(i,is,m,5)*f(i,it,5)
  405 continue
  415 continue
  400 continue
      return
      end
