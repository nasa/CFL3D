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
      subroutine tdq(n,r1,r2,r3,r4,r5,kx,ky,kz,lx,ly,lz,mx,my,mz,c,
     .               ub,rho,u,v,w,max,eig2,eig3,xm2a)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Multiply the inverse of the diagonalizing matrix
c     T times the change in characteristic combination of variables.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension  r1(max),r2(max),r3(max),r4(max),r5(max)
      dimension  c(max),ub(max), u(max), v(max), w(max),rho(max)
      dimension  xm2a(max),eig2(max), eig3(max)
      dimension  kx(max),ky(max),kz(max)
      dimension  lx(max),ly(max),lz(max)
      dimension  mx(max),my(max),mz(max)
c
#   ifdef CMPLX
      complex kx,ky,kz,lx,ly,lz,mx,my,mz
#   else
      real kx,ky,kz,lx,ly,lz,mx,my,mz
#   endif
c
      common /precond/ cprec,uref,avn
c
c     M(inverse)*T*R
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1000 m=1,n
         t1      = 1.0/rho(m)
         t2      = t1*r2(m)
         t3      = t1*r3(m)
         t5      = t1* c(m)*(r4(m)-r5(m))
c
         r5(m) = r4(m)+r5(m)
         r1(m) = r1(m)+r5(m)
         r5(m) = c(m)*c(m)*r5(m)
c
         r2(m) = lx(m)*t2+mx(m)*t3+kx(m)*t5
         r3(m) = ly(m)*t2+my(m)*t3+ky(m)*t5
         r4(m) = lz(m)*t2+mz(m)*t3+kz(m)*t5
 1000    continue
      else
cdir$ ivdep
         do 10001 m=1,n
c
c        modifications for preconditioning
c
         xm2ar = 1.0/xm2a(m)
         fplus = (eig2(m)-ub(m))*xm2ar
         fmins = -(eig3(m)-ub(m))*xm2ar
         fsum = 2.0/(fplus+fmins)
c
         t1      = 1.0/rho(m)
         t2      = t1*r2(m)
         t3      = t1*r3(m)
         t5      = t1* c(m)*(fmins*r4(m)-fplus*r5(m))*fsum
c
         r5(m) = (r4(m)+r5(m))*fsum
         r1(m) = r1(m)+r5(m)
         r5(m) = c(m)*c(m)*r5(m)
c
         r2(m) = lx(m)*t2+mx(m)*t3+kx(m)*t5
         r3(m) = ly(m)*t2+my(m)*t3+ky(m)*t5
         r4(m) = lz(m)*t2+mz(m)*t3+kz(m)*t5
10001    continue
      end if
      return
      end
