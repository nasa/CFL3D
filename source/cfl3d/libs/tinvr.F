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
      subroutine tinvr(n,r1,r2,r3,r4,r5,kx,ky,kz,lx,ly,lz,mx,my,mz,c,ub,
     .                 rho,u,v,w,max,itinv,eig2,eig3,xm2a)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Multiply the inverse of the diagonalizing matrix
c     T times the residual contribution.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension r1(max),r2(max),r3(max),r4(max),r5(max)
      dimension  c(max),ub(max), u(max), v(max), w(max),rho(max)
      dimension  xm2a(max),eig2(max),eig3(max)
      dimension kx(max),ky(max),kz(max)
      dimension lx(max),ly(max),lz(max)
      dimension mx(max),my(max),mz(max)
c
#   ifdef CMPLX
      complex kx,ky,kz,lx,ly,lz,mx,my,mz,lxi,lyi,lzi,mxi
#   else
      real kx,ky,kz,lx,ly,lz,mx,my,mz,lxi,lyi,lzi,mxi
#   endif
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /precond/ cprec,uref,avn
c
cdir$ ivdep
      do 500 m=1,n
      lxi   =  ky(m)-kz(m) 
      lyi   =  kz(m)-kx(m) 
      lzi   =  kx(m)-ky(m) 
      mxi   =  lxi*lxi+lyi*lyi+lzi*lzi 
      mxi   =  1.e0/sqrt(mxi)
c
      lx(m) =  lxi*mxi 
      ly(m) =  lyi*mxi 
      lz(m) =  lzi*mxi 
      mx(m) =  ky(m)*lz(m)-ly(m)*kz(m)
      my(m) = -kx(m)*lz(m)+lx(m)*kz(m)
      mz(m) =  kx(m)*ly(m)-lx(m)*ky(m)
  500 continue
c
      if (real(cprec) .eq. 0.) then
c
      if (itinv.eq.0) then
c
c     T(inverse)*R
c
cdir$ ivdep
      do 1000 m=1,n
      vb    =  u(m)*lx(m)+v(m)*ly(m)+w(m)*lz(m)
      wb    =  u(m)*mx(m)+v(m)*my(m)+w(m)*mz(m)
c
      q2    =  0.5*(u(m)*u(m)+v(m)*v(m)+w(m)*w(m))
c
      t3    =  1.e0/c(m)
      t1    = -q2*r1(m)-r5(m)+u(m)*r2(m)+v(m)*r3(m)+w(m)*r4(m)
      t1    =  gm1*t1*t3*t3
      t2    = -ub(m)*r1(m)+kx(m)*r2(m)+ky(m)*r3(m)+kz(m)*r4(m)
      t2    =  t2*t3
      t3    = -vb*r1(m)+lx(m)*r2(m)+ly(m)*r3(m)+lz(m)*r4(m)
c
      r3(m) = -wb*r1(m)+mx(m)*r2(m)+my(m)*r3(m)+mz(m)*r4(m)
      r1(m) =  r1(m)+t1
      r2(m) =  t3
      r4(m) =  0.5e0*(t2-t1)
      r5(m) =  r4(m)-t2
 1000 continue
c
      else
c
c     T(inverse)*M*R
c
cdir$ ivdep
      do 2000 m=1,n
      t1    =  1.e0/c(m)
c
      r5(m) =  r5(m)*t1*t1
      r1(m) =  r1(m)-r5(m)
c
      t2    =  rho(m)*r2(m)
      t3    =  rho(m)*r3(m)
      t4    =  rho(m)*r4(m)
c
      r2(m) =           lx(m)*t2+ly(m)*t3+lz(m)*t4
      r3(m) =           mx(m)*t2+my(m)*t3+mz(m)*t4
      r4(m) =  0.5*(t1*(kx(m)*t2+ky(m)*t3+kz(m)*t4)+r5(m))
      r5(m) = -r4(m)+r5(m)
 2000 continue
      end if
c
      else
c
      if (itinv.eq.0) then
c
c     T(inverse)*M*R
c
cdir$ ivdep
      do 10001 m=1,n

      t1    =  1.e0/c(m)
      rrho = 1.e0/rho(m)
      xm2 = xm2a(m)*t1
      xm2ar = 1.0/xm2a(m)
      fplus = (eig2(m)-ub(m))*xm2ar
      fmins = -(eig3(m)-ub(m))*xm2ar
 
      r11 = r1(m)
      r21 = r2(m)
      r31 = r3(m)
      r41 = r4(m)
      r51 = r5(m)

      vmag1 = u(m)**2 + v(m)**2 + w(m)**2
      r5t = gm1*(0.5*vmag1*r11 
     .    - (u(m)*r21 + v(m)*r31 + w(m)*r41) + r51) 
c
c ---- multiplication by inverse of precond. matrix
c
      r1(m) = r11 - (1.-xm2)*r5t*t1*t1 
      r2(m) = rrho*(-u(m)*r11 + r21)
      r3(m) = rrho*(-v(m)*r11 + r31)
      r4(m) = rrho*(-w(m)*r11 + r41)
      r5(m) = xm2*r5t
c
c ---- multiplication by T(inverse)
c
      r5t =  r5(m)*t1*t1
      r1(m) =  r1(m)-r5t
c
      t2    =  rho(m)*r2(m)
      t3    =  rho(m)*r3(m)
      t4    =  rho(m)*r4(m)
c
      r2(m) =           lx(m)*t2+ly(m)*t3+lz(m)*t4
      r3(m) =           mx(m)*t2+my(m)*t3+mz(m)*t4
      r4(m) =  0.5*(t1*(kx(m)*t2+ky(m)*t3+kz(m)*t4)
     .       + r5t*fplus)
      r5(m) =  -0.5*(t1*(kx(m)*t2+ky(m)*t3+kz(m)*t4)
     .       - r5t*fmins)

10001 continue
c
      else
c
c     T(inverse)*R
c
cdir$ ivdep
      do 20001 m=1,n
c
      t1    =  1.e0/c(m)
      xm2ar = 1.0/xm2a(m)
      fplus = (eig2(m)-ub(m))*xm2ar
      fmins = -(eig3(m)-ub(m))*xm2ar
c
      r5t =  r5(m)*t1*t1
      r1(m) =  r1(m)-r5t
c
      t2    =  rho(m)*r2(m)
      t3    =  rho(m)*r3(m)
      t4    =  rho(m)*r4(m)
c
      r2(m) =           lx(m)*t2+ly(m)*t3+lz(m)*t4
      r3(m) =           mx(m)*t2+my(m)*t3+mz(m)*t4
      r4(m) =  0.5*(t1*(kx(m)*t2+ky(m)*t3+kz(m)*t4)
     .       + r5t*fplus)
      r5(m) =  -0.5*(t1*(kx(m)*t2+ky(m)*t3+kz(m)*t4)
     .       - r5t*fmins)
20001 continue
      end if
c
      end if
c
      return
      end
