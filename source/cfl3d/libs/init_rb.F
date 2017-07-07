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
      subroutine init_rb(maxaes,maxbl,zorig,aesrfdat) 
c
c     $Id$
c
c***********************************************************************
c     Purpose: initialize rigid body dynamic data, and calculate the 
c              stm, stmi, and bmat arrays used for the solution to the
c              dynamic equations, where stm is the state transition
c              matrix, stmi is the integral of the stm, and bmat
c              is the array containing the generalized masses. These
c              arrays depend only on  uinf and ref. length  and time 
c              step, so need only to be calculated once, at the start
c              of a calculation, for constant time step. Note that 
c              ainf=uinf/xmach is used for the non-dimensionalization 
c              of time in the rigid body equations of motion.
c              
c        Reference: Vinh, Lam-Son, Edwards, J.W., Seidel, D. A., Batina,
c                   J. T., "Transonic Stability and Control of Aircraft
c                   Using CFD Methods," AIAA Paper 88-4374.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension temp(4,4)
      dimension aesrfdat(5,maxaes),zorig(maxbl)
c
      common /rbstmt1/ bmatrb(4,4),gforcnrb(4),gforcnmrb(4),gforcsrb(4),
     .                 stmrb(4,4),stmirb(4,4),xsrb(4),xxnrb(4),xnmrb(4),
     .                 x0rb(4)
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat 
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
c
      xcg      = xmc
      crefrb   = cref
      xtmref   = xcg - xmc
      areat    = sref
      grefl    = aesrfdat(2,1)
      uinfrb   = aesrfdat(3,1)
      qinfrb   = aesrfdat(4,1)
      ainf     = uinfrb/xmach
      timesc   = grefl/ainf
      dts      = dt*timesc
      zrg0     = zorig(1)
c
      do m=1,4
         do n = 1,4
            stmrb(n,m)   = 0.
            stmirb(n,m)  = 0.
            bmatrb(n,m)  = 0.
         end do
         gforcnrb(m)  = 0.
         gforcnmrb(m) = 0.
         gforcsrb(m)  = 0.
         x0rb(m)      = 0.
         xxnrb(m)     = x0rb(m)
         xsrb(m)      = xxnrb(m)
      end do
c
      do i = 1,4
         stmrb(i,i)  = 1.
         stmirb(i,i) = dts
      enddo
c
      stmrb(1,2)   = dts
      stmrb(1,4)   = 0.
      stmrb(3,4)   = dts
      stmrb(2,4)   = 0.
      stmirb(1,2)  = 0.5*dts*dts
      stmirb(1,4)  = 0.
      stmirb(2,4)  = 0.
      stmirb(3,4)  = 0.5*dts*dts
      bmatrb(2,2)  = 1./tmass
      bmatrb(4,4)  = 1./yinert
c
c     overwrite stmi with matrix product stmi*bmat
c     (theta*B in the references's notation)
c
      do j=1,4
         do i=1,4
            temp(i,j) = 0.
            do k=1,4
               temp(i,j) = temp(i,j) + stmirb(i,k)*bmatrb(k,j)
            end do
         end do
      end do
      do j=1,4
         do i=1,4
            stmirb(i,j) = temp(i,j)
         end do
      end do
c
      return
      end
