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
      subroutine diagk(i,npl,jdim,kdim,idim,q,res,dtj,sk,t,vol,vist3d,
     .                 blank,iover)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve scalar tridiagonal equations to approximate the
c     spatially-split factor in the K-direction of the 3-d spatially-
c     split algorithm.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5),t(npl*(jdim-1)*kdim,35)
      dimension sk(jdim,kdim,idim-1,5),vol(jdim,kdim,idim-1)
      dimension q(jdim,kdim,idim,5),dtj(jdim,kdim,idim-1)
      dimension vist3d(jdim,kdim,idim),blank(jdim,kdim,idim)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
      common /entfix/ epsa_l,epsa_r
      common /zero/iexp
c
c     10.**(-iexp) is machine zero
      zero    = 10.**(-iexp)
      epsa_l  = 2.*epsa_r
c
c     k-implicit j-sweep line inversions af
c
      kdim1 = kdim-1
      jdim1 = jdim-1
      xmre  = 8.e0*xmach/reue
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
c     load rhs (delta q * dtj) into f
c
c      temporary variables are as below:  
c
c      t(1-4)  :  kx,ky,kz,cell-face area
c      t(5)    :  c = sound speed
c      t(6)    :  rho
c      t(7)    :  u
c      t(8)    :  v
c      t(9)    :  w
c      t(10)   :  p
c      t(11)   :  ubar = u*kx+v*ky+w*kz+kt
c      t(12)   :  cell volume
c      t(14-16):  lx,ly,lz  (unit vectors in plane)
c      t(17-19):  mx,my,mz  (unit vectors in plane)
c      t(21)   :  dtj = 1./( dt * j ) = vol / dt
c      t(22-24):  a,b,c tridiagonal coefficients
c      t(25)   :  diagonal viscous term
c      t(26-30):  rhs vectors
c      t(33,34):  preconditioned acoustic evals - u'+/-a'
c      t(35)   :  preconditioned reference Mach number-squared *
c                 sound speed
c
      jv  = npl*jdim1
      do 1009 ipl=1,npl
      ii  = i+ipl-1
      jv0 = (ipl-1)*jdim1+1
      do 1005 k=1,kdim1
      jk  = (k-1)*jv + jv0
      do 1004 l=1,5
cdir$ ivdep
      do 1000 izz=1,jdim1
      t(izz+jk-1,25+l) = res(izz,k,ii,l)
      t(izz+jk-1,l+5)  = q(izz,k,ii,l)
 1000 continue
 1004 continue
cdir$ ivdep
      do 1001 izz=1,jdim1
      t(izz+jk-1,16) = sk(izz,k,ii,1)
      t(izz+jk-1,17) = sk(izz,k,ii,2)
      t(izz+jk-1,18) = sk(izz,k,ii,3)
      t(izz+jk-1,19) = sk(izz,k,ii,4)
      t(izz+jk-1,20) = sk(izz,k,ii,5)
      t(izz+jk-1,21) = tfacp1*dtj(izz,k,ii)
 1001 continue
      if (ivisc(3).gt.1) then
cdir$ ivdep
         do 1002 izz=1,jdim1
         t(izz+jk-1,12) = vol(izz,k,ii)
         t(izz+jk-1,31) = vist3d(izz,k,ii)
 1002    continue
      else if(ivisc(3).gt.0) then
cdir$ ivdep
         do 1003 izz=1,jdim1
         t(izz+jk-1,12) = vol(izz,k,ii)
 1003    continue
      end if
 1005 continue
      jk = kdim1*jv+jv0
cdir$ ivdep
      do 1006 izz=1,jdim1
      t(izz+jk-1,16) = sk(izz,kdim,ii,1)
      t(izz+jk-1,17) = sk(izz,kdim,ii,2)
      t(izz+jk-1,18) = sk(izz,kdim,ii,3)
      t(izz+jk-1,19) = sk(izz,kdim,ii,4)
      t(izz+jk-1,20) = sk(izz,kdim,ii,5)
 1006 continue
 1009 continue
c
      n = jv*kdim1
      do 1010 l=1,5
cdir$ ivdep
      do 1007 izz=1,n
      t(izz,25+l) = t(izz,25+l)*t(izz,21)
 1007 continue
 1010 continue
c
c      average metric
c
cdir$ ivdep
      do 1008 izz=1,n
      t1       = t(izz,16)+t(izz+jv,16) 
      t2       = t(izz,17)+t(izz+jv,17) 
      t3       = t(izz,18)+t(izz+jv,18) 
      t4       = t1*t1+t2*t2+t3*t3
      t4       = 1.e0/sqrt(t4)
      t(izz,1) = t1*t4
      t(izz,2) = t2*t4
      t(izz,3) = t3*t4
      t(izz,13)= 0.5*(t(izz,20)+t(izz+jv,20))
 1008 continue
cdir$ ivdep
      do 1011 izz=1,n+jv
      t(izz,4) = 0.50*t(izz,19)
 1011 continue
c
c      recover primitives
c
cdir$ ivdep
      do 1012 izz=1,n
      t(izz,16) = 1.e0/t(izz,6)
 1012 continue
c
c      viscous term
c
      if (ivisc(3).gt.0) then
      if (ivisc(3).gt.1) then
cdir$ ivdep
         do 1013 izz=1,n
         t(izz,32) = t(izz,16)*(1.e0+t(izz,31))
 1013    continue
      else
cdir$ ivdep
         do 1014 izz=1,n
         t(izz,32) = t(izz,16)
 1014    continue
      end if
cdir$ ivdep
      do 1015 izz=1,n+jv
      t(izz,25) = xmre*t(izz,4)*t(izz,4)
 1015 continue
cdir$ ivdep
      do 1016 izz=1,jv
      t(izz,25) = t(izz,25)*t(izz,32)/t(izz,12)
 1016 continue
      ns = n-jv
cdir$ ivdep
      do 1017 izz=1,ns
      t(izz+jv,25) = t(izz+jv,25)*(t(izz,32)
     .              +t(izz+jv,32))/(t(izz,12)+t(izz+jv,12))
 1017 continue
cdir$ ivdep
      do 1018 izz=1,jv
      t(izz+n,25) = t(izz+n,25)*t(izz+ns,32)/t(izz+ns,12)
 1018 continue
      else
cdir$ ivdep
      do 1019 izz=1,n+jv
      t(izz,25) = 0.e0
 1019 continue
      end if
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1020 izz=1,n
         t(izz,5)  = sqrt(gamma*t(izz,10)*t(izz,16))
         t(izz,11) = t(izz,1)*t(izz,7)+t(izz,2)*t(izz,8)
     .             + t(izz,3)*t(izz,9)+t(izz,13)
 1020    continue
      else
cdir$ ivdep
         do 10201 izz=1,n
         t(izz,5)  = sqrt(gamma*t(izz,10)*t(izz,16))
         t(izz,11) = t(izz,1)*t(izz,7)+t(izz,2)*t(izz,8)
     .             + t(izz,3)*t(izz,9)+t(izz,13)
c
c -----  calculation of preconditioning quantities
c
         vmag1 =  t(izz,7)*t(izz,7) + t(izz,8)*t(izz,8)
     .         +  t(izz,9)*t(izz,9)
         vel2 = ccmax(vmag1,avn*uref**2)
         vel = sqrt(ccmin(t(izz,5)*t(izz,5),vel2))
         vel = cprec*vel + (1.-cprec)*t(izz,5)
         xm2 = (vel/t(izz,5))**2
         xmave = t(izz,11)/t(izz,5)
         t11 = 0.5*(1.+xm2)
         t21 = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
         t(izz,33) = t11*t(izz,11) + t21*t(izz,5)
         t(izz,34) = t11*t(izz,11) - t21*t(izz,5)
         t(izz,35) = xm2*t(izz,5)
10201    continue
      end if
c
c     t(inverse) r
c
      maxf = jv*kdim
      call tinvr(n,t(1,26),t(1,27),t(1,28),t(1,29),t(1,30),t(1,1),
     .             t(1,2), t(1,3), t(1,14),t(1,15),t(1,16),t(1,17),
     .             t(1,18),t(1,19),t(1,5), t(1,11),t(1,6), t(1,7),
     .             t(1,8), t(1,9), maxf,1, t(1,33),t(1,34),t(1,35)) 
c
c     assemble and solve decoupled matrix equations
c
      il   = 1
      iu   = kdim1
c
      epsi = 0.
cdir$ ivdep
      do 1021 izz=1,n
      t(izz,31) = t(izz,11)
      t(izz,32) = ccabs(t(izz,31))
c
c     limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
c
      if (real(epsa_l) .gt. 0.) then
         cc    = ccabs(t(izz,5))
         uu    = ccabs(t(izz,7))
         vv    = ccabs(t(izz,8))
         ww    = ccabs(t(izz,9))
         epsaa = epsa_l*(cc + uu + vv + ww)
         epsbb = 0.25/ccmax(epsaa,zero)
         epscc = 2.00*epsaa
         if (real(t(izz,32)).lt.real(epscc))
     .       t(izz,32) = t(izz,32)*t(izz,32)*epsbb + epsaa
      end if
c
      t(izz,24) = t(izz,31)+t(izz,32)
      t(izz,31) = t(izz,31)-t(izz,32)
      t(izz,23) = t(izz,21)+t(izz+jv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .           +t(izz+jv,25)+t(izz,25)
 1021 continue
cdir$ ivdep
      do 1022 izz=1,n-jv
      t(izz+jv,22) = -t(izz,24)*t(izz+jv,4)-t(izz+jv,25)
      t(izz,24)    =  t(izz+jv,31)*t(izz+jv,4)-t(izz+jv,25)
 1022 continue
c
      if (iover.eq.1)
     . call dabckz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      call dlutr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24))
      call dfbtr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,26))
      call dfbtr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,27))
      call dfbtr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,28))
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1023 izz=1,n
         t(izz,31) = t(izz,11)+t(izz,5)
         t(izz,32) = ccabs(t(izz,31))
c
c        limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_l) .gt. 0.) then
            cc    = ccabs(t(izz,5))
            uu    = ccabs(t(izz,7))
            vv    = ccabs(t(izz,8))
            ww    = ccabs(t(izz,9))
            epsaa = epsa_l*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t(izz,32)).lt.real(epscc))
     .          t(izz,32) = t(izz,32)*t(izz,32)*epsbb + epsaa
         end if
c
         t(izz,24) = t(izz,31)+t(izz,32)
         t(izz,31) = t(izz,31)-t(izz,32)
         t(izz,23) = t(izz,21)+t(izz+jv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+jv,25)+t(izz,25)
 1023    continue
      else
cdir$ ivdep
         do 10231 izz=1,n
         t(izz,31) = t(izz,33)
         t(izz,32) = ccabs(t(izz,31))
c
c        limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_l) .gt. 0.) then
            cc    = ccabs(t(izz,5))
            uu    = ccabs(t(izz,7))
            vv    = ccabs(t(izz,8))
            ww    = ccabs(t(izz,9))
            epsaa = epsa_l*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t(izz,32)).lt.real(epscc))
     .          t(izz,32) = t(izz,32)*t(izz,32)*epsbb + epsaa
         end if
c
         t(izz,24) = t(izz,31)+t(izz,32)
         t(izz,31) = t(izz,31)-t(izz,32)
         t(izz,23) = t(izz,21)+t(izz+jv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+jv,25)+t(izz,25)
10231    continue
      end if
cdir$ ivdep
      do 1024 izz=1,n-jv
      t(izz+jv,22) = -t(izz,24)*t(izz+jv,4)-t(izz+jv,25)
      t(izz,24)    =  t(izz+jv,31)*t(izz+jv,4)-t(izz+jv,25)
 1024 continue
c
      if (iover.eq.1)
     . call dabckz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      call dlutr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24))
      call dfbtr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,29))
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1025 izz=1,n
         t(izz,31) = t(izz,11)-t(izz,5)
         t(izz,32) = ccabs(t(izz,31))
c
c        limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_l) .gt. 0.) then
            cc    = ccabs(t(izz,5))
            uu    = ccabs(t(izz,7))
            vv    = ccabs(t(izz,8))
            ww    = ccabs(t(izz,9))
            epsaa = epsa_l*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t(izz,32)).lt.real(epscc))
     .          t(izz,32) = t(izz,32)*t(izz,32)*epsbb + epsaa
         end if
c
         t(izz,24) = t(izz,31)+t(izz,32)
         t(izz,31) = t(izz,31)-t(izz,32)
         t(izz,23) = t(izz,21)+t(izz+jv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+jv,25)+t(izz,25)
 1025    continue
      else
cdir$ ivdep
         do 10251 izz=1,n
         t(izz,31) = t(izz,34)
         t(izz,32) = ccabs(t(izz,31))
c
c        limit eigenvalue a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_l) .gt. 0.) then
            cc    = ccabs(t(izz,5))
            uu    = ccabs(t(izz,7))
            vv    = ccabs(t(izz,8))
            ww    = ccabs(t(izz,9))
            epsaa = epsa_l*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t(izz,32)).lt.real(epscc))
     .          t(izz,32) = t(izz,32)*t(izz,32)*epsbb + epsaa
         end if
c
         t(izz,24) = t(izz,31)+t(izz,32)
         t(izz,31) = t(izz,31)-t(izz,32)
         t(izz,23) = t(izz,21)+t(izz+jv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+jv,25)+t(izz,25)
10251    continue
      end if
cdir$ ivdep
      do 1026 izz=1,n-jv
      t(izz+jv,22) = -t(izz,24)*t(izz+jv,4)-t(izz+jv,25)
      t(izz,24)    =  t(izz+jv,31)*t(izz+jv,4)-t(izz+jv,25)
 1026 continue
c
      if (iover.eq.1)
     . call dabckz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      call dlutr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24))
      call dfbtr(jv,jv,kdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,30))
c
c      t * delta q
c
      call tdq  (n,t(1,26),t(1,27),t(1,28),t(1,29),t(1,30),t(1,1),
     .             t(1,2), t(1,3), t(1,14),t(1,15),t(1,16),t(1,17),
     .             t(1,18),t(1,19),t(1,5), t(1,11),t(1,6), t(1,7),
     .             t(1,8), t(1,9), maxf, t(1,33), t(1,34), t(1,35))
c
c     update delta q
c
      do 1300 ipl=1,npl
      ii  = i+ipl-1
      jv0 = (ipl-1)*jdim1 + 1
      do 1300 l=1,5
      do 1300 k=1,kdim1
      jk  = jv0+(k-1)*jv
cdir$ ivdep
      do 1027 izz=1,jdim1
      res(izz,k,ii,l) = t(izz+jk-1,25+l)
 1027 continue
 1300 continue
      return
      end
