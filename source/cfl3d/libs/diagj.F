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
      subroutine diagj(i,npl,jdim,kdim,idim,q,res,dtj,sj,t,iperd,
     .                 vol,vist3d,blank,iover)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve scalar tridiagonal equations to approximate the
c     spatially-split factor in the J-direction of the 3-d spatially-
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
      dimension res(jdim,kdim,idim-1,5),t(npl*kdim*jdim,35)
      dimension sj(jdim,kdim,idim-1,5),vol(jdim,kdim,idim-1)
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
c     j-implicit k-sweep line inversions af
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      xmre  = 8.e0*xmach/reue
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
c     load rhs (-residual) into f
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
c      t(14-16):  lx,ly,lz  (unit vectors in plane)
c      t(17-19):  mx,my,mz  (unit vectors in plane)
c      t(26-30):  rhs vectors
c      t(21)   :  dtj = 1./( dt * j ) = vol / dt
c      t(22-24):  a,b,c tridiagonal coefficients
c      t(33,34):  preconditioned acoustic evals - u'+/-a'
c      t(35)   :  preconditioned reference Mach number-squared *
c                 sound speed
c
      kv = npl*kdim
      do 1009 j=1,jdim1
      kj = (j-1)*kv+1
      do 1004 l=1,5
c
      jj = 1-jdim
      do 8466 ii=1,kv
      jj = jj+jdim
 8466 t(kj+ii-1,25+l) = -res(j+jj-1,1,i,l)
c      call q8vgathp(kv,res(j,1,i,l),jdim,kv,kv,t(kj,25+l))
c
      jj = 1-jdim
      do 8467 ii=1,kv
      jj = jj+jdim
 8467 t(kj+ii-1,l+5) = q(j+jj-1,1,i,l)
c      call q8vgathp(kv,q(j,1,i,l),jdim,kv,kv,t(kj,l+15))
 1004 continue
c
      jj = 1-jdim
      do 8458 ii=1,kv
      jj = jj+jdim
 8458 t(kj+ii-1,21) = tfacp1*dtj(j+jj-1,1,i)
c      call q8vgathp(kv,dtj(j,1,i),jdim,kv,kv,t(kj,21))
      if(ivisc(2) .gt. 0) then
        jj=1-jdim
        do 9458 ii=1,kv
        jj=jj+jdim
 9458   t(kj+ii-1,12)=vol(j+jj-1,1,i)
      end if
 1009 continue
c
      do 1119 j=1,jdim
      kj = (j-1)*kv+1
      do 1119 l=1,5
c
      jj = 1-jdim
      do 8459 ii=1,kv
      jj = jj+jdim
 8459 t(kj+ii-1,15+l) = sj(j+jj-1,1,i,l)
c      call q8vgathp(kv,sj(j,1,i,l),jdim,kv,kv,t(kj,5+l))
 1119 continue
      if(ivisc(2) .gt. 1) then
        ic=0
        do 8558 ipl=1,npl
          ii=i+ipl-1
          do 8558 k=1,kdim
            ic=ic+1
            if(k .ne. kdim) then
              do 1118 j=1,jdim1
                kj=(j-1)*kv
                t(kj+ic,31)=vist3d(j,k,ii)
 1118         continue
            else
              do 1120 j=1,jdim1
                kj=(j-1)*kv
                t(kj+ic,31)=vist3d(j,kdim1,ii)
 1120         continue
            end if
 8558   continue
      end if
c
      n = kv*jdim1
c
c      average metric
c
cdir$ ivdep
      do 1001 izz=1,n
      t1       = t(izz,16)+t(izz+kv,16) 
      t2       = t(izz,17)+t(izz+kv,17) 
      t3       = t(izz,18)+t(izz+kv,18) 
      t4       = t1*t1+t2*t2+t3*t3
      t4       = 1.e0/sqrt(t4)
      t(izz,1) = t1*t4
      t(izz,2) = t2*t4
      t(izz,3) = t3*t4
      t(izz,13)= 0.5*(t(izz,20)+t(izz+kv,20))
 1001 continue
cdir$ ivdep
      do 1002 izz=1,n+kv
      t(izz,4) = 0.50*t(izz,19)
 1002 continue
c
c     recover primitives
c
c   viscous term
      if(ivisc(2) .gt. 0) then
      if(ivisc(2) .gt. 1) then
cdir$ ivdep
      do 7013 izz=1,n
        t(izz,32)=(1.e0+t(izz,31))/t(izz,6)
 7013 continue
      else
cdir$ ivdep
      do 7014 izz=1,n
        t(izz,32)=1.e0/t(izz,6)
 7014 continue
      end if
cdir$ ivdep
      do 7015 izz=1,n+kv
        t(izz,25)=xmre*t(izz,4)*t(izz,4)
 7015 continue
cdir$ ivdep
      do 7016 izz=1,kv
        t(izz,25)=t(izz,25)*t(izz,32)/t(izz,12)
 7016 continue
      ns=n-kv
cdir$ ivdep
      do 7017 izz=1,ns
        t(izz+kv,25)=t(izz+kv,25)*(t(izz,32)
     .              +t(izz+kv,32))/(t(izz,12)+t(izz+kv,12))
 7017 continue
cdir$ ivdep
      do 7018 izz=1,kv
        t(izz+n,25)=t(izz+n,25)*t(izz+ns,32)/t(izz+ns,12)
 7018 continue
      else
cdir$ ivdep
      do 7019 izz=1,n+kv
        t(izz,25)=0.e0
 7019 continue
      end if
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1003 izz=1,n
         t(izz,5)  = sqrt(gamma*t(izz,10)/t(izz,6))
         t(izz,11) = t(izz,1)*t(izz,7)+t(izz,2)*t(izz,8)
     .             + t(izz,3)*t(izz,9)+t(izz,13)
 1003    continue
      else
cdir$ ivdep
         do 10031 izz=1,n
         t(izz,5)  = sqrt(gamma*t(izz,10)/t(izz,6))
         t(izz,11) = t(izz,1)*t(izz,7)+t(izz,2)*t(izz,8)
     .             + t(izz,3)*t(izz,9)+t(izz,13)
c
c -----  calculation of preconditioning quantities
c
         vmag1 =  t(izz,7)*t(izz,7) + t(izz,8)*t(izz,8)
     .          + t(izz,9)*t(izz,9)
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
10031    continue
      end if
c
c     t(inverse) r
c
      maxf = kv*jdim
      call tinvr(n,t(1,26),t(1,27),t(1,28),t(1,29),t(1,30),t(1,1),
     .             t(1,2), t(1,3), t(1,14),t(1,15),t(1,16),t(1,17),
     .             t(1,18),t(1,19),t(1,5), t(1,11),t(1,6), t(1,7),
     .             t(1,8), t(1,9), maxf,0, t(1,33),t(1,34),t(1,35))
c
c     assemble and solve decoupled matrix equations
c
      il   = 1
      iu   = jdim1
c
      epsi = 0.
cdir$ ivdep
      do 1005 izz=1,n
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
      t(izz,23) = t(izz,21)+t(izz+kv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .           +t(izz+kv,25)+t(izz,25)
 1005 continue
      if (iperd.eq.1) then
cdir$ ivdep
         do 1006 izz=1,kv
         t(izz,22)      = -t(izz+n-kv,24)*t(izz,4) - t(izz,25)
         t(izz+n-kv,24) = t(izz,31)*t(izz+n,4) - t(izz+n,25)
 1006    continue
      end if
cdir$ ivdep
      do 1007 izz=1,n-kv
      t(izz+kv,22) = -t(izz,24)*t(izz+kv,4) - t(izz+kv,25)
      t(izz,24)    =  t(izz+kv,31)*t(izz+kv,4) - t(izz+kv,25)
 1007 continue
c
      if (iover.eq.1)
     . call dabcjz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      if (iperd.eq.0) then
         call dlutr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24))
         call dfbtr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,26))
         call dfbtr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,27))
         call dfbtr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,28))
      else
         call dlutrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),
     .               t(1,31),t(1,32))
         call dfbtrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,26),
     .               t(1,31),t(1,32))
         call dfbtrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,27),
     .               t(1,31),t(1,32))
         call dfbtrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,28),
     .               t(1,31),t(1,32))
      end if
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1008 izz=1,n
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
         t(izz,23) = t(izz,21)+t(izz+kv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+kv,25)+t(izz,25)
 1008    continue
      else
cdir$ ivdep
         do 10081 izz=1,n
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
         t(izz,23) = t(izz,21)+t(izz+kv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+kv,25)+t(izz,25)
10081    continue
      end if
      if (iperd.eq.1) then
cdir$ ivdep
         do 1011 izz=1,kv
         t(izz,22)      = -t(izz+n-kv,24)*t(izz,4) - t(izz,25)
         t(izz+n-kv,24) =  t(izz,31)*t(izz+n,4) - t(izz+n,25)
 1011    continue
      end if
cdir$ ivdep
      do 1012 izz=1,n-kv
      t(izz+kv,22) = -t(izz,24)*t(izz+kv,4) - t(izz+kv,25)
      t(izz,24)    =  t(izz+kv,31)*t(izz+kv,4) - t(izz+kv,25)
 1012 continue
c
      if (iover.eq.1)
     . call dabcjz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      if (iperd.eq.0) then
         call dlutr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24))
         call dfbtr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,29))
      else
         call dlutrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),
     .               t(1,31),t(1,32))
         call dfbtrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,29),
     .               t(1,31),t(1,32))
      end if
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1013 izz=1,n
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
         t(izz,23) = t(izz,21)+t(izz+kv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+kv,25)+t(izz,25)
 1013    continue
      else
cdir$ ivdep
         do 10131 izz=1,n
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
         t(izz,23) = t(izz,21)+t(izz+kv,4)*t(izz,24)-t(izz,4)*t(izz,31)
     .              +t(izz+kv,25)+t(izz,25)
10131    continue
      end if
      if (iperd.eq.1) then
cdir$ ivdep
         do 1014 izz=1,kv
         t(izz,22)      = -t(izz+n-kv,24)*t(izz,4) - t(izz,25)
         t(izz+n-kv,24) =  t(izz,31)*t(izz+n,4) - t(izz+n,25)
 1014    continue
      end if
cdir$ ivdep
      do 1015 izz=1,n-kv
      t(izz+kv,22) = -t(izz,24)*t(izz+kv,4) - t(izz+kv,25)
      t(izz,24)    =  t(izz+kv,31)*t(izz+kv,4) - t(izz+kv,25)
 1015 continue
c
      if (iover.eq.1)
     . call dabcjz(i,npl,jdim,kdim,idim,t(1,22),t(1,23),t(1,24),blank)
c
      if (iperd.eq.0) then
         call dlutr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24))
         call dfbtr(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,30))
      else
         call dlutrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),
     .               t(1,31),t(1,32))
         call dfbtrp(kv,kv,jdim,il,iu,t(1,22),t(1,23),t(1,24),t(1,30),
     .               t(1,31),t(1,32))
      end if
c
c      t * delta q
c
      call tdq  (n,t(1,26),t(1,27),t(1,28),t(1,29),t(1,30),t(1,1),
     .             t(1,2), t(1,3), t(1,14),t(1,15),t(1,16),t(1,17),
     .             t(1,18),t(1,19),t(1,5), t(1,11),t(1,6), t(1,7),
     .             t(1,8), t(1,9), maxf, t(1,33),t(1,34),t(1,35))
c
c     update delta q
c
      do 1300 j=1,jdim1
      kj = (j-1)*kv+1
      do 1300 l=1,5
c
      jj = 1-jdim
      do 8445 ii=1,kv
      jj = jj+jdim
 8445 res(j+jj-1,1,i,l) = t(kj+ii-1,25+l)
c      call q8vscatp(kv,t(kj,25+l),jdim,kv,kv,res(j,1,i,l))
 1300 continue
c
      do 1301 ipl=1,npl
      ii = i+ipl-1
      do 1301 l=1,5
cdir$ ivdep
      do 1016 izz=1,jdim1
      res(izz,kdim,ii,l) = 0.e0
 1016 continue
 1301 continue
      return
      end
