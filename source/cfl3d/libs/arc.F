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
      subroutine arc(jdim1,kdim1,msub1,msub2,
     .               jjmax1,kkmax1,lmax1,x1,y1,z1,
     .               limit0,jjmax2,kkmax2,
     .               x2,y2,z2,xie2,eta2,mblkpt,
     .               ifit,itmax,jcorr,kcorr,
     .               sxie,seta,sxie2,seta2,xie2s,eta2s,
     .               intmx,icheck,nblkj,nblkk,jmm,kmm,j21,j22,k21,k22,
     .               npt,xif1,xif2,etf1,etf2,nou,bou,nbuf,
     .               ibufdim,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform arc-length correction to generalized coordinates
c     near boundary if required when shearing correction has failed.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x1(jdim1,kdim1,msub1),y1(jdim1,kdim1,msub1),
     .          z1(jdim1,kdim1,msub1)
      dimension x2(jdim1,kdim1,msub2),y2(jdim1,kdim1,msub2),
     .          z2(jdim1,kdim1,msub2)
      dimension xie2(npt),eta2(npt),mblkpt(npt)
      dimension mblk2nd(maxbl)
      dimension jjmax1(msub1),kkmax1(msub1),jjmax2(msub2),kkmax2(msub2)
      dimension sxie(jdim1,kdim1,msub1),seta(jdim1,kdim1,msub1)
      dimension sxie2(jdim1,kdim1,msub2),seta2(jdim1,kdim1,msub2)
      dimension xie2s(jdim1,kdim1),eta2s(jdim1,kdim1),nblkj(jdim1),
     .          nblkk(kdim1),jmm(kdim1),kmm(jdim1)
c
      integer   xif1(msub1),xif2(msub1),etf1(msub1),etf2(msub1)
c
      l2    = 1
      jmax2 = jjmax2(l2)
      kmax2 = kkmax2(l2)
c
c     arc lengths for (unexpanded) "from" grid(s)
c
      do 601 l=1,lmax1
      jmax1 = jjmax1(l)
      kmax1 = kkmax1(l) 
      do 603 j=2,jmax1-1
      seta(j-1,1,l) = 0.
      do 604 k=3,kmax1-1
      seta(j-1,k-1,l) = seta(j-1,k-2,l)+sqrt((x1(j,k,l)-x1(j,k-1,l))**2+
     .                                       (y1(j,k,l)-y1(j,k-1,l))**2+
     .                                       (z1(j,k,l)-z1(j,k-1,l))**2)
  604 continue
  603 continue
      do 503 k=2,kmax1-1
      sxie(1,k-1,l) = 0.
      do 504 j=3,jmax1-1
      sxie(j-1,k-1,l) = sxie(j-2,k-1,l)+sqrt((x1(j,k,l)-x1(j-1,k,l))**2+
     .                                       (y1(j,k,l)-y1(j-1,k,l))**2+
     .                                       (z1(j,k,l)-z1(j-1,k,l))**2)
  504 continue
  503 continue
  601 continue
c
c     arc lengths for "to" grid
c
      do 613 j=1,jmax2
      seta2(j,1,1) = 0.
      do 614 k=2,kmax2
      seta2(j,k,1) = seta2(j,k-1,1)+sqrt((x2(j,k,1)-x2(j,k-1,1))**2+
     .                                   (y2(j,k,1)-y2(j,k-1,1))**2+
     .                                   (z2(j,k,1)-z2(j,k-1,1))**2)
  614 continue
  613 continue
      do 513 k=1,kmax2
      sxie2(1,k,1) = 0.
      do 514 j=2,jmax2
      sxie2(j,k,1) = sxie2(j-1,k,1)+sqrt((x2(j,k,1)-x2(j-1,k,1))**2+
     .                                   (y2(j,k,1)-y2(j-1,k,1))**2+
     .                                   (z2(j,k,1)-z2(j-1,k,1))**2)
  514 continue
  513 continue
c
      if (kcorr.eq.1) then
c
c     try arc length correction near eta=0
c
c     Arc length in eta direction now replaces the 
c     physical coordinates x,y,z from which computational
c     coordinates away from the boundaries are to be found.
c     The current xie values are assumed correct.
c
c
      do 6901 j=j21,j22-1
      km   = kmm(j)
      if (km.lt.k21) go to 6901
c
      test1 = .5
      test2 = .75
      kup  = test1*kmm(j)
      kup2 = test2*kmm(j)
      if (kup2.le.k21) go to 6901
c
      do 6900 k=k21,kup2
c
c     center of "to" grid cell j,k
c
      sc = .25*(seta2(j,k,1)+seta2(j+1,k,1)+seta2(j+1,k+1,1)
     .         +seta2(j,k+1,1))
c
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      jp   = xie2(ll)
      xie  = xie2(ll)-jp
      lp   = mblkpt(ll)
      jmax = jjmax1(lp)
      kmax = kkmax1(lp)
      if (k.eq.k21) then
         kp = 1
      else 
         kp = eta2s(j,k-1)+1
      end if
      eta   = kp
      lsrch = 2
      call topol2(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,xie,eta,seta,
     .            limit0,sc,jp,kp,lp,lsrch,itmax,xiet,etat,xif1,xif2,
     .            etf1,etf2,nou,bou,nbuf,ibufdim,myid,mblk2nd,
     .            maxbl)
      eta2s(j,k) = etat
 6900 continue
c
c     Blend arc length corrected and uncorrected computational 
c     coordinates of "to" grid...arc length correction for k .le. kup,
c     uncorrected for k .ge. kup2, and linear interpolation of corrected
c     and uncorrected values in between.
c
      do 555 k=k21,kup
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      eta2(ll) = eta2s(j,k)
  555 continue
      if(kup2 .gt. kup) then
        do 556 k=kup,kup2
        phi1      = 1.-float(k-kup)/float(kup2-kup)
        phi2      = 1.-phi1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        eta2(ll) = eta2(ll)*phi2+eta2s(j,k)*phi1
  556   continue
      end if
c
c     check for monotonicity of corrected coordinates
c     
      do 557 k=k21+1,km
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      ll1 = (j22-j21)*(k-1-k21) + (j-j21+1)
      if (real(eta2(ll)).le.real(eta2(ll1))) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),'(''  WARNING: corrected eta values are'',
     .   '' not monotonic at j,k = '',i5,'','',i5)') j,k
      end if
  557 continue
 6901 continue
      end if
c
      if (jcorr.eq.1) then
c
c     try arc length correction near xie=0
c
c     Arc length in xie direction now replaces the 
c     physical coordinates x,y,z from which computational
c     coordinates away from the boundaries are to be found.
c     The current eta values are assumed correct.
c
c
      do 7901 k=k21,k22-1
      jm = jmm(k)
      if (jm.lt.j21) go to 7901
c
      test1 = .5
      test2 = .75
      jup  = test1*jmm(k)
      jup2 = test2*jmm(k)
      if (jup2.le.j21) go to 7901
c
      do 7900 j=j21,jup2
c
c     center of "to" grid cell j,k
c
      sc = .25*(sxie2(j,k,1)+sxie2(j+1,k,1)+sxie2(j+1,k+1,1)
     .         +sxie2(j,k+1,1))
c
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      kp   = eta2(ll)
      eta  = eta2(ll)-kp
      lp   = mblkpt(ll)
      jmax = jjmax1(lp)
      kmax = kkmax1(lp)
      if (j.eq.j21) then
         jp = 1
      else 
         jp = xie2s(j-1,k)+1
      end if
      xie   = jp
      lsrch = 1
      call topol2(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,xie,eta,sxie,
     .            limit0,sc,jp,kp,lp,lsrch,itmax,xiet,etat,xif1,xif2,
     .            etf1,etf2,nou,bou,nbuf,ibufdim,myid,mblk2nd,
     .            maxbl)
      xie2s(j,k) = xiet
 7900 continue
c
c     Blend arc length corrected and uncorrected computational 
c     coordinates of "to" grid...arc length correction for j .le. jup,
c     uncorrected for j .ge. jup2, and linear interpolation of corrected
c     and uncorrected values in between.
c
      do 2555 j=j21,jup
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      xie2(ll) = xie2s(j,k)
 2555 continue
      do 2556 j=jup,jup2
      phi1      = 1.-float(j-jup)/float(jup2-jup)
      phi2      = 1.-phi1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      xie2(ll) = xie2(ll)*phi2+xie2s(j,k)*phi1
 2556 continue
c
c     check for monotonicity of corrected coordinates
c     
      do 3557 j=j21+1,jm
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      ll1 = (j22-j21)*(k-k21) + (j-1-j21+1)
      if (real(xie2(ll)).le.real(xie2(ll1))) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),'(''  WARNING: corrected xie values are'',
     .   '' not monotonic at j,k = '',i5,'','',i5)') j,k
      end if
 3557 continue
 7901 continue
      end if
      return
      end
