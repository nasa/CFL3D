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
      subroutine topol2(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,xie,eta,
     .                  sarc,limit0,sc,jp,kp,lp,lsrch,itmax,xiet,etat,
     .                  xif1,xif2,etf1,etf2,nou,bou,nbuf,
     .                  ibufdim,myid,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  One dimensional equivalent to subroutine topol. Performs
c     search in one direction only.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension sarc(jdim1,kdim1,msub1),jjmax1(msub1),kkmax1(msub1)
      dimension mblk2nd(maxbl)
      integer   xif1(msub1),xif2(msub1),etf1(msub1),etf2(msub1)
c
      common /tol/ epsc,epsc0,epsreen,epscoll
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
c      
c     lsrch=1...search in xie direction, eta value is given
c     lsrch=2...search in eta direction, xie value is given
c
      jmax   = jjmax1(lp)
      kmax   = kkmax1(lp)
      js     = xif1(lp)
      je     = xif2(lp)
      ks     = etf1(lp)
      ke     = etf2(lp)
      intern = 0
      limit  = limit0
 6901 continue
      intern = intern+1
      if (intern.gt.itmax) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' program terminated in dynamic'',
     .         '' patching routines - see file '',a60)') grdmov
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)' stopping...intern > itmax '
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     bi-linear fit in target cell jp,kp in "from" block lp
c
      dy2 = sarc(jp+1,kp,lp)-sarc(jp,kp,lp)
      dy3 = sarc(jp+1,kp+1,lp)-sarc(jp,kp,lp)
      dy4 = sarc(jp,kp+1,lp)-sarc(jp,kp,lp)
      a1  = sarc(jp,kp,lp)
      a2  = dy2
      a3  = dy4
      a4  = dy3-a2-a3
c
c     for fixed eta
c
      if (lsrch.eq.1) then
      xie = (sc-a1-a3*eta)/(a2+a4*eta)
c
c     check to ensure point is inside cell
c
      imiss = 0
      if (real(xie).lt.-real(epsc) .or. real(xie).gt.1.+real(epsc)) then
         imiss = 1
      end if
c
      if (imiss.eq.0) go to 6902
c
c     try new target cell
c
      jpc = jp
      if (real(xie).ge.0) jinc = abs(xie)
      if (real(xie).lt.0) jinc = abs(xie-1)
      if (limit.gt.3.and.intern.gt.5)  limit=3
      if (limit.gt.1.and.intern.gt.10) limit=1
      jinc = min(jinc,limit)
      if (real(xie).gt.1.0) then
         jpc = jp + jinc
      else if (real(xie).lt.0.) then
         jpc = jp - jinc 
      end if
      jp = min(jpc,je-1)
      jp = max(js,jp)
      go to 6901
c
c     for fixed xie
c
      else if (lsrch.eq.2) then
      eta = (sc-a1-a2*xie)/(a3+a4*xie)
c
      imiss = 0
      if (real(eta).lt.-real(epsc) .or. real(eta).gt.1.+real(epsc)) then
         imiss = 1
      end if
      if (imiss.eq.0) go to 6902
c
c     try new target cell
c
      kpc = kp
      kpl = kp
      if (real(eta).ge.0) kinc = abs(eta)
      if (real(eta).lt.0) kinc = abs(eta-1)
      if (limit.gt.3 .and. intern.gt.5)  limit = 3
      if (limit.gt.1 .and. intern.gt.10) limit = 1
      kinc = min( kinc , limit )
      if (real(eta).gt.1.0) then
         kpc = kp + kinc
      else if (real(eta).lt.0.) then
         kpc = kp - kinc 
      end if
      kp = min(kpc,ke-1)
      kp = max(ks,kp)
      go to 6901
      end if
 6902 continue
      xiet = jp+xie
      etat = kp+eta
      return
      end
