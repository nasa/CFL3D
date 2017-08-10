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
      subroutine gfluxr(i,npl,xkap,jdim,kdim,idim,res,q,qj0,sj,t,nvtq,
     .                  nv,nfa,wfa,iwfa,jbctyp,isf,nbl,bcj,nou,bou,nbuf,
     .                  ibufdim,myid,mblk2nd,maxbl,idef)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute residual contributions for the 
c     right-hand-side in the J-direction from the inviscid terms.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension wfa(1),iwfa(maxbl*7*3)
      dimension jbctyp(2)
      dimension sj(jdim*kdim,idim-1,5)
      dimension q(jdim,kdim,idim,5),qj0(kdim,idim-1,5,4)
      dimension res(jdim,kdim,idim-1,5),bcj(kdim,idim-1,2)
      dimension t(nvtq,nv),mblk2nd(maxbl)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /cpurate/ rate(5),ratesub(5),ncell(20)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /twod/ i2d
      common /chk/ ichk
      common /fvfds/ rkap0(3),ifds(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /sklton/ isklton
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      kg    = kdim*npl
      n     = jdim*kdim*npl
      nr    = n-1
c
      if (real(xkap).lt.-2.e0) then
c
c      first order
c
c     fill temp arrays with interior values
c
      do 200 l=1,5
cdir$ ivdep
      do 50 izz=1,nr
      t(izz+1,20+l) = q(izz,1,i,l)
      t(izz,25+l)   = q(izz,1,i,l)
   50 continue
c
c     fill temp arrays with boundary values
c
      jj  = 1-jdim
cdir$ ivdep
      do 100 kk=1,kg 
      jj  = jj+jdim
      jj2 = jj+jdim1
      t(jj,20+l)  = qj0(kk,i,l,1)
      t(jj2,25+l) = qj0(kk,i,l,3)
  100 continue
  200 continue
c
      else
c
c      higher order
c
      do 1000 l=1,5
c
c     calculate interior gradients
c
cdir$ ivdep
      do 300 izz=1,nr
      t(izz+1,1) = q(izz+1,1,i,l)-q(izz,1,i,l)
  300 continue
c
c     edge gradients - left boundary
c
      do 400 ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*kdim*jdim+1-jdim
      do 400 kk=1,kdim1
      jc      = jc + jdim
      t(jc,1) = (1.0-bcj(kk,ii,1))*(q(1,kk,ii,l)-qj0(kk,ii,l,1))
     .              +bcj(kk,ii,1) * qj0(kk,ii,l,2)
  400 continue
c
c      edge gradients - right boundary
c
      do 500 ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*kdim*jdim
      do 500 kk=1,kdim1
      jc      = jc+jdim
      t(jc,1) = (1.0-bcj(kk,ii,2))*(qj0(kk,ii,l,3)-q(jdim1,kk,ii,l))
     .              +bcj(kk,ii,2) * qj0(kk,ii,l,4)
  500 continue
c
c     zero out gradients of dummy interfaces at k=kdim
c
      do 600 ipl=1,npl
      kvl = (ipl-1)*jdim*kdim+kdim1*jdim
      do 600 j=1,jdim
      t(j+kvl,1) = 0.0
  600 continue
c 
c      gradient limiting - cell interface interpolations
c
cdir$ ivdep
      do 700 izz=1,nr
      t(izz,2) = t(izz+1,1)
  700 continue
c
      ifl = iflim(2)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = jdim1
      end if
      call xlim(xkap,nr,t(1,1),t(1,2),q(1,1,i,l),ifl,ncells,l)
c
cdir$ ivdep
      do 800 izz=1,nr
      t(izz+1,20+l) = q(izz,1,i,l)+t(izz,2)
      t(izz  ,25+l) = q(izz,1,i,l)-t(izz,1)
  800 continue 
 1000 continue
      end if
c
      do 2000 l=1,5
c
c      edge values - left boundary
c
      do 1300 ipl=1,npl
      ii = i+ipl-1
      kl = 0
      do 1100 kk=1,kdim1
      kl = kl + 1
      t(kl,1) = qj0(kk,ii,l,1) - qj0(kk,ii,l,2)
      t(kl,2) = q(1,kk,ii,l)   - qj0(kk,ii,l,1)
      t(kl,3) = qj0(kk,ii,l,1)
 1100 continue
      ifl = iflim(2)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = jdim1
      end if
      call xlim(xkap,kdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = (ipl-1)*(kdim*jdim) + 1 - jdim
      kl = 0
      do 1200 kk=1,kdim1
      kl = kl + 1
      jc = jc + jdim
      t(jc,20+l) = (1.0-bcj(kk,ii,1))*(qj0(kk,ii,l,1)+t(kl,2))
     .                 +bcj(kk,ii,1) * qj0(kk,ii,l,1)
      t(jc,25+l) = (1.0-bcj(kk,ii,1))* t(jc,25+l)
     .                 +bcj(kk,ii,1) * qj0(kk,ii,l,1)
 1200 continue
 1300 continue
c
c      edge values - right boundary
c
      do 1800 ipl=1,npl
      ii = i+ipl-1
      kl = 0
      do 1600 kk=1,kdim1
      kl = kl + 1
      t(kl,2) = qj0(kk,ii,l,4) - qj0(kk,ii,l,3)
      t(kl,1) = qj0(kk,ii,l,3) - q(jdim1,kk,ii,l)
      t(kl,3) = qj0(kk,ii,l,3)
 1600 continue
      ifl = iflim(2)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = jdim1
      end if
      call xlim(xkap,kdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = (ipl-1)*(kdim*jdim)
      kl = 0
      do 1700 kk=1,kdim1
      kl = kl + 1
      jc = jc + jdim
      t(jc,20+l) = (1.0-bcj(kk,ii,2))* t(jc,20+l)
     .                 +bcj(kk,ii,2) * qj0(kk,ii,l,3)
      t(jc,25+l) = (1.0-bcj(kk,ii,2))*(qj0(kk,ii,l,3)-t(kl,1))
     .                 +bcj(kk,ii,2) * qj0(kk,ii,l,3)
 1700 continue
 1800 continue
c
c     fill end point for safety
      t(n,25+l) = t(nr,25+l)
c
 2000 continue
c
      if (ichk.eq.1) then
         epsz = 1.0e-03
         epss = 1.0e+03
         do 5432 ipl=1,npl
         ii = i+ipl-1
         do 5432 k=1,kdim
         do 5432 j=1,jdim
         ic = jdim*kdim*(ipl-1) + (k-1)*jdim + j
         if (real(t(ic,21)).lt.real(epsz) .or. 
     .       real(t(ic,25)).lt.real(epsz) .or.
     .       real(t(ic,21)).gt.real(epss) .or. 
     .       real(t(ic,25)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in gflux left - small ',
     .                 '(or large) density and/or pressure at '
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' j,k,i,t(21),t(25) = ',j,k,ii,
     .      real(t(ic,21)),real(t(ic,25))
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (real(t(ic,26)).lt.real(epsz) .or.
     .       real(t(ic,30)).lt.real(epsz) .or.
     .       real(t(ic,26)).gt.real(epss) .or.
     .       real(t(ic,30)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in gflux right - small ',
     .                 ' (or large) density and/or pressure at '
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' j,k,i,t(26),t(30) = ',j,k,ii,
     .      real(t(ic,26)),real(t(ic,30))
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
 5432    continue
      end if
c
      jkpro = jdim*kdim*npl-jdim
      if (ifds(2).eq.0) then
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),185) nbl
      end if
  185 format(49h   computing inviscid fluxes, J-direction - flux-,
     .24hvector splitting - block,i4)
c
      call fluxp(sj(1,i,1),sj(1,i,2),sj(1,i,3),sj(1,i,4),sj(1,i,5),
     .           t(1,21),t(1,31),jkpro,t,jkpro,nvtq,nou,bou,nbuf,
     .           ibufdim)
c
      do 1150 l=1,5
cdir$ ivdep
      do 1011 izz=1,jkpro
      t(izz,20+l) = t(izz,30+l)
 1011 continue
 1150 continue
c
      call fluxm(sj(1,i,1),sj(1,i,2),sj(1,i,3),sj(1,i,4),sj(1,i,5),
     .           t(1,26),t(1,31),jkpro,t,jkpro,nvtq,nou,bou,nbuf,
     .           ibufdim)
c
      do 1400 l=1,5
cdir$ ivdep
      do 1013 izz=1,jkpro
      t(izz,30+l) = t(izz,30+l) + t(izz,20+l)
 1013 continue
 1400 continue
c
      else if (ifds(2).eq.1) then
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),184)
      end if
  184 format(49h   computing inviscid fluxes, J-direction - flux-,
     .20hdifference splitting)
c
      call fhat(sj(1,i,1),sj(1,i,2),sj(1,i,3),sj(1,i,4),sj(1,i,5),
     .          t(1,31),t(1,26),t(1,21),jkpro,nvtq) 
c
      else
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),183)
      end if
  183 format(49h   computing inviscid fluxes, J-direction - MAPS+,
     .15h flux splitting)
c
      call fmaps(sj(1,i,1),sj(1,i,2),sj(1,i,3),sj(1,i,4),sj(1,i,5),
     .          t(1,31),t(1,26),t(1,21),jkpro,nvtq)
c
      end if
c
c      conservative at embedded boundaries - wfa array contains fluxes
c      from a finer grid
c
      if (nfa.gt.0) then
      iis  = i
      iie  = i+npl-1
      lfcc = iwfa(7)
      do 7055 ifa=1,nfa
      ic   = (ifa-1)*7
      jsb  = iwfa(ic+1)
      ksb  = iwfa(ic+2)
      isb  = iwfa(ic+3)
      jeb  = iwfa(ic+4)
      keb  = iwfa(ic+5)
      ieb  = iwfa(ic+6)
      ifts = iwfa(ic+7)
c
      if (ifts.ne.lfcc) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*) ' inconsistent summations - stopping'
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      if (isklton.gt.0)then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),269) jsb,ksb,keb,isb,ieb
      end if
  269 format(2x,33h installing accumulated fluxes on,
     .3h j=,i3,16h at ks,ke,is,ie=,4i4)
c
c     loop over all planes in I-direction
c
      lfcc = ifts
      do 754 ii=1,idim1
c
c     skip planes not in embedded region
c
      if (ii.ge.isb .and. ii.lt.ieb) then
c
c     check for ii in region to be updated on this pass (npl planes of data)
c
      if (iis.le.ii .and. ii.le.iie) then
         ipl = ii-i+1
         do 705 l=1,5
         lfc = ifts+(l-1)*(keb-ksb)*(ieb-isb)+(ii-isb)*(keb-ksb)
         loc = (ipl-1)*jdim*kdim+(ksb-1)*jdim+jsb
         do 705 k=ksb,keb-1
         t(loc,30+l) = wfa(lfc)
         loc  = loc+jdim
         lfc  = lfc+1
         lfcc = lfcc+1
  705    continue
      else
c
c     increment counter for planes not updated on this pass
c
      lfcc = lfcc+5*(keb-ksb)
      end if
      end if
  754 continue
 7055 continue
      end if
c
      do 450 l=1,5
cdir$ ivdep
      do 1014 izz=1,nr
      res(izz,1,i,l) = t(izz+1,30+l)-t(izz,30+l)
 1014 continue
  450 continue
c
c     geometric conservation law terms for deforming grids
c
      if (idef.gt.0) then
         oogmo = 1./(gamma-1.)
cdir$ ivdep
         do 1016 izz=2,nr
         t(izz,41) = -(sj(izz+1,i,5)*sj(izz+1,i,4)
     .             - sj(izz,i,5)*sj(izz,i,4))                     
 1016    continue
cdir$ ivdep
         do 1017 izz=1,nr,jdim   
         t(izz+jdim1,41) = -(sj(izz+jdim1,i,5)*sj(izz+jdim1,i,4)
     .                   - sj(izz+jdim1-1,i,5)*sj(izz+jdim1-1,i,4))
         t(izz,41)       = -(sj(izz+1,i,5)*sj(izz+1,i,4)
     .                   - sj(izz,i,5)*sj(izz,i,4))                     
 1017    continue
cdir$ ivdep
         do 1018 izz=1,nr  
         t(izz,36) = q(izz,1,i,1) 
         t(izz,37) = q(izz,1,i,1)*q(izz,1,i,2) 
         t(izz,38) = q(izz,1,i,1)*q(izz,1,i,3) 
         t(izz,39) = q(izz,1,i,1)*q(izz,1,i,4) 
         t(izz,40) = q(izz,1,i,5)*oogmo
     .             + 0.5*q(izz,1,i,1)*(q(izz,1,i,2)*q(izz,1,i,2)
     .             + q(izz,1,i,3)*q(izz,1,i,3)
     .             + q(izz,1,i,4)*q(izz,1,i,4)) 
 1018    continue 
         do 452 l=1,5
cdir$ ivdep
         do 1019 izz=1,nr
         res(izz,1,i,l) = res(izz,1,i,l) + t(izz,35+l)*t(izz,41) 
 1019    continue
  452    continue
      end if  
c
c      store finer-grid fluxes for enforcing conservation on coarser meshes
c
      if (isf.eq.1) then
      if (jbctyp(1).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),379) nbl
      end if
  379 format(2x,42h installing j fluxes into qj0 for edge j=0,
     .10h for block,i3)
c
c      left boundary
c
      do 333 ipl=1,npl
      ii = i+ipl-1
      do 333 l=1,5
      jk = (ipl-1)*jdim*kdim+1-jdim
      do 333 k=1,kdim1
      jk = jk+jdim
      qj0(k,ii,l,2) = t(jk,30+l)
  333 continue
      end if      
      if (jbctyp(2).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),389) nbl
      end if
  389 format(2x,45h installing j fluxes into qj0 for edge j=jdim,
     .10h for block,i3)
c
c      right boundary
c
      do 444 ipl=1,npl
      ii = i+ipl-1
      do 444 l=1,5
      jk = (ipl-1)*jdim*kdim
      do 444 k=1,kdim1
      jk = jk+jdim
      qj0(k,ii,l,4) = t(jk,30+l)
  444 continue
      end if
      end if
      return
      end
