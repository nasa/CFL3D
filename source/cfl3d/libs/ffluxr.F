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
      subroutine ffluxr(k,npl,xkap,jdim,kdim,idim,res,q,qi0,si,t,nvtq,
     .                  nv,nfa,wfa,iwfa,ibctyp,isf,nbl,bci,nou,bou,nbuf,
     .                  ibufdim,myid,mblk2nd,maxbl,idef)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute residual contributions for the 
c     right-hand-side in the I-direction from the inviscid terms.
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
      dimension ibctyp(2)
      dimension si(jdim*kdim,idim,5)
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),bci(jdim,kdim,2)
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
      jdim1 = jdim-1
      idim1 = idim-1
      kdim1 = kdim-1
c
      if (npl.eq.kdim1 .and. nvtq.ge.jdim*kdim*idim) npl = kdim
c
      jv = npl*jdim
      n  = jv*idim
      nr = n-jv
c
      if (npl.ne.kdim) then
c
c        fill temp arrays with metrics if all k planes are not done at once
c
         do 40 i=1,idim
         jk  = (i-1)*jv+1
         jk1 = (k-1)*jdim+1
         do 40 l=1,5
cdir$ ivdep
         do 10 izz=1,jv
         t(izz+jk-1,14+l) = si(izz+jk1-1,i,l)
   10    continue
   40    continue
      end if
c
c     all orders
c
c     fill temp arrays with interior values
c
      do 100 l=1,5
      if (npl.eq.kdim) then
cdir$ ivdep
         do 50 izz=1,nr
         t(izz+jv,20+l) = q(izz,1,1,l)
         t(izz,25+l)    = q(izz,1,1,l)
   50    continue
      else
         do 70 i=1,idim1
         jk = (i-1)*jv+1
cdir$ ivdep
         do 60 izz=1,jv
         t(izz+jv+jk-1,20+l) = q(izz,k,i,l)
         t(izz+jk-1,25+l)    = q(izz,k,i,l)
   60    continue
   70    continue
      end if
c
c      fill temp array with boundary values
c
cdir$ ivdep
      do 80 izz=1,jv
      t(izz,20+l)    = qi0(izz,k,l,1)
      t(izz+nr,25+l) = qi0(izz,k,l,3)
   80 continue
  100 continue
c
      if (real(xkap).gt.-2.e0) then
c
c     higher order
c
      do 1000 l=1,5
c
c     gradients on interior planes
c
cdir$ ivdep
      do 300 izz=1,nr-jv
      t(izz+jv,1) = t(izz+jv,25+l)-t(izz,25+l)
  300 continue
c
c     initialize gradients to zero on planes i=1 and i=idim
c
      kv = jv*idim1
      do 350 izz=1,jv
      t(izz,1)    = 0.0
      t(kv+izz,1) = 0.0
  350 continue
c
c     edge gradients - left boundary
c
      do 400 kpl=1,npl
      kk = k+kpl-1
      jc = (kpl-1)*jdim
      do 400 jj=1,jdim1
      jc = jc + 1
      t(jc,1) = (1.0-bci(jj,kk,1))*(q(jj,kk,1,l)-qi0(jj,kk,l,1))
     .              +bci(jj,kk,1) * qi0(jj,kk,l,2)
  400 continue
c
c     edge gradients - right boundary
c
      do 500 kpl=1,npl
      kk = k+kpl-1
      jc = idim1*jv + (kpl-1)*jdim
      do 500 jj=1,jdim1
      jc = jc+1
      t(jc,1) = (1.0-bci(jj,kk,2))*(qi0(jj,kk,l,3)-q(jj,kk,idim1,l))
     .              +bci(jj,kk,2) * qi0(jj,kk,l,4)
  500 continue
c
c      gradient limiting - cell interface interpolations
c
cdir$ ivdep
      do 700 izz=1,nr
      t(izz,2) = t(izz+jv,1)
  700 continue
c
      ifl = iflim(1)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = idim1
      end if
      call xlim(xkap,nr,t(1,1),t(1,2),t(1,25+l),ifl,ncells,l)
c
cdir$ ivdep
      do 800 izz=1,nr
      t(izz+jv,20+l) = t(izz+jv,20+l)+t(izz,2)
      t(izz,   25+l) = t(izz,   25+l)-t(izz,1)
  800 continue
 1000 continue
      end if
c
      do 2000 l=1,5
c
c     edge values - left boundary
c
      do 1300 kpl=1,npl
      kk = k+kpl-1
      jl = 0
      do 1100 jj=1,jdim1
      jl = jl + 1
      t(jl,1) = qi0(jj,kk,l,1) - qi0(jj,kk,l,2)
      t(jl,2) = q(jj,kk,1,l)   - qi0(jj,kk,l,1)
      t(jl,3) = qi0(jj,kk,l,1)
 1100 continue
      ifl = iflim(1)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = idim1
      end if
      call xlim(xkap,jdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = (kpl-1)*jdim
      jl = 0
      do 1200 jj=1,jdim1
      jl = jl + 1
      jc = jc + 1
      t(jc,20+l) = (1.0-bci(jj,kk,1))*(qi0(jj,kk,l,1)+t(jl,2))
     .                 +bci(jj,kk,1) * qi0(jj,kk,l,1)
      t(jc,25+l) = (1.0-bci(jj,kk,1))* t(jc,25+l)
     .                 +bci(jj,kk,1) * qi0(jj,kk,l,1)
 1200 continue
 1300 continue
c
c     edge values - right boundary
c
      do 1800 kpl=1,npl
      kk = k+kpl-1
      jl = 0
      do 1600 jj=1,jdim1
      jl = jl + 1
      t(jl,2) = qi0(jj,kk,l,4) - qi0(jj,kk,l,3)
      t(jl,1) = qi0(jj,kk,l,3) - q(jj,kk,idim1,l)
      t(jl,3) = qi0(jj,kk,l,3)
 1600 continue
      ifl = iflim(1)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = idim1
      end if
      call xlim(xkap,jdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = idim1*jv + (kpl-1)*jdim
      jl = 0
      do 1700 jj=1,jdim1
      jl = jl + 1
      jc = jc + 1
      t(jc,20+l) = (1.0-bci(jj,kk,2))* t(jc,20+l)
     .                 +bci(jj,kk,2) * qi0(jj,kk,l,3)
      t(jc,25+l) = (1.0-bci(jj,kk,2))*(qi0(jj,kk,l,3)-t(jl,1))
     .                 +bci(jj,kk,2) * qi0(jj,kk,l,3)
 1700 continue
 1800 continue
c
c     fill end points for safety
      do 802 kpl=1,npl
      jc = (kpl-1)*jdim+jdim
      jl = idim1*jv + (kpl-1)*jdim + jdim
      t(jc,20+l)  = t(jc-1,20+l)
      t(jl,25+l)  = t(jl-1,25+l)
  802 continue
c
 2000 continue
c
      if (ichk.eq.1) then
         epsz = 1.0e-03
         epss = 1.0e+03
         do 5432 ipl=1,npl
         kk = k+ipl-1
         do 5432 i=1,idim
         do 5432 j=1,jdim
         kc = jv*(i-1) + (ipl-1)*jdim + j
         if (real(t(kc,21)).lt.real(epsz) .or.
     .       real(t(kc,25)).lt.real(epsz) .or.
     .       real(t(kc,21)).gt.real(epss) .or.
     .       real(t(kc,25)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in fflux left - small ',
     .                 '(or large) density and/or pressure at '
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' j,k,i,t(21),t(25) = ',j,kk,i,
     .      real(t(kc,21)),real(t(kc,25))
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (real(t(kc,26)).lt.real(epsz) .or.
     .       real(t(kc,30)).lt.real(epsz) .or.
     .       real(t(kc,26)).gt.real(epss) .or.
     .       real(t(kc,30)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in fflux right - small',
     .                 ' (or large) density and/or pressure at '
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' j,k,i,t(26),t(30) = ',j,kk,i,
     .      real(t(kc,26)),real(t(kc,30))
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
 5432    continue
      end if
c
      if (ifds(1).eq.0) then
c
      if (isklton.gt.0 .and. k.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),185) nbl
      end if
  185 format(49h   computing inviscid fluxes, I-direction - flux-,
     .24hvector splitting - block,i4)
c
      if (npl.eq.kdim) then
         call fluxp(si(1,1,1),si(1,1,2),si(1,1,3),si(1,1,4),si(1,1,5),
     .              t(1,21),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
      else
         call fluxp(t(1,15),t(1,16),t(1,17),t(1,18),t(1,19),
     .              t(1,21),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
      end if
c
      do 150 l=1,5
cdir$ ivdep
      do 1015 izz=1,n
      t(izz,20+l) = t(izz,30+l)
 1015 continue
  150 continue
c
      if (npl.eq.kdim) then
         call fluxm(si(1,1,1),si(1,1,2),si(1,1,3),si(1,1,4),si(1,1,5),
     .              t(1,26),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
      else
         call fluxm(t(1,15),t(1,16),t(1,17),t(1,18),t(1,19),
     .              t(1,26),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
      end if
c
      do 1400 l=1,5
cdir$ ivdep
      do 1018 izz=1,n
      t(izz,30+l) = t(izz,30+l) + t(izz,20+l)
 1018 continue
 1400 continue
c
      else if (ifds(1).eq.1) then
c
      if (isklton.gt.0 .and. k.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),184)
      end if
  184 format(49h   computing inviscid fluxes, I-direction - flux-,
     .20hdifference splitting)
      if (npl.eq.kdim) then
         call fhat(si(1,1,1),si(1,1,2),si(1,1,3),si(1,1,4),si(1,1,5),
     .             t(1,31),t(1,26),t(1,21),n,nvtq) 
      else
         call fhat(t(1,15),t(1,16),t(1,17),t(1,18),t(1,19),
     .             t(1,31),t(1,26),t(1,21),n,nvtq) 
      end if
c
      else
c
      if (isklton.gt.0 .and. k.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),183)
      end if
  183 format(49h   computing inviscid fluxes, I-direction - MAPS+,
     .15h flux splitting)
      if (npl.eq.kdim) then
         call fmaps(si(1,1,1),si(1,1,2),si(1,1,3),si(1,1,4),si(1,1,5),
     .             t(1,31),t(1,26),t(1,21),n,nvtq)
      else
         call fmaps(t(1,15),t(1,16),t(1,17),t(1,18),t(1,19),
     .             t(1,31),t(1,26),t(1,21),n,nvtq)
      end if
c
      end if
c
c      conservative at embedded boundaries - wfa array contains fluxes
c      from a finer grid
c
      if (nfa.gt.0) then
      kks  = k
      kke  = k+npl-1
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
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),269) isb,jsb,jeb,ksb,keb
      end if
  269 format(2x,33h installing accumulated fluxes on,
     .3h i=,i3,16h at js,je,ks,ke=,4i4)
c
c     loop over all planes in K-direction
c
      lfcc = ifts
      do 754 kk=1,kdim1
c
c     skip planes not in embedded region
c
      if (kk.ge.ksb .and. kk.lt.keb) then
c
c     check for kk in region to be updated on this pass (npl planes of data)
c
      if (kks.le.kk .and. kk.le.kke) then
         kpl = kk-k+1
         do 705 l=1,5
         lfc = ifts+(l-1)*(keb-ksb)*(jeb-jsb)+(kk-ksb)*(jeb-jsb)
         loc = (kpl-1)*jdim+(isb-1)*jdim*npl+jsb
         do 705 j=jsb,jeb-1
         t(loc,30+l) = wfa(lfc)
         loc  = loc+1
         lfc  = lfc+1
         lfcc = lfcc+1
  705    continue
      else
c
c     increment counter for planes not updated on this pass
c
      lfcc = lfcc+5*(jeb-jsb)
      end if
      end if
  754 continue
 7055 continue
      end if
c
      if (npl.eq.kdim) then
         do 450 l=1,5
cdir$ ivdep
         do 1020 izz=1,nr
         res(izz,1,1,l) = res(izz,1,1,l)+t(izz+jv,30+l)-t(izz,30+l)
 1020    continue
  450    continue
c
c        geometric conservation law terms for deforming grids
c
         if (idef.gt.0) then 
            oogmo = 1./(gamma-1.)
cdir$ ivdep
            do 10160 izz=1,nr
            t(izz,36) = q(izz,1,1,1)
            t(izz,37) = q(izz,1,1,1)*q(izz,1,1,2)
            t(izz,38) = q(izz,1,1,1)*q(izz,1,1,3)
            t(izz,39) = q(izz,1,1,1)*q(izz,1,1,4)
            t(izz,40) = q(izz,1,1,5)*oogmo
     .                + 0.5*q(izz,1,1,1)*(q(izz,1,1,2)*q(izz,1,1,2)
     .                + q(izz,1,1,3)*q(izz,1,1,3)
     .                + q(izz,1,1,4)*q(izz,1,1,4))
            t(izz,41) = -(si(izz+jv,1,5)*si(izz+jv,1,4)
     .                - si(izz,1,5)*si(izz,1,4))
10160       continue
            do 4510 l=1,5
cdir$ ivdep
            do 10220 izz=1,nr
            res(izz,1,1,l) = res(izz,1,1,l) + t(izz,35+l)*t(izz,41)
10220       continue
 4510       continue
         end if  
       else
         do 451 i=1,idim1
         jk = (i-1)*jv
         do 451 l=1,5
cdir$ ivdep
         do 1021 izz=1,jv
         res(izz,k,i,l) = res(izz,k,i,l)+t(izz+jk+jv,30+l)
     .                                  -t(izz+jk,30+l)
 1021    continue
  451    continue
c
c        geometric conservation law terms for deforming grids
c
         if (idef.gt.0) then
            oogmo = 1./(gamma-1.)
            jk1   = (k-1)*jdim
            do 4520 i=1,idim1
cdir$ ivdep
            do 10230 izz=1,jv
            izz1 = izz + jk1
            t(izz,36) = q(izz,k,i,1) 
            t(izz,37) = q(izz,k,i,1)*q(izz,k,i,2) 
            t(izz,38) = q(izz,k,i,1)*q(izz,k,i,3) 
            t(izz,39) = q(izz,k,i,1)*q(izz,k,i,4) 
            t(izz,40) = q(izz,k,i,5)*oogmo
     .                + 0.5*q(izz,k,i,1)*(q(izz,k,i,2)*q(izz,k,i,2)
     .                + q(izz,k,i,3)*q(izz,k,i,3) 
     .                + q(izz,k,i,4)*q(izz,k,i,4)) 
            t(izz,41) = -(si(izz1,i+1,5)*si(izz1,i+1,4)
     .                - si(izz1,i,5)*si(izz1,i,4))
10230       continue
c
            do 4550 l=1,5
cdir$ ivdep
            do 10240 izz=1,jv
            res(izz,k,i,l) = res(izz,k,i,l)
     .                     + t(izz,35+l)*t(izz,41)
10240       continue
 4550       continue
c
 4520       continue  
         end if 
      end if
c
c      store finer-grid fluxes for enforcing conservation on coarser meshes
c
      if (isf.eq.1) then
      if (ibctyp(1).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),379) nbl
      end if
  379 format(2x,42h installing i fluxes into qi0 for edge i=0,
     .10h for block,i3)
c
c      left boundary
c
      do 333 kpl=1,npl
      kk = k+kpl-1
      do 333 l=1,5
      jk = (kpl-1)*jdim
      do 333 j=1,jdim1
      jk = jk+1
      qi0(j,kk,l,2) = t(jk,30+l)
  333 continue
      end if      
      if (ibctyp(2).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),389) nbl
      end if
  389 format(2x,45h installing i fluxes into qi0 for edge i=idim,
     .10h for block,i3)
c
c      right boundary
c
      do 444 kpl=1,npl
      kk = k+kpl-1
      do 444 l=1,5
      jk = idim1*jv+(kpl-1)*jdim
      do 444 j=1,jdim1
      jk = jk+1
      qi0(j,kk,l,4) = t(jk,30+l)
  444 continue
      end if
      end if
      return
      end
