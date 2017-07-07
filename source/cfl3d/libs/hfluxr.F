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
      subroutine hfluxr(i,npl,xkap,jdim,kdim,idim,res,q,qk0,sk,t,nvtq,
     .                  nv,nfa,wfa,iwfa,kbctyp,isf,nbl,bck,nou,bou,nbuf,
     .                  ibufdim,myid,mblk2nd,maxbl,idef)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute residual contributions for the 
c     right-hand-side in the K-direction from the inviscid terms.
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
      dimension kbctyp(2)
      dimension sk(jdim*kdim,idim-1,5)
      dimension q(jdim,kdim,idim,5),qk0(jdim,idim-1,5,4)
      dimension res(jdim,kdim,idim-1,5),bck(jdim,idim-1,2)
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
      n     = npl*jdim*kdim
      nr    = n-jdim
      nr2   = n-2*jdim
c
      if (real(xkap).lt.-2.e0) then
c
c     first order
c
c     fill temp arrays with interior values
c
      do 100 l=1,5
cdir$ ivdep
      do 50 izz=1,nr
      t(izz+jdim,20+l) = q(izz,1,i,l)
      t(izz,25+l)      = q(izz,1,i,l)
   50 continue
  100 continue
c
c     fill temp arrays with boundary values
c
      do 200 ipl=1,npl
      ii  = i+ipl-1
      jk  = 1+(ipl-1)*jdim*kdim
      jk2 = 1+(ipl-1)*jdim*kdim + kdim1*jdim
      do 150 l=1,5
cdir$ ivdep
      do 120 izz=1,jdim
      t(izz+jk-1,20+l)  = qk0(izz,ii,l,1)
      t(izz+jk2-1,25+l) = qk0(izz,ii,l,3)
  120 continue
  150 continue
  200 continue
c
      else
c
c     higher order
c
      do 1000 l=1,5
c
c     calculate interior gradients
c
cdir$ ivdep
      do 300 izz=1,nr2
      t(izz+jdim,1) = q(izz+jdim,1,i,l)-q(izz,1,i,l)
  300 continue
c
c      edge gradients - left boundary
c
      do 400 ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*kdim*jdim
      do 400 jj=1,jdim1
      jc      = jc + 1
      t(jc,1) = (1.0-bck(jj,ii,1))*(q(jj,1,ii,l)-qk0(jj,ii,l,1))
     .              +bck(jj,ii,1) * qk0(jj,ii,l,2)
  400 continue
c
c      edge gradients - right boundary
c
      do 500 ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*kdim*jdim + kdim1*jdim
      do 500 jj=1,jdim1
      jc      = jc + 1
      t(jc,1) = (1.0-bck(jj,ii,2))*(qk0(jj,ii,l,3)-q(jj,kdim1,ii,l))
     .              +bck(jj,ii,2) * qk0(jj,ii,l,4)
  500 continue
c
c     zero out gradients for dummy interfaces at j=jdim
c
      do 600 ipl=1,npl
      jvl = (ipl-1)*jdim*kdim
      do 600 k=1,kdim
      t(k*jdim+jvl,1) = 0.0
  600 continue
c
c      gradient limiting - cell interface interpolations
c
cdir$ ivdep
      do 700 izz=1,nr
      t(izz,2) = t(izz+jdim,1)
  700 continue
c
      ifl = iflim(3)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = kdim1
      end if
      call xlim(xkap,nr,t(1,1),t(1,2),q(1,1,i,l),ifl,ncells,l)
c
cdir$ ivdep
      do 800 izz=1,nr
      t(izz+jdim,20+l) = q(izz,1,i,l)+t(izz,2)
      t(izz,     25+l) = q(izz,1,i,l)-t(izz,1)
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
      jl = 0
      do 1100 jj=1,jdim1
      jl = jl + 1
      t(jl,1) = qk0(jj,ii,l,1) - qk0(jj,ii,l,2)
      t(jl,2) = q(jj,1,ii,l)   - qk0(jj,ii,l,1)
      t(jl,3) = qk0(jj,ii,l,1)
 1100 continue
      ifl = iflim(3)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = kdim1
      end if
      call xlim(xkap,jdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = (ipl-1)*(kdim*jdim)
      jl = 0
      do 1200 jj=1,jdim1
      jl = jl + 1
      jc = jc + 1
      t(jc,20+l) = (1.0-bck(jj,ii,1))*(qk0(jj,ii,l,1)+t(jl,2))
     .                 +bck(jj,ii,1) * qk0(jj,ii,l,1)
      t(jc,25+l) = (1.0-bck(jj,ii,1))* t(jc,25+l)
     .                 +bck(jj,ii,1) * qk0(jj,ii,l,1)
 1200 continue
 1300 continue
c
c      edge values - right boundary
c
      do 1800 ipl=1,npl
      ii = i+ipl-1
      jl = 0
      do 1600 jj=1,jdim1
      jl = jl + 1
      t(jl,2) = qk0(jj,ii,l,4) - qk0(jj,ii,l,3)
      t(jl,1) = qk0(jj,ii,l,3) - q(jj,kdim1,ii,l)
      t(jl,3) = qk0(jj,ii,l,3)
 1600 continue
      ifl = iflim(3)
      if (ifl.eq.4) then
        if (i2d.eq.1) then
          ncells = sqrt(float(ncell(level)))
        else
          ncells = float(ncell(level))**(1./3.)
        end if
      else
        ncells = kdim1
      end if
      call xlim(xkap,jdim1,t(1,1),t(1,2),t(1,3),ifl,ncells,l)
      jc = (ipl-1)*(kdim*jdim) + kdim1*jdim
      jl = 0
      do 1700 jj=1,jdim1
      jl = jl + 1
      jc = jc + 1
      t(jc,20+l) = (1.0-bck(jj,ii,2))* t(jc,20+l)
     .                 +bck(jj,ii,2) * qk0(jj,ii,l,3)
      t(jc,25+l) = (1.0-bck(jj,ii,2))*(qk0(jj,ii,l,3)-t(jl,1))
     .                 +bck(jj,ii,2) * qk0(jj,ii,l,3)
 1700 continue
 1800 continue
c
c     fill end points for safety
      t(jdim,20+l)  = t(jdim-1,20+l)
      t(jdim*kdim*npl,25+l)  = t(jdim*kdim*npl-1,25+l)
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
         ic = (ipl-1)*jdim*kdim + (k-1)*jdim + j
         if (real(t(ic,21)).lt.real(epsz) .or.
     .       real(t(ic,25)).lt.real(epsz) .or.
     .       real(t(ic,21)).gt.real(epss) .or.
     .       real(t(ic,25)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in hflux left - small ',
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
            write(bou(nou(1),1),*)' stopping in hflux right - small ',
     .                 ' (or large) density and/or pressure at '
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' j,k,i,t(26),t(30) = ',j,k,ii,
     .      real(t(ic,26)),real(t(ic,30))
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
 5432    continue
      end if
c
      if (ifds(3).eq.0) then
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),185) nbl
      end if
  185 format(49h   computing inviscid fluxes, K-direction - flux-,
     .24hvector splitting - block,i4)
c
      call fluxp(sk(1,i,1),sk(1,i,2),sk(1,i,3),sk(1,i,4),sk(1,i,5),
     .           t(1,21),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
c
      do 1150 l=1,5
cdir$ ivdep
      do 1011 izz=1,n
      t(izz,20+l) = t(izz,30+l)
 1011 continue
 1150 continue
c
      call fluxm(sk(1,i,1),sk(1,i,2),sk(1,i,3),sk(1,i,4),sk(1,i,5),
     .           t(1,26),t(1,31),n,t,n,nvtq,nou,bou,nbuf,ibufdim)
c
      do 1400 l=1,5
cdir$ ivdep
      do 1013 izz=1,n
      t(izz,30+l) = t(izz,30+l)+t(izz,20+l)
 1013 continue
 1400 continue
c
      else if (ifds(3).eq.1) then
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),184)
      end if
  184 format(49h   computing inviscid fluxes, K-direction - flux-,
     .20hdifference splitting)
c
      call fhat(sk(1,i,1),sk(1,i,2),sk(1,i,3),sk(1,i,4),sk(1,i,5),
     .          t(1,31),t(1,26),t(1,21),n,nvtq) 
c
      else
c
      if (isklton.gt.0 .and. i.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),183)
      end if
  183 format(49h   computing inviscid fluxes, K-direction - MAPS+,
     .15h flux splitting)
c
      call fmaps(sk(1,i,1),sk(1,i,2),sk(1,i,3),sk(1,i,4),sk(1,i,5),
     .          t(1,31),t(1,26),t(1,21),n,nvtq)
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
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),269) ksb,jsb,jeb,isb,ieb
      end if
  269 format(2x,33h installing accumulated fluxes on,
     .3h k=,i3,16h at js,je,is,ie=,4i4)
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
         lfc  = ifts+(l-1)*(jeb-jsb)*(ieb-isb)+(ii-isb)*(jeb-jsb)
         loc  = (ipl-1)*jdim*kdim+(ksb-1)*jdim+jsb
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
      do 450 l=1,5
cdir$ ivdep
      do 1014 izz=1,nr
      res(izz,1,i,l) = res(izz,1,i,l)+t(izz+jdim,30+l)-t(izz,30+l)
 1014 continue
  450 continue
c
c     geometric conservation law terms for deforming grids
c
      if (idef.gt.0) then
         oogmo = 1./(gamma-1.)
cdir$ ivdep
         do 1016 izz=1,nr
         t(izz,41) = -(sk(izz+jdim,i,5)*sk(izz+jdim,i,4)
     .             - sk(izz,i,5)*sk(izz,i,4))                     
 1016    continue
         do 1018 ii = i,i+npl-1
         ijk2  = (ii-i+1)*jdim*kdim-jdim  
cdir$ ivdep
         do 1017 izz=1,jdim  
           t(izz+ijk2,41) = -(sk(izz+ijk2,i,5)*sk(izz+ijk2,i,4)
     .            - sk(izz+ijk2-jdim,i,5)*sk(izz+ijk2-jdim,i,4))                     
 1017    continue
 1018    continue 
         do 1019 izz=1,nr    
         t(izz,36) = q(izz,1,i,1) 
         t(izz,37) = q(izz,1,i,1)*q(izz,1,i,2) 
         t(izz,38) = q(izz,1,i,1)*q(izz,1,i,3) 
         t(izz,39) = q(izz,1,i,1)*q(izz,1,i,4) 
         t(izz,40) = q(izz,1,i,5)*oogmo
     .             + 0.5*q(izz,1,i,1)*(q(izz,1,i,2)*q(izz,1,i,2)
     .             + q(izz,1,i,3)*q(izz,1,i,3)
     .             + q(izz,1,i,4)*q(izz,1,i,4)) 
 1019    continue 
         do 452 l=1,5
cdir$ ivdep
         do 1020 izz=1,nr
         res(izz,1,i,l) = res(izz,1,i,l) + t(izz,35+l)*t(izz,41) 
 1020    continue
  452    continue
      end if 
c
c      store finer-grid fluxes for enforcing conservation on coarser meshes
c
      if (isf.eq.1) then
      if (kbctyp(1).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),379) nbl
      end if
  379 format(2x,42h installing k fluxes into qk0 for edge k=0,
     .10h for block,i3)
c
c      left boundary
c
      do 333 ipl=1,npl
      ii = i+ipl-1
      do 333 l=1,5
      jk = (ipl-1)*jdim*kdim
      do 333 j=1,jdim1
      jk = jk+1
      qk0(j,ii,l,2) = t(jk,30+l)
  333 continue
      end if      
      if (kbctyp(2).eq.21) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),389) nbl
      end if
  389 format(2x,45h installing k fluxes into qk0 for edge k=kdim,
     .10h for block,i3)
c
c      right boundary
c
      do 444 ipl=1,npl
      ii = i+ipl-1
      do 444 l=1,5
      jk = (ipl-1)*jdim*kdim+jdim*kdim1
      do 444 j=1,jdim1
      jk = jk+1
      qk0(j,ii,l,4) = t(jk,30+l)
  444 continue
      end if
      end if
      return
      end
