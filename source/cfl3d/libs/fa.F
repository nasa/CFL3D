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
      subroutine fa(nbl,lw,w,mgwk,wk,nwork,iw,ifamax,nwfa,nifa,nfajki,
     .              maxbl,maxseg,jdimg,kdimg,idimg,jsg,ksg,isg,jeg,
     .              keg,ieg,jbcinfo,kbcinfo,ibcinfo,nblock,nblcg,
     .              nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Accumulate fluxes to insure conservation.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
c     maxbl - maximum number of blocks
c
      dimension nou(nbuf)
      dimension w(mgwk)
      dimension lw(65,maxbl)
      dimension wk(nwork)
      dimension iw(ifamax),nfajki(3)
      dimension jbctyp(2),kbctyp(2),ibctyp(2)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl),nblcg(maxbl),
     .          jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
c
      common /sklton/ isklton
c
      nfajki(1) = 0
      nfajki(2) = 0
      nfajki(3) = 0
c
c      nwfa  = total number of coarser grid fluxes accumulated from
c              finer grids
c      nifa  = number of integers stored in iw array, equal
c              to sum of neta in 3 directions times seven
c      neta  = number of finer grid boundary accumulations in either the
c              j-/k-/i-directions, viz.
c              nfajki(1) = neta in j-direction
c              nfajki(2) = neta in k-direction
c              nfajki(3) = neta in i-direction
c
      nwfa = 0
      nifa = 0
      neta = 0
c
c      accumulate fluxes to ensure conservation - J direction
c
      do 6500 nblc=1,nblock
      if (nbl.eq.nblc) go to 6500
c
      nblcc    = nblcg(nblc)
      if (nblcc.eq.nbl) then
         jf       = jdimg(nblc)
         kf       = kdimg(nblc)
         if       = idimg(nblc)
         js       = jsg(nblc)
         ks       = ksg(nblc)
         is       = isg(nblc)
         je       = jeg(nblc)
         ke       = keg(nblc)
         ie       = ieg(nblc)
         jbctyp(1) = jbcinfo(nblc,1,1,1)
         jbctyp(2) = jbcinfo(nblc,1,1,2)
c
         if ((isklton.gt.0 .and. jbctyp(1).eq.21) .or.
     .       (isklton.gt.0 .and. jbctyp(2).eq.21)) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),*)
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),7)nblc,nbl
         end if
    7    format(1x,33hflux accumulation in j from block,i3,
     .          1x,16hfor use on block,i3)
c
         call fa2xj(jf,kf,if,w(lw(2,nblc)),js,ks,is,je,ke,ie,
     .              jbctyp,wk,iw,nwfa,nifa,neta,nou,bou,nbuf,ibufdim)
      end if
 6500 continue
c
      nfajki(1) = neta
      neta = 0
c
c      K-direction
c
      do 6502 nblc=1,nblock
      if (nbl.eq.nblc) go to 6502
c
      nblcc    = nblcg(nblc)
      if (nblcc.eq.nbl) then
         jf       = jdimg(nblc)
         kf       = kdimg(nblc)
         if       = idimg(nblc)
         js       = jsg(nblc)
         ks       = ksg(nblc)
         is       = isg(nblc)
         je       = jeg(nblc)
         ke       = keg(nblc)
         ie       = ieg(nblc)
         kbctyp(1) = kbcinfo(nblc,1,1,1)
         kbctyp(2) = kbcinfo(nblc,1,1,2)
c
         if ((isklton.gt.0 .and. kbctyp(1).eq.21) .or.
     .       (isklton.gt.0 .and. kbctyp(2).eq.21))  then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),*)
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),17)nblc,nbl
         end if
   17    format(1x,33hflux accumulation in k from block,i3,
     .          1x,16hfor use on block,i3)
c
         call fa2xk(jf,kf,if,w(lw(3,nblc)),js,ks,is,je,ke,ie,
     .              kbctyp,wk,iw,nwfa,nifa,neta,nou,bou,nbuf,ibufdim)
      end if
 6502 continue
c
      nfajki(2) = neta
      neta      = 0
c
c      I-direction
c
      do 6504 nblc=1,nblock
      if (nbl.eq.nblc) go to 6504
c
      nblcc    = nblcg(nblc)
      if (nblcc.eq.nbl) then
         jf       = jdimg(nblc)
         kf       = kdimg(nblc)
         if       = idimg(nblc)
         js       = jsg(nblc)
         ks       = ksg(nblc)
         is       = isg(nblc)
         je       = jeg(nblc)
         ke       = keg(nblc)
         ie       = ieg(nblc)
         ibctyp(1) = ibcinfo(nblc,1,1,1)
         ibctyp(2) = ibcinfo(nblc,1,1,2)
c
         if ((isklton.gt.0 .and. ibctyp(1).eq.21) .or.
     .       (isklton.gt.0 .and. ibctyp(2).eq.21))  then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),*)
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),71)nblc,nbl
         end if
   71    format(1x,33hflux accumulation in i from block,i3,
     .          1x,16hfor use on block,i3)
c
         call fa2xi(jf,kf,if,w(lw(4,nblc)),js,ks,is,je,ke,ie,
     .              ibctyp,wk,iw,nwfa,nifa,neta,nou,bou,nbuf,ibufdim)
      end if
 6504 continue
      nfajki(3) = neta
      return
      end
