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
      subroutine pre_embed(nbl,lw,lw2,icount,maxbl,maxseg,lbcemb,
     .                     nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                     nbcidim,jbcinfo,kbcinfo,ibcinfo,iemg,
     .                     igridg,jdimg,kdimg,idimg,isav_emb,
     .                     is_emb,ie_emb,nbcemb,nou,bou,nbuf,ibufdim,
     .                     myid,maxgr,ierrflg)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Set up data arrays for embedded bc message passing
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension lw(65,maxbl),lw2(43,maxbl)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension igridg(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          iemg(maxgr)
      dimension isav_emb(lbcemb,12),is_emb(5),ie_emb(5)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
c
      call lead(nbl,lw,lw2,maxbl)
c
c     check to see if this block is embedded
c
      igrid = igridg(nbl)
      iem   = iemg(igrid)
      if (iem .le. 0) return
c
c     nblc is the coarser block in which the current block is embedded
c     (note: nblc,js,je,...ke are set by the call to lead)
c
      jc = jdimg(nblc)
      kc = kdimg(nblc)
      ic = idimg(nblc) 
c
      nsi = (idim-1)/(ie-is)
c
      if (jdim.ne.(je-js)*2+1 .or. je.gt.jc .or. js.lt.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)' embedding error:  jdim,je,js=',
     .   jdim,je,js
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      else
     .   if (kdim.ne.(ke-ks)*2+1 .or. ke.gt.kc .or. ks.lt.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' embedding error:  kdim,ke,ks=',
     .      kdim,ke,ks
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      else
     .   if (ie.gt.ic .or. is.lt.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' embedding error:  idim,ie,is=',
     .      idim,ie,is
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (nsi.lt.1 .or. nsi.gt.2) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)' embedding error:  idim,ie,is,nsi=',
     .   idim,ie,is,nsi
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if 
c
c     determine all the embeded type bcs for this block
c
      do 802 nseg=1,nbci0(nbl)
      if (ibcinfo(nbl,nseg,1,1).eq.21) then
         nface  = 1
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js 
         isav_emb(icount,6)  = je 
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  802 continue
c
      do 803 nseg=1,nbcidim(nbl)
      if (ibcinfo(nbl,nseg,1,2).eq.21) then
         nface  = 2 
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js
         isav_emb(icount,6)  = je
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  803 continue
c
      do 804 nseg=1,nbcj0(nbl)
      if (jbcinfo(nbl,nseg,1,1).eq.21) then
         nface  = 3
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js
         isav_emb(icount,6)  = je
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  804 continue
c
      do 805 nseg=1,nbcjdim(nbl)
      if (jbcinfo(nbl,nseg,1,2).eq.21) then
         nface  = 4
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js
         isav_emb(icount,6)  = je
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  805 continue
c
      do 806 nseg=1,nbck0(nbl)
      if (kbcinfo(nbl,nseg,1,1).eq.21) then
         nface  = 5
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js
         isav_emb(icount,6)  = je
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  806 continue
c
      do 807 nseg=1,nbckdim(nbl)
      if (kbcinfo(nbl,nseg,1,2).eq.21) then
         nface  = 6
         icount = icount + 1
         isav_emb(icount,1)  = nbl
         isav_emb(icount,2)  = nface
         isav_emb(icount,3)  = is
         isav_emb(icount,4)  = ie
         isav_emb(icount,5)  = js
         isav_emb(icount,6)  = je
         isav_emb(icount,7)  = ks
         isav_emb(icount,8)  = ke
         isav_emb(icount,9)  = nblc
         isav_emb(icount,10) = nsi
      end if
  807 continue
c
      nbcemb = icount
c
      if (nbcemb.gt.lbcemb) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'  Stopping in pre_embed:'
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'  nbcemb too small'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      return
      end
