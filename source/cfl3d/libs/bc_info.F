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
      subroutine bc_info(ntime,nbl,w,mgwk,ibcinfo,jbcinfo,
     .                   kbcinfo,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                   nbckdim,bcfilei,bcfilej,bcfilek,igridg,itrans,
     .                   irotat,idefrm,idimg,jdimg,kdimg,nblg,iibg,
     .                   kkbg,jjbg,ibcg,ibpntsg,iipntsg,lig,lbg,
     .                   isav_blk,isav_prd,iemg,nbli,
     .                   geom_miss,period_miss,ninter,iindex,nou,bou,
     .                   nbuf,myid,mblk2nd,ibufdim,maxbl,maxseg,iitot,
     .                   intmax,nsub1,mxbli,lbcprd,epsc0,epsrot,lwdat,
     .                   maxgr,iovrlp,bcfiles,mxbcfil)
c
c     $Id$
c
c***********************************************************************
c      Purpose: print out boundary data/conditions to main output file.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 filname
      character*80 bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
c
      integer bcfilei,bcfilej,bcfilek
c
      dimension nou(nbuf)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension w(mgwk),mblk2nd(maxbl),lwdat(maxbl,maxseg,6)
      dimension iibg(iitot),kkbg(iitot),jjbg(iitot),ibcg(iitot),
     .          ibpntsg(maxbl,4),iipntsg(maxbl),lig(maxbl),lbg(maxbl),
     .          iovrlp(maxbl)
      dimension isav_blk(2*mxbli,17),isav_prd(lbcprd,12)
      dimension geom_miss(2*mxbli),period_miss(lbcprd)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2),
     .          idimg(maxbl),jdimg(maxbl),kdimg(maxbl)
      dimension igridg(maxbl),nblg(maxgr),iemg(maxgr),
     .          itrans(maxbl),irotat(maxbl),idefrm(nbl)
      dimension iindex(intmax,6*nsub1+9)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /sklton/ isklton
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /is_perbc/ is_prd(5),ie_prd(5),nbcprd
c
      if (isklton.ne.1) return
c
      iuns = max(irotat(nbl),itrans(nbl),idefrm(nbl))
c
      if (iuns.eq.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl,igridg(nbl)
      end if
      if (iuns.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),8) nbl,igridg(nbl)
      end if
c
    7 format(30h boundary conditions for block,i6,1x,6h(grid ,i6,1h))
    8 format(30h boundary conditions for block,i6,1x,6h(grid ,i6,1h),
     .       1x,14h- dynamic mesh)
c
c***********************************************************************
c      physical boundary conditions.
c***********************************************************************
c
c********************
c     i=0 boundary
c********************
c
      do 802 nseg=1,nbci0(nbl)
      ista = 1
      iend = 1
      jsta = ibcinfo(nbl,nseg,2,1)
      jend = ibcinfo(nbl,nseg,3,1)
      ksta = ibcinfo(nbl,nseg,4,1)
      kend = ibcinfo(nbl,nseg,5,1)
      ldata= lwdat(nbl,nseg,1)
      mdim = jend-jsta
      ndim = kend-ksta
      filname = bcfiles(bcfilei(nbl,nseg,1))
c
      if (ibcinfo(nbl,nseg,1,1).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1012)
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2005) then
c       ngnew is the block number of the block periodic with block nbl
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 1
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk 
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (ibcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (ibcinfo(nbl,nseg,1,1).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,1,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  802 continue
c
c********************
c     i=idim boundary
c********************
c
      do 803 nseg=1,nbcidim(nbl)
      ista = idim
      iend = idim
      jsta = ibcinfo(nbl,nseg,2,2)
      jend = ibcinfo(nbl,nseg,3,2)
      ksta = ibcinfo(nbl,nseg,4,2)
      kend = ibcinfo(nbl,nseg,5,2)
      ldata= lwdat(nbl,nseg,2)
      mdim = jend-jsta
      ndim = kend-ksta
      filname = bcfiles(bcfilei(nbl,nseg,2))
c
      if (ibcinfo(nbl,nseg,1,2).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1012) 
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .               nou,bou,nbuf,ibufdim)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2005) then
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 2
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (ibcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (ibcinfo(nbl,nseg,1,2).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,2,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  803 continue
c
c********************
c     j=0 boundary
c********************
c
      do 804 nseg=1,nbcj0(nbl)
      ista = jbcinfo(nbl,nseg,2,1)
      iend = jbcinfo(nbl,nseg,3,1)
      jsta = 1
      jend = 1
      ksta = jbcinfo(nbl,nseg,4,1)
      kend = jbcinfo(nbl,nseg,5,1)
      ldata= lwdat(nbl,nseg,3)
      mdim = kend-ksta
      ndim = iend-ista
      filname = bcfiles(bcfilej(nbl,nseg,1))
c
      if (jbcinfo(nbl,nseg,1,1).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1012)
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2005) then
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 3
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (jbcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (jbcinfo(nbl,nseg,1,1).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,3,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  804 continue
c
c********************
c     j=jdim boundary
c********************
c
      do 805 nseg=1,nbcjdim(nbl)
      ista = jbcinfo(nbl,nseg,2,2)
      iend = jbcinfo(nbl,nseg,3,2)
      jsta = jdim
      jend = jdim 
      ksta = jbcinfo(nbl,nseg,4,2)
      kend = jbcinfo(nbl,nseg,5,2)
      ldata= lwdat(nbl,nseg,4)
      mdim = kend-ksta
      ndim = iend-ista
      filname = bcfiles(bcfilej(nbl,nseg,2))
c
      if (jbcinfo(nbl,nseg,1,2).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1012)
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .               nou,bou,nbuf,ibufdim)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2005) then
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 4
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (jbcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (jbcinfo(nbl,nseg,1,2).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,4,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  805 continue
c
c********************
c     k=0 boundary
c********************
c
      do 806 nseg=1,nbck0(nbl)
      ista = kbcinfo(nbl,nseg,2,1)
      iend = kbcinfo(nbl,nseg,3,1)
      jsta = kbcinfo(nbl,nseg,4,1)
      jend = kbcinfo(nbl,nseg,5,1)
      ksta = 1
      kend = 1
      ldata= lwdat(nbl,nseg,5)
      mdim = jend-jsta
      ndim = iend-ista
      filname = bcfiles(bcfilek(nbl,nseg,1))
c
      if (kbcinfo(nbl,nseg,1,1).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1012)
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2005) then
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 5
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (kbcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (kbcinfo(nbl,nseg,1,1).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,5,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  806 continue
c
c********************
c     k=kdim boundary
c********************
c
      do 807 nseg=1,nbckdim(nbl)
      ista = kbcinfo(nbl,nseg,2,2)
      iend = kbcinfo(nbl,nseg,3,2)
      jsta = kbcinfo(nbl,nseg,4,2)
      jend = kbcinfo(nbl,nseg,5,2)
      ksta = kdim
      kend = kdim
      ldata= lwdat(nbl,nseg,6)
      mdim = jend-jsta
      ndim = iend-ista
      filname = bcfiles(bcfilek(nbl,nseg,2))
c
      if (kbcinfo(nbl,nseg,1,2).eq.9999)
     .  call out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1000)
     .  call out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1001)
     .  call out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1002)
     .  call out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1003)
     .  call out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim,myid)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1005)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim,0)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1006)
     .  call out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim,1)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1008)
     .  call out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1011)
     .  call out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1012)
     .  call out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1013)
     .  call out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .               nou,bou,nbuf,ibufdim)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2002)
     .  call out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2003)
     .  call out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2009)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2010)
     .  call out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2019)
     .  call out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2004)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2014)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2024)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2034)
     .  call out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2005) then
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        jdimp = jdimg(ngnew)
        kdimp = kdimg(ngnew)
        idimp = idimg(ngnew)
        nblp  = ngnew
        call out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),nbl,jdimp,kdimp,idimp,
     .  nblp,filname,myid,mblk2nd,maxbl)
        nface = 6
        do lcntp=is_prd(level),ie_prd(level)
           nbl_chk   = isav_prd(lcntp,1)
           nface_chk = isav_prd(lcntp,2)
           nseg_chk  = isav_prd(lcntp,11)
           if (nbl.eq.nbl_chk .and. nface.eq.nface_chk
     .        .and. nseg.eq.nseg_chk) then
              if (real(period_miss(lcntp)).gt.real(epsrot)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           blocking '',
     .           ''check....periodic mismatch = '',e14.7)')
     .           real(period_miss(lcntp))
              end if
           end if
        end do
      end if
c
      if (kbcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ngnew = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
        end if
        call out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,ngnew,myid,
     .  mblk2nd,maxbl)
      end if
c
      if (kbcinfo(nbl,nseg,1,2).eq.2007)
     .  call out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2008)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2016)
     .  call out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,5)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2026)
     .  call out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2018)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2028)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,2)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2038)
     .  call out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,3)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2102)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,0)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2103)
     .  call out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,6,
     .  nou,bou,nbuf,ibufdim,mdim,ndim,w(ldata),filname,myid,mblk2nd,
     .  maxbl,1)
c
  807 continue
c
c*****************************
c      embedded grid boundares
c*****************************
c
      igrid = igridg(nbl)
      iem = iemg(igrid)
      if (iem.gt.0) then
c
         jc  = jdimg(nblc)
         kc  = kdimg(nblc)
         ic  = idimg(nblc)
c
         nsi = (idim-1)/(ie-is)
c
         if (nsi.eq.2) then
c
c          full coarsening in i-direction
c
           if (jbcinfo(nbl,1,1,1).eq.21) then
              i21       = jbcinfo(nbl,1,2,1)
              i22       = jbcinfo(nbl,1,3,1)
              k21       = jbcinfo(nbl,1,4,1)
              k22       = jbcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1351)i21,i22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1557)js,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1557)js,nblc,is,ie,ks,ke
              end if
           end if
           if (jbcinfo(nbl,1,1,2).eq.21) then
              i21       = jbcinfo(nbl,1,2,2)
              i22       = jbcinfo(nbl,1,3,2)
              k21       = jbcinfo(nbl,1,4,2)
              k22       = jbcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1352)i21,i22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1557)je,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1557)je,nblc,is,ie,ks,ke
              end if
           end if
           if (kbcinfo(nbl,1,1,1).eq.21) then
              i21       = kbcinfo(nbl,1,2,1)
              i22       = kbcinfo(nbl,1,3,1)
              j21       = kbcinfo(nbl,1,4,1)
              j22       = kbcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1353)i21,i22,j21,j22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1558)ks,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1558)ks,nblc,is,ie,js,je
              end if
           end if
           if (kbcinfo(nbl,1,1,2).eq.21) then
              i21       = kbcinfo(nbl,1,2,2)
              i22       = kbcinfo(nbl,1,3,2)
              j21       = kbcinfo(nbl,1,4,2)
              j22       = kbcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1354)i21,i22,j21,j22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1558)ke,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1558)ke,nblc,is,ie,js,je
              end if
           end if
           if (ibcinfo(nbl,1,1,1).eq.21) then
              j21       = ibcinfo(nbl,1,2,1)
              j22       = ibcinfo(nbl,1,3,1)
              k21       = ibcinfo(nbl,1,4,1)
              k22       = ibcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1355)j21,j22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1559)is,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1559)is,nblc,is,ie,js,je
              end if
           end if
           if (ibcinfo(nbl,1,1,2).eq.21) then
              j21       = ibcinfo(nbl,1,2,2)
              j22       = ibcinfo(nbl,1,3,2)
              k21       = ibcinfo(nbl,1,4,2)
              k22       = ibcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1356)j21,j22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1559)ie,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1559)ie,nblc,is,ie,js,je
              end if
           end if
c
         else if (nsi.eq.1) then
c
c          semi coarsening in i-direction
c
           if (jbcinfo(nbl,1,1,1).eq.21) then
              i21       = jbcinfo(nbl,1,2,1)
              i22       = jbcinfo(nbl,1,3,1)
              k21       = jbcinfo(nbl,1,4,1)
              k22       = jbcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1451)i21,i22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1557)js,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1557)js,nblc,is,ie,ks,ke
              end if
           end if
           if (jbcinfo(nbl,1,1,2).eq.21) then
              i21       = jbcinfo(nbl,1,2,2)
              i22       = jbcinfo(nbl,1,3,2)
              k21       = jbcinfo(nbl,1,4,2)
              k22       = jbcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1452)i21,i22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1557)je,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1557)je,nblc,is,ie,ks,ke
              end if
           end if
           if (kbcinfo(nbl,1,1,1).eq.21) then
              i21       = kbcinfo(nbl,1,2,1)
              i22       = kbcinfo(nbl,1,3,1)
              j21       = kbcinfo(nbl,1,4,1)
              j22       = kbcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1453)i21,i22,j21,j22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1558)ks,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1558)ks,nblc,is,ie,js,je
              end if
           end if
           if (kbcinfo(nbl,1,1,2).eq.21) then
              i21       = kbcinfo(nbl,1,2,2)
              i22       = kbcinfo(nbl,1,3,2)
              j21       = kbcinfo(nbl,1,4,2)
              j22       = kbcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1454)i21,i22,j21,j22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1558)ke,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1558)ke,nblc,is,ie,js,je
              end if
           end if
           if (ibcinfo(nbl,1,1,1).eq.21) then
              j21       = ibcinfo(nbl,1,2,1)
              j22       = ibcinfo(nbl,1,3,1)
              k21       = ibcinfo(nbl,1,4,1)
              k22       = ibcinfo(nbl,1,5,1)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1455)j21,j22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1559)is,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1559)is,nblc,is,ie,js,je
              end if
           end if
           if (ibcinfo(nbl,1,1,2).eq.21) then
              j21       = ibcinfo(nbl,1,2,2)
              j22       = ibcinfo(nbl,1,3,2)
              k21       = ibcinfo(nbl,1,4,2)
              k22       = ibcinfo(nbl,1,5,2)
              if (isklton.eq.1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1456)j21,j22,k21,k22
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),1559)ie,nblc
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),1559)ie,nblc,is,ie,js,je
              end if
           end if
c
         end if
      end if
c
 1351 format(' ','  j=   1  embedded, full coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1352 format(' ','  j=jdim  embedded, full coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1353 format(' ','  k=   1  embedded, full coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1354 format(' ','  k=kdim  embedded, full coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1355 format(' ','  i=   1  embedded, full coarsening - i  type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1356 format(' ','  i=idim  embedded, full coarsening - i  type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1451 format(' ','  j=   1  embedded, semi coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1452 format(' ','  j=jdim  embedded, semi coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1453 format(' ','  k=   1  embedded, semi coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1454 format(' ','  k=kdim  embedded, semi coarsening - i  type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1455 format(' ','  i=   1  embedded, semi coarsening - i  type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1456 format(' ','  i=idim  embedded, semi coarsening - i  type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1557 format(' ',9x,' connects to j=',i5,' of block',i5)
 1558 format(' ',9x,' connects to k=',i5,' of block',i5)
 1559 format(' ',9x,' connects to i=',i5,' of block',i5)
c
c*************************
c     1-1 block boundaries
c*************************
c
      if (abs(nbli).gt.0) then 
c
c     loop over all 1-1 interfaces
c
      do lcnt = is_blk(level),ie_blk(level)
         ic_blk = isav_blk(lcnt,4)
         if (nbl.eq.ic_blk) then
            n      = isav_blk(lcnt,1)
            it     = isav_blk(lcnt,2)
            ir     = isav_blk(lcnt,3)
            in_blk = isav_blk(lcnt,5)
            jface  = isav_blk(lcnt,6)
            jedge  = isav_blk(lcnt,7)
            jside  = isav_blk(lcnt,8)
            lwt    = isav_blk(lcnt,9)
            iedge  = isav_blk(lcnt,10)
            iss    = isav_blk(lcnt,11)
            ise    = isav_blk(lcnt,12)
            jss    = isav_blk(lcnt,13)
            jse    = isav_blk(lcnt,14)
            kss    = isav_blk(lcnt,15)
            kse    = isav_blk(lcnt,16)
            iti    = isav_blk(lcnt,17)
c
c           set dimensions of current blocks involved
c
            idimn = idimg(in_blk)
            jdimn = jdimg(in_blk)
            kdimn = kdimg(in_blk)
            idimc = idimg(ic_blk)
            jdimc = jdimg(ic_blk)
            kdimc = kdimg(ic_blk)
c
c           k= constant interface  (qk0)
c
            if (jside.eq.3) then
               if (isklton.eq.1) then
                  if (iedge.eq.1 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1645)iss,ise,jss,jse
                  end if
                  if (iedge.eq.2 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1646)iss,ise,jss,jse
                  end if
                  if (jface.eq.1 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1857)jedge,in_blk
                  end if
                  if (jface.eq.2 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1858)jedge,in_blk
                  end if
                  if (jface.eq.3 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1859)jedge,in_blk
                  end if
               end if
c
c           j = constant interface  (qj0)
c
            else if (jside.eq.4) then
                if (isklton.eq.1) then
                   if (iedge.eq.1 .and. nbl.eq.ic_blk) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1),1643)iss,ise,kss,kse
                   end if
                   if (iedge.eq.2 .and. nbl.eq.ic_blk) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1),1644)iss,ise,kss,kse
                   end if
                   if (jface.eq.1 .and. nbl.eq.ic_blk) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1),1857)jedge,in_blk
                   end if
                   if (jface.eq.2 .and. nbl.eq.ic_blk) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1),1858)jedge,in_blk
                   end if
                   if (jface.eq.3 .and. nbl.eq.ic_blk) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1),1859)jedge,in_blk
                   end if
               end if
c
c           i = constant interface  (qi0)
c
            else if (jside.eq.5) then
               if (isklton.eq.1) then
                  if (iedge.eq.1 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1641)jss,jse,kss,kse
                  end if
                  if (iedge.eq.2 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1642)jss,jse,kss,kse
                  end if
                  if (jface.eq.1 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1857)jedge,in_blk
                  end if
                  if (jface.eq.2 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1858) jedge,in_blk
                  end if
                  if (jface.eq.3 .and. nbl.eq.ic_blk) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),1859)jedge,in_blk
                  end if
               end if
c
            end if
c
            if (real(geom_miss(lcnt)).gt.real(epsc0)) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''           '',
     .         ''blocking check....geometric mismatch = '',e14.7)')
     .         real(geom_miss(lcnt))
            end if
c
         end if
c
         end do
c
      end if
c
 1641 format(' ','  i=   1  1-1 blocking                   type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1642 format(' ','  i=idim  1-1 blocking                   type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1643 format(' ','  j=   1  1-1 blocking                   type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1644 format(' ','  j=jdim  1-1 blocking                   type    0',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1645 format(' ','  k=   1  1-1 blocking                   type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1646 format(' ','  k=kdim  1-1 blocking                   type    0',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1857 format(' ',9x,' connects to i =',i5,' of block',i5)
 1858 format(' ',9x,' connects to j =',i5,' of block',i5)
 1859 format(' ',9x,' connects to k =',i5,' of block',i5)
c
c****************************
c     patched grid boundaries
c****************************
c
      if(abs(ninter).gt.0) then
c
c     loop over all patch interfaces
c
      do 8001 icheck=1,abs(ninter)
      lmax1 =  iindex(icheck,1)
      nblcc =  iindex(icheck,lmax1+2)
      if (nblcc.ne.nbl) go to 8001
c
c     set range of points requiring interpolation on "to" side of
c     patch interface
c
      j21       = iindex(icheck,2*lmax1+6)
      j22       = iindex(icheck,2*lmax1+7)
      k21       = iindex(icheck,2*lmax1+8)
      k22       = iindex(icheck,2*lmax1+9)
c
c     patch surface to be interpolated to is an i=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.1) then
          if (iindex(icheck,2*lmax1+3).eq.11) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),1741)j21,j22,k21,k22
          else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),1742)j21,j22,k21,k22
          end if
      end if
c
c     patch surface to be interpolated to is a j=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.2) then
         if (iindex(icheck,2*lmax1+3).eq.21) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1743)j21,j22,k21,k22
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1744)j21,j22,k21,k22
         end if
      end if
c
c     patch surface to be interpolated to is a k=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.3) then
         if (iindex(icheck,2*lmax1+3).eq.31) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1745)j21,j22,k21,k22
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1746)j21,j22,k21,k22
         end if
      end if
c
      do 1705 l=1,lmax1
      mbl = iindex(icheck,l+1)
      mtype  = iindex(icheck,l+lmax1+2)
c
c     patch surface to be interpolated from is an i=constant surface
c
      if (mtype/10.eq.1) then
         if (mtype.eq.11) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1957) 1,mbl
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1957) idimg(mbl),mbl
         end if
      end if
c
c     patch surface to be interpolated from is a j=constant surface
c
      if (mtype/10.eq.2) then
         if (mtype.eq.21) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1958) 1,mbl
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1958) jdimg(mbl),mbl
         end if
      end if

c
c     patch surface to be interpolated from is a k=constant surface
c
      if (mtype/10.eq.3) then
         if (mtype.eq.31) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1959) 1,mbl
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1959) kdimg(mbl),mbl
         end if
      end if
c
 1705 continue
c
 8001 continue
c
      end if
c
 1741 format(' ','  i=   1  patched-grid interpolation     type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1742 format(' ','  i=idim  patched-grid interpolation     type    0',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1743 format(' ','  j=   1  patched-grid interpolation     type    0',
     .       '  k=',i5,',',i5,'  i=',i5,',',i5)
 1744 format(' ','  j=jdim  patched-grid interpolation     type    0',
     .       '  k=',i5,',',i5,'  i=',i5,',',i5)
 1745 format(' ','  k=   1  patched-grid interpolation     type    0',
     .       '  j=',i5,',',i5,'  i=',i5,',',i5)
 1746 format(' ','  k=kdim  patched-grid interpolation     type    0',
     .       '  j=',i5,',',i5,'  i=',i5,',',i5)
 1957 format(' ',9x,' interpolated from i =',i5,' of block',i5)
 1958 format(' ',9x,' interpolated from j =',i5,' of block',i5)
 1959 format(' ',9x,' interpolated from k =',i5,' of block',i5)
c
c****************************
c     chimera grid boundaries
c****************************
c
      if (iovrlp(nbl).gt.0) then
c
         idim1 = idim-1
         jdim1 = jdim-1
         kdim1 = kdim-1
c
         lsta = lbg(nbl)
         lend = lsta-1
c
         if (ibpntsg(nbl,1).gt.0) then
            lend = lsta+ibpntsg(nbl,1)-1
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1600)
         end if
c
         if (ibpntsg(nbl,2).gt.0) then
            mp0   = 0
            mpdim = 0
            lsta = lend+1
            lend = lsta+ibpntsg(nbl,2)-1
            do l=lsta,lend
               mp = max(0,1+jjbg(l)-jdim1) + 1 + max(0,-jjbg(l))
               mpp = max(0,3-mp)
               mp0 = mp0 + min(1,mpp)
               mpp = max(0,mp-2)
               mpdim = mpdim + min(1,mpp)
            end do
            if (mp0.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1603)
            end if
            if (mpdim.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1604)
            end if
         end if
c
         if (ibpntsg(nbl,3).gt.0) then
            mp0   = 0
            mpdim = 0
            lsta = lend+1
            lend = lsta+ibpntsg(nbl,3)-1
            do l=lsta,lend
               mp = max(0,1+kkbg(l)-kdim1) + 1 + max(0,-kkbg(l))
               mpp = max(0,3-mp)
               mp0 = mp0 + min(1,mpp)
               mpp = max(0,mp-2)
               mpdim = mpdim + min(1,mpp)
            end do
            if (mp0.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1605)
            end if
            if (mpdim.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1606)
            end if
         end if
c
         if (ibpntsg(nbl,4).gt.0) then
            mp0   = 0
            mpdim = 0
            lsta = lend+1
            lend = lsta+ibpntsg(nbl,4)-1
            do l=lsta,lend
               mp = max(0,1+iibg(l)-idim1) + 1 + max(0,-iibg(l))
               mpp = max(0,3-mp)
               mp0 = mp0 + min(1,mpp)
               mpp = max(0,mp-2)
               mpdim = mpdim + min(1,mpp)
            end do
            if (mp0.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1601)
            end if
            if (mpdim.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1602)
            end if
         end if
c
      end if
c
 1600 format(' ','  hole    chimera grid interpolation     type    0')
 1601 format(' ','  i=   1  chimera grid interpolation     type    0')
 1602 format(' ','  i=idim  chimera grid interpolation     type    0')
 1603 format(' ','  j=   1  chimera grid interpolation     type    0')
 1604 format(' ','  j=jdim  chimera grid interpolation     type    0')
 1605 format(' ','  k=   1  chimera grid interpolation     type    0')
 1606 format(' ','  k=kdim  chimera grid interpolation     type    0')
c
      return
      end
      subroutine out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001) ista,iend,ksta,kend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2001) ista,iend,ksta,kend
      end if
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002) ista,iend,ksta,kend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2002) ista,iend,ksta,kend
      end if
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003) ista,iend,jsta,jend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2003) ista,iend,jsta,jend
      end if
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004) ista,iend,jsta,jend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2004) ista,iend,jsta,jend
      end if
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005) jsta,jend,ksta,kend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2005) jsta,jend,ksta,kend
      end if
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006) jsta,jend,ksta,kend
      if (noninflag .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2006) jsta,jend,ksta,kend
      end if
      end if
c
 1001 format(' ','  j=   1  freestream                     type 1000',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  freestream                     type 1000',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  freestream                     type 1000',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  freestream                     type 1000',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  freestream                     type 1000',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  freestream                     type 1000',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
 2001 format(' ','  j=   1  NONINERTIAL freestream         type 1000',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  NONINERTIAL freestream         type 1000',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  NONINERTIAL freestream         type 1000',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  NONINERTIAL freestream         type 1000',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  NONINERTIAL freestream         type 1000',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  NONINERTIAL freestream         type 1000',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  symmetry plane                 type 1001',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  symmetry plane                 type 1001',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  symmetry plane                 type 1001',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  symmetry plane                 type 1001',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  symmetry plane                 type 1001',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  symmetry plane                 type 1001',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c 
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004) ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  1-pt extrapolation             type 1002',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  1-pt extrapolation             type 1002',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  1-pt extrapolation             type 1002',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  1-pt extrapolation             type 1002',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  1-pt extrapolation             type 1002',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  1-pt extrapolation             type 1002',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
      return 
      end
      subroutine out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,myid)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
c
      if (nface.eq.3) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1001)ista,iend,ksta,kend
        if (iipv .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write (bou(nou(1),1),3001)
        end if
      end if
c
      if (nface.eq.4) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1002)ista,iend,ksta,kend
        if (iipv .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write (bou(nou(1),1),3001)
        end if
      end if
c
      if (nface.eq.5) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1003)ista,iend,jsta,jend
        if (iipv .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write (bou(nou(1),1),3001)
        end if
      end if
c
      if (nface.eq.6) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1004)ista,iend,jsta,jend
        if (iipv .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write (bou(nou(1),1),3001)
        end if
      end if
c
      if (nface.eq.1) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1005)jsta,jend,ksta,kend
          if (iipv .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' i=1 bc1003 not appropriate for'',
     .       '' i2d=-1 (farfield point vortex correction)'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (nface.eq.2) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1006)jsta,jend,ksta,kend
          if (iipv .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' i=1 bc1003 not appropriate for'',
     .       '' i2d=-1 (farfield point vortex correction)'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
      end if
c
 1001 format(' ','  j=   1  characteristic inflow/outflow  type 1003',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  characteristic inflow/outflow  type 1003',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  characteristic inflow/outflow  type 1003',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  characteristic inflow/outflow  type 1003',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  characteristic inflow/outflow  type 1003',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  characteristic inflow/outflow  type 1003',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
 3001 format(' ',10x,'With farfield point vortex correction')
      return
      end
      subroutine out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,inormmom)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1001)ista,iend,ksta,kend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2001)ista,iend,ksta,kend
         end if
      end if
c
      if (nface.eq.4) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1002)ista,iend,ksta,kend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2002)ista,iend,ksta,kend
         end if
      end if
c
      if (nface.eq.5) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1003)ista,iend,jsta,jend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2003)ista,iend,jsta,jend
         end if
      end if
c
      if (nface.eq.6) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1004)ista,iend,jsta,jend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2004)ista,iend,jsta,jend
         end if
      end if
c
      if (nface.eq.1) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1005)jsta,jend,ksta,kend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2005)jsta,jend,ksta,kend
         end if
      end if
c
      if (nface.eq.2) then
         if (inormmom .eq. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1006)jsta,jend,ksta,kend
         else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),2006)jsta,jend,ksta,kend
         end if
      end if
c
 1001 format(' ','  j=   1  inviscid surface               type 1005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  inviscid surface               type 1005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  inviscid surface               type 1005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  inviscid surface               type 1005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  inviscid surface               type 1005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  inviscid surface               type 1005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2001 format(' ','  j=   1  inviscid surface w norm mom eq type 1006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  inviscid surface w norm mom eq type 1006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  inviscid surface w norm mom eq type 1006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  inviscid surface w norm mom eq type 1006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  inviscid surface w norm mom eq type 1006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  inviscid surface w norm mom eq type 1006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  tunnel inflow                  type 1008',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  tunnel inflow                  type 1008',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  tunnel inflow                  type 1008',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  tunnel inflow                  type 1008',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  tunnel inflow                  type 1008',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  tunnel inflow                  type 1008',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  singular axis - half plane     type 1011',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  singular axis - half plane     type 1011',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  singular axis - half plane     type 1011',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  singular axis - half plane     type 1011',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  singular axis - half plane     type 1011',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  singular axis - half plane     type 1011',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  singular axis - full plane     type 1012',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  singular axis - full plane     type 1012',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  singular axis - full plane     type 1012',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  singular axis - full plane     type 1012',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  singular axis - full plane     type 1012',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  singular axis - full plane     type 1012',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
      subroutine out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c 
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  singular axis - extrapolation  type 1013',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  singular axis - extrapolation  type 1013',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  singular axis - extrapolation  type 1013',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  singular axis - extrapolation  type 1013',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  singular axis - extrapolation  type 1013',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  singular axis - extrapolation  type 1013',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return 
      end
      subroutine out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,1)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
      end if
c
      if (nface.eq.4) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,1)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
      end if
c
      if (nface.eq.5) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,1)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
      end if
c
      if (nface.eq.6) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,1)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
c
      end if
c
      if (nface.eq.1) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,1)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
      end if
c
      if (nface.eq.2) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,1)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2002:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(pratio)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) filname
         end if
      end if
c
 1001 format(' ','  j=   1  set p, extrapolate rho,u,v,w   type 2002',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  set p, extrapolate rho,u,v,w   type 2002',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  set p, extrapolate rho,u,v,w   type 2002',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  set p, extrapolate rho,u,v,w   type 2002',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  set p, extrapolate rho,u,v,w   type 2002',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  set p, extrapolate rho,u,v,w   type 2002',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'P/Pinf     = ',f8.4)
 1008 format(11x,'P/Pinf set from file:')
 1009 format('           ',a60)
c
      return 
      end
      subroutine out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /maxiv/ ivmx
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,5
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,5
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,5
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,5
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,5
         do 5 ipp=1,2
         do 5 j=jsta,jend1
         jj = j-jsta+1
         do 5 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,5
         do 6 ipp=1,2
         do 6 j=jsta,jend1
         jj = j-jsta+1
         do 6 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2003:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  engine inflow data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2003...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            xme   = bcdata(1,1,1,1)
            pte   = bcdata(1,1,1,2)
            tte   = bcdata(1,1,1,3)
            alpe  = bcdata(1,1,1,4)
            betae = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(xme)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
 1001 format(' ','  j=   1  engine inflow                  type 2003',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  engine inflow                  type 2003',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  engine inflow                  type 2003',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  engine inflow                  type 2003',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  engine inflow                  type 2003',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  engine inflow                  type 2003',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'inlet mach     = ',f8.4)
 1008 format(11x,'Pt/Pinf        = ',f8.4)
 1009 format(11x,'Tt/Tinf        = ',f8.4)
 1010 format(11x,'alphe          = ',f8.4)
 1011 format(11x,'betae          = ',f8.4)
 1012 format(11x,'inflow data set from file:')
 1013 format('           ',a60)
 1014 format(11x,'turb1 (nondim) = ',e12.4)
 1015 format(11x,'turb2 (nondim) = ',e12.4)
c
      return
      end
      subroutine out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl,iout)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /wallfun/ iwf(3)
      common /wallvel_2034/ uub_2034,vvb_2034,wwb_2034,vvb_xaxisrot2034
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      iextra=0
      if (iout .eq. 1 .or. iout .eq. 2) iextra=1
c
      if (nface.eq.3) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,2+iextra
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1001)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1101)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1201)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1301)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1001)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1101)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1201)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1301)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1001)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1101)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1201)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1301)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1101)ista,iend,ksta,kend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1201)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),1301)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(2) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
      end if
c 
      if (nface.eq.4) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,2+iextra
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1002)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1102)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1202)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1302)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1002)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1102)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1202)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1302)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1002)ista,iend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1102)ista,iend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1202)ista,iend,ksta,kend
              else
                write(bou(nou(1),1),1302)ista,iend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1102)ista,iend,ksta,kend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1202)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),1302)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(2) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
      end if
c 
      if (nface.eq.5) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,2+iextra
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1003)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1103)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1203)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1303)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1003)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1103)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1203)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1303)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1003)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1103)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1203)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1303)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1003)ista,iend,jsta,jend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1103)ista,iend,jsta,jend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1203)ista,iend,jsta,jend
            else
              write(bou(nou(1),1),1303)ista,iend,jsta,jend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(3) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
      end if
c 
      if (nface.eq.6) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,2+iextra
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1004)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1104)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1204)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1304)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1004)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1104)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1204)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1304)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1004)ista,iend,jsta,jend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1104)ista,iend,jsta,jend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1204)ista,iend,jsta,jend
              else
                write(bou(nou(1),1),1304)ista,iend,jsta,jend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1004)ista,iend,jsta,jend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1104)ista,iend,jsta,jend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1204)ista,iend,jsta,jend
            else
              write(bou(nou(1),1),1304)ista,iend,jsta,jend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(3) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
      end if
c 
      if (nface.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,2+iextra
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1005)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1105)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1205)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1305)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1005)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1105)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1205)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1305)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1005)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1105)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1205)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1305)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1105)jsta,jend,ksta,kend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1205)jsta,jend,ksta,kend
            else
              write(bou(nou(1),1),1305)jsta,jend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(1) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
      end if
c 
      if (nface.eq.2) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,2+iextra
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),*)'  Stopping in bc2004:'
            else if (iout .eq. 1) then
              write(bou(nou(1),1),*)'  Stopping in bc2014:'
            else if (iout .eq. 2) then
              write(bou(nou(1),1),*)'  Stopping in bc2024:'
            else
              write(bou(nou(1),1),*)'  Stopping in bc2034:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            if (iout .eq. 1) then
              stopindex = bcdata(1,1,1,3)
            end if
            if (iout .eq. 2) then
              gammawall = bcdata(1,1,1,3)
            end if
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1006)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1106)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1206)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1306)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1006)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1106)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1206)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1306)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              if (iout .eq. 0) then
                write(bou(nou(1),1),1006)jsta,jend,ksta,kend
              else if (iout .eq. 1) then
                write(bou(nou(1),1),1106)jsta,jend,ksta,kend
              else if (iout .eq. 2) then
                write(bou(nou(1),1),1206)jsta,jend,ksta,kend
              else
                write(bou(nou(1),1),1306)jsta,jend,ksta,kend
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
              end if
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
            end if
            if (iout .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) int(stopindex)
            end if
            if (iout .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1113) gammawall
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
              write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            else if (iout .eq. 1) then
              write(bou(nou(1),1),1106)jsta,jend,ksta,kend
            else if (iout .eq. 2) then
              write(bou(nou(1),1),1206)jsta,jend,ksta,kend
            else
              write(bou(nou(1),1),1306)jsta,jend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1999)uub_2034,vvb_2034,wwb_2034,
     +           vvb_xaxisrot2034
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(1) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
      end if
c
 1001 format(' ','  j=   1  viscous wall                   type 2004',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  viscous wall                   type 2004',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  viscous wall                   type 2004',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  viscous wall                   type 2004',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  viscous wall                   type 2004',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  viscous wall                   type 2004',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1101 format(' ','  j=   1  viscous laminar wall           type 2014',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1102 format(' ','  j=jdim  viscous laminar wall           type 2014',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1103 format(' ','  k=   1  viscous laminar wall           type 2014',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1104 format(' ','  k=kdim  viscous laminar wall           type 2014',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1105 format(' ','  i=   1  viscous laminar wall           type 2014',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1106 format(' ','  i=idim  viscous laminar wall           type 2014',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1201 format(' ','  j=   1  viscous laminar wall           type 2024',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1202 format(' ','  j=jdim  viscous laminar wall           type 2024',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1203 format(' ','  k=   1  viscous laminar wall           type 2024',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1204 format(' ','  k=kdim  viscous laminar wall           type 2024',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1205 format(' ','  i=   1  viscous laminar wall           type 2024',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1206 format(' ','  i=idim  viscous laminar wall           type 2024',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1301 format(' ','  j=   1  viscous wall + imposed vel     type 2034',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1302 format(' ','  j=jdim  viscous wall + imposed vel     type 2034',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1303 format(' ','  k=   1  viscous wall + imposed vel     type 2034',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1304 format(' ','  k=kdim  viscous wall + imposed vel     type 2034',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1305 format(' ','  i=   1  viscous wall + imposed vel     type 2034',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1306 format(' ','  i=idim  viscous wall + imposed vel     type 2034',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'Tw/Tinf    = ',f8.4)
 1008 format(11x,'C_q        = ',f8.4)
 1013 format(11x,'stopindex  = ',i8)
 1113 format(11x,'intermtbc  = ',f8.4)
 1009 format(11x,'Tw/Tinf    = adiabatic wall')
 1010 format(11x,'Tw/Tinf    = stagnation')
 1011 format(11x,'Twtype, C_q set from file:')
 1012 format('           ',a60)
 1999 format(11x,'Imposed u,v,w,vrot_about_x_axis = ',4f10.5)
c
 3001 format(' ',10x,'Wall function employed in j-direction')
 3002 format(' ',10x,'Wall function employed in k-direction')
 3003 format(' ',10x,'Wall function employed in i-direction')
c
      return
      end
      subroutine out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,
     .                  jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,
     .                  maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,4
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 2737 l=2,4
         do 2737 i=ista,iend1
           ii = i-ista+1
           do 2737 k=ksta,kend1
             kk = k-ksta+1
             if (bcdata(kk,ii,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 2737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (kdim .ne. kdimp .or. idim .ne. idimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' idim,kdim,idimp,kdimp='',4i5)')
     .      idim,kdim,idimp,kdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. jdimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      ''dimension 2 in j'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      if (nface.eq.4) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,4
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 3737 l=2,4
         do 3737 i=ista,iend1
           ii = i-ista+1
           do 3737 k=ksta,kend1
             kk = k-ksta+1
             if (bcdata(kk,ii,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 3737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (kdim .ne. kdimp .or. idim .ne. idimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' idim,kdim,idimp,kdimp='',4i5)')
     .      idim,kdim,idimp,kdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. jdimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      ''dimension 2 in j'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      if (nface.eq.5) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,4
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 4737 l=2,4
         do 4737 i=ista,iend1
           ii = i-ista+1
           do 4737 j=jsta,jend1
             jj = j-jsta+1
             if (bcdata(jj,ii,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 4737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (jdim .ne. jdimp .or. idim .ne. idimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' idim,jdim,idimp,jdimp='',4i5)')
     .      idim,jdim,idimp,jdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. kdimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      '' dimension 2 in k'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      if (nface.eq.6) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,4
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 5737 l=2,4
         do 5737 i=ista,iend1
           ii = i-ista+1
           do 5737 j=jsta,jend1
             jj = j-jsta+1
             if (bcdata(jj,ii,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 5737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (jdim .ne. jdimp .or. idim .ne. idimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' idim,jdim,idimp,jdimp='',4i5)')
     .      idim,jdim,idimp,jdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. kdimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      '' dimension 2 in k'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      if (nface.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,4
         do 5 ipp=1,2
         do 5 j=jsta,jend1
         jj = j-jsta+1
         do 5 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 6737 l=2,4
         do 6737 k=ksta,kend1
           kk = k-ksta+1
           do 6737 j=jsta,jend1
             jj = j-jsta+1
             if (bcdata(jj,kk,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 6737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (jdim .ne. jdimp .or. kdim .ne. kdimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' jdim,kdim,jdimp,kdimp='',4i5)') 
     .      jdim,kdim,jdimp,kdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. idimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      '' dimension 2 in i'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      if (nface.eq.2) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,4
         do 6 ipp=1,2
         do 6 j=jsta,jend1
         jj = j-jsta+1
         do 6 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2005:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  periodic BC data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            blnb   = bcdata(1,1,1,1)
            dthtx  = bcdata(1,1,1,2)
            dthty  = bcdata(1,1,1,3)
            dthtz  = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) nblp 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) real(dthtx)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1014) real(dthty)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1015) real(dthtz)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
         ip = 1
c
c   perform check on dthtx,y,z inputs:
         do 7737 l=2,4
         do 7737 k=ksta,kend1
           kk = k-ksta+1
           do 7737 j=jsta,jend1
             jj = j-jsta+1
             if (bcdata(jj,kk,ip,l) .ne. bcdata(1,1,ip,l)) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' dthtx,y,z currently must be'',
     .          '' const over whole face'')')
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
 7737    continue
         if (bcdata(1,1,ip,2).ne.0. .and. (bcdata(1,1,ip,3).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,3).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,4).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (bcdata(1,1,ip,4).ne.0. .and. (bcdata(1,1,ip,2).ne.0. .or.
     .       bcdata(1,1,ip,3).ne.0.) ) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'('' 2 of the 3 dtht values'',
     .       '' currently must = 0'')')
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c   perform check on periodic block:
         if (jdim .ne. jdimp .or. kdim .ne. kdimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' jdim,kdim,jdimp,kdimp='',4i5)')
     .      jdim,kdim,jdimp,kdimp
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Block you are periodic with'',
     .      '' currently must be of same 2 dimensions'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' (and orientation) on the face as'',
     .      '' the BC block'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if ((nbl .ne. nblp) .and. idimp .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Periodic block cannot be'',
     .      '' dimension 2 in i'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  (on ANY multigrid level)'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
 1001 format(' ','  j=   1  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  periodic - angular rotation    type 2005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  periodic - angular rotation    type 2005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'periodic with j=jdim, block',i5)
 1008 format(11x,'periodic with j=   1, block',i5)
 1009 format(11x,'periodic with k=kdim, block',i5)
 1010 format(11x,'periodic with k=   1, block',i5)
 1011 format(11x,'periodic with i=idim, block',i5)
 1012 format(11x,'periodic with i=   1, block',i5)
 1013 format(11x,'dthetax    = ',f8.4)
 1014 format(11x,'dthetay    = ',f8.4)
 1015 format(11x,'dthetaz    = ',f8.4)
c
 2001 format(' ','  j=   1  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  periodic - angular rotation    type 2005',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  periodic - angular rotation    type 2005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  periodic - angular rotation    type 2005',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2007 format(11x,'periodicity set from file:')
 2008 format('           ',a60)
c
      return
      end
      subroutine out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,nblc,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      lijk = abs(int(real(bcdata(1,1,1,3))))
      ldir = int(real(bcdata(1,1,1,3)))/abs(int(real(bcdata(1,1,1,3))))
c
      xfact = 1.0
      yfact = 1.0
      zfact = 1.0
      if (bcdata(1,1,1,4) .eq. 1) then
         xfact = 0.0
      else if (bcdata(1,1,1,4) .eq. 2) then
         yfact = 0.0
      else if (bcdata(1,1,1,4) .eq. 3) then
         zfact = 0.0
      else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)' stopping...bcdata(4) = ',
     .   bcdata(1,1,1,4),' should be +1, +2, or +3' 
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      if (nface.eq.3) then
c
         j = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,4
         do 1 ippp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(kk,ii,ippp,3)))) .gt. 3) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c 
         if (lijk .eq. 1 .or. lijk .eq. 3) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 1 or 3 on a j=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1001)ista,iend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4001)ista,iend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 1 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3005)
         end if
         if (lijk .eq. 1 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3006)
         end if
         if (lijk .eq. 3 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
         if (lijk .eq. 3 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3004)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ista.ne.1 .or. iend.ne.idim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
      if (nface.eq.4) then
c
         j = jdim - 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,2
         do 2 ippp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(kk,ii,ippp,3)))) .gt. 3) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (lijk .eq. 1 .or. lijk .eq. 3) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 1 or 3 on a j=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1002)ista,iend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4002)ista,iend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 1 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3005)
         end if
         if (lijk .eq. 1 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3006)
         end if
         if (lijk .eq. 3 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
         if (lijk .eq. 3 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3004)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ista.ne.1 .or. iend.ne.idim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
c
      if (nface.eq.5) then
c
         k = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,4
         do 3 ippp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(jj,ii,ippp,3)))) .gt. 3) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (lijk .eq. 1 .or. lijk .eq. 2) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 1 or 2 on a k=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1003)ista,iend,jsta,jend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4003)ista,iend,jsta,jend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 2 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
         if (lijk .eq. 2 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
         if (lijk .eq. 1 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3005)
         end if
         if (lijk .eq. 1 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3006)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ista.ne.1 .or. iend.ne.idim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
      if (nface.eq.6) then
c
         k = kdim - 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,4
         do 4 ippp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(jj,ii,ippp,3)))) .gt. 3) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (lijk .eq. 1 .or. lijk .eq. 2) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 1 or 2 on a k=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1004)ista,iend,jsta,jend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4004)ista,iend,jsta,jend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 2 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
         if (lijk .eq. 2 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
         if (lijk .eq. 1 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3005)
         end if
         if (lijk .eq. 1 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3006)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ista.ne.1 .or. iend.ne.idim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
      if (nface.eq.1) then
c
         i = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,4
         do 5 ippp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(jj,kk,ippp,3)))) .gt. 3) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (lijk .eq. 2 .or. lijk .eq. 3) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 2 or 3 on an i=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1005)jsta,jend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4005)jsta,jend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 2 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
         if (lijk .eq. 2 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
         if (lijk .eq. 3 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
         if (lijk .eq. 3 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3004)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (jsta.ne.1 .or. jend.ne.jdim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
      if (nface.eq.2) then
c
         i = idim - 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,4
         do 6 ippp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ippp,l)) .lt. -1.e10) iflg = 1
         if (int(abs(real(bcdata(jj,kk,ippp,3)))) .gt. 3) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2006:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  radial equilib. press. data',
     . ' incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (lijk .eq. 2 .or. lijk .eq. 3) then
c           lijk ok
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping...abs(bcdata(3)) must be',
     .                 ' 2 or 3 on an i=constant face'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            if (real(bcdata(1,1,1,1)) .le. 0) then 
               pratio = bcdata(1,1,1,2)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1006)jsta,jend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1007) real(pratio)
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4006)jsta,jend,ksta,kend
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),4007)nblc
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
         if (lijk .eq. 2 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
         if (lijk .eq. 2 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
         if (lijk .eq. 3 .and. ldir .eq. +1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
         if (lijk .eq. 3 .and. ldir .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3004)
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (jsta.ne.1 .or. jend.ne.jdim) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping...bc segment must',
     .         ' span the entire block face in the radial direction'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
 1001 format(' ','  j=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'P/Pinf     = ',f8.4)
c
 2001 format(' ','  j=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2007 format(11x,'P/Pinf set from file:')
 2008 format('           ',a60)
c
 3001 format(11x,'integrated in +j direction')
 3002 format(11x,'integrated in -j direction')
 3003 format(11x,'integrated in +k direction')
 3004 format(11x,'integrated in -k direction')
 3005 format(11x,'integrated in +i direction')
 3006 format(11x,'integrated in -i direction')
c
 4001 format(' ','  j=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 4002 format(' ','  j=jdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 4003 format(' ','  k=   1  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 4004 format(' ','  k=kdim  set p via radial equilibrium   type 2006',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 4005 format(' ','  i=   1  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 4006 format(' ','  i=idim  set p via radial equilibrium   type 2006',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 4007 format(11x,'continued from block',i5)
c
      return 
      end
      subroutine out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /maxiv/ ivmx
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,5
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c 
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c 
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,5
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,5
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,5
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,5
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,7)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,5
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2007:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  primative variable data ',
     .                         'incorrectly set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2007...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (filname.eq.'null') then
            rho1 = bcdata(1,1,1,1)
            u1   = bcdata(1,1,1,2)
            v1   = bcdata(1,1,1,3)
            w1   = bcdata(1,1,1,4)
            p1   = bcdata(1,1,1,5)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(rho1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(w1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(p1)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,7)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
 1001 format(' ','  j=   1  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  set all primative variables    type 2007',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  set all primative variables    type 2007',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'rho/rhoinf       = ',f8.4)
 1008 format(11x,'u/ainf           = ',f8.4)
 1009 format(11x,'v/ainf           = ',f8.4)
 1010 format(11x,'w/ainf           = ',f8.4)
 1011 format(11x,'p/rhoinf/ainf**2 = ',f8.4)
 1014 format(11x,'turb1 (nondim)   = ',e12.4)
 1015 format(11x,'turb2 (nondim)   = ',e12.4)
c
 2001 format(' ','  j=   1  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  set all primative variables    type 2007',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  set all primative variables    type 2007',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  set all primative variables    type 2007',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2007 format(11x,'primative variables set from file:')
 2008 format('           ',a60)
c
      return 
      end
      subroutine out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl,iout)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,4
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c 
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3001)ista,iend,ksta,kend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4001)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,4
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3002)ista,iend,ksta,kend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4002)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,4
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3003)ista,iend,jsta,jend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4003)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,4
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3004)ista,iend,jsta,jend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4004)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,4
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3005)jsta,jend,ksta,kend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4005)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,4
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),*)'  Stopping in bc2008:'
            else if (iout .eq. 1) then
            write(bou(nou(1),1),*)'  Stopping in bc2018:'
            else if (iout .eq. 2) then
            write(bou(nou(1),1),*)'  Stopping in bc2028:'
            else if (iout .eq. 3) then
            write(bou(nou(1),1),*)'  Stopping in bc2038:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            val1  = bcdata(1,1,1,1)
            u1    = bcdata(1,1,1,2)
            v1    = bcdata(1,1,1,3)
            w1    = bcdata(1,1,1,4)
            if (iout .eq. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            else if (iout .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2011) real(w1)
            else if (iout .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3011) real(w1)
            else if (iout .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(val1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(u1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(v1)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(w1)
            end if
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iout .eq. 0) then
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            else if (iout .eq. 1) then
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            else if (iout .eq. 2) then
            write(bou(nou(1),1),3006)jsta,jend,ksta,kend
            else if (iout .eq. 3) then
            write(bou(nou(1),1),4006)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
 1001 format(' ','  j=   1  set rho,u,v,w, extrapolate p   type 2008',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  set rho,u,v,w, extrapolate p   type 2008',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  set rho,u,v,w, extrapolate p   type 2008',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  set rho,u,v,w, extrapolate p   type 2008',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  set rho,u,v,w, extrapolate p   type 2008',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  set rho,u,v,w, extrapolate p   type 2008',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1008 format(11x,'rho/rhoinf     = ',f8.4)
 1009 format(11x,'u/ainf         = ',f8.4)
 1010 format(11x,'v/ainf         = ',f8.4)
 1011 format(11x,'w/ainf         = ',f8.4)
 1012 format(11x,'inflow data set from file:')
 1013 format('           ',a60)
 1014 format(11x,'turb1 (nondim) = ',e12.4)
 1015 format(11x,'turb2 (nondim) = ',e12.4)
 2001 format(' ','  j=   1  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  set T,rhou,rhov,rhow, extrap p type 2018',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2008 format(11x,'T/Tinf         = ',f8.4)
 2009 format(11x,'rhou/rhoinfainf= ',f8.4)
 2010 format(11x,'rhov/rhoinfainf= ',f8.4)
 2011 format(11x,'rhow/rhoinfainf= ',f8.4)
 3001 format(' ','  j=   1 set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 3002 format(' ','  j=jdim set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 3003 format(' ','  k=   1 set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 3004 format(' ','  k=kdim set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 3005 format(' ','  i=   1 set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 3006 format(' ','  i=idim set f,rhou,rhov,rhow, ext rho,p type 2028',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 3008 format(11x,'freq*lref/aref    = ',f8.4)
 3009 format(11x,'rhoumax/rhoinfainf= ',f8.4)
 3010 format(11x,'rhovmax/rhoinfainf= ',f8.4)
 3011 format(11x,'rhowmax/rhoinfainf= ',f8.4)
 4001 format(' ','  j=   1  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 4002 format(' ','  j=jdim  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 4003 format(' ','  k=   1  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 4004 format(' ','  k=kdim  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 4005 format(' ','  i=   1  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 4006 format(' ','  i=idim  set rho,u,v,w+rndm, extrap p   type 2038',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return 
      end
      subroutine out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl,iflgg)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /maxiv/ ivmx
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,4
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2001)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2001)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,4
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2002)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2002)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,4
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2003)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2003)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,4
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2004)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2004)ista,iend,jsta,jend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,4
         do 5 ipp=1,2
         do 5 j=jsta,jend1
         jj = j-jsta+1
         do 5 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,5)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,6)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,4
         do 6 ipp=1,2
         do 6 j=jsta,jend1
         jj = j-jsta+1
         do 6 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2009/2010:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2009/2010...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            alpe  = bcdata(1,1,1,3)
            betae = bcdata(1,1,1,4)
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(alpe)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(betae)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,5)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,6)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            if (iflgg .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
 1001 format(' ','  j=   1  nozzle total BC                type 2009',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  nozzle total BC                type 2009',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  nozzle total BC                type 2009',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  nozzle total BC                type 2009',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  nozzle total BC                type 2009',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  nozzle total BC                type 2009',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2001 format(' ','  j=   1  nozzle total BC                type 2010',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2002 format(' ','  j=jdim  nozzle total BC                type 2010',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 2003 format(' ','  k=   1  nozzle total BC                type 2010',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2004 format(' ','  k=kdim  nozzle total BC                type 2010',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 2005 format(' ','  i=   1  nozzle total BC                type 2010',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 2006 format(' ','  i=idim  nozzle total BC                type 2010',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1008 format(11x,'Pt/Pinf        = ',f8.4)
 1009 format(11x,'Tt/Tinf        = ',f8.4)
 1010 format(11x,'alphe          = ',f8.4)
 1011 format(11x,'betae          = ',f8.4)
 1012 format(11x,'inflow data set from file:')
 1013 format('           ',a60)
 1014 format(11x,'turb1 (nondim) = ',e12.4)
 1015 format(11x,'turb2 (nondim) = ',e12.4)
c
      return
      end
      subroutine out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /maxiv/ ivmx
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,2
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,2
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,2
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,2
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,2
         do 5 ipp=1,2
         do 5 j=jsta,jend1
         jj = j-jsta+1
         do 5 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,3)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,4)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,2
         do 6 ipp=1,2
         do 6 j=jsta,jend1
         jj = j-jsta+1
         do 6 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2019:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  nozzle total BC data incorrectly',
     .      ' set'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (itrflg1.gt.0 .or. itrflg2.gt.0) then
            if (ivmx.gt.5 .and. itrflg1*itrflg2.eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  stopping in bc2019...must',
     .         ' set both turbulence'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'  quantities for 2-eq turb.',
     .         ' models...only one is set'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
         if (filname.eq.'null') then
            pte   = bcdata(1,1,1,1)
            tte   = bcdata(1,1,1,2)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pte)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(tte)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,3)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,4)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
      end if
c
 1001 format(' ','  j=   1  OVERFLOW-type nozzle total BC  type 2019',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  OVERFLOW-type nozzle total BC  type 2019',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  OVERFLOW-type nozzle total BC  type 2019',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  OVERFLOW-type nozzle total BC  type 2019',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  OVERFLOW-type nozzle total BC  type 2019',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  OVERFLOW-type nozzle total BC  type 2019',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1008 format(12x,'Pt/Ptinf        = ',f8.4)
 1009 format(12x,'Tt/Ttinf        = ',f8.4)
 1012 format(11x,'inflow data set from file:')
 1013 format('           ',a60)
 1014 format(11x,'turb1 (nondim) = ',e12.4)
 1015 format(11x,'turb2 (nondim) = ',e12.4)
c
      return
      end
      subroutine out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl,iout)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /wallfun/ iwf(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,2+iout
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(kk,ii,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(2) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
      end if
c 
      if (nface.eq.4) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,2+iout
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(kk,ii,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(2) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3001)
         end if
      end if
c 
      if (nface.eq.5) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,2+iout
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(jj,ii,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,jsta,jend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(3) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
      end if
c 
      if (nface.eq.6) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,2+iout
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(jj,ii,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(3) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3002)
         end if
      end if
c 
      if (nface.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,2+iout
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(jj,kk,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(1) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
      end if
c 
      if (nface.eq.2) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,2+iout
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
         if ((real(bcdata(jj,kk,ipp,l)).ne.0.0).and.
     +       (real(dt).lt.0.).and.(l.eq.3)) iflg = 2
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2016:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  viscous wall data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
           if (iflg.eq.2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2016:'
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
     +        '  unsteady jet not allowed with dt < 0'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
           end if

         if (filname.eq.'null') then
            ctemp = bcdata(1,1,1,1)
            cq    = bcdata(1,1,1,2)
            cqu   = bcdata(1,1,1,3)
            sjetx = bcdata(1,1,1,4)
            sjety = bcdata(1,1,1,5)
            sjetz = bcdata(1,1,1,6)
            rfreq = bcdata(1,1,1,7)
            if (real(ctemp) .gt. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007) real(ctemp)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else if (ctemp .eq. 0.) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1009)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)ista,iend,ksta,kend
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1010)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1008) real(cq)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1013) real(cqu)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1014) real(sjetx)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1015) real(sjety)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1016) real(sjetz)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1017) real(rfreq)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012) filname
         end if
         if (iwf(1) .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),3003)
         end if
      end if
c
 1001 format(' ','  j=   1  viscous wall                   type 2016',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  viscous wall                   type 2016',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  viscous wall                   type 2016',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  viscous wall                   type 2016',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  viscous wall                   type 2016',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  viscous wall                   type 2016',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'Tw/Tinf    = ',f8.4)
 1008 format(11x,'cq         = ',f8.4)
 1009 format(11x,'Tw/Tinf    = adiabatic wall')
 1010 format(11x,'Tw/Tinf    = stagnation')
 1011 format(11x,'Twtype, C_q set from file:')
 1012 format('           ',a60)
 1013 format(11x,'cqu        = ',f8.4)
 1014 format(11x,'sjetx      = ',f8.4)
 1015 format(11x,'sjetx      = ',f8.4)
 1016 format(11x,'sjetx      = ',f8.4)
 1017 format(11x,'rfreq      = ',f8.4)
c
 3001 format(' ',10x,'Wall function employed in j-direction')
 3002 format(' ',10x,'Wall function employed in k-direction')
 3003 format(' ',10x,'Wall function employed in i-direction')
c
      return
      end
      subroutine out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (nface.eq.3) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,9
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do i=ista,iend1
         ii = i-ista+1
         do k=ksta,kend1
         kk = k-ksta+1
           if(real(bcdata(kk,ii,ipp,4)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,5)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(kk,ii,ipp,7)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,8)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c 
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.4) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,9
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do i=ista,iend1
         ii = i-ista+1
         do k=ksta,kend1
         kk = k-ksta+1
           if(real(bcdata(kk,ii,ipp,4)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,5)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(kk,ii,ipp,7)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,8)) .eq. 0. .and.
     +        real(bcdata(kk,ii,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.5) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,9
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do i=ista,iend1
         ii = i-ista+1
         do j=jsta,jend1
         jj = j-jsta+1
           if(real(bcdata(jj,ii,ipp,4)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,5)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(jj,ii,ipp,7)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,8)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.6) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,9
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do i=ista,iend1
         ii = i-ista+1
         do j=jsta,jend1
         jj = j-jsta+1
           if(real(bcdata(jj,ii,ipp,4)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,5)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(jj,ii,ipp,7)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,8)) .eq. 0. .and.
     +        real(bcdata(jj,ii,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.1) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,9
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do k=ksta,kend1
         kk = k-ksta+1
         do j=jsta,jend1
         jj = j-jsta+1
           if(real(bcdata(jj,kk,ipp,4)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,5)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(jj,kk,ipp,7)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,8)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
      if (nface.eq.2) then
c
c        check to see if turbulence data is input (itrflg1/2 = 1) or
c        if freestream values are to be used (itrflg1/2 = 0); the check
c        assumes if the first point has been set, all points have been
c
         ipp     = 1
         itrflg1 = 0
         itrflg2 = 0
         if (real(bcdata(1,1,ipp,10)) .ge. 0.) itrflg1 = 1
         if (real(bcdata(1,1,ipp,11)) .ge. 0.) itrflg2 = 1
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,9
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  variable data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         iflg = 0
         do ipp=1,2
         do k=ksta,kend1
         kk = k-ksta+1
         do j=jsta,jend1
         jj = j-jsta+1
           if(real(bcdata(jj,kk,ipp,4)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,5)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,6)) .eq. 0.) iflg = 1
           if(real(bcdata(jj,kk,ipp,7)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,8)) .eq. 0. .and.
     +        real(bcdata(jj,kk,ipp,9)) .eq. 0.) iflg = 1
         enddo
         enddo
         enddo
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2026:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  normal vectors cannot be zero'
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            vmag       = bcdata(1,1,1,1)
            rfreq      = bcdata(1,1,1,2)
            sideangj   = bcdata(1,1,1,3)
            sxa        = bcdata(1,1,1,4)
            sya        = bcdata(1,1,1,5)
            sza        = bcdata(1,1,1,6)
            sxb        = bcdata(1,1,1,7)
            syb        = bcdata(1,1,1,8)
            szb        = bcdata(1,1,1,9)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(vmag)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(rfreq)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) real(sideangj)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1011) real(sxa),real(sya),real(sza)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007) real(sxb),real(syb),real(szb)
            if (itrflg1.gt.0) then
               t1 = bcdata(1,1,1,10)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1014) real(t1)
            end if
            if (itrflg2.gt.0) then
               t2 = bcdata(1,1,1,11)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1015) real(t2)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1012)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1013) filname
         end if
c
      end if
c
 1001 format(' ','  j=   1  sweeping jet   type 2026',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  sweeping jet   type 2026',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  sweeping jet   type 2026',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  sweeping jet   type 2026',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  sweeping jet   type 2026',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  sweeping jet   type 2026',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1008 format(11x,'vmag/ainf      = ',f8.4)
 1009 format(11x,'rfreq*lref/ainf= ',f8.4)
 1010 format(11x,'sideangj (deg) = ',f8.4)
 1011 format(11x,'sxa,sya,sza    = ',3f8.4)
 1007 format(11x,'sxb,syb,szb    = ',3f8.4)
 1012 format(11x,'inflow data set from file:')
 1013 format('           ',a60)
 1014 format(11x,'turb1 (nondim) = ',e12.4)
 1015 format(11x,'turb2 (nondim) = ',e12.4)
c
      return 
      end
      subroutine out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,
     .                  filname,myid,mblk2nd,maxbl,iflag)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
#   ifdef CMPLX
      complex lref
#   else
      real lref
#   endif
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      twopi = 8.0*atan(1.0)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (iflag .eq. 0) then
        l_set=4
      else
        l_set=5
      end if
c
      if (nface.eq.3) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,l_set
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c 
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1001)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3001)ista,iend,ksta,kend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) 
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2001)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4001)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c 
      if (nface.eq.4) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,l_set
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1002)ista,iend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3002)ista,iend,ksta,kend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2002)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4002)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.5) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,l_set
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1003)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3003)ista,iend,jsta,jend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2003)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4003)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.6) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,l_set
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1004)ista,iend,jsta,jend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3004)ista,iend,jsta,jend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2004)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4004)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,l_set
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3005)jsta,jend,ksta,kend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2005)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4005)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      if (nface.eq.2) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,l_set
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2102:'
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)'  Stopping in bc2103:'
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),3006)jsta,jend,ksta,kend
            end if
            if (iflag .eq. 0) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1007)
            else
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1012)
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) 
     .            real(twopi)*real(rfreqp)/real(lref)
            if (iflag .eq. 1) then
              phioff = bcdata(1,1,1,5)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),1011) real(phioff)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            if (iflag .eq. 0) then
              write(bou(nou(1),1),2006)ista,iend,ksta,kend
            else
              write(bou(nou(1),1),4006)ista,iend,ksta,kend
            end if
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
 1001 format(' ','  j=   1  time-varying pressure          type 2102',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  time-varying pressure          type 2102',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  time-varying pressure          type 2102',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  time-varying pressure          type 2102',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  time-varying pressure          type 2102',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  time-varying pressure          type 2102',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1007 format(11x,'P(t)/Pinf  = P1 + eps*sin(k*t)')
 1008 format(11x,'P1         = ',f8.4)
 1009 format(11x,'eps        = ',f8.4)
 1010 format(11x,'k          = ',f8.4)
 1011 format(11x,'phioff(deg)= ',f8.4)
 1012 format(11x,'P(t)/Pinf  = P1 + eps*sin(k*t+phioff)')
c
 2001 format(' ','  j=   1  time-varying pressure          type 2102')
 2002 format(' ','  j=jdim  time-varying pressure          type 2102')
 2003 format(' ','  k=   1  time-varying pressure          type 2102')
 2004 format(' ','  k=kdim  time-varying pressure          type 2102')
 2005 format(' ','  i=   1  time-varying pressure          type 2102')
 2006 format(' ','  i=idim  time-varying pressure          type 2102')
 2007 format(11x,'P(t)/Pinf set from file:')
 2008 format('           ',a60)
c
 3001 format(' ','  j=   1  time-varying pressure          type 2103',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 3002 format(' ','  j=jdim  time-varying pressure          type 2103',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 3003 format(' ','  k=   1  time-varying pressure          type 2103',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 3004 format(' ','  k=kdim  time-varying pressure          type 2103',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 3005 format(' ','  i=   1  time-varying pressure          type 2103',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 3006 format(' ','  i=idim  time-varying pressure          type 2103',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
 4001 format(' ','  j=   1  time-varying pressure          type 2103')
 4002 format(' ','  j=jdim  time-varying pressure          type 2103')
 4003 format(' ','  k=   1  time-varying pressure          type 2103')
 4004 format(' ','  k=kdim  time-varying pressure          type 2103')
 4005 format(' ','  i=   1  time-varying pressure          type 2103')
 4006 format(' ','  i=idim  time-varying pressure          type 2103')
c
      return 
      end
      subroutine out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,
     .                  nface,nou,bou,nbuf,ibufdim)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      if (nface.eq.3) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1001) ista,iend,ksta,kend
      end if
c
      if (nface.eq.4) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1002) ista,iend,ksta,kend
      end if
c
      if (nface.eq.5) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1003) ista,iend,jsta,jend
      end if
c
      if (nface.eq.6) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1004) ista,iend,jsta,jend
      end if
c
      if (nface.eq.1) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1005) jsta,jend,ksta,kend
      end if
c
      if (nface.eq.2) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1006) jsta,jend,ksta,kend
      end if
c
 1001 format(' ','  j=   1  exact soln (MMS)               type 9999',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1002 format(' ','  j=jdim  exact soln (MMS)               type 9999',
     .       '  i=',i5,',',i5,'  k=',i5,',',i5)
 1003 format(' ','  k=   1  exact soln (MMS)               type 9999',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1004 format(' ','  k=kdim  exact soln (MMS)               type 9999',
     .       '  i=',i5,',',i5,'  j=',i5,',',i5)
 1005 format(' ','  i=   1  exact soln (MMS)               type 9999',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
 1006 format(' ','  i=idim  exact soln (MMS)               type 9999',
     .       '  j=',i5,',',i5,'  k=',i5,',',i5)
c
      return
      end
