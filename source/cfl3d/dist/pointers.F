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
      subroutine pointers (lw,lw2,maxl,lembed,nstart,nwork,mwork,maxbl,
     .                     maxgr,maxseg,lwdat,levelg,igridg,iviscg,
     .                     idimg,jdimg,kdimg,nblcg,itrans,irotat,idefrm,
     .                     nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .                     ibcinfo,jbcinfo,kbcinfo,ngrid,ncgg,nblg,
     .                     iemg,nblock,myhost,myid,mblk2nd,nou,bou,nbuf,
     .                     ibufdim,nblfine,ilamlog,jlamlog,
     .                     klamlog,ilamhig,jlamhig,klamhig,idegg,iwfg,
     .                     idiagg,iflimg,ifdsg,rkap0g,jsg,ksg,isg,jeg,
     .                     keg,ieg,memblock,icall,nmds,maxaes,mpihost)
#   ifdef CMPLX
#   else
      use module_kwstm, only:kws_get_nummem
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To set up pointers for main storage array w
c
c     mpihost...flag to determine memory requirements for the host node:
c               = 0, get memory requirements for a sequential host
c               = 1, get memory requirements for a parallel host
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#endif
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension nblfine(maxbl),mblk2nd(maxbl),memblock(maxbl)
      dimension lw(65,maxbl),lwdat(maxbl,maxseg,6),lw2(43,maxbl)
      dimension levelg(maxbl),igridg(maxbl),iviscg(maxbl,3),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),nblcg(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl), 
     .          nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr)
      dimension ilamlog(maxbl),ilamhig(maxbl),jlamlog(maxbl),
     .          jlamhig(maxbl),klamlog(maxbl),klamhig(maxbl),
     .          iwfg(maxbl,3),idegg(maxbl,3),idiagg(maxbl,3),
     .          iflimg(maxbl,3),ifdsg(maxbl,3),rkap0g(maxbl,3),
     .          jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /elastic/ ndefrm,naesrf
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /des/ cdes,ides,cddes
      common /plot3dtyp/ ifunct
      common /memry/ lowmem_ux
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
c
      if (icall .eq. 0) then
         ierrflg = -99
      else
         ierrflg = -1
      end if
c
      lglobal = 1
      lembed  = 0
      do 103 igrid=1,ngrid
         lglobal = max0(lglobal,ncgg(igrid)+1)
         lembed  = max0(lembed,iemg(igrid))
  103 continue
      maxl    = lglobal+lembed
c
c      maxl    : level of finest grids
c      lembed  : levels of embedment above finest global grids
c      lglobal : level of finest global grids
c            1 : level of coarsest global grids
c
c     starting locations for data in w are stored in 
c     lw(i,n) where i is type of data and n is block number
c
c     Note: for the parallel code, the pointers for the host node
c     are sized based on the corresponding size in the largest block,
c     and the pointers for all blocks start at 1 (i.e. lw(1,nbl) = 1
c     for all nbl).
c
c     i correspondence to data is as follows
c     i=1  q
c       2  qj0
c       3  qk0
c       4  qi0
c       5  sj
c       6  sk
c       7  si
c       8  vol
c       9  dtj
c      10  x
c      11  y
c      12  z
c      13  vist3d
c      14  snk0 or smin
c      15  snj0/sni0 or xjb
c      16  q1
c      17  qr
c      18  blank
c      19  xib or tursav
c      20  sig
c      21  sqtq
c      22  g
c      23  tj0
c      24  tk0
c      25  ti0
c      26  xkb
c      27  blnum
c      28  vj0
c      29  vk0
c      30  vi0
c      31  bcj
c      32  bck
c      33  bci
c      34  qc0
c      35  dqc0
c      36  xtbj
c      37  xtbk
c      38  xtbi
c      39  atbj
c      40  atbk
c      41  atbi
c      42  bcdataj
c      43  bcdatak
c      44  bcdatai
c      45  xib2
c      46  qavg
c      47  ux
c      48  cmuv
c      49  volj0
c      50  volk0
c      51  voli0
c      52  xmdj
c      53  xmdk
c      54  xmdi
c      55  velg
c      56  xnm2
c      57  ynm2
c      58  znm2
c      59  deltj
c      60  deltk
c      61  delti
c      62  xnm1
c      63  ynm1
c      64  znm1
c      65  q2avg
c
c
      if (myid.eq.myhost .and. mpihost.eq.1) then
c
         jmxx = jdimg(nblg(1))
         kmxx = kdimg(nblg(1))
         imxx = idimg(nblg(1))
         jmkmim    = jmxx*kmxx*imxx
         jm1km1im1 = (jmxx-1)*(kmxx-1)*(imxx-1)
         jmkmim1   = jmxx*kmxx*(imxx-1)
         kmim1     = kmxx*(imxx-1)
         jmim1     = jmxx*(imxx-1)
         jmkm      = jmxx*kmxx
         kmim      = kmxx*imxx
         jmim      = jmxx*imxx
c
         do igrid = 2,ngrid
            nbl = nblg(igrid)
            jmxx = jdimg(nblg(igrid))
            kmxx = kdimg(nblg(igrid))
            imxx = idimg(nblg(igrid))
            jmkmim_ig    = jmxx*kmxx*imxx
            jm1km1im1_ig = (jmxx-1)*(kmxx-1)*(imxx-1)
            jmkmim1_ig   = jmxx*kmxx*(imxx-1)
            kmim1_ig     = kmxx*(imxx-1)
            jmim1_ig     = jmxx*(imxx-1)
            jmkm_ig      = jmxx*kmxx
            kmim_ig      = kmxx*imxx
            jmim_ig      = jmxx*imxx
            if (jmkmim_ig    .gt. jmkmim)    jmkmim    = jmkmim_ig
            if (jm1km1im1_ig .gt. jm1km1im1) jm1km1im1 = jm1km1im1_ig
            if (jmkmim1_ig   .gt. jmkmim1)   jmkmim1   = jmkmim1_ig
            if (kmim1_ig     .gt. kmim1)     kmim1     = kmim1_ig
            if (jmim1_ig     .gt. jmim1)     jmim1     = jmim1_ig
            if (jmkm_ig      .gt. jmkm)      jmkm      = jmkm_ig
            if (kmim_ig      .gt. kmim)      kmim      = kmim_ig
            if (jmim_ig      .gt. jmim)      jmim      = jmim_ig
         end do
c
c        determine max boundary data array dimensions
c
         jbcmax = 0
         kbcmax = 0
         ibcmax = 0
         do igrid = 1,ngrid
            nbl = nblg(igrid)
c
            jbc = 0
            do nseg=1,nbcj0(nbl)
              if(jbcinfo(nbl,nseg,7,1).ne.0) jbc = jbc +
     .          (jbcinfo(nbl,nseg,3,1)-jbcinfo(nbl,nseg,2,1))*
     .          (jbcinfo(nbl,nseg,5,1)-jbcinfo(nbl,nseg,4,1))
            end do
            do nseg=1,nbcjdim(nbl)
              if(jbcinfo(nbl,nseg,7,2).ne.0) jbc = jbc +
     .          (jbcinfo(nbl,nseg,3,2)-jbcinfo(nbl,nseg,2,2))*
     .          (jbcinfo(nbl,nseg,5,2)-jbcinfo(nbl,nseg,4,2))
            end do
            if (jbc .gt. jbcmax) jbcmax = jbc
c
            kbc = 0
            do nseg=1,nbck0(nbl)
              if(kbcinfo(nbl,nseg,7,1).ne.0) kbc = kbc +
     .          (kbcinfo(nbl,nseg,3,1)-kbcinfo(nbl,nseg,2,1))*
     .          (kbcinfo(nbl,nseg,5,1)-kbcinfo(nbl,nseg,4,1))
            end do
            do nseg=1,nbckdim(nbl)
              if(kbcinfo(nbl,nseg,7,2).ne.0) kbc = kbc +
     .          (kbcinfo(nbl,nseg,3,2)-kbcinfo(nbl,nseg,2,2))*
     .          (kbcinfo(nbl,nseg,5,2)-kbcinfo(nbl,nseg,4,2))
            end do
            if (kbc .gt. kbcmax) kbcmax = kbc
c
            ibc = 0
            do nseg=1,nbci0(nbl)
              if(ibcinfo(nbl,nseg,7,1).ne.0) ibc = ibc +
     .          (ibcinfo(nbl,nseg,3,1)-ibcinfo(nbl,nseg,2,1))*
     .          (ibcinfo(nbl,nseg,5,1)-ibcinfo(nbl,nseg,4,1))
            end do
            do nseg=1,nbcidim(nbl)
              if(ibcinfo(nbl,nseg,7,2).ne.0) ibc = ibc +
     .          (ibcinfo(nbl,nseg,3,2)-ibcinfo(nbl,nseg,2,2))*
     .          (ibcinfo(nbl,nseg,5,2)-ibcinfo(nbl,nseg,4,2))
            end do
           if (ibc .gt. ibcmax) ibcmax = ibc
c
         end do
c
      end if
c
      nstart = 1
c
      do 1000 nbl=1,nblock
c
      ivisc(1) = iviscg(nbl,1)
      ivisc(2) = iviscg(nbl,2)
      ivisc(3) = iviscg(nbl,3)
      igrid    = igridg(nbl)
      level    = levelg(nbl)
      iuns     = max(itrans(nbl),irotat(nbl),idefrm(nbl))
      iem      = iemg(igrid)
c
c     determine if there are any finer embedded blocks for
c     global block nbl
c
      nfiner = 0
      do mbl = 1,nblock
         if (levelg(mbl) .gt. lglobal) then
            if (nblcg(mbl).eq.nbl) then
               nfiner = nfiner + 1
               nblfine(nfiner) = mbl
            end if
         end if
      end do
c
      j  = jdimg(nbl)
      k  = kdimg(nbl)
      i  = idimg(nbl)
      j1 = j-1
      k1 = k-1
      i1 = i-1
c
      if (myid.ne.myhost .or. (myid.eq.myhost .and. mpihost.eq.0)) then
         jmkmim    = j*k*i
         jm1km1im1 = j1*k1*i1
         jmkmim1   = j*k*i1
         kmim1     = k*i1
         jmim1     = j*i1
         jmkm      = j*k
         kmim      = k*i
         jmim      = j*i
      end if
c
c
c     augment pointers if host or if block is local to this node,
c     or if this block needs embedded grid data from a non-local block
c
      nblc = nblcg(nbl)
      iaug = 0
      if (myid.eq.myhost .or. myid.eq.mblk2nd(nbl)) iaug = 1
      if (myid.eq.mblk2nd(nblc)) iaug = 1
      if (nfiner.gt.0) then
         do nf = 1,nfiner
            if (myid.eq.mblk2nd(nblfine(nf))) iaug = 1
         end do
      end if
c
      if (iaug.gt.0) then
c 
c........q(jdim,kdim,idim,5)
c        flow variables at cell centers
c
         lw(1,nbl) = nstart
         ns = jmkmim*5
c
c........qj0(kdim,idim-1,5,4)
c        flow variables at j-boundary
c
         lw(2,nbl) = lw(1,nbl) + ns
         ns = kmim1*5*4
c
c........qk0(jdim,idim-1,5,4)
c        flow variables at k-boundary
c
         lw(3,nbl)  = lw(2,nbl) + ns
         ns = jmim1*5*4
c
c........qi0(jdim,kdim,5,4)
c        flow variables at i-boundary
c
         lw(4,nbl) = lw(3,nbl) + ns
         ns = jmkm*5*4
c
c........sj(jdim,kdim,idim-1,5)
c        j-face metrics
c
         lw(5,nbl)  = lw(4,nbl) + ns
         ns = jmkmim1*5
c
c........sk(jdim,kdim,idim-1,5)
c        k-face metrics
c
         lw(6,nbl) = lw(5,nbl) + ns
         ns = jmkmim1*5
c
c........si(jdim,kdim,idim,5)
c        i-face metrics
c
         lw(7,nbl) = lw(6,nbl) + ns
         ns = jmkmim*5
c
c........vol(jdim,kdim,idim-1)
c        cell volume
c
         lw(8,nbl) = lw(7,nbl) + ns
         ns = jmkmim1
c
c........dtj(jdim,kdim,idim-1)    
c        cell time step
c
         lw(9,nbl) = lw(8,nbl) + ns
         ns = jmkmim1
c
c........x(jdim,kdim,idim)
c        grid point x coordinate
c
         lw(10,nbl) = lw(9,nbl) + ns
         ns = jmkmim
c
c........y(jdim,kdim,idim)
c        grid point y coordinate
c
         lw(11,nbl) = lw(10,nbl) + ns
         ns = jmkmim
c
c........z(jdim,kdim,idim)
c        grid point x coordinate
c
         lw(12,nbl) = lw(11,nbl) + ns
         ns = jmkmim
c
c........vist3d(jdim,kdim,idim)
c        eddy viscosity at cell center
c
         lw(13,nbl) = lw(12,nbl) + ns
         ns = 0
         if (ivmx.gt.1) ns = jmkmim
c
c........snk0(jdim-1,kdim-1,idim-1)
c        distance to nearest solid surface (field-eq turb models);
c        directed distance to nearest k-surface in current block
c        (B-L turb model)
c
         lw(14,nbl) = lw(13,nbl) +  ns
         ns = 0
         if (ivisc(3).gt.1 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4)
     .   ns = jm1km1im1
c
c........sni0/snj0/xjb(jdim-1,kdim-1,idim-1)
c        directed distance to nearest i-surface or j-surface in 
c        current block (B-L turb model); j-location of the nearest
c        surface point (B-B turb model or LES#25)
c
c        sni0 and snj0 share same storage location because
c        both can not be used at the same time for any given block
c        (a limitation of the currernt implementaion of the B-L model)
c
         lw(15,nbl) = lw(14,nbl) + ns
         ns = 0
         if (ivisc(1).eq.2 .or. ivisc(2).eq.2 .or. ivisc(3).eq.4 .or.
     .       ivisc(2).eq.4 .or. ivisc(1).eq.4 .or. ivisc(3).eq.25.or.
     .       ivisc(2).eq.25.or. ivisc(1).eq.25)
     .   ns = jm1km1im1
c
c........q1(jdim,kdim,idim,5)
c        cell center flow variables restricted from finer grid
c
c        global grids with multigrid:
c          q1 level i stored at level i except 
c          for finest level of global grids
c
c        embedded grids with multigrid :
c          q1 for level i stored at level i except
c          for finest level of embedded grids
c
         lw(16,nbl) = lw(15,nbl) + ns
         ns = jmkmim*5
         if (mgflag.eq.0) ns = 0
         if (iem.eq.0) then
c           global grid
            if (mgflag.eq.1 .and. level.eq.lglobal) ns = 0
         else
c           embedded grid
            if (mgflag.ne.2 .or. level.eq.maxl) ns = 0
         end if
c
c........qr(jdim,kdim,idim-1,5)
c        cell center residuals restricted from finer grid
c
c        global grids with multigrid:
c          qr for level i stored at level i except
c          for finest level of global grids
c
c        embedded grids with/without multigrid:
c          qr for level i stored at level i+1
c
         lw(17,nbl) = lw(16,nbl) + ns
         ns = jmkmim1*5
c
         if (iem.eq.0) then
c           global grid
            if (mgflag.eq.0)      ns = 0
            if (level.eq.lglobal) ns = 0
         end if
c
c........blank(jdim,kdim,idim)
c        cell center blank values (primarily overset applications)
c
         lw(18,nbl) = lw(17,nbl) + ns
         ns = jmkmim
c
c........xib/tursav(jdim,kdim,idim,2)
c        Baldwin-Barth, Spalart, k-omega or SST turbulence model 
c        parameters (for B-B & LES#25, the second half of storage is used
c        for k-index of nearest surface point; for Spalart, the
c        second half of storage is allocated but unused)
c
         lw(19,nbl) = lw(18,nbl) + ns
         ns = 0
         if(ivisc(3) .ge. 4 .or. ivisc(2) .ge. 4 .or. ivisc(1) .ge. 4)
     .      ns = jmkmim*2
c        extra memory needed for ivisc=30 (3-eqn model)
         if(ivisc(3).eq.30 .or. ivisc(2).eq.30 .or. ivisc(1).eq.30)
     .      ns = jmkmim*3
c        extra memory needed for ivisc=40 (4-eqn model)
         if(ivisc(3).eq.40 .or. ivisc(2).eq.40 .or. ivisc(1).eq.40)
     .      ns = jmkmim*4
c	 extra memory needed for ivisc=72 (stress-omega model,Wilcox 3rd Edition)
#   ifdef CMPLX
c        RSM not working in complex mode
#   else
         if(ivisc(3).eq.72 .or. ivisc(2).eq.72 .or. ivisc(1).eq.72)
     .      ns = jmkmim*kws_get_nummem()
#   endif
c
c........Johnson-King data no longer stored
c
         lw(20,nbl) = lw(19,nbl) + ns
         ns = 0
         lw(21,nbl) = lw(20,nbl) + ns
         ns = 0
         lw(22,nbl) = lw(21,nbl) + ns
         ns = 0
c
c........tj0(kdim,idim-1,2,4)
c        turbulence model variable at j-boundary
c
         lw(23,nbl) = lw(22,nbl) + ns
         ns = 0
         if(ivisc(3) .ge. 4 .or. ivisc(2) .ge. 4 .or. ivisc(1) .ge. 4)
     .      ns = kmim1*2*4
c        extra memory needed for ivisc=30 (3-eqn model)
         if(ivisc(3).eq.30 .or. ivisc(2).eq.30 .or. ivisc(1).eq.30)
     .      ns = kmim1*3*4
c        extra memory needed for ivisc=40 (4-eqn model)
         if(ivisc(3).eq.40 .or. ivisc(2).eq.40 .or. ivisc(1).eq.40)
     .      ns = kmim1*4*4
c        extra memory needed for ivisc=72 (stress-omega  model)
#   ifdef CMPLX
c        RSM not working in complex mode
#   else
         if(ivisc(3).eq.72 .or. ivisc(2).eq.72 .or. ivisc(1).eq.72)
     .      ns = kmim1*kws_get_nummem()*4
#   endif

c
c........tk0(jdim,idim-1,2,4)
c        turbulence model variable at k-boundary
c
         lw(24,nbl) = lw(23,nbl) + ns
         ns = 0
         if(ivisc(3) .ge. 4 .or. ivisc(2) .ge. 4 .or. ivisc(1) .ge. 4)
     .      ns = jmim1*2*4
c        extra memory needed for ivisc=30 (3-eqn model)
         if(ivisc(3).eq.30 .or. ivisc(2).eq.30 .or. ivisc(1).eq.30)
     .      ns = jmim1*3*4
c        extra memory needed for ivisc=40 (4-eqn model)
         if(ivisc(3).eq.40 .or. ivisc(2).eq.40 .or. ivisc(1).eq.40)
     .      ns = jmim1*4*4
c        extra memory needed for ivisc=72 (4-eqn model)
#   ifdef CMPLX
c        RSM not working in complex mode
#   else
         if(ivisc(3).eq.72 .or. ivisc(2).eq.72 .or. ivisc(1).eq.72)
     .      ns = jmim1*kws_get_nummem()*4
#   endif
c
c........ti0(jdim,jdim,2,4)
c        turbulence model variable at i-boundary
c
         lw(25,nbl) = lw(24,nbl) + ns
         ns = 0
         if(ivisc(3) .ge. 4 .or. ivisc(2) .ge. 4 .or. ivisc(1) .ge. 4)
     .      ns = jmkm*2*4
c        extra memory needed for ivisc=30 (3-eqn model)
         if(ivisc(3).eq.30 .or. ivisc(2).eq.30 .or. ivisc(1).eq.30)
     .      ns = jmkm*3*4
c        extra memory needed for ivisc=40 (4-eqn model)
         if(ivisc(3).eq.40 .or. ivisc(2).eq.40 .or. ivisc(1).eq.40)
     .      ns = jmkm*4*4
c        extra memory needed for ivisc=72 (stress-omega model)
#   ifdef CMPLX
c        RSM not working in complex mode
#   else
         if(ivisc(3).eq.72 .or. ivisc(2).eq.72 .or. ivisc(1).eq.72)
     .      ns = jmkm*kws_get_nummem()*4
#   endif
c
c........xkb(jdim-1,kdim-1,idim-1)
c        k-index of nearest surface point (B-B turb model or LES#25)
c        directed distance from MAX on k (B-L turb model)
c
         lw(26,nbl) = lw(25,nbl) + ns
         ns = 0
         if(ivisc(3).eq.2 .or.
     .      ivisc(3).eq.4 .or. ivisc(2).eq.4 .or. ivisc(1).eq.4 .or.
     .      ivisc(3).eq.25.or. ivisc(2).eq.25.or. ivisc(1).eq.25)
     .      ns = jm1km1im1
c
c........blnum(jdim-1,kdim-1,idim-1)
c        block number containing nearest surface point (B-B turb model
c        or LES#25)
c        directed distance from MAX on i or j (B-L turb model)
c
         lw(27,nbl) = lw(26,nbl) + ns
         ns = 0
         if(ivisc(2).eq.2 .or. ivisc(1).eq.2 .or.
     .      ivisc(3).eq.4 .or. ivisc(2).eq.4 .or. ivisc(1).eq.4 .or.
     .      ivisc(3).eq.25.or. ivisc(2).eq.25.or. ivisc(1).eq.25)
     .      ns = jm1km1im1
c
c........vj0(kdim,idim-1,1,4)
c        eddy viscosity at j-boundary
c
         lw(28,nbl) = lw(27,nbl) + ns
         ns = 0
         if(ivmx.gt.1) ns = kmim1*1*4
c
c........vk0(jdim,idim-1,1,4)
c        eddy viscosity at k-boundary
c
         lw(29,nbl) = lw(28,nbl) + ns
         ns = 0
         if(ivmx.gt.1) ns = jmim1*1*4
c
c........vi0(jdim,kdim,1,4)
c        eddy viscosity at i-boundary
c
         lw(30,nbl) = lw(29,nbl) + ns
         ns = 0
         if(ivmx.gt.1) ns =jmkm*1*4
c
c........bcj(kdim,idim-1,2)
c        boundary data type flag at j-boundary
c
         lw(31,nbl) = lw(30,nbl) + ns
         ns         = kmim1*2
c
c........bck(jdim,idim-1,2)
c        boundary data type flag at k-boundary
c
         lw(32,nbl) = lw(31,nbl) + ns
         ns         = jmim1*2
c
c........bci(jdim,kdim,2)
c        boundary data type flag at i-boundary
c
         lw(33,nbl) = lw(32,nbl) + ns
         ns         = jmkm*2
c
c........qc0(jdim,kdim,idim-1,5)
c        flow variables at previous time step/subiteration
c
         lw(34,nbl) = lw(33,nbl) + ns
         ns = 0
         if (real(dt).gt.0.) ns = jmkmim1*5
c
c........dqc0(jdim,kdim,idim-1,5)
c        change in flow variables from previous time step/subiteration
c
         lw(35,nbl) = lw(34,nbl) + ns
         ns = 0
         if (real(dt).gt.0.) ns=jmkmim1*5
c
c........xtbj(kdim,idim-1,3,2)
c        j-boundary velcity compoments on moving grids
c
         lw(36,nbl) = lw(35,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=kmim1*3*2
c
c........xtbk(jdim,idim-1,3,2)
c        k-boundary velcity compoments on moving grids
c
         lw(37,nbl) = lw(36,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=jmim1*3*2
c
c........xtbi(jdim,kdim,3,2)
c        i-boundary velcity compoments on moving grids
c
         lw(38,nbl) = lw(37,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=jmkm*3*2
c
c........atbj(kdim,idim-1,3,2)
c        j-boundary acceleration compoments on moving grids
c
         lw(39,nbl) = lw(38,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=kmim1*3*2
c
c........atbk(jdim,idim-1,3,2)
c        k-boundary acceleration compoments on moving grids
c
         lw(40,nbl) = lw(39,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=jmim1*3*2
c
c........atbi(jdim,kdim,3,2)
c        i-boundary acceleration compoments on moving grids
c
         lw(41,nbl) = lw(40,nbl) + ns
         ns = 0
         if (iuns.gt.0) ns=jmkm*3*2
c
c........bcdataj(npts,mpts,2,12)
c        2000 series bc data for j-boundaries
c        lwdat further subdivides the pointers starting at lw(42,nbl)
c        to indicate a starting location for each boundary segment
c
         lw(42,nbl) = lw(41,nbl) + ns
c
         ns  = 0
         ns1 = 0
         do 2600 nseg=1,nbcj0(nbl)
         lwdat(nbl,nseg,3) = lw(42,nbl) + ns1
         if(jbcinfo(nbl,nseg,7,1).ne.0) ns1 = ns1 +
     .       (jbcinfo(nbl,nseg,3,1)-jbcinfo(nbl,nseg,2,1))*
     .       (jbcinfo(nbl,nseg,5,1)-jbcinfo(nbl,nseg,4,1))*12*2
 2600    continue
         do 2610 nseg=1,nbcjdim(nbl)
         lwdat(nbl,nseg,4) = lw(42,nbl) + ns1
         if(jbcinfo(nbl,nseg,7,2).ne.0) ns1 = ns1 +
     .       (jbcinfo(nbl,nseg,3,2)-jbcinfo(nbl,nseg,2,2))*
     .       (jbcinfo(nbl,nseg,5,2)-jbcinfo(nbl,nseg,4,2))*12*2
 2610    continue
c
         ns = ns + ns1
         if (myid.eq.myhost .and. mpihost.eq.1) then
            ns = jbcmax*12*2
         end if
c
c........bcdatak(npts,mpts,2,12)
c        2000 series bc data for k-boundaries
c        lwdat further subdivides the pointers starting at lw(43,nbl)
c        to indicate a starting location for each boundary segment
c
         lw(43,nbl) = lw(42,nbl) + ns
c
         ns  = 0
         ns1 = 0
         do 2620 nseg=1,nbck0(nbl)
         lwdat(nbl,nseg,5) = lw(43,nbl) + ns1
         if(kbcinfo(nbl,nseg,7,1).ne.0) ns1 = ns1 +
     .       (kbcinfo(nbl,nseg,3,1)-kbcinfo(nbl,nseg,2,1))*
     .       (kbcinfo(nbl,nseg,5,1)-kbcinfo(nbl,nseg,4,1))*12*2
 2620    continue
         do 2630 nseg=1,nbckdim(nbl)
         lwdat(nbl,nseg,6) = lw(43,nbl) + ns1
         if(kbcinfo(nbl,nseg,7,2).ne.0) ns1 = ns1 +
     .       (kbcinfo(nbl,nseg,3,2)-kbcinfo(nbl,nseg,2,2))*
     .       (kbcinfo(nbl,nseg,5,2)-kbcinfo(nbl,nseg,4,2))*12*2
 2630    continue
c
         ns = ns + ns1
         if (myid.eq.myhost .and. mpihost.eq.1) then
            ns = kbcmax*12*2
         end if
c
c........bcdatai(npts,mpts,2,12)
c        2000 series bc data for i-boundaries
c        lwdat further subdivides the pointers starting at lw(44,nbl)
c        to indicate a starting location for each boundary segment
c
         lw(44,nbl) = lw(43,nbl) + ns
c
         ns  = 0
         ns1 = 0
         do 2640 nseg=1,nbci0(nbl)
         lwdat(nbl,nseg,1) = lw(44,nbl)+ns1
         if(ibcinfo(nbl,nseg,7,1).ne.0) ns1 = ns1 +
     .       (ibcinfo(nbl,nseg,3,1)-ibcinfo(nbl,nseg,2,1))*
     .       (ibcinfo(nbl,nseg,5,1)-ibcinfo(nbl,nseg,4,1))*12*2
 2640    continue
         do 2650 nseg=1,nbcidim(nbl)
         lwdat(nbl,nseg,2) = lw(44,nbl)+ns1
         if(ibcinfo(nbl,nseg,7,2).ne.0) ns1 = ns1 +
     .       (ibcinfo(nbl,nseg,3,2)-ibcinfo(nbl,nseg,2,2))*
     .       (ibcinfo(nbl,nseg,5,2)-ibcinfo(nbl,nseg,4,2))*12*2
 2650    continue
c
         ns = ns + ns1
         if (myid.eq.myhost .and. mpihost.eq.1) then
            ns = ibcmax*12*2
         end if
c
c........xib2(jdim,kdim,idim,4)
c        extra storage for subiterations for time-accurate turb. models
c
         lw(45,nbl) = lw(44,nbl) + ns
         ns = 0
         if((ivisc(3) .ge. 4 .or. ivisc(2) .ge. 4 .or. ivisc(1) .ge. 4)
     .      .and. real(dt) .gt. 0.) ns=jmkmim*4
c        extra memory needed for ivisc=30 (3-eqn model)
         if(ivisc(3).eq.30 .or. ivisc(2).eq.30 .or. ivisc(1).eq.30)
     .      ns = jmkmim*6
c        extra memory needed for ivisc=40 (4-eqn model)
         if(ivisc(3).eq.40 .or. ivisc(2).eq.40 .or. ivisc(1).eq.40)
     .      ns = jmkmim*8
c        extra memory needed for ivisc=72 (stress-omega model)
         if(ivisc(3).eq.72 .or. ivisc(2).eq.72 .or. ivisc(1).eq.72)
     .      ns = jmkmim*14
c
c........qavg(jdim,kdim,idim,5) iteration-averaged data
c
         lw(46,nbl) = lw(45,nbl) + ns
         ns = 0

         if (iteravg.eq.1 .or. iteravg.eq.2) ns = jmkmim*5
c
c........ux(jdim-1,kdim-1,idim-1,9)

c        velocity gradients for EASM k-omega or k-epsilon turb models
c        or for k-o/SST/k-e model when isstprod or isstdenom are active
c        or for any ivisc >= 11 or for SARC or SSTRC or for Wilcox06
c
         lw(47,nbl) = lw(46,nbl) + ns
         ns = 0
         if (lowmem_ux .eq. 1) then
           if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.ge.11 .or.
     .         ivmx.eq.16 .or.
     .        ((ivmx.eq.6 .or. ivmx.eq.7 .or. ivmx.eq.10) .and.
     .        ikoprod.eq.1) .or.
     .        (ivmx.eq.7 .and. isstdenom.eq.1) .or.
     .        (ivmx.eq.7 .and. isst2003.eq.1) .or.
     .        ((ivmx.eq.6 .or. ivmx.eq.7) .and. isstrc.gt.0) .or.
     .        ((ivmx.eq.6 .or. ivmx.eq.7) .and. isstsf.eq.1) .or.
     .        ((ivmx.eq.5) .and. isarc2d.eq.1) .or.
     .        ((ivmx.eq.5) .and. isarc3d.eq.1) .or.
     .        ((ivmx.eq.5) .and. isar.eq.1) .or.
     .        (ivmx.eq.6 .and. i_wilcox06.eq.1) .or.
     .        (ivmx.eq.6 .and. i_wilcox98.eq.1) .or. ivmax==72 .or.
     .        ides.ge.2 .or. i_nonlin.ne.0 .or. i_tauijs.ne.0)
     .        ns = jm1km1im1*9
         else
           ns = jm1km1im1*9
         end if
c
c        cmuv(jdim-1,kdim-1,idim-1)
c........EASM, var-g (now considered the "Standard" EASM)
c
         lw(48,nbl) = lw(47,nbl) + ns
         ns = 0
         if(ivisc(3).eq.8 .or. ivisc(2).eq.8 .or. ivisc(1).eq.8 .or.
     .      ivisc(3).eq.9 .or. ivisc(2).eq.9 .or. ivisc(1).eq.9 .or.
     .      ivisc(3).eq.13 .or. ivisc(2).eq.13 .or. ivisc(1).eq.13 .or.
     .      ivisc(3).eq.14 .or. ivisc(2).eq.14 .or. ivisc(1).eq.14)
     .      ns = jm1km1im1
c
c        volj0(kdim,idim-1,1,4)
c........cell volumes at j-boundary
c
         lw(49,nbl) = lw(48,nbl) + ns
         ns = 0
         if(ivmx.ge.1) ns =kmim1*1*4
         
c
c        volk0(jdim,idim-1,1,4)
c........cell volumes at k-boundary
c
         lw(50,nbl) = lw(49,nbl) + ns
         ns = 0
         if(ivmx.ge.1) ns =jmim1*1*4

c
c        voli0(jdim,kdim,1,4)
c........cell volumes at i-boundary
c
         lw(51,nbl) = lw(50,nbl) + ns
         ns = 0
         if(ivmx.ge.1) ns =jmkm*1*4
c
c........aeroelatic j-surface data
c
         lw(52,nbl) = lw(51,nbl) + ns
         ns = 0
         if(naesrf.gt.0) ns =kmim*6*maxaes*nmds
c
c........aeroelatic k-surface data
c
         lw(53,nbl) = lw(52,nbl) + ns
         ns = 0
         if(naesrf.gt.0) ns =jmim*6*maxaes*nmds
c
c........aeroelatic i-surface data
c
         lw(54,nbl) = lw(53,nbl) + ns
         ns = 0
         if(naesrf.gt.0) ns =jmkm*6*maxaes*nmds
c
c........intermediate grid velocity for blocks undergoing BOTH
c        rigid motion and deforming motion
c        Note: this data no longer stored
c
         lw(55,nbl) = lw(54,nbl) + ns
         ns = 0
c
c........x-coordinate at time step n-2 (only for deforming
c        grids, when second order time derivatives are needed)
c
         lw(56,nbl) = lw(55,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0 .and. abs(ita).gt.1) ns = jmkmim
c
c........y-coordinate at time step n-2 (only for deforming
c        grids, when second order time derivatives are needed)
c
         lw(57,nbl) = lw(56,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0 .and. abs(ita).gt.1) ns = jmkmim
c
c........z-coordinate at time step n-2 (only for deforming
c        grids, when second order time derivatives are needed)
c
         lw(58,nbl) = lw(57,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0 .and. abs(ita).gt.1) ns = jmkmim
c
c........displacements of j=const surfaces for deforming grids
c
         lw(59,nbl) = lw(58,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = kmim*3*2
c
c........displacements of k=const surfaces for deforming grids
c
         lw(60,nbl) = lw(59,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = jmim*3*2
c
c........displacements of i=const surfaces for deforming grids
c
         lw(61,nbl) = lw(60,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = jmkm*3*2
c
c........x-coordinate at time step n-1
c
         lw(62,nbl) = lw(61,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = jmkmim
c
c........y-coordinate at time step n-1
c
         lw(63,nbl) = lw(62,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = jmkmim
c
c........z-coordinate at time step n-1
c
         lw(64,nbl) = lw(63,nbl) + ns
         ns = 0
         if (idefrm(nbl).gt.0) ns = jmkmim
c
c
c........q2avg(jdim,kdim,idim,5) Average of the square of q
c
         lw(65,nbl) = lw(64,nbl) + ns
         ns = 0
         if (ipertavg.eq.1 .or. ipertavg.eq.2) ns = jmkmim*5
c
         nstart = lw(65,nbl) + ns
c
         if (myid.eq.myhost .and. mpihost.eq.1) then
            if (nbl.lt.nblock) then
               nstart = 1
            end if
         end if
c
c        memblock is the memory need to store the w array for this 
c        block; it is used in precfl3d to determine the memory 
c        requirement for each node
c
         memblock(nbl) = lw(65,nbl) + ns - lw(1,nbl) + 1
c
      else
c
c        if block is not local to this node, set pointers to previous
c        blocks starting location
c
         do ii = 1,65
            lw(ii,nbl) = nstart
         end do
         do ii = 1, maxseg
            do jj = 1,6
               lwdat(nbl,ii,jj) = nstart
            end do
         end do
c
         nstart = lw(65,nbl)
c
      end if
c
c     set up lw2 array
c
      jdim        = jdimg(nbl)
      kdim        = kdimg(nbl)
      idim        = idimg(nbl)
      lw2(1, nbl) = jdimg(nbl)
      lw2(2, nbl) = kdimg(nbl)
      lw2(3, nbl) = idimg(nbl)
      lw2(4, nbl) = nblcg(nbl)
      lw2(5, nbl) = jdim/2+1
      lw2(6, nbl) = kdim/2+1
      lw2(7, nbl) = idim/2+1
      lw2(8, nbl) = ilamlog(nbl)
      lw2(9, nbl) = ilamhig(nbl)
      lw2(10,nbl) = jlamlog(nbl)
      lw2(11,nbl) = jlamhig(nbl)
      lw2(12,nbl) = klamlog(nbl)
      lw2(13,nbl) = klamhig(nbl)
      lw2(14,nbl) = iviscg(nbl,1)
      lw2(15,nbl) = iviscg(nbl,2)
      lw2(16,nbl) = iviscg(nbl,3)
      lw2(17,nbl) = idegg(nbl,1)
      lw2(18,nbl) = idegg(nbl,2)
      lw2(19,nbl) = idegg(nbl,3)
      lw2(20,nbl) = iwfg(nbl,1)
      lw2(21,nbl) = iwfg(nbl,2)
      lw2(22,nbl) = iwfg(nbl,3)
      lw2(23,nbl) = idiagg(nbl,1)
      lw2(24,nbl) = idiagg(nbl,2)
      lw2(25,nbl) = idiagg(nbl,3)
      lw2(26,nbl) = iflimg(nbl,1)
      lw2(27,nbl) = iflimg(nbl,2)
      lw2(28,nbl) = iflimg(nbl,3)
      lw2(29,nbl) = ifdsg(nbl,1)
      lw2(30,nbl) = ifdsg(nbl,2)
      lw2(31,nbl) = ifdsg(nbl,3)
c     factor multiplying rkap0g must be consistant with subroutine lead!
      lw2(32,nbl) = int(rkap0g(nbl,1)*1e6)
      lw2(33,nbl) = int(rkap0g(nbl,2)*1e6)
      lw2(34,nbl) = int(rkap0g(nbl,3)*1e6)
      lw2(35,nbl) = jsg(nbl)
      lw2(36,nbl) = ksg(nbl)
      lw2(37,nbl) = isg(nbl)
      lw2(38,nbl) = jeg(nbl)
      lw2(39,nbl) = keg(nbl)
      lw2(40,nbl) = ieg(nbl)
      if (iemg(igridg(nbl)).gt.0) then
         nblc       = nblcg(nbl)
         lw2(5,nbl) = jdimg(nblc)
         lw2(6,nbl) = kdimg(nblc)
         lw2(7,nbl) = idimg(nblc)
      end if
c
 1000 continue
c
      ftot_w  = float(nstart)/1.e+06
      nwork   = mwork-nstart
      ftot_wk = float(nwork)/1.e+06
c
      if (icall .gt. 0) then
c
#if defined DIST_MPI
         if (myid.eq.myhost) then
            write(11,*)
            write(11,1997)
            write(11,1999)
            write(11,3500) myid,real(ftot_w),real(ftot_wk),
     .                     real(ftot_w+ftot_wk)
         end if
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),3500) myid,real(ftot_w),real(ftot_wk),
     .                             real(ftot_w+ftot_wk)
         if (real(ftot_wk).le.0.) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),4000)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
#else
         write(11,*)
         write(11,1997)
         write(11,2000) real(ftot_w)
         write(11,3000) real(ftot_wk)
         if (real(ftot_wk).le.0.) then
            write(11,4000)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
#endif
      end if
c
    4 format(5i6)
   37 format(2x,4h lw(,i2,1h,,i3,3h)= ,i10)
   38 format(2x,7h lwdat(,i3,1h,,i2,1h,,i2,3h)= ,i10)
  223 format(/1x,5higrid,1x,5hblock,2x,4hjdim,2x,4hkdim,2x,4hidim)
  910 format(/52h SUMMARY OF STARTING LOCATIONS FOR BLOCK INFORMATION,
     .8h STORAGE)
  911 format(/52h SUMMARY OF STARTING LOCATIONS FOR BLOCK INFORMATION,
     .8h ON HOST)
 1997 format(1x,38hSUMMARY OF PRIMARY MEMORY REQUIREMENTS)
 1999 format(2x,39hnode   w storage (mw)   wk storage (mw),
     .          17h       total (mw))
 2000 format(2x,28hmemory (mw) for w  storage =,e12.5)
 3000 format(2x,28hmemory (mw) for wk storage =,e12.5)
 3500 format(2x,i4,5x,e12.5,6x,e12.5,5x,e12.5)
 4000 format(2x,26hstopping...mwork too small)
c
      return
      end
