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
      subroutine tlns3d_to_cfl3d(iform,iver,ipatch,ironver,ngrid,
     .           iln,jln,kln,ipar,nnodes,isd,ifrom,xif1,xif2,etf1,
     .           etf2,ito,xit1,xit2,ett1,ett2,nfb,iredundant,iconcat,
     .           xmap,msegn,mbloc,iovrlp,nseg,imap,ivisb,twotref,
     .           itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,iturbb,
     .           cq,cporous,cpchamber,cradiation,
     .           llimit,iitmax,mmcxie,mmceta,iifit,iic0,iiorph,
     .           bcvali,bcvalj,bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,
     .           nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,nblk,limblk,
     .           isva,nblon,rkap0g,levelg,igridg,iflimg,ifdsg,iviscg,
     .           jdimg,kdimg,idimg,idiagg,nblcg,idegg,jsg,ksg,isg,
     .           jeg,keg,ieg,mit,ilamlog,ilamhig,jlamlog,jlamhig,
     .           klamlog,klamhig,ncgg,nblg,iemg,inewgg,iadvance,iforce,
     .           intmax,nsub1,msegt,maxseg,mxbli,nrotat,ntrans,tlref,
     .           rlref,ioflag,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c
c     Converts TLNS3D input and map files into a CFL3D input file.
c
c     NOTES:
c
c       "TRUST BUT VERIFY" -  R. Reagan
c
c        due to differences between tlns3d and cfl3d, it is not always
c        possible to translate between the two codes without some
c        "interpretation", so it is strongly suggested that the user
c        check over the cfl3d input file (and ronnie input files for
c        patched grids).
c
c***********************************************************************
c
c     mbloc/maxgr...maximum number of grids
c     maxseg........maximum number of segments/grid face
c     mxbli.........maximum number if 1-1 zonal interfaces
c     intmax........maximum number of patched interfaces
c     nsub1.........maximum number of source ("from") zones
c                   that may patch to a target ("to") zone
c     msegt.........maximum number of entries in tlns3d map
c                   file (currently only 20 entries are actually used) 
c
      character*10 datahdr(10)
      character*80 gridin,plt3dg,plt3dq,output,resid,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             gridout,sdgridin,sdgridout,dovrlap,dpatch,dresid,
     .             rout
      character*120 bou(ibufdim,nbuf)
c
      integer ifrom(intmax,nsub1),xif1(intmax,nsub1),
     .        xif2(intmax,nsub1),etf1(intmax,nsub1),
     .        etf2(intmax,nsub1),ito(intmax),xit1(intmax),
     .        xit2(intmax),ett1(intmax),ett2(intmax),nfb(intmax)
c
      dimension nou(nbuf)
      dimension iredundant(mbloc,msegn),iconcat(intmax)
      dimension xmap(msegt,msegn,mbloc),iovrlp(mbloc)
      dimension nseg(mbloc)
      dimension imap(msegt,msegn,mbloc)
      dimension ivisb(msegt,mbloc),twotref(msegn,mbloc)
      dimension itrb1(mbloc),itrb2(mbloc),
     .          jtrb1(mbloc),jtrb2(mbloc),
     .          ktrb1(mbloc),ktrb2(mbloc),
     .          iturbb(mbloc)
      dimension cq(msegn,mbloc),cporous(msegn,mbloc),
     .          cpchamber(msegn,mbloc),cradiation(msegn,mbloc)
      dimension iln(mbloc),jln(mbloc),kln(mbloc)
      dimension llimit(intmax),iitmax(intmax),mmcxie(intmax),
     .          mmceta(intmax),iifit(intmax),iic0(intmax),
     .          iiorph(intmax)
      dimension bcvali(mbloc,maxseg,7,2),
     .          bcvalj(mbloc,maxseg,7,2),bcvalk(mbloc,maxseg,7,2),
     .          nbci0(mbloc),nbcidim(mbloc),nbcj0(mbloc),nbcjdim(mbloc),
     .          nbck0(mbloc),nbckdim(mbloc),ibcinfo(mbloc,maxseg,7,2),
     .          jbcinfo(mbloc,maxseg,7,2),kbcinfo(mbloc,maxseg,7,2)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli),nblon(mxbli)
      dimension rkap0g(mbloc,3),levelg(mbloc),igridg(mbloc),
     .          iflimg(mbloc,3),ifdsg(mbloc,3),iviscg(mbloc,3),
     .          jdimg(mbloc),kdimg(mbloc),idimg(mbloc),idiagg(mbloc,3),
     .          nblcg(mbloc),idegg(mbloc,3),
     .          jsg(mbloc),ksg(mbloc),isg(mbloc),jeg(mbloc),keg(mbloc),
     .          ieg(mbloc),mit(5,mbloc),
     .          ilamlog(mbloc),ilamhig(mbloc),jlamlog(mbloc),
     .          jlamhig(mbloc),klamlog(mbloc),klamhig(mbloc)
      dimension ncgg(mbloc),nblg(mbloc),iemg(mbloc),
     .          inewgg(mbloc)
      dimension iadvance(mbloc),iforce(mbloc)
c
      common /unit5/ iunit5
      common /cflfiles/gridin,plt3dg,plt3dq,output,resid,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 gridout,sdgridin,sdgridout,dovrlap,dpatch,dresid
      common /ronfiles/rout
      common /ron/ titleron(20)
      common /chk1/ ncgmax
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /info2/ tinf,reue,c2spe,bref,cref,sref,xmc,ymc,zmc,
     .               rfreq,alphau,cloc,cfltau
      common /info3/ ncg,isnd,ialph,irest,iunst,ntstep,
     .               ita,ihstry,nplot3d,nprint,nwrest,ichk,i2d,
     .               mglev(5),nem(5),mitl(5,5),mtt,ndv,isdform,
     .               ip3dgrd,iplt3dtyp
      common /mgv/ epsssc(3),epsssr(3),issc,issr
      common /alphait/ cltarg,resupdt,cltol,dalim,dcldal,alphalast,
     .                 ialphit,nttlast,icycupdt
      common /precon1/ cprec,uref,avn,iprecon
c
c     set idimg,jdimg,kdim to the new (split) values
c
      do n=1,ngrid
         idimg(n) = iln(n)
         jdimg(n) = jln(n)
         kdimg(n) = kln(n)
      end do
      nbloc = ngrid
c
      nbloc = abs(nbloc)
      if (nbloc .gt. mbloc) then
         write(6,*) 'stopping: parameter mbloc must be ',
     .              'at least ',nbloc
         call termn8(0,-3,ibufdim,nbuf,bou,nou)
      end if
c
c     always use map file with third line (iold=1)
c
      iold = 1
c
c***********************************************************************
c
c     read in tlns3d map file
c
c***********************************************************************
c
      rewind(8)
      read(8,*)
      read(8,*) nbloc
      do 600 ibloc = 1,nbloc
c
      read(8,*)
      read(8,*)
      read(8,*)  nseg(ibloc),(ivisb(n,ibloc),n=1,3),
     .           iturbb(ibloc),itrb1(ibloc),itrb2(ibloc),jtrb1(ibloc),
     .           jtrb2(ibloc),ktrb1(ibloc),ktrb2(ibloc)
      read(8,*)
      read(8,*) (ivisb(n,ibloc),n=4,15)
      read(8,*) (ivisb(n,ibloc),n=16,27)
      read(8,*) (ivisb(n,ibloc),n=28,34)
      read(8,*)
c
      read(8,*)
      read(8,*)
c
c     imap array
c
      do 600 iseg  = 1,nseg(ibloc)
c
      if (iold.eq.0) then
         read(8,*)
         read(8,*) ibloc1,isegg,(imap(n,iseg,ibloc),n= 1,12)
         read(8,*) (imap(n,iseg,ibloc),n=13,20),twotref(iseg,ibloc)
      else
         read(8,*)
         read(8,*) ibloc1,isegg,(imap(n,iseg,ibloc),n= 1,12)
         read(8,*) (imap(n,iseg,ibloc),n=13,20)
      end if
c
      if (iold.gt.0) then
        read(8,*) twotref(iseg,ibloc),cq(iseg,ibloc),
     .            cporous(iseg,ibloc),cpchamber(iseg,ibloc),
     .            cradiation(iseg,ibloc)
      end if
c
c     additional data for cfl3d
c
      read(8,*)
      read(8,*) (xmap(n,iseg,ibloc),n=1,5)
      read(8,*) (xmap(n,iseg,ibloc),n=6,10)
      read(8,*) (xmap(n,iseg,ibloc),n=11,15)
      read(8,*) (xmap(n,iseg,ibloc),n=16,20)
      read(8,*) (xmap(n,iseg,ibloc),n=21,25)
      read(8,*) (xmap(n,iseg,ibloc),n=26,26)
      read(8,*) (imap(n,iseg,ibloc),n=21,23)
      read(8,*)
c
  600 continue
c***********************************************************************
c
c     translate tlns3d mapfile data to cfl3d data
c
c***********************************************************************
c
      do 900 n=1,ngrid
c
      iviscg(n,1) = ivisb(1,n)
      iviscg(n,2) = ivisb(2,n)
      iviscg(n,3) = ivisb(3,n)
c
c     ivv = 0 if there are no viscous blocks
      ivv = 1
      if(iviscg(n,1).eq.0 .and. iviscg(n,2).eq.0 .and. iviscg(n,3).eq.0)
     .ivv = 0
c
  900 continue
c
      iflag1 = 0
      nbli   = 0
      ninter = 0
c
      do n=1,ngrid
         do iseg=1,nseg(n)
            iredundant(n,iseg) = 0
         end do
      end do
c
      do 1000 n=1,ngrid
c
c     some items carried in the ivisb array are independent
c     of whether the original block gets split or not (e.g.
c     iadvance, idiag), but some are not (i.e. embedded grid
c     info, illamlo,etc). those that are not independent of
c     splitting are now almost certainly in error. all transfers
c     from the ivisb array are placed here for completeness
c     in case they are ever split correctly in future versions.
c     data that is believed to be incorrect is subsequently
c     overwritten with generic data
c
      ncgg(n)     = ivisb( 4,n)
      iemg(n)     = ivisb( 5,n)
      iadvance(n) = ivisb( 6,n)
      iforce(n)   = ivisb( 7,n)
      ilamlog(n)  = ivisb( 8,n)
      ilamhig(n)  = ivisb( 9,n)
      jlamlog(n)  = ivisb(10,n)
      jlamhig(n)  = ivisb(11,n)
      klamlog(n)  = ivisb(12,n)
      klamhig(n)  = ivisb(13,n)
      inewgg(n)   = ivisb(14,n)
      igridg(n)   = ivisb(15,n)
      isg(n)      = ivisb(16,n)
      jsg(n)      = ivisb(17,n)
      ksg(n)      = ivisb(18,n)
      ieg(n)      = ivisb(19,n)
      jeg(n)      = ivisb(20,n)
      keg(n)      = ivisb(21,n)
      idiagg(n,1) = ivisb(22,n)
      idiagg(n,2) = ivisb(23,n)
      idiagg(n,3) = ivisb(24,n)
      iflimg(n,1) = ivisb(25,n)
      iflimg(n,2) = ivisb(26,n) 
      iflimg(n,3) = ivisb(27,n)
      ifdsg(n,1)  = ivisb(28,n)
      ifdsg(n,2)  = ivisb(29,n)
      ifdsg(n,3)  = ivisb(30,n)
      irkapi      = ivisb(31,n)
      irkapj      = ivisb(32,n)
      irkapk      = ivisb(33,n)
      if(irkapi .eq. -1) rkapi = -1.0
      if(irkapi .eq.  0) rkapi = 0.
      if(irkapi .eq.  1) rkapi = 1.0
      if(irkapi .eq.  3) rkapi = 0.3333
      if(irkapj .eq. -1) rkapj = -1.0
      if(irkapj .eq.  0) rkapj = 0.
      if(irkapj .eq.  1) rkapj = 1.0
      if(irkapj .eq.  3) rkapj = 0.3333
      if(irkapk .eq. -1) rkapk = -1.0
      if(irkapk .eq.  0) rkapk = 0.
      if(irkapk .eq.  1) rkapk = 1.0
      if(irkapk .eq.  3) rkapk = 0.3333
      rkap0g(n,1) = rkapi
      rkap0g(n,2) = rkapj
      rkap0g(n,3) = rkapk
      iovrlp(n)   = ivisb(34,n)
c
c     questionable data now  overwritten
c     note that iforce(n) is later reset using
c     bc segment information, which is correct after
c     splitting
c 
      iforce(n)  = 0
      iemg(n)    = 0
      isg(n)     = 0
      jsg(n)     = 0
      ksg(n)     = 0
      ieg(n)     = 0
      jeg(n)     = 0
      keg(n)     = 0
      igridg(n)  = 0
      ilamlog(n) = 0
      ilamhig(n) = 0
      jlamlog(n) = 0
      jlamhig(n) = 0
      klamlog(n) = 0
      klamhig(n) = 0
c
c     initialize number of bc segments and force flags
c
      nbci0(n)   = 0
      nbcidim(n) = 0 
      nbcj0(n)   = 0
      nbcjdim(n) = 0 
      nbck0(n)   = 0
      nbckdim(n) = 0 
      ifo1       = 0
      ifo2       = 0
      jfo1       = 0
      jfo2       = 0
      kfo1       = 0
      kfo2       = 0
c
      do 1010 iseg=1,nseg(n)
      ibctt   = imap(1,iseg,n)
      if (ibctt.eq.0 .or. ibctt.eq.1 .or. ibctt.eq.-1) then
         ibctyp = 0
      else
         ibctyp = ibctt
      end if
c
      ift   = imap(2,iseg,n)
      i1tmn = imap(3,iseg,n)
      i1tmx = imap(4,iseg,n)
      i2tmn = imap(5,iseg,n)
      i2tmx = imap(6,iseg,n)
c
c     i0 boundary data
c
      if (ift.eq.1) then
c 
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbci0(n).gt.0 .and. ibctyp.eq.0) then
            m = 1
            if ((ibctyp.eq.ibcinfo(n,nbci0(n),1,m)) .and.
     .          (i1tmn .eq.ibcinfo(n,nbci0(n),2,m)) .and.
     .          (i1tmx .eq.ibcinfo(n,nbci0(n),3,m)) .and.
     .          (i2tmn .eq.ibcinfo(n,nbci0(n),4,m)) .and.
     .          (i2tmx .eq.ibcinfo(n,nbci0(n),5,m))) then
                go to 6001
            end if
         end if
         nbci0(n) = nbci0(n) + 1
         ns = nbci0(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 1
         ibcinfo(n,ns,1,m) = ibctyp
         ibcinfo(n,ns,2,m) = i1tmn
         ibcinfo(n,ns,3,m) = i1tmx
         ibcinfo(n,ns,4,m) = i2tmn
         ibcinfo(n,ns,5,m) = i2tmx
         ibcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               ifo1 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               ifo1 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         ibcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0)  then
            do mm=1,abs(ndata)
               bcvali(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6001 continue
c     
c     idim boundary data
c
      if (ift.eq.2) then
c
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbcidim(n).gt.0 .and. ibctyp.eq.0) then
            m = 2
            if ((ibctyp.eq.ibcinfo(n,nbcidim(n),1,m)) .and.
     .          (i1tmn .eq.ibcinfo(n,nbcidim(n),2,m)) .and.
     .          (i1tmx .eq.ibcinfo(n,nbcidim(n),3,m)) .and.
     .          (i2tmn .eq.ibcinfo(n,nbcidim(n),4,m)) .and.
     .          (i2tmx .eq.ibcinfo(n,nbcidim(n),5,m))) then
                go to 6002
            end if
         end if
         nbcidim(n) = nbcidim(n) + 1
         ns = nbcidim(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 2
         ibcinfo(n,ns,1,m) = ibctyp
         ibcinfo(n,ns,2,m) = i1tmn
         ibcinfo(n,ns,3,m) = i1tmx
         ibcinfo(n,ns,4,m) = i2tmn
         ibcinfo(n,ns,5,m) = i2tmx
         ibcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               ifo2 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               ifo2 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         ibcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0) then
            do mm=1,abs(ndata)
               bcvali(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6002 continue
c     
c     j0 boundary data
c
      if (ift.eq.3) then
c
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbcj0(n).gt.0 .and. ibctyp.eq.0) then
            m = 1
            if ((ibctyp.eq.jbcinfo(n,nbcj0(n),1,m)) .and.
     .          (i1tmn .eq.jbcinfo(n,nbcj0(n),2,m)) .and.
     .          (i1tmx .eq.jbcinfo(n,nbcj0(n),3,m)) .and.
     .          (i2tmn .eq.jbcinfo(n,nbcj0(n),4,m)) .and.
     .          (i2tmx .eq.jbcinfo(n,nbcj0(n),5,m))) then
                go to 6003
            end if
         end if
         nbcj0(n) = nbcj0(n) + 1
         ns = nbcj0(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 1
         jbcinfo(n,ns,1,m) = ibctyp
         jbcinfo(n,ns,2,m) = i2tmn
         jbcinfo(n,ns,3,m) = i2tmx
         jbcinfo(n,ns,4,m) = i1tmn
         jbcinfo(n,ns,5,m) = i1tmx
         jbcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               jfo1 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               jfo1 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         jbcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0) then
            do mm=1,abs(ndata)
               bcvalj(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6003 continue
c     
c     jdim boundary data
c
      if (ift.eq.4) then
c
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbcjdim(n).gt.0 .and. ibctyp.eq.0) then
            m = 2
            if ((ibctyp.eq.jbcinfo(n,nbcjdim(n),1,m)) .and.
     .          (i1tmn .eq.jbcinfo(n,nbcjdim(n),2,m)) .and.
     .          (i1tmx .eq.jbcinfo(n,nbcjdim(n),3,m)) .and.
     .          (i2tmn .eq.jbcinfo(n,nbcjdim(n),4,m)) .and.
     .          (i2tmx .eq.jbcinfo(n,nbcjdim(n),5,m))) then
                go to 6004
            end if
         end if
         nbcjdim(n) = nbcjdim(n) + 1
         ns = nbcjdim(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 2
         jbcinfo(n,ns,1,m) = ibctyp
         jbcinfo(n,ns,2,m) = i2tmn
         jbcinfo(n,ns,3,m) = i2tmx
         jbcinfo(n,ns,4,m) = i1tmn
         jbcinfo(n,ns,5,m) = i1tmx
         jbcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               jfo2 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               jfo2 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         jbcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0) then
            do mm=1,abs(ndata)
               bcvalj(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6004 continue
c     
c     k0 boundary data
c
      if (ift.eq.5) then
c
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbck0(n).gt.0 .and. ibctyp.eq.0) then
            m = 1
            if ((ibctyp.eq.kbcinfo(n,nbck0(n),1,m)) .and.
     .          (i1tmn .eq.kbcinfo(n,nbck0(n),2,m)) .and.
     .          (i1tmx .eq.kbcinfo(n,nbck0(n),3,m)) .and.
     .          (i2tmn .eq.kbcinfo(n,nbck0(n),4,m)) .and.
     .          (i2tmx .eq.kbcinfo(n,nbck0(n),5,m))) then
                go to 6005
            end if
         end if
         nbck0(n) = nbck0(n) + 1
         ns = nbck0(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 1
         kbcinfo(n,ns,1,m) = ibctyp
         kbcinfo(n,ns,2,m) = i1tmn
         kbcinfo(n,ns,3,m) = i1tmx
         kbcinfo(n,ns,4,m) = i2tmn
         kbcinfo(n,ns,5,m) = i2tmx
         kbcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               kfo1 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               kfo1 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         kbcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0) then
            do mm=1,abs(ndata)
               bcvalk(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6005 continue
c     
c     kdim boundary data
c
      if (ift.eq.6) then
c
c        don't include duplicate bc data - duplicate data can occur
c        when patch data is split
c
         if (nbckdim(n).gt.0 .and. ibctyp.eq.0) then
            m = 2
            if ((ibctyp.eq.kbcinfo(n,nbckdim(n),1,m)) .and.
     .          (i1tmn .eq.kbcinfo(n,nbckdim(n),2,m)) .and.
     .          (i1tmx .eq.kbcinfo(n,nbckdim(n),3,m)) .and.
     .          (i2tmn .eq.kbcinfo(n,nbckdim(n),4,m)) .and.
     .          (i2tmx .eq.kbcinfo(n,nbckdim(n),5,m))) then
                go to 6006
            end if
         end if
         nbckdim(n) = nbckdim(n) + 1
         ns = nbckdim(n)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
         m  = 2
         kbcinfo(n,ns,1,m) = ibctyp
         kbcinfo(n,ns,2,m) = i1tmn
         kbcinfo(n,ns,3,m) = i1tmx
         kbcinfo(n,ns,4,m) = i2tmn
         kbcinfo(n,ns,5,m) = i2tmx
         kbcinfo(n,ns,6,m) = imap(22,iseg,n)
         if (ivv.eq.0) then
            if (ibctyp.eq.1005 .or. ibctyp.eq.1006) then
               kfo2 = 1
            end if
         else
            if (ibctyp.eq.1004 .or. abs(ibctyp).eq.2004 .or.
     .          abs(ibctyp).eq.2014) then
               kfo2 = 1
            end if
         end if
         ndata = imap(21,iseg,n)
         kbcinfo(n,ns,7,m) = ndata
         if (abs(ndata) .gt. 0) then
            do mm=1,abs(ndata)
               bcvalk(n,ns,mm,m) = xmap(mm,iseg,n)
            end do
         end if
      end if
c
 6006 continue
c
 1010 continue 
c
c     set iforce for this block based on solid surface bc information
c
      ifo = 0
      jfo = 0
      kfo = 0
      if (ifo1.eq.1 .and. ifo2.eq.0) ifo = 1
      if (ifo1.eq.0 .and. ifo2.eq.1) ifo = 2
      if (ifo1.eq.1 .and. ifo2.eq.1) ifo = 3
      if (jfo1.eq.1 .and. jfo2.eq.0) jfo = 1
      if (jfo1.eq.0 .and. jfo2.eq.1) jfo = 2
      if (jfo1.eq.1 .and. jfo2.eq.1) jfo = 3
      if (kfo1.eq.1 .and. kfo2.eq.0) kfo = 1
      if (kfo1.eq.0 .and. kfo2.eq.1) kfo = 2
      if (kfo1.eq.1 .and. kfo2.eq.1) kfo = 3
      iforce(n) = ifo*100 + jfo*10 + kfo
c
c     extract 1-1 data 
c 
c     on the first pass, check for matching partners of 
c     1-1 segments and mark one as redundant
c
      do iseg=1,nseg(n)
c        t = target ("to") s = source ("from")
         if (iredundant(n,iseg).eq.0) then
         if (imap(1,iseg,n).eq.0 .or. imap(1,iseg,n).eq.1) then
            ift   = imap(2,iseg,n)
            i1tmn = imap(3,iseg,n)
            i1tmx = imap(4,iseg,n)
            i2tmn = imap(5,iseg,n)
            i2tmx = imap(6,iseg,n)
            ifs   = imap(8,iseg,n)
            i1smn = imap(9,iseg,n)
            i1smx = imap(10,iseg,n)
            i2smn = imap(11,iseg,n)
            i2smx = imap(12,iseg,n)
            ibt = n
            ibs = imap(7,iseg,n)
            do mseg=1,nseg(ibs)
c              following goto avoids a segment flagging itself
c              as redundant
               if (ibs.eq.ibt .and. mseg.eq.iseg) go to 1616
               if (imap(1,mseg,ibs).eq.0 .or. imap(1,mseg,ibs).eq.1
     .            .and. iredundant(ibs,mseg).eq.0) then
                  mft   = imap(2,mseg,ibs)
                  m1tmn = imap(3,mseg,ibs)
                  m1tmx = imap(4,mseg,ibs)
                  m2tmn = imap(5,mseg,ibs)
                  m2tmx = imap(6,mseg,ibs)
                  mfs   = imap(8,mseg,ibs)
                  m1smn = imap(9,mseg,ibs)
                  m1smx = imap(10,mseg,ibs)
                  m2smn = imap(11,mseg,ibs)
                  m2smx = imap(12,mseg,ibs)
                  mbt = ibs
                  mbs = imap(7,mseg,ibs)
c                 check if block numbers and face types are the same
c                 note: source face numbers may be negative
                  if (mbt.eq.ibs .and. mbs.eq.ibt .and.
     .                mft.eq.abs(ifs) .and. ift.eq.abs(mfs)) then
c                     rearrange indicies to run from lo to hi if needed
                      i1tmn_1 = i1tmn
                      i1tmx_1 = i1tmx
                      if (i1tmn.gt.i1tmx) then
                         i1tmx_1 = i1tmn
                         i1tmn_1 = i1tmx
                      end if
                      i2tmn_1 = i2tmn
                      i2tmx_1 = i2tmx
                      if (i2tmn.gt.i2tmx) then
                         i2tmx_1 = i2tmn
                         i2tmn_1 = i2tmx
                      end if
                      i1smn_1 = i1smn
                      i1smx_1 = i1smx
                      if (i1smn.gt.i1smx) then
                         i1smx_1 = i1smn
                         i1smn_1 = i1smx
                      end if
                      i2smn_1 = i2smn
                      i2smx_1 = i2smx
                      if (i2smn.gt.i2smx) then
                         i2smx_1 = i2smn
                         i2smn_1 = i2smx
                      end if
c
                      m1tmn_1 = m1tmn
                      m1tmx_1 = m1tmx
                      if (m1tmn.gt.m1tmx) then
                         m1tmx_1 = m1tmn
                         m1tmn_1 = m1tmx
                      end if
                      m2tmn_1 = m2tmn
                      m2tmx_1 = m2tmx
                      if (m2tmn.gt.m2tmx) then
                         m2tmx_1 = m2tmn
                         m2tmn_1 = m2tmx
                      end if
                      m1smn_1 = m1smn
                      m1smx_1 = m1smx
                      if (m1smn.gt.m1smx) then
                         m1smx_1 = m1smn
                         m1smn_1 = m1smx
                      end if
                      m2smn_1 = m2smn
                      m2smx_1 = m2smx
                      if (m2smn.gt.m2smx) then
                         m2smx_1 = m2smn
                         m2smn_1 = m2smx
                      end if
c
                      iflag = 0
                      iflag = iflag + (i1tmn_1 - m1smn_1)
                      iflag = iflag + (i1tmx_1 - m1smx_1)
                      iflag = iflag + (i2tmn_1 - m2smn_1)
                      iflag = iflag + (i2tmx_1 - m2smx_1)
                      if (iflag .eq. 0) then
                         iredundant(ibs,mseg) = 1
c                        write(6,*)'segment, block ',iseg,n,
c    .                   ' is redundant with'
c                        write(6,*)'segment, block ',mseg,ibs
                      end if
                  end if
               end if
 1616          continue
            end do
         end if
         end if
      end do
c
c     on this pass, the cfl3d 1-1 interface data is set from
c     the target and source sections of the non-redundant map data
c
      do 1020 iseg=1,nseg(n)
c     t = target ("to") s = source ("from")
      if (imap(1,iseg,n).eq.0 .or. imap(1,iseg,n).eq.1) then
      if (iredundant(n,iseg).eq.0) then
         nbli = nbli + 1
         nblon(nbli) = 0
         ift   = imap(2,iseg,n)
         i1tmn = imap(3,iseg,n)
         i1tmx = imap(4,iseg,n)
         i2tmn = imap(5,iseg,n)
         i2tmx = imap(6,iseg,n)
         ifs   = imap(8,iseg,n)
         i1smn = imap(9,iseg,n)
         i1smx = imap(10,iseg,n)
         i2smn = imap(11,iseg,n)
         i2smx = imap(12,iseg,n)
         nbt = n
         nbs = imap(7,iseg,n)
         nblk(1,nbli) = nbt
         nblk(2,nbli) = nbs
         if (ift.eq.1) then
            limblk(1,1,nbli) = 1
            limblk(1,4,nbli) = 1
            limblk(1,2,nbli) = i1tmn
            limblk(1,5,nbli) = i1tmx
            limblk(1,3,nbli) = i2tmn
            limblk(1,6,nbli) = i2tmx
            isva(1,1,nbli)   = 2
            isva(1,2,nbli)   = 3
         end if
         if (ift.eq.2) then
            limblk(1,1,nbli) = idimg(nbt)
            limblk(1,4,nbli) = idimg(nbt)
            limblk(1,2,nbli) = i1tmn
            limblk(1,5,nbli) = i1tmx
            limblk(1,3,nbli) = i2tmn
            limblk(1,6,nbli) = i2tmx
            isva(1,1,nbli)   = 2
            isva(1,2,nbli)   = 3
         end if
         if (ift.eq.3) then
            limblk(1,1,nbli) = i2tmn
            limblk(1,4,nbli) = i2tmx
            limblk(1,2,nbli) = 1
            limblk(1,5,nbli) = 1
            limblk(1,3,nbli) = i1tmn
            limblk(1,6,nbli) = i1tmx
            isva(1,1,nbli)   = 1
            isva(1,2,nbli)   = 3
         end if
         if (ift.eq.4) then
            limblk(1,1,nbli) = i2tmn
            limblk(1,4,nbli) = i2tmx
            limblk(1,2,nbli) = jdimg(nbt)
            limblk(1,5,nbli) = jdimg(nbt)
            limblk(1,3,nbli) = i1tmn
            limblk(1,6,nbli) = i1tmx
            isva(1,1,nbli)   = 1
            isva(1,2,nbli)   = 3
         end if
         if (ift.eq.5) then
            limblk(1,1,nbli) = i1tmn
            limblk(1,4,nbli) = i1tmx
            limblk(1,2,nbli) = i2tmn
            limblk(1,5,nbli) = i2tmx
            limblk(1,3,nbli) = 1
            limblk(1,6,nbli) = 1
            isva(1,1,nbli)   = 1
            isva(1,2,nbli)   = 2
         end if
         if (ift.eq.6) then
            limblk(1,1,nbli) = i1tmn
            limblk(1,4,nbli) = i1tmx
            limblk(1,2,nbli) = i2tmn
            limblk(1,5,nbli) = i2tmx
            limblk(1,3,nbli) = kdimg(nbt)
            limblk(1,6,nbli) = kdimg(nbt)
            isva(1,1,nbli)   = 1
            isva(1,2,nbli)   = 2
         end if
         if (abs(ifs).eq.1) then
            limblk(2,1,nbli) = 1
            limblk(2,4,nbli) = 1
            limblk(2,2,nbli) = i1smn
            limblk(2,5,nbli) = i1smx
            limblk(2,3,nbli) = i2smn
            limblk(2,6,nbli) = i2smx
            isva(2,1,nbli)   = 2
            isva(2,2,nbli)   = 3
         end if
         if (abs(ifs).eq.2) then
            limblk(2,1,nbli) = idimg(nbs)
            limblk(2,4,nbli) = idimg(nbs)
            limblk(2,2,nbli) = i1smn
            limblk(2,5,nbli) = i1smx
            limblk(2,3,nbli) = i2smn
            limblk(2,6,nbli) = i2smx
            isva(2,1,nbli)   = 2
            isva(2,2,nbli)   = 3 
         end if
         if (abs(ifs).eq.3) then
            limblk(2,1,nbli) = i2smn
            limblk(2,4,nbli) = i2smx
            limblk(2,2,nbli) = 1
            limblk(2,5,nbli) = 1
            limblk(2,3,nbli) = i1smn
            limblk(2,6,nbli) = i1smx
            isva(2,1,nbli)   = 1
            isva(2,2,nbli)   = 3
         end if
         if (abs(ifs).eq.4) then
            limblk(2,1,nbli) = i2smn
            limblk(2,4,nbli) = i2smx
            limblk(2,2,nbli) = jdimg(nbs)
            limblk(2,5,nbli) = jdimg(nbs)
            limblk(2,3,nbli) = i1smn
            limblk(2,6,nbli) = i1smx
            isva(2,1,nbli)   = 1
            isva(2,2,nbli)   = 3
         end if
         if (abs(ifs).eq.5) then
            limblk(2,1,nbli) = i1smn
            limblk(2,4,nbli) = i1smx
            limblk(2,2,nbli) = i2smn
            limblk(2,5,nbli) = i2smx
            limblk(2,3,nbli) = 1
            limblk(2,6,nbli) = 1
            isva(2,1,nbli)   = 1
            isva(2,2,nbli)   = 2
         end if
         if (abs(ifs).eq.6) then
            limblk(2,1,nbli) = i1smn
            limblk(2,4,nbli) = i1smx
            limblk(2,2,nbli) = i2smn
            limblk(2,5,nbli) = i2smx
            limblk(2,3,nbli) = kdimg(nbs)
            limblk(2,6,nbli) = kdimg(nbs)
            isva(2,1,nbli)   = 1
            isva(2,2,nbli)   = 2
         end if
c
         if (ifs.lt.0) then
            itemp = isva(2,1,nbli)
            isva(2,1,nbli) = isva(2,2,nbli)
            isva(2,2,nbli) = itemp
         end if
c
c        check order of side 1 against the order of side 2,
c        set via isva array
c
         irange1 = abs(limblk(1,4,nbli) - limblk(1,1,nbli))
         irange2 = abs(limblk(2,4,nbli) - limblk(2,1,nbli))
         jrange1 = abs(limblk(1,5,nbli) - limblk(1,2,nbli))
         jrange2 = abs(limblk(2,5,nbli) - limblk(2,2,nbli))
         krange1 = abs(limblk(1,6,nbli) - limblk(1,3,nbli))
         krange2 = abs(limblk(2,6,nbli) - limblk(2,3,nbli))
         if (isva(1,1,nbli) .eq. 1) then
            llrange11 = irange1
         else if (isva(1,1,nbli) .eq. 2) then
            llrange11 = jrange1
         else if (isva(1,1,nbli) .eq. 3) then
            llrange11 = krange1
         end if
         if (isva(1,2,nbli) .eq. 1) then
            llrange12 = irange1
         else if (isva(1,2,nbli) .eq. 2) then
            llrange12 = jrange1
         else if (isva(1,2,nbli) .eq. 3) then
            llrange12 = krange1
         end if
         if (isva(2,1,nbli) .eq. 1) then
            llrange21 = irange2
         else if (isva(2,1,nbli) .eq. 2) then
            llrange21 = jrange2
         else if (isva(2,1,nbli) .eq. 3) then
            llrange21 = krange2
         end if
         if (isva(2,2,nbli) .eq. 1) then
            llrange22 = irange2
         else if (isva(2,2,nbli) .eq. 2) then
            llrange22 = jrange2
         else if (isva(2,2,nbli) .eq. 3) then
            llrange22 = krange2
         end if
         if (llrange11.ne.llrange21 .or. llrange12.ne.llrange22) then
            itemp = isva(2,1,nbli)
            isva(2,1,nbli) = isva(2,2,nbli)
            isva(2,2,nbli) = itemp
         end if
      end if
      end if
 1020 continue
c
c     extract patched interface data
c
      nbtlast = 0
      iftlast = 0
      i1tmnlast = 0
      i1tmxlast = 0
      i2tmnlast = 0
      i2tmxlast = 0
      do 1030 iseg=1,nseg(n)
c     t = target ("to") s = source ("from")
      if (imap(1,iseg,n) .eq. -1) then
         ift    = imap(2,iseg,n)
         i1tmn  = imap(3,iseg,n)
         i1tmx  = imap(4,iseg,n)
         i2tmn  = imap(5,iseg,n)
         i2tmx  = imap(6,iseg,n)
         ifs    = abs(imap(8,iseg,n))
         i1smn  = imap(9,iseg,n)
         i1smx  = imap(10,iseg,n)
         i2smn  = imap(11,iseg,n)
         i2smx  = imap(12,iseg,n)
c
c        reverse order on source zone indicies to go from
c        low to high, if needed
c
         if (i1smn.gt.i1smx) then
            itemp = i1smn
            i1smn = i1smx
            i1smx = itemp
         end if
         if (i2smn.gt.i2smx) then
            itemp = i2smn
            i2smn = i2smx
            i2smx = itemp
         end if
c
         nbt    = n
         nbs    = imap(7,iseg,n)
c
c        check to see if this is the same target ("to") face as the
c        last segment; if it is, then do not increment the ninter
c        counter (no. of patched interfaces), but just add to the
c        list of source ("from") faces
c
         ichk = 0
         ichk = ichk + abs(nbt   - nbtlast)
         ichk = ichk + abs(ift   - iftlast)
         ichk = ichk + abs(i1tmn - i1tmnlast) 
         ichk = ichk + abs(i1tmx - i1tmxlast)
         ichk = ichk + abs(i2tmn - i2tmnlast)
         ichk = ichk + abs(i2tmx - i2tmxlast)
c
c        target ("to") data
c
         if (ichk.ne.0) then
            ninter      = ninter + 1
            nfb(ninter) = 0
            ito1        = nbt
            if (ift.eq.1 .or.ift.eq.2) then
               ito2 = 11
               if (ift.eq.2) ito2 = 12
               xit1(ninter) = i1tmn
               xit2(ninter) = i1tmx
               ett1(ninter) = i2tmn
               ett2(ninter) = i2tmx
               ito(ninter)  = 100*ito1 + ito2
            end if
            if (ift.eq.3 .or.ift.eq.4) then
               ito2 = 21
               if (ift.eq.4) ito2 = 22
               xit1(ninter) = i1tmn
               xit2(ninter) = i1tmx
               ett1(ninter) = i2tmn
               ett2(ninter) = i2tmx
               ito(ninter)  = 100*ito1 + ito2
            end if
            if (ift.eq.5 .or.ift.eq.6) then
               ito2 = 31
               if (ift.eq.6) ito2 = 32
               xit1(ninter) = i2tmn
               xit2(ninter) = i2tmx
               ett1(ninter) = i1tmn
               ett2(ninter) = i1tmx
               ito(ninter)  = 100*ito1 + ito2
            end if
         end if
c
c        source ("from") data
c
         nfb(ninter) = nfb(ninter) + 1
         nfbb        = nfb(ninter)
         ifrom1 = nbs
         if (ifs.eq.1 .or.ifs.eq.2) then
            ifrom2 = 11
            if (ifs.eq.2) ifrom2 = 12
c
c           check bounds of source zone to make sure they lie in
c           the range of dimensions of the source zone - an early  
c           version of a tlns3d block splitter code would introduce
c           errors in the tlns3d map file wherein some source zone
c           interface ranges actually were larger than the source
c           zone dimensions
c
            if (i1smx.gt.jdimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i1smx>jdimg(nbs): i1smx,jdimg(nbs) = ',
     .         i1smx,jdimg(nbs) 
               i1smx = jdimg(nbs)
            end if
            if (i2smx.gt.kdimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i2smx>kdimg(nbs): i2smx,kdimg(nbs) = ',
     .         i2smx,kdimg(nbs)
               i2smx = kdimg(nbs)
            end if
c
            xif1(ninter,nfbb) = i1smn
            xif2(ninter,nfbb) = i1smx
            etf1(ninter,nfbb) = i2smn
            etf2(ninter,nfbb) = i2smx
            ifrom(ninter,nfbb)  = 100*ifrom1 + ifrom2
         end if
         if (ifs.eq.3 .or.ifs.eq.4) then
            ifrom2 = 21
c
c           check bounds of source zone to make sure they lie in
c           the range of dimensions of the source zone - an early
c           version of a tlns3d block splitter code would introduce
c           errors in the tlns3d map file wherein some source zone
c           interface ranges actually were larger than the source
c           zone dimensions
c
            if (i1smx.gt.kdimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i1smx>kdimg(nbs): i1smx,kdimg(nbs) = ',
     .         i1smx,kdimg(nbs)
               i1smx = kdimg(nbs)
            end if
            if (i2smx.gt.idimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i2smx>idimg(nbs): i2smx,idimg(nbs) = ',
     .         i2smx,idimg(nbs)
               i2smx = idimg(nbs)
            end if
c
            if (ifs.eq.4) ifrom2 = 22
            xif1(ninter,nfbb) = i1smn
            xif2(ninter,nfbb) = i1smx
            etf1(ninter,nfbb) = i2smn
            etf2(ninter,nfbb) = i2smx
            ifrom(ninter,nfbb)  = 100*ifrom1 + ifrom2
         end if
         if (ifs.eq.5 .or.ifs.eq.6) then
            ifrom2 = 31
c
c           check bounds of source zone to make sure they lie in
c           the range of dimensions of the source zone - an early
c           version of a tlns3d block splitter code would introduce
c           errors in the tlns3d map file wherein some source zone
c           interface ranges actually were larger than the source
c           zone dimensions
c
            if (i2smx.gt.jdimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i2smx>jdimg(nbs): i2smx,jdimg(nbs) = ',
     .         i2smx,jdimg(nbs)
               i2smx = jdimg(nbs)
            end if
            if (i1smx.gt.idimg(nbs)) then
               write(6,*)'correcting error in tlns3d map file:'
               write(6,*)'i1smx>idimg(nbs): i1smx,idimg(nbs) = ',
     .         i1smx,idimg(nbs)
               i1smx = idimg(nbs)
            end if
c
            if (ifs.eq.6) ifrom2 = 32
            xif1(ninter,nfbb) = i2smn
            xif2(ninter,nfbb) = i2smx
            etf1(ninter,nfbb) = i1smn
            etf2(ninter,nfbb) = i1smx
            ifrom(ninter,nfbb)  = 100*ifrom1 + ifrom2
         end if
         nbtlast = nbt
         iftlast = ift
         i1tmnlast = i1tmn
         i1tmxlast = i1tmx
         i2tmnlast = i2tmn
         i2tmxlast = i2tmx
      end if
 1030 continue
c
 1000 continue
c
c     further concatenate any applicable patches. The logic above 
c     will concatenate patches to any target segment that have the
c     same target ("to") index range *and* occur sequentially
c     in the tlns3d map file. The problem arises if the patches
c     are not in sequence in the mapfile. In that case, the
c     target is split into multiple segments, and individual segments
c     may not be completely covered by source ("from") segments
c
c
c     idebug > 0 for debugging purposes to write out patch 
c     data before elimination of split patches and to indicate
c     which patches were determined to need concatenation
c
      idebug = 0
c
      if(idebug.gt.0) then
         write(21,*)'first pass through patch interfaces'
         do 2416 int=1,ninter
         write(21,201) int,ito(int),xit1(int),xit2(int),ett1(int),
     .   ett2(int),nfb(int)
         do 2418 l = 1,nfb(int)
         write(21,203) ifrom(int,l),xif1(int,l),xif2(int,l),
     .   etf1(int,l),etf2(int,l)
2418     continue
2416     continue
      end if
c
      ninter0 = ninter
c
      do 1990 n = 1,ninter0
      iconcat(n) = 0
 1990 continue
c
      ninter = 0
      do 2000 n=1,ninter0
      if (iconcat(n).eq.0) then
         ninter = ninter + 1
         xit1n = xit1(n)
         xit2n = xit2(n)
         ett1n = ett1(n)
         ett2n = ett2(n)
         iton  = ito(n) 
         do 2100 m=1,ninter0
         if (n.ne.m) then
            xit1m = xit1(m)
            xit2m = xit2(m)
            ett1m = ett1(m)
            ett2m = ett2(m)
            itom  = ito(m)
            ichk = 0
            ichk = ichk + abs(xit1m-xit1n)
            ichk = ichk + abs(xit2m-xit2n)
            ichk = ichk + abs(ett1m-ett1n)
            ichk = ichk + abs(ett2m-ett2n)
            ichk = ichk + abs(itom-iton)
            if (ichk .eq.0) then
c              patches m and n can be concatenated
               nfbm   = nfb(m)
               do 2110 nf = 1,nfbm
               nfb(n) = nfb(n)+1
               nfbb   = nfb(n)
               xif1(n,nfbb)  = xif1(m,nf) 
               xif2(n,nfbb)  = xif2(m,nf)
               etf1(n,nfbb)  = etf1(m,nf)
               etf2(n,nfbb)  = etf2(m,nf)
               ifrom(n,nfbb) = ifrom(m,nf)
 2110          continue 
               iconcat(m) = 1
               if (idebug.gt.0) write(21,*)' patch ',n,
     .            ' will concatenate with patch ',m
            end if
         end if
 2100    continue
      end if
 2000 continue
      if (idebug.gt.0) then
         write(21,*)'old ninter = ',ninter0
         write(21,*)'new ninter = ',ninter
      end if
c
c     set default patching parameters iifit,llimit,iitmax,iico,iiorph
c
      if (ninter .gt. 0) then
         do 998 n=1,ninter
         iifit(n)  = 1
         llimit(n) = 1
         iitmax(n) = 100
         iic0(n)   = 0
         iiorph(n) = 0
c        for now, set mmcxie and mmceta = 50
c
         mmcxie(n) = 50
         mmceta(n) = 50
  998    continue
      end if
c
c     check multigridability of boundary-condition segments
c
      ncgmax2 = 0
      iflg = 1
      do nnnn = 1,ncgmax+1
         nnn = nnnn-1
         do igrid=1,ngrid
            do m=1,2
               if (m.eq.1) then
                  ns = nbci0(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               else
                  ns = nbcidim(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               end if
               do iseg=1,ns
                  jtest1 = (ibcinfo(igrid,iseg,2,m)-1)/2**nnn + 1
                  jtest2 = (ibcinfo(igrid,iseg,3,m)-1)/2**nnn + 1
                  ktest1 = (ibcinfo(igrid,iseg,4,m)-1)/2**nnn + 1 
                  ktest2 = (ibcinfo(igrid,iseg,5,m)-1)/2**nnn + 1
                  if ((jtest1/2*2 .eq. jtest1) .or.
     .                (ktest1/2*2 .eq. ktest1) .or.
     .                (jtest2/2*2 .eq. jtest2) .or.
     .                (ktest2/2*2 .eq. ktest2)) then
                      iflg = 0
                  end if
               end do
            end do
c
            do m=1,2
               if (m.eq.1) then
                  ns = nbcj0(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               else
                  ns = nbcjdim(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               end if
               do iseg=1,ns
                  itest1 = (jbcinfo(igrid,iseg,2,m)-1)/2**nnn + 1
                  if (idimg(igrid) .eq. 2) then
                     itest2 = 1
                  else
                     itest2 = (jbcinfo(igrid,iseg,3,m)-1)/2**nnn + 1
                  end if
                  ktest1 = (jbcinfo(igrid,iseg,4,m)-1)/2**nnn + 1
                  ktest2 = (jbcinfo(igrid,iseg,5,m)-1)/2**nnn + 1
                  if ((itest1/2*2 .eq. itest1) .or.
     .                (ktest1/2*2 .eq. ktest1) .or.
     .                (itest2/2*2 .eq. itest2) .or.
     .                (ktest2/2*2 .eq. ktest2)) then
                      iflg = 0
                  end if
               end do
            end do
c
            do m=1,2
               if (m.eq.1) then
                  ns = nbck0(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               else
                  ns = nbckdim(igrid)
                  if (ns.gt.maxseg) then
                     write(6,*)'stopping...must increase ',
     .               'parameter maxseg to at least ',ns
                     call termn8(0,-3,ibufdim,nbuf,bou,nou)
                  end if
               end if
               do iseg=1,ns
                  itest1 = (kbcinfo(igrid,iseg,2,m)-1)/2**nnn + 1
                  if (idimg(igrid) .eq. 2) then
                     itest2 = 1
                  else
                     itest2 = (kbcinfo(igrid,iseg,3,m)-1)/2**nnn + 1
                  end if
                  ktest1 = (kbcinfo(igrid,iseg,4,m)-1)/2**nnn + 1
                  ktest2 = (kbcinfo(igrid,iseg,5,m)-1)/2**nnn + 1
                  if ((itest1/2*2 .eq. itest1) .or.
     .                (ktest1/2*2 .eq. ktest1) .or.
     .                (itest2/2*2 .eq. itest2) .or.
     .                (ktest2/2*2 .eq. ktest2)) then
                      iflg = 0
                  end if
               end do
            end do
         end do
         if (iflg .ne. 0) ncgmax2 = ncgmax2 + 1
      end do
c
      write(*,'(1x)')
      write(*,'("bc segment dimensions are multigridable ",
     .          "to ncg = ",i2)') ncgmax2
      write(*,'(1x)')
c
      ncgmax = min(ncgmax,ncgmax2)
c
      write(*,'("from block segment dimensions and ",
     .          "overall block dimensions, ncgmax = ",i2)') 
     .          ncgmax
      icgflg = 0
      do igrd=1,ngrid
         if(ncgmax.lt.ncgg(igrd)) icgflg = 1
      end do
      if (icgflg.gt.0) then
         write(*,'("this value will be used for ncg in the",
     .   " split input file")')
         write(*,'(1x)')
         write(*,'("   NOTE: this is a smaller ncg than used in",
     .   " the unsplit input file")')
        write(*,'("   you may want to consider an alternative",
     .   " splitting to maintain")')
        write(*,'("   the same level of multigridability")')
         do igrd=1,ngrid
            ncgg(igrd) = ncgmax
         end do
      end if
      icgflg2 = 0
      do igrd=1,ngrid
         if(ncgmax.gt.ncgg(igrd)) icgflg2 = 1
      end do
      if (icgflg2.gt.0) then
         write(*,'(1x)')
         write(*,'("   NOTE: this is a larger ncg than used in ",
     .   "the unsplit input file")')
         write(*,'("   the unsplit ncg value will be maintained ",
     .   "but you might want to")')
         write(*,'("   increase the value of ncg and alter the ",
     .   " mulitgrid parameters")')
         write(*,'("   to reflect this level of multigridability")')
      end if
      write(*,'(1x)')
c
      ncg = ncgmax
c
c***********************************************************************
c
c     begin generation of cfl3d input file
c
c***********************************************************************
c
c
c     set wall temp data for bc2004 in version 5 based on isnd and c2spe
c
      if (isnd.eq.0) then
         twall = 0
      else
         if (c2spe.gt.0) then
            twall = c2spe
         else
            twall = -1
         end if
      end if
c
      write(7,18)
   18 format(19hinput/output files:)
c
c     grid file name already read in is used
c     - use generic names for other files
c
      write(7,'(a60)')gridout
      write(7,'(a60)')plt3dg
      write(7,'(a60)')plt3dq
      write(7,'(a60)')output
      write(7,'(a60)')resid
      write(7,'(a60)')turbres
      write(7,'(a60)')blomx
      write(7,'(a60)')output2
      write(7,'(a60)')printout
      write(7,'(a60)')pplunge
      write(7,'(a60)')ovrlap
      write(7,'(a60)')patch
      write(7,'(a60)')restrt
c
c     check for keyword input
c
      nskip = 14
      do n=1,nskip
         read(10,*)
      end do
c
      iunit5sav = iunit5
      iunit5    = 10
      call readkey(ititr,myid,ibufdim,nbuf,bou,nou,7,-1)
      iunit5    = iunit5sav
c
      write(7,11)(title(i),i=1,20)
   11 format(1h ,20a4)
c
c     versions 4.1 and 6/5 have different ways of triggering the
c     Cl-specified option...v4.1 uses the sign on Mach no. as a
c     trigger, whereas v6 uses the sign on Sref (in v6, the sign
c     on Mach no. is used to trigger the preconditioning option;
c     there is no preconditioning option in v4.1). The location
c     of the extra data for the specification of Cl also differs
c     between the two versions.
c
      if (ialphit .gt. 0) then
         if (iver .eq. 4) then
            xmach = -abs(xmach)
         else
            sref  = -abs(sref)
         end if
      end if
      if (iprecon .gt. 0) then
         if (iver .gt. 4) then
            xmach = -abs(xmach)
         end if
      end if
c
      if (iver.eq.4) then
         if (ipar.eq.0) then
         write(7,20)
   20    format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .   6x,4hisnd,5x,5hc2spe)
         write(7,21) xmach,alpha,beta,reue,tinf,isnd,c2spe
         else
         write(7,520)
  520    format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .   6x,4hisnd,5x,5hc2spe,5x,5hnodes)
         write(7,21) xmach,alpha,beta,reue,tinf,isnd,c2spe,nnodes
         end if
   21    format(3f10.5,e10.3,f10.5,i10,f10.5,i10)
      else
         write(7,22)
   22    format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .   5x,5hialph,4x,6hihstry) 
         write(7,23) xmach,alpha,beta,reue,tinf,ialph,ihstry
   23    format(3f10.5,e10.3,f10.5,i10,i10,i10)
      end if
c
      xmach = abs(xmach)
c
      if (iver .eq. 4) then
         if (ialphit .eq. 0) then
            write(7,25)
   25       format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,
     .      7x,3hzmc)
            write(7,26)sref,cref,bref,xmc,ymc,zmc
         else
            write (7,251)
  251       format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,
     .      7x,3hzmc,4x,6hcltarg,3x,7hresupdt)
            write (7,26) sref,cref,bref,xmc,ymc,zmc,cltarg,resupdt
         end if
   26    format(f10.3,7f10.4)
      else
         write(7,25)
         write(7,26)sref,cref,bref,xmc,ymc,zmc
      end if
c
      if (iver.eq.4) then
         write(7,31)
   31    format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,5hiunst,
     .   5x,5hrfreq,4x,6halphau,6x,4hcloc)
         write(7,32)dt,irest,iflagts,fmax,iunst,rfreq,alphau,cloc
   32    format(f10.5,2i10,f10.5,i10,3f10.5)
      else
         write(7,33)
   33    format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,5hiunst,
     .   3x,7hcfl_tau)
         write(7,34)dt,irest,iflagts,fmax,iunst,cfltau
   34    format(f10.5,2i10,f10.5,i10,3f10.5)
      end if
c
      ngridout = ngrid
      if (ip3dgrd .gt. 0) ngridout = -ngrid
c
c     in version 4.1, assume no plot3d files are output; 
c     in version 6/5 and 4.1hp, output solid sufaces to the
c     plot3d files; don't output print files in any version 
c
      if (iver.eq.4 .and. ipar.eq.0) then
         nplot3d = 0
         nprint  = 0
      else
         if(abs(nplot3d) .gt. 0) nplot3d = -1
         if(abs(nprint)  .gt. 0) nprint  = 0
      end if
c
      nwrest = 9999
      ichk   = 0
c
      write(7,35)
   35 format(5x,5hngrid,3x,7hnplot3d,4x,6hnprint,4x,6hnwrest,
     .       6x,4hichk,7x,3hi2d,3x,7h ntstep,4x,6h   ita)
      write(7,36) ngridout,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
   36 format(10i10)
c
      write(7,61)
   61 format(7x,3hncg,7x,3hiem,2x,8hiadvance,4x,6hiforce,
     .       2x,8hivisc(i),2x,8hivisc(j),2x,8hivisc(k))
      do 71 igrid=1,ngrid
   71 write(7,36) ncgg(igrid),iemg(igrid),iadvance(igrid),
     .            iforce(igrid),(iviscg(igrid,l),l=1,3)
c
      write(7,73)
   73 format(6x,4hidim,6x,4hjdim,6x,4hkdim)
      i2dd = 0
      do 7001 igrid=1,ngrid
      if (i2d.gt.0 .and. idimg(igrid).ne.2) i2dd = 1
      write(7,78) idimg(igrid),jdimg(igrid),kdimg(igrid)
   78 format(3i10)
 7001 continue
c
      write(7,80)
   80 format(4x,6hilamlo,4x,6hilamhi,4x,6hjlamlo,
     .       4x,6hjlamhi,4x,6hklamlo,4x,6hklamhi)
      do 7201 igrid=1,ngrid
      write(7,36) ilamlog(igrid),ilamhig(igrid),jlamlog(igrid),
     .            jlamhig(igrid),klamlog(igrid),klamhig(igrid)
 7201 continue
c
      write(7,82)
   82 format(5x,5hinewg,4x,6higridc,
     .       8x,2his,8x,2hjs,8x,2hks,8x,2hie,8x,2hje,8x,2hke)
      do 7002 igrid=1,ngrid
      inewg  = inewgg(igrid)
      igridc = igridg(igrid)
      is     = isg(igrid)
      js     = jsg(igrid)
      ks     = ksg(igrid)
      ie     = ieg(igrid)
      je     = jeg(igrid)
      ke     = keg(igrid)
      write(7,36) inewg,igridc,is,js,ks,ie,je,ke
 7002 continue
c
      write(7,84)
   84 format(2x,8hidiag(i),2x,8hidiag(j),2x,8hidiag(k),
     .       2x,8hiflim(i),2x,8hiflim(j),2x,8hiflim(k))
      do 9002 igrid=1,ngrid
      write(7,36) idiagg(igrid,1),idiagg(igrid,2),idiagg(igrid,3),
     .            iflimg(igrid,1),iflimg(igrid,2),iflimg(igrid,3)
 9002 continue
c
      write(7,86)
   86 format(3x,7hifds(i),3x,7hifds(j),3x,7hifds(k),
     .       2x,8hrkap0(i),2x,8hrkap0(j),2x,8hrkap0(k))
      do 7012 igrid=1,ngrid
      write(7,88) ifdsg(igrid,1),ifdsg(igrid,2),ifdsg(igrid,3),
     .            rkap0g(igrid,1),rkap0g(igrid,2),rkap0g(igrid,3)
   88 format(3i10,3f10.4)
 7012 continue
c
      write(7,90)
   90 format(6x,4hgrid,5x,5hnbci0,3x,7hnbcidim,5x,5hnbcj0,3x,7hnbcjdim,
     .                 5x,5hnbck0,3x,7hnbckdim,4x,6hiovrlp)
      do 9007 igrid=1,ngrid
      write(7,36) igrid,nbci0(igrid),nbcidim(igrid),nbcj0(igrid),
     .            nbcjdim(igrid),nbck0(igrid),nbckdim(igrid),
     .            iovrlp(igrid)
 9007 continue
c
c     I0/IDIM Boundary Data
c
      do 1209 m=1,2
      if (m.eq.1) then
         write(7,102)
      else
         write(7,103)
      end if
      do 1207 igrid=1,ngrid
      if (m.eq.1) then
         ns = nbci0(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      else
         ns = nbcidim(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      end if
      do 1201 iseg=1,ns
      mseg = iseg
      if (ibcinfo(igrid,iseg,6,m) .eq. 0) then
         mseg = -iseg
      end if
      if (iver.ne.4 .and. ibcinfo(igrid,iseg,1,m).eq.1004) then
         ibcinfo(igrid,iseg,1,m) = 2004
         ibcinfo(igrid,iseg,7,m) = 2
         bcvali(igrid,iseg,1,m)  = twall
         bcvali(igrid,iseg,2,m)  = 0. 
      end if
      write(7,36) igrid,mseg,ibcinfo(igrid,iseg,1,m),
     .                       ibcinfo(igrid,iseg,2,m),
     .                       ibcinfo(igrid,iseg,3,m),
     .                       ibcinfo(igrid,iseg,4,m),
     .                       ibcinfo(igrid,iseg,5,m),
     .                       ibcinfo(igrid,iseg,7,m)
      if (ibcinfo(igrid,iseg,7,m).gt.0) then
         ndata  = ibcinfo(igrid,iseg,7,m)
         ibctyp = ibcinfo(igrid,iseg,1,m)
         call getdhdr(datahdr,ibctyp,iver,ndata)
         write(7,104) (datahdr(mm),mm=1,ndata)
         write(7,105) (bcvali(igrid,iseg,mm,m),mm=1,ndata)
      end if
      if (ibcinfo(igrid,iseg,7,m).lt.0) then
         write(7,1213)
         ibcfile = bcvali(igrid,iseg,1,m)
         write(7,'(''bc_file'',i2)') ibcfile
      end if
 1213 format(17h     bc data file)

c
 1201 continue
 1207 continue
 1209 continue
c
  102 format(3hi0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .             6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
  103 format(5hidim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .               6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
  104 format(10a10)
  105 format(10f10.5)
c
c     J0/JDIM Boundary Data
c
      do 1309 m=1,2
      if (m.eq.1) then
         write(7,112)
      else
         write(7,113)
      end if
      do 1307 igrid=1,ngrid
      if (m.eq.1) then
         ns = nbcj0(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      else
         ns = nbcjdim(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      end if
      do 1301 iseg=1,ns
      mseg = iseg
      if (jbcinfo(igrid,iseg,6,m) .eq. 0) then
         mseg = -iseg
      end if
      if (iver.ne.4 .and. jbcinfo(igrid,iseg,1,m).eq.1004) then
         jbcinfo(igrid,iseg,1,m) = 2004
         jbcinfo(igrid,iseg,7,m) = 2
         bcvalj(igrid,iseg,1,m)  = twall
         bcvalj(igrid,iseg,2,m)  = 0.
      end if
      write(7,36) igrid,mseg,jbcinfo(igrid,iseg,1,m),
     .                       jbcinfo(igrid,iseg,2,m),
     .                       jbcinfo(igrid,iseg,3,m),
     .                       jbcinfo(igrid,iseg,4,m),
     .                       jbcinfo(igrid,iseg,5,m),
     .                       jbcinfo(igrid,iseg,7,m)
      if (jbcinfo(igrid,iseg,7,m).gt.0) then
         ndata = jbcinfo(igrid,iseg,7,m)
         ibctyp = jbcinfo(igrid,iseg,1,m)
         call getdhdr(datahdr,ibctyp,iver,ndata)
         write(7,104) (datahdr(mm),mm=1,ndata)
         write(7,105) (bcvalj(igrid,iseg,mm,m),mm=1,ndata)
      end if
      if (jbcinfo(igrid,iseg,7,m).lt.0) then
         write(7,1213)
         ibcfile = bcvalj(igrid,iseg,1,m)
         write(7,'(''bc_file'',i2)') ibcfile
      end if
c
 1301 continue
 1307 continue
 1309 continue
c
  112 format(3hj0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
  113 format(5hjdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .               6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
c
c     K0/KDIM Boundary Data
c
      do 1409 m=1,2
      if (m.eq.1) then
         write(7,122)
      else
         write(7,123)
      end if
      do 1407 igrid=1,ngrid
      if (m.eq.1) then
         ns = nbck0(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      else
         ns = nbckdim(igrid)
         if (ns.gt.maxseg) then
            write(6,*)'stopping...must increase ',
     .      'parameter maxseg to at least ',ns
            call termn8(0,-3,ibufdim,nbuf,bou,nou)
         end if
      end if
      do 1401 iseg=1,ns
      mseg = iseg
      if (kbcinfo(igrid,iseg,6,m) .eq. 0) then
         mseg = -iseg
      end if
      if (iver.ne.4 .and. kbcinfo(igrid,iseg,1,m).eq.1004) then
         kbcinfo(igrid,iseg,1,m) = 2004
         kbcinfo(igrid,iseg,7,m) = 2
         bcvalk(igrid,iseg,1,m)  = twall
         bcvalk(igrid,iseg,2,m)  = 0. 
      end if
      write(7,36) igrid,mseg,kbcinfo(igrid,iseg,1,m),
     .                       kbcinfo(igrid,iseg,2,m),
     .                       kbcinfo(igrid,iseg,3,m),
     .                       kbcinfo(igrid,iseg,4,m),
     .                       kbcinfo(igrid,iseg,5,m),
     .                       kbcinfo(igrid,iseg,7,m)
      if (kbcinfo(igrid,iseg,7,m).gt.0) then
         ndata = kbcinfo(igrid,iseg,7,m)
         ibctyp = kbcinfo(igrid,iseg,1,m)
         call getdhdr(datahdr,ibctyp,iver,ndata)
         write(7,104) (datahdr(mm),mm=1,ndata)
         write(7,105) (bcvalk(igrid,iseg,mm,m),mm=1,ndata)
      end if
      if (kbcinfo(igrid,iseg,7,m).lt.0) then
         write(7,1213)
         ibcfile = bcvalk(igrid,iseg,1,m)
         write(7,'(''bc_file'',i2)') ibcfile
      end if
c
 1401 continue
 1407 continue
 1409 continue
c
  122 format(3hk0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
  123 format(5hkdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .               6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
c
c     check to see if ncg has been lowered due to the splitting
c     process, and, if so, see if mseq, mgflag, etc. need to be 
c     altered
c
      if (icgflg .ne. 0) then
        if (ncgmax .eq. 0) then
           mgflag = 0
           mseq   = 1
        end if
        if (mseq .eq. 1) then
           if (mgflag .ne. 0) then
c             set mglevg to upper limit consistant with new ncg
              mglev(1) = ncgmax + 1
              do ii = 1,mglev(1)
                 mitl(ii,1) = 1
              end do
           else
              mglev(1) = 1
              mitl(1,1) = 1
           end if
        else
c          set mseq to upper limit consistant with new ncg
           if (mseq .gt. ncgmax + 1) mseq = ncgmax + 1
           if (mgflag .ne. 0) then
c             set mglevg to upper limit consistant with new ncg
              do mm = 1,mseq
                 mglev(mm) = mm
                 do ii = 1,mglev(mm)
                    mitl(ii,mm) = 1
                 end do
              end do
           else
              do mm = 1,mseq
                 mglev(mm) = 1
                 mitl(1,mm) = 1 
              end do
           end if
        end if
      end if
c
      if (iver .eq. 4) then
         write(7,130)
  130    format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam)     
         write(7,36) mseq,mgflag,iconsf,mtt,ngam
      else
         if (iprecon .eq. 0) then
            write(7,130)
            write(7,36) mseq,mgflag,iconsf,mtt,ngam
         else
            write(7,131)
  131       format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam,
     .      5x,5hcprec,6x,4huref,7x,3havn)
            write(7,361) mseq,mgflag,iconsf,mtt,ngam,cprec,uref,avn
  361       format(5i10,3f10.4)
         end if
      end if
c
      write(7,132)
  132 format(6x,4hissc,1x,9hepsssc(1),1x,9hepsssc(2),1x,9hepsssc(3),
     .       6x,4hissr,1x,9hepsssr(1),1x,9hepsssr(2),1x,9hepsssr(3))
      write(7,360) issc,(epsssc(n),n=1,3),issr,(epsssr(n),n=1,3)
  360 format(2(i10,3f10.4))
c
      write(7,136)
  136 format(6x,4hncyc,4x,6hmglevg,5x,5hnemgl,5x,5hnitfo)
      do 138 m=1,mseq
      write(7,36)ncyc1(m),mglev(m),nem(m),nitfo1(m)
  138 continue
c
      write(7,140)
  140 format(6x,4hmit1,6x,4hmit2,6x,4hmit3,6x,4hmit4,6x,9hmit5  ...)
      do 142 m=1,mseq
      write(7,36)(mitl(i,m),i=1,(mglev(m)+nem(m)))
  142 continue
c
c     1-1 blocking data - version 4.x
c
      if (iver.le.4) then
c
         write(7,144)
  144    format(19h 1-1 blocking data:)
         write(7,146)
  146    format(6x,4hnbli)
         write(7,36) nbli
         write(7,148)
  148    format(5x,5hnblon,3x,7hnblk(1),3x,7hnblk(2))
         if (nbli.gt.0) then
            do 150 n=1,nbli
            write(7,36)nblon(n),(nblk(ib,n),ib=1,2)
  150       continue
         end if 
c
         write(7,152)
  152    format(6x,4hista,6x,4hjsta,6x,4hksta,
     .          6x,4hiend,6x,4hjend,6x,4hkend,
     .          1x,9hisva(1,1),1x,9hisva(1,2))
         if (nbli.gt.0) then
            do 154 n=1,nbli
            write(7,36)(limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
  154       continue
         end if
c
         write(7,156)
  156    format(6x,4hista,6x,4hjsta,6x,4hksta,
     .          6x,4hiend,6x,4hjend,6x,4hkend,
     .          1x,9hisva(2,1),1x,9hisva(2,2))
         if (nbli.gt.0) then
            do 157 n=1,nbli
            write(7,36)(limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
  157       continue
         end if
c
      else
c
c     1-1 blocking data - version 5.x
c
         write(7,144)
         write(7,146)
         write(7,37) nbli
c
         write(7,1521)
 1521    format(2x,6hnumber,4x,4hgrid,4x,4hista,4x,4hjsta,4x,4hksta,
     .          4x,4hiend,4x,4hjend,4x,4hkend,
     .          3x,5hisva1,3x,5hisva2)
         if (nbli.gt.0) then
            mm = 0
            do 1541 n=1,nbli
               mm = mm + 1
               write(7,37)mm,nblk(1,n),
     .         (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
   37          format(10i8)
 1541       continue
         end if
c
         write(7,1521)
         if (nbli.gt.0) then
            mm = 0
            do 1571 n=1,nbli
               mm = mm + 1
               write(7,37)mm,nblk(2,n),
     .         (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
 1571       continue
         end if
c
      end if
c
c     patched grid flag
c
      ninter1 = 0
      if (ninter.gt.0) ninter1 = -1
      write(7,158)
  158 format(22h patch interface data:)
      write(7,160)
  160 format(4x,6hninter)
      write(7,36)ninter1
c
c     plot3d output data - just set to output surface data for
c     version 4.1hp or 5.0; no plot3d output for version 4.1
c
      write(7,*) 'plot3d output:'
      write(7,162)
  162 format(2x,4hgrid,1x,5hiptyp,2x,4hista,
     .2x,4hiend,2x,4hiinc,2x,4hjsta,2x,4hjend,2x,
     .4hjinc,2x,4hksta,2x,4hkend,2x,4hkinc)
      if (iver.eq.4 .and. ipar.gt.0 .and. 
     .   abs(nplot3d).gt.0) then
         write(7,76) 1,iplt3dtyp,0,0,0,0,0,0,0,0,0
   76    format(11i6)
      end if
      if (iver.ge.5) then
         write(7,76) 1,iplt3dtyp,0,0,0,0,0,0,0,0,0
      end if
c
c     no movie output - just output a placeholder in the input file
c
      movie = 0
      write(7,164)
  164 format(1x,5hmovie)
      if ((iver.eq.4 .and. nplot3d.ne.0) .or. iver.ge.5) then
         write(7,166) movie
      end if
c
c     print out data - just set to output surface data for
c     version 4.1hp or 5.0; no plot3d output for version 4.1
c
      write(7,*) 'print out:'
      write(7,162)
      if (iver.eq.4 .and. ipar.gt.0 .and.
     .   abs(nprint).gt.0) then
         write(7,76) 1,0,0,0,0,0,0,0,0,0,0
      end if
c
c     no contol surface data - just output a placeholder in the input file
c
      ncs = 0
c
      if (iver.gt.4) then
         write(7,*) 'control surfaces:'
         write(7,*) '  ncs'
         write(7,166)ncs
  166    format(i6)
         write(7,168)
  168    format(2x,4hgrid,2x,4hista,2x,4hiend,2x,4hjsta,2x,
     .          4hjend,2x,4hksta,2x,4hkend,1x,5hiwall,1x,5hinorm)
      end if
c
c        dynamic mesh data
c
c        count number of translating/rotating zones in the split grid
c
         ntrans = 0
         nrotat = 0
         do igd=1,ngrid
            itrans = xmap( 8,1,igd)
            irotat = xmap(16,1,igd)
            if (itrans.gt.0) then
                ntrans = ntrans + 1
            end if
            if (irotat.gt.0) then
                nrotat = nrotat + 1
            end if
         end do
c
         write(7,*) 'moving grid data - translation'
         write(7,*) 'ntrans'
         write(7,'(i7)') ntrans
         write(7,*) '  lref'
         if (ntrans.gt.0) write(7,'(f9.3)') tlref
         write(7,*) '  grid   itrans    rfreq     xmag',
     .   '     ymag     zmag'
         if (ntrans.gt.0) then
            do igd=1,ngrid
               itrans = xmap( 8,1,igd)
               rfreqt = xmap( 9,1,igd)
               xmag   = xmap(10,1,igd) 
               ymag   = xmap(11,1,igd)
               zmag   = xmap(12,1,igd)
               if (itrans.gt.0) then
                  write(7,'(i7,i9,4f9.3)') igd,itrans,rfreqt,xmag,ymag,
     .            zmag
               end if
            end do
         end if
         write(7,*) '  grid    dxmax    dymax    dzmax'
         if (ntrans.gt.0) then
            do igd=1,ngrid
               dxmax  = xmap(13,1,igd)
               dymax  = xmap(14,1,igd)
               dzmax  = xmap(15,1,igd)
               if (itrans.gt.0) then
                  write(7,'(i7,4f9.3)') igd,dxmax,dymax,dzmax
               end if
            end do
         end if
c
         write(7,*) 'moving grid data - rotation'
         write(7,*) 'nrotat'
         write(7,'(i7)') nrotat
         write(7,*) '     lref'
         if (nrotat.gt.0) write(7,'(f10.4)') rlref
         write(7,*) '  grid   irotat    rfreq   thxmag   thymag',
     .              '   thzmag   xorig   yorig   zorig'
         if (nrotat.gt.0) then
            do igd=1,ngrid
               irotat = xmap(16,1,igd)
               rfreqr = xmap(17,1,igd)
               thxmag = xmap(18,1,igd)
               thymag = xmap(19,1,igd)
               thzmag = xmap(20,1,igd)
               xorig  = xmap(22,1,igd)
               yorig  = xmap(23,1,igd)
               zorig  = xmap(23,1,igd)
               if (irotat.gt.0) then
                  write(7,'(i7,i9,7f9.3)') igd,irotat,rfreqr,thxmag,
     .            thymag,thzmag,xorig,yorig,zorig
               end if
            end do
         end if
         write(7,*) '  grid   thxmax   thymax   thzmax'
         if (nrotat.gt.0) then
            do igd=1,ngrid
               irotat = xmap(16,1,igd)
               thxmax = xmap(24,1,igd)
               thymax = xmap(25,1,igd)
               thzmax = xmap(26,1,igd)
               if (irotat.gt.0) then
                  write(7,'(i7,3f9.3)') igd,thxmax,thymax,thzmax
               end if
            end do
         end if
c
c        ignore dynamic mesh data for now!
c
c     sensitivity data (version 4 only)
c
      if (iver.eq.4 .and. isd.gt.0) then
         write(7,*) 'sensitivity data:'
         write(7,*) '       ndv   isdform'
         write(7,177) ndv,isdform
  177    format(2i10)
         write(7,*) 'sensitivity i/o files'
         write(7, '(a60)') sdgridout
         write(7, '(a60)') dovrlap
         write(7, '(a60)') dpatch
         write(7, '(a60)') dresid
       end if
c
c***********************************************************************
c
c     begin generation of ronnie patched-grid input file
c
c***********************************************************************
c
      if (ninter.gt.0) then
c
      write(9,18)
      write(9,'(a60)')gridout
      write(9,'(a60)')rout
      write(9,'(a60)')patch
c
      if (ironver .eq. 1) then
          itrace = -1
          write(9,'(''ioflag  itrace'')')
          write(9,'(i6,2x,i6)') ioflag,itrace
      else
          ioflag = 0
      end if
      write(9,11)(titleron(i),i=1,20)
c
      write(9,200)
  200 format(1x,5hngrid)
      write(9,201) ngridout
  201 format(13i6)
c
      write(9,202)
  202 format(3x,3hncg,3x,3hiem,2x,4hidim,2x,4hjdim,2x,4hkdim)
      do 400 n=1,ngrid
      write(9,201) ncg,iemg(n),idimg(n),jdimg(n),kdimg(n)
  400 continue
c
      write(9,206)
  206 format(1x,6hninter)
      write(9,201)ninter
      write(9,208)
  208 format(3x,3hint,1x,5hiifit,1x,5hlimit,1x,5hitmax,1x,
     .5hmcxie,1x,5hmceta,3x,3hc-0,1x,5hiorph)
      do 410 int=1,ninter
      write(9,201) int,iifit(int),llimit(int),iitmax(int),mmcxie(int),
     .mmceta(int),iic0(int),iiorph(int)
  410 continue 
c
      if (ioflag .eq. 0) then
         write(9,210)
  210    format(3x,3hint,4x,2hto,2x,4hxie1,2x,4hxie2,2x,4heta1,2x,
     .   4heta2,3x,3hnfb,4(2x,4hfrom),3x,3h...)
         in = 0
         do 412 int=1,ninter0
         if (iconcat(int).eq.0) then
            in = in + 1
            write(9,201) in,ito(int),xit1(int),xit2(int),ett1(int),
     .      ett2(int),nfb(int),(ifrom(int,l),l=1,nfb(int))
         end if
  412    continue
      else
         write(9,212)
  212    format(3x,3hint,4x,2hto,2x,4hxie1,2x,4hxie2,2x,4heta1,2x,
     .   4heta2,3x,3hnfb,/,8x,4hfrom,2x,4hxie1,2x,4hxie2,2x,4heta1,
     .   2x,4heta2)
         in = 0
         do 416 int=1,ninter0
         if (iconcat(int).eq.0) then
            in = in + 1
            write(9,201) in,ito(int),xit1(int),xit2(int),ett1(int),
     .      ett2(int),nfb(int)
            do 418 l = 1,nfb(int)
            write(9,203) ifrom(int,l),xif1(int,l),xif2(int,l),
     .      etf1(int,l),etf2(int,l)
 418        continue
         end if
 416     continue
 203     format(6x,13i6)
      end if
c
      end if
      return
      end
