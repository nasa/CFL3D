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
      subroutine global2(maxbl,maxgr,msub1,nintr,intmx,ngrid,idimg,
     .                   jdimg,kdimg,levelg,ncgg,nblg,iindx,llimit,
     .                   iitmax,mmcxie,mmceta,ncheck,iifit,iic0,
     .                   iiorph,iitoss,ifiner,dx,dy,dz,dthetx,
     .                   dthety,dthetz,myid,mptch,mxxe,icall,iunit25,
     .                   nou,bou,ibufdim,nbuf,ifrom,xif1,etf1,xif2,
     .                   etf2,igridg,iemg,nblock,ioflag,imode)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read dynamic patch input parameters
c
c     ioflag  - flag for input format
c             = 0 old style ronnie input
c             = 1 new style ronnie input
c             = 2 dynamic patching style from cfl3d input file
c
c     imode   - flag for indicating stand-alone ronnie or part of cfl3d
c             = 0 stand-alone ronnie patched-grid preprocessor
c             = 1 dynamic patching within cfl3d (also within precfl3d)
c
c***********************************************************************
c
c     maxgr   - maximum number of grids
c     maxbl   - maximum number of blocks
c     intmx   - maximum number of block interpolations
c     mptch   - maximum dim. of any block face involved in a dynamic
c               patch interface
c     mxxe    - size of 1-d array used to store interpolation 
c               coefficients for all dynamic patch interfaces,
c               including those on coarser blocks
c     msub1   - maximum number of blocks a single dynamic patch face 
c               may be interpolated from
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
      integer xi1,xi2,et1,et2,xif11,xif22,etf11,etf22
      integer ifrom(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .        etf2(msub1)
c
      real realval(20)
c
      dimension nou(nbuf),iemg(maxgr),igridg(maxbl)
      dimension iindx(intmx,6*msub1+9),llimit(intmx),
     .          iitmax(intmx),mmcxie(intmx),mmceta(intmx),
     .          ncheck(maxbl),iifit(intmx),iic0(intmx),iiorph(intmx),
     .          iitoss(intmx),ifiner(intmx)
      dimension dx(intmx,msub1),dy(intmx,msub1),dz(intmx,msub1),
     .          dthetx(intmx,msub1),dthety(intmx,msub1),
     .          dthetz(intmx,msub1)
      dimension idimg(maxbl),jdimg(maxbl),kdimg(maxbl),levelg(maxbl),
     .          ncgg(maxgr),nblg(maxgr)
c
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /conversion/ radtodeg
      common /unit5/ iunit5
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /igrdtyp/ ip3dgrd,ialph
c
   36 format(13i6)
   38 format(i6,i6,10i9)
   39 format(3x,6f9.4)
   40 format(6x,i6,10i9)
c
      if (icall .eq. 0) then
         ierrflg = -99
      else
         ierrflg = -1
      end if
c
      nintr = 0
c
      if (imode.eq.1) then
         read(iunit5,*,end=999)
      end if
      read(iunit5,*)
      read(iunit5,*)nintr
c
      if (nintr .eq. 0) return
c
      lnbuf = lnbuf + 2
c
      if (ioflag.eq.2) then
         write(iunit25,1260)
 1260    format(/,1x,24hdynamic patch input data)
         write(iunit25,1261)
 1261    format(1x,5hnintr)
      else
         write(iunit25,1262)
 1262    format(/,16hpatch input data)
         write(iunit25,1263)
 1263    format(1x,6hninter)
      end if
c
      do 680 igrid = 1,ngrid
      ncg = ncgg(igrid)
      if(igrid.eq.1) then
        ncgmax = ncg
      else
        ncgmax = max(ncgmax,ncg)
      end if
  680 continue
      ntest = (ncgmax+1)*nintr
      if (ntest.gt.intmx) then
         write(iunit25,9220)ntest
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 9220 format(51h too many blocks to interpolate; intmx should be at,
     .        7h least ,i4)
      write(iunit25,36)nintr
c
c     NOTE:  A patched grid with two-way communciation between
c            the grids on either side of the patch requires at 
c            least 2 block interpolations.
c
c      int   - patch interpolation number
c      iifit - type of fit
c            = 1 bi-linear fit
c            = 2 serendipity (degenerate) bi-quadratic fit
c            = 3 quadratic fit in xie, linear fit in eta 
c            = 4 linear fit in xie, quadratic fit in eta
c     llimit - maximum step size to find cell centers (Newton iteration)
c     iitmax - maximum number of iterations for the above
c     mmcxie - percent (0-100) of the total number of points on xie=0
c              "to" boundary which may be found to be not coincident 
c              with the "from" xie=0 boundary before the entire
c              boundary is assumed not to be coincident. 
c              If mmcxie > 100, then the boundaries are assumed to be
c              coincident regardless of the number of individual points
c              found not to be coincident, and steps are taken to 
c              render them coincident in computational space
c     mmceta - same as mmcxie, but for eta=0 boundary
c
c        ito - block number/topology of "to" surface
c              ("to" refers to the block being interpolated)
c              See note below.
c
c     nfblks - number of block boundaries which make up the "from"
c              side of the patch surface ("from" refers to the
c              the block(s) from which the interpolations are made)
c
c      ifrom - block number/topology of "from" surface(s)
c
c     NOTE:  The general form is ifrom (or ito) = Nmn,
c            where N  indicates the global block number of the block
c                  m  indicates the generalized coodinate which is 
c                     constant along the patch surface (m=1 implies
c                     i=const, m=2 implies j=const,
c                     and m=3 implies k=const)
c                  n  indicates on which of the two possible 
c                     m=constant surfaces the patch occurs on
c                     (n=1 implies the patch occurs on m=1,
c                     while n=2 implies the patch occurs on m=mmax)
c            for example
c
c            = 0111 patching occurs along the i=1 boundary  of block 1
c
c            = 0112 patching occurs along the i=imax boundary of block 1
c
c            = 0121 patching occurs along the j=1 boundary of block 1
c
c            = 0122 patching occurs along the j=jmax boundary of block 1
c
c            = 0131 patching occurs along the k=1 boundary of block 1
c
c            = 0132 patching occurs along the k=kmax boundary of block 1
c
c
c
      read(iunit5,*)
      if (ioflag.eq.2) then
         write(iunit25,1252)
 1252    format(3x,3hint,1x,5hiifit,4x,5hlimit,4x,5hitmax,4x,
     .   5hmcxie,4x,5hmceta,6x,3hc-0,4x,5hiorph,4x,5hitoss)
      else
         write(iunit25,1253)
 1253    format(3x,3hint,1x,5hiifit,4x,5hlimit,4x,5hitmax,4x,
     .   5hmcxie,4x,5hmceta,6x,3hc-0,4x,5hiorph)
      end if
      if (nintr.gt.0) then
         do 7030 n=1,nintr
         iifit(n) = -99
7030     continue
         do 7031 n=1,nintr
         if (ioflag.eq.2) then
            read(iunit5,*) int,iif,llimit(int),iitmax(int),mmcxie(int),
     .      mmceta(int),iic0(int),iiorph(int),iitoss(int)
         else
            read(iunit5,*) int,iif,llimit(int),iitmax(int),mmcxie(int),
     .      mmceta(int),iic0(int),iiorph(int)
            iitoss(int) = 0
         end if
      if (iifit(int) .ne. -99) then
         write(11,'('' program terminated in dynamic patching '',
     .         ''routines - see file '',a60)') grdmov
         write(iunit25,*) ' stopping...attempting to set data for int ',
     .   int,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      iifit(int) = iif
         if (iifit(int).eq.0) iifit(int) = 1
         if (iic0(int).gt.0) then
            if (iifit(int) .ge. 0) then
               iifit(int)  = 1
            else 
               iifit(int)  = -1
            end if
            mmcxie(int) = 0
            mmceta(int) = 0
         end if
         if (ioflag.eq.2) then
            write(iunit25,38) int,iifit(int),llimit(int),iitmax(int),
     .      mmcxie(int),mmceta(int),iic0(int),iiorph(int),iitoss(int)
         else
            write(iunit25,38) int,iifit(int),llimit(int),iitmax(int),
     .      mmcxie(int),mmceta(int),iic0(int),iiorph(int)
         end if
 7031    continue
      end if 
c
      nfbmax = 0
      do 1702 n=1,nintr
      iindx(n,1) = -99 
1702  continue
c
 8146 format(3x,3hint,4x,2hto,2x,4hxie1,2x,4hxie2,2x,4heta1,
     .2x,4heta2,3x,3hnfb,4(2x,4hfrom),3x,3h...)
 8147 format(3x,3hint,4x,2hto,5x,4hxie1,5x,4hxie2,5x,4heta1,5x,4heta2,
     .6x,3hnfb)
 8148 format(10x,2hdx,7x,2hdy,7x,2hdz,3x,6hdthetx,3x,6hdthety,
     .3x,6hdthetz)
 8149 format(8x,4hfrom,5x,4hxie1,5x,4hxie2,5x,4heta1,5x,4heta2)
c
      if (ioflag.eq.0) then
          read(iunit5,*)
          write(iunit25,8146)
      else if (ioflag.eq.1) then
          read(iunit5,*)
          read(iunit5,*)
          write(iunit25,8147)
          write(iunit25,8149)
      end if
c 
      do 1703 n=1,nintr
      if (ioflag.eq.2) then
         read(iunit5,*)
         read(iunit5,*) int,ito,xi1,xi2,et1,et2,nfb
      else if (ioflag.eq.1) then
         read(iunit5,*) int,ito,xi1,xi2,et1,et2,nfb
      else if (ioflag.eq.0) then
         read(iunit5,*) int,ito,xi1,xi2,et1,et2,nfb,(ifrom(l),l=1,nfb)
         do l=1,nfb
            xif1(l)       = 0
            xif2(l)       = 0
            etf1(l)       = 0
            etf2(l)       = 0
            dx(int,l)     = 0
            dy(int,l)     = 0
            dz(int,l)     = 0
            dthetx(int,l) = 0
            dthety(int,l) = 0
            dthetz(int,l) = 0
         end do
      end if
c     make sure parameter msub1 is big enough
      if (nfb.gt.msub1) then
         write(iunit25,9019) nfb
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      nfbmax = max(nfbmax,nfb)
 9019 format(1x,37hmsub1 is too small; should be least: ,i2)
      if (ioflag.eq.1) then
         do l=1,nfb
            read(iunit5,*) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
            dx(int,l)     = 0
            dy(int,l)     = 0
            dz(int,l)     = 0
            dthetx(int,l) = 0
            dthety(int,l) = 0
            dthetz(int,l) = 0
         end do
      else if (ioflag.eq.2) then
         do l=1,nfb
            read(iunit5,*)
            read(iunit5,*) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
            read(iunit5,*)
            read(iunit5,*) (realval(i),i=1,6)
            dx(int,l)     = realval(1)
            dthetx(int,l) = realval(4)
            if(ialph > 0) then
              dy(int,l)     = -realval(3)
              dthety(int,l) = -realval(6)
              dz(int,l)     = realval(2)
              dthetz(int,l) = realval(5)
            else
              dy(int,l)     = realval(2)
              dthety(int,l) = realval(5)
              dz(int,l)     = realval(3)
              dthetz(int,l) = realval(6)
            end if
         end do
      end if
c
      if (iindx(int,1) .ne. -99) then
         write(11,'('' program terminated in dynamic patching '',
     .         ''routines - see file '',a60)') grdmov
        write(iunit25,*)
     .  ' stopping...attempting to set data for int ',int,
     .  ' more than once'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      iindx(int,1)       = nfb
      do 1701 l=1,nfb
      iindx(int,l+1)     = ifrom(l)/100
      iindx(int,l+nfb+2) = ifrom(l)-iindx(int,l+1)*100
 1701 continue
c
      iindx(int,nfb+2)   = ito/100
      nbl                 = nblg(iindx(int,nfb+2))
      iindx(int,2*nfb+3) = ito-iindx(int,nfb+2)*100
c
c     set full ranges if zeros are input
      if (iindx(int,2*nfb+3)/10.eq.1) then
         if (xi1.eq.0 .and. xi2.eq.0) then
            xi1 = 1
            xi2 = jdimg(nbl)
         end if
         if (et1.eq.0 .and. et2.eq.0) then
            et1 = 1
            et2 = kdimg(nbl)
         end if
      end if
      if (iindx(int,2*nfb+3)/10.eq.2) then
         if (xi1.eq.0 .and. xi2.eq.0) then
            xi1 = 1
            xi2 = kdimg(nbl)
         end if
         if (et1.eq.0 .and. et2.eq.0) then
            et1 = 1
            et2 = idimg(nbl)
         end if
      end if
      if (iindx(int,2*nfb+3)/10.eq.3) then
         if (xi1.eq.0 .and. xi2.eq.0) then
            xi1 = 1
            xi2 = jdimg(nbl)
         end if
         if (et1.eq.0 .and. et2.eq.0) then
            et1 = 1
            et2 = idimg(nbl)
         end if
      end if
c
c     set full ranges if zeros are input for search range
c     in "from" block(s)
c
      do 1705 l=1,nfb
      mbl = nblg(iindx(int,l+1))
      if (iindx(int,l+nfb+2)/10.eq.1) then
         if (xif1(l).eq.0 .and. xif2(l).eq.0) then
            xif1(l) = 1
            xif2(l) = jdimg(mbl)
         end if
         if (etf1(l).eq.0 .and. etf2(l).eq.0) then
            etf1(l) = 1
            etf2(l) = kdimg(mbl)
         end if
      end if
      if (iindx(int,l+nfb+2)/10.eq.2) then
         if (xif1(l).eq.0 .and. xif2(l).eq.0) then
            xif1(l) = 1
            xif2(l) = kdimg(mbl)
         end if
         if (etf1(l).eq.0 .and. etf2(l).eq.0) then
            etf1(l) = 1
            etf2(l) = idimg(mbl)
         end if
      end if
      if (iindx(int,l+nfb+2)/10.eq.3) then
         if (xif1(l).eq.0 .and. xif2(l).eq.0) then
            xif1(l) = 1
            xif2(l) = jdimg(mbl)
         end if
         if (etf1(l).eq.0 .and. etf2(l).eq.0) then
            etf1(l) = 1
            etf2(l) = idimg(mbl)
         end if
      end if
 1705 continue
c
      if (ioflag.eq.0 .or. ioflag.eq.1) then
         write(iunit25,38) int,ito,xi1,xi2,et1,et2,nfb
         do l=1,nfb
            write(iunit25,40) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
         end do
      else
         write(iunit25,8147)
         write(iunit25,38) int,ito,xi1,xi2,et1,et2,nfb
         do l=1,nfb
            write(iunit25,8149)
            write(iunit25,40) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
            write(iunit25,8148)
            write(iunit25,39) real(dx(int,l)),real(dy(int,l)),
     .                        real(dz(int,l)),real(dthetx(int,l)),
     .                        real(dthety(int,l)),real(dthetz(int,l))
            dthetx(int,l) = dthetx(int,l)/radtodeg
            dthety(int,l) = dthety(int,l)/radtodeg
            dthetz(int,l) = dthetz(int,l)/radtodeg
         end do
      end if
c
c     check input ranges of xie and eta
c
      if (iindx(int,2*nfb+3)/10.eq.1) then
         if (xi1.lt.1 .or. xi1.gt.jdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if 
         if (xi2.lt.1 .or. xi2.gt.jdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xi2 .le. xi1) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if 
         if (et1.lt.1 .or. et1.gt.kdimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (et2.lt.1 .or. et2.gt.kdimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (et2 .le. et1) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this i=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (iindx(int,2*nfb+3)/10.eq.2) then
c
         if (xi1.lt.1 .or. xi1.gt.kdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if   
         if (xi2.lt.1 .or. xi2.gt.kdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xi2 .le. xi1) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (et1.lt.1 .or. et1.gt.idimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (et2.lt.1 .or. et2.gt.idimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (et2 .le. et1) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this j=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (iindx(int,2*nfb+3)/10.eq.3) then
c
         if (xi1.lt.1 .or. xi1.gt.jdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if   
         if (xi2.lt.1 .or. xi2.gt.jdimg(nbl)) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xi2 .le. xi1) then
            write(iunit25,*)' incorrect xie range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (et1.lt.1 .or. et1.gt.idimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (et2.lt.1 .or. et2.gt.idimg(nbl)) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (et2 .le. et1) then
            write(iunit25,*)' incorrect eta range input on ',
     .      'this k=constant "to" face'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
c     check input search ranges of xie and eta - "from" block(s)
c
      do 1716 l=1,nfb
      mbl = nblg(iindx(int,l+1))
      if (iindx(int,l+nfb+2)/10.eq.1) then
         if (xif1(l).lt.1 .or. xif1(l).gt.jdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l).lt.1 .or. xif2(l).gt.jdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l) .le. xif1(l)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf1(l).lt.1 .or. etf1(l).gt.kdimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf2(l).lt.1 .or. etf2(l).gt.kdimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (etf2(l) .le. etf1(l)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the i=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (iindx(int,l+nfb+2)/10.eq.2) then
c
         if (xif1(l).lt.1 .or. xif1(l).gt.kdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l).lt.1 .or. xif2(l).gt.kdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l) .le. xif1(l)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf1(l).lt.1 .or. etf1(l).gt.idimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf2(l).lt.1 .or. etf2(l).gt.idimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (etf2(l) .le. etf1(l)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the j=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (iindx(int,l+nfb+2)/10.eq.3) then
c
         if (xif1(l).lt.1 .or. xif1(l).gt.jdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l).lt.1 .or. xif2(l).gt.jdimg(mbl)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (xif2(l) .le. xif1(l)) then
            write(iunit25,*)' incorrect xie search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf1(l).lt.1 .or. etf1(l).gt.idimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (etf2(l).lt.1 .or. etf2(l).gt.idimg(mbl)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
          if (etf2(l) .le. etf1(l)) then
            write(iunit25,*)' incorrect eta search range input on ',
     .      'the k=constant "from" face',l
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
1716  continue
c
      iindx(int,2*nfb+6) = xi1
      iindx(int,2*nfb+7) = xi2
      iindx(int,2*nfb+8) = et1
      iindx(int,2*nfb+9) = et2
      do 1707 l=1,nfb
      iindx(int,2*nfb+9+l) = xif1(l)
      iindx(int,3*nfb+9+l) = xif2(l)
      iindx(int,4*nfb+9+l) = etf1(l)
      iindx(int,5*nfb+9+l) = etf2(l)
1707  continue
c
      iindx(int,2*nfb+4) = (xi2-xi1)*(et2-et1)
      if (int .eq. 1) then
         iindx(int,2*nfb+5) = 1
      else
         nfb1 = iindx(int-1,1)
         iindx(int,2*nfb+5) = iindx(int-1,2*nfb1+5)
     .                      + iindx(int-1,2*nfb1+4)
      end if
 1703 continue
c
      if (imode .eq. 0) then
         write(iunit25,1919)
 1919    format(/1x,16hSUMMARY BY GRIDS)
         write(iunit25,117)
  117    format(3x,4hgrid,2x,5hlevel,2x,5hblock,
     .          3x,4hjdim,3x,4hkdim,3x,4hidim,3x,9hgrid pts.)
c
         lfem   = 0
         lfgm   = ncgmax+1
         lcgm   = 1
         igptot = 0
         do 105 igrid=1,ngrid
         iem         = iemg(igrid)
         nbl         = nblg(igrid)
         ncg         = ncgg(igrid)
         igridg(nbl) = igrid
         levelg(nbl) = iem+lfgm
         if (iem.eq.0) lfgm = max(lfgm,levelg(nbl))
         if (iem.gt.0) lfem = max(lfem,levelg(nbl))
         igpts       = jdimg(nbl)*kdimg(nbl)*idimg(nbl)
         igptot      = igptot+igpts
         write(iunit25,7)igrid,levelg(nbl),nbl,
     .                   jdimg(nbl),kdimg(nbl),idimg(nbl),igpts
         if (ncg.gt.0 .and. iem.eq.0) then
         do 104 n=1,ncg
         nbl         = nbl+1
         igridg(nbl) = igrid
         levelg(nbl) = levelg(nbl-1)-1
         igpts       = jdimg(nbl)*kdimg(nbl)*idimg(nbl)
         igptot      = igptot+igpts
         write(iunit25,7)igrid,levelg(nbl),nbl,
     .                   jdimg(nbl),kdimg(nbl),idimg(nbl),igpts
  104    continue
         end if
  105    continue
    7    format(6i7,i12)
         write(iunit25,2021)igptot
 2021    format(/,37x,5hTOTAL,i12)
c
         write(iunit25,2304)
 2304    format(/1x,17hSUMMARY BY LEVELS)
         write(iunit25,110)
  110    format(2x,5hlevel,3x,4hgrid,2x,5hblock)
         lf = lfem
         if (lfem.eq.0) lf = lfgm
         do 115 levelc=lf,lcgm,-1
         do 114 nbl=1,nblock
         if (levelc.ne.levelg(nbl)) go to 114
         igrid = igridg(nbl)
         write(iunit25,7)levelg(nbl),igrid,nbl
  114    continue
  115    continue
c
         if (mseq.gt.(ncgmax+1)) then
         write(iunit25,153)mseq,ncgmax
  153    format(1x,12hmseq,ncgmax=,2i4)
         write(iunit25,155)mseq,ncgmax
  155    format(1x,28herror in input, mseq, ncgmax,2i5)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         write(iunit25,4008)
 4008    format(/)
         if (lfem.ne.0) write(iunit25,4009)lfem
 4009    format(31h level of finest embedded mesh=,i3)
         write(iunit25,4019)lcgm
 4019    format(31h level of coarsest global mesh=,i3)
         write(iunit25,4029)lfgm
 4029    format(31h level of finest global mesh  =,i3)
c
         write(iunit25,7227)
 7227    format(/1x,25hSUMMARY OF GRID SEQUENCES)
         write(iunit25,9)
    9    format(1x,8hsequence,2x,14hstarting level,4x,
     .          12hending level)
         do 19 m=1,mseq
         write(iunit25,29)m,levelt(m),levelb(m)
   19    continue
   29    format(2x,i7,2i16)
c
      end if
c
c  generate additional patch data for coarser blocks
c
      ntemp  = nintr
      ifirst = 1
      do 8400 n=1,nintr
      nfb    = iindx(n,1)
      ifiner(n) = 0
      if (ncgg(iindx(n,nfb+2)).gt.0) then
         if (ifirst.eq.1) then
            if (ioflag.eq.2) then
               write(iunit25,8601)
            else
               write(iunit25,8602)
            end if
            ifirst = 0
         end if
         ito = (iindx(n,nfb+2)*100) + iindx(n,2*nfb+3)
         itonew = nblg(iindx(n,nfb+2))*100 + iindx(n,2*nfb+3)
         xi1 = iindx(n,2*nfb+6)
         xi2 = iindx(n,2*nfb+7)
         et1 = iindx(n,2*nfb+8)
         et2 = iindx(n,2*nfb+9)
         do 9101 i=1,nfb
         ifrom(i) = nblg(iindx(n,i+1))*100 + iindx(n,i+nfb+2)
 9101    continue
         write(iunit25,8147)
         write(iunit25,38) n,itonew,xi1,xi2,et1,et2,nfb
         do 9102 l=1,nfb
         write(iunit25,8149)
         write(iunit25,40) ifrom(l),xif1(l),xif2(l),etf1(l),etf2(l)
         if (ioflag.eq.2) then
            write(iunit25,8148)
            write(iunit25,39) real(dx(n,l)),real(dy(n,l)),
     .                        real(dz(n,l)),
     .                        real(dthetx(n,l))*real(radtodeg),
     .                        real(dthety(n,l))*real(radtodeg),
     .                        real(dthetz(n,l))*real(radtodeg)
         end if
 9102    continue
c
         do 8200 m=1,ncgg(iindx(n,nfb+2))
         ntemp = ntemp + 1
         lprev                = ntemp - 1
         if (m.eq.1) then
            ifiner(ntemp)      = n
         else
            ifiner(ntemp)      = ntemp - 1
         end if
         iifit(ntemp)         = iifit(n)
         iic0(ntemp)          = iic0(n)
         iiorph(ntemp)        = iiorph(n)
         iitoss(ntemp)        = iitoss(n)
         llimit(ntemp)        = llimit(n)
         iitmax(ntemp)        = iitmax(n)
         mmcxie(ntemp)        = mmcxie(n)
         mmceta(ntemp)        = mmceta(n)
         do 8201 l=1,nfb
         dx(ntemp,l)          = dx(n,l) 
         dy(ntemp,l)          = dy(n,l)
         dz(ntemp,l)          = dz(n,l)
         dthetx(ntemp,l)      = dthetx(n,l)
         dthety(ntemp,l)      = dthety(n,l)
         dthetz(ntemp,l)      = dthetz(n,l)
 8201    continue
         iindx(ntemp,1)       = nfb
         iindx(ntemp,nfb+2)   = nblg(iindx(n,nfb+2)) + m
         iindx(ntemp,nfb+3)   = iindx(n,nfb+3)
         iindx(ntemp,2*nfb+3) = iindx(n,2*nfb+3)
         nfb1                  = iindx(lprev,1)
         nbl                   = iindx(ntemp,nfb+2)      
c
         xi1 = iindx(n,2*nfb+6)
         xi2 = iindx(n,2*nfb+7)
         et1 = iindx(n,2*nfb+8)
         et2 = iindx(n,2*nfb+9)
c
         iindx(ntemp,2*nfb+6) = xi1/(2*m) + 1
         iindx(ntemp,2*nfb+7) = xi2/(2*m) + 1
         iindx(ntemp,2*nfb+8) = et1/(2*m) + 1
         iindx(ntemp,2*nfb+9) = et2/(2*m) + 1
c
c        2-d meshes
         if (iindx(n,2*nfb+6).le. 2)
     .      iindx(ntemp,2*nfb+6) = iindx(n,2*nfb+6)
         if (iindx(n,2*nfb+7).le. 2)
     .      iindx(ntemp,2*nfb+7) = iindx(n,2*nfb+7)
         if (iindx(n,2*nfb+8).le. 2) 
     .      iindx(ntemp,2*nfb+8) = iindx(n,2*nfb+8)
         if (iindx(n,2*nfb+9).le. 2)
     .      iindx(ntemp,2*nfb+9) = iindx(n,2*nfb+9)
c
         xi1 = iindx(ntemp,2*nfb+6)
         xi2 = iindx(ntemp,2*nfb+7)
         et1 = iindx(ntemp,2*nfb+8)
         et2 = iindx(ntemp,2*nfb+9)
c
         do 8202 l=1,nfb
         xif11 = iindx(n,2*nfb+9+l)
         xif22 = iindx(n,3*nfb+9+l)
         etf11 = iindx(n,4*nfb+9+l)
         etf22 = iindx(n,5*nfb+9+l)
c
         iindx(ntemp,2*nfb+9+l) = xif11/(2*m) + 1
         iindx(ntemp,3*nfb+9+l) = xif22/(2*m) + 1
         iindx(ntemp,4*nfb+9+l) = etf11/(2*m) + 1
         iindx(ntemp,5*nfb+9+l) = etf22/(2*m) + 1
c
c        2-d meshes
         if (iindx(n,2*nfb+9+l).le. 2)
     .      iindx(ntemp,2*nfb+9+l) = iindx(n,2*nfb+9+l)
         if (iindx(n,3*nfb+9+l).le. 2)
     .      iindx(ntemp,3*nfb+9+l) = iindx(n,3*nfb+9+l)
         if (iindx(n,4*nfb+9+l).le. 2)
     .      iindx(ntemp,4*nfb+9+l) = iindx(n,4*nfb+9+l)
         if (iindx(n,5*nfb+9+l).le. 2)
     .      iindx(ntemp,5*nfb+9+l) = iindx(n,5*nfb+9+l)
8202     continue
c
         iindx(ntemp,2*nfb+4) = (xi2-xi1)*(et2-et1) 
         iindx(ntemp,2*nfb+5) = iindx(lprev,2*nfb1+4) +
     .                           iindx(lprev,2*nfb1+5) 
         mxxchk = iindx(ntemp,2*nfb+4)+iindx(ntemp,2*nfb+5)-1
c
         if (icall .gt. 0) then
            if (mxxchk.gt.mxxe) then
               write(11,'('' program terminated in dynamic patching '',
     .         ''routines - see file '',a60)') grdmov
               write(iunit25,*)'  stopping...parameter mxxe too small'
               write(iunit25,*)'  must make  mxxe at least ',mxxchk
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         do 8100 i=1,nfb
         iindx(ntemp,i+1)     = nblg(iindx(n,i+1)) + m
         iindx(ntemp,nfb+i+2) = iindx(n,nfb+i+2)
         if (ncgg(iindx(n,i+1)) .lt. m) 
     .       iindx(ntemp,i+1) = nblg(iindx(n,i+1)) +
     .                           ncgg(iindx(n,i+1))
             ifrom(i) = (iindx(ntemp,i+1)*100) + iindx(ntemp,i+nfb+2)
 8100    continue
         ito = (iindx(ntemp,nfb+2)*100) + iindx(ntemp,2*nfb+3)
         write(iunit25,8147)
         write(iunit25,37) ntemp,ito,xi1,xi2,et1,et2,nfb
   37    format(2h *,i4,i6,10i9)
         do 9302 l=1,nfb
         write(iunit25,8149)
         write(iunit25,40) ifrom(l),
     .               iindx(ntemp,2*nfb+9+l),iindx(ntemp,3*nfb+9+l),
     .               iindx(ntemp,4*nfb+9+l),iindx(ntemp,5*nfb+9+l)
         if (ioflag.eq.2) then
            write(iunit25,8148)
            write(iunit25,39) real(dx(n,l)),real(dy(n,l)),
     .                        real(dz(n,l)),
     .                        real(dthetx(n,l))*real(radtodeg),
     .                        real(dthety(n,l))*real(radtodeg),
     .                        real(dthetz(n,l))*real(radtodeg)
         end if
 9302    continue
 8200    continue
      end if
c
c     renumber "to" and "from" block numbers on the finest global levels
c     and embeded levels to reflect CFL3D block numbering system
c
      iindx(n,nfb+2) = nblg(iindx(n,nfb+2))
      do 8300 i=1,nfb
      iindx(n,i+1)   = nblg(iindx(n,i+1))
 8300 continue
c
 8400 continue
 8601   format(/,1x,47hREARRANGED BLOCK NUMBERS AND DYNAMIC PATCH DATA,
     .              29h GENERATED FOR COARSER MESHES,
     .         /,11x,47h    block references to CFL3D internal ordering,
     .         /,11x,43h         * ---> coarser level interpolation,/)
 8602   format(/,1x,39hREARRANGED BLOCK NUMBERS AND PATCH DATA,
     .              29h GENERATED FOR COARSER MESHES,
     .         /,11x,47h    block references to CFL3D internal ordering,
     .         /,11x,43h         * ---> coarser level interpolation,/)
      nintr = ntemp
c
c     check parameter sizes, multigridability 
c
      mpa = 0
      do 1713 n=1,nintr
      if (n.gt.1) nfb1 = iindx(n-1,1)
      nfb = iindx(n,1)
      nbl = iindx(n,nfb+2)
      xi1 = iindx(n,2*nfb+6)
      xi2 = iindx(n,2*nfb+7)
      et1 = iindx(n,2*nfb+8)
      et2 = iindx(n,2*nfb+9)
      id2  = idimg(nbl)
      jd2  = jdimg(nbl)
      kd2  = kdimg(nbl)
      itop = iindx(n,2*nfb+3)/10
      if (itop .eq. 1) then
         mpa = max(jd2,kd2,mpa)
      end if
      if (itop .eq. 2) then
         mpa = max(id2,kd2,mpa)
      end if
      if (itop .eq. 3) then
         mpa = max(id2,jd2,mpa)
      end if
      iindx(n,2*nfb+4) = (xi2-xi1)*(et2-et1) 
      do 1714 l=1,nfb
      mbl  = iindx(n,l+1)
      mbl  = iindx(n,l+1)
      id1  = idimg(mbl)
      jd1  = jdimg(mbl)
      kd1  = kdimg(mbl)
      itop = (iindx(n,l+nfb+2)-iindx(n,l+nfb+2)/100*100)/10
      if (itop .eq. 1) then
         mpa = max(jd1,kd1,mpa)
      end if
      if (itop .eq. 2) then
         mpa = max(id1,kd1,mpa)
      end if
      if (itop .eq. 3) then
         mpa = max(id1,jd1,mpa)
      end if
c
      lmptch = max(lmptch,mpa)
c
 1714 continue
c
      if (icall .gt. 0) then
         if    (mpa.gt.mptch) then
            write(iunit25,1929)mpa
 1929       format(1x,31hmptch is too small; should be: ,i4)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
      if (n.eq.1) then
         iindx(n,2*nfb+5)  = 1
      else
         iindx(n,2*nfb+5)  = iindx(n-1,2*nfb1+4)+iindx(n-1,2*nfb1+5)    
      end if
      mxxchk = iindx(nintr,2*nfb+4)+iindx(nintr,2*nfb+5)-1
 1713 continue
c
      if (icall .gt. 0) then
         if (mxxchk.gt.mxxe) then
            write(11,'('' program terminated in dynamic patching '',
     .      ''routines - see file '',a60)') grdmov
            write(iunit25,*)'  stopping...parameter mxxe too small'
            write(iunit25,*)'  must make mxxe at least ',mxxchk
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
c
c     check multigridability - only check "to" side...
c     "from" side need not always need to be multigridable
c
      istop = 0
      do n=1,nintr
         nfb = iindx(n,1)
         nbl = iindx(n,nfb+2)
         xi1 = iindx(n,2*nfb+6)
         xi2 = iindx(n,2*nfb+7)
         et1 = iindx(n,2*nfb+8)
         et2 = iindx(n,2*nfb+9)
         if (levelg(nbl) .gt. 1) then
            if (float(xi1/2) .eq. float(xi1)/2.) then
               write(iunit25,
     .         '('' begining "to" xie index not multigridable'',
     .         '' for interpolation '',i4)') n
              istop=1
            end if
            if (float(xi2/2) .eq. float(xi2)/2.) then
               write(iunit25,'('' ending "to" xie index not'',
     .         '' multigridable for interpolation '',i4)') n
              istop=1
            end if
c           skip eta multigrid checks if 2d (but check to make sure
c           et1=1 and et2=2 if 2d)
            if (idimg(nbl).eq.2) then
               if (et1.ne.1 .or. et2.ne.2) then
                  write(iunit25,'('' for 2D cases must have eta1 = 1'',
     .            '' and eta2 = 2'')')
                  istop=1
               end if
            else
               if (float(et1/2) .eq. float(et1)/2.) then
                  write(iunit25,'('' begining "to" eta index not'',
     .            '' multigridable for interpolation '',i4)') n
                  istop=1
               end if
               if (float(et2/2) .eq. float(et2)/2.) then
                  write(iunit25,'('' ending "to" eta index not'',
     .            '' multigridable for interpolation '',i4)') n
                 istop=1
               end if
            end if
         end if
      end do
      if (istop.gt.0) then
         write(iunit25,'('' stopping'')')
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c 
      ichk = 0
      if(ichk.gt.0) then
c
c     check iindx array
      write(8,*)'nintr: ',nintr
      write(8,*)
      do 6969 n=1,nintr
      write(8,*)
      write(8,*)
      write(8,*)' interp. no.: ',n
      write(8,*)
      write(8,*)' nfb: ',iindx(n,1)
      nfb = iindx(n,1)
      write(8,*)' from block(s): ',(iindx(n,l+1),l=1,nfb)
      write(8,*)'  to: ',iindx(n,nfb+2)
      write(8,*)' topology (from): ',(iindx(n,l+nfb+2),l=1,nfb)
      write(8,*)' topology (  to): ',iindx(n,2*nfb+3)
      write(8,*)' number of points: ',iindx(n,2*nfb+4)
      write(8,*)' starting index: ',iindx(n,2*nfb+5)
      write(8,*)' xie range: ',iindx(n,2*nfb+6),iindx(n,2*nfb+7)
      write(8,*)' eta range: ',iindx(n,2*nfb+8),iindx(n,2*nfb+9)
      nfb = iindx(n,1)
      do 6968 l=1,nfb
      write(8,*)' xie search range in from block ',l,' :',
     .            iindx(n,2*nfb+9+l),iindx(n,3*nfb+9+l)
      write(8,*)' eta search range, from block ',l,' :',
     .            iindx(n,4*nfb+9+l),iindx(n,5*nfb+9+l)
 6968 continue
 6969 continue
c
      end if
c
c     adjust patch parameters so dynamic patch data can be appended to
c     regular patch data
c
      if (icall .eq. 0) then
         lmsub1 = nfbmax
         lnsub1 = max(lnsub1,lmsub1)
         lintmx = nintr
         if (lintmax .gt. 1) then
            lintmax = lintmax + lintmx
         else
            lintmax = lintmx
         end if
         lmxxe = 0
         do int=1,nintr
            nfbl   = iindx(int,1)
            npts   = iindx(int,2*nfbl+4)
            nst    = iindx(int,2*nfbl+5)
            ntot   = nst + npts - 1
         end do
         lmxxe  = ntot
         if (lmaxxe .gt. 1) then
            lmaxxe  = lmaxxe + lmxxe
         else
            lmaxxe  = lmxxe
         end if
      end if
c
  999 continue
c
      return
      end
