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
      subroutine dynptch(lw,lw2,work,mgwk,wk,nwork,ncall,
     .                   maxgr,maxbl,msub1,intmx,mxxe,mptch,jdimg,
     .                   kdimg,idimg,xorig,yorig,zorig,nblock,ngrid,
     .                   levelg,ncgg,nblg,windex,ninter,iindex,nblkpt,
     .                   windx,nintr,iindx,mblkpt,llimit,iitmax,
     .                   mmcxie,mmceta,ncheck,iifit,iic0,
     .                   iiorph,iitoss,ifiner,factjlo,factjhi,
     .                   factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .                   dthetz,dthetxx,dthetyy,dthetzz,
     .                   isav_dpat,isav_dpat_b,intmax,maxxe,nsub1,
     .                   lw_temp,ireq_ar,lout,ifrom,xif1,xif2,etf1,
     .                   etf2,jjmax1,kkmax1,iiint1,iiint2,nblk1,
     .                   nblk2,jimage,kimage,jte,kte,jmm,kmm,
     .                   xte,yte,zte,xmi,ymi,zmi,xmie,ymie,
     .                   zmie,sxie,seta,sxie2,seta2,xie2s,
     .                   eta2s,temp,x2,y2,z2,x1,y1,z1,
     .                   myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .                   ibufdim,igridg,iemg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Establish zone-to-zone communication for block inter-
c     faces that move relative to one another, using a patched-grid 
c     technique (non-conservative)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      integer xif1(msub1),xif2(msub1),etf1(msub1),etf2(msub1),
     .        ifrom(msub1)
c
      dimension work(mgwk),wk(nwork)
      dimension lw(65,maxbl),lw2(43,maxbl),lw_temp(65,maxbl)
      dimension jjmax1(msub1),kkmax1(msub1),iiint1(msub1),iiint2(msub1)
      dimension jimage(msub1,mptch+2,mptch+2),
     .          kimage(msub1,mptch+2,mptch+2)
      dimension jte(msub1),kte(msub1),lout(msub1)
      dimension jmm(mptch+2),kmm(mptch+2)
      dimension nblk1(mptch+2),nblk2(mptch+2)
      dimension xte(mptch+2,mptch+2,msub1),
     .          yte(mptch+2,mptch+2,msub1),
     .          zte(mptch+2,mptch+2,msub1)
      dimension xmi(mptch+2,mptch+2,msub1),
     .          ymi(mptch+2,mptch+2,msub1),
     .          zmi(mptch+2,mptch+2,msub1)
      dimension xmie(mptch+2,mptch+2,msub1),
     .          ymie(mptch+2,mptch+2,msub1),
     .          zmie(mptch+2,mptch+2,msub1)
      dimension sxie(mptch+2,mptch+2,msub1),
     .          seta(mptch+2,mptch+2,msub1),
     .          sxie2(mptch+2,mptch+2),
     .          seta2(mptch+2,mptch+2)
      dimension xie2s(mptch+2,mptch+2),eta2s(mptch+2,mptch+2)
      dimension temp((mptch+2)*(mptch+2))
      dimension x2(mptch+2,mptch+2),y2(mptch+2,mptch+2),
     .          z2(mptch+2,mptch+2)
      dimension x1(mptch+2,mptch+2),y1(mptch+2,mptch+2),
     .          z1(mptch+2,mptch+2)
      dimension ireq_ar(intmx*3),mblk2nd(maxbl)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension windx(mxxe,2),iindx(intmx,6*msub1+9),
     .          windex(maxxe,2),iindex(intmax,6*nsub1+9), 
     .          llimit(intmx),iitmax(intmx),mmcxie(intmx),
     .          mmceta(intmx),ncheck(maxbl),iifit(intmx),
     .          mblkpt(mxxe),iic0(intmx),iiorph(intmx),iitoss(intmx),
     .          ifiner(intmx),nblkpt(maxxe)
      dimension factjlo(intmx,msub1),factjhi(intmx,msub1),
     .          factklo(intmx,msub1),factkhi(intmx,msub1)
      dimension dx(intmx,msub1),dy(intmx,msub1),dz(intmx,msub1),
     .          dthetx(intmx,msub1),dthety(intmx,msub1),
     .          dthetz(intmx,msub1)
      dimension dthetxx(intmax,nsub1),dthetyy(intmax,nsub1),
     .          dthetzz(intmax,nsub1)
      dimension isav_dpat(intmx,17),isav_dpat_b(intmx,msub1,6)
      dimension ncgg(maxgr),nblg(maxgr),levelg(maxbl),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl)
      dimension iemg(maxgr),igridg(maxbl)
c
      common /savnint/ninter0
      common /sklton/ isklton
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /is_dpatch/ maxdcnt
      common /tracer/ itrace
c
c     read in the dynamic patch data at the bottom of the input file
c
      if (myid .eq. myhost) then
         if (ncall .eq. 1) then
            icallgl = 1
            iunit   = 25
            ioflag  = 2
            imode   = 1
            call global2(maxbl,maxgr,msub1,nintr,intmx,ngrid,idimg,
     .                   jdimg,kdimg,levelg,ncgg,nblg,iindx,llimit,
     .                   iitmax,mmcxie,mmceta,ncheck,iifit,iic0,
     .                   iiorph,iitoss,ifiner,dx,dy,dz,dthetx,
     .                   dthety,dthetz,myid,mptch,mxxe,icallgl,iunit,
     .                   nou,bou,ibufdim,nbuf,ifrom,xif1,etf1,xif2,
     .                   etf2,igridg,iemg,nblock,ioflag,imode)
         end if
      end if
c
c     set tracing flag for debugging problem patch interfaces; note:
c     1) itrace should be set < 0 unless debugging is required
c     2) trace output is found in fort.7
c     3) currently trace output is NOT supported in parallel mode
c
c       itrace < 0, do not write search history for current "to" cell
c       itrace = 0, overwrite history from previous "to" cell with current
c       itrace = 1, retain the search history for ALL cells (may get huge file)
c
      itrace = -1
c
#if defined DIST_MPI
      if (ncall .eq. 1) then
         call MPI_Bcast(nintr, 1, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
      end if
#endif
      if (nintr .eq. 0) return
c
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset1   = maxbl
      ioffset2   = intmx
      itag_lw    = 1
      itag_nblkp = itag_lw    + ioffset1
      itag_windx = itag_nblkp + ioffset2
c
      if (ncall .eq. 1) then
         msglen = intmx*(6*msub1+9)
         call MPI_Bcast(iindx, msglen, MPI_INTEGER, myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(llimit, intmx, MPI_INTEGER, myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(iitmax, intmx, MPI_INTEGER, myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(mmcxie, intmx, MPI_INTEGER, myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(mmceta, intmx, MPI_INTEGER, myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(iifit, intmx, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
         call MPI_Bcast(iic0, intmx, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
         call MPI_Bcast(iiorph, intmx, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
         call MPI_Bcast(iitoss, intmx, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
         call MPI_Bcast(ifiner, intmx, MPI_INTEGER, myhost,
     .                  mycomm, ierr)
         call MPI_Bcast(dx, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dy, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dz, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthetx, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthety, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthetz, intmx*msub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
      end if
#endif
c     ioutpt controls output for patch diagnostics...output patch
c     diagnostics only on first call to patcher
c
      ioutpt = 0
      if (ncall .eq. 1) ioutpt = 1
c
c     save number of regular (not dynamic) patch interpolations
c
      if (ncall .eq. 1) then
         ninter0 = ninter
      end if
c
      if (myid.eq.myhost) then
         if (isklton.eq.1) write(11,405) 
      end if
  405 format(/,52h calculating dynamic grid interpolation coefficients)
c
c     zero out dynamic patch interpolation coefficients (first call to
c     patcher only - for subsequent calls, the previous values are used
c     as starting points for the search routine
c
      if (ncall .eq. 1) then
         do 5 ll=1,2
         do 5 mm=1,mxxe
         windx(mm,ll) = 0.
    5    continue
      end if
c
c      set up pre_patch data for dynamic grids
c
      ierrflg = -1
c
      if (ncall .eq. 1) then
         icount_dpat = 0
         do 6809 nbl=1,nblock
#if defined DIST_MPI
            nd_srce = mblk2nd(nbl)
            mytag   = itag_lw + nbl
            if (myid.eq.nd_srce) then
               call MPI_Send(lw,65*maxbl,MPI_INTEGER,myhost,
     .                       mytag,mycomm,ierr)
            else
               if (myid.eq.myhost) then
                  call MPI_Recv(lw_temp,65*maxbl,MPI_INTEGER,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  call pre_patch(nbl, lw_temp, icount_dpat, nintr,
     .                           iindx,intmx,msub1,isav_dpat,
     .                           isav_dpat_b,jjmax1,kkmax1,
     .                           iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                           ierrflg)
               end if
            end if
#else
            call pre_patch(nbl,lw,icount_dpat,nintr,
     .                     iindx,intmx,msub1,isav_dpat,
     .                     isav_dpat_b,jjmax1,kkmax1,
     .                     iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                     ierrflg)
#endif
 6809    continue
c
         if (myid.eq.myhost) then
            maxdcnt = icount_dpat
         end if
c
#if defined DIST_MPI
c        share pre_patch data with nodes
c
         call MPI_Bcast(maxdcnt, 1, MPI_INTEGER, 
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(isav_dpat, intmx*17, MPI_INTEGER, 
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(isav_dpat_b, intmx*msub1*6, MPI_INTEGER,
     .                  myhost, mycomm, ierr)
#endif
      end if
c
c     cycle through blocks, check for and establishing connection
c     information for those blocks with patching
c
      if (myid.eq.myhost) then
         if (ioutpt .gt. 0) write(25,91)
      end if
   91 format(/,/,1x,46hBEGINNING GENERALIZED-COORDINATE INTERPOLATION)
c
      it_thro = 0
      do 7000 nbl=1,nblock
#if defined DIST_MPI
      if (myid.eq.mblk2nd(nbl)) then
#endif
      call lead(nbl,lw,lw2,maxbl)
      it_thro = it_thro + 1
      call patcher(nbl,lw,work,mgwk,wk,nwork,ncall,ioutpt,it_thro,
     .             maxbl,msub1,intmx,mxxe,mptch,jdimg,kdimg,
     .             idimg,windx,nintr,iindx,llimit,iitmax,
     .             mmcxie,mmceta,ncheck,iifit,mblkpt,iic0,
     .             iiorph,iitoss,ifiner,factjlo,factjhi,
     .             factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .             dthetz,isav_dpat,isav_dpat_b,
     .             xte,yte,zte,xmi,ymi,zmi,xmie,ymie,zmie,
     .             jjmax1,kkmax1,jimage,kimage,xorig,yorig,zorig,
     .             jte,kte,sxie,seta,sxie2,seta2,xie2s,eta2s,
     .             temp,x2,y2,z2,nblk1,nblk2,jmm,kmm,x1,y1,z1,
     .             lout,xif1,xif2,etf1,etf2,ireq_ar,
     .             myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .             ibufdim)
#if defined DIST_MPI
      end if
#endif
      if (ioutpt .gt. 0) then
#   ifdef FASTIO
         call writ_buffast(nbl,25,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,33)
         call writ_buffast(nbl,9,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,33)
#   else
         call writ_buf(nbl,25,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
         call writ_buf(nbl,9,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
7000  continue
c
c     append dynamic patch data to regular patch data
c     (see top of subroutine patcher for definition of iindex array)
c
      len0 = 0
      if (ninter0 .gt. 0) then
         nfb0 = iindex(ninter0,1)
         len0 = iindex(ninter0,2*nfb0+4)
     .        + iindex(ninter0,2*nfb0+5) - 1
      end if
      do 1500 n=1,nintr
      iindex(n+ninter0,1)        = iindx(n,1)
      nfb = iindx(n,1)
      if (myid.eq.myhost) then
         do 1501 ll=1,nfb
         iindex(n+ninter0,1+ll)     = iindx(n,1+ll)
 1501    continue
         iindex(n+ninter0,nfb+2)    = iindx(n,nfb+2)
         do 1502 ll=1,nfb
         iindex(n+ninter0,nfb+2+ll) = iindx(n,nfb+2+ll)
 1502    continue
         iindex(n+ninter0,2*nfb+3)  = iindx(n,2*nfb+3)
         iindex(n+ninter0,2*nfb+4)  = iindx(n,2*nfb+4)
         if (n .eq. 1) then
            iindex(n+ninter0,2*nfb+5) = len0 + 1
         else 
            nfb1 = iindex(n+ninter0-1,1)
            iindex(n+ninter0,2*nfb+5)  = iindex(n+ninter0-1,2*nfb1+4)
     .                                 + iindex(n+ninter0-1,2*nfb1+5)
        end if
         iindex(n+ninter0,2*nfb+6)  = iindx(n,2*nfb+6)
         iindex(n+ninter0,2*nfb+7)  = iindx(n,2*nfb+7)
         iindex(n+ninter0,2*nfb+8)  = iindx(n,2*nfb+8)
         iindex(n+ninter0,2*nfb+9)  = iindx(n,2*nfb+9)
      end if
c
      lst  = iindx(n,2*nfb+5)
      len  = lst + iindx(n,2*nfb+4) - 1
c
#if defined DIST_MPI
      nbl_src = iindx(n,nfb+2)
      nd_src  = mblk2nd(nbl_src)
      mytag   = itag_nblkp + n
      msglen  = len - lst + 1
c
      if (myid.eq.myhost) then
         call MPI_Recv(nblkpt(len0+lst), msglen, MPI_INTEGER,
     .                 nd_src, mytag, mycomm, istat, ierr)
      else if (myid.eq.nd_src) then
         call MPI_Send(mblkpt(lst), msglen, MPI_INTEGER, myhost,
     .                 mytag, mycomm, ierr)
      end if
c
      mytag = itag_windx + n
      if (myid.eq.myhost) then
         call MPI_Recv(windex(len0+lst,1), msglen, MY_MPI_REAL,
     .                 nd_src, mytag, mycomm, istat, ierr)
      else if (myid.eq.nd_src) then
         call MPI_Send(windx(lst,1), msglen, MY_MPI_REAL,
     .                 myhost, mytag, mycomm, ierr)
      end if
c
      mytag = itag_windx + maxbl + n
      if (myid.eq.myhost) then
      call MPI_Recv(windex(len0+lst,2), msglen, MY_MPI_REAL,
     .              nd_src, mytag, mycomm, istat, ierr)
      else if (myid.eq.nd_src) then
         call MPI_Send(windx(lst,2), msglen, MY_MPI_REAL,
     .                 myhost, mytag, mycomm, ierr)
      end if
#else
      do 1503 nnn=lst,len
      nblkpt(len0+nnn) = mblkpt(nnn)
 1503 continue
      do 1504 ll=1,2
      do 1505 nnn=lst,len
      windex(len0+nnn,ll) =  windx(nnn,ll)
 1505 continue
 1504 continue
#endif
c
      if (myid.eq.myhost) then
         if (n.eq.nintr .and. (len0+len-1).gt.maxxe) then
            write(11,'('' program terminated in dynamic patching '',
     .            ''routines - see file '',a60)') grdmov
            nou(4) = min(nou(4)+1,ibufdim)
            write(bou(nou(4),4),*)
            nou(4) = min(nou(4)+1,ibufdim)
            write(bou(nou(4),4),1) len0+len-1
    1       format(1x,42hstopping...parameter maxxe is too small...,
     .                  16hmust be at least,i6)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do 1506 ll=1,nfb
         dthetxx(n+ninter0,ll) = dthetx(n,ll)
         dthetyy(n+ninter0,ll) = dthety(n,ll)
         dthetzz(n+ninter0,ll) = dthetz(n,ll)
         iindex(n+ninter0,2*nfb+ll+9) = iindx(n,2*nfb+ll+9)
         iindex(n+ninter0,3*nfb+ll+9) = iindx(n,3*nfb+ll+9)
         iindex(n+ninter0,4*nfb+ll+9) = iindx(n,4*nfb+ll+9)
         iindex(n+ninter0,5*nfb+ll+9) = iindx(n,5*nfb+ll+9)
1506     continue
      end if
1500  continue
c
      if (myid.eq.myhost) then
c
         ninter = ninter0+nintr
c
         if (ioutpt .gt. 0) then
c
c           check iindex array
c
            write(25,*)
            write(25,*)
            write(25,*)' ***** IINDEX ARRAY *****'
            write(25,*)
            write(25,*)' ninter: ',ninter
            write(25,*)
            do 6969 n=1,ninter
            write(25,*)
            write(25,*)' interp. no.: ',n
            write(25,*)
            write(25,*)' nfb: ',iindex(n,1)
            nfb = iindex(n,1)
            write(25,*)' from block(s): ',(iindex(n,l+1),l=1,nfb)
            write(25,*)' to: ',iindex(n,nfb+2)
            write(25,*)' topology (from): ',(iindex(n,l+nfb+2),l=1,nfb)
            write(25,*)' topology (  to): ',iindex(n,2*nfb+3)
            write(25,*)' number of points: ',iindex(n,2*nfb+4)
            write(25,*)' starting index: ',iindex(n,2*nfb+5)
            write(25,*)' xie range: ',iindex(n,2*nfb+6),
     .                   iindex(n,2*nfb+7)
            write(25,*)' eta range: ',iindex(n,2*nfb+8),
     .                   iindex(n,2*nfb+9)
            do 6968 l=1,nfb
            write(25,*)' xie search range in from block ',l,' :',
     .                 iindex(n,2*nfb+9+l),iindex(n,3*nfb+9+l)
            write(25,*)' eta search range, from block ',l,' :',
     .                 iindex(n,4*nfb+9+l),iindex(n,5*nfb+9+l)
 6968       continue
 6969       continue
            write(25,*)
            nfb = iindx(nintr,1)
            lst = iindx(nintr,2*nfb+5)
            len = lst + iindx(nintr,2*nfb+4) - 1
            write(25,9997) len
9997        format('  minimum dimension for parameter MAXXE in CFL3D:',
     .      i6)
            call my_flush(25)
         end if
c
      end if
c
#if defined DIST_MPI
      if (ncall .eq. 1) then 
         call MPI_Bcast(iindex, intmax*(6*nsub1+9), MPI_INTEGER, 
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthetxx, intmax*nsub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthetyy, intmax*nsub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
         call MPI_Bcast(dthetzz, intmax*nsub1, MY_MPI_REAL,
     .                  myhost, mycomm, ierr)
      end if
      call MPI_Bcast(nblkpt, maxxe, MPI_INTEGER,
     .               myhost, mycomm, ierr)
      call MPI_Bcast(windex, 2*maxxe, MY_MPI_REAL,
     .               myhost, mycomm, ierr)
#endif
      return
      end
