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
      subroutine qinter(iseq,lembed,maxl,lw,lw2,w,mgwk,wk,nwork,maxbl,
     .                  maxgr,levelg,igridg,idimg,jdimg,kdimg,ngrid,
     .                  nblg,iemg,inewgg,itrans,irotat,idefrm,xorig,
     .                  yorig,zorig,xorig0,yorig0,zorig0,thetax,thetay,
     .                  thetaz,time2,nou,bou,nbuf,ibufdim,myid,myhost,
     .                  mycomm,mblk2nd,nsegdfrm,idfrmseg,xorgae,
     .                  yorgae,zorgae,thtxae,thtyae,thtzae,maxsegdg,
     .                  nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Interpolate the solution from a coarser mesh to a finer
c     mesh. The finer mesh can be either a global mesh or an embedded
c     mesh. Also updates grid position of finer mesh if meshes are in
c     motion.
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
c
#endif
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension mblk2nd(maxbl),levelg(maxbl),igridg(maxbl),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),nblg(maxgr),
     .          iemg(maxgr),inewgg(maxgr),itrans(maxbl),irotat(maxbl),
     .          xorig(maxbl),yorig(maxbl),zorig(maxbl),xorig0(maxbl),
     .          yorig0(maxbl),zorig0(maxbl),thetax(maxbl),
     .          thetay(maxbl),thetaz(maxbl),time2(maxbl),idefrm(maxbl)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg) 
      dimension xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /maxiv/ ivmx
c
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset  = maxbl
      itag_q   = 1
      itag_qv  = itag_q  + ioffset
      itag_qt  = itag_qv + ioffset
      itag_dat = itag_qt + ioffset
#endif
      mode = 0
      iflg = 0
c
c     interpolate solution to finer mesh - global meshes
c
      if (iseq.le.mseq .and. iseq.ne.1) then
c
         if (myid.eq.myhost) then
            write(11,200)
         end if
         iflg = 1
c
         do 1000 igrid=1,ngrid
c
         nbl  = nblg(igrid)
         iem  = iemg(igrid)
         if (iem.gt.0) go to 1000
c
c        finer mesh nblf, coarser mesh nblz
         nblf = nblg(igrid)+(mseq-iseq)
         nblz = nblf+1
c
         if(mblk2nd(nblf) .eq. myid) then
c
         call lead(nblf,lw,lw2,maxbl)
         iwk1  = 1
         iwk2  = iwk1+jdim*kdim*idim*7
         iwk3  = iwk2+jj2*kk2*ii2*7
         lqc   = lw( 1,nblz)
         lq1c  = lw(16,nblz)
         lturc = lw(19,nblz)
         lvisc = lw(13,nblz)
         nroom = nwork-iwk3
         mdim  = jdim*kk2*ii2*7
         if (nroom.lt.mdim) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),220)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),230)nblz,nblf,igrid
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),240)jdim,kdim,idim
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),250)jj2,kk2,ii2
c
         call addx(w(lq),w(lqc),jdim,kdim,idim,jj2,kk2,ii2,w(lq1c),
     .             wk(iwk1),wk(iwk2),wk(iwk3),nblz,w(lblk),
     .             nou,bou,nbuf,ibufdim,5,myid)
         if (ivmx.ge.4) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),260)
           call addx(w(lxib),w(lturc),jdim,kdim,idim,jj2,kk2,ii2,
     .             w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),nblz,w(lblk),
     .             nou,bou,nbuf,ibufdim,nummem,myid)
           call addx(w(lvis),w(lvisc),jdim,kdim,idim,jj2,kk2,ii2,
     .             w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),nblz,w(lblk),
     .             nou,bou,nbuf,ibufdim,1,myid)
         end if
c
         iuns = max(itrans(nblf),irotat(nblf),idefrm(nblf))
         if (iuns.gt.0) then
c
c           for deforming mesh, interpolate x,y,z (old and new) and
c           boundary deltas up to finer level
c
            if (idefrm(nblf) .gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),265)
               do iseg=1,nsegdfrm(nblf)
                  xorgae(nblf,iseg) = xorgae(nblz,iseg)
                  yorgae(nblf,iseg) = yorgae(nblz,iseg)
                  zorgae(nblf,iseg) = zorgae(nblz,iseg)
                  thtxae(nblf,iseg) = thtxae(nblz,iseg)
                  thtyae(nblf,iseg) = thtyae(nblz,iseg)
                  thtzae(nblf,iseg) = thtzae(nblz,iseg)
               end do
               lxc     = lw(10,nblz)
               lyc     = lw(11,nblz)
               lzc     = lw(12,nblz)
               lxnm2c  = lw(56,nblz)
               lynm2c  = lw(57,nblz)
               lznm2c  = lw(58,nblz)
               ldeltjc = lw(59,nblz) 
               ldeltkc = lw(60,nblz)
               ldeltic = lw(61,nblz)
               call xyzintr(jdim,kdim,idim,jj2,kk2,ii2,w(lx),
     .                      w(ly),w(lz),w(lxc),w(lyc),w(lzc))
               call xyzintr(jdim,kdim,idim,jj2,kk2,ii2,w(lxnm2),
     .                      w(lynm2),w(lznm2),w(lxnm2c),w(lynm2c),
     .                      w(lznm2c))
               call delintr(jdim,kdim,idim,jj2,kk2,ii2,w(ldeltj),
     .                      w(ldeltk),w(ldelti),w(ldeltjc),
     .                      w(ldeltkc),w(ldeltic))
            end if
c
c           update position of finer mesh to current position of 
c           coarser mesh to account for rigid grid motion
c
            if (itrans(nblf).gt.0 .or. irotat(nblf).gt.0) then
               xorig(nblf)  = xorig(nblz)
               yorig(nblf)  = yorig(nblz)
               zorig(nblf)  = zorig(nblz)
               thetax(nblf) = thetax(nblz)
               thetay(nblf) = thetay(nblz)
               thetaz(nblf) = thetaz(nblz)
               time2(nblf)  = time2(nblz)
               call grdmove(nblf,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                      xorig0(nblf),yorig0(nblf),zorig0(nblf),
     .                      xorig(nblf),yorig(nblf),zorig(nblf),
     .                      thetax(nblf),thetay(nblf),thetaz(nblf))
            end if
         end if
c
         end if
c
         nblout = nblg(igrid)
#   ifdef FASTIO
         call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,13)
#   else
         call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
c
 1000    continue
c
         if (myid.eq.myhost) then
            if (lembed.eq.0 .or. iseq.ne.mseq) then
               write(11,210)
               return
            end if
         end if
c
      end if
c
c     embedded grids - either all new  (iseq>1)
c                      or some new     (iseq=1   inewg=1)
c
c     cycle through levels  coarsest to finest
c
      if (lembed.gt.0) then
c
         do 2000 level=mseq+1,maxl
         do 2100 igrid=1,ngrid
         nbl   = nblg(igrid)
         iem   = iemg(igrid)
         inewg = inewgg(igrid)
c
         if (iem.eq.0) go to 2100
         if (iseq.eq.1 .and. inewg.eq.0) go to 2100
         if (levelg(nbl).ne.level) go to 2100
c
#if defined DIST_MPI
         if (myid.ne.myhost) then
#endif
c        embeded mesh nbl; coarser mesh nblc
         call lead(nbl,lw,lw2,maxbl)
         jj2   = jdimg(nblc)
         kk2   = kdimg(nblc)
         ii2   = idimg(nblc)
         lqc   = lw( 1,nblc)
         lq1c  = lw(16,nblc)
         lturc = lw(19,nblc)
         lvisc = lw(13,nblc)
         iwk1  = 1
         iwk2  = iwk1+jdim*kdim*idim*7
         iwk3  = iwk2+jj2*kk2*ii2*7
         iwk4  = iwk3+jdim*kk2*ii2
         nroom = nwork-iwk4
         mdim  = jdim*kdim*ii2
         if (nroom.lt.mdim) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),280)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
#if defined DIST_MPI
         nd_dest = mblk2nd(nbl)
         nd_srce = mblk2nd(nblc)
c
         if (nd_dest.ne.nd_srce) then
c        
c        embedded and coarser block reside on different processors
c
         jki   = jj2*kk2*ii2
         jki2  = jj2*kk2*ii2*nummem
         jki5  = jj2*kk2*ii2*5
c
         mytag   = itag_q + nbl
         if (myid.eq.nd_srce) then
            call MPI_Send(w(lqc),jki5,MY_MPI_REAL,
     .                    nd_dest,mytag,mycomm,ierr)
         else if (myid .eq. nd_dest) then
            iwk5  = iwk4 + jdim*kdim*ii2
            nroom = nwork-iwk5
            if (nroom.lt.jki5) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),280)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            call MPI_Recv(wk(iwk5),jki5,MY_MPI_REAL,
     .                    nd_srce,mytag,mycomm,istat,ierr)
         end if
         iwk6  = iwk5
         iwk7  = iwk5
         if (ivmx.ge.4) then
            mytag = itag_qv + nbl
            if (myid.eq.nd_srce) then
               call MPI_Send(w(lvisc),jki,MY_MPI_REAL,
     .                       nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. nd_dest) then
               iwk6  = iwk5 + jki5
               nroom = nwork-iwk6
               if (nroom.lt.jki) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),280)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               call MPI_Recv(wk(iwk6),jki,MY_MPI_REAL,
     .                       nd_srce,mytag,mycomm,istat,ierr)
            end if
            mytag = itag_qt + nbl
            if (myid.eq.nd_srce) then
               call MPI_Send(w(lturc),jki2,MY_MPI_REAL,
     .                       nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. nd_dest) then
               iwk7  = iwk6 + jki2
               nroom = nwork-iwk7
               if (nroom.lt.jki2) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),280)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               call MPI_Recv(wk(iwk6),jki2,MY_MPI_REAL,
     .                       nd_srce,mytag,mycomm,istat,ierr)
            end if
         end if
c
         if (myid .eq. mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),290)nblc,nbl,igridg(nbl)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),240)jdim,kdim,idim
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),250)jj2,kk2,ii2
            ipass = 1
            call add2x(w(lq),wk(iwk5),jdim,kdim,idim,jj2,kk2,ii2,
     .                 w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .                 js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .                 nou,bou,nbuf,ibufdim,5,myid)
            if (ivmx.ge.4) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),300)
               call add2x(w(lxib),wk(iwk7),jdim,kdim,idim,jj2,kk2,ii2,
     .                 w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .                 js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .                 nou,bou,nbuf,ibufdim,nummem,myid)
               call add2x(w(lvis),wk(iwk6),jdim,kdim,idim,jj2,kk2,ii2,
     .                 w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .                 js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .                 nou,bou,nbuf,ibufdim,1,myid)
            end if
         end if
c
         else
c
c        embedded and coarser block reside on same processor
c
         if (myid .eq. mblk2nd(nbl)) then
c
#endif
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),290)nblc,nbl,igridg(nbl)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),240)jdim,kdim,idim
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),250)jj2,kk2,ii2
         ipass = 1
         call add2x(w(lq),w(lqc),jdim,kdim,idim,jj2,kk2,ii2,
     .              w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .              js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .              nou,bou,nbuf,ibufdim,5,myid)
         if (ivmx.ge.4) then
            nou(1) = min(nou(1)+1,ibufdim)  
            write(bou(nou(1),1),300)
            call add2x(w(lxib),w(lturc),jdim,kdim,idim,jj2,kk2,ii2,
     .              w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .              js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .              nou,bou,nbuf,ibufdim,nummem,myid)
            call add2x(w(lvis),w(lvisc),jdim,kdim,idim,jj2,kk2,ii2,
     .              w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),wk(iwk4),
     .              js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .              nou,bou,nbuf,ibufdim,1,myid)
         end if
#if defined DIST_MPI
c
         end if
c
         end if
#endif
c
         iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
         if (iuns.gt.0) then
c
c           update position of embeded mesh to current position of 
c           coarser (global) mesh
c
#if defined DIST_MPI
            nd_dest = mblk2nd(nbl)
            nd_srce = mblk2nd(nblc)
c
            if (nd_dest.ne.nd_srce) then
c
c           coarse and embedded blocks reside on different processors
c
            mytag = itag_dat + 1
            if (myid.eq.nd_srce) then
               call MPI_Send(xorig,3,MY_MPI_REAL,
     .                       nd_dest,mytag,mycomm,ierr)
            else if (myid.eq.nd_dest) then
               call MPI_Recv(xorig,3,MY_MPI_REAL,
     .                       nd_srce,mytag,mycomm,istat,ierr)
            end if
            mytag = itag_dat + 2
            if (myid.eq.nd_srce) then
               call MPI_Send(thetax,3,MY_MPI_REAL,
     .                       nd_dest,mytag,mycomm,ierr)
            else if (myid.eq.nd_dest) then
               call MPI_Recv(thetax,3,MY_MPI_REAL,
     .                       nd_srce,mytag,mycomm,istat,ierr)
            end if
            mytag = itag_dat + 3
            if (myid.eq.nd_srce) then
               call MPI_Send(time2,1,MY_MPI_REAL,
     .                       nd_dest,mytag,mycomm,ierr)
            else if (myid.eq.nd_dest) then
               call MPI_Recv(time2,1,MY_MPI_REAL,
     .                       nd_srce,mytag,mycomm,istat,ierr)
            end if
c
            else
c
c           coarse and embedded blocks reside on same processor
c
            if (myid .eq. mblk2nd(nbl)) then
c
#endif
            xorig(nbl)  = xorig(nblc)
            yorig(nbl)  = yorig(nblc)
            zorig(nbl)  = zorig(nblc)
            thetax(nbl) = thetax(nblc)
            thetay(nbl) = thetay(nblc)
            thetaz(nbl) = thetaz(nblc)
            time2(nbl)  = time2(nblc)
c
            if (idefrm(nbl) .eq. 0) then
               call grdmove(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                      xorig0(nbl),yorig0(nbl),zorig0(nbl),
     .                      xorig(nbl),yorig(nbl),zorig(nbl),
     .                      thetax(nbl),thetay(nbl),thetaz(nbl))
            end if
#if defined DIST_MPI
c
            end if
c
            end if
#endif
         end if
#if defined DIST_MPI
         end if
#endif
c
         nblout = nblg(igrid)
#   ifdef FASTIO
         call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,14)
#   else
         call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
c
 2100    continue
 2000    continue
c
         if (myid.eq.myhost) then
            if (iflg.gt.0) write(11,210)
         end if
c
      end if
c
  200 format(47h***** BEGINNING SEQUENCING TO FINER LEVEL *****/)
  210 format(/44h***** ENDING SEQUENCING TO FINER LEVEL *****)
  220 format(42h not enough work space for subroutine addx)
  230 format(40h interpolating solution on coarser block,i4,
     .       18h  to   finer block,i4,6h (grid,i4,1h))
  240 format(31h   jdim,kdim,idim (finer grid)=,3i5)
  250 format(31h   jj2,kk2,ii2  (coarser grid)=,3i5)
  260 format(45h   interpolating turb quantities from coarser,
     .       15h to finer block)
  265 format(49h   interpolating deforming mesh data from coarser,
     .       15h to finer block)
  270 format(43h stopping...no way to mesh sequence up from,
     .       34h starting level for deforming mesh)
  280 format(43h not enough work space for subroutine add2x)
  290 format(40h interpolating solution on coarser block,i4,
     .       18h  to embeded block,i4,6h (grid,i4,1h))
  300 format(45h   interpolating turb quantities from coarser,
     .       17h to embeded block)
c
      return
      end
