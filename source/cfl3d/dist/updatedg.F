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
      subroutine updatedg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                    maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                    nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                    nblock,levelg,igridg,idefrm,ncgg,iadvance,nou,
     .                    bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                    utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                    xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                    rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,itrans,irotat,islavept,
     .                    nslave,iskip,jskip,kskip,xs,xxn,nsegdfrm,
     .                    idfrmseg,iaesurf,maxsegdg,iwk,nmaster,nt,
     .                    xorig,yorig,zorig,xorgae0,yorgae0,zorgae0,
     .                    icouple,ireq,nnodes,nblelst,iskmax,jskmax,
     .                    kskmax,ue)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Update a deforming grid to a new shape, and obtain 
c     corresponding grid-boundary velocities for use in the boundary
c     conditions. Also collocate new grid position to coarser levels
c     levels and obtain grid-boundary velocities on coarser levels.
c
c     surface shapes are updated on a segment basis, according to the
c     general rule:
c
c     idfrmseg(nbl,is) <  99...deformation via prescribed motion
c                      =  99...deformation via aeroelastic motion
c
c     where nbl denotes the block and is denotes the segment
c
c     currently, the following specfic motions are supported:
c
c     idfrmseg(nbl,is) =   1...translation only
c                          2...rotation (+/- symmetric) only
c                         99...aeroelastic only
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
      character*120 bou(ibufdim,nbuf)
      character*35 cmd1,cmd2
      character*3 char1
c
      dimension nou(nbuf)
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),wk(nwork),iwk(maxbl)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension levelg(maxbl),igridg(maxbl),idefrm(maxbl)
      dimension itrans(maxbl),irotat(maxbl)
      dimension ncgg(maxgr),iadvance(maxbl),mblk2nd(maxbl)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4),xs(2*nmds,maxaes),
     .          xxn(2*nmds,maxaes)
      dimension aesrfdat(5,maxaes),islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave),ue(3*nslave)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension xorgae0(maxbl,maxsegdg),yorgae0(maxbl,maxsegdg),
     .          zorgae0(maxbl,maxsegdg),icouple(maxbl,maxsegdg)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension ireq(maxbl)
      dimension nblelst(maxbl,2)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
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
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /elastic_ss/ idef_ss
c
      if (isklton.eq.1) then
#if defined DIST_MPI
         if (myid.eq.1) then
#endif
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),400)
  400    format(1x,1h )
#if defined DIST_MPI
         end if
#endif
#   ifdef FASTIO
         call writ_buffast(1,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,29)
#   else
         call writ_buf(1,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset  = maxbl
      itag_ns  = 1
      itag_wk  = itag_ns + ioffset
      itag_org = itag_wk + ioffset
      itag     = 0
      nreq     = 1
c
c     make sure all current block origins are available on this
c     processor
c
      if (myid .ne. myhost) then
         ist = 1
         do nbl=1,nblock
            itag = itag + 1
            if (nbl.gt.1) ist  = ist + 1
            if (myid .ne. mblk2nd(nbl)) then

c              recieve the data on non-local nodes
c
               nd_srce = mblk2nd(nbl)
               mytag   = itag_org + itag
               call MPI_Irecv (wk(ist),3,MY_MPI_REAL,
     .                         nd_srce,mytag,mycomm,ireq(nreq),ierr)
               nreq = nreq + 1
            else
c
c              send the data to non-local nodes
c
               do inode = 1,nnodes
                  if (inode .ne. myid) then
                     nd_dest = inode
                     mytag   = itag_org + itag
                     wk(ist)   = xorig(nbl)
                     wk(ist+1) = yorig(nbl)
                     wk(ist+2) = zorig(nbl)
                     call MPI_Send (wk(ist),3,MY_MPI_REAL,
     .                              nd_dest,mytag,mycomm,ierr)
                  end if
               end do
            end if
            if (nreq-1 .gt. 0) then
               call MPI_Wait(ireq(nreq-1),istat,ierr)
            end if
            xorig(nbl) = wk(ist)
            yorig(nbl) = wk(ist+1)
            zorig(nbl) = wk(ist+2)
         end do
      end if
#endif
c
c     update surface positions due to forced motion and/or
c     aeroelastic motion
c
      icnt  = 0
      do 100 nbl = 1,nblock
c
      if (myid.eq.mblk2nd(nbl) .and.  (levelg(nbl).ge.lglobal .and.
     .    levelg(nbl).le.levelt(iseqr))) then
c
         call lead(nbl,lw,lw2,maxbl)
c
c        zero out deltj, deltk, delti from previous step
c        (except if this is a static deformation case)
c
         if (idefrm(nbl).gt.0 .and. idef_ss .eq. 0) then
            do lll=1,kdim*idim*3*2
               w(ldeltj+lll-1) = 0.
            end do
            do lll=1,jdim*idim*3*2
               w(ldeltk+lll-1) = 0.
            end do
            do lll=1,kdim*jdim*3*2
               w(ldelti+lll-1) = 0.
            end do
         end if
c
         if (idefrm(nbl).gt.0 .and. idefrm(nbl).lt.999) then
c
            if (isklton.eq.1)then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),440) nbl
            end if
  440       format(1x,28hdeforming surface of block  ,i4,
     .             13h to new shape)
c
            lwk1 = 1
            lwk2 = lwk1 + kdim*idim*2
            lwk3 = lwk2 + jdim*idim*2
            lwk4 = lwk3 + jdim*kdim*2
            if (nwork.lt.lwk4) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),444)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),443) lwk4 - nwork
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
  444       format(38h not enough work space for subroutines,
     .             23h trnsurf/rotsurf/aesurf)
  443       format(' Increase keyword memadd allocation by ',i8)
c
c           wk(lwk1-lwk4) contain a flag that prevents points common
c           to multiple segments from being updated more than once
c           per motion type (e.g. translation and rotation). this
c           array must be set to 1 before each call to a specific
c           surface motion type
c
            do lll=lwk1,lwk4
               wk(lll) = 1.
            end do
c
c           update segment origin in case forced deformation is 
c           coupled to either forced motion of block origin or
c           if deforming rotation is coupled to deforming translation
c
            do is = 1,nsegdfrm(nbl)
               if (icouple(nbl,is) .ne. 0) then
                  if (idfrmseg(nbl,is).gt.0 .and. 
     .               idfrmseg(nbl,is).lt.99) then
                     if (icouple(nbl,is) .ne. 0) then
                        nblmast = icouple(nbl,is)
                        xorgae(nbl,is) = xorig(nblmast)
                        yorgae(nbl,is) = yorig(nblmast)
                        zorgae(nbl,is) = zorig(nblmast)
                     end if
                   end if
               end if
            end do
c
c           surface translation
c
            do is = 1,nsegdfrm(nbl)
               if (idfrmseg(nbl,is) .eq. 1) then
                  call trnsurf(jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                         w(ldeltj),w(ldeltk),w(ldelti),nbl,
     .                         idfrmseg(nbl,is),xorgae(nbl,is),
     .                         yorgae(nbl,is),zorgae(nbl,is),
     .                         utrnsae(nbl,is),vtrnsae(nbl,is),
     .                         wtrnsae(nbl,is),rfrqtae(nbl,is),
     .                         icsi(nbl,is),icsf(nbl,is),jcsi(nbl,is),
     .                         jcsf(nbl,is),kcsi(nbl,is),kcsf(nbl,is),
     .                         time,nou,bou,nbuf,ibufdim,myid,wk(lwk1),
     .                         wk(lwk2),wk(lwk3),xorgae0(nbl,is),
     .                         yorgae0(nbl,is),zorgae0(nbl,is))
               end if
            end do
c
            do lll=lwk1,lwk4
               wk(lll) = 1.
            end do
c
c           surface rotation
c
            do is = 1,nsegdfrm(nbl)
               if (idfrmseg(nbl,is) .eq. 2) then
                  call rotsurf(jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                         w(ldeltj),w(ldeltk),w(ldelti),nbl,
     .                         idfrmseg(nbl,is),xorgae(nbl,is),
     .                         yorgae(nbl,is),zorgae(nbl,is),
     .                         omgxae(nbl,is),omgyae(nbl,is),
     .                         omgzae(nbl,is),thtxae(nbl,is),
     .                         thtyae(nbl,is),thtzae(nbl,is),
     .                         rfrqrae(nbl,is),icsi(nbl,is),
     .                         icsf(nbl,is),jcsi(nbl,is),jcsf(nbl,is),
     .                         kcsi(nbl,is),kcsf(nbl,is),time,nou,bou,
     .                         nbuf,ibufdim,myid,wk(lwk1),wk(lwk2),
     .                         wk(lwk3))
               end if
            end do
c
            do lll=lwk1,lwk4
               wk(lll) = 1.
            end do
c
c           aeroelastic surface 
c
            call aesurf(nbl,jdim,kdim,idim,w(ldeltj),w(ldeltk),
     .                  w(ldelti),w(lxmdj),w(lxmdk),w(lxmdi),
     .                  wk(lwk1),wk(lwk2),wk(lwk3),maxbl,maxseg,
     .                  nmds,maxaes,aesrfdat,xs,xxn,icsi,icsf,
     .                  jcsi,jcsf,kcsi,kcsf,nsegdfrm,idfrmseg,
     .                  iaesurf,maxsegdg)
c
            do lll=lwk1,lwk4
               wk(lll) = 1.
            end do
c
         end if
      end if
c
      if (isklton.eq.1) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,30)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
  100 continue
c
#if defined DIST_MPI
      itag     = 0
      nreq     = 1
c
c     make sure all current segment origins are available on this
c     processor
c
      if (myid .ne. myhost) then
         ist = 1
         do nbl=1,nblock
            itag = itag + 1
            if (nbl .gt. 1) ist  = ist + 3*nsegdfrm(nbl-1)
            if (myid .ne. mblk2nd(nbl)) then

c              recieve the data on non-local nodes
c
               nd_srce = mblk2nd(nbl)
               mytag   = itag_org + itag
               nval    = 3*nsegdfrm(nbl)
               call MPI_Irecv (wk(ist),nval,MY_MPI_REAL,
     .                         nd_srce,mytag,mycomm,ireq(nreq),ierr)
               nreq = nreq + 1
            else
c
c              send the data to non-local nodes
c
               do inode = 1,nnodes
                  if (inode .ne. myid) then
                     nd_dest   = inode
                     mytag     = itag_org + itag
                     nval      = 3*nsegdfrm(nbl)
                     do is=1,nsegdfrm(nbl)
                        ist2 = 3*(is-1)
                        wk(ist+ist2)   = xorgae(nbl,is)
                        wk(ist+ist2+1) = yorgae(nbl,is)
                        wk(ist+ist2+2) = zorgae(nbl,is)
                     end do
                     call MPI_Send (wk(ist),nval,MY_MPI_REAL,
     .                              nd_dest,mytag,mycomm,ierr)
                  end if
               end do
            end if
            if (nreq-1 .gt. 0) then
               call MPI_Wait(ireq(nreq-1),istat,ierr)
            end if
            do is=1,nsegdfrm(nbl)
               ist2 = 3*(is-1)
               xorgae(nbl,is) = wk(ist+ist2)
               yorgae(nbl,is) = wk(ist+ist2+1)
               zorgae(nbl,is) = wk(ist+ist2+2)
            end do
         end do
      end if
#endif
c
c     update block origin in case motion of block origin is
c     coupled to forced deformation. note: only translational
c     motion updates block origin - rotational motion leaves
c     origin unchanged (i.e origin is center of roation)
c
      do nbl=1,nblock
         if (myid.eq.mblk2nd(nbl) .and. (levelg(nbl).ge.lglobal .and.
     .       levelg(nbl).le.levelt(iseqr))) then
             do is = 1,nsegdfrm(nbl)
                if (icouple(nbl,is) .ne. 0) then
                   nblmast = icouple(nbl,is)
                   if (idfrmseg(nblmast,is) .eq. 1) then
                      xorig(nbl) = xorgae(nblmast,is)
                      yorig(nbl) = yorgae(nblmast,is)
                      zorig(nbl) = zorgae(nblmast,is)
                   end if
                end if
             end do
         end if
      end do
c
c     get a list of all points on deforming solid surfaces. the
c     list structure is as folllows: each successive 9 entries
c     in the list give, in order, the x, y, z, deltx, delty, deltz,
c     xnm1, ynm1, znm1 values for the surface point (where xnm1, 
c     etc. are the points at time n-1); these 9 data are repeated
c     for each solid surface point that undergoes deformation. Thus,
c     if there are a total of nsurf solid surface points that undergo
c     deformation, the list will be of dimension 9*nsurf
c 
      nsurf = 0
      do nbl = 1,nblock
         iwk(nbl) = 0
         if ((levelg(nbl).ge.lglobal .and.
     .       levelg(nbl).le.levelt(iseqr))) then
             if (idefrm(nbl).gt.0 .and. idefrm(nbl).lt.999) then
                call lead(nbl,lw,lw2,maxbl)
                if (myid.eq.mblk2nd(nbl)) then
                   call getsurf(w(lx),w(ly),w(lz),w(ldeltj),
     .                          w(ldeltk),w(ldelti),w(lxnm1),
     .                          w(lynm1),w(lznm1),icsi,icsf,
     .                          jcsi,jcsf,kcsi,kcsf,wk,nwork,
     .                          nbl,idim,jdim,kdim,nsurf,iwk(nbl),
     .                          nsegdfrm,maxbl,idfrmseg,maxsegdg)
                end if
#if defined DIST_MPI
                mytag = itag_ns + nbl
                nd_srce = mblk2nd(nbl)
                if (myid.eq.mblk2nd(nbl)) then
                   call MPI_Send(iwk(nbl),1,MPI_INTEGER,myhost,mytag,
     .                           mycomm,ierr)
                else if (myid.eq.myhost) then
                   call MPI_Recv(iwk(nbl),1,MPI_INTEGER,nd_srce,mytag,
     .                           mycomm,istat,ierr)
                end if
                if (iwk(nbl).gt.0) then
                   mytag   = itag_wk + nbl
                   nd_srce = mblk2nd(nbl)
                   numdat  = 9*iwk(nbl)
                   if (myid.eq.mblk2nd(nbl)) then
                      call MPI_Send(wk,numdat,MY_MPI_REAL,
     .                              myhost,mytag,mycomm,ierr) 
                   else if (myid.eq.myhost) then
                      ns = 9*nsurf+1
                      call MPI_Recv(wk(ns),numdat,MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,istat,ierr)
                   end if
                end if
#endif
                if (myid.eq.myhost) then
                   nsurf = nsurf + iwk(nbl)
                end if
             end if
         end if
      end do
c
#if defined DIST_MPI
      call MPI_Bcast (nsurf,1,MPI_INTEGER,myhost,mycomm,ierr) 
      call MPI_Bcast (wk,nsurf*9,MY_MPI_REAL,myhost,
     .                mycomm,ierr)
#endif
      if (iunst.gt.1 .or. idef_ss.gt.0) then
 
         do 150 nbl = 1,nblock
         if (myid.eq.mblk2nd(nbl) .and. (levelg(nbl).ge.lglobal .and.
     .       levelg(nbl).le.levelt(iseqr))) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            if (idefrm(nbl) .gt. 0) then
c
c              temporary storage locations:
c              lvel  = start of grid point velocity array
c              lacci = start of i-boundary point acceleration array
c              laccj = start of j-boundary point acceleration array
c              lacck = start of k-boundary point acceleration array
c              lt1wk = start of work array for subroutine metric
c              lt2wk = start of work array for subroutine metric
c              lt3wk = start of work array for subroutine metric
c              lt4wk = start of work array for subroutine cellvol
c
               lvel  = 9*nsurf + 1
               lacci = jdim*kdim*idim*3+lvel
               laccj = jdim*kdim*3*2+lacci
               lacck = kdim*idim*3*2+laccj
               lt1wk = jdim*idim*3*2+lacck
               lt2wk = jdim*kdim*idim*3+lt1wk
               lt3wk = jdim*kdim*6+lt2wk
               lt4wk = jdim*kdim*idim*5+lt3wk 
c
               mdim  = jdim*kdim*idim*3
               if (nwork.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),445)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),446) mdim - nwork
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
c
               if (isklton.eq.1)then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'('' deforming block  '',i4,
     .                  '' to new shape'')') nbl
               end if
c

               nflag = 0
               irst = 0
               if (idef_ss .gt. 0) irst = 1
               call deform_surf(nbl,idim,jdim,kdim,w(ldeltj),w(ldeltk),
     .                          w(ldelti),lw,lw2,icsi,icsf,jcsi,jcsf,
     .                          kcsi,kcsf,maxbl,mseq,time,dt,ita,nou,
     .                          bou,nbuf,ibufdim,myid,idefrm,nbci0,
     .                          nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                          ibcinfo,jbcinfo,kbcinfo,maxseg,wk,ue,
     .                          nsurf,irst,nflag,islavept,nslave,
     .                          nsegdfrm,idfrmseg,iaesurf,
     .                          maxsegdg,nmaster,iseqr)

c
            end if
         end if
150      continue
      end if
c
#if defined DIST_MPI
      do nn = 1,maxbl
        if(nblelst(nn,1).ne.0) then
          nbl   = nblelst(nn,1)
          n1    = nblelst(nn,2)
          n2    = nblelst(nn+1,2)
          nstrt = 3*(n1-1) + 1  
          nend  = 3*(n2-1) + 1  
          ndata = nend - nstrt
          if (myid .eq. mblk2nd(nbl)) then 
            mytag   = mblk2nd(nbl)
            nd_dest = myhost
            call MPI_Send (ue(nstrt),ndata,MY_MPI_REAL,
     .                     nd_dest,mytag,mycomm,ierr)
          else if (myid .eq. myhost) then 
            mytag   = mblk2nd(nbl)
            nd_srce = mblk2nd(nbl)
            call MPI_Recv (ue(nstrt),ndata,MY_MPI_REAL,
     .                     nd_srce,mytag,mycomm,istat,ierr)
          end if
        end if
      enddo
      call MPI_Bcast (ue,3*nslave,MY_MPI_REAL,myhost,
     .                mycomm,ierr)
#endif
c
c
c     update faces and interior of deforming meshes, and finite
c     difference old and new grid positions to get new grid velocities
c     (velocities stored in work array, at position lvel)
c
      if (iunst.gt.1 .or. idef_ss.gt.0) then
 
           do n = 1,nslave
             nbl = islavept(n,9,iseqr)
             ll  = islavept(n,1,iseqr)
             if (myid .eq. mblk2nd(nbl)) then
               call lead(nbl,lw,lw2,maxbl)
               xst(n) = w(lx+ll)  
               yst(n) = w(ly+ll)  
               zst(n) = w(lz+ll)  
             end if
           enddo
#if defined DIST_MPI
           do nn = 1,maxbl
             if(nblelst(nn,1).ne.0) then
               nbl   = nblelst(nn,1)
               nstrt = nblelst(nn,2) 
               nend  = nblelst(nn+1,2) 
               ndata = nend - nstrt
               if (myid .eq. mblk2nd(nbl)) then
                   mytag   = mblk2nd(nbl)
                   nd_dest = myhost
                   call MPI_Send (xst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_dest,
     .                            mytag,mycomm,ierr)
               else if (myid .eq. myhost) then
                   mytag   = mblk2nd(nbl)
                   nd_srce = mblk2nd(nbl)
                   call MPI_Recv (xst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_srce,mytag,
     .                            mycomm,istat,ierr)
               end if
             end if
           enddo
           do nn = 1,maxbl
             if(nblelst(nn,1).ne.0) then
               nbl   = nblelst(nn,1)
               nstrt = nblelst(nn,2) 
               nend  = nblelst(nn+1,2) 
               ndata = nend - nstrt
               if (myid .eq. mblk2nd(nbl)) then
                   mytag   = mblk2nd(nbl)
                   nd_dest = myhost
                   call MPI_Send (yst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_dest,
     .                            mytag,mycomm,ierr)
               else if (myid .eq. myhost) then
                   mytag   = mblk2nd(nbl)
                   nd_srce = mblk2nd(nbl)
                   call MPI_Recv (yst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_srce,mytag,
     .                            mycomm,istat,ierr)
               end if
             end if
           enddo
           do nn = 1,maxbl
             if(nblelst(nn,1).ne.0) then
               nbl   = nblelst(nn,1)
               nstrt = nblelst(nn,2) 
               nend  = nblelst(nn+1,2) 
               ndata = nend - nstrt
               if (myid .eq. mblk2nd(nbl)) then
                   mytag   = mblk2nd(nbl)
                   nd_dest = myhost
                   call MPI_Send (zst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_dest,
     .                            mytag,mycomm,ierr)
               else if (myid .eq. myhost) then
                   mytag   = mblk2nd(nbl)
                   nd_srce = mblk2nd(nbl)
                   call MPI_Recv (zst(nstrt),ndata,
     .                            MY_MPI_REAL,nd_srce,mytag,
     .                            mycomm,istat,ierr)
               end if
             end if
           enddo
           call MPI_Bcast (xst,nslave,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
           call MPI_Bcast (yst,nslave,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
           call MPI_Bcast (zst,nslave,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
#endif
         call deform_el(islavept,nslave,nmaster,ue,xst,yst,zst,
     .                  nt,myhost,mycomm,myid,nnodes,mblk2nd,nblelst,
     .                  maxbl,iseqr)
      end if
      if (iunst.gt.1 .or. idef_ss.gt.0) then
         ivert = 0
         do 200 nbl = 1,nblock
         if (myid.eq.mblk2nd(nbl) .and. (levelg(nbl).ge.lglobal .and.
     .       levelg(nbl).le.levelt(iseqr))) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            if (idefrm(nbl) .gt. 0) then
c
c              temporary storage locations:
c              lvel  = start of grid point velocity array
c              lacci = start of i-boundary point acceleration array
c              laccj = start of j-boundary point acceleration array
c              lacck = start of k-boundary point acceleration array
c              lt1wk = start of work array for subroutine metric
c              lt2wk = start of work array for subroutine metric
c              lt3wk = start of work array for subroutine metric
c              lt4wk = start of work array for subroutine cellvol
c
               lvel  = 9*nsurf + 1
               lacci = jdim*kdim*idim*3+lvel
               laccj = jdim*kdim*3*2+lacci
               lacck = kdim*idim*3*2+laccj
               lt1wk = jdim*idim*3*2+lacck
               lt2wk = jdim*kdim*idim*3+lt1wk
               lt3wk = jdim*kdim*6+lt2wk
               lt4wk = jdim*kdim*idim*5+lt3wk 
c
               mdim  = jdim*kdim*idim*3
               if (nwork.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),445)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),446) mdim - nwork
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  445          format(37h not enough work space for subroutine,
     .                 7h deform)
  446       format(' Increase keyword memadd allocation by ',i8)
c
               if (isklton.eq.1)then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'('' deforming block  '',i4,
     .                  '' to new shape'')') nbl
               end if
c
               nflag = 0
               irst = 0
               if (idef_ss .gt. 0) irst = 1
               call deform(nbl,idim,jdim,kdim,w(lx),w(ly),w(lz),
     .                     w(lxnm2),w(lynm2),w(lznm2),w(lxnm1),
     .                     w(lynm1),w(lznm1),w(ldeltj),w(ldeltk),
     .                     w(ldelti),ue,wk(lvel),icsi,icsf,jcsi,
     .                     jcsf,kcsi,kcsf,maxbl,dt,nou,bou,
     .                     nbuf,ibufdim,myid,idefrm,nbci0,nbcidim,
     .                     nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,
     .                     kbcinfo,maxseg,wk,nsurf,irst,nflag,
     .                     islavept,nslave,iskip,jskip,kskip,nsegdfrm,
     .                     idfrmseg,iaesurf,maxsegdg,nmaster,iseqr,
     .                     iskmax,jskmax,kskmax,nt)
               if (iunst .gt. 1) then
c
c              calculate face-average values of velocity and acceleration
c              on block boundaries and place in permanent storage for use
c              in boundary condition routines
c
c              set boundary accelerations to zero...this is done for
c              rigid rotation/translation too
c
               do lll=lacci,lt1wk-1
                  wk(lll) = 0.
               end do
c
               call xtbatb(jdim,kdim,idim,w(lxtbj),w(lxtbk),w(lxtbi),
     .                     w(latbj),w(latbk),w(latbi),wk(lvel),
     .                     wk(lacci),wk(laccj),wk(lacck))
c
               nroom = nwork-lt3wk
               mdim  = jdim*kdim*idim*5
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),425)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),426) mdim - nroom
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  425          format(45h not enough work space for metric subroutines)
  426          format(' Increase keyword memadd allocation by ',i8)
c
c              calculate spatial metrics for updated grid
c
               iflag = -1
               call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                     w(lsk),w(lsi),wk(lt2wk),wk(lt3wk),nbl,
     .                     iflag,icnt,nbci0,nbcj0,nbck0,nbcidim,
     .                     nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                     maxbl,maxseg,nou,bou,nbuf,ibufdim,myid,
     .                     mblk2nd)
c
c              calculate temporal metrics for updated grid
c
               call tmetric(jdim,kdim,idim,w(lsj),w(lsk),w(lsi),
     .                      w(lx),w(ly),w(lz),wk(lvel),wk(lt1wk),
     .                      wk(lt2wk),wk(lt3wk),nbl)
c
c              calculate volumes for updated grid
c
               nroom = nwork-lt4wk
               mdim = jdim*kdim*15
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),430)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),431) mdim - nroom 
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  430          format(45h not enough work space for subroutine cellvol)
  431       format(' Increase keyword memadd allocation by ',i8)
c
               iflagv = 0
               call cellvol(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                      w(lsk),w(lsi),w(lvol),wk(lt4wk),nou,bou,
     .                      nbuf,ibufdim,myid,mblk2nd,maxbl,nbl,
     .                      1,iflagv,imin,imax,jmin,jmax,kmin,kmax)
c    .                      0,iflagv,imin,imax,jmin,jmax,kmin,kmax)
               iflagv = 0
               if(iflagv.eq.1) then
                 idim1=imax-imin+1
                 jdim1=jmax-jmin+1
                 kdim1=kmax-kmin+1
                 call s_gwrite(nbl,idim,jdim,kdim,imin,imax,jmin,jmax,
     .                         kmin,kmax,idim1,jdim1,kdim1,w(lx),w(ly),
     .                         w(lz))
                 call gridgen_glyf1(nbl)
                 call gridgen_glyf2(nbl)
                 cmd1='gridgen -b dgplot3d_gridgen0001.glf'
                 cmd2='gridgen -b dgplot3d_gridgen0002.glf'
                 if(nbl.lt.10) then
                   write(char1,'(i1)') nbl
                   cmd1(30:30) = char1(1:1)
                   cmd2(30:30) = char1(1:1)
                 else if(nbl.lt.100) then
                   write(char1,'(i2)') nbl
                   cmd1(29:30) = char1(1:2)
                   cmd2(29:30) = char1(1:2)
                 else
                   write(char1,'(i3)') nbl
                   cmd1(28:30) = char1(1:3)
                   cmd2(28:30) = char1(1:3)
                 end if
                 Isys=system(cmd1)
                 Isys=system(cmd2)
                 call s_gread(nbl,idim,jdim,kdim,imin,imax,jmin,jmax,
     .                        kmin,kmax,idim1,jdim1,kdim1,w(lx),w(ly),
     .                        w(lz))
c
                 iflag = -1
                 call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                       w(lsk),w(lsi),wk(lt2wk),wk(lt3wk),nbl,
     .                       iflag,icnt,nbci0,nbcj0,nbck0,nbcidim,
     .                       nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                       maxbl,maxseg,nou,bou,nbuf,ibufdim,myid,
     .                       mblk2nd)
c
c                calculate temporal metrics for updated grid
c
                 call tmetric(jdim,kdim,idim,w(lsj),w(lsk),w(lsi),
     .                        w(lx),w(ly),w(lz),wk(lvel),wk(lt1wk),
     .                        wk(lt2wk),wk(lt3wk),nbl)
c
c                calculate volumes for updated grid
c
                 iflagv = 0
                 call cellvol(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                        w(lsk),w(lsi),w(lvol),wk(lt4wk),nou,bou,
     .                        nbuf,ibufdim,myid,mblk2nd,maxbl,nbl,
     .                        1,iflagv,imin,imax,jmin,jmax,kmin,kmax)
               end if 
c
c              coarser levels
c
               ncg = ncgg(igridg(nbl)) - (mseq-iseqr)
               if (ncg.gt.0 .and. mgflag.gt.0) then
                  nbll = nbl
                  do 1820 m=1,ncg
                  nbll = nbll+1
                  nbllm1 = nbll - 1
                  do is=1,nsegdfrm(nbl)
                     xorgae(nbll,is)  = xorgae(nbl,is)
                     yorgae(nbll,is)  = yorgae(nbl,is)
                     zorgae(nbll,is)  = zorgae(nbl,is)
                     thtxae(nbll,is) = thtxae(nbl,is)
                     thtyae(nbll,is) = thtyae(nbl,is)
                     thtzae(nbll,is) = thtzae(nbl,is)
                  end do
                  lvolc  = lw( 8,nbll)
                  lxc    = lw(10,nbll)
                  lyc    = lw(11,nbll)
                  lzc    = lw(12,nbll)
                  lxtbjc = lw(36,nbll)
                  lxtbkc = lw(37,nbll)
                  lxtbic = lw(38,nbll)
                  latbjc = lw(39,nbll)
                  latbkc = lw(40,nbll)
                  latbic = lw(41,nbll)
                  lvelc  = lt1wk
c
                  if(isklton.eq.1)then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),850) nbll,ii2,jj2,kk2
                  end if
  850             format(1x,24h  creating coarser block,i4,
     .            24h of dimensions (I/J/K) :,3i4)
c
c                 collocate xyz
                  call collx(w(lx),w(ly),w(lz),w(lxc),w(lyc),w(lzc),
     .                       jdim,kdim,idim,jj2,kk2,ii2)
c
c                 collocate grid point velocity
                  call collxt(wk(lvel),wk(lvelc),jdim,kdim,idim,
     .                        jj2,kk2,ii2,nbllm1,nou,bou,nbuf,ibufdim)
                  nv = jj2*kk2*ii2*3
                  do 1825 izz = 1,nv
                  wk(lvel+izz-1) = wk(lvelc+izz-1)
 1825             continue
c
c                 collocate volumes
                  call collv(w(lvol),w(lvolc),jdim,kdim,idim,jj2,
     .                       kk2,ii2)
c
c                 collocate i0/idim boundary velocity/acceleration
                  call collxtb(w(lxtbi),w(lxtbic),jdim,kdim,
     .                        jj2,kk2,nbllm1)
                  call collxtb(w(latbi),w(latbic),jdim,kdim,
     .                        jj2,kk2,nbllm1)
c
c                 collocate j0/jdim boundary velocity/acceleration
                  call collxtb(w(lxtbj),w(lxtbjc),kdim,idim-1,
     .                        kk2,ii2-1,nbllm1)
                  call collxtb(w(latbj),w(latbjc),kdim,idim-1,
     .                        kk2,ii2-1,nbllm1)
c
c                 collocate k0/kdim boundary velocity/acceleration
                  call collxtb(w(lxtbk),w(lxtbkc),jdim,idim-1,
     .                        jj2,ii2-1,nbllm1)
                  call collxtb(w(latbk),w(latbkc),jdim,idim-1,
     .                        jj2,ii2-1,nbllm1)
c
c                 calculate spatial metrics for updated coarser grid
c
                  call lead(nbll,lw,lw2,maxbl)
c
                  lvel  = 1 + 9*nsurf
                  lt1wk = jdim*kdim*idim*3+lvel
                  lt2wk = jdim*kdim*idim*3+lt1wk
                  lt3wk = jdim*kdim*6+lt2wk
c
                  iflag = -1
                  call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                        w(lsk),w(lsi),wk(lt2wk),wk(lt3wk),nbll,
     .                        iflag,icnt,nbci0,nbcj0,nbck0,nbcidim,
     .                        nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                        maxbl,maxseg,nou,bou,nbuf,ibufdim,myid,
     .                        mblk2nd)
c
c                 calculate temporal metrics for updated coarser grid
c
                  call tmetric(jdim,kdim,idim,w(lsj),w(lsk),w(lsi),
     .                         w(lx),w(ly),w(lz),wk(lvel),wk(lt1wk),
     .                         wk(lt2wk),wk(lt3wk),nbll)
c
 1820             continue
c
                  call lead(nbl,lw,lw2,maxbl)
               end if
c
               end if
c
            end if
         end if
c
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,31)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
c
  200    continue
      end if
c
      return
      end
      subroutine elglobfe(sa,xst,yst,zst,xix,xiy,xiz,etax,etay,
     .                    etaz,zetax,zetay,zetaz,ei,ej,ek,gij,gjk,
     .                    gik,ooj,eps,stiffl,islavept,ija,nslave,
     .                    nmaster,nnodes,myhost,myid,mycomm,mblk2nd,
     .                    nblelst,maxbl,iseqr)
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
       dimension xst(nslave),yst(nslave),zst(nslave),xix(nslave),
     .           xiy(nslave),xiz(nslave),etax(nslave),etay(nslave),
     .           etaz(nslave),zetax(nslave),zetay(nslave),zetaz(nslave),
     .           ei(nslave),ej(nslave),ek(nslave),gij(nslave),
     .           gjk(nslave),gik(nslave),tinv(6,6),t(6,6),ul(6,6),
     .           c(6,6),cjp1(6,6),cjm1(6,6),cim1(6,6),cip1(6,6),
     .           ckm1(6,6),ckp1(6,6),ooj(nslave)
      dimension sa(245*nslave+2),stiffl(24,24)
      dimension ija(245*nslave+2)
      dimension islavept(nslave,nmaster,5)
      dimension n1(0:7),n2(0:7),n11(0:7),n11i(20,0:7),n33i(20,0:7),
     .          iimax(0:7),n44i(0:7)
      dimension mblk2nd(maxbl)
      dimension nblelst(maxbl,2)
      common /twod/ i2d
c
c
c
      isl1 = 8
      isl2 = 2
      isl3 = 3
      isl4 = 4
      isl5 = 5
      if(i2d.ne.0) then
        isl6 = 4
        isl7 = 5
      else
        isl6 = 6
        isl7 = 7
      end if
      n2(0) = 1
      n2(1) = 4
      n2(2) = 16
      n2(3) = 13
      n2(4) = 7
      n2(5) = 19
      n2(6) = 10
      n2(7) = 22
c
       do 200 n = 1,nslave
c
        n1(0)    = n
        n1(1)    = islavept(n,3,iseqr)
        n1(2)    = islavept(n,5,iseqr)
        n1(3)    = islavept(n1(2),3,iseqr)
        n1(4)    = islavept(n,7,iseqr)
        n1(5)    = islavept(n1(2),7,iseqr)
        n1(6)    = islavept(n1(1),7,iseqr)
        n1(7)    = islavept(n1(6),5,iseqr)
        nbl = islavept(n,9,iseqr)
        if(n1(1).ne.n.and.n1(2).ne.n.and.n1(4).ne.n.and.
     .     islavept(n1(1),9,iseqr).eq.nbl.and.
     .     islavept(n1(2),9,iseqr).eq.nbl.and.
     .     islavept(n1(4),9,iseqr).eq.nbl) then
c
        xi =.25*(xst(n1(4))-xst(n1(0))+xst(n1(5))-xst(n1(2))
     .          +xst(n1(6))-xst(n1(1))+xst(n1(7))-xst(n1(3)))
        zi =.25*(zst(n1(4))-zst(n1(0))+zst(n1(5))-zst(n1(2))
     .          +zst(n1(6))-zst(n1(1))+zst(n1(7))-zst(n1(3)))
        yi =.25*(yst(n1(4))-yst(n1(0))+yst(n1(5))-yst(n1(2))
     .          +yst(n1(6))-yst(n1(1))+yst(n1(7))-yst(n1(3)))
        xj =.25*(xst(n1(1))-xst(n1(0))+xst(n1(3))-xst(n1(2))
     .          +xst(n1(6))-xst(n1(4))+xst(n1(7))-xst(n1(5)))
        yj =.25*(yst(n1(1))-yst(n1(0))+yst(n1(3))-yst(n1(2))
     .          +yst(n1(6))-yst(n1(4))+yst(n1(7))-yst(n1(5)))
        zj =.25*(zst(n1(1))-zst(n1(0))+zst(n1(3))-zst(n1(2))
     .          +zst(n1(6))-zst(n1(4))+zst(n1(7))-zst(n1(5)))
        denomi= sqrt(xi*xi+yi*yi+zi*zi)
        xi = xi/denomi
        yi = yi/denomi
        zi = zi/denomi
        denomj= sqrt(xj*xj+yj*yj+zj*zj)
        xj = xj/denomj
        yj = yj/denomj
        zj = zj/denomj
        xk = yi*zj - zi*yj
        yk = zi*xj - xi*zj
        zk = xi*yj - yi*xj
        xi = yj*zk - zj*yk
        yi = zj*xk - xj*zk
        zi = xj*yk - yj*xk
c
c  Note the row-column ordering here because
c  the transformation is from local coord to
c  the x,y,z - coordinate system
c
        t11 = xi
        t21 = yi
        t31 = zi
        t12 = xj
        t22 = yj
        t32 = zj
        t13 = xk
        t23 = yk
        t33 = zk
        t(1,1) = t11*t11
        t(1,2) = t12*t12
        t(1,3) = t13*t13
        t(1,4) = 2.*t11*t12
        t(1,5) = 2.*t12*t13
        t(1,6) = 2.*t11*t13
        t(2,1) = t21*t21
        t(2,2) = t22*t22
        t(2,3) = t23*t23
        t(2,4) = 2.*t21*t22
        t(2,5) = 2.*t22*t23
        t(2,6) = 2.*t21*t23
        t(3,1) = t31*t31
        t(3,2) = t32*t32
        t(3,3) = t33*t33
        t(3,4) = 2.*t31*t32
        t(3,5) = 2.*t32*t33
        t(3,6) = 2.*t31*t33
        t(4,1) = t11*t21
        t(4,2) = t12*t22
        t(4,3) = t13*t23
        t(4,4) = t11*t22+t21*t12
        t(4,5) = t12*t23+t13*t22
        t(4,6) = t11*t23+t13*t21
        t(5,1) = t31*t21
        t(5,2) = t22*t32
        t(5,3) = t23*t33
        t(5,4) = t21*t32+t22*t31
        t(5,5) = t22*t33+t32*t23
        t(5,6) = t21*t33+t23*t31
        t(6,1) = t11*t31
        t(6,2) = t12*t32
        t(6,3) = t13*t33
        t(6,4) = t11*t32+t12*t31
        t(6,5) = t12*t33+t13*t32
        t(6,6) = t11*t33+t31*t13
 
c
c  The inverse must be used:
c
        call INVDET(ul,t,6,DTNRM,DETM)
c
        do jj = 1,6
         tinv(1,jj) =  ei(n) *ul(1,jj)
         tinv(2,jj) =  ej(n) *ul(2,jj)
         tinv(3,jj) =  ek(n) *ul(3,jj)
         tinv(4,jj) =  gij(n)*ul(4,jj)
         tinv(5,jj) =  gjk(n)*ul(5,jj)
         tinv(6,jj) =  gik(n)*ul(6,jj)
        enddo
        do jj = 1,6
          do ii = 1,6
            sum = 0.
            do kk = 1,6
              sum = sum + t(ii,kk)*tinv(kk,jj)
            enddo
            c(ii,jj) = sum
          enddo
        enddo
c
      c0 = 1./2.
      c1 = 1./6.
      c2 = 1./9.
      c3 = 1./12.
      c4 = 1./18.
      c5 = 1./24.
      c6 = 1./36.
c
c    Rows 1-3
c
       term1 = xix(n)  *xix(n)  +etax(n) *etax(n) +zetax(n)*zetax(n)
       term2 = etax(n) *xix(n)  +zetax(n)*xix(n)  +zetax(n)* etax(n)
       term3 = xiy(n)  *xiy(n)  +etay(n) *etay(n) +zetay(n)*zetay(n)
       term4 = etay(n) *xiy(n)  +zetay(n)*xiy(n)  +zetay(n)* etay(n)
       term5 = xiz(n)  *xiz(n)  +etaz(n) *etaz(n) +zetaz(n)*zetaz(n)
       term6 = etaz(n) *xiz(n)  +zetaz(n)*xiz(n)  +zetaz(n)* etaz(n)
       term7 = xix(n)  *xiy(n)  +etay(n) *etax(n) +zetay(n)*zetax(n)
       term8 = xix(n)  *etay(n) +xix(n)  *zetay(n)+xiy(n)  *etax(n)
     .        +zetay(n)*etax(n) +xiy(n)  *zetax(n)+etay(n) *zetax(n)
       term9 = xiz(n)  *xix(n)  +etaz(n) *etax(n) +zetaz(n)*zetax(n)
       term10= xix(n)  *etaz(n) +xix(n)  *zetaz(n)+xiz(n)  *etax(n)
     .        +zetaz(n)*etax(n) +xiz(n)  *zetax(n)+etaz(n) *zetax(n)
       term11= xiz(n)  *xiy(n)  +etaz(n) *etay(n) +zetaz(n)*zetay(n)
       term12= xiy(n)  *etaz(n) +xiy(n)  *zetaz(n)+xiz(n)  *etay(n)
     .        +zetaz(n)*etay(n) +xiz(n)  *zetay(n)+etaz(n) *zetay(n)
 
 
       term13 = c2*term1+c1*term2
       term14 = c2*term3+c1*term4
       term15 = c2*term5+c1*term6
       term16 = c2*term7 +c3*term8
       term17 = c2*term9 +c3*term10
       term18 = c2*term11+c3*term12
       stiffl(1,1)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16
     .             +(c(1,6)+c(6,1))*term17
     .             +(c(4,6)+c(6,4))*term18
       stiffl(1,2)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(1,2)+c(4,4))*term16
     .             +(c(1,5)+c(6,4))*term17
     .             +(c(4,5)+c(6,2))*term18
       stiffl(1,3)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(1,5)+c(4,6))*term16
     .             +(c(1,3)+c(6,6))*term17
     .             +(c(4,3)+c(6,5))*term18
       stiffl(2,2)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16
     .             +(c(2,5)+c(5,2))*term18
     .             +(c(4,5)+c(5,4))*term17
       stiffl(2,3)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16
     .             +(c(2,3)+c(5,5))*term18
     .             +(c(4,3)+c(5,6))*term17
       stiffl(3,3)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(6,5)+c(5,6))*term16
     .             +(c(3,5)+c(5,3))*term18
     .             +(c(3,6)+c(6,3))*term17
  
       term13= -c2*xix(n)  *  xix(n)+c4* etax(n) *etax(n)
     .         +c3*zetax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14= -c2*xiy(n)  *  xiy(n)+c4* etay(n) *etay(n)
     .         +c3*zetay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15= -c2*xiz(n)  *  xiz(n)+c4* etaz(n) *etaz(n)
     .         +c3*zetaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16= -c2*xiy(n)  *  xix(n)+c4* etay(n) *etax(n)
     .         +c5*zetay(n)* etax(n)+c5* etay(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*( etax(n)*  xiy(n)+   zetax(n)*  xiy(n)
     .            -  xix(n)* etay(n)-     xix(n)*zetay(n))
       term18= -c2*xiz(n)  *  xix(n)+c4* etaz(n) *etax(n)
     .         +c5*zetaz(n)* etax(n)+c5* etaz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*( etax(n)*  xiz(n)+   zetax(n)*  xiz(n)
     .            -  xix(n)* etaz(n)-     xix(n)*zetaz(n))
       term20= -c2*xiz(n)  *  xiy(n)+c4* etaz(n) *etay(n)
     .         +c5*zetaz(n)* etay(n)+c5* etaz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*( etay(n)*  xiz(n)+   zetay(n)*  xiz(n)
     .            -  xiy(n)* etaz(n)-     xiy(n)*zetaz(n))
       stiffl(1,4)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(4,1)-c(1,4))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(6,1)-c(1,6))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(6,4)-c(4,6))*term21
       stiffl(2,4)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(2,1)-c(4,4))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(5,1)-c(4,6))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(5,4)-c(2,6))*term21
       stiffl(3,4)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(5,1)-c(6,4))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(3,1)-c(6,6))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(3,4)-c(5,6))*term21
       stiffl(3,5)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(5,4)-c(6,2))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(3,4)-c(6,5))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(3,2)-c(5,5))*term21
       stiffl(2,5)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(2,4)-c(4,2))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(5,4)-c(4,5))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(5,2)-c(2,5))*term21
       stiffl(1,5)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(4,4)-c(1,2))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(6,4)-c(1,5))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(6,2)-c(4,5))*term21
       stiffl(3,6)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(5,6)-c(6,5))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(3,6)-c(6,3))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(3,5)-c(5,3))*term21
       stiffl(2,6)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(2,6)-c(4,5))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(5,6)-c(4,3))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(5,5)-c(2,3))*term21
       stiffl(1,6)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(4,6)-c(1,5))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(6,6)-c(1,3))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(6,5)-c(4,3))*term21
  
       term13=  c4*xix(n)  *  xix(n)+c4* etax(n) *etax(n)
     .         +c3*  xix(n)* etax(n)-c2*zetax(n)*zetax(n)
       term14=  c4*xiy(n)  *  xiy(n)+c4* etay(n) *etay(n)
     .         +c3*  xiy(n)* etay(n)-c2*zetay(n)*zetay(n)
       term15=  c4*xiz(n)  *  xiz(n)+c4* etaz(n) *etaz(n)
     .         +c3*  xiz(n)* etaz(n)-c2*zetaz(n)*zetaz(n)
       term16=  c4*xiy(n)  *  xix(n)+c5* etay(n) * xix(n)
     .         +c5*  xiy(n)* etax(n)+c4* etay(n)* etax(n)
     .         -c2*zetay(n)*zetax(n)
       term17= c3*(zetax(n)*  xiy(n)+   zetax(n)* etay(n)
     .            -  xix(n)*zetay(n)-    etax(n)*zetay(n))
       term18=  c4*xiz(n)  *  xix(n)+c5* etaz(n) * xix(n)
     .         +c5*  xiz(n)* etax(n)+c4* etaz(n)* etax(n)
     .         -c2*zetaz(n)*zetax(n)
       term19= c3*(zetax(n)*  xiz(n)+   zetax(n)* etaz(n)
     .            -  xix(n)*zetaz(n)-    etax(n)*zetaz(n))
       term20=  c4*xiz(n)  *  xiy(n)+c5* etaz(n) * xiy(n)
     .         +c5*  xiz(n)* etay(n)+c4* etaz(n)* etay(n)
     .         -c2*zetaz(n)*zetay(n)
       term21= c3*(zetay(n)*  xiz(n)+   zetay(n)* etaz(n)
     .            -  xiy(n)*zetaz(n)-    etay(n)*zetaz(n))
       stiffl(1,7)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(2,7)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(3,7)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(3,8)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(2,8)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(1,8)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(3,9)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(2,9)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(1,9)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
       term13= -c4*xix(n)  *  xix(n)+c6* etax(n) *etax(n)
     .         -c1*  xix(n)*zetax(n)-c4*zetax(n)*zetax(n)
       term14= -c4*xiy(n)  *  xiy(n)+c6* etay(n) *etay(n)
     .         -c1*  xiy(n)*zetay(n)-c4*zetay(n)*zetay(n)
       term15= -c4*xiz(n)  *  xiz(n)+c6* etaz(n) *etaz(n)
     .         -c1*  xiz(n)*zetaz(n)-c4*zetaz(n)*zetaz(n)
       term16= -c4*xiy(n)  *  xix(n)-c3*zetay(n) * xix(n)
     .         +c6* etay(n)* etax(n)-c3*  xiy(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*(-etax(n)*  xiy(n)+     xix(n)* etay(n)
     .            +zetax(n)* etay(n)-    etax(n)*zetay(n))
       term18= -c4*xiz(n)  *  xix(n)-c3*zetaz(n) * xix(n)
     .         +c6* etaz(n)* etax(n)-c3*  xiz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*(-etax(n)*  xiz(n)+     xix(n)* etaz(n)
     .            +zetax(n)* etaz(n)-    etax(n)*zetaz(n))
       term20= -c4*xiz(n)  *  xiy(n)-c3*zetaz(n) * xiy(n)
     .         +c6* etaz(n)* etay(n)-c3*  xiz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*(-etay(n)*  xiz(n)+     xiy(n)* etaz(n)
     .            +zetay(n)* etaz(n)-    etay(n)*zetaz(n))
       stiffl(1,10)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(2,10)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(3,10)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(3,11)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(2,11)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(1,11)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(3,12)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(2,12)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(1,12)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*xix(n)  *  xix(n)-c4* etax(n) *etax(n)
     .         -c1*  xix(n)* etax(n)+c6*zetax(n)*zetax(n)
       term14= -c4*xiy(n)  *  xiy(n)-c4* etay(n) *etay(n)
     .         -c1*  xiy(n)* etay(n)+c6*zetay(n)*zetay(n)
       term15= -c4*xiz(n)  *  xiz(n)-c4* etaz(n) *etaz(n)
     .         -c1*  xiz(n)* etaz(n)+c6*zetaz(n)*zetaz(n)
       term16= -c4*xiy(n)  *  xix(n)-c3* etay(n) * xix(n)
     .         -c3*  xiy(n)* etax(n)-c4* etay(n)* etax(n)
     .         +c6*zetay(n)*zetax(n)
       term17= c5*(-zetax(n)*  xiy(n) -  zetax(n)* etay(n)
     .             +zetay(n)*  xix(n) +   etax(n)*zetay(n))
       term18= -c4*xiz(n)  *  xix(n)-c3* etaz(n) * xix(n)
     .         -c3*  xiz(n)* etax(n)-c4* etaz(n)* etax(n)
     .         +c6*zetaz(n)*zetax(n)
       term19= c5*(-zetax(n)*  xiz(n) -  zetax(n)* etaz(n)
     .             +zetaz(n)*  xix(n) +   etax(n)*zetaz(n))
       term20= -c4*xiz(n)  *  xiy(n)-c3* etaz(n) * xiy(n)
     .         -c3*  xiz(n)* etay(n)-c4* etaz(n)* etay(n)
     .         +c6*zetaz(n)*zetay(n)
       term21= c5*(-zetay(n)*  xiz(n) -  zetay(n)* etaz(n)
     .             +zetaz(n)*  xiy(n) +   etay(n)*zetaz(n))
       stiffl(1,13)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(2,13)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(3,13)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(3,14)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(2,14)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(1,14)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(3,15)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(2,15)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(1,15)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*xix(n)  *  xix(n)-c2* etax(n) *etax(n)
     .         +c3*  xix(n)*zetax(n)+c4*zetax(n)*zetax(n)
       term14=  c4*xiy(n)  *  xiy(n)-c2* etay(n) *etay(n)
     .         +c3*  xiy(n)*zetay(n)+c4*zetay(n)*zetay(n)
       term15=  c4*xiz(n)  *  xiz(n)-c2* etaz(n) *etaz(n)
     .         +c3*  xiz(n)*zetaz(n)+c4*zetaz(n)*zetaz(n)
       term16=  c4*xiy(n)  *  xix(n)+c5*zetay(n) * xix(n)
     .         +c5*  xiy(n)*zetax(n)-c2* etay(n)* etax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*(-etay(n)*  xix(n) +   etax(n)*  xiy(n)
     .             +etax(n)*zetay(n) -  zetax(n)* etay(n))
       term18=  c4*xiz(n)  *  xix(n)+c5*zetaz(n) * xix(n)
     .         +c5*  xiz(n)*zetax(n)-c2* etaz(n)* etax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*(-etaz(n)*  xix(n) +   etax(n)*  xiz(n)
     .             +etax(n)*zetaz(n) -  zetax(n)* etaz(n))
       term20=  c4*xiz(n)  *  xiy(n)+c5*zetaz(n) * xiy(n)
     .         +c5*  xiz(n)*zetay(n)-c2* etaz(n)* etay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*(-etaz(n)*  xiy(n) +   etay(n)*  xiz(n)
     .             +etay(n)*zetaz(n) -  zetay(n)* etaz(n))
       stiffl(1,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(2,16)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(3,16)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(3,17)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(2,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(1,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(3,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(2,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(1,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c6*xix(n)  *  xix(n)-c4* etax(n) *etax(n)
     .         -c1*zetax(n)* etax(n)-c4*zetax(n)*zetax(n)
       term14=  c6*xiy(n)  *  xiy(n)-c4* etay(n) *etay(n)
     .         -c1*zetay(n)* etay(n)-c4*zetay(n)*zetay(n)
       term15=  c6*xiz(n)  *  xiz(n)-c4* etaz(n) *etaz(n)
     .         -c1*zetaz(n)* etaz(n)-c4*zetaz(n)*zetaz(n)
       term16=  c6*xiy(n)  *  xix(n)-c4* etay(n) *etax(n)
     .         -c3*zetay(n)* etax(n)-c3* etay(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*( etax(n)*  xiy(n) +  zetax(n)*  xiy(n)
     .            -  xix(n)* etay(n) -    xix(n)*zetay(n))
       term18=  c6*xiz(n)  *  xix(n)-c4* etaz(n) *etax(n)
     .         -c3*zetaz(n)* etax(n)-c3* etaz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*( etax(n)*  xiz(n) +  zetax(n)*  xiz(n)
     .            -  xix(n)* etaz(n) -    xix(n)*zetaz(n))
       term20=  c6*xiz(n)  *  xiy(n)-c4* etaz(n) *etay(n)
     .         -c3*zetaz(n)* etay(n)-c3* etaz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*( etay(n)*  xiz(n) +  zetay(n)*  xiz(n)
     .            -  xiy(n)* etaz(n) -    xiy(n)*zetaz(n))
       stiffl(1,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(2,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(3,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(3,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(2,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(1,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(3,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(2,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(1,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c6*xix(n)  *  xix(n)-c3* etax(n) * xix(n)
     .         -c3*zetax(n)*  xix(n)-c6* etax(n)* etax(n)
     .         -c3*zetax(n)* etax(n)-c6*zetax(n)*zetax(n)
       term14= -c6*xiy(n)  *  xiy(n)-c3* etay(n) * xiy(n)
     .         -c3*zetay(n)*  xiy(n)-c6* etay(n)* etay(n)
     .         -c3*zetay(n)* etay(n)-c6*zetay(n)*zetay(n)
       term15= -c6*xiz(n)  *  xiz(n)-c3* etaz(n) * xiz(n)
     .         -c3*zetaz(n)*  xiz(n)-c6* etaz(n)* etaz(n)
     .         -c3*zetaz(n)* etaz(n)-c6*zetaz(n)*zetaz(n)
       term16= -c6*xiy(n)  *  xix(n)-c6* etay(n) *etax(n)
     .         -c6*zetay(n)*zetax(n)
       term17= c5*(-etay(n)*  xix(n) -  zetay(n)*  xix(n)
     .            -  xiy(n)* etax(n) -  zetay(n)* etax(n)
     .            -  xiy(n)*zetax(n) -   etay(n)*zetax(n))
       term18= -c6*xiz(n)  *  xix(n)-c6* etaz(n) *etax(n)
     .         -c6*zetaz(n)*zetax(n)
       term19= c5*(-etaz(n)*  xix(n) -  zetaz(n)*  xix(n)
     .            -  xiz(n)* etax(n) -  zetaz(n)* etax(n)
     .            -  xiz(n)*zetax(n) -   etaz(n)*zetax(n))
       term20= -c6*xiz(n)  *  xiy(n)-c6* etaz(n) *etay(n)
     .         -c6*zetaz(n)*zetay(n)
       term21= c5*(-etaz(n)*  xiy(n) -  zetaz(n)*  xiy(n)
     .            -  xiz(n)* etay(n) -  zetaz(n)* etay(n)
     .            -  xiz(n)*zetay(n) -   etaz(n)*zetay(n))
       stiffl(1,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)+c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)+c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)+c(6,4))*term21
       stiffl(2,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)+c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)+c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)+c(5,4))*term21
       stiffl(3,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)+c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)+c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)+c(3,4))*term21
       stiffl(3,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)+c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)+c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)+c(3,2))*term21
       stiffl(2,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)+c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)+c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)+c(5,2))*term21
       stiffl(1,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)+c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)+c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)+c(6,2))*term21
       stiffl(3,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)+c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)+c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)+c(3,5))*term21
       stiffl(2,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)+c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)+c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)+c(5,5))*term21
       stiffl(1,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)+c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)+c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)+c(6,5))*term21
c
c
c    Rows 4-6
c
       term1 = xix(n)  *xix(n)  +etax(n) *etax(n) +zetax(n)*zetax(n)
       term2 =-etax(n) *xix(n)  -zetax(n)*xix(n)  +zetax(n)* etax(n)
       term3 = xiy(n)  *xiy(n)  +etay(n) *etay(n) +zetay(n)*zetay(n)
       term4 =-etay(n) *xiy(n)  -zetay(n)*xiy(n)  +zetay(n)* etay(n)
       term5 = xiz(n)  *xiz(n)  +etaz(n) *etaz(n) +zetaz(n)*zetaz(n)
       term6 =-etaz(n) *xiz(n)  -zetaz(n)*xiz(n)  +zetaz(n)* etaz(n)
       term7 = xix(n)  *xiy(n)  +etay(n) *etax(n) +zetay(n)*zetax(n)
       term8 =-xix(n)  *etay(n) -xix(n)  *zetay(n)-xiy(n)  *etax(n)
     .        +zetay(n)*etax(n) -xiy(n)  *zetax(n)+etay(n) *zetax(n)
       term9 = xiz(n)  *xix(n)  +etaz(n) *etax(n) +zetaz(n)*zetax(n)
       term10=-xix(n)  *etaz(n) -xix(n)  *zetaz(n)-xiz(n)  *etax(n)
     .        +zetaz(n)*etax(n) -xiz(n)  *zetax(n)+etaz(n) *zetax(n)
       term11= xiz(n)  *xiy(n)  +etaz(n) *etay(n) +zetaz(n)*zetay(n)
       term12=-xiz(n)  *etay(n) -xiz(n)  *zetay(n)-xiy(n)  *etaz(n)
     .        +zetay(n)*etaz(n) -xiy(n)  *zetaz(n)+etay(n) *zetaz(n)
       term13 = c2*term1+c1*term2
       term14 = c2*term3+c1*term4
       term15 = c2*term5+c1*term6
       term16 = c2*term7 +c3*term8
       term17 = c2*term9 +c3*term10
       term18 = c2*term11+c3*term12
       stiffl(4,4)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16
     .             +(c(1,6)+c(6,1))*term17
     .             +(c(4,6)+c(6,4))*term18
       stiffl(4,5)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(1,2)+c(4,4))*term16
     .             +(c(1,5)+c(6,4))*term17
     .             +(c(4,5)+c(6,2))*term18
       stiffl(4,6)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(1,5)+c(4,6))*term16
     .             +(c(1,3)+c(6,6))*term17
     .             +(c(4,3)+c(6,5))*term18
       stiffl(5,5)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16
     .             +(c(4,5)+c(5,4))*term17
     .             +(c(2,5)+c(5,2))*term18
       stiffl(5,6)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16
     .             +(c(4,3)+c(5,6))*term17
     .             +(c(2,3)+c(5,5))*term18
       stiffl(6,6)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(6,5)+c(5,6))*term16
     .             +(c(3,6)+c(6,3))*term17
     .             +(c(3,5)+c(5,3))*term18
  
       term13= -c4*xix(n)  *  xix(n)+c6* etax(n) *etax(n)
     .         +c1*  xix(n)*zetax(n)-c4*zetax(n)*zetax(n)
       term14= -c4*xiy(n)  *  xiy(n)+c6* etay(n) *etay(n)
     .         +c1*  xiy(n)*zetay(n)-c4*zetay(n)*zetay(n)
       term15= -c4*xiz(n)  *  xiz(n)+c6* etaz(n) *etaz(n)
     .         +c1*  xiz(n)*zetaz(n)-c4*zetaz(n)*zetaz(n)
       term16= -c4*xiy(n)  *  xix(n)+c3*zetay(n) * xix(n)
     .         +c3*  xiy(n)*zetax(n)+c6* etay(n)* etax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*( etax(n)*  xiy(n)-     xix(n)* etay(n)
     .            +zetax(n)* etay(n)-    etax(n)*zetay(n))
       term18= -c4*xiz(n)  *  xix(n)+c3*zetaz(n) * xix(n)
     .         +c3*  xiz(n)*zetax(n)+c6* etaz(n)* etax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*( etax(n)*  xiz(n)-     xix(n)* etaz(n)
     .            +zetax(n)* etaz(n)-    etax(n)*zetaz(n))
       term20= -c4*xiz(n)  *  xiy(n)+c3*zetaz(n) * xiy(n)
     .         +c3*  xiz(n)*zetay(n)+c6* etaz(n)* etay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*( etay(n)*  xiz(n)-     xiy(n)* etaz(n)
     .            +zetay(n)* etaz(n)-    etay(n)*zetaz(n))
       stiffl(4,7)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(5,7)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(6,7)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(6,8)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(5,8)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(4,8)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(6,9)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(5,9)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(4,9)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
 
       term13=  c4*xix(n)  *  xix(n)+c4* etax(n) *etax(n)
     .         -c3*  xix(n)* etax(n)-c2*zetax(n)*zetax(n)
       term14=  c4*xiy(n)  *  xiy(n)+c4* etay(n) *etay(n)
     .         -c3*  xiy(n)* etay(n)-c2*zetay(n)*zetay(n)
       term15=  c4*xiz(n)  *  xiz(n)+c4* etaz(n) *etaz(n)
     .         -c3*  xiz(n)* etaz(n)-c2*zetaz(n)*zetaz(n)
       term16=  c4*xiy(n)  *  xix(n)-c5* etax(n) * xiy(n)
     .         +c4* etay(n)* etax(n)-c5*  xix(n)* etay(n)
     .         -c2*zetay(n)*zetax(n)
       term17= c3*(zetay(n)*  xix(n)-   zetay(n)* etax(n)
     .            -zetax(n)*  xiy(n)+    etay(n)*zetax(n))
       term18=  c4*xiz(n)  *  xix(n)-c5* etax(n) * xiz(n)
     .         +c4* etaz(n)* etax(n)-c5*  xix(n)* etaz(n)
     .         -c2*zetaz(n)*zetax(n)
       term19= c3*(zetaz(n)*  xix(n)-   zetaz(n)* etax(n)
     .            -zetax(n)*  xiz(n)+    etaz(n)*zetax(n))
       term20=  c4*xiz(n)  *  xiy(n)-c5* etay(n) * xiz(n)
     .         +c4* etaz(n)* etay(n)-c5*  xiy(n)* etaz(n)
     .         -c2*zetaz(n)*zetay(n)
       term21= c3*(zetaz(n)*  xiy(n)-   zetaz(n)* etay(n)
     .            -zetay(n)*  xiz(n)+    etaz(n)*zetay(n))
       stiffl(4,10)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(5,10)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(6,10)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(6,11)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(5,11)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(4,11)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(6,12)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(5,12)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(4,12)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*xix(n)  *  xix(n)-c3*zetax(n) * xix(n)
     .         -c2* etax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14=  c4*xiy(n)  *  xiy(n)-c3*zetay(n) * xiy(n)
     .         -c2* etay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15=  c4*xiz(n)  *  xiz(n)-c3*zetaz(n) * xiz(n)
     .         -c2* etaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16=  c4*xiy(n)  *  xix(n)-c5*zetay(n) * xix(n)
     .         -c2* etay(n)* etax(n)-c5*  xiy(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*(- etax(n)*  xiy(n) +    xix(n)* etay(n)
     .             -zetax(n)* etay(n) +   etax(n)*zetay(n))
       term18=  c4*xiz(n)  *  xix(n)-c5*zetaz(n) * xix(n)
     .         -c2* etaz(n)* etax(n)-c5*  xiz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*(- etax(n)*  xiz(n) +    xix(n)* etaz(n)
     .             -zetax(n)* etaz(n) +   etax(n)*zetaz(n))
       term20=  c4*xiz(n)  *  xiy(n)-c5*zetaz(n) * xiy(n)
     .         -c2* etaz(n)* etay(n)-c5*  xiz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*(- etay(n)*  xiz(n) +    xiy(n)* etaz(n)
     .             -zetay(n)* etaz(n) +   etay(n)*zetaz(n))
       stiffl(4,13)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(5,13)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(6,13)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(6,14)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(5,14)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(4,14)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(6,15)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(5,15)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(4,15)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*xix(n)  *  xix(n)+c1* etax(n) * xix(n)
     .         -c4* etax(n)* etax(n)+c6*zetax(n)*zetax(n)
       term14= -c4*xiy(n)  *  xiy(n)+c1* etay(n) * xiy(n)
     .         -c4* etay(n)* etay(n)+c6*zetay(n)*zetay(n)
       term15= -c4*xiz(n)  *  xiz(n)+c1* etaz(n) * xiz(n)
     .         -c4* etaz(n)* etaz(n)+c6*zetaz(n)*zetaz(n)
       term16= -c4*xiy(n)  *  xix(n)+c3* etay(n) * xix(n)
     .         +c3*  xiy(n)* etax(n)-c4* etay(n)* etax(n)
     .         +c6*zetay(n)*zetax(n)
       term17= c5*(zetax(n)*  xiy(n) -  zetax(n)* etay(n)
     .             - xix(n)*zetay(n) +   etax(n)*zetay(n))
       term18= -c4*xiz(n)  *  xix(n)+c3* etaz(n) * xix(n)
     .         +c3*  xiz(n)* etax(n)-c4* etaz(n)* etax(n)
     .         +c6*zetaz(n)*zetax(n)
       term19= c5*(zetax(n)*  xiz(n) -  zetax(n)* etaz(n)
     .             - xix(n)*zetaz(n) +   etax(n)*zetaz(n))
       term20= -c4*xiz(n)  *  xiy(n)+c3* etaz(n) * xiy(n)
     .         +c3*  xiz(n)* etay(n)-c4* etaz(n)* etay(n)
     .         +c6*zetaz(n)*zetay(n)
       term21= c5*(zetay(n)*  xiz(n) -  zetay(n)* etaz(n)
     .             - xiy(n)*zetaz(n) +   etay(n)*zetaz(n))
       stiffl(4,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(5,16)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(6,16)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(6,17)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(5,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(4,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(6,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(5,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(4,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c6*(xix(n)  *  xix(n)+ etax(n) *etax(n)
     .           +zetax(n)  *zetax(n))
     .        + c3*(etax(n) *  xix(n)+zetax(n) * xix(n)
     .        -     zetax(n)* etax(n))
       term14= -c6*(xiy(n)  *  xiy(n)+ etay(n) *etay(n)
     .           +zetay(n)  *zetay(n))
     .        + c3*(etay(n) *  xiy(n)+zetay(n) * xiy(n)
     .        -     zetay(n)* etay(n))
       term15= -c6*(xiz(n)  *  xiz(n)+ etaz(n) *etaz(n)
     .           +zetaz(n)  *zetaz(n))
     .        + c3*(etaz(n) *  xiz(n)+zetaz(n) * xiz(n)
     .        -     zetaz(n)* etaz(n))
       term16= -c6*xiy(n)  *  xix(n)-c4* etay(n) *etax(n)
     .         -c6*zetay(n)*zetax(n)+c5* etay(n)*  xix(n)
     .         +c5*zetay(n)*  xix(n)
       term17= c5*( etax(n)*  xiy(n) -   etax(n)*zetay(n)
     .            +  xiy(n)*zetax(n) -   etay(n)*zetax(n))
       term18= -c6*xiz(n)  *  xix(n)-c4* etaz(n) *etax(n)
     .         -c6*zetaz(n)*zetax(n)+c5* etaz(n)*  xix(n)
     .         +c5*zetaz(n)*  xix(n)
       term19= c5*( etax(n)*  xiz(n) -   etax(n)*zetaz(n)
     .            +  xiz(n)*zetax(n) -   etaz(n)*zetax(n))
       term20= -c6*xiz(n)  *  xiy(n)-c4* etaz(n) *etay(n)
     .         -c6*zetaz(n)*zetay(n)+c5* etay(n)*  xiz(n)
     .         +c5*zetay(n)*  xiz(n)
       term21= c5*( etaz(n)*  xiy(n) -   etaz(n)*zetay(n)
     .            +  xiy(n)*zetaz(n) -   etay(n)*zetaz(n))
       stiffl(4,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)+c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)+c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)+c(6,4))*term21
       stiffl(5,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)+c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)+c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)+c(5,4))*term21
       stiffl(6,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)+c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)+c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)+c(3,4))*term21
       stiffl(6,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)+c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)+c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)+c(3,2))*term21
       stiffl(5,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)+c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)+c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)+c(5,2))*term21
       stiffl(4,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)+c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)+c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)+c(6,2))*term21
       stiffl(6,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)+c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)+c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)+c(3,5))*term21
       stiffl(5,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)+c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)+c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)+c(5,5))*term21
       stiffl(4,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)+c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)+c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)+c(6,5))*term21
c
       term13=  c6*xix(n)  *  xix(n)-c4* etax(n) *etax(n)
     .         -c1*zetax(n)* etax(n)-c4*zetax(n)*zetax(n)
       term14=  c6*xiy(n)  *  xiy(n)-c4* etay(n) *etay(n)
     .         -c1*zetay(n)* etay(n)-c4*zetay(n)*zetay(n)
       term15=  c6*xiz(n)  *  xiz(n)-c4* etaz(n) *etaz(n)
     .         -c1*zetaz(n)* etaz(n)-c4*zetaz(n)*zetaz(n)
       term16=  c6*xiy(n)  *  xix(n)-c4* etay(n) *etax(n)
     .         -c3*zetay(n)* etax(n)-c3* etay(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*(-etax(n)*  xiy(n) -  zetax(n)*  xiy(n)
     .            +  xix(n)* etay(n) +  zetay(n)*  xix(n))
       term18=  c6*xiz(n)  *  xix(n)-c4* etaz(n) *etax(n)
     .         -c3*zetaz(n)* etax(n)-c3* etaz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*(-etax(n)*  xiz(n) -  zetax(n)*  xiz(n)
     .            +  xix(n)* etaz(n) +  zetaz(n)*  xix(n))
       term20=  c6*xiz(n)  *  xiy(n)-c4* etaz(n) *etay(n)
     .         -c3*zetaz(n)* etay(n)-c3* etaz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*(-etay(n)*  xiz(n) -  zetay(n)*  xiz(n)
     .            +  xiy(n)* etaz(n) +  zetaz(n)*  xiy(n))
       stiffl(4,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(5,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(6,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(6,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(5,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(4,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(6,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(5,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(4,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
c    Rows 7-9
c
       term13=  c2*xix(n)  *  xix(n)+c2* etax(n) *etax(n)
     .                              +c2*zetax(n)*zetax(n)
     .         +c1* etax(n)*  xix(n)-c1*zetax(n) * xix(n)
     .         -c1*zetax(n)* etax(n)
       term14=  c2*xiy(n)  *  xiy(n)+c2* etay(n) *etay(n)
     .                              +c2*zetay(n)*zetay(n)
     .         +c1* etay(n)*  xiy(n)-c1*zetay(n) * xiy(n)
     .         -c1*zetay(n)* etay(n)
       term15=  c2*xiz(n)  *  xiz(n)+c2* etaz(n) *etaz(n)
     .                              +c2*zetaz(n)*zetaz(n)
     .         +c1* etaz(n)*  xiz(n)-c1*zetaz(n) * xiz(n)
     .         -c1*zetaz(n)* etaz(n)
       term16=  c2*xiy(n)  *  xix(n)+c2* etay(n)* etax(n)
     .         +c2*zetay(n)*zetax(n)+c3* etay(n)*  xix(n)
     .         -c3*zetay(n)*  xix(n)+c3*  xiy(n)* etax(n)
     .         -c3*zetay(n)* etax(n)-c3*  xiy(n)*zetax(n)
     .         -c3* etay(n)*zetax(n)
       term18=  c2*xiz(n)  *  xix(n)+c2* etaz(n)* etax(n)
     .         +c2*zetaz(n)*zetax(n)+c3* etaz(n)*  xix(n)
     .         -c3*zetaz(n)*  xix(n)+c3*  xiz(n)* etax(n)
     .         -c3*zetaz(n)* etax(n)-c3*  xiz(n)*zetax(n)
     .         -c3* etaz(n)*zetax(n)
       term20=  c2*xiz(n)  *  xiy(n)+c2* etaz(n)* etay(n)
     .         +c2*zetaz(n)*zetay(n)+c3* etaz(n)*  xiy(n)
     .         -c3*zetaz(n)*  xiy(n)+c3*  xiz(n)* etay(n)
     .         -c3*zetaz(n)* etay(n)-c3*  xiz(n)*zetay(n)
     .         -c3* etaz(n)*zetay(n)
       stiffl(7,7)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16
     .             +(c(1,6)+c(6,1))*term18
     .             +(c(4,6)+c(6,4))*term20
       stiffl(8,8)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16
     .             +(c(5,4)+c(4,5))*term18
     .             +(c(5,2)+c(2,5))*term20
       stiffl(7,8)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16
     .             +(c(6,4)+c(1,5))*term18
     .             +(c(6,2)+c(4,5))*term20
       stiffl(9,9)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16
     .             +(c(3,6)+c(6,3))*term18
     .             +(c(3,5)+c(5,3))*term20
       stiffl(8,9)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16
     .             +(c(5,6)+c(4,3))*term18
     .             +(c(5,5)+c(2,3))*term20
       stiffl(7,9)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16
     .             +(c(6,6)+c(1,3))*term18
     .             +(c(6,5)+c(4,3))*term20
 
       term13= -c2*  xix(n)*  xix(n)+c4* etax(n) *etax(n)
     .         -c3*zetax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14= -c2*  xiy(n)*  xiy(n)+c4* etay(n) *etay(n)
     .         -c3*zetay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15= -c2*  xiz(n)*  xiz(n)+c4* etaz(n) *etaz(n)
     .         -c3*zetaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16= -c2*  xiy(n)*  xix(n)+c4* etay(n)* etax(n)
     .         -c5*zetay(n)* etax(n)-c5* etay(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= -c3* etax(n)*  xiy(n)+c3*zetax(n)*  xiy(n)
     .         +c3*  xix(n)* etay(n)-c3*  xix(n)*zetay(n)
       term18= -c2*  xiz(n)*  xix(n)+c4* etaz(n)* etax(n)
     .         -c5*zetaz(n)* etax(n)-c5* etaz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= -c3* etax(n)*  xiz(n)+c3*zetax(n)*  xiz(n)
     .         +c3*  xix(n)* etaz(n)-c3*  xix(n)*zetaz(n)
       term20= -c2*  xiz(n)*  xiy(n)+c4* etaz(n)* etay(n)
     .         -c5*zetaz(n)* etay(n)-c5* etaz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= -c3* etay(n)*  xiz(n)+c3*zetay(n)*  xiz(n)
     .         +c3*  xiy(n)* etaz(n)-c3*  xiy(n)*zetaz(n)
       stiffl(7,10)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(8,10)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(9,10)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(9,11)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(8,11)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(7,11)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(9,12)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(8,12)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(7,12)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c6*xix(n)  *  xix(n)-c6* etax(n) *etax(n)
     .                              -c6*zetax(n)*zetax(n)
     .         -c3* etax(n)*  xix(n)+c3*zetax(n) * xix(n)
     .         +c3*zetax(n)* etax(n)
       term14= -c6*xiy(n)  *  xiy(n)-c6* etay(n) *etay(n)
     .                              -c6*zetay(n)*zetay(n)
     .         -c3* etay(n)*  xiy(n)+c3*zetay(n) * xiy(n)
     .         +c3*zetay(n)* etay(n)
       term15= -c6*xiz(n)  *  xiz(n)-c6* etaz(n) *etaz(n)
     .                              -c6*zetaz(n)*zetaz(n)
     .         -c3* etaz(n)*  xiz(n)+c3*zetaz(n) * xiz(n)
     .         +c3*zetaz(n)* etaz(n)
       term16= -c6*xiy(n)  *  xix(n)-c6* etay(n)* etax(n)
     .         -c6*zetay(n)*zetax(n)-c5* etay(n)*  xix(n)
     .         +c5*zetay(n)*  xix(n)-c5*  xiy(n)* etax(n)
     .         +c5*zetay(n)* etax(n)+c5*  xiy(n)*zetax(n)
     .         +c5* etay(n)*zetax(n)
       term18= -c6*xiz(n)  *  xix(n)-c6* etaz(n)* etax(n)
     .         -c6*zetaz(n)*zetax(n)-c5* etaz(n)*  xix(n)
     .         +c5*zetaz(n)*  xix(n)-c5*  xiz(n)* etax(n)
     .         +c5*zetaz(n)* etax(n)+c5*  xiz(n)*zetax(n)
     .         +c5* etaz(n)*zetax(n)
       term20= -c6*xiz(n)  *  xiy(n)-c6* etaz(n)* etay(n)
     .         -c6*zetaz(n)*zetay(n)-c5* etaz(n)*  xiy(n)
     .         +c5*zetaz(n)*  xiy(n)-c5*  xiz(n)* etay(n)
     .         +c5*zetaz(n)* etay(n)+c5*  xiz(n)*zetay(n)
     .         +c5* etaz(n)*zetay(n)
       stiffl(7,13)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16
     .             +(c(1,6)+c(6,1))*term18
     .             +(c(4,6)+c(6,4))*term20
       stiffl(8,13)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16
     .             +(c(5,1)+c(4,6))*term18
     .             +(c(5,4)+c(2,6))*term20
       stiffl(9,13)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16
     .             +(c(3,1)+c(6,6))*term18
     .             +(c(3,4)+c(5,6))*term20
       stiffl(9,14)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16
     .             +(c(3,4)+c(6,5))*term18
     .             +(c(3,2)+c(5,5))*term20
       stiffl(8,14)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16
     .             +(c(5,4)+c(4,5))*term18
     .             +(c(5,2)+c(2,5))*term20
       stiffl(7,14)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16
     .             +(c(6,4)+c(1,5))*term18
     .             +(c(6,2)+c(4,5))*term20
       stiffl(9,15)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16
     .             +(c(3,6)+c(6,3))*term18
     .             +(c(3,5)+c(5,3))*term20
       stiffl(8,15)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16
     .             +(c(5,6)+c(4,3))*term18
     .             +(c(5,5)+c(2,3))*term20
       stiffl(7,15)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16
     .             +(c(6,6)+c(1,3))*term18
     .             +(c(6,5)+c(4,3))*term20
c
       term13=  c6*  xix(n)*  xix(n)-c4* etax(n) *etax(n)
     .         +c1*zetax(n)* etax(n)-c4*zetax(n)*zetax(n)
       term14=  c6*  xiy(n)*  xiy(n)-c4* etay(n) *etay(n)
     .         +c1*zetay(n)* etay(n)-c4*zetay(n)*zetay(n)
       term15=  c6*  xiz(n)*  xiz(n)-c4* etaz(n) *etaz(n)
     .         +c1*zetaz(n)* etaz(n)-c4*zetaz(n)*zetaz(n)
       term16=  c6*  xiy(n)*  xix(n)-c4* etay(n)* etax(n)
     .         +c3*zetay(n)* etax(n)+c3* etay(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= -c5* etay(n)*  xix(n)+c5*zetay(n)*  xix(n)
     .         +c5*  xiy(n)* etax(n)-c5*  xiy(n)*zetax(n)
       term18=  c6*  xiz(n)*  xix(n)-c4* etaz(n)* etax(n)
     .         +c3*zetaz(n)* etax(n)+c3* etaz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= -c5* etaz(n)*  xix(n)+c5*zetaz(n)*  xix(n)
     .         +c5*  xiz(n)* etax(n)-c5*  xiz(n)*zetax(n)
       term20=  c6*  xiz(n)*  xiy(n)-c4* etaz(n)* etay(n)
     .         +c3*zetaz(n)* etay(n)+c3* etaz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= -c5* etaz(n)*  xiy(n)+c5*zetaz(n)*  xiy(n)
     .         +c5*  xiz(n)* etay(n)-c5*  xiz(n)*zetay(n)
       stiffl(7,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(8,16)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(9,16)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(9,17)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(8,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(7,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(9,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(8,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(7,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*  xix(n)*  xix(n)-c3*zetax(n)*  xix(n)
     .         -c2* etax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14=  c4*  xiy(n)*  xiy(n)-c3*zetay(n)*  xiy(n)
     .         -c2* etay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15=  c4*  xiz(n)*  xiz(n)-c3*zetaz(n)*  xiz(n)
     .         -c2* etaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16=  c4*  xiy(n)*  xix(n)-c5*zetay(n)*  xix(n)
     .         -c2* etay(n)* etax(n)-c5*  xiy(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*(  xiy(n)* etax(n) -    xix(n)* etay(n)
     .            +zetax(n)* etay(n) -   etax(n)*zetay(n))
       term18=  c4*  xiz(n)*  xix(n)-c5*zetaz(n)*  xix(n)
     .         -c2* etaz(n)* etax(n)-c5*  xiz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*(  xiz(n)* etax(n) -    xix(n)* etaz(n)
     .            +zetax(n)* etaz(n) -   etax(n)*zetaz(n))
       term20=  c4*  xiz(n)*  xiy(n)-c5*zetaz(n)*  xiy(n)
     .         -c2* etaz(n)* etay(n)-c5*  xiz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*(  xiz(n)* etay(n) -    xiy(n)* etaz(n)
     .            +zetay(n)* etaz(n) -   etay(n)*zetaz(n))
       stiffl(7,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(8,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(9,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(9,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(8,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(7,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(9,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(8,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(7,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*  xix(n)*  xix(n)-c4* etax(n)* etax(n)
     .         -c1* etax(n)*  xix(n)+c6*zetax(n)*zetax(n)
       term14= -c4*  xiy(n)*  xiy(n)-c4* etay(n)* etay(n)
     .         -c1* etay(n)*  xiy(n)+c6*zetay(n)*zetay(n)
       term15= -c4*  xiz(n)*  xiz(n)-c4* etaz(n)* etaz(n)
     .         -c1* etaz(n)*  xiz(n)+c6*zetaz(n)*zetaz(n)
       term16= -c4*  xiy(n)*  xix(n)-c3* etay(n)*  xix(n)
     .         -c3*  xiy(n)* etax(n)-c4* etay(n)* etax(n)
     .         +c6*zetay(n)*zetax(n)
       term17= c5*(- xix(n)*zetay(n) -   etax(n)*zetay(n)
     .            +zetax(n)*  xiy(n) +  zetax(n)* etay(n))
       term18= -c4*  xiz(n)*  xix(n)-c3* etaz(n)*  xix(n)
     .         -c3*  xiz(n)* etax(n)-c4* etaz(n)* etax(n)
     .         +c6*zetaz(n)*zetax(n)
       term19= c5*(- xix(n)*zetaz(n) -   etax(n)*zetaz(n)
     .            +zetax(n)*  xiz(n) +  zetax(n)* etaz(n))
       term20= -c4*  xiz(n)*  xiy(n)-c3* etaz(n)*  xiy(n)
     .         -c3*  xiz(n)* etay(n)-c4* etaz(n)* etay(n)
     .         +c6*zetaz(n)*zetay(n)
       term21= c5*(- xiy(n)*zetaz(n) -   etay(n)*zetaz(n)
     .            +zetay(n)*  xiz(n) +  zetay(n)* etaz(n))
       stiffl(7,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(8,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(9,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(9,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(8,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(7,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(9,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(8,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(7,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c     Rows 10-12
c
       term13=  c2*  xix(n)*  xix(n)+c2* etax(n)* etax(n)
     .         -c1* etax(n)*  xix(n)+c1*zetax(n)*  xix(n)
     .         -c1*zetax(n)* etax(n)+c2*zetax(n)*zetax(n)
       term14=  c2*  xiy(n)*  xiy(n)+c2* etay(n)* etay(n)
     .         -c1* etay(n)*  xiy(n)+c1*zetay(n)*  xiy(n)
     .         -c1*zetay(n)* etay(n)+c2*zetay(n)*zetay(n)
       term15=  c2*  xiz(n)*  xiz(n)+c2* etaz(n)* etaz(n)
     .         -c1* etaz(n)*  xiz(n)+c1*zetaz(n)*  xiz(n)
     .         -c1*zetaz(n)* etaz(n)+c2*zetaz(n)*zetaz(n)
       term16=  c2*  xiy(n)*  xix(n)-c3* etay(n)*  xix(n)
     .         +c3*zetay(n)*  xix(n)-c3*  xiy(n)* etax(n)
     .         +c2* etay(n)* etax(n)-c3*zetay(n)* etax(n)
     .         +c3*  xiy(n)*zetax(n)-c3* etay(n)*zetax(n)
     .         +c2*zetay(n)*zetax(n)
       term17= 0.
       term18=  c2*  xiz(n)*  xix(n)-c3* etaz(n)*  xix(n)
     .         +c3*zetaz(n)*  xix(n)-c3*  xiz(n)* etax(n)
     .         +c2* etaz(n)* etax(n)-c3*zetaz(n)* etax(n)
     .         +c3*  xiz(n)*zetax(n)-c3* etaz(n)*zetax(n)
     .         +c2*zetaz(n)*zetax(n)
       term19= 0.
       term20=  c2*  xiz(n)*  xiy(n)-c3* etaz(n)*  xiy(n)
     .         +c3*zetaz(n)*  xiy(n)-c3*  xiz(n)* etay(n)
     .         +c2* etaz(n)* etay(n)-c3*zetaz(n)* etay(n)
     .         +c3*  xiz(n)*zetay(n)-c3* etaz(n)*zetay(n)
     .         +c2*zetaz(n)*zetay(n)
       term21= 0.
       stiffl(10,10)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(11,11)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(10,11)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(12,12)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(11,12)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(10,12)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c6*xix(n)  *  xix(n)-c4* etax(n) *etax(n)
     .                              -c4*zetax(n)*zetax(n)
     .         +c1*zetax(n)* etax(n)
       term14=  c6*xiy(n)  *  xiy(n)-c4* etay(n) *etay(n)
     .                              -c4*zetay(n)*zetay(n)
     .         +c1*zetay(n)* etay(n)
       term15=  c6*xiz(n)  *  xiz(n)-c4* etaz(n) *etaz(n)
     .                              -c4*zetaz(n)*zetaz(n)
     .         +c1*zetaz(n)* etaz(n)
       term16=  c6*  xiy(n)*  xix(n)-c4* etay(n)* etax(n)
     .         +c3*zetay(n)* etax(n)+c3* etay(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= -c5* etax(n)*  xiy(n)+c5*zetax(n)*  xiy(n)
     .         +c5*  xix(n)* etay(n)-c5*  xix(n)*zetay(n)
       term18=  c6*  xiz(n)*  xix(n)-c4* etaz(n)* etax(n)
     .         +c3*zetaz(n)* etax(n)+c3* etaz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= -c5* etax(n)*  xiz(n)+c5*zetax(n)*  xiz(n)
     .         +c5*  xix(n)* etaz(n)-c5*  xix(n)*zetaz(n)
       term20=  c6*  xiz(n)*  xiy(n)-c4* etaz(n)* etay(n)
     .         +c3*zetaz(n)* etay(n)+c3* etaz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= -c5* etay(n)*  xiz(n)+c5*zetay(n)*  xiz(n)
     .         +c5*  xiy(n)* etaz(n)-c5*  xiy(n)*zetaz(n)
       stiffl(10,13)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(11,13)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(4,4)+c(2,1))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(4,6)+c(5,1))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(2,6)+c(5,4))*term20+(c(2,6)-c(5,4))*term21
       stiffl(12,13)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(6,4)+c(5,1))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(6,6)+c(3,1))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(5,6)+c(3,4))*term20+(c(5,6)-c(3,4))*term21
       stiffl(12,14)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(6,2)+c(5,4))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(6,5)+c(3,4))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(5,5)+c(3,2))*term20+(c(5,5)-c(3,2))*term21
       stiffl(11,14)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(4,2)+c(2,4))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(4,5)+c(5,4))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(2,5)+c(5,2))*term20+(c(2,5)-c(5,2))*term21
       stiffl(10,14)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(1,2)+c(4,4))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(1,5)+c(6,4))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(4,5)+c(6,2))*term20+(c(4,5)-c(6,2))*term21
       stiffl(12,15)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(6,5)+c(5,6))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(6,3)+c(3,6))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(5,3)+c(3,5))*term20+(c(5,3)-c(3,5))*term20
       stiffl(11,15)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(4,5)+c(2,6))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(4,3)+c(5,6))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(2,3)+c(5,5))*term20+(c(2,3)-c(5,5))*term21
       stiffl(10,15)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(1,5)+c(4,6))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(1,3)+c(6,6))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(4,3)+c(6,5))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c6*  xix(n)*  xix(n)+c3* etax(n) * xix(n)
     .         -c3*zetax(n)*  xix(n)-c6* etax(n)* etax(n)
     .         +c3*zetax(n)* etax(n)-c6*zetax(n)*zetax(n)
       term14= -c6*  xiy(n)*  xiy(n)+c3* etay(n) * xiy(n)
     .         -c3*zetay(n)*  xiy(n)-c6* etay(n)* etay(n)
     .         +c3*zetay(n)* etay(n)-c6*zetay(n)*zetay(n)
       term15= -c6*  xiz(n)*  xiz(n)+c3* etaz(n) * xiz(n)
     .         -c3*zetaz(n)*  xiz(n)-c6* etaz(n)* etaz(n)
     .         +c3*zetaz(n)* etaz(n)-c6*zetaz(n)*zetaz(n)
       term16= -c6*  xiy(n)*  xix(n)+c5* etay(n)*  xix(n)
     .         -c5*zetay(n)*  xix(n)+c5*  xiy(n)* etax(n)
     .         -c6* etay(n)* etax(n)+c5*zetay(n)* etax(n)
     .         -c5*  xiy(n)*zetax(n)+c5* etay(n)*zetax(n)
     .         -c6*zetay(n)*zetax(n)
       term17= 0.
       term18= -c6*  xiz(n)*  xix(n)+c5* etaz(n)*  xix(n)
     .         -c5*zetaz(n)*  xix(n)+c5*  xiz(n)* etax(n)
     .         -c6* etaz(n)* etax(n)+c5*zetaz(n)* etax(n)
     .         -c5*  xiz(n)*zetax(n)+c5* etaz(n)*zetax(n)
     .         -c6*zetaz(n)*zetax(n)
       term19= 0.
       term20= -c6*  xiz(n)*  xiy(n)+c5* etaz(n)*  xiy(n)
     .         -c5*zetaz(n)*  xiy(n)+c5*  xiz(n)* etay(n)
     .         -c6* etaz(n)* etay(n)+c5*zetaz(n)* etay(n)
     .         -c5*  xiz(n)*zetay(n)+c5* etaz(n)*zetay(n)
     .         -c6*zetaz(n)*zetay(n)
       term21= 0.
       stiffl(10,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(11,16)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(12,16)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(12,17)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(11,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(10,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(12,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(11,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(10,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*  xix(n)*  xix(n)+c1* etax(n)*  xix(n)
     .         -c4* etax(n)* etax(n)+c6*zetax(n)*zetax(n)
       term14= -c4*  xiy(n)*  xiy(n)+c1* etay(n)*  xiy(n)
     .         -c4* etay(n)* etay(n)+c6*zetay(n)*zetay(n)
       term15= -c4*  xiz(n)*  xiz(n)+c1* etaz(n)*  xiz(n)
     .         -c4* etaz(n)* etaz(n)+c6*zetaz(n)*zetaz(n)
       term16= -c4*  xiy(n)*  xix(n)+c3* etay(n)*  xix(n)
     .         +c3*  xiy(n)* etax(n)-c4* etay(n)* etax(n)
     .         +c6*zetay(n)*zetax(n)
       term17= c5*(-zetax(n)*  xiy(n) + zetax(n)* etay(n)
     .             +  xix(n)*zetay(n) -  etax(n)*zetay(n))
       term18= -c4*  xiz(n)*  xix(n)+c3* etaz(n)*  xix(n)
     .         +c3*  xiz(n)* etax(n)-c4* etaz(n)* etax(n)
     .         +c6*zetaz(n)*zetax(n)
       term19= c5*(-zetax(n)*  xiz(n) + zetax(n)* etaz(n)
     .             +  xix(n)*zetaz(n) -  etax(n)*zetaz(n))
       term20= -c4*  xiz(n)*  xiy(n)+c3* etaz(n)*  xiy(n)
     .         +c3*  xiz(n)* etay(n)-c4* etaz(n)* etay(n)
     .         +c6*zetaz(n)*zetay(n)
       term21= c5*(-zetay(n)*  xiz(n) + zetay(n)* etaz(n)
     .             +  xiy(n)*zetaz(n) -  etay(n)*zetaz(n))
       stiffl(10,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(11,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(12,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(12,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(11,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(10,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(12,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(11,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(10,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*  xix(n)*  xix(n)+c3*zetax(n)*  xix(n)
     .         -c2* etax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14=  c4*  xiy(n)*  xiy(n)+c3*zetay(n)*  xiy(n)
     .         -c2* etay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15=  c4*  xiz(n)*  xiz(n)+c3*zetaz(n)*  xiz(n)
     .         -c2* etaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16=  c4*  xiy(n)*  xix(n)+c5*zetay(n)*  xix(n)
     .         -c2* etay(n)* etax(n)+c5*  xiy(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*(-etax(n)*  xiy(n) +    xix(n)* etay(n)
     .            +zetax(n)* etay(n) -   etax(n)*zetay(n))
       term18=  c4*  xiz(n)*  xix(n)+c5*zetaz(n)*  xix(n)
     .         -c2* etaz(n)* etax(n)+c5*  xiz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*(-etax(n)*  xiz(n) +    xix(n)* etaz(n)
     .            +zetax(n)* etaz(n) -   etax(n)*zetaz(n))
       term20=  c4*  xiz(n)*  xiy(n)+c5*zetaz(n)*  xiy(n)
     .         -c2* etaz(n)* etay(n)+c5*  xiz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*(-etay(n)*  xiz(n) +    xiy(n)* etaz(n)
     .            +zetay(n)* etaz(n) -   etay(n)*zetaz(n))
       stiffl(10,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(11,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(12,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(12,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(11,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(10,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(12,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(11,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(10,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
c     Rows 13-15
c
c
       term13=  c2*  xix(n)*  xix(n)+c2* etax(n) *etax(n)
     .                              +c2*zetax(n)*zetax(n)
     .         +c1* etax(n)*  xix(n)-c1*zetax(n)*  xix(n)
     .         -c1*zetax(n)* etax(n)
       term14=  c2*  xiy(n)*  xiy(n)+c2* etay(n) *etay(n)
     .                              +c2*zetay(n)*zetay(n)
     .         +c1* etay(n)*  xiy(n)-c1*zetay(n)*  xiy(n)
     .         -c1*zetay(n)* etay(n)
       term15=  c2*  xiz(n)*  xiz(n)+c2* etaz(n) *etaz(n)
     .                              +c2*zetaz(n)*zetaz(n)
     .         +c1* etaz(n)*  xiz(n)-c1*zetaz(n)*  xiz(n)
     .         -c1*zetaz(n)* etaz(n)
       term16=  c2*  xiy(n)*  xix(n)+c3* etay(n)*  xix(n)
     .         -c3*zetay(n)*  xix(n)+c3*  xiy(n)* etax(n)
     .         +c2* etay(n)* etax(n)-c3*zetay(n)* etax(n)
     .         -c3*  xiy(n)*zetax(n)-c3* etay(n)*zetax(n)
     .         +c2*zetay(n)*zetax(n)
       term17= 0.
       term18=  c2*  xiz(n)*  xix(n)+c3* etaz(n)*  xix(n)
     .         -c3*zetaz(n)*  xix(n)+c3*  xiz(n)* etax(n)
     .         +c2* etaz(n)* etax(n)-c3*zetaz(n)* etax(n)
     .         -c3*  xiz(n)*zetax(n)-c3* etaz(n)*zetax(n)
     .         +c2*zetaz(n)*zetax(n)
       term19= 0.
       term20=  c2*  xiz(n)*  xiy(n)+c3* etaz(n)*  xiy(n)
     .         -c3*zetaz(n)*  xiy(n)+c3*  xiz(n)* etay(n)
     .         +c2* etaz(n)* etay(n)-c3*zetaz(n)* etay(n)
     .         -c3*  xiz(n)*zetay(n)-c3* etaz(n)*zetay(n)
     .         +c2*zetaz(n)*zetay(n)
       term21= 0.
       stiffl(13,13)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(14,14)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(4,2)+c(2,4))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(4,5)+c(5,4))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(2,5)+c(5,2))*term20+(c(2,5)-c(5,2))*term21
       stiffl(13,14)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(1,2)+c(4,4))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(1,5)+c(6,4))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(4,5)+c(6,2))*term20+(c(4,5)-c(6,2))*term21
       stiffl(15,15)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(6,5)+c(5,6))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(6,3)+c(3,6))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(5,3)+c(3,5))*term20+(c(5,3)-c(3,5))*term20
       stiffl(14,15)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(4,5)+c(2,6))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(4,3)+c(5,6))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(2,3)+c(5,5))*term20+(c(2,3)-c(5,5))*term21
       stiffl(13,15)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(1,5)+c(4,6))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(1,3)+c(6,6))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(4,3)+c(6,5))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c2*  xix(n)*  xix(n)+c4* etax(n) *etax(n)
     .         -c3*zetax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14= -c2*  xiy(n)*  xiy(n)+c4* etay(n) *etay(n)
     .         -c3*zetay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15= -c2*  xiz(n)*  xiz(n)+c4* etaz(n) *etaz(n)
     .         -c3*zetaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16= -c2*  xiy(n)*  xix(n)+c4* etay(n)* etax(n)
     .         -c5*zetay(n)* etax(n)-c5* etay(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= -c3* etax(n)*  xiy(n)+c3*zetax(n)*  xiy(n)
     .         +c3*  xix(n)* etay(n)-c3*  xix(n)*zetay(n)
       term18= -c2*  xiz(n)*  xix(n)+c4* etaz(n)* etax(n)
     .         -c5*zetaz(n)* etax(n)-c5* etaz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= -c3* etax(n)*  xiz(n)+c3*zetax(n)*  xiz(n)
     .         +c3*  xix(n)* etaz(n)-c3*  xix(n)*zetaz(n)
       term20= -c2*  xiz(n)*  xiy(n)+c4* etaz(n)* etay(n)
     .         -c5*zetaz(n)* etay(n)-c5* etaz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= -c3* etay(n)*  xiz(n)+c3*zetay(n)*  xiz(n)
     .         +c3*  xiy(n)* etaz(n)-c3*  xiy(n)*zetaz(n)
       stiffl(13,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(14,16)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(15,16)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(15,17)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(14,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(13,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(15,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(14,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(13,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*  xix(n)*  xix(n)+c1*zetax(n)*  xix(n)
     .         +c6* etax(n)* etax(n)-c4*zetax(n)*zetax(n)
       term14= -c4*  xiy(n)*  xiy(n)+c1*zetay(n)*  xiy(n)
     .         +c6* etay(n)* etay(n)-c4*zetay(n)*zetay(n)
       term15= -c4*  xiz(n)*  xiz(n)+c1*zetaz(n)*  xiz(n)
     .         +c6* etaz(n)* etaz(n)-c4*zetaz(n)*zetaz(n)
       term16= -c4*  xiy(n)*  xix(n)+c3*zetay(n)*  xix(n)
     .         +c6* etay(n)* etax(n)+c3*  xiy(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*(- etax(n)*  xiy(n) +   xix(n)* etay(n)
     .             -zetax(n)* etay(n) +  etax(n)*zetay(n))
       term18= -c4*  xiz(n)*  xix(n)+c3*zetaz(n)*  xix(n)
     .         +c6* etaz(n)* etax(n)+c3*  xiz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*(- etax(n)*  xiz(n) +   xix(n)* etaz(n)
     .             -zetax(n)* etaz(n) +  etax(n)*zetaz(n))
       term20= -c4*  xiz(n)*  xiy(n)+c3*zetaz(n)*  xiy(n)
     .         +c6* etaz(n)* etay(n)+c3*  xiz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*(- etay(n)*  xiz(n) +   xiy(n)* etaz(n)
     .             -zetay(n)* etaz(n) +  etay(n)*zetaz(n))
       stiffl(13,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(14,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(15,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(15,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(14,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(13,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(15,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(14,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(13,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*  xix(n)*  xix(n)+c3* etax(n)*  xix(n)
     .         +c4* etax(n)* etax(n)-c2*zetax(n)*zetax(n)
       term14=  c4*  xiy(n)*  xiy(n)+c3* etay(n)*  xiy(n)
     .         +c4* etay(n)* etay(n)-c2*zetay(n)*zetay(n)
       term15=  c4*  xiz(n)*  xiz(n)+c3* etaz(n)*  xiz(n)
     .         +c4* etaz(n)* etaz(n)-c2*zetaz(n)*zetaz(n)
       term16=  c4*  xiy(n)*  xix(n)+c5* etay(n)*  xix(n)
     .         +c5*  xiy(n)* etax(n)+c4* etay(n)* etax(n)
     .         -c2*zetay(n)*zetax(n)
       term17= c3*(-zetax(n)* xiy(n) -  zetax(n)* etay(n)
     .            +  xix(n)*zetay(n) +   etax(n)*zetay(n))
       term18=  c4*  xiz(n)*  xix(n)+c5* etaz(n)*  xix(n)
     .         +c5*  xiz(n)* etax(n)+c4* etaz(n)* etax(n)
     .         -c2*zetaz(n)*zetax(n)
       term19= c3*(-zetax(n)* xiz(n) -  zetax(n)* etaz(n)
     .            +  xix(n)*zetaz(n) +   etax(n)*zetaz(n))
       term20=  c4*  xiz(n)*  xiy(n)+c5* etaz(n)*  xiy(n)
     .         +c5*  xiz(n)* etay(n)+c4* etaz(n)* etay(n)
     .         -c2*zetaz(n)*zetay(n)
       term21= c3*(-zetay(n)* xiz(n) -  zetay(n)* etaz(n)
     .            +  xiy(n)*zetaz(n) +   etay(n)*zetaz(n))
       stiffl(13,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(14,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(15,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(15,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(14,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(13,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(15,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(14,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(13,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
c
c     Rows 16-18
c
c
c
       term13=  c2*  xix(n)*  xix(n)-c1* etax(n) * xix(n)
     .         +c1*zetax(n)*  xix(n)+c2* etax(n)* etax(n)
     .         -c1*zetax(n)* etax(n)+c2*zetax(n)*zetax(n)
       term14=  c2*  xiy(n)*  xiy(n)-c1* etay(n) * xiy(n)
     .         +c1*zetay(n)*  xiy(n)+c2* etay(n)* etay(n)
     .         -c1*zetay(n)* etay(n)+c2*zetay(n)*zetay(n)
       term15=  c2*  xiz(n)*  xiz(n)-c1* etaz(n) * xiz(n)
     .         +c1*zetaz(n)*  xiz(n)+c2* etaz(n)* etaz(n)
     .         -c1*zetaz(n)* etaz(n)+c2*zetaz(n)*zetaz(n)
       term16=  c2*  xiy(n)*  xix(n)-c3* etay(n)*  xix(n)
     .         +c3*zetay(n)*  xix(n)-c3*  xiy(n)* etax(n)
     .         +c2* etay(n)* etax(n)-c3*zetay(n)* etax(n)
     .         +c3*  xiy(n)*zetax(n)-c3* etay(n)*zetax(n)
     .         +c2*zetay(n)*zetax(n)
       term17= 0.
       term18=  c2*  xiz(n)*  xix(n)-c3* etaz(n)*  xix(n)
     .         +c3*zetaz(n)*  xix(n)-c3*  xiz(n)* etax(n)
     .         +c2* etaz(n)* etax(n)-c3*zetaz(n)* etax(n)
     .         +c3*  xiz(n)*zetax(n)-c3* etaz(n)*zetax(n)
     .         +c2*zetaz(n)*zetax(n)
       term19= 0.
       term20=  c2*  xiy(n)*  xiz(n)-c3* etay(n)*  xiz(n)
     .         +c3*zetay(n)*  xiz(n)-c3*  xiy(n)* etaz(n)
     .         +c2* etay(n)* etaz(n)-c3*zetay(n)* etaz(n)
     .         +c3*  xiy(n)*zetaz(n)-c3* etay(n)*zetaz(n)
     .         +c2*zetay(n)*zetaz(n)
       term21= 0.
       stiffl(16,16)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(17,17)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(16,17)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(18,18)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(17,18)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(16,18)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13=  c4*  xix(n)*  xix(n)-c3* etax(n)*  xix(n)
     .         +c4* etax(n)* etax(n)-c2*zetax(n)*zetax(n)
       term14=  c4*  xiy(n)*  xiy(n)-c3* etay(n)*  xiy(n)
     .         +c4* etay(n)* etay(n)-c2*zetay(n)*zetay(n)
       term15=  c4*  xiz(n)*  xiz(n)-c3* etaz(n)*  xiz(n)
     .         +c4* etaz(n)* etaz(n)-c2*zetaz(n)*zetaz(n)
       term16=  c4*  xiy(n)*  xix(n)-c5* etay(n)*  xix(n)
     .         -c5*  xiy(n)* etax(n)+c4* etay(n)* etax(n)
     .         -c2*zetay(n)*zetax(n)
       term17= c3*( zetax(n)*  xiy(n) - zetax(n)* etay(n)
     .             -  xix(n)*zetay(n) +  etax(n)*zetay(n))
       term18=  c4*  xiz(n)*  xix(n)-c5* etaz(n)*  xix(n)
     .         -c5*  xiz(n)* etax(n)+c4* etaz(n)* etax(n)
     .         -c2*zetaz(n)*zetax(n)
       term19= c3*( zetax(n)*  xiz(n) - zetax(n)* etaz(n)
     .             -  xix(n)*zetaz(n) +  etax(n)*zetaz(n))
       term20=  c4*  xiz(n)*  xiy(n)-c5* etaz(n)*  xiy(n)
     .         -c5*  xiz(n)* etay(n)+c4* etaz(n)* etay(n)
     .         -c2*zetaz(n)*zetay(n)
       term21= c3*( zetay(n)*  xiz(n) - zetay(n)* etaz(n)
     .             -  xiy(n)*zetaz(n) +  etay(n)*zetaz(n))
       stiffl(16,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(17,19)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(18,19)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(18,20)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(17,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(16,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(18,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(17,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(16,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
       term13= -c4*  xix(n)*  xix(n)-c1*zetax(n)*  xix(n)
     .         +c6* etax(n)* etax(n)-c4*zetax(n)*zetax(n)
       term14= -c4*  xiy(n)*  xiy(n)-c1*zetay(n)*  xiy(n)
     .         +c6* etay(n)* etay(n)-c4*zetay(n)*zetay(n)
       term15= -c4*  xiz(n)*  xiz(n)-c1*zetaz(n)*  xiz(n)
     .         +c6* etaz(n)* etaz(n)-c4*zetaz(n)*zetaz(n)
       term16= -c4*  xiy(n)*  xix(n)-c3*zetay(n)*  xix(n)
     .         +c6* etay(n)* etax(n)-c3*  xiy(n)*zetax(n)
     .         -c4*zetay(n)*zetax(n)
       term17= c5*(  etax(n)* xiy(n) -    xix(n)* etay(n)
     .            -zetax(n)* etay(n) +   etax(n)*zetay(n))
       term18= -c4*  xiz(n)*  xix(n)-c3*zetaz(n)*  xix(n)
     .         +c6* etaz(n)* etax(n)-c3*  xiz(n)*zetax(n)
     .         -c4*zetaz(n)*zetax(n)
       term19= c5*(  etax(n)* xiz(n) -    xix(n)* etaz(n)
     .            -zetax(n)* etaz(n) +   etax(n)*zetaz(n))
       term20= -c4*  xiz(n)*  xiy(n)-c3*zetaz(n)*  xiy(n)
     .         +c6* etaz(n)* etay(n)-c3*  xiz(n)*zetay(n)
     .         -c4*zetaz(n)*zetay(n)
       term21= c5*(  etay(n)* xiz(n) -    xiy(n)* etaz(n)
     .            -zetay(n)* etaz(n) +   etay(n)*zetaz(n))
       stiffl(16,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(17,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(18,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(18,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(17,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(16,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(18,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(17,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(16,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
c
c
c
c     Rows 19-21
c
c
c
       term13=  c2*  xix(n)*  xix(n)-c1* etax(n)*  xix(n)
     .         -c1*zetax(n)*  xix(n)+c2* etax(n)* etax(n)
     .         +c1*zetax(n)* etax(n)+c2*zetax(n)*zetax(n)
       term14=  c2*  xiy(n)*  xiy(n)-c1* etay(n)*  xiy(n)
     .         -c1*zetay(n)*  xiy(n)+c2* etay(n)* etay(n)
     .         +c1*zetay(n)* etay(n)+c2*zetay(n)*zetay(n)
       term15=  c2*  xiz(n)*  xiz(n)-c1* etaz(n)*  xiz(n)
     .         -c1*zetaz(n)*  xiz(n)+c2* etaz(n)* etaz(n)
     .         +c1*zetaz(n)* etaz(n)+c2*zetaz(n)*zetaz(n)
       term16=  c2*  xiy(n)*  xix(n)-c3* etay(n)*  xix(n)
     .         -c3*zetay(n)*  xix(n)-c3*  xiy(n)* etax(n)
     .         +c2* etay(n)* etax(n)+c3*zetay(n)* etax(n)
     .         -c3*  xiy(n)*zetax(n)+c3* etay(n)*zetax(n)
     .         +c2*zetay(n)*zetax(n)
       term17= 0.
       term18=  c2*  xiz(n)*  xix(n)-c3* etaz(n)*  xix(n)
     .         -c3*zetaz(n)*  xix(n)-c3*  xiz(n)* etax(n)
     .         +c2* etaz(n)* etax(n)+c3*zetaz(n)* etax(n)
     .         -c3*  xiz(n)*zetax(n)+c3* etaz(n)*zetax(n)
     .         +c2*zetaz(n)*zetax(n)
       term19= 0.
       term20=  c2*  xiz(n)*  xiy(n)-c3* etaz(n)*  xiy(n)
     .         -c3*zetaz(n)*  xiy(n)-c3*  xiz(n)* etay(n)
     .         +c2* etaz(n)* etay(n)+c3*zetaz(n)* etay(n)
     .         -c3*  xiz(n)*zetay(n)+c3* etaz(n)*zetay(n)
     .         +c2*zetaz(n)*zetay(n)
       term21= 0.
       stiffl(19,19)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(20,20)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(19,20)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(21,21)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(20,21)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(19,21)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
       term13= -c2*  xix(n)*  xix(n)+c4* etax(n)* etax(n)
     .         +c3*zetax(n)* etax(n)+c4*zetax(n)*zetax(n)
       term14= -c2*  xiy(n)*  xiy(n)+c4* etay(n)* etay(n)
     .         +c3*zetay(n)* etay(n)+c4*zetay(n)*zetay(n)
       term15= -c2*  xiz(n)*  xiz(n)+c4* etaz(n)* etaz(n)
     .         +c3*zetaz(n)* etaz(n)+c4*zetaz(n)*zetaz(n)
       term16= -c2*  xiy(n)*  xix(n)+c4* etay(n)* etax(n)
     .         +c5*zetay(n)* etax(n)+c5* etay(n)*zetax(n)
     .         +c4*zetay(n)*zetax(n)
       term17= c3*(  etax(n)* xiy(n) +  zetax(n)*  xiy(n)
     .            -  xix(n)* etay(n) -    xix(n)*zetay(n))
       term18= -c2*  xiz(n)*  xix(n)+c4* etaz(n)* etax(n)
     .         +c5*zetaz(n)* etax(n)+c5* etaz(n)*zetax(n)
     .         +c4*zetaz(n)*zetax(n)
       term19= c3*(  etax(n)* xiz(n) +  zetax(n)*  xiz(n)
     .            -  xix(n)* etaz(n) -    xix(n)*zetaz(n))
       term20= -c2*  xiz(n)*  xiy(n)+c4* etaz(n)* etay(n)
     .         +c5*zetaz(n)* etay(n)+c5* etaz(n)*zetay(n)
     .         +c4*zetaz(n)*zetay(n)
       term21= c3*(  etay(n)* xiz(n) +  zetay(n)*  xiz(n)
     .            -  xiy(n)* etaz(n) -    xiy(n)*zetaz(n))
       stiffl(19,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(20,22)=c(4,1)*term13+c(2,4)*term14+c(5,6)*term15
     .             +(c(2,1)+c(4,4))*term16+(c(4,4)-c(2,1))*term17
     .             +(c(5,1)+c(4,6))*term18+(c(4,6)-c(5,1))*term19
     .             +(c(5,4)+c(2,6))*term20+(c(2,6)-c(5,4))*term21
       stiffl(21,22)=c(6,1)*term13+c(5,4)*term14+c(3,6)*term15
     .             +(c(5,1)+c(6,4))*term16+(c(6,4)-c(5,1))*term17
     .             +(c(3,1)+c(6,6))*term18+(c(6,6)-c(3,1))*term19
     .             +(c(3,4)+c(5,6))*term20+(c(5,6)-c(3,4))*term21
       stiffl(21,23)=c(6,4)*term13+c(5,2)*term14+c(3,5)*term15
     .             +(c(5,4)+c(6,2))*term16+(c(6,2)-c(5,4))*term17
     .             +(c(3,4)+c(6,5))*term18+(c(6,5)-c(3,4))*term19
     .             +(c(3,2)+c(5,5))*term20+(c(5,5)-c(3,2))*term21
       stiffl(20,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(19,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(21,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(20,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(19,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
c
c
c
c     Rows 22-24
c
c
c
c
       term13=  c2*  xix(n)*  xix(n)+c1* etax(n)*  xix(n)
     .         +c1*zetax(n)*  xix(n)+c2* etax(n)* etax(n)
     .         +c1*zetax(n)* etax(n)+c2*zetax(n)*zetax(n)
       term14=  c2*  xiy(n)*  xiy(n)+c1* etay(n)*  xiy(n)
     .         +c1*zetay(n)*  xiy(n)+c2* etay(n)* etay(n)
     .         +c1*zetay(n)* etay(n)+c2*zetay(n)*zetay(n)
       term15=  c2*  xiz(n)*  xiz(n)+c1* etaz(n)*  xiz(n)
     .         +c1*zetaz(n)*  xiz(n)+c2* etaz(n)* etaz(n)
     .         +c1*zetaz(n)* etaz(n)+c2*zetaz(n)*zetaz(n)
       term16=  c2*  xiy(n)*  xix(n)+c3* etay(n)*  xix(n)
     .         +c3*zetay(n)*  xix(n)+c3*  xiy(n)* etax(n)
     .         +c2* etay(n)* etax(n)+c3*zetay(n)* etax(n)
     .         +c3*  xiy(n)*zetax(n)+c3* etay(n)*zetax(n)
     .         +c2*zetay(n)*zetax(n)
       term17= 0.
       term18=  c2*  xiz(n)*  xix(n)+c3* etaz(n)*  xix(n)
     .         +c3*zetaz(n)*  xix(n)+c3*  xiz(n)* etax(n)
     .         +c2* etaz(n)* etax(n)+c3*zetaz(n)* etax(n)
     .         +c3*  xiz(n)*zetax(n)+c3* etaz(n)*zetax(n)
     .         +c2*zetaz(n)*zetax(n)
       term19= 0.
       term20=  c2*  xiz(n)*  xiy(n)+c3* etaz(n)*  xiy(n)
     .         +c3*zetaz(n)*  xiy(n)+c3*  xiz(n)* etay(n)
     .         +c2* etaz(n)* etay(n)+c3*zetaz(n)* etay(n)
     .         +c3*  xiz(n)*zetay(n)+c3* etaz(n)*zetay(n)
     .         +c2*zetaz(n)*zetay(n)
       term21= 0.
       stiffl(22,22)=c(1,1)*term13+c(4,4)*term14+c(6,6)*term15
     .             +(c(1,4)+c(4,1))*term16+(c(1,4)-c(4,1))*term17
     .             +(c(1,6)+c(6,1))*term18+(c(1,6)-c(6,1))*term19
     .             +(c(4,6)+c(6,4))*term20+(c(4,6)-c(6,4))*term21
       stiffl(23,23)=c(4,4)*term13+c(2,2)*term14+c(5,5)*term15
     .             +(c(2,4)+c(4,2))*term16+(c(4,2)-c(2,4))*term17
     .             +(c(5,4)+c(4,5))*term18+(c(4,5)-c(5,4))*term19
     .             +(c(5,2)+c(2,5))*term20+(c(2,5)-c(5,2))*term21
       stiffl(22,23)=c(1,4)*term13+c(4,2)*term14+c(6,5)*term15
     .             +(c(4,4)+c(1,2))*term16+(c(1,2)-c(4,4))*term17
     .             +(c(6,4)+c(1,5))*term18+(c(1,5)-c(6,4))*term19
     .             +(c(6,2)+c(4,5))*term20+(c(4,5)-c(6,2))*term21
       stiffl(24,24)=c(6,6)*term13+c(5,5)*term14+c(3,3)*term15
     .             +(c(5,6)+c(6,5))*term16+(c(6,5)-c(5,6))*term17
     .             +(c(3,6)+c(6,3))*term18+(c(6,3)-c(3,6))*term19
     .             +(c(3,5)+c(5,3))*term20+(c(5,3)-c(3,5))*term21
       stiffl(23,24)=c(4,6)*term13+c(2,5)*term14+c(5,3)*term15
     .             +(c(2,6)+c(4,5))*term16+(c(4,5)-c(2,6))*term17
     .             +(c(5,6)+c(4,3))*term18+(c(4,3)-c(5,6))*term19
     .             +(c(5,5)+c(2,3))*term20+(c(2,3)-c(5,5))*term21
       stiffl(22,24)=c(1,6)*term13+c(4,5)*term14+c(6,3)*term15
     .             +(c(4,6)+c(1,5))*term16+(c(1,5)-c(4,6))*term17
     .             +(c(6,6)+c(1,3))*term18+(c(1,3)-c(6,6))*term19
     .             +(c(6,5)+c(4,3))*term20+(c(4,3)-c(6,5))*term21
      do ii = 1,24
        do jj = ii+1,24
          stiffl(jj,ii) = stiffl(ii,jj)
        enddo
      enddo
c
c     n1(0) = n
c     n1(1) = j + 1
c     n1(2) = k + 1
c     n1(3) = j + 1,k + 1
c     n1(4) = i + 1
c     n1(5) = i + 1,k + 1
c     n1(6) = j + 1,i + 1
c     n1(7) = j + 1,i + 1,k + 1
c
           n11i    = 0
           n33i    = 0
           n44i    = 1
           do ii = 0,7
             n11i(1,ii) = n1(ii)
             iimax(ii)  = 1
             if(islavept(n1(ii),8,iseqr).eq.0) n44i(ii) = 0
             if(islavept(n1(ii),11,iseqr).gt.1) then
               iimax(ii)  = islavept(n1(ii),11,iseqr)
               do ii4 = 2,iimax(ii)
                 n11i(ii4,ii) = islavept(n1(ii),12+ii4-2,iseqr)
                 if(islavept(n11i(ii4,ii),8,iseqr).eq.0) then
                   n44i(ii) = 0
                 end if
               enddo
             end if
             do ii4 = 1,iimax(ii)
               n33i(ii4,ii) = 3*(n11i(ii4,ii)-1)
             enddo
           enddo
        oj  = 1./ooj(n)
c
        do ii = 0,7
         if(n44i(ii).ne.0) then
              do ii3 = 1,iimax(ii)
                  n22i = n11i(ii3,ii)
                     do j = 1,3
                      irow = n2(ii)+j-1
                      itst1= n33i(ii3,ii)+j
                      itst2= ija(itst1)
                      sa (itst1)=sa(itst1)+stiffl(irow,irow)*oj
                      i1 = 0
                      do j1= 1,3
                       if(j1.ne.j) then
                        sa (itst2+i1) = sa(itst2+i1)+
     .                               stiffl(irow,n2(ii)+j1-1)*oj
                        ija(itst2+i1) = n33i(ii3,ii)+j1
                        i1 = i1 + 1
                       end if
                      enddo
                      j2sta =  ija(itst1)+2
                      j2end =  ija(itst1+1)-1
                      do kk = 0,7
                       if(kk.ne.ii) then
                          do j2 = j2sta,j2end,3
                            if(ija(j2).eq.0) then
                             do j1 = 1,3
                               icol  = n2(kk)+j1-1
                               j3 = j2 + j1 - 1
                               sa (j3) = sa(j3)+stiffl(irow,icol)*oj
                               ija(j3) = n33i(1,kk)+j1
                             enddo
                             goto 2200
                            else
                              do ii4 = 1,iimax(kk)
                               itst1 = n33i(ii4,kk)+1
                               if(ija(j2).eq.itst1) then
                                do j1 = 1,3
                                  icol  = n2(kk)+j1-1
                                  j3 = j2 + j1 - 1
                                  sa (j3)=sa(j3)+stiffl(irow,icol)*oj
                                enddo
                                goto 2200
                               end if
                              enddo
                            end if
                          enddo
2200                      continue
                       end if
                      enddo
                     enddo
              enddo
3000       continue
         else
            do j = 1,3
              sa (3*(n1(ii)-1)+j) = 1.d0
            enddo
         end if
        enddo
        end if
200    continue
 
       return
       end
       subroutine s_gwrite(nbl,idim,jdim,kdim,ista,iend,jsta,jend,ksta,
     .                     kend,idim1,jdim1,kdim1,x,y,z)
       dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      open(5000+nbl,status='unknown',form='formatted')
      write(5000+nbl,'(i8)') 1
      write(5000+nbl,'(3i8)')idim1,jdim1,kdim1
      write(5000+nbl,'(5e21.12)') (((x(j,k,i),i=ista,iend),j=jsta,jend),
     .                               k=ksta,kend)
      write(5000+nbl,'(5e21.12)') (((y(j,k,i),i=ista,iend),j=jsta,jend),
     .                               k=ksta,kend)
      write(5000+nbl,'(5e21.12)') (((z(j,k,i),i=ista,iend),j=jsta,jend),
     .                               k=ksta,kend)
       close(5000+nbl)
 
       return
       end
      subroutine s_gread(nbl,idim,jdim,kdim,ista,iend,jsta,jend,ksta,
     .                   kend,idim1,jdim1,kdim1,x,y,z)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      open(5000+nbl,status='old',form='formatted')
      read(5000+nbl,*) nblocks
      read(5000+nbl,*)idim1,jdim1,kdim1
      read(5000+nbl,*)
     .          (((x(j,k,i),i=ista,iend),j=jsta,jend),k=ksta,kend)
     .         ,(((y(j,k,i),i=ista,iend),j=jsta,jend),k=ksta,kend)
     .         ,(((z(j,k,i),i=ista,iend),j=jsta,jend),k=ksta,kend)
      close(5000+nbl)
      return
      end
      subroutine gridgen_glyf1(nbl)
      character*8 bgcontrl,fgcontl
      character*5 filnm
      character*3 char1
      character*24 filnm2
c
      open(9700,file='gridgen.options1',status='old',
     .    form='formatted')
      read(9700,12030) bgcontrl,fgcontrl
12030 format(2a8)
      read(9700,*) xmgrlx,xmgprlng,iter
      close(9700)
      filnm2 = 'dgplot3d_gridgen0001.glf'
      if(nbl.lt.10) then
        write(char1,'(i1)') nbl
        filnm2(19:19) = char1(1:1)
      else if(nbl.lt.100) then
        write(char1,'(i2)') nbl
        filnm2(18:19) = char1(1:2)
      else
        write(char1,'(i3)') nbl
        filnm2(17:19) = char1(1:3)
      end if
      open(9900,file=filnm2,status='unknown',form='formatted')
      filnm(1:5) = 'fort.'
      iflnum = 5000+nbl
      write(9900,11021) filnm,iflnum
      write(9900,11022) 1,1,0
      write(9900,11025)
      write(9900,11026)  1
      write(9900,11039) bgcontrl,fgcontrl,xmgrlx,xmgprlng,iter
     .                 ,filnm,iflnum
11021 format('# Gridgen Journal File V1 (Gridgen 15.09 REL 3)',
     .   /,'# Created Sun Aug 13 14:47:02 2006',/,
     .     '#',/,'package require PWI_Glyph 1.6.9',/,/,
     .     '# Delete any existing grids and database entities.  ',
     .     'Reset AS/W, defaults, and',/,'# tolerances.',/,
     .     'gg::memClear',/,'gg::aswDeleteBC -glob "*"',/,
     .     'gg::aswDeleteVC -glob "*"',/,
     .     'gg::aswSet GENERIC -dim 3',/,'gg::defReset',/,
     .     'gg::tolReset',/,'# Check for user input after ',
     .     'each function, but delay screen updates until',/,
     .     '# script is finished.',/,'gg::updatePolicy INPUT_ONLY',/,/,
     .     'set _ggTemp_(1) [gg::blkImport "',a5,i4,'" \ ',/,
     .     '   -style PLOT3D -form ASCII -precision DOUBLE]')
11022 format('set _BL(',i1,') [lindex $_ggTemp_(',i1,') ',i1,']')
11025 format('unset _ggTemp_(1)')
11026 format('set _ggTemp_(2) [list',' $_BL(',i1,')]')
11039 format('gg::blkEllSolverBegin $_ggTemp_(2)',/,
     .'# gg::blkTFISolverRun $_ggTemp_(2)',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -bg_control ',a8,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -fg_control ',a8,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_calc INTERPOLATE',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -angle_calc INTERPOLATE',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_blend 2',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -angle_blend 2',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -mg_relax ',f3.2,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -mg_prolongation ',f3.2,/,
     .'  gg::blkEllSolverStep -iterations ',i5,' -nodisplay',/,
     .'gg::blkEllSolverEnd',/,'unset _ggTemp_(2)',/,
     .'gg::blkExport ALL "',a5,i4,'" \ ',/,
     .'   -style PLOT3D -form ASCII -precision DOUBLE',/,
     .'gg::terminate')
       close(9900)
      return
      end
      subroutine gridgen_glyf2(nbl)
      character*8 bgcontrl,fgcontl
      character*5 filnm
      character*3 char1
      character*24 filnm2
c
      open(9700,file='gridgen.options2',status='old',
     .    form='formatted')
      read(9700,12030) bgcontrl,fgcontrl
12030 format(2a8)
      read(9700,*) xmgrlx,xmgprlng,iter
      close(9700)
      filnm2 = 'dgplot3d_gridgen0002.glf'
      if(nbl.lt.10) then
        write(char1,'(i1)') nbl
        filnm2(19:19) = char1(1:1)
      else if(nbl.lt.100) then
        write(char1,'(i2)') nbl
        filnm2(18:19) = char1(1:2)
      else
        write(char1,'(i3)') nbl
        filnm2(17:19) = char1(1:3)
      end if
      open(9900,file=filnm2,status='unknown',form='formatted')
      filnm(1:5) = 'fort.'
      iflnum = 5000+nbl
      write(9900,11021) filnm,iflnum
      write(9900,11022) 1,1,0
      write(9900,11025)
      write(9900,11026)  1
      write(9900,11039) bgcontrl,fgcontrl,xmgrlx,xmgprlng,iter
     .                 ,filnm,iflnum
11021 format('# Gridgen Journal File V1 (Gridgen 15.09 REL 3)',
     .   /,'# Created Sun Aug 13 14:47:02 2006',/,
     .     '#',/,'package require PWI_Glyph 1.6.9',/,/,
     .     '# Delete any existing grids and database entities.  ',
     .     'Reset AS/W, defaults, and',/,'# tolerances.',/,
     .     'gg::memClear',/,'gg::aswDeleteBC -glob "*"',/,
     .     'gg::aswDeleteVC -glob "*"',/,
     .     'gg::aswSet GENERIC -dim 3',/,'gg::defReset',/,
     .     'gg::tolReset',/,'# Check for user input after ',
     .     'each function, but delay screen updates until',/,
     .     '# script is finished.',/,'gg::updatePolicy INPUT_ONLY',/,/,
     .     'set _ggTemp_(1) [gg::blkImport "',a5,i4,'" \ ',/,
     .     '   -style PLOT3D -form ASCII -precision DOUBLE]')
11022 format('set _BL(',i1,') [lindex $_ggTemp_(',i1,') ',i1,']')
11025 format('unset _ggTemp_(1)')
11026 format('set _ggTemp_(2) [list',' $_BL(',i1,')]')
11039 format('gg::blkEllSolverBegin $_ggTemp_(2)',/,
     .'# gg::blkTFISolverRun $_ggTemp_(2)',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -bg_control ',a8,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -fg_control ',a8,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_calc INTERPOLATE',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -angle_calc INTERPOLATE',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -spacing_blend 2',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -angle_blend 2',/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -mg_relax ',f3.2,/,
     .'  gg::blkEllSolverAtt $_ggTemp_(2) -mg_prolongation ',f3.2,/,
     .'  gg::blkEllSolverStep -iterations ',i5,' -nodisplay',/,
     .'gg::blkEllSolverEnd',/,'unset _ggTemp_(2)',/,
     .'gg::blkExport ALL "',a5,i4,'" \ ',/,
     .'   -style PLOT3D -form ASCII -precision DOUBLE',/,
     .'gg::terminate')
       close(9900)
      return
      end
