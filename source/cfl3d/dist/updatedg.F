      subroutine updatedg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                    maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                    nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                    nblock,levelg,igridg,idefrm,ncgg,iadvance,nou,
     .                    bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                    utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                    xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                    rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,itrans,irotat,slavept,nslave,
     .                    iskip,jskip,kskip,xs,xxn,nsegdfrm,idfrmseg,
     .                    iaesurf,maxsegdg,iwk,nmaster,nt,xorig,
     .                    yorig,zorig,xorgae0,yorgae0,zorgae0,icouple,
     .                    ireq,nnodes)
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
      dimension aesrfdat(5,maxaes),slavept(nslave,3,nmaster,5)
      dimension iskip(maxbl),jskip(maxbl),kskip(maxbl)
      dimension xorgae0(maxbl,maxsegdg),yorgae0(maxbl,maxsegdg),
     .          zorgae0(maxbl,maxsegdg),icouple(maxbl,maxsegdg)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension ireq(maxbl)
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
      common /unst/ time,cfltau,ntstep,ita,iunst
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
         call writ_buf(1,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
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
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
  444       format(38h not enough work space for subroutines,
     .             23h trnsurf/rotsurf/aesurf)
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
     .                  maxbl,maxseg,nmds,maxaes,aesrfdat,
     .                  xs,xxn,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                  nsegdfrm,idfrmseg,iaesurf,maxsegdg)
c
         end if
      end if
c
      if (isklton.eq.1) then
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
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
               call MPI_Irecv (wk(ist),3,MY_MPI_REAL,
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
c
c     update faces and interior of deforming meshes, and finite
c     difference old and new grid positions to get new grid velocities
c     (velocities stored in work array, at position lvel)
c
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
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  445          format(37h not enough work space for subroutine,
     .                 7h deform)
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
     .                     w(lxnm2),w(lynm2),w(lznm2),
     .                     w(lxnm1),w(lynm1),w(lznm1),
     .                     w(ldeltj),w(ldeltk),w(ldelti),wk(lvel),
     .                     icsi,icsf,jcsi,jcsf,kcsi,kcsf,maxbl,time,
     .                     dt,ita,nou,bou,nbuf,ibufdim,myid,idefrm,
     .                     nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                     ibcinfo,jbcinfo,kbcinfo,maxseg,wk,nsurf,irst,
     .                     nflag,ivert,slavept,nslave,iskip,jskip,
     .                     kskip,nsegdfrm,idfrmseg,iaesurf,maxsegdg,
     .                     nmaster,iseqr)
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
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  425          format(45h not enough work space for metric subroutines)
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
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  430          format(45h not enough work space for subroutine cellvol)
c
               call cellvol(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                      w(lsk),w(lsi),w(lvol),wk(lt4wk),nou,bou,
     .                      nbuf,ibufdim,myid,mblk2nd,maxbl,nbl)
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
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
c
  200    continue
      end if
c
      return
      end
