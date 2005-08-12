      subroutine setslave(lw,lw2,w,mgwk,wk,nwork,maxbl,maxgr,maxseg,
     .                    nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .                    ibcinfo,jbcinfo,kbcinfo,nblock,idefrm,
     .                    iadvance,nou,bou,nbuf,ibufdim,myid,myhost,
     .                    mycomm,mblk2nd,icsi,icsf,jcsi,jcsf,kcsi,
     .                    kcsf,slavept,nslave,nsegdfrm,idfrmseg,
     .                    maxsegdg,iwk,iwork,nmaster,ngrid,jskip,kskip,
     .                    iskip,nblg,levelg,lfgm)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set up slave point data array for deforming meshes
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
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),wk(nwork),iwk(iwork)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension idefrm(maxbl),nblg(maxgr)
      dimension iadvance(maxbl),mblk2nd(maxbl)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg)
      dimension iskip(maxbl),jskip(maxbl),kskip(maxbl)
      dimension slavept(nslave,3,nmaster,5),levelg(maxbl)
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
      common /elastic_ss/ idef_ss
      common /twod/ i2d
      common /deformz/ beta1,alpha1,isktyp,negvol,meshdef,nsprgit
c
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset  = maxbl
      itag_ns  = 1
      itag_wk  = itag_ns + ioffset
      itag_org = itag_wk + ioffset
#endif
c
      do iseq=1,mseq
c
      lglobal = lfgm-(mseq-iseq)
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
     .      levelg(nbl).le.levelt(iseq))) then
            if (idefrm(nbl).gt.0 .and. idefrm(nbl).lt.999) then
               call lead(nbl,lw,lw2,maxbl)
               if (myid.eq.mblk2nd(nbl)) then
                  call getsurf(w(lx),w(ly),w(lz),w(ldeltj),
     .                         w(ldeltk),w(ldelti),w(lxnm1),
     .                         w(lynm1),w(lznm1),icsi,icsf,
     .                         jcsi,jcsf,kcsi,kcsf,wk,nwork,
     .                         nbl,idim,jdim,kdim,nsurf,iwk(nbl),
     .                         nsegdfrm,maxbl,idfrmseg,maxsegdg)
               end if
#if defined DIST_MPI
               mytag = itag_ns + nbl
               nd_srce = mblk2nd(nbl)
               if (myid.eq.mblk2nd(nbl)) then
                  call MPI_Send(iwk(nbl),1,MPI_INTEGER,myhost,mytag,
     .                          mycomm,ierr)
               else if (myid.eq.myhost) then
                  call MPI_Recv(iwk(nbl),1,MPI_INTEGER,nd_srce,mytag,
     .                          mycomm,istat,ierr)
               end if
               if (iwk(nbl).gt.0) then
                  mytag   = itag_wk + nbl
                  nd_srce = mblk2nd(nbl)
                  numdat  = 9*iwk(nbl)
                  if (myid.eq.mblk2nd(nbl)) then
                     call MPI_Send(wk,numdat,MY_MPI_REAL,
     .                             myhost,mytag,mycomm,ierr) 
                  else if (myid.eq.myhost) then
                     ns = 9*nsurf+1
                     call MPI_Recv(wk(ns),numdat,MY_MPI_REAL,
     .                             nd_srce,mytag,mycomm,istat,ierr)
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
      call MPI_Bcast (wk,nsurf*9,MY_MPI_REAL,
     .                myhost,mycomm,ierr)
#endif
c
      if (nsurf.gt.iwork) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping in setslave...parameter '',
     .                         ''iwork too small'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      if (iunst.gt.1 .or. idef_ss.gt.0) then
         ivert = 0
         do nbl=1,nblock
         if ((levelg(nbl).ge.lglobal .and.
     .      levelg(nbl).le.levelt(iseq))) then
            if (myid.eq.mblk2nd(nbl) .and. idefrm(nbl) .gt. 0) then
               call lead(nbl,lw,lw2,maxbl)
c
c              determine the master surface point for the 4 corner
c              points of each block subface that is NOT on a solid
c              surface.
c
               iskp = iskip(nbl)
               jskp = jskip(nbl)
               kskp = kskip(nbl)
c
c              i=1 subfaces
c
               do nseg=1,nbci0(nbl)
                  ii     = 1
                  ibctyp = ibcinfo(nbl,nseg,1,ii)
                  if (abs(ibctyp).ne.2004 .and.
     .                abs(ibctyp).ne.2014 .and.
     .                abs(ibctyp).ne.2016 .and.
     .                abs(ibctyp).ne.1005 .and.
     .                abs(ibctyp).ne.1006 .or.
     .                real(beta1).lt.0.) then
                      js = ibcinfo(nbl,nseg,2,ii)
                      je = ibcinfo(nbl,nseg,3,ii)
                      ks = ibcinfo(nbl,nseg,4,ii)
                      ke = ibcinfo(nbl,nseg,5,ii)
                      do j=js,je,jskp
                         do k=ks,ke,kskp
                            call setcorner(j,k,1,w(lx),w(ly),
     .                                     w(lz),w(lx),w(ly),w(lz),
     .                                     jdim,kdim,idim,wk,nsurf,1,
     .                                     ivert,slavept,nslave,nou,
     .                                     bou,ibufdim,nbuf,myid,
     .                                     iwk,nmaster,iseq)
                         end do
                      end do
                  end if
               end do
c
c              i=idim subfaces
c
               do nseg=1,nbcidim(nbl)
                  ii     = 2
                  ibctyp = ibcinfo(nbl,nseg,1,ii)
                  if (abs(ibctyp).ne.2004 .and.
     .                abs(ibctyp).ne.2014 .and.
     .                abs(ibctyp).ne.2016 .and.
     .                abs(ibctyp).ne.1005 .and.
     .                abs(ibctyp).ne.1006 .or.
     .                real(beta1).lt.0.) then
                      js = ibcinfo(nbl,nseg,2,ii)
                      je = ibcinfo(nbl,nseg,3,ii)
                      ks = ibcinfo(nbl,nseg,4,ii)
                      ke = ibcinfo(nbl,nseg,5,ii)
                      do j=js,je,jskp
                         do k=ks,ke,kskp
                            call setcorner(j,k,idim,w(lx),w(ly),
     .                                     w(lz),w(lx),w(ly),w(lz),
     .                                     jdim,kdim,idim,wk,nsurf,1,
     .                                     ivert,slavept,nslave,nou,
     .                                     bou,ibufdim,nbuf,myid,
     .                                     iwk,nmaster,iseq)
                         end do
                      end do
                  end if
               end do
c
               if (i2d .eq. 0) then
c
c                 j=1 subfaces
c
                  do nseg=1,nbcj0(nbl)
                     jj     = 1
                     jbctyp = jbcinfo(nbl,nseg,1,jj)
                     if (abs(jbctyp).ne.2004 .and.
     .                   abs(jbctyp).ne.2014 .and.
     .                   abs(jbctyp).ne.2016 .and.
     .                   abs(jbctyp).ne.1005 .and.
     .                   abs(jbctyp).ne.1006 .or.
     .                   real(beta1).lt.0.) then
                         is = jbcinfo(nbl,nseg,2,jj)
                         ie = jbcinfo(nbl,nseg,3,jj)
                         ks = jbcinfo(nbl,nseg,4,jj)
                         ke = jbcinfo(nbl,nseg,5,jj)
                         do i=is,ie,iskp
                            do k=ks,ke,kskp
                               call setcorner(1,k,i,w(lx),w(ly),
     .                                        w(lz),w(lx),w(ly),w(lz),
     .                                        jdim,kdim,idim,wk,nsurf,1,
     .                                        ivert,slavept,nslave,nou,
     .                                        bou,ibufdim,nbuf,myid,
     .                                        iwk,nmaster,iseq)
                            end do
                         end do
                     end if
                  end do
c
c                 j=jdim subfaces
c
                  do nseg=1,nbcjdim(nbl)
                     jj     = 2
                     jbctyp = jbcinfo(nbl,nseg,1,jj)
                     if (abs(jbctyp).ne.2004 .and.
     .                   abs(jbctyp).ne.2014 .and.
     .                   abs(jbctyp).ne.2016 .and.
     .                   abs(jbctyp).ne.1005 .and.
     .                   abs(jbctyp).ne.1006 .or.
     .                   real(beta1).lt.0.) then
                         is = jbcinfo(nbl,nseg,2,jj)
                         ie = jbcinfo(nbl,nseg,3,jj)
                         ks = jbcinfo(nbl,nseg,4,jj)
                         ke = jbcinfo(nbl,nseg,5,jj)
                         do i=is,ie,iskp
                            do k=ks,ke,kskp
                               call setcorner(jdim,k,i,w(lx),w(ly),
     .                                        w(lz),w(lx),w(ly),w(lz),
     .                                        jdim,kdim,idim,wk,nsurf,1,
     .                                        ivert,slavept,nslave,nou,
     .                                        bou,ibufdim,nbuf,myid,
     .                                        iwk,nmaster,iseq)
                            end do
                         end do
                     end if
                  end do
c
c                 k=1 subfaces
c
                  do nseg=1,nbck0(nbl)
                     kk     = 1
                     kbctyp = kbcinfo(nbl,nseg,1,kk)
                     if (abs(kbctyp).ne.2004 .and.
     .                   abs(kbctyp).ne.2014 .and.
     .                   abs(kbctyp).ne.2016 .and.
     .                   abs(kbctyp).ne.1005 .and.
     .                   abs(kbctyp).ne.1006 .or.
     .                   real(beta1).lt.0.) then
                         is = kbcinfo(nbl,nseg,2,kk)
                         ie = kbcinfo(nbl,nseg,3,kk)
                         js = kbcinfo(nbl,nseg,4,kk)
                         je = kbcinfo(nbl,nseg,5,kk)
                         do i=is,ie,iskp
                            do j=js,je,jskp
                               call setcorner(j,1,i,w(lx),w(ly),
     .                                        w(lz),w(lx),w(ly),w(lz),
     .                                        jdim,kdim,idim,wk,nsurf,1,
     .                                        ivert,slavept,nslave,nou,
     .                                        bou,ibufdim,nbuf,myid,
     .                                        iwk,nmaster,iseq)
                            end do
                         end do
                     end if
                  end do
c
c                 k=kdim subfaces
c
                  do nseg=1,nbckdim(nbl)
                     kk     = 2
                     kbctyp = kbcinfo(nbl,nseg,1,kk)
                     if (abs(kbctyp).ne.2004 .and.
     .                   abs(kbctyp).ne.2014 .and.
     .                   abs(kbctyp).ne.2016 .and.
     .                   abs(kbctyp).ne.1005 .and.
     .                   abs(kbctyp).ne.1006 .or.
     .                   real(beta1).lt.0.) then
                         is = kbcinfo(nbl,nseg,2,kk)
                         ie = kbcinfo(nbl,nseg,3,kk)
                         js = kbcinfo(nbl,nseg,4,kk)
                         je = kbcinfo(nbl,nseg,5,kk)
                         do i=is,ie,iskp
                            do j=js,je,jskp
                               call setcorner(j,kdim,i,w(lx),w(ly),
     .                                        w(lz),w(lx),w(ly),w(lz),
     .                                        jdim,kdim,idim,wk,nsurf,1,
     .                                        ivert,slavept,nslave,nou,
     .                                        bou,ibufdim,nbuf,myid,
     .                                        iwk,nmaster,iseq)
                            end do
                         end do
                     end if
                  end do
c
               end if
c
            end if
         end if
         end do 
      end if
c
      end do
c
      return
      end
