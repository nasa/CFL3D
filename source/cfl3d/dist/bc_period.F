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
      subroutine bc_period(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,maxgr,
     .                     maxseg,iadvance,bcfilei,bcfilej,bcfilek,
     .                     lwdat,xorig,yorig,zorig,jdimg,kdimg,idimg,
     .                     lbcprd,isav_prd,
     .                     period_miss,epsrot,ireq_ar,index_ar,
     .                     ireq_snd,keep_trac,keep_trac2,myid,myhost,
     .                     mycomm,mblk2nd,nou,bou,nbuf,ibufdim,istat2,
     .                     istat_size,bcfiles,mxbcfil,nummem)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Update periodic boundary conditions.
c      note: certain arrays, e.g. ireq_ar(lbcprd*5) are dimensioned
c            large enough to allow future use of boundary volume
c            data, as is done in bc_blkint
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
#   ifdef BUILD_MPE
#     include "mpef.h"
#   endif
      dimension mp(4)
#endif
c
      character*80 filname
      character*80  bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
c
      integer bcfilei,bcfilej,bcfilek
c
      dimension nou(nbuf)
      dimension istat2(istat_size,lbcprd*5)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl),
     .          lwdat(maxbl,maxseg,6),jdimg(maxbl),kdimg(maxbl),
     .          idimg(maxbl),iadvance(maxbl),mblk2nd(maxbl)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension period_miss(lbcprd),ireq_ar(lbcprd*5),
     .          index_ar(lbcprd*5),ireq_snd(lbcprd*5),
     .          keep_trac(lbcprd,10),keep_trac2(lbcprd*5)
      dimension isav_prd(lbcprd,12)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbll,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /sklton/ isklton
      common /maxiv/ ivmx
      common /zero/ iexp
      common /is_perbc/ is_prd(5),ie_prd(5),nbcprd
c
c     epsrot...tolerence on 1-1 geometric mismatch after rotation of
c              periodic face (10**(-iexp) is machine zero)
c              if mismatch < epsrot, no warning
c              if mismatch > epsrot, print warning but do not stop!!!!
c
      epsrot  = max(1.e-09,10.**(-iexp+1))
c
      if (ntime.gt.0 .and. nbcprd.gt.0) then 
c
c***********************************************************************
c        First Case: all data needed to set periodic bc lies on the
c                    current processor
c***********************************************************************
c
         iwk_indx = 1
c
         do lcnt = is_prd(level),ie_prd(level)
c           nbll is current (to) block
c           nblp is periodic (from) block
            nbll = isav_prd(lcnt,1)
            nblp = isav_prd(lcnt,12)
            nd_recv = mblk2nd(nbll)
            nd_srce = mblk2nd(nblp)
c
            if (nd_srce.eq.myid .and. nd_recv.eq.myid) then
c
               if (iadvance(nbll).ge.0) then
c
                  call lead(nbll,lw,lw2,maxbl)
c
                  n      = lcnt
                  nface  = isav_prd(lcnt,2)
                  ista   = isav_prd(lcnt,3)
                  iend   = isav_prd(lcnt,4)
                  jsta   = isav_prd(lcnt,5)
                  jend   = isav_prd(lcnt,6)
                  ksta   = isav_prd(lcnt,7)
                  kend   = isav_prd(lcnt,8)
                  mdim   = isav_prd(lcnt,9)
                  ndim   = isav_prd(lcnt,10)
                  nseg   = isav_prd(lcnt,11)
                  ldata  = lwdat(nbll,nseg,nface)
                  if (nface.eq.1) then
                     filname = bcfiles(bcfilei(nbll,nseg,1))
                  else if (nface.eq.2) then
                     filname = bcfiles(bcfilei(nbll,nseg,2))
                  else if (nface.eq.3) then
                     filname = bcfiles(bcfilej(nbll,nseg,1))
                  else if (nface.eq.4) then
                     filname = bcfiles(bcfilej(nbll,nseg,2))
                  else if (nface.eq.5) then
                     filname = bcfiles(bcfilek(nbll,nseg,1))
                  else if (nface.eq.6) then
                     filname = bcfiles(bcfilek(nbll,nseg,2))
                  end if
                  idimp = idimg(nblp)
                  jdimp = jdimg(nblp)
                  kdimp = kdimg(nblp)
                  if (nface.eq.1 .or. nface.eq.2) maxdims = jdimp*kdimp
                  if (nface.eq.3 .or. nface.eq.4) maxdims = kdimp*idimp
                  if (nface.eq.5 .or. nface.eq.6) maxdims = jdimp*idimp
c
c                 k = constant interface
c
                  if (nface.eq.5 .or. nface.eq.6) then
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*9 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for chkrot'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        call chkrot(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                              w(lz),nblp,jdimp,kdimp,idimp,
     .                              w(lw(10,nblp)),w(lw(11,nblp)),
     .                              w(lw(12,nblp)),nface,w(ldata),
     .                              wk(iwk_indx),wk(iwk_indx+maxdims*3),
     .                              wk(iwk_indx+maxdims*6),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              mdim,ndim,lcnt,xorig,yorig,zorig,
     .                              maxbl,period_miss,lbcprd,nou,bou,
     .                              nbuf,ibufdim,myid)
                     end if
c
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*36 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for bctype 2005'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                     end if
                     iskltsav = isklton
                     isklton  = 0
                     call bc2005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .                           w(lqi0),ista,iend,jsta,jend,ksta,kend,
     .                           nface,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .                           w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,
     .                           ndim,w(ldata),filname,
     .                           w(lw(1,nblp)),w(lw(13,nblp)),
     .                           w(lw(19,nblp)),jdimp,kdimp,idimp,
     .                           wk(iwk_indx),wk(iwk_indx+maxdims*10),
     .                           wk(iwk_indx+maxdims*20),nbll,nblp,
     .                           nou,bou,nbuf,ibufdim,myid,mblk2nd,
     .                           maxbl,nummem)
                     isklton = iskltsav
                  end if
c
c                 j = constant interface
c
                  if (nface.eq.3 .or. nface.eq.4) then
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*9 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for chkrot'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        call chkrot(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                              w(lz),nblp,jdimp,kdimp,idimp,
     .                              w(lw(10,nblp)),w(lw(11,nblp)),
     .                              w(lw(12,nblp)),nface,w(ldata),
     .                              wk(iwk_indx),wk(iwk_indx+maxdims*3),
     .                              wk(iwk_indx+maxdims*6),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              mdim,ndim,lcnt,xorig,yorig,zorig,
     .                              maxbl,period_miss,lbcprd,nou,bou,
     .                              nbuf,ibufdim,myid)
                     end if
c
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*36 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for bctype 2005'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                     end if
                     iskltsav = isklton
                     isklton  = 0
                     call bc2005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .                           w(lqi0),ista,iend,jsta,jend,ksta,kend,
     .                           nface,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .                           w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,
     .                           ndim,w(ldata),filname,
     .                           w(lw(1,nblp)),w(lw(13,nblp)),
     .                           w(lw(19,nblp)),jdimp,kdimp,idimp,
     .                           wk(iwk_indx),wk(iwk_indx+maxdims*10),
     .                           wk(iwk_indx+maxdims*20),nbll,nblp,
     .                           nou,bou,nbuf,ibufdim,myid,mblk2nd,
     .                           maxbl,nummem)
                     isklton = iskltsav
                  end if
c
c                 i = constant interface
c
                  if (nface.eq.1 .or. nface.eq.2) then
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*9 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for chkrot'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        call chkrot(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                              w(lz),nblp,jdimp,kdimp,idimp,
     .                              w(lw(10,nblp)),w(lw(11,nblp)),
     .                              w(lw(12,nblp)),nface,w(ldata),
     .                              wk(iwk_indx),wk(iwk_indx+maxdims*3),
     .                              wk(iwk_indx+maxdims*6),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              mdim,ndim,lcnt,xorig,yorig,zorig,
     .                              maxbl,period_miss,lbcprd,nou,bou,
     .                              nbuf,ibufdim,myid)
                     end if
c
                     if (isklton.eq.1) then
                        kcheck = iwk_indx + maxdims*36 - 1
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)' stop in bc_period',
     .                     '...insufficient wk storage for bctype 2005'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                     end if
                     iskltsav = isklton
                     isklton  = 0
                     call bc2005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .                           w(lqi0),ista,iend,jsta,jend,ksta,kend,
     .                           nface,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .                           w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,
     .                           ndim,w(ldata),filname,
     .                           w(lw(1,nblp)),w(lw(13,nblp)),
     .                           w(lw(19,nblp)),jdimp,kdimp,idimp,
     .                           wk(iwk_indx),wk(iwk_indx+maxdims*10),
     .                           wk(iwk_indx+maxdims*20),nbll,nblp,
     .                           nou,bou,nbuf,ibufdim,myid,mblk2nd,
     .                           maxbl,nummem)
                     isklton = iskltsav
                  end if
c
               end if
c
            end if
c
         end do
c
c***********************************************************************
c        Second Case: data needed to set periodic bc lies on another
c                     processor
c***********************************************************************
#if defined DIST_MPI
#        ifdef BUILD_MPE
c
c        begin monitoring message passing
c
         call MPE_Log_event (30, 0, "Start BC_PERIOD")
#        endif
c
c        set baseline tag values
c
         ioffset = lbcprd
         itag_x  = 1
         itag_q  = itag_x + ioffset
         itag_v  = itag_q + ioffset
         itag_t  = itag_v + ioffset
c
c        post a bunch of receives first (for non-buffering implementations)
c        set the request index and index for wk
c
         kqintl = 1
         ireq   = 0
c
         do lcnt = is_prd(level),ie_prd(level)
c           nbll is current (to) block
c           nblp is periodic (from) block
            nbll = isav_prd(lcnt,1)
            if (iadvance(nbll).ge.0) then
            nblp    = isav_prd(lcnt,12)
            nd_recv = mblk2nd(nbll)
            nd_srce = mblk2nd(nblp)
            if (nd_recv.eq.myid) then
               if (nd_srce.ne.myid) then
                  n     = lcnt
                  nseg  = isav_prd(lcnt,11)
                  nface = isav_prd(lcnt,2)
                  ldata = lwdat(nbll,nseg,nface)
                  idimp = idimg(nblp)
                  jdimp = jdimg(nblp)
                  kdimp = kdimg(nblp)
                  if (nface.eq.1 .or. nface.eq.2) maxdims = jdimp*kdimp
                  if (nface.eq.3 .or. nface.eq.4) maxdims = kdimp*idimp
                  if (nface.eq.5 .or. nface.eq.6) maxdims = jdimp*idimp
c
c                 receive xyz data for checking mismatch
c 
                  if (isklton.eq.1) then
                     ldim = 3
                     np   = 1
                     kcheck = kqintl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_x + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,1)   = kqintl
                     keep_trac(n,2)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
c
c                 receive q data
c
                  ldim = 5
                  np   = 2
                  kcheck = kqintl + maxdims*ldim*np
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stop in bc_period',
     .               '....work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  mytag = itag_q + n
                  ireq  = ireq + 1
                  call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                           MY_MPI_REAL,
     .                           nd_srce,mytag,mycomm,
     .                           ireq_ar(ireq),ierr)
                  keep_trac(n,3)   = kqintl
                  keep_trac(n,4)   = ireq
                  keep_trac2(ireq) = lcnt
                  kqintl = kcheck
c
c                 receive vist3d data
c
                  if (ivmx.ge.2) then
                     ldim = 1
                     np   = 2
                     kcheck = kqintl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_v + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,5)   = kqintl
                     keep_trac(n,6)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
c
c                 receive turb. data 
c
                  if (ivmx.ge.4) then
                     ldim = nummem
                     np   = 2
                     kcheck = kqintl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_t + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,7)   = kqintl
                     keep_trac(n,8)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
               end if
            end if 
            end if
         end do
c
         if (myid.ne.myhost) then
            if (ireq.gt.lbcprd*4) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),999) ireq,lbcprd*4
 999           format(' problem in bc_period...ireq = ',i4,
     .         ' but max allowable value = lbcprd*4 = ',i4)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
        end if
c
c        loop over all periodic boundaries looking for blocks that
c        need to send out info to other processors
c
         ktl   = kqintl
         ireq2 = 0
c
         do lcnt = is_prd(level),ie_prd(level)
c           nbll is current (to) block
c           nblp is periodic (from) block
            nbll    = isav_prd(lcnt,1)
            if (iadvance(nbll).ge.0) then
            nblp    = isav_prd(lcnt,12)
            nd_recv = mblk2nd(nbll)
            nd_srce = mblk2nd(nblp)
            if (nd_srce.eq.myid) then
               if (nd_recv.ne.myid) then
                  n      = lcnt
                  nseg   = isav_prd(lcnt,11)
                  nface  = isav_prd(lcnt,2)
                  ldata  = lwdat(nbll,nseg,nface)
                  nface  = isav_prd(lcnt,2)
                  ista   = isav_prd(lcnt,3)
                  iend   = isav_prd(lcnt,4)
                  jsta   = isav_prd(lcnt,5)
                  jend   = isav_prd(lcnt,6)
                  ksta   = isav_prd(lcnt,7)
                  kend   = isav_prd(lcnt,8)
                  idimp  = idimg(nblp)
                  jdimp  = jdimg(nblp)
                  kdimp  = kdimg(nblp)
                  if (nface.eq.1 .or. nface.eq.2) maxdims = jdimp*kdimp
                  if (nface.eq.3 .or. nface.eq.4) maxdims = kdimp*idimp
                  if (nface.eq.5 .or. nface.eq.6) maxdims = jdimp*idimp
c
                  if (isklton.eq.1) then
c
c                    set up mp array for one plane of grid-point data; 
c                    mp indicates the planes to be loaded into the work
c                    array for transfer to another processor
c
                     np   = 1
                     if (nface.eq.1) then
                        mp(1) = idimp
                     else if (nface.eq.2) then
                        mp(1) = 1
                     else if (nface.eq.3) then
                        mp(1) = jdimp
                     else if (nface.eq.4) then
                        mp(1) = 1
                     else if (nface.eq.5) then
                        mp(1) = kdimp
                     else if (nface.eq.6) then
                        mp(1) = 1
                     end if
c
c                    load 1 plane of xyz data from full 3D periodic block
c                    to a work array and send to the appropriate processor
c
                     lwxr = lw(10,nblp)
                     ldim = 3
                     kcheck = ktl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (nface.eq.1 .or. nface.eq.2) then
                        call ld_dati(w(lwxr),jdimp,kdimp,idimp,wk(ktl),
     .                               ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                               nou,bou,nbuf,ibufdim,myid)
                     else if (nface.eq.3 .or. nface.eq.4) then
                        call ld_datj(w(lwxr),jdimp,kdimp,idimp,wk(ktl),
     .                               ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                               nou,bou,nbuf,ibufdim,myid)
                     else
                        call ld_datk(w(lwxr),jdimp,kdimp,idimp,wk(ktl),
     .                               ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                               nou,bou,nbuf,ibufdim,myid)
                     end if
                     mytag = itag_x + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*np,
     .                             MY_MPI_REAL,
     .                             nd_recv, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
c
c                 set up mp array for 2 planes of cell-center data;
c                 mp indicates the planes to be loaded into the work
c                 array for transfer to another processor
c
                  np   = 2
                  if (nface.eq.1) then
                     mp(1) = idimp-1
                     mp(2) = max(idimp-2,1)
                  else if (nface.eq.2) then
                     mp(1) = 1
                     mp(2) = min(idimp-1,2)
                  else if (nface.eq.3) then
                     mp(1) = jdimp-1
                     mp(2) = max(jdimp-2,1)
                  else if (nface.eq.4) then
                     mp(1) = 1
                     mp(2) = min(jdimp-1,2)
                  else if (nface.eq.5) then
                     mp(1) = kdimp-1
                     mp(2) = max(kdimp-2,1)
                  else if (nface.eq.6) then
                     mp(1) = 1
                     mp(2) = min(kdimp-1,2)
                  end if
c
c                 load 2 planes of q data from full 3D periodic block
c                 to a work array and send to the appropriate processor
c
                  lws  = lw( 1,nblp)
                  ldim = 5
                  kcheck = ktl + maxdims*ldim*np
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stop in bc_period',
     .               '....work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lws),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                  else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lws),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                  else
                     call ld_datk(w(lws),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                  end if
                  mytag = itag_q + n
                  ireq2 = ireq2 + 1
                  call MPI_ISend(wk(ktl), maxdims*ldim*np,
     .                          MY_MPI_REAL,
     .                          nd_recv, mytag, mycomm,
     .                          ireq_snd(ireq2), ierr)
                  ktl = kcheck
c
c                 load 2 planes of vist3d data from full 3D periodic block
c                 to a work array and send to the appropriate processor
c
                  if (ivmx.ge.2) then
                     lwst = lw(13,nblp)
                     ldim = 1
                     kcheck = ktl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else
                     call ld_datk(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     end if
                     mytag = itag_v + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*np,
     .                             MY_MPI_REAL,
     .                             nd_recv, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
c
c                 load 2 planes of turb data from full 3D periodic block
c                 to a work array and send to the appropriate processor
c
                  if (ivmx.ge.4) then
                     lwst = lw(19,nblp)
                     ldim = nummem
                     kcheck = ktl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stop in bc_period',
     .                  '....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else
                     call ld_datk(w(lwst),jdimp,kdimp,idimp,wk(ktl),
     .                            ldim,mp,np,1,idimp,1,jdimp,1,kdimp,
     .                            nou,bou,nbuf,ibufdim,myid)
                     end if
                     mytag = itag_t + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*np,
     .                             MY_MPI_REAL,
     .                             nd_recv, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
               end if
            end if
            end if
         end do
c
c        set periodic bc's 
c
         ndone  = 0
c
         do while (ndone.lt.ireq)
c
         call MPI_Waitsome(ireq,ireq_ar,nrecvd,index_ar,
     .   istat2,ierr)
c
         if (nrecvd.gt.0) then
            ndone = ndone + nrecvd
            do nnn=1,nrecvd
               lcnt    = keep_trac2(index_ar(nnn))
               n       = lcnt
               nbll    = isav_prd(lcnt,1)
               nd_recv = mblk2nd(nbll)
c
               if (iadvance(nbll).ge.0) then
c
               call lead(nbll,lw,lw2,maxbl)
c
               nface  = isav_prd(lcnt,2)
               ista   = isav_prd(lcnt,3)
               iend   = isav_prd(lcnt,4)
               jsta   = isav_prd(lcnt,5)
               jend   = isav_prd(lcnt,6)
               ksta   = isav_prd(lcnt,7)
               kend   = isav_prd(lcnt,8)
               mdim   = isav_prd(lcnt,9)
               ndim   = isav_prd(lcnt,10)
               nseg   = isav_prd(lcnt,11)
               nblp   = isav_prd(lcnt,12)
               nd_srce = mblk2nd(nblp)
               ldata  = lwdat(nbll,nseg,nface)
               if (nface.eq.1) then
                  filname = bcfiles(bcfilei(nbll,nseg,1))
               else if (nface.eq.2) then
                  filname = bcfiles(bcfilei(nbll,nseg,2))
               else if (nface.eq.3) then
                  filname = bcfiles(bcfilej(nbll,nseg,1))
               else if (nface.eq.4) then
                  filname = bcfiles(bcfilej(nbll,nseg,2))
               else if (nface.eq.5) then
                  filname = bcfiles(bcfilek(nbll,nseg,1))
               else if (nface.eq.6) then
                  filname = bcfiles(bcfilek(nbll,nseg,2))
               end if
               idimp = idimg(nblp)
               jdimp = jdimg(nblp)
               kdimp = kdimg(nblp)
               if (nface.eq.1 .or. nface.eq.2) maxdims = jdimp*kdimp
               if (nface.eq.3 .or. nface.eq.4) maxdims = kdimp*idimp
               if (nface.eq.5 .or. nface.eq.6) maxdims = jdimp*idimp
c
c              k = constant interface
c
               if (nface.eq.5 .or. nface.eq.6) then
                  if (isklton.eq.1) then
c
c                    check periodic mismatch
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call chkrotk_d(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                                 w(lz),nblp,jdimp,kdimp,idimp,
     .                                 nface,w(ldata),wk(kqintl),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 mdim,ndim,lcnt,xorig,yorig,zorig,
     .                                 maxbl,period_miss,lbcprd,nou,bou,
     .                                 nbuf,ibufdim,myid)
                     end if
                  end if
c
c
c                 set q data
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     iflag=0
                     call bc2005k_d(jdim,kdim,idim,w(lqk0),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              nface,mdim,ndim,w(ldata),
     .                              filname,jdimp,kdimp,idimp,
     .                              wk(kqintl),nbll,nblp,ldim,nou,
     .                              bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .                              iflag)
                  end if
c
c                 set vist3d data
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        kqintl = keep_trac(n,5)
                        iflag=1
                        call bc2005k_d(jdim,kdim,idim,w(lvk0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
c
c                 set turb. data
c
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        kqintl = keep_trac(n,7)
                        iflag=2
                        call bc2005k_d(jdim,kdim,idim,w(ltk0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
               end if
c
c              j = constant interface
c
               if (nface.eq.3 .or. nface.eq.4) then
                  if (isklton.eq.1) then
c
c                    check periodic mismatch
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call chkrotj_d(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                                 w(lz),nblp,jdimp,kdimp,idimp,
     .                                 nface,w(ldata),wk(kqintl),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 mdim,ndim,lcnt,xorig,yorig,zorig,
     .                                 maxbl,period_miss,lbcprd,nou,bou,
     .                                 nbuf,ibufdim,myid)
                     end if
                  end if
c
c                 set q data
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     iflag=0
                     call bc2005j_d(jdim,kdim,idim,w(lqj0),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              nface,mdim,ndim,w(ldata),
     .                              filname,jdimp,kdimp,idimp,
     .                              wk(kqintl),nbll,nblp,ldim,nou,
     .                              bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .                              iflag)
                  end if
c
c                 set vist3d data
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        kqintl = keep_trac(n,5)
                        iflag=1
                        call bc2005j_d(jdim,kdim,idim,w(lvj0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
c
c                 set turb. data
c
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        kqintl = keep_trac(n,7)
                        iflag=2
                        call bc2005j_d(jdim,kdim,idim,w(ltj0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
               end if
c
c              i = constant interface
c
               if (nface.eq.1 .or. nface.eq.2) then
                  if (isklton.eq.1) then
c
c                    check periodic mismatch
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call chkroti_d(nbll,jdim,kdim,idim,w(lx),w(ly),
     .                                 w(lz),nblp,jdimp,kdimp,idimp,
     .                                 nface,w(ldata),wk(kqintl),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 mdim,ndim,lcnt,xorig,yorig,zorig,
     .                                 maxbl,period_miss,lbcprd,nou,bou,
     .                                 nbuf,ibufdim,myid)
                     end if
                  end if
c
c                 set q data
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     iflag=0
                     call bc2005i_d(jdim,kdim,idim,w(lqi0),
     .                              ista,iend,jsta,jend,ksta,kend,
     .                              nface,mdim,ndim,w(ldata),
     .                              filname,jdimp,kdimp,idimp,
     .                              wk(kqintl),nbll,nblp,ldim,nou,
     .                              bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .                              iflag)
                  end if
c
c                 set vist3d data
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        kqintl = keep_trac(n,5)
                        iflag=1
                        call bc2005i_d(jdim,kdim,idim,w(lvi0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
c
c                 set turb. data
c
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        kqintl = keep_trac(n,7)
                        iflag=2
                        call bc2005i_d(jdim,kdim,idim,w(lti0),
     .                                 ista,iend,jsta,jend,ksta,kend,
     .                                 nface,mdim,ndim,w(ldata),
     .                                 filname,jdimp,kdimp,idimp,
     .                                 wk(kqintl),nbll,nblp,ldim,nou,
     .                                 bou,nbuf,ibufdim,myid,mblk2nd,
     .                                 maxbl,iflag)
                     end if
                  end if
               end if
c
               end if
c
            end do
c
         end if
c
         end do
c
c        make sure all sends are completed before exiting
c
         if (ireq2.gt.0) then
            call MPI_Waitall (ireq2, ireq_snd, istat2, ierr)
         end if
c
#        ifdef BUILD_MPE
c        end monitoring message passing
c
         call MPE_Log_event (39, 0, "End BC_PERIOD")
#        endif
#endif
c
      end if
c
      return
      end
