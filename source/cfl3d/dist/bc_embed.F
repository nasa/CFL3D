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
      subroutine bc_embed(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,maxgr,
     .                    lbcemb,iadvance,idimg,jdimg,kdimg,
     .                    isav_emb,ireq_ar,ireq_snd,
     .                    index_ar,keep_trac,keep_trac2,myid,myhost,
     .                    mycomm,mblk2nd,nou,bou,nbuf,ibufdim,iviscg,
     .                    istat2,istat_size,nummem)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Update embedded grid boundaries
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
      dimension mp(4)
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension istat2(istat_size,lbcemb*3)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl),iadvance(maxbl),
     .          mblk2nd(maxbl),iviscg(maxbl,3)
      dimension isav_emb(lbcemb,12)
      dimension ireq_ar(lbcemb*3),index_ar(lbcemb*3),
     .          ireq_snd(lbcemb*3),keep_trac(lbcemb,6),
     .          keep_trac2(lbcemb*3)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /maxiv/ ivmx
      common /is_embedbc/ is_emb(5),ie_emb(5),nbcemb
c
      if (ie_emb(level).lt.is_emb(level)) return
c
      if (ntime.gt.0 .and. nbcemb.gt.0) then
c
c***********************************************************************
c        First Case: all data needed to set embeded bc lies on the
c                    current processor
c***********************************************************************
c
         do lcnt = is_emb(level),ie_emb(level)
c           nblf is finer (embeded) block
c           nblc is coarser block in which block nblf is embedded
            nblf = isav_emb(lcnt,1)
            nblc = isav_emb(lcnt,9)         
            nd_recv = mblk2nd(nblf)
            nd_srce = mblk2nd(nblc)
c
            if (nd_srce.eq.myid .and. nd_recv.eq.myid) then
c
               if (iadvance(nblf).ge.0) then
c
               call lead(nblf,lw,lw2,maxbl)
c
               n       = lcnt
               nface   = isav_emb(lcnt,2)
               is      = isav_emb(lcnt,3)
               ie      = isav_emb(lcnt,4)
               js      = isav_emb(lcnt,5)
               je      = isav_emb(lcnt,6)
               ks      = isav_emb(lcnt,7)
               ke      = isav_emb(lcnt,8)
               nblc    = isav_emb(lcnt,9)
               nsi     = isav_emb(lcnt,10)
               nd_srce = mblk2nd(nblc)
               idimc   = idimg(nblc)
               jdimc   = jdimg(nblc)
               kdimc   = kdimg(nblc)
               if (nface.eq.1 .or. nface.eq.2) maxdims = jdimc*kdimc
               if (nface.eq.3 .or. nface.eq.4) maxdims = kdimc*idimc
               if (nface.eq.5 .or. nface.eq.6) maxdims = jdimc*idimc
c
               lqc  = lw(1,nblc)
               lqct = lw(19,nblc)
               lqcv = lw(13,nblc)
c
               if (nsi.eq.2) then
c
c                 full coarsening in i-direction
c
c                 interpolate q
c
                  ldim = 5
                  call i2x(jdimc,kdimc,idimc,w(lqc),jdim,kdim,
     .                     idim,w(lqj0),w(lqk0),w(lqi0),
     .                     js,ks,is,je,ke,ie,nblc,ldim,
     .                     nblf,w(lbcj),w(lbck),w(lbci),nface)
c
c                 interpolate vist3d
c
                  if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                .or. iviscg(nblc,3).ge.2) then
                     ldim = 1
                     call i2x(jdimc,kdimc,idimc,w(lqcv),jdim,kdim,
     .                        idim,w(lvj0),w(lvk0),w(lvi0),
     .                        js,ks,is,je,ke,ie,nblc,ldim,
     .                        nblf,w(lbcj),w(lbck),w(lbci),nface)
                  end if
c
c                 interpolate turb. data
c
                  if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                .or. iviscg(nblc,3).ge.4) then
                     ldim = nummem
                     call i2x(jdimc,kdimc,idimc,w(lqct),jdim,kdim,
     .                        idim,w(lqj0),w(lqk0),w(lqi0),
     .                        js,ks,is,je,ke,ie,nblc,ldim,
     .                        nblf,w(lbcj),w(lbck),w(lbci),nface)
                  end if
c
               else if (nsi.eq.1) then
c
c                 semi coarsening in i-direction
c
c                 interpolate q
c
                  ldim = 5
                  call i2xs(jdimc,kdimc,idimc,w(lqc),jdim,kdim,
     .                      idim,w(lqj0),w(lqk0),w(lqi0),
     .                      js,ks,is,je,ke,ie,nblc,ldim,
     .                      nblf,w(lbcj),w(lbck),w(lbci),nface)
c
c                 interpolate vist3d
c
                  if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                .or. iviscg(nblc,3).ge.2) then
                     ldim = 1
                     call i2xs(jdimc,kdimc,idimc,w(lqcv),jdim,kdim,
     .                         idim,w(lvj0),w(lvk0),w(lvi0),
     .                         js,ks,is,je,ke,ie,nblc,ldim,
     .                         nblf,w(lbcj),w(lbck),w(lbci),nface)
                  end if
c
c                 interpolate turb. data
c
                  if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                .or. iviscg(nblc,3).ge.4) then
                     ldim = nummem
                     call i2xs(jdimc,kdimc,idimc,w(lqct),jdim,kdim,
     .                         idim,w(lqj0),w(lqk0),w(lqi0),
     .                         js,ks,is,je,ke,ie,nblc,ldim,
     .                         nblf,w(lbcj),w(lbck),w(lbci),nface)
                  end if
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
         call MPE_Log_event (40, 0, "Start BC_EMBED")
#        endif
c
c        set baseline tag values
c
         ioffset = lbcemb
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
         do lcnt = is_emb(level),ie_emb(level)
c           nblf is finer (embeded) block
c           nblc is coarser block in which block nblf is embedded
            nblf    = isav_emb(lcnt,1)
            if (iadvance(nblf).ge.0) then
            nface   = isav_emb(lcnt,2)
            nblc    = isav_emb(lcnt,9)
            nd_recv = mblk2nd(nblf)
            nd_srce = mblk2nd(nblc)
            if (nd_recv.eq.myid) then
               if (nd_srce.ne.myid) then
                  n     = lcnt
                  idimc = idimg(nblc)
                  jdimc = jdimg(nblc)
                  kdimc = kdimg(nblc)
                  if (nface.eq.1 .or. nface.eq.2) maxdims = jdimc*kdimc
                  if (nface.eq.3 .or. nface.eq.4) maxdims = kdimc*idimc
                  if (nface.eq.5 .or. nface.eq.6) maxdims = jdimc*idimc
c
c                 receive q data
c
                  ldim = 5
                  np   = 3
                  kcheck = kqintl + maxdims*ldim*np
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in ',
     .               'bc_embed....work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  mytag = itag_q + n
                  ireq  = ireq + 1
                  call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                           MY_MPI_REAL,
     .                           nd_srce,mytag,mycomm,
     .                           ireq_ar(ireq),ierr)
                  keep_trac(n,1)  = kqintl
                  keep_trac(n,2)  = ireq
                  keep_trac2(ireq) = lcnt
                  kqintl = kcheck
c
c                 receive vist3d data
c
                  if (ivmx.ge.2) then
                     ldim = 1
                     np   = 3
                     kcheck = kqintl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in ',
     .                  'bc_embed....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_v + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,3)  = kqintl
                     keep_trac(n,4)  = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
c
c                 receive turb. data
c
                  if (ivmx.ge.4) then
                     ldim = nummem
                     np   = 3
                     kcheck = kqintl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in ',
     .                  'bc_embed....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_t + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*np,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,5)  = kqintl
                     keep_trac(n,6)  = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
               end if
            end if
            end if
         end do
c
         if (myid.ne.myhost) then
            if (ireq.gt.lbcemb*4) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),999) ireq,lbcemb*4
 999           format(' problem in bc_embed...ireq = ',i4,
     .         ' but max allowable value = lbcemb*4 = ',i4)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
        end if
c
c        loop over all embedded boundaries looking for blocks that
c        need to send out info to other processors
c
         ktl   = kqintl
         ireq2 = 0
c
         do lcnt = is_emb(level),ie_emb(level)
c           nblf is finer (embeded) block
c           nblc is coarser block in which block nblf is embedded
            nblf    = isav_emb(lcnt,1)
            if (iadvance(nblf).ge.0) then
            nblc    = isav_emb(lcnt,9)
            nface   = isav_emb(lcnt,2)
            nd_recv = mblk2nd(nblf)
            nd_srce = mblk2nd(nblc)
            if (nd_srce.eq.myid) then
               if (nd_recv.ne.myid) then
                  n      = lcnt
                  is     = isav_emb(lcnt,3)
                  ie     = isav_emb(lcnt,4)
                  js     = isav_emb(lcnt,5)
                  je     = isav_emb(lcnt,6)
                  ks     = isav_emb(lcnt,7)
                  ke     = isav_emb(lcnt,8)
                  idimc  = idimg(nblc)
                  jdimc  = jdimg(nblc)
                  kdimc  = kdimg(nblc)
                  if (nface.eq.1 .or. nface.eq.2) maxdims = jdimc*kdimc
                  if (nface.eq.3 .or. nface.eq.4) maxdims = kdimc*idimc
                  if (nface.eq.5 .or. nface.eq.6) maxdims = jdimc*idimc
c
c                 set up mp array for 3 planes of cell-center data;
c                 mp indicates the planes to be loaded into the work
c                 array for transfer to another processor
c
                  np   = 3
                  if (nface.eq.1) then
                     mp(1) = is
                     mp(2) = is-1
                     mp(3) = is-2
                  else if (nface.eq.2) then
                     mp(1) = ie-1
                     mp(2) = ie
                     mp(3) = ie+1
                  else if (nface.eq.3) then
                     mp(1) = js
                     mp(2) = js-1
                     mp(3) = js-2
                  else if (nface.eq.4) then
                     mp(1) = je-1
                     mp(2) = je
                     mp(3) = je+1
                  else if (nface.eq.5) then
                     mp(1) = ks
                     mp(2) = ks-1
                     mp(3) = ks-2
                  else if (nface.eq.6) then
                     mp(1) = ke-1
                     mp(2) = ke
                     mp(3) = ke+1
                  end if
c
c                 load 3 planes of q data from full 3D embeded block
c                 to a work array and send to the appropriate processor
c
                  lws  = lw( 1,nblc)
                  ldim = 5
                  kcheck = ktl + maxdims*ldim*np
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in ',
     .               'bc_embed....work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lws),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                  else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lws),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                  else
                     call ld_datk(w(lws),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
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
c                 load 3 planes of vist3d data from full 3D embeded block
c                 to a work array and send to the appropriate processor
c
                  if (ivmx.ge.2) then
                     lwst = lw(13,nblc)
                     ldim = 1
                     kcheck = ktl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in ',
     .                  'bc_embed....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else
                     call ld_datk(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
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
c                 load 3 planes of turb data from full 3D embeded block
c                 to a work array and send to the appropriate processor
c
                  if (ivmx.ge.4) then
                     lwst = lw(19,nblc)
                     ldim = nummem
                     kcheck = ktl + maxdims*ldim*np
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in ',
     .                  'bc_embed....work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (nface.eq.1 .or. nface.eq.2) then
                     call ld_dati(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else if (nface.eq.3 .or. nface.eq.4) then
                     call ld_datj(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
     .                            nou,bou,nbuf,ibufdim,myid)
                     else
                     call ld_datk(w(lwst),jdimc,kdimc,idimc,wk(ktl),
     .                            ldim,mp,np,1,idimc,1,jdimc,1,kdimc,
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
c        set embeded-grid bcs
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
               nblf = isav_emb(lcnt,1)
               nd_recv = mblk2nd(nblf)
c
               if (iadvance(nblf).ge.0) then
c
               call lead(nblf,lw,lw2,maxbl)
c
               nface   = isav_emb(lcnt,2)
               is      = isav_emb(lcnt,3)
               ie      = isav_emb(lcnt,4)
               js      = isav_emb(lcnt,5)
               je      = isav_emb(lcnt,6)
               ks      = isav_emb(lcnt,7)
               ke      = isav_emb(lcnt,8)
               nblc    = isav_emb(lcnt,9)
               nsi     = isav_emb(lcnt,10)
               nd_srce = mblk2nd(nblc)
               idimc   = idimg(nblc)
               jdimc   = jdimg(nblc)
               kdimc   = kdimg(nblc)
               if (nface.eq.1 .or. nface.eq.2) maxdims = jdimc*kdimc
               if (nface.eq.3 .or. nface.eq.4) maxdims = kdimc*idimc
               if (nface.eq.5 .or. nface.eq.6) maxdims = jdimc*idimc
c
c              k = constant interface
c
               if (nface.eq.5 .or. nface.eq.6) then
c
                  if (nsi.eq.2) then
c
c                    full coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)
                        call i2xk_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                              jdim,kdim,idim,w(lqk0),js,ks,is,
     .                              je,ke,ie,nblc,ldim,nbl,w(lbck),
     .                              nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xk_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                 jdim,kdim,idim,w(lvk0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbck),
     .                                 nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                   .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xk_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                 jdim,kdim,idim,w(lqk0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbck),
     .                                 nface)
                        end if
                     end if
c
                  else if (nsi.eq.1) then
c
c                    semi coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)

                        call i2xsk_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                               jdim,kdim,idim,w(lqk0),js,ks,is,
     .                               je,ke,ie,nblc,ldim,nbl,w(lbck),
     c                               nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xsk_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                  jdim,kdim,idim,w(lvk0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbck),
     .                                  nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                   .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xsk_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                  jdim,kdim,idim,w(lqk0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbck),
     .                                  nface)
                        end if
                     end if
c
                  end if
c
               end if
c
c              j = constant interface
c
               if (nface.eq.3 .or. nface.eq.4) then
c
                  if (nsi.eq.2) then
c
c                    full coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)
                        call i2xj_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                              jdim,kdim,idim,w(lqj0),js,ks,is,
     .                              je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                              nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xj_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                 jdim,kdim,idim,w(lvj0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                                 nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                  .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xj_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                 jdim,kdim,idim,w(lqj0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                                 nface)
                        end if
                     end if
c
                  else if (nsi.eq.1) then
c
c                    semi coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)
                        call i2xsj_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                               jdim,kdim,idim,w(lqj0),js,ks,is,
     .                               je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                               nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xsj_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                  jdim,kdim,idim,w(lvj0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                                  nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                  .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xsj_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                  jdim,kdim,idim,w(lqj0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbcj),
     .                                  nface)
                        end if
                     end if
c
                  end if
c
               end if
c
c              i = constant interface
c
               if (nface.eq.1 .or. nface.eq.2) then
c
                  if (nsi.eq.2) then
c
c                    full coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)
                        call i2xi_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                              jdim,kdim,idim,w(lqi0),js,ks,is,
     .                              je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                              nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xi_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                 jdim,kdim,idim,w(lvi0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                                 nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                   .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xi_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                 jdim,kdim,idim,w(lqi0),js,ks,is,
     .                                 je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                                 nface)
                        end if
                     end if
c
                  else if (nsi.eq.1) then
c
c                    semi coarsening in i-direction
c
c                    interpolate q
c
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        ldim    = 5
                        kqintlq = keep_trac(n,1)
                        call i2xsi_d(jdimc,kdimc,idimc,wk(kqintlq),
     .                               jdim,kdim,idim,w(lqi0),js,ks,is,
     .                               je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                               nface)
                     end if
c
c                    interpolate vist3d
c
                     if (iviscg(nblc,1).ge.2 .or. iviscg(nblc,2).ge.2
     .                  .or. iviscg(nblc,3).ge.2) then
                        if (index_ar(nnn) .eq. keep_trac(n,4)) then
                           ldim = 1
                           kqintlv = keep_trac(n,3)
                           call i2xsi_d(jdimc,kdimc,idimc,wk(kqintlv),
     .                                  jdim,kdim,idim,w(lvi0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                                  nface)
                        end if
                     end if
c
c                    interpolate turb. data
c
                     if (iviscg(nblc,1).ge.4 .or. iviscg(nblc,2).ge.4
     .                  .or. iviscg(nblc,3).ge.4) then
                        if (index_ar(nnn) .eq. keep_trac(n,6)) then
                           ldim = nummem
                           kqintlt = keep_trac(n,5)
                           call i2xsi_d(jdimc,kdimc,idimc,wk(kqintlt),
     .                                  jdim,kdim,idim,w(lqi0),js,ks,is,
     .                                  je,ke,ie,nblc,ldim,nbl,w(lbci),
     .                                  nface)
                        end if
                     end if
c
                  end if
c
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
         call MPE_Log_event (49, 0, "End BC_EMBED")
#        endif
#endif
c
      end if
c
      return
      end
