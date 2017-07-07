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
      subroutine bc_blkint(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,maxgr,
     .                     mxbli,iadvance,geom_miss,epsc0,nblk,nbli,
     .                     limblk,isva,nblon,jdimg,kdimg,idimg,
     .                     mblk2nd,isav_blk,ireq_ar,
     .                     index_ar,ireq_snd,keep_trac,keep_trac2,
     .                     nou,bou,nbuf,ibufdim,myid,myhost,mycomm,
     .                     istat2,istat_size,nummem)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Update 1-1 block interface boundary conditions.
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
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension istat2(istat_size,mxbli*5)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl),
     .          iadvance(maxbl)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),isva(2,2,mxbli),
     .          nblon(mxbli),geom_miss(2*mxbli),isav_blk(2*mxbli,17)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl),mblk2nd(maxbl)
      dimension ireq_ar(mxbli*5),index_ar(mxbli*5),ireq_snd(mxbli*5),
     .          keep_trac(mxbli,10),keep_trac2(mxbli*5)
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
      common /sklton/ isklton
      common /maxiv/ ivmx
      common /zero/ iexp
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
c
c     epsc0....tolerence on 1-1 geometric mismatch
c              if mismatch < epsc0, no warning
c              if mismatch > epsc0, print warning but do not stop!!!!
c
      epsc0 = max(1.e-09,10.**(-iexp+1))
c
c     ivolint...flag for boundary-volume interpolation
c               = 0, do not set boundary volumes from neighboring block
c               = 1, set boundary volumes from neighboring block
c                 (0 is same as previous versions)
c
      if (ntime.gt.0 .and. abs(nbli).gt.0) then 
c
c***********************************************************************
c        First Case: all data needed to set 1-1 bc lies on the
c                    current processor
c***********************************************************************
c
         do lcnt = is_blk(level),ie_blk(level)
c           ic_blk is current (to) block
c           in_blk is neighbor (from) block
            ic_blk  = isav_blk(lcnt,4)
            in_blk  = isav_blk(lcnt,5)
            nd_dest = mblk2nd(ic_blk)
            nd_srce = mblk2nd(in_blk)
c
            if (nd_srce.eq.myid .and. nd_dest.eq.myid) then
c
               if (iadvance(ic_blk).ge.0) then
c
                  call lead(ic_blk,lw,lw2,maxbl)
c
                  n      = isav_blk(lcnt,1)
                  it     = isav_blk(lcnt,2)
                  ir     = isav_blk(lcnt,3)
                  jface  = isav_blk(lcnt,6)
                  jedge  = isav_blk(lcnt,7)
                  jside  = isav_blk(lcnt,8)
                  lwt    = isav_blk(lcnt,9)
                  iedge  = isav_blk(lcnt,10)
                  iss    = isav_blk(lcnt,11)
                  ise    = isav_blk(lcnt,12)
                  jss    = isav_blk(lcnt,13)
                  jse    = isav_blk(lcnt,14)
                  kss    = isav_blk(lcnt,15)
                  kse    = isav_blk(lcnt,16)
                  iti    = isav_blk(lcnt,17)
c
c                 set dimensions of blocks involved
c
                  idimn = idimg(in_blk)
                  jdimn = jdimg(in_blk)
                  kdimn = kdimg(in_blk)
                  idimc = idimg(ic_blk)
                  jdimc = jdimg(ic_blk)
                  kdimc = kdimg(ic_blk)
c
c                 set pointers for neighboring block
c
                  lws  = lw( 1,in_blk)
                  lwxr = lw(10,in_blk)
                  lwyr = lw(11,in_blk)
                  lwzr = lw(12,in_blk)
c
c                 set pointers for current block (note: lwt has
c                 already been set above, and corresponds to 
c                 qi0, qj0, or qk0, depending on the interface)
c
                  lwxt = lw(10,ic_blk)
                  lwyt = lw(11,ic_blk)
                  lwzt = lw(12,ic_blk)
c
c                 k = constant interface 
c
                  if (jside.eq.3) then
c
c                    check geometric mismatch
c
                     if (isklton.eq.1) then
                        call cblkk(nbli,idimn,jdimn,kdimn,idimc,jdimc,
     .                             kdimc,limblk(1,1,n),isva(1,1,n),
     .                             it,ir,iedge,w(lwxr),w(lwyr),w(lwzr),
     .                             w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                             geom_miss,mxbli)
                     end if
c
c                    interpolate q
c
                     ldim = 5
                     call blockk (w(lws),w(lwt),idimn,jdimn,kdimn,
     .                           idimc,jdimc,limblk(1,1,n),
     .                           isva(1,1,n),it,ir,ldim,
     .                           w(lbck),iedge,0)
c
c                    interpolate turbulent quantities
c
                     if (ivmx.ge.2) then
                        ldim = 1
                        lwst = lw(13,in_blk)
                        lwtt = lw(29,ic_blk)
                        if (limblk(it,3,n).ne.1) then
                           lwtt   = lwtt + jdimc*(idimc-1)*1*2
                        end if
                        call blockk (w(lwst),w(lwtt),idimn,jdimn,
     .                              kdimn,idimc,jdimc,limblk(1,1,n),
     .                              isva(1,1,n),it,ir,ldim,
     .                              w(lbck),iedge,0)
                     end if
                     if (ivmx.ge.4) then
                        ldim = nummem
                        lwst = lw(19,in_blk)    
                        lwtt = lw(24,ic_blk)    
                        if (limblk(it,3,n).ne.1) then
                           lwtt   = lwtt + jdimc*(idimc-1)*2*nummem
                        end if
                        call blockk (w(lwst),w(lwtt),idimn,jdimn,
     .                               kdimn,idimc,jdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbck),iedge,0)
                     end if
c
c                    interpolate cell volumes
c
                     if (ivmx.ge.1 .and. ivolint.gt.0) then
                        ldim   = 1
                        lvol   = lw(8,in_blk)
                        lvolk0 = lw(50,ic_blk)
                        if (limblk(it,3,n).ne.1) then
                           lvolk0 = lvolk0 + jdimc*(idimc-1)*1*2
                        end if
                        call blockk (w(lvol),w(lvolk0),idimn-1,jdimn,
     .                               kdimn,idimc,jdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbck),iedge,1)
                     end if
c
c                 j = constant interface
c
                  else if (jside.eq.4) then
c
c                    check geometric mismatch
c
                     if (isklton.eq.1) then
                        call cblkj(nbli,idimn,jdimn,kdimn,idimc,jdimc,
     .                             kdimc,limblk(1,1,n),isva(1,1,n),
     .                             it,ir,iedge,w(lwxr),w(lwyr),w(lwzr),
     .                             w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                             geom_miss,mxbli)
                     end if
c
c                    interpolate q
c
                     ldim = 5
                     call blockj (w(lws),w(lwt),idimn,jdimn,kdimn,
     .                            idimc,kdimc,limblk(1,1,n),
     .                            isva(1,1,n),it,ir,ldim,
     .                            w(lbcj),iedge,0)
c
c                    interpolate turbulent quantities
c
                     if (ivmx.ge.2) then
                        ldim = 1
                        lwst = lw(13,in_blk)    
                        lwtt = lw(28,ic_blk)     
                        if (limblk(it,2,n).ne.1) then
                           lwtt   = lwtt + kdimc*(idimc-1)*1*2
                        end if
                        call blockj (w(lwst),w(lwtt),idimn,jdimn,
     .                              kdimn,idimc,kdimc,limblk(1,1,n),
     .                              isva(1,1,n),it,ir,ldim,
     .                              w(lbcj),iedge,0)
                     end if
                     if (ivmx.ge.4) then
                        ldim = nummem
                        lwst = lw(19,in_blk)    
                        lwtt = lw(23,ic_blk)     
                        if (limblk(it,2,n).ne.1) then
                           lwtt   = lwtt + kdimc*(idimc-1)*2*nummem
                        end if
                        call blockj (w(lwst),w(lwtt),idimn,jdimn,
     .                               kdimn,idimc,kdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbcj),iedge,0)
                     end if
c
c                    interpolate cell volumes
c
                     if (ivmx.ge.1 .and. ivolint.gt.0) then
                        ldim   = 1
                        lvol   = lw(8,in_blk)
                        lvolj0 = lw(49,ic_blk)
                        if (limblk(it,2,n).ne.1) then
                           lvolj0 = lvolj0 + kdimc*(idimc-1)*1*2
                        end if
                        call blockj (w(lvol),w(lvolj0),idimn-1,jdimn,
     .                               kdimn,idimc,kdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbcj),iedge,1)
                     end if
c
c                 i = constant interface
c
                  else if (jside.eq.5) then
c
c                    check geometric mismatch
c
                     if (isklton.eq.1) then
                        call cblki(nbli,idimn,jdimn,kdimn,idimc,jdimc,
     .                             kdimc,limblk(1,1,n),isva(1,1,n),
     .                             it,ir,iedge,w(lwxr),w(lwyr),w(lwzr),
     .                             w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                             geom_miss,mxbli)
                     end if
c
c                    interpolate q
c
                     ldim = 5
                     call blocki (w(lws),w(lwt),idimn,jdimn,kdimn,
     .                            jdimc,kdimc,limblk(1,1,n),
     .                            isva(1,1,n),it,ir,ldim,
     .                            w(lbci),iedge,0)
c
c                    interpolate turbulent quantities
c
                     if (ivmx.ge.2) then
                        ldim = 1
                        lwst = lw(13,in_blk)    
                        lwtt = lw(30,ic_blk)     
                        if (limblk(it,1,n).ne.1) then
                           lwtt   = lwtt + jdimc*kdimc*1*2
                        end if
                        call blocki (w(lwst),w(lwtt),idimn,jdimn,
     .                               kdimn,jdimc,kdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbci),iedge,0)
                     end if
                     if (ivmx.ge.4) then
                        ldim = nummem
                        lwst = lw(19,in_blk)    
                        lwtt = lw(25,ic_blk)     
                        if (limblk(it,1,n).ne.1) then
                           lwtt   = lwtt + jdimc*kdimc*2*nummem
                        end if
                        call blocki (w(lwst),w(lwtt),idimn,jdimn,
     .                               kdimn,jdimc,kdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbci),iedge,0)
                     end if
c
c                    interpolate cell volumes
c
                     if (ivmx.ge.1 .and. ivolint.gt.0) then
                        ldim   = 1
                        lvol   = lw(8,in_blk)
                        lvoli0 = lw(51,ic_blk)
                        if (limblk(it,1,n).ne.1) then
                           lvoli0 = lvoli0 + jdimc*kdimc*1*2
                        end if
                        call blocki (w(lvol),w(lvoli0),idimn-1,jdimn,
     .                               kdimn,jdimc,kdimc,limblk(1,1,n),
     .                               isva(1,1,n),it,ir,ldim,
     .                               w(lbci),iedge,1)
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
c***********************************************************************
c        Second Case: data needed to set 1-1 bc lies on another 
c                     processor
c***********************************************************************
#if defined DIST_MPI
#        ifdef BUILD_MPE
c
c        begin monitoring message passing
c
         call MPE_Log_event (20, 0, "Start BC_BLKINT")
#        endif
c
c        set baseline tag values
c
         ioffset = mxbli
         itag_x   = 1
         itag_q   = itag_x + ioffset
         itag_v   = itag_q + ioffset
         itag_t   = itag_v + ioffset
         itag_vol = itag_t + ioffset
c
c        post a bunch of receives first (for non-buffering implementations)
c        set the request index and index for wk
c
         kqintl = 1
         ireq   = 0
c
         do lcnt = is_blk(level),ie_blk(level)
c           ic_blk is current (to) block
c           in_blk is neighbor (from) block
            ic_blk  = isav_blk(lcnt,4)
            if (iadvance(ic_blk).ge.0) then
            in_blk  = isav_blk(lcnt,5)
            nd_dest = mblk2nd(ic_blk)
            nd_srce = mblk2nd(in_blk)
            if (nd_srce.ne.myid) then
               if (nd_dest.eq.myid) then
                  n = isav_blk(lcnt,1)
                  jface = isav_blk(lcnt,6)
                  idimn = idimg(in_blk)
                  jdimn = jdimg(in_blk)
                  kdimn = kdimg(in_blk)
                  if (jface.eq.1) maxdims = jdimn*kdimn
                  if (jface.eq.2) maxdims = kdimn*idimn
                  if (jface.eq.3) maxdims = jdimn*idimn
                  if (isklton.eq.1) then
                     kcheck = kqintl + maxdims*3
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_x + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*3,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,1)   = kqintl
                     keep_trac(n,2)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
c
                  ldim = 5
                  kcheck = kqintl + maxdims*ldim*2
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in bc_blkint',
     .               '...work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  mytag = itag_q + n
                  ireq  = ireq + 1
                  call MPI_IRecv (wk(kqintl), maxdims*ldim*2,
     .                           MY_MPI_REAL,
     .                           nd_srce,mytag,mycomm,
     .                           ireq_ar(ireq),ierr)
                  keep_trac(n,3)   = kqintl
                  keep_trac(n,4)   = ireq
                  keep_trac2(ireq) = lcnt
                  kqintl = kcheck
                  if (ivmx.ge.2) then
                     ldim = 1
                     kcheck = kqintl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_v + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*2,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,5)   = kqintl
                     keep_trac(n,6)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
                  if (ivmx.ge.4) then
                     ldim = nummem
                     kcheck = kqintl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_t + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*2,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,7)   = kqintl
                     keep_trac(n,8)   = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
                  if (ivmx.ge.1 .and. ivolint.gt.0) then
                     ldim = 1
                     kcheck = kqintl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_vol + n
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl), maxdims*ldim*2,
     .                              MY_MPI_REAL,
     .                              nd_srce,mytag,mycomm,
     .                              ireq_ar(ireq),ierr)
                     keep_trac(n,9)   = kqintl
                     keep_trac(n,10)  = ireq
                     keep_trac2(ireq) = lcnt
                     kqintl = kcheck
                  end if
               end if
            end if 
            end if
         end do
c
         if (myid.ne.myhost) then
            if (ireq.gt.mxbli*5) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),999) ireq,mxbli*5
 999           format(' problem in bc_blkint...ireq = ',i4,
     .         ' but max allowable value = mxbli*5 = ',i4)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
        end if
c
c        loop over all blocks looking for blocks that need to send out
c        info to other processors
c
         ktl   = kqintl
         ireq2 = 0
c
         do lcnt = is_blk(level),ie_blk(level)
c           ic_blk is current (to) block
c           in_blk is neighbor (from) block
            ic_blk  = isav_blk(lcnt,4)
            if (iadvance(ic_blk).ge.0) then
            in_blk  = isav_blk(lcnt,5)
            nd_dest = mblk2nd(ic_blk)
            nd_srce = mblk2nd(in_blk)
            if (nd_dest.ne.myid) then
               if (nd_srce.eq.myid) then
                  n     = isav_blk(lcnt,1)
                  it    = isav_blk(lcnt,2)
                  ir    = isav_blk(lcnt,3)
                  jface = isav_blk(lcnt,6)
                  jside = isav_blk(lcnt,8)
                  iedge = isav_blk(lcnt,10)
                  idimn = idimg(in_blk)
                  jdimn = jdimg(in_blk)
                  kdimn = kdimg(in_blk)
                  if (jface.eq.1) maxdims = jdimn*kdimn
                  if (jface.eq.2) maxdims = kdimn*idimn
                  if (jface.eq.3) maxdims = jdimn*idimn
c
                  lws  = lw( 1,in_blk)
                  lwxr = lw(10,in_blk)
                  lwyr = lw(11,in_blk)
                  lwzr = lw(12,in_blk)
                  if (isklton.eq.1) then
                     kcheck = ktl + maxdims*3
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (jside.eq.3) then
                        call pre_cblkk (idimn, jdimn, kdimn,
     .                                 limblk(1,1,n), isva(1,1,n), it, 
     .                                 ir, w(lwxr), w(lwyr), w(lwzr), 
     .                                 wk(ktl), maxdims, nval1,myid,
     .                                 ibufdim,nbuf,bou,nou)
                     else if (jside.eq.4) then
                        call pre_cblkj (idimn, jdimn, kdimn,
     .                                 limblk(1,1,n), isva(1,1,n), it,
     .                                 ir, w(lwxr), w(lwyr), w(lwzr),
     .                                 wk(ktl), maxdims, nval1,myid,
     .                                 ibufdim,nbuf,bou,nou)
                     else
                        call pre_cblki (idimn, jdimn, kdimn,
     .                                 limblk(1,1,n), isva(1,1,n), it,
     .                                 ir, w(lwxr), w(lwyr), w(lwzr),
     .                                 wk(ktl), maxdims, nval1,myid,
     .                                 ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_x + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*3,
     .                             MY_MPI_REAL,
     .                             nd_dest, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
                  ldim = 5
                  kcheck = ktl + maxdims*ldim*2
                  if (kcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in bc_blkint',
     .               '...work array insufficient',kcheck
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  if (jside.eq.3) then
                     call pre_blockk (idimn, jdimn, kdimn,
     .                               limblk(1,1,n), isva(1,1,n), it,
     .                               ir, w(lws), wk(ktl), maxdims, ldim,
     .                               nval2,0,myid,mblk2nd,maxbl,
     .                               bou,ibufdim,nbuf,nou)
                  else if (jside.eq.4) then
                     call pre_blockj (idimn, jdimn, kdimn,
     .                               limblk(1,1,n), isva(1,1,n), it,
     .                               ir, w(lws), wk(ktl), maxdims, ldim,
     .                               nval2,0,myid,mblk2nd,maxbl,
     .                               bou,ibufdim,nbuf,nou)
                  else
                     call pre_blocki (idimn, jdimn, kdimn,
     .                               limblk(1,1,n), isva(1,1,n), it,
     .                               ir, w(lws), wk(ktl), maxdims, ldim,
     .                               nval2,0,myid,mblk2nd,maxbl,
     .                               bou,ibufdim,nbuf,nou)
                  end if
                  mytag = itag_q + n
                  ireq2 = ireq2 + 1
                  call MPI_ISend(wk(ktl), maxdims*ldim*2,
     .                          MY_MPI_REAL,
     .                          nd_dest, mytag, mycomm,
     .                          ireq_snd(ireq2), ierr)
                  ktl = kcheck
                  if (ivmx.ge.2) then
                     lwst = lw(13,in_blk)
                     ldim = 1
                     kcheck = ktl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (jside.eq.3) then
                        call pre_blockk (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n), 
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else if (jside.eq.4) then
                        call pre_blockj (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n), 
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else
                        call pre_blocki (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n), 
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     end if
                     mytag = itag_v + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*2,
     .                             MY_MPI_REAL,
     .                             nd_dest, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
                  if (ivmx.ge.4) then
                     lwst = lw(19,in_blk)
                     ldim = nummem
                     kcheck = ktl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck             
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (jside.eq.3) then
                        call pre_blockk (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n), 
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else if (jside.eq.4) then
                        call pre_blockj (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n),
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else
                        call pre_blocki (idimn, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n), 
     .                                  it, ir, w(lwst), wk(ktl), 
     .                                  maxdims, ldim, nval2,0,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     end if
                     mytag = itag_t + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*2,
     .                             MY_MPI_REAL,
     .                             nd_dest, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
                  if (ivmx.ge.1 .and. ivolint.gt.0) then
                     lwst = lw(8,in_blk)
                     ldim = 1
                     kcheck = ktl + maxdims*ldim*2
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in bc_blkint',
     .                  '...work array insufficient',kcheck
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     if (jside.eq.3) then
                        call pre_blockk (idimn-1, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n),
     .                                  it, ir, w(lwst), wk(ktl),
     .                                  maxdims, ldim, nval2,1,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else if (jside.eq.4) then
                        call pre_blockj (idimn-1, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n),
     .                                  it, ir, w(lwst), wk(ktl),
     .                                  maxdims, ldim, nval2,1,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     else
                        call pre_blocki (idimn-1, jdimn, kdimn,
     .                                  limblk(1,1,n), isva(1,1,n),
     .                                  it, ir, w(lwst), wk(ktl),
     .                                  maxdims, ldim, nval2,1,
     .                                  myid,mblk2nd,maxbl,
     .                                  bou,ibufdim,nbuf,nou)
                     end if
                     mytag = itag_vol + n
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), maxdims*ldim*2,
     .                             MY_MPI_REAL,
     .                             nd_dest, mytag, mycomm,
     .                             ireq_snd(ireq2), ierr)
                     ktl = kcheck
                  end if
               end if
            end if
            end if
         end do
c
c        set 1-1 interface bc's 
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
               lcnt   = keep_trac2(index_ar(nnn))
               n      = isav_blk(lcnt,1)
               it     = isav_blk(lcnt,2)
               ir     = isav_blk(lcnt,3)
               ic_blk = isav_blk(lcnt,4)
               in_blk = isav_blk(lcnt,5)
               jface  = isav_blk(lcnt,6)
               jedge  = isav_blk(lcnt,7)
               jside  = isav_blk(lcnt,8)
               lwt    = isav_blk(lcnt,9)
               iedge  = isav_blk(lcnt,10)
               iss    = isav_blk(lcnt,11)
               ise    = isav_blk(lcnt,12)
               jss    = isav_blk(lcnt,13)
               jse    = isav_blk(lcnt,14)
               kss    = isav_blk(lcnt,15)
               kse    = isav_blk(lcnt,16)
               iti    = isav_blk(lcnt,17)
               nd_dest = mblk2nd(ic_blk)
               nd_srce = mblk2nd(in_blk)
c
               if (iadvance(ic_blk).ge.0) then
c
               call lead(ic_blk,lw,lw2,maxbl)
c
c              set dimensions of blocks involved
c
               idimn = idimg(in_blk)
               jdimn = jdimg(in_blk)
               kdimn = kdimg(in_blk)
               idimc = idimg(ic_blk)
               jdimc = jdimg(ic_blk)
               kdimc = kdimg(ic_blk)
               if (jface.eq.1) maxdims = jdimn*kdimn
               if (jface.eq.2) maxdims = kdimn*idimn
               if (jface.eq.3) maxdims = jdimn*idimn
c
c              set pointers for neighboring block
c
               lws  = lw( 1,in_blk)
               lwxr = lw(10,in_blk)
               lwyr = lw(11,in_blk)
               lwzr = lw(12,in_blk)
c
c              set pointers for current block (note: lwt has
c              already been set above, and corresponds to 
c              qi0, qj0, or qk0, depending on the interface)
c
               lwxt = lw(10,ic_blk)
               lwyt = lw(11,ic_blk)
               lwzt = lw(12,ic_blk)
c
c              k = constant interface 
c
               if (jside.eq.3) then
c
c                 check geometric mismatch
c
                  if (isklton.eq.1) then
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call cblkk_d(nbli,idimn,jdimn,kdimn,idimc,jdimc,
     .                               kdimc,limblk(1,1,n),isva(1,1,n),
     .                               it,ir,iedge,wk(kqintl),maxdims,
     .                               w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                               geom_miss,mxbli)
                     end if
                  end if
c
c                 interpolate q
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     call blockk_d (wk(kqintl),w(lwt),idimn,jdimn,kdimn,
     .                             idimc,jdimc,limblk(1,1,n),
     .                             isva(1,1,n),it,ir,maxdims,ldim,
     .                             w(lbck),iedge)
                  end if
c
c                 interpolate turbulent quantities
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        lwtt = lw(29,ic_blk)
                        if (limblk(it,3,n).ne.1) then
                           lwtt   = lwtt + jdimc*(idimc-1)*1*2
                        end if
                        kqintl = keep_trac(n,5)
                        call blockk_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                kdimn,idimc,jdimc,limblk(1,1,n),
     .                                isva(1,1,n),it,ir,maxdims,ldim,
     .                                w(lbck),iedge)
                     end if
                  end if
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        lwtt = lw(24,ic_blk)
                        if (limblk(it,3,n).ne.1) then
                           lwtt   = lwtt + jdimc*(idimc-1)*2*nummem
                        end if
                        kqintl = keep_trac(n,7)
                        call blockk_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                kdimn,idimc,jdimc,limblk(1,1,n),
     .                                isva(1,1,n),it,ir,maxdims,ldim,
     .                                w(lbck),iedge)
                     end if
                  end if
c
c                 interpolate cell volumes
c
                  if (ivmx.ge.1 .and. ivolint.gt.0) then
                     if (index_ar(nnn) .eq. keep_trac(n,10)) then
                        ldim   = 1
                        lvolk0 = lw(50,ic_blk)
                        if (limblk(it,3,n).ne.1) then
                           lvolk0 = lvolk0 + jdimc*(idimc-1)*1*2
                        end if
                        kqintl = keep_trac(n,9)
                        call blockk_d (wk(kqintl),w(lvolk0),idimn,jdimn,
     .                                 kdimn,idimc,jdimc,limblk(1,1,n),
     .                                 isva(1,1,n),it,ir,maxdims,ldim,
     .                                 w(lbck),iedge)
                     end if
                  end if
c
c              j = constant interface
c
               else if (jside.eq.4) then
c
c                 check geometric mismatch
c
                  if (isklton.eq.1) then
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call cblkj_d(nbli,idimc,jdimc,kdimc,
     .                               limblk(1,1,n),isva(1,1,n),
     .                               it,ir,iedge,wk(kqintl),maxdims,
     .                               w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                               geom_miss,mxbli)
                     end if
                  end if
c
c                 interpolate q
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     call blockj_d (wk(kqintl),w(lwt),idimn,jdimn,kdimn,
     .                             idimc,kdimc,limblk(1,1,n),
     .                             isva(1,1,n),it,ir,maxdims,ldim,
     .                             w(lbcj),iedge)
                  end if
c
c                 interpolate turbulent quantities
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        lwtt = lw(28,ic_blk)
                        if (limblk(it,2,n).ne.1) then
                           lwtt   = lwtt + kdimc*(idimc-1)*1*2
                        end if
                        kqintl = keep_trac(n,5)
                        call blockj_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                kdimn,idimc,kdimc,limblk(1,1,n),
     .                                isva(1,1,n),it,ir,maxdims,ldim,
     .                                w(lbcj),iedge)
                     end if
                  end if
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        lwtt = lw(23,ic_blk)
                        if (limblk(it,2,n).ne.1) then
                           lwtt   = lwtt + kdimc*(idimc-1)*2*nummem
                        end if
                        kqintl = keep_trac(n,7)
                        call blockj_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                 kdimn,idimc,kdimc,limblk(1,1,n),
     .                                 isva(1,1,n),it,ir,maxdims,ldim,
     .                                 w(lbcj),iedge)
                     end if
                  end if
c
c                 interpolate cell volumes
c
                  if (ivmx.ge.1 .and. ivolint.gt.0) then
                     if (index_ar(nnn) .eq. keep_trac(n,10)) then
                        ldim   = 1
                        lvolj0 = lw(49,ic_blk)
                        if (limblk(it,2,n).ne.1) then
                           lvolj0 = lvolj0 + kdimc*(idimc-1)*1*2
                        end if
                        kqintl = keep_trac(n,9)
                        call blockj_d (wk(kqintl),w(lvolj0),idimn,jdimn,
     .                                 kdimn,idimc,kdimc,limblk(1,1,n),
     .                                 isva(1,1,n),it,ir,maxdims,ldim,
     .                                 w(lbcj),iedge)
                     end if
                  end if
c
c              i = constant interface
c
               else if (jside.eq.5) then
c
c                 check geometric mismatch
c
                  if (isklton.eq.1) then
                     if (index_ar(nnn) .eq. keep_trac(n,2)) then
                        kqintl = keep_trac(n,1)
                        call cblki_d(nbli,idimn,jdimn,kdimn,idimc,jdimc,     
     .                               kdimc,limblk(1,1,n),isva(1,1,n),        
     .                               it,ir,iedge,wk(kqintl),maxdims,
     .                               w(lwxt),w(lwyt),w(lwzt),iti,lcnt,
     .                               geom_miss,mxbli)
                     end if
                  end if
c
c                 interpolate q
c
                  if (index_ar(nnn) .eq. keep_trac(n,4)) then
                     ldim = 5
                     kqintl = keep_trac(n,3)
                     call blocki_d (wk(kqintl),w(lwt),idimn,jdimn,kdimn,
     .                              jdimc,kdimc,limblk(1,1,n),
     .                              isva(1,1,n),it,ir,maxdims,ldim,
     .                              w(lbci),iedge)
                  end if
c
c                 interpolate turbulent quantities
c
                  if (ivmx.ge.2) then
                     if (index_ar(nnn) .eq. keep_trac(n,6)) then
                        ldim = 1
                        lwtt = lw(30,ic_blk)
                        if (limblk(it,1,n).ne.1) then
                           lwtt   = lwtt + jdimc*kdimc*1*2
                        end if
                        kqintl = keep_trac(n,5)
                        call blocki_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                kdimn,jdimc,kdimc,limblk(1,1,n),
     .                                isva(1,1,n),it,ir,maxdims,ldim,
     .                                w(lbci),iedge)
                     end if
                  end if
                  if (ivmx.ge.4) then
                     if (index_ar(nnn) .eq. keep_trac(n,8)) then
                        ldim = nummem
                        lwtt = lw(25,ic_blk)
                        if (limblk(it,1,n).ne.1) then
                           lwtt   = lwtt + jdimc*kdimc*2*nummem
                        end if
                        kqintl = keep_trac(n,7)
                        call blocki_d (wk(kqintl),w(lwtt),idimn,jdimn,
     .                                kdimn,jdimc,kdimc,limblk(1,1,n),
     .                                isva(1,1,n),it,ir,maxdims,ldim,
     .                                w(lbci),iedge)
                     end if
                  end if
c
c                 interpolate cell volumes
c
                  if (ivmx.ge.1 .and. ivolint.gt.0) then
                     if (index_ar(nnn) .eq. keep_trac(n,10)) then
                        ldim   = 1
                        lvoli0 = lw(51,ic_blk)
                        if (limblk(it,1,n).ne.1) then
                           lvoli0 = lvoli0 + jdimc*kdimc*1*2
                        end if
                        kqintl = keep_trac(n,9)
                        call blocki_d (wk(kqintl),w(lvoli0),idimn,jdimn,
     .                                 kdimn,jdimc,kdimc,limblk(1,1,n),
     .                                 isva(1,1,n),it,ir,maxdims,ldim,
     .                                 w(lbci),iedge)
                     end if
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
         call MPE_Log_event (29, 0, "End BC_BLKINT")
#        endif
#endif
c
      end if
c
      return
      end
