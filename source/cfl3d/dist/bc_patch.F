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
      subroutine bc_patch(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,maxgr,
     .                    intmax,nsub1,maxxe,iadvance,jdimg,kdimg,
     .                    idimg,ninter,windex,iindex,nblkpt,dthetxx,
     .                    dthetyy,dthetzz,isav_pat,isav_pat_b,
     .                    ireq_ar,ireq_snd,index_ar,
     .                    keep_trac,keep_trac2,myid,myhost,mycomm,
     .                    mblk2nd,nou,bou,nbuf,ibufdim,istat2,
     .                    istat_size,nummem)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Update patched-grid interface boundary conditions.
c      Interface conditions are set by linearly interpolating q
c      from one grid to ghost cells of another grid using generalized
c      coordinates. Interpolation may be done with or without limiting:
c
c     lim_ptch = 1...limiter employed for interpolation
c                0...no limiter employed
c
c     unlimited is standard, but the use of limiting can sometimes
c     cure problems at interfaces
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
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension istat2(istat_size,intmax*nsub1*3)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension iadvance(maxbl),mblk2nd(maxbl)
      dimension ireq_ar(intmax*nsub1*3),index_ar(intmax*nsub1*3),
     .          ireq_snd(intmax*nsub1*3),keep_trac(intmax,nsub1,6),
     .          keep_trac2(intmax*nsub1*3,2)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),nblkpt(maxxe),
     .          dthetxx(intmax,nsub1),dthetyy(intmax,nsub1),
     .          dthetzz(intmax,nsub1)
      dimension isav_pat(intmax,17),isav_pat_b(intmax,nsub1,6)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /maxiv/ ivmx
      common /is_patch/ is_pat(5),ie_pat(5),ipatch1st
c
c     lim_ptch = 0 for patching without limiting; = 1 with limiting
c
      lim_ptch = 0
c
      if (ntime.gt.0 .and. abs(ninter).gt.0) then
c
c***********************************************************************
c        First Case: some or all data needed to set patch bc lies on
c                    the current processor
c***********************************************************************
c
c        first order interpolation for q if icyc .le. nitfo
c        or if ipatch1st is specified as 1
c
         ifo = 0
         if (icyc.gt.nitfo .and. ipatch1st.ne.1) ifo = 1
c
         do lcnt = is_pat(level),ie_pat(level)
c           nb_to is current (to) block
c           nb_from is neighbor (from) block
            nb_to   = isav_pat(lcnt,1)
            nfb     = isav_pat(lcnt,2)
            nd_recv = mblk2nd(nb_to)
            if (iadvance(nb_to).ge.0) then
               do l = 1,nfb
                  nb_from = isav_pat_b(lcnt,l,1)
                  nd_srce = mblk2nd(nb_from)
c
                  if (nd_srce.eq.myid .and. nd_recv.eq.myid) then
c
                  call lead(nb_to,lw,lw2,maxbl)
c
                  j21    = isav_pat(lcnt,3)
                  j22    = isav_pat(lcnt,4)
                  k21    = isav_pat(lcnt,5)
                  k22    = isav_pat(lcnt,6)
                  lqedge = isav_pat(lcnt,7)
                  lqedgt = isav_pat(lcnt,8) 
                  lqedgv = isav_pat(lcnt,9)
                  lqedgb = isav_pat(lcnt,10)
                  jmax2  = isav_pat(lcnt,11)
                  kmax2  = isav_pat(lcnt,12)
                  mint1  = isav_pat(lcnt,13)
                  mint2  = isav_pat(lcnt,14)
                  lst    = isav_pat(lcnt,15)
                  npt    = isav_pat(lcnt,16)
                  icheck = isav_pat(lcnt,17)
                  mtype  = isav_pat_b(lcnt,l,2)
	          jmax1  = isav_pat_b(lcnt,l,3)
	          kmax1  = isav_pat_b(lcnt,l,4)
                  idim1  = idimg(nb_from)
                  jdim1  = jdimg(nb_from)
                  kdim1  = kdimg(nb_from)
                  iint1  = isav_pat_b(lcnt,l,5)
                  iint2  = isav_pat_b(lcnt,l,6)
                  lqintl = lw(1,nb_from)
                  lqintt = lw(19,nb_from)
                  lqintv = lw(13,nb_from)
                  dthtx  = dthetxx(icheck,l)
                  dthty  = dthetyy(icheck,l)
                  dthtz  = dthetzz(icheck,l)
                  lwk1   = 1
                  lwk2   = lwk1 + jmax1*kmax1
                  lwk3   = lwk2 + jmax1*kmax1
                  lwk4   = lwk3 + jmax1*kmax1
                  klimit = lwk4 + jdim1*kdim1*idim1*5 - 1
                  if (klimit.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in bc_patch...',
     .                             ' insufficient wk storage'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
c
c                 interpolate q data
c
                  ldim = 5
                  iq = 1
                  jkldim = ldim*jmax1*kmax1
c
c                 interpolate to first layer of ghost cells
c
                  call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .            nsub1,jmax1,kmax1,l,iint1,w(lqintl),w(lqedgb),
     .            mint1,w(lqedge),windex(lst,1),windex(lst,2),
     .            wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .            icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,k22,
     .            wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
c                 interpolate to second layer of ghost cells
c
                  call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .            nsub1,jmax1,kmax1,l,iint2,w(lqintl),w(lqedgb),
     .            mint2,w(lqedge),windex(lst,1),windex(lst,2),
     .            wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .            icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,k22,
     .            wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
c                 interpolate vist3d data
c
                  if (ivmx.ge.2) then
                     ldim = 1
                     iq = 0
                     jkldim = ldim*jmax1*kmax1
c
c                    interpolate to first layer of ghost cells
c
                     call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint1,w(lqintv),w(lqedgb),
     .               mint1,w(lqedgv),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .               k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
c                    interpolate to second layer of ghost cells
c
                     call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint2,w(lqintv),w(lqedgb),
     .               mint2,w(lqedgv),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .               k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
                  end if
c
c                 interpolate turb. data
c
                  if (ivmx.ge.4) then
                     ldim = nummem
                     iq = 0
                     jkldim = ldim*jmax1*kmax1
c
c                    interpolate to first layer of ghost cells
c
                     call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint1,w(lqintt),w(lqedgb),
     .               mint1,w(lqedgt),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .               k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
c                    interpolate to second layer of ghost cells
c
                     call int2(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint2,w(lqintt),w(lqedgb),
     .               mint2,w(lqedgt),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .               k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
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
c***********************************************************************
c        Second Case: some or all data needed to set patch bc lies on
c                     another processor
c***********************************************************************
#if defined(DIST_MPI)
#        ifdef BUILD_MPE
c
c        begin monitoring message passing
c
         call MPE_Log_event (10, 0, "Start BC_PATCH")
#        endif
c
c        first order interpolation for q if icyc .le. nitfo
c        or if ipatch1st is specified as 1
c
         ifo = 0
         if (icyc.gt.nitfo .and. ipatch1st.ne.1) ifo = 1
c
c        set baseline tag values
c
         ioffset = intmax*nsub1
         itag_q  = 1
         itag_v  = itag_q + ioffset
         itag_t  = itag_v + ioffset
c
c        post a bunch of receives first (for non-buffering implementations)
c        set the request index and index for wk
c
         kqintl = 1
         ireq   = 0
         itag   = 0
c
         do lcnt = is_pat(level),ie_pat(level)
c           nb_to is current (to) block
c           nb_from is neighbor (from) block
            nfb = isav_pat(lcnt,2)
            n   = isav_pat(lcnt,17)
            do ll = 1, nfb
               itag    = itag + 1
               jmax1   = isav_pat_b(lcnt,ll,3)
               kmax1   = isav_pat_b(lcnt,ll,4)
               nb_from = isav_pat_b(lcnt,ll,1)
               nb_to   = isav_pat(lcnt,1)
               if (iadvance(nb_to).ge.0) then
               nd_srce = mblk2nd(nb_from)
               nd_recv = mblk2nd(nb_to)
               if (nd_srce.ne.myid) then
                  if (nd_recv.eq.myid) then
c
c                    receive q data
c
                     ldim = 5
                     jkldim = ldim*jmax1*kmax1
                     kcheck = kqintl + 2*jkldim
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in',
     .                  ' bc_patch....work array insufficient'
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     mytag = itag_q + itag
                     ireq  = ireq + 1
                     call MPI_IRecv (wk(kqintl),2*jkldim,
     .                             MY_MPI_REAL,
     .                             nd_srce,mytag,mycomm,
     .                             ireq_ar(ireq),ierr)
                     keep_trac(n,ll,1)  = kqintl
                     keep_trac(n,ll,2)  = ireq
                     keep_trac2(ireq,1) = lcnt
                     keep_trac2(ireq,2) = ll
                     kqintl = kcheck
c
c                    receive vist3d data
c
                     if (ivmx.ge.2) then
                        ldim = 1
                        jkldim = ldim*jmax1*kmax1
                        kcheck = kqintl + 2*jkldim
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)'stopping in',
     .                     ' bc_patch....work array insufficient'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        mytag = itag_v + itag
                        ireq  = ireq + 1
                        call MPI_IRecv (wk(kqintl),2*jkldim,
     .                                 MY_MPI_REAL,
     .                                 nd_srce,mytag,mycomm,
     .                                 ireq_ar(ireq),ierr)
                        keep_trac(n,ll,3)  = kqintl
                        keep_trac(n,ll,4)  = ireq
                        keep_trac2(ireq,1) = lcnt
                        keep_trac2(ireq,2) = ll
                        kqintl = kcheck
                     end if
c
c                    receive turb. data
c
                     if (ivmx.ge.4) then
                        ldim = nummem
                        jkldim = ldim*jmax1*kmax1
                        kcheck = kqintl + 2*jkldim
                        if (kcheck.gt.nwork) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),*)'stopping in',
     .                     ' bc_patch....work array insufficient'
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        mytag = itag_t + itag
                        ireq  = ireq + 1
                        call MPI_IRecv (wk(kqintl),2*jkldim,
     .                                 MY_MPI_REAL,
     .                                 nd_srce,mytag,mycomm,
     .                                 ireq_ar(ireq),ierr)
                        keep_trac(n,ll,5)  = kqintl
                        keep_trac(n,ll,6)  = ireq
                        keep_trac2(ireq,1) = lcnt
                        keep_trac2(ireq,2) = ll
                        kqintl = kcheck
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
         end do
c
c        check tags
c
         if (myid.ne.myhost) then
            if (itag_q + itag .gt. itag_v) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'not enough tags between ',
     .                       'itag_q and itag_v'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'itag,itag_q,itag_v = ',
     .                        itag,itag_q,itag_v
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            if (itag_v + itag .gt. itag_t) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'not enough tags between ',
     .                       'itag_v and itag_t'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'itag,itag_v,itag_t = ',
     .                       itag,itag_v,itag_t
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
c
            if (ireq.gt.intmax*nsub1*3) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),999) ireq,intmax*nsub1*3
 999           format(' problem in bc_patch...ireq = ',i4,
     .         ' but max allowable value = intmax*nsub1*3 = ',i4)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
c        loop over all patch boundaries looking for blocks that 
c        need to send out info to other processors
c
         ktl   = kqintl
         ireq2 = 0
         itag  = 0
c
         do lcnt = is_pat(level),ie_pat(level)
c           nb_to is current (to) block
c           nb_from is neighbor (from) block
            nfb = isav_pat(lcnt,2)
            do ll = 1, nfb
               itag = itag + 1
               nb_from = isav_pat_b(lcnt,ll,1)
               nb_to   = isav_pat(lcnt,1)
               if (iadvance(nb_to).ge.0) then
               nd_srce = mblk2nd(nb_from)
               nd_recv = mblk2nd(nb_to)
               if (nd_srce.eq.myid) then
                  if (nd_recv.ne.myid) then
c
                     jdim1 = jdimg(nb_from)
                     kdim1 = kdimg(nb_from)
                     idim1 = idimg(nb_from)
                     mtype = isav_pat_b(lcnt,ll,2)
                     jmax1 = isav_pat_b(lcnt,ll,3)
                     kmax1 = isav_pat_b(lcnt,ll,4)
                     iint1 = isav_pat_b(lcnt,ll,5)
                     iint2 = isav_pat_b(lcnt,ll,6)
c
c                    load 2 planes of q data from full "from" block
c                    to a work array and send to the appropriate processor
c
                     ldim = 5
                     lqintl = lw(1,nb_from)
                     call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                              kmax1, ldim, w(lqintl), wk(ktl), 
     .                              mtype, iint1)
                     nvals = jmax1*kmax1*ldim
                     ktla = ktl + nvals
                     kcheck = ktla + nvals
                     if (kcheck.gt.nwork) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*)' stopping in',
     .                  ' bc_patch....work array insufficient'
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                     call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                              kmax1, ldim, w(lqintl), wk(ktla),
     .                              mtype, iint2)
                     mytag = itag_q + itag
                     ireq2 = ireq2 + 1
                     call MPI_ISend(wk(ktl), 2*nvals, 
     .                            MY_MPI_REAL,
     .                            nd_recv, mytag, mycomm,
     .                            ireq_snd(ireq2), ierr)
                     ktl = kcheck
c
c                    load 2 planes of vist3d data from full "from" block
c                    to a work array and send to the appropriate processor
c
                     if (ivmx.ge.2) then
                         lqintv = lw(13,nb_from)
                         ldim = 1
                         call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                                  kmax1, ldim, w(lqintv), wk(ktl),
     .                                  mtype, iint1)
                         nvals = jmax1*kmax1*ldim
                         ktla = ktl + nvals
                         kcheck = ktla + nvals
                         if (kcheck.gt.nwork) then
                            nou(1) = min(nou(1)+1,ibufdim)
                            write(bou(nou(1),1),*)' stopping in',
     .                      ' bc_patch....work array insufficient'
                            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                         end if
                         call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                                  kmax1, ldim, w(lqintv),wk(ktla),
     .                                  mtype, iint2)
                         mytag = itag_v + itag
                         ireq2 = ireq2 + 1
                         call MPI_ISend(wk(ktl), 2*nvals, 
     .                                MY_MPI_REAL,
     .                                nd_recv, mytag, mycomm,
     .                                ireq_snd(ireq2), ierr)
                         ktl = kcheck
                     end if
c
c                    load 2 planes of turb. data from full "from" block
c                    to a work array and send to the appropriate processor
c
                     if (ivmx.ge.4) then
                         lqintt = lw(19,nb_from)
                         ldim = nummem
                         call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                                  kmax1, ldim, w(lqintt), wk(ktl),
     .                                  mtype, iint1)        
                         nvals = jmax1*kmax1*ldim
                         ktla = ktl + nvals
                         kcheck = ktla + nvals
                         if (kcheck.gt.nwork) then
                            nou(1) = min(nou(1)+1,ibufdim)
                            write(bou(nou(1),1),*)' stopping in',
     .                      ' bc_patch....work array insufficient'
                            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                         end if
                         call get_bvals(jdim1, kdim1, idim1, jmax1, 
     .                                  kmax1, ldim, w(lqintt),wk(ktla),
     .                                  mtype, iint2)
                         mytag = itag_t + itag
                         ireq2 = ireq2 + 1
                         call MPI_ISend(wk(ktl), 2*nvals, 
     .                                MY_MPI_REAL,
     .                                nd_recv, mytag, mycomm,
     .                                ireq_snd(ireq2), ierr)
                         ktl = kcheck
                     end if
                  end if 
               end if
               end if
            end do
         end do
c
c        check tags
c
         if (myid.ne.myhost) then
            if (itag_q + itag .gt. itag_v) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'not enough tags between ',
     .                       'itag_q and itag_v'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'itag,itag_q,itag_v = ',
     .                        itag,itag_q,itag_v
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            if (itag_v + itag .gt. itag_t) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'not enough tags between ',
     .                       'itag_v and itag_t'
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'itag,itag_v,itag_t = ',
     .                       itag,itag_v,itag_t
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
c        set patch interface bc's
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
               lcnt   = keep_trac2(index_ar(nnn),1)
               nb_to  = isav_pat(lcnt,1)        
               nfb    = isav_pat(lcnt,2)
               j21    = isav_pat(lcnt,3)
               j22    = isav_pat(lcnt,4)
               k21    = isav_pat(lcnt,5)
               k22    = isav_pat(lcnt,6)
               lqedge = isav_pat(lcnt,7)
               lqedgt = isav_pat(lcnt,8) 
               lqedgv = isav_pat(lcnt,9)
               lqedgb = isav_pat(lcnt,10)
               jmax2  = isav_pat(lcnt,11)
               kmax2  = isav_pat(lcnt,12)
               mint1  = isav_pat(lcnt,13)
               mint2  = isav_pat(lcnt,14)
               lst    = isav_pat(lcnt,15)
               npt    = isav_pat(lcnt,16)
               icheck = isav_pat(lcnt,17)
               n      = icheck
c
               if (iadvance(nb_to).ge.0) then
c
                  call lead(nb_to,lw,lw2,maxbl)
c
                  l = keep_trac2(index_ar(nnn),2)
                  nb_from = isav_pat_b(lcnt,l,1)
                  nd_srce = mblk2nd(nb_from)
                  mtype   = isav_pat_b(lcnt,l,2)
	          jmax1   = isav_pat_b(lcnt,l,3)
	          kmax1   = isav_pat_b(lcnt,l,4)
                  idim1   = idimg(nb_from)
                  jdim1   = jdimg(nb_from)
                  kdim1   = kdimg(nb_from)
                  iint1   = isav_pat_b(lcnt,l,5)
                  iint2   = isav_pat_b(lcnt,l,6)
                  lqintl  = lw(1,nb_from)
                  lqintt  = lw(19,nb_from)
                  lqintv  = lw(13,nb_from)
                  dthtx   = dthetxx(icheck,l)
                  dthty   = dthetyy(icheck,l)
                  dthtz   = dthetzz(icheck,l)
                  lwk1    = ktl
                  lwk2    = lwk1 + jmax1*kmax1
                  lwk3    = lwk2 + jmax1*kmax1
                  lwk4    = lwk3 + jmax1*kmax1
                  klimit  = lwk4 + jmax1*kmax1*5 - 1 
                  if (klimit.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in bc_patch...',
     .                             ' insufficient wk storage'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
c
c                 interpolate q data
c
                  if (index_ar(nnn) .eq. keep_trac(n,l,2)) then
                     ldim = 5
                     iq = 1
                     jkldim = ldim*jmax1*kmax1
                     kqintl = keep_trac(n,l,1)
c
c                    interpolate to first layer of ghost cells
c
                     call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint1,wk(kqintl),w(lqedgb),
     .               mint1,w(lqedge),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,k22,
     .               wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
                     kqintl = kqintl + jkldim
c
c                    interpolate to second layer of ghost cells
c
                     call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .               nsub1,jmax1,kmax1,l,iint2,wk(kqintl),w(lqedgb),
     .               mint2,w(lqedge),windex(lst,1),windex(lst,2),
     .               wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .               icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,k22,
     .               wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
                  end if
c
c                 interpolate vist3d data
c
                  if (ivmx.ge.2) then
                      if (index_ar(nnn) .eq. keep_trac(n,l,4)) then
                         ldim = 1
                         iq = 0
                         jkldim = ldim*jmax1*kmax1
                         kqintl = keep_trac(n,l,3)
c
c                        interpolate to first layer of ghost cells
c
                         call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .                   nsub1,jmax1,kmax1,l,iint1,wk(kqintl),w(lqedgb),
     .                   mint1,w(lqedgv),windex(lst,1),windex(lst,2),
     .                   wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .                   icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .                   k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
                         kqintl = kqintl + jkldim
c
c                        interpolate to second layer of ghost cells
c
                         call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .                   nsub1,jmax1,kmax1,l,iint2,wk(kqintl),w(lqedgb),
     .                   mint2,w(lqedgv),windex(lst,1),windex(lst,2),
     .                   wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .                   icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .                   k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
                      end if
                  end if
c
c                 interpolate turb. data
c
                  if (ivmx.ge.4) then
                      if (index_ar(nnn) .eq. keep_trac(n,l,6)) then
                         ldim = nummem
                         iq = 0
                         jkldim = ldim*jmax1*kmax1
                         kqintl = keep_trac(n,l,5)
c
c                        interpolate to first layer of ghost cells
c
                         call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .                   nsub1,jmax1,kmax1,l,iint1,wk(kqintl),w(lqedgb),
     .                   mint1,w(lqedgt),windex(lst,1),windex(lst,2),
     .                   wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .                   icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .                   k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
c
                         kqintl = kqintl + jkldim
c
c                        interpolate to second layer of ghost cells
c
                         call int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,
     .                   nsub1,jmax1,kmax1,l,iint2,wk(kqintl),w(lqedgb),
     .                   mint2,w(lqedgt),windex(lst,1),windex(lst,2),
     .                   wk(lwk1),wk(lwk2),wk(lwk3),nblkpt(lst),intmax,
     .                   icheck,mtype,iindex,ifo,ldim,npt,j21,j22,k21,
     .                   k22,wk(lwk4),dthtx,dthty,dthtz,lim_ptch)
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
         call MPE_Log_event (19, 0, "End BC_PATCH")
#        endif
#endif
c
      end if
c
      return
      end
