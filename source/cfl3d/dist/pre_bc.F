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
      subroutine pre_bc(lw,lw2,lw_temp,maxbl,maxgr,maxseg,ninter,intmax,
     .                  nsub1,iindex,isav_pat,jjmax1,
     .                  kkmax1,iiint1,iiint2,mxbli,nbli,limblk,isva,
     .                  nblon,nblk,lbcprd,isav_prd,
     .                  bcvali,bcvalj,bcvalk,nblg,lbcemb,
     .                  iemg,igridg,isav_emb,
     .                  iviscg,jdimg,kdimg,idimg,nbci0,nbcj0,nbck0,
     .                  nbcjdim,nbckdim,nbcidim,ibcinfo,jbcinfo,
     .                  kbcinfo,iadvance,myid,myhost,mycomm,
     .                  mblk2nd,nou,bou,nbuf,ibufdim,isav_pat_b,levelg,
     .                  nblock,isav_blk)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set up the data arrays needed for asynchronous message
c     passing of boundary condition data. 
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
      dimension lw(65,maxbl),lw2(43,maxbl),lw_temp(65,maxbl)
      dimension jjmax1(nsub1),kkmax1(nsub1),iiint1(nsub1),iiint2(nsub1)
      dimension igridg(maxbl),nblg(maxgr),levelg(maxbl),iviscg(maxbl,3)
      dimension bcvali(maxbl,maxseg,12,2),bcvalj(maxbl,maxseg,12,2),
     .          bcvalk(maxbl,maxseg,12,2),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension iadvance(maxbl),mblk2nd(maxbl),iemg(maxgr),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),isva(2,2,mxbli),
     .          nblon(mxbli)
      dimension iindex(intmax,6*nsub1+9),isav_pat(intmax,17),
     .          isav_pat_b(intmax,nsub1,6)
      dimension isav_blk(2*mxbli,17)
      dimension isav_prd(lbcprd,12)
      dimension isav_emb(lbcemb,12)
c
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /is_perbc/ is_prd(5),ie_prd(5),nbcprd
      common /is_patch/ is_pat(5),ie_pat(5),ipatch1st
      common /is_embedbc/ is_emb(5),ie_emb(5),nbcemb
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
      ierrflg = -1
c
c     initialize arrays required for patch boundaries
c
      do nnn = 1,intmax
         do lk = 1,17
            isav_pat(nnn,lk) = 0
         end do
         do lk = 1,nsub1
            do lm = 1,6
               isav_pat_b(nnn,lk,lm) = 0
            end do
         end do
      end do
c
c     initialize arrays required for 1:1 block boundaries
c
      do nnn = 1,2*mxbli
         do lk = 1,17
            isav_blk(nnn,lk) = 0
         end do
      end do
c
c     initialize arrays required for periodic block boundaries
c
      do nnn = 1,lbcprd
         do lk = 1,12
            isav_prd(nnn,lk) = 0
         end do
      end do
c
c     initialize arrays required for embedded block boundaries
c
      do nnn = 1,lbcemb
         do lk = 1,12
            isav_emb(nnn,lk) = 0
         end do
      end do
c
      icount_pat = 0
      icount_blk = 0
      icount_prd = 0
      icount_emb = 0
c
#if defined DIST_MPI
c
c     set baseline tag values
c
      itag_lw = 1
c
#endif
      do levl = 1,levt
         is_blk(levl) = icount_blk + 1
         is_pat(levl) = icount_pat + 1
         is_prd(levl) = icount_prd + 1
         is_emb(levl) = icount_emb + 1
         ie_blk(levl) = is_blk(levl) - 1
         ie_pat(levl) = is_pat(levl) - 1
         ie_prd(levl) = is_prd(levl) - 1
         ie_emb(levl) = is_emb(levl) - 1
         do 6909 nbl=1,nblock
            if (levl.ne.levelg(nbl)) go to 6909
#if defined DIST_MPI
            nd_srce = mblk2nd(nbl)
c
c           note: pre_patch, pre_blockbc, etc., are called
c           from the host, but need the lw array that is
c           set on the nodes (and is different from that on
c           the host). hence, the need to pass lw/lw_temp
c           back to the host
c
            if (myid.eq.nd_srce) then
               mytag = itag_lw + nbl
               call MPI_Send(lw, 65*maxbl, MPI_INTEGER, myhost,
     .                        mytag, mycomm, ierr)
            else
               if (myid.eq.myhost) then
                  mytag = itag_lw + nbl
                  call MPI_Recv(lw_temp, 65*maxbl, MPI_INTEGER,
     .                           nd_srce, mytag, mycomm, istat, ierr)
                  icount_pat1 = icount_pat
                  call pre_patch(nbl,lw_temp,icount_pat,ninter,
     .                            iindex,intmax,nsub1,isav_pat,
     .                            isav_pat_b,jjmax1,kkmax1,
     .                            iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                            ierrflg)
                  icount_blk1 = icount_blk
                  call pre_blockbc(nbl,lw_temp,icount_blk,idimg,
     .                             jdimg,kdimg,isav_blk,nblk,nbli,
     .                             limblk,isva,nblon,mxbli,nou,bou,nbuf,
     .                             ibufdim,myid,maxbl,ierrflg)
                  icount_prd1 = icount_prd
                  call pre_period(nbl,lw_temp,lw2,icount_prd,
     .                            maxbl,maxseg,lbcprd,
     .                            nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                            nbcidim,jbcinfo,kbcinfo,ibcinfo,
     .                            igridg,jdimg,kdimg,idimg,isav_prd,
     .                            is_prd,ie_prd,nbcprd,nou,bou,nbuf,
     .                            ibufdim,bcvali,bcvalj,bcvalk,myid,
     .                            nblg,maxgr,ierrflg)
                  icount_emb1 = icount_emb
                  call pre_embed(nbl,lw_temp,lw2,icount_emb,
     .                           maxbl,maxseg,lbcemb,
     .                           nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                           nbcidim,jbcinfo,kbcinfo,ibcinfo,iemg,
     .                           igridg,jdimg,kdimg,idimg,isav_emb,
     .                           is_emb,ie_emb,nbcemb,nou,bou,nbuf,
     .                           ibufdim,myid,maxgr,ierrflg)
                  if (icount_blk .gt. icount_blk1) then
                     ie_blk(levl) = icount_blk
                  end if
                  if (icount_pat .gt. icount_pat1) then
                     ie_pat(levl) = icount_pat
                  end if
                  if (icount_prd .gt. icount_prd1) then
                     ie_prd(levl) = icount_prd
                  end if
                  if (icount_emb .gt. icount_emb1) then
                     ie_emb(levl) = icount_emb
                  end if
               endif
            endif
#else
            icount_pat1 = icount_pat
            call pre_patch(nbl,lw,icount_pat,ninter,
     .                      iindex,intmax,nsub1,isav_pat,
     .                      isav_pat_b,jjmax1,kkmax1,
     .                      iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                      ierrflg)
            icount_blk1 = icount_blk
            call pre_blockbc(nbl,lw,icount_blk,idimg,
     .                       jdimg,kdimg,isav_blk,nblk,nbli,
     .                       limblk,isva,nblon,mxbli,nou,bou,nbuf,
     .                       ibufdim,myid,maxbl,ierrflg)
            icount_prd1 = icount_prd
            call pre_period(nbl,lw,lw2,icount_prd,
     .                      maxbl,maxseg,lbcprd,
     .                      nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                      nbcidim,jbcinfo,kbcinfo,ibcinfo,
     .                      igridg,jdimg,kdimg,idimg,isav_prd,
     .                      is_prd,ie_prd,nbcprd,nou,bou,nbuf,ibufdim,
     .                      bcvali,bcvalj,bcvalk,myid,nblg,maxgr,
     .                      ierrflg)
            icount_emb1 = icount_emb
            call pre_embed(nbl,lw,lw2,icount_emb,
     .                     maxbl,maxseg,lbcemb,
     .                     nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                     nbcidim,jbcinfo,kbcinfo,ibcinfo,iemg,
     .                     igridg,jdimg,kdimg,idimg,isav_emb,
     .                     is_emb,ie_emb,nbcemb,nou,bou,nbuf,ibufdim,
     .                     myid,maxgr,ierrflg)
            if (icount_blk .gt. icount_blk1) then
               ie_blk(levl) = icount_blk
            end if
            if (icount_pat .gt. icount_pat1) then
               ie_pat(levl) = icount_pat
            end if
            if (icount_prd .gt. icount_prd1) then
               ie_prd(levl) = icount_prd
            end if
            if (icount_emb .gt. icount_emb1) then
               ie_emb(levl) = icount_emb
            end if
#endif
 6909    continue
      end do
c
#if defined DIST_MPI
c
c     share information with the processors
c
      call MPI_Bcast (isav_pat, intmax*17, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (isav_pat_b, intmax*nsub1*6, MPI_INTEGER,
     .                myhost, mycomm, ierr)
      call MPI_Bcast (isav_blk, 2*mxbli*17, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (isav_prd, lbcprd*12, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (isav_emb, lbcemb*12, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (is_blk, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (ie_blk, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (is_pat, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (ie_pat, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (is_prd, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (ie_prd, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (nbcprd, 1, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (is_emb, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (ie_emb, 5, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (nbcemb, 1, MPI_INTEGER, myhost,
     .                mycomm, ierr)
#endif
c
      return
      end
