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
      subroutine prntcp(jdim,kdim,idim,cpi,cpj,cpk,q,nbl,maxbl,maxseg,
     .                  ibcinfo,jbcinfo,kbcinfo,nbci0,nbcj0,nbck0,
     .                  nbcidim,nbcjdim,nbckdim,thetay,mblk2nd,myid,
     .                  myhost,mycomm,irite)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Writes pressures on body (actually cell centers closest
c     to body) to output file for unsteady flow
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
      dimension itrnsfr(2)
#endif
c
      dimension q(jdim,kdim,idim,5),cpi(jdim*kdim),cpj(kdim*idim),
     .          cpk(jdim*idim),thetay(maxbl),mblk2nd(maxbl)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /conversion/ radtodeg
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset   = maxseg
      itag_size = 1
      itag_cp   = itag_size + ioffset
c
#endif
      alphd = radtodeg*(alpha+thetay(nbl))
c
      cpc   = 2.e0/(gamma*xmach*xmach)
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c 
c     pressures on i0 segments
c
      if (nbci0(nbl).gt.0) then
         do 10 ns=1,nbci0(nbl)
         if (abs(ibcinfo(nbl,ns,1,1)).eq.1005 .or. 
     .       abs(ibcinfo(nbl,ns,1,1)).eq.1006 .or.
     .       abs(ibcinfo(nbl,ns,1,1)).eq.2004 .or.
     .       abs(ibcinfo(nbl,ns,1,1)).eq.2014 .or.
     .       abs(ibcinfo(nbl,ns,1,1)).eq.2024 .or.
     .       abs(ibcinfo(nbl,ns,1,1)).eq.2034 .or.
     .       abs(ibcinfo(nbl,ns,1,1)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                jst  = ibcinfo(nbl,ns,2,1)
                jfn  = ibcinfo(nbl,ns,3,1)-1
                kst  = ibcinfo(nbl,ns,4,1)
                kfn  = ibcinfo(nbl,ns,5,1)-1
                jj   = 0
                jlen = jfn-jst+1
                klen = kfn-kst+1
                do 11 j=jst,jfn
                jj   = jj+1
                kk   = 0
                do 12 k=kst,kfn 
                kk   = kk+1
                izz  = (kk-1)*jlen + jj
                cpi(izz) = (q(j,k,1,5)/p0-1.e0)*cpc
   12           continue
   11          continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = jlen
                itrnsfr(2) = klen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr) 
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                jlen = itrnsfr(1)
                klen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpi,jlen*klen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpi,jlen*klen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,1,ns
                write(20,4)
                write(20,1)real(time),real(alphd),jlen,klen
                do 13 j=1,jlen
                write(20,2)(real(cpi((k-1)*jlen + j)),k=1,klen)
   13           continue
             end if
         end if
   10    continue
      end if
c 
c     pressures on idim segments
c
      if (nbcidim(nbl).gt.0) then
         do 20 ns=1,nbcidim(nbl)
         if (abs(ibcinfo(nbl,ns,1,2)).eq.1005 .or. 
     .       abs(ibcinfo(nbl,ns,1,2)).eq.1006 .or.
     .       abs(ibcinfo(nbl,ns,1,2)).eq.2004 .or.
     .       abs(ibcinfo(nbl,ns,1,2)).eq.2014 .or.
     .       abs(ibcinfo(nbl,ns,1,2)).eq.2024 .or.
     .       abs(ibcinfo(nbl,ns,1,2)).eq.2034 .or.
     .       abs(ibcinfo(nbl,ns,1,2)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                jst  = ibcinfo(nbl,ns,2,2)
                jfn  = ibcinfo(nbl,ns,3,2)-1
                kst  = ibcinfo(nbl,ns,4,2)
                kfn  = ibcinfo(nbl,ns,5,2)-1
                jj   = 0
                jlen = jfn-jst+1
                klen = kfn-kst+1
                do 21 j=jst,jfn
                jj   = jj+1
                kk   = 0
                do 22 k=kst,kfn
                kk   = kk+1
                izz  = (kk-1)*jlen + jj
                cpi(izz) = (q(j,k,idim1,5)/p0-1.e0)*cpc
   22           continue
   21           continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = jlen
                itrnsfr(2) = klen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                jlen = itrnsfr(1)
                klen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpi,jlen*klen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpi,jlen*klen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,2,ns
                write(20,4)
                write(20,1)real(time),real(alphd),jlen,klen
                do 23 j=1,jlen
                write(20,2)(real(cpi((k-1)*jlen + j)),k=1,klen)
   23           continue
             end if
         end if
   20    continue
      end if
c
c     pressures on j0 segments
c
      if (nbcj0(nbl).gt.0) then
         do 30 ns=1,nbcj0(nbl)
         if (abs(jbcinfo(nbl,ns,1,1)).eq.1005 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.1006 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.2004 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.2014 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.2024 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.2034 .or.
     .       abs(jbcinfo(nbl,ns,1,1)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                ist  = jbcinfo(nbl,ns,2,1)
                ifn  = jbcinfo(nbl,ns,3,1)-1
                kst  = jbcinfo(nbl,ns,4,1)
                kfn  = jbcinfo(nbl,ns,5,1)-1
                ii   = 0
                ilen = ifn-ist+1
                klen = kfn-kst+1
                do 31 i=ist,ifn
                ii   = ii+1
                kk   = 0
                do 32 k=kst,kfn
                kk   = kk+1
                izz  = (ii-1)*klen + kk
                cpj(izz) = (q(1,k,i,5)/p0-1.e0)*cpc
   32           continue
   31           continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = ilen
                itrnsfr(2) = klen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                ilen = itrnsfr(1)
                klen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpj,ilen*klen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpj,ilen*klen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,3,ns
                write(20,5)
                write(20,1)real(time),real(alphd),klen,ilen
                do 33 k=1,klen
                write(20,2)(real(cpj((i-1)*klen + k)),i=1,ilen)
   33           continue
             end if
         end if
   30    continue
      end if
c
c     pressures on jdim segments
c
      if (nbcjdim(nbl).gt.0) then
         do 40 ns=1,nbcjdim(nbl)
         if (abs(jbcinfo(nbl,ns,1,2)).eq.1005 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.1006 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.2004 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.2014 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.2024 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.2034 .or.
     .       abs(jbcinfo(nbl,ns,1,2)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                ist  = jbcinfo(nbl,ns,2,2)
                ifn  = jbcinfo(nbl,ns,3,2)-1
                kst  = jbcinfo(nbl,ns,4,2)
                kfn  = jbcinfo(nbl,ns,5,2)-1
                ii   = 0
                ilen = ifn-ist+1
                klen = kfn-kst+1
                do 41 i=ist,ifn
                ii   = ii+1
                kk   = 0
                do 42 k=kst,kfn
                kk   = kk+1
                izz  = (ii-1)*klen + kk
                cpj(izz) = (q(jdim1,k,i,5)/p0-1.e0)*cpc
   42           continue
   41           continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = ilen
                itrnsfr(2) = klen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                ilen = itrnsfr(1)
                klen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpj,ilen*klen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpj,ilen*klen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,4,ns
                write(20,5)
                write(20,1)real(time),real(alphd),klen,ilen
                do 43 k=1,klen
                write(20,2)(real(cpj((i-1)*klen + k)),i=1,ilen)
   43           continue
             end if
         end if
   40    continue
      end if
c
c     pressures on k0 segments
c
      if (nbck0(nbl).gt.0) then
         do 50 ns=1,nbck0(nbl)
         if (abs(kbcinfo(nbl,ns,1,1)).eq.1005 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.1006 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.2004 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.2014 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.2024 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.2034 .or.
     .       abs(kbcinfo(nbl,ns,1,1)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                ist  = kbcinfo(nbl,ns,2,1)
                ifn  = kbcinfo(nbl,ns,3,1)-1
                jst  = kbcinfo(nbl,ns,4,1)
                jfn  = kbcinfo(nbl,ns,5,1)-1
                ii   = 0
                ilen = ifn-ist+1
                jlen = jfn-jst+1
                do 51 i=ist,ifn
                ii   = ii+1
                jj   = 0
                do 52 j=jst,jfn
                jj   = jj+1
                izz  = (ii-1)*jlen + jj
                cpk(izz) = (q(j,1,i,5)/p0-1.e0)*cpc
   52           continue
   51           continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = ilen
                itrnsfr(2) = jlen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                ilen = itrnsfr(1)
                jlen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpk,ilen*jlen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpk,ilen*jlen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,5,ns
                write(20,6)
                write(20,1)real(time),real(alphd),jlen,ilen
                do 53 i=1,ilen
                write(20,2)(real(cpk((i-1)*jlen + j)),j=1,jlen)
   53           continue
             end if
         end if
   50    continue
      end if
c
c     pressures on kdim segments
c
      if (nbckdim(nbl).gt.0) then
         do 60 ns=1,nbckdim(nbl)
         if (abs(kbcinfo(nbl,ns,1,2)).eq.1005 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.1006 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.2004 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.2014 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.2024 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.2034 .or.
     .       abs(kbcinfo(nbl,ns,1,2)).eq.2016) then
             if (myid.eq.mblk2nd(nbl)) then
                ist  = kbcinfo(nbl,ns,2,2)
                ifn  = kbcinfo(nbl,ns,3,2)-1
                jst  = kbcinfo(nbl,ns,4,2)
                jfn  = kbcinfo(nbl,ns,5,2)-1
                ii   = 0
                ilen = ifn-ist+1
                jlen = jfn-jst+1
                do 61 i=ist,ifn
                ii   = ii+1
                jj   = 0
                do 62 j=jst,jfn
                jj   = jj+1
                izz  = (ii-1)*jlen + jj
                cpk(izz) = (q(j,kdim1,i,5)/p0-1.e0)*cpc
   62           continue
   61           continue
             end if
#if defined DIST_MPI
c
             nd_srce = mblk2nd(nbl)
             mytag1  = itag_size + ns
             mytag2  = itag_cp   + ns
c
             if (myid.eq.nd_srce) then
                itrnsfr(1) = ilen
                itrnsfr(2) = jlen
                call MPI_Send(itrnsfr,2,MPI_INTEGER,myhost,
     .                        mytag1,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(itrnsfr,2,MPI_INTEGER,nd_srce,
     .                        mytag1,mycomm,istat,ierr)
                ilen = itrnsfr(1)
                jlen = itrnsfr(2)
             end if
c
             if (myid.eq.nd_srce) then
                call MPI_Send(cpk,ilen*jlen,MY_MPI_REAL,
     .                        myhost,mytag2,mycomm,ierr)
             end if
             if (myid.eq.myhost) then
                call MPI_Recv(cpk,ilen*jlen,MY_MPI_REAL,
     .                        nd_srce,mytag2,mycomm,istat,ierr)
             end if
c
#endif
             if (myid.eq.myhost) then
                write(20,3) nbl,6,ns
                write(20,6)
                write(20,1)real(time),real(alphd),jlen,ilen
                do 63 i=1,ilen
                write(20,2)(real(cpk((i-1)*jlen + j)),j=1,jlen)
   63           continue
             end if
         end if
   60    continue
      end if
c
    1 format(2f10.5,2i10)
    2 format(10f10.6)
    3 format(10h     block,10h      face,10h   segment,/,3i10)
    4 format(10h      time,10h     alpha,10h      jlen,10h      klen)
    5 format(10h      time,10h     alpha,10h      klen,10h      ilen)
    6 format(10h      time,10h     alpha,10h      jlen,10h      ilen)
c
      return
      end
