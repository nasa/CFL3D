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
      subroutine rrest(nbl,jdim,kdim,idim,q,qj0,qk0,qi0,ncycmax,ntr,rms,
     .     clw,cdw,cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,
     .     n_clcd, clcd, nblocks_clcd, blocks_clcd,
     .     fmdotw,cftmomw,cftpw,cftvw,cfttotw,rmstr,
     .     nneg,iskip,vist3d,tursav,
     .     smin,xjb,xkb,blnum,cmuv,maxbl,mblk2nd,
     .     myid,myhost,mycomm,nou,bou,nbuf,ibufdim,
     .     igrid,wk,idima,jdima,kdima,vj0,vk0,vi0,
     .     tj0,tk0,ti0,qavg,q2avg,nummem)
c     
c     $Id$
c     
c***********************************************************************
c     Purpose:  Read the restart file for a block. iskip = 1 means this
c     is the first block to be read in from the restart file;
c     all other blocks have iskip = 0
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
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4),
     .     qj0(kdim,idim-1,5,4), qk0(jdim,idim-1,5,4),
     .     vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .     vi0(jdim,kdim,1,4),tj0(kdim,idim-1,nummem,4),
     .     tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
      integer blocks_clcd
      integer stats
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),cdpw(ncycmax),
     .     cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .     cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .     clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd),
     .     fmdotw(ncycmax),cftmomw(ncycmax),cftpw(ncycmax),
     .     cftvw(ncycmax),cfttotw(ncycmax),
     .     rmstr(ncycmax,nummem),nneg(ncycmax,nummem)
      dimension titlw(20),mblk2nd(maxbl)
      dimension tursav(jdim,kdim,idim,nummem),vist3d(jdim,kdim,idim),
     .     smin(jdim-1,kdim-1,idim-1),xjb(jdim-1,kdim-1,idim-1),
     .     xkb(jdim-1,kdim-1,idim-1),blnum(jdim-1,kdim-1,idim-1),
     .     cmuv(jdim-1,kdim-1,idim-1)
      dimension qavg(jdim,kdim,idim,5)
      dimension q2avg(jdim,kdim,idim,5)
      dimension wk(idima+1,jdima+1,kdima+1)
      dimension irinddata(6)
c     
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /twod/ i2d
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .     idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .     iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wrestq/ irest,irest2
      common /reyue/ reue,tinf,ivisc(3)
      common /sminn/ isminc,ismincforce
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .     tur10(7)
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /conversion/ radtodeg
      common /ghost/ irghost,iwghost
      common /igrdtyp/ ip3dgrd,ialph

      common /random_input/ randomize
      real, allocatable :: harvest(:)
c     
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c     
#if defined DIST_MPI
c     set baseline tag values
c     
      ioffset      = maxbl
      itag_alpha   = 1
      itag_ntt     = itag_alpha   + ioffset
      itag_rms     = itag_ntt     + ioffset
      itag_clw     = itag_rms     + ioffset
      itag_cdw     = itag_clw     + ioffset
      itag_cdpw    = itag_cdw     + ioffset
      itag_cdvw    = itag_cdpw    + ioffset
      itag_cxw     = itag_cdvw    + ioffset
      itag_cyw     = itag_cxw     + ioffset
      itag_czw     = itag_cyw     + ioffset
      itag_cmxw    = itag_czw     + ioffset
      itag_cmyw    = itag_cmxw    + ioffset
      itag_cmzw    = itag_cmyw    + ioffset
      itag_fmdotw  = itag_cmzw    + ioffset
      itag_cftmomw = itag_fmdotw  + ioffset
      itag_cftpw   = itag_cftmomw + ioffset
      itag_cftvw   = itag_cftpw   + ioffset
      itag_cfttotw = itag_cftvw   + ioffset
      itag_q       = itag_cfttotw + ioffset
      itag_xmachw  = itag_q       + ioffset
      itag_reuew   = itag_xmachw  + ioffset
      itag_time    = itag_reuew   + ioffset
      itag_iv1     = itag_time    + ioffset
      itag_iv2     = itag_iv1     + ioffset
      itag_iv3     = itag_iv2     + ioffset
      itag_rmstr   = itag_iv3     + ioffset
      itag_nneg    = itag_rmstr   + ioffset
      itag_qv      = itag_nneg    + ioffset
      itag_qt      = itag_qv      + ioffset
      itag_smin    = itag_qt      + ioffset
      itag_xjb     = itag_smin    + ioffset
      itag_xkb     = itag_xjb     + ioffset
      itag_blnum   = itag_xkb     + ioffset
      itag_cmuv    = itag_blnum   + ioffset
      itag_qj0     = itag_cmuv    + ioffset
      itag_qk0     = itag_qj0     + ioffset
      itag_qi0     = itag_qk0     + ioffset
      itag_vj0     = itag_qi0     + ioffset
      itag_vk0     = itag_vj0     + ioffset
      itag_vi0     = itag_vk0     + ioffset
      itag_tj0     = itag_vi0     + ioffset
      itag_tk0     = itag_tj0     + ioffset
      itag_ti0     = itag_tk0     + ioffset
      itag_qavg    = itag_ti0     + ioffset
      itag_xnum    = itag_qavg    + ioffset
      itag_q2avg   = itag_xnum    + ioffset
      itag_xnum2   = itag_q2avg   + ioffset
      itag_clcd    = itag_xnum2   + ioffset
c     
      jki   = jdim*kdim*idim
      jkim  = jdim1*kdim1*idim1
      jki5  = 5*jki
      jkim5 = 5*jkim
      jki2  = 2*jki
      jk20  = jdim*kdim*20
      ki20  = kdim*idim1*20
      ji20  = jdim*idim1*20
      jk4   = jdim*kdim*4
      ki4   = kdim*idim1*4
      ji4   = jdim*idim1*4
      jk8   = jdim*kdim*8
      ki8   = kdim*idim1*8
      ji8   = jdim*idim1*8
c     
      if (myid.eq.myhost) then
c     
#endif
         if (icgns .eq. 1) then
#if defined CGNS
            write(901,'('' cgns read in rrest'')')
            write(901,'(''  ...reading conserved variables and'',
     +           '' translating to primitive'')')
            write(901,'(''(  this may result in machine-order '',
     +           '' differences compared to a non-cgns run)'')')
            call getsolinfo(iccg,ibase,igrid,nsoluse,irinddata)
            irind=1
            jrind=1
            krind=1
c     Note: the following assumes irinddata(2)=(1), (4)=(3),
c     (6)=(5) - (currently checked for in getsolinfo):
            if (irinddata(1) .eq. 1) irind=0
            if (irinddata(3) .eq. 1) jrind=0
            if (irinddata(5) .eq. 1) krind=0
            call readsoln(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +         jdim,kdim,wk,nsoluse,irind,jrind,krind,i2d,ialph,q,iprim)
            if (iprim .eq. 0) then
c     convert q to primitive (they are in CGNS file as conserved)
               do i=1,idim1
                  do k=1,kdim1
                     do j=1,jdim1
                        q(j,k,i,2)=q(j,k,i,2)/q(j,k,i,1)
                        q(j,k,i,3)=q(j,k,i,3)/q(j,k,i,1)
                        q(j,k,i,4)=q(j,k,i,4)/q(j,k,i,1)
                        q(j,k,i,5)=gm1*(q(j,k,i,5)
     $                       -0.5*q(j,k,i,1)*(q(j,k,i,2)**2+
     +                       q(j,k,i,3)**2+q(j,k,i,4)**2))
                     enddo
                  enddo
               enddo
            end if
            if (irghost .ne. 0) then
               write(901,'('' cgns reading specific BC values'')')
               call readbcs(iccg,ibase,igrid,idim,jdim,kdim,i2d,
     +              qj0,qk0,qi0,vj0,vk0,vi0,tj0,tk0,ti0,nummem)
            end if
            if (iskip .gt. 0) then
c             readhist called with all 7 turb histories, but dummy
c             variables inserted for those not used
              if (nummem .eq. 2) then
               call readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     +              cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     +              cftpw,cftvw,cfttotw,rmstr(1,1),rmstr(1,2),
     +              rmstr(1,1),rmstr(1,1),rmstr(1,1),rmstr(1,1),
     +              rmstr(1,1),nneg(1,1),nneg(1,2),nneg(1,1),
     +              nneg(1,1),nneg(1,1),nneg(1,1),nneg(1,1),nummem)
              else if (nummem .eq. 3) then
               call readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     +              cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     +              cftpw,cftvw,cfttotw,rmstr(1,1),rmstr(1,2),
     +              rmstr(1,3),rmstr(1,1),rmstr(1,1),rmstr(1,1),
     +              rmstr(1,1),nneg(1,1),nneg(1,2),nneg(1,3),
     +              nneg(1,1),nneg(1,1),nneg(1,1),nneg(1,1),nummem)
              else
               call readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     +              cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     +              cftpw,cftvw,cfttotw,rmstr(1,1),rmstr(1,2),
     +              rmstr(1,3),rmstr(1,4),rmstr(1,5),rmstr(1,6),
     +              rmstr(1,7),nneg(1,1),nneg(1,2),nneg(1,3),
     +              nneg(1,4),nneg(1,5),nneg(1,6),nneg(1,7),nummem)
              end if
            end if
            call readref(iccg,ibase,xmachw,reuew,
     +           rho0w,c0w,vk0w,xlength0w,tinfw,alphw)
            if (rho0w.ne.1. .or. c0w.ne.1. .or. xlength0w.ne.1.) then
               write(901,'('' ERROR, cgns restart file has different'',
     +              '' nondimensionalization'')')
               write(901,'('' rho0w,c0w,xlength0w='',3e15.5)')
     +              rho0w,c0w,xlength0w
               write(901,'('' Stopping.'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            if (abs(vk0w-1.) .gt. 1.e-4) then
               write(901,'('' ERROR, cgns restart file has different'',
     +              '' nondimensionalization'')')
               write(901,'('' vk0w='',e15.5)') vk0w
               write(901,'('' Stopping.'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            if(tinfw .ne. tinf) then
               write(901,'('' WARNING:  changing tinf!'')')
               write(901,'(''    old='',f12.5,'' new='',f12.5)')
     +              tinfw,tinf
            end if
            call readtime(iccg,ibase,time)
#endif
         else
            read(2) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
c     write(15,97)jt,kt,it
c     97 format(1x,15hrestarting on a,i5,3h x ,i5,3h x ,i5,5h grid)
            if (jt.ne.jdim .or. kt.ne.kdim .or. it.ne.idim) then
               write(11,*) ' stopping.inconsistent restart grid indices'
               write(11,*) ' block indices  idim,jdim,kdim =',idim,jdim,
     $              kdim
               write(11,*) ' restart indices  idim,jdim,kdim =',it,jt,kt
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c     
c     write(15,100) real(titlw),real(xmachw),jt,kt,it,real(alphw),
c     real(reuew),ntr,real(time)
c     100 format(20a4/1x,e12.5,3i5,2e12.5,i5,e12.5)
c     
         ntt = ntr
c     
c     use previous alpha as starting point if alpha-iteration is used
c     
         if (ialphit.gt.0) then
            alpha = alphw/radtodeg
c     reset freestream for new alpha on host
            call init_mast
         end if
c     
         if (icgns .ne. 1) then
            if (ntr.le.ncycmax) then
c     
c     Convergence data (residual,force coefficients, mass flow, etc.)
c     
               if (iskip.gt.0) then
                  read(2) (rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .                 (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .                 (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .                 (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .                 (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .                 (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .                 (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .                 (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
               end if
c     
            else
c     
               if (irest2 .eq. 1) then
                  if (iskip.gt.0) then
                     read(2) (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr),
     .                    (dum,     n=1,ntr),(dum,     n=1,ntr)
                  end if
               else
c     
                  write(11,1239)
 1239             format(/,1x,11hstopping...,
     .                 40hprevious number of iterations computed >,
     .                 1x,18h dimension ncycmax)
                  write(11,*)' ntr,ncycmax = ',ntr,ncycmax
                  write(11,*)' increase value of ncycmax to at LEAST ',
     .                 ntr+ncycmax
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
c     
               end if
            end if
         end if
#if defined DIST_MPI
c     
      end if
c     
c     send/receive data to/on the appropriate node
c     
      nd_dest = mblk2nd(nbl)
c     
      if (ialphit .gt. 0) then
         mytag = itag_alpha + nbl
         if (myid .eq. myhost) then
            call MPI_Send(alpha,1,MY_MPI_REAL,
     .           nd_dest,mytag,mycomm,ierr)
         else if (myid .eq. mblk2nd(nbl)) then
            call MPI_Recv(alpha,1,MY_MPI_REAL,
     .           myhost,mytag,mycomm,istat,ierr)
c     reset freestream for new alpha on compute  processors
            call init_mast
         end if
      end if
c     
      mytag = itag_ntt + nbl
      if (myid .eq. myhost) then
         call MPI_Send(ntt,1,MPI_INTEGER,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(ntt,1,MPI_INTEGER,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_rms + nbl
      if (myid .eq. myhost) then
         call MPI_Send(rms,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(rms,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_clw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(clw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(clw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cdw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cdw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cdw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cdpw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cdpw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cdpw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cdvw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cdvw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cdvw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cxw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cxw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cxw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cyw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cyw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cyw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_czw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(czw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(czw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cmxw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cmxw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cmxw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cmyw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cmyw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cmyw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cmzw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cmzw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cmzw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_fmdotw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(fmdotw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(fmdotw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cftmomw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cftmomw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cftmomw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cftpw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cftpw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cftpw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cftvw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cftvw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cftvw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_cfttotw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(cfttotw,ntt,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(cfttotw,ntt,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
#endif
c     
c     Primative variables (rho,u,v,w,p)
c     
#if defined DIST_MPI
      if (myid.eq.myhost) then
c     
#endif
         if (icgns .ne. 1) then
            read(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),
     $           i=1,idim1),l=1,5)
            if (irghost .ne. 0)
     .           read(2) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),
     $           l=1,5),m=1,4),
     .           ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     .           ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
         end if
c     
#if defined DIST_MPI
c     
c     fill in edge values of q array for safety before passing
c     the data to the appropriate node
c     
         do l=1,5
            do i=1,idim1
               k = kdim
               do j=1,jdim1
                  q(j,k,i,l) = qiv(l)
               end do
               j=jdim
               do k=1,kdim
                  q(j,k,i,l) = qiv(l)
               end do
            end do
            i = idim
            do j=1,jdim
               do k=1,kdim
                  q(j,k,i,l) = qiv(l)
               end do
            end do
         end do
c     
      end if
c     
c     send/receive data to/on the appropriate node
c     
      nd_dest = mblk2nd(nbl)
c     
      mytag = itag_q + nbl
      if (myid .eq. myhost) then
         call MPI_Send(q,jki5,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(q,jki5,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      if (irghost .ne.0) then
         mytag = itag_qj0 + nbl
         if (myid .eq. myhost) then
            call MPI_Send(qj0,ki20,MY_MPI_REAL,
     .           nd_dest,mytag,mycomm,ierr)
         else if (myid .eq. mblk2nd(nbl)) then
            call MPI_Recv(qj0,ki20,MY_MPI_REAL,
     .           myhost,mytag,mycomm,istat,ierr)
         end if
c     
         mytag = itag_qk0 + nbl
         if (myid .eq. myhost) then
            call MPI_Send(qk0,ji20,MY_MPI_REAL,
     .           nd_dest,mytag,mycomm,ierr)
         else if (myid .eq. mblk2nd(nbl)) then
            call MPI_Recv(qk0,ji20,MY_MPI_REAL,
     .           myhost,mytag,mycomm,istat,ierr)
         end if
c     
         mytag = itag_qi0 + nbl
         if (myid .eq. myhost) then
            call MPI_Send(qi0,jk20,MY_MPI_REAL,
     .           nd_dest,mytag,mycomm,ierr)
         else if (myid .eq. mblk2nd(nbl)) then
            call MPI_Recv(qi0,jk20,MY_MPI_REAL,
     .           myhost,mytag,mycomm,istat,ierr)
         end if
      end if
c     
      mytag = itag_xmachw + nbl
      if (myid .eq. myhost) then
         call MPI_Send(xmachw,1,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(xmachw,1,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_reuew + nbl
      if (myid .eq. myhost) then
         call MPI_Send(reuew,1,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(reuew,1,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      mytag = itag_time + nbl
      if (myid .eq. myhost) then
         call MPI_Send(time,1,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(time,1,MY_MPI_REAL,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
      if (myid.eq.mblk2nd(nbl)) then
#endif
         n   = jdim*kdim
         if (xmachw.ne.xmach) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''WARNING: fix-up taken for Mach '',
     .           ''number inconsistancy in restart file'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  scaling u,v,w by M(new)/M(old);'',
     .           '' cp held fixed: M(new), M(old)= '',e12.4,e12.4)')
     .           real(xmach),real(xmachw)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            do 1011 i=1,idim
c     
c     u/cref,v/cref,w/cref
c     
c     dir$ ivdep
               do 1000 izz=1,n
                  q(izz,1,i,2) = q(izz,1,i,2)*xmach/xmachw
                  q(izz,1,i,3) = q(izz,1,i,3)*xmach/xmachw
                  q(izz,1,i,4) = q(izz,1,i,4)*xmach/xmachw
c     
c     cp=constant
c     
                  q(izz,1,i,5) =  q(izz,1,i,5)*gamma
                  q(izz,1,i,5) = (q(izz,1,i,5)-1.e0)*xmach*
     $                 xmach/xmachw/xmachw
                  q(izz,1,i,5) = (q(izz,1,i,5)+1.e0)/gamma
 1000          continue
c     
c     check for negative density,pressure
c     
               do 3235 j=1,jdim
                  do 3234 k=1,kdim
                     if (real(q(j,k,i,1)).le.0. .or.
     $                    real(q(j,k,i,5)).le.0.) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*) ' neg. density,pressure'
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*) real(q(j,k,i,1)),
     $                       real(q(j,k,i,5))
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),*) i,j,k
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
 3234             continue
 3235          continue
c     
 1011       continue
         end if
         if (reuew.ne.reue) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''WARNING:  ReUe inconsistency in'',
     .           '' restart file data Solution will proceed with '',
     $           ''new value!'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  ReUe(new), ReUe(old)= '',
     .           e12.4,e12.4)') real(reue),real(reuew)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
         end if
c
         if (real(randomize) .gt. 1.e-12) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'('' adding random component to'',
     .           '' restart values'')')
           jki_size=(jdim-1)*(kdim-1)*(idim-1)
           allocate( harvest(jki_size),stat=stats )
           call umalloc(jki_size,0,'harvest',memuse,stats)
           call random_seed
           call random_number(harvest)
           irnd_count=0
           do i=1,idim-1
             do k=1,kdim-1
               do j=1,jdim-1
                 irnd_count=irnd_count+1
                 rnd_nbr=harvest(irnd_count)-0.5
                 q(j,k,i,1)=q(j,k,i,1)*(1.0+rnd_nbr*randomize)
                 q(j,k,i,2)=q(j,k,i,2)*(1.0+rnd_nbr*randomize)
                 q(j,k,i,3)=q(j,k,i,3)*(1.0+rnd_nbr*randomize)
                 q(j,k,i,4)=q(j,k,i,4)*(1.0+rnd_nbr*randomize)
                 q(j,k,i,5)=q(j,k,i,5)*(1.0+rnd_nbr*randomize)
               enddo
             enddo
           enddo
         end if
c
#if defined DIST_MPI
c     
      end if
#endif
c     
c     Turbulence quantities
c     
#if defined DIST_MPI
      if (myid .eq. myhost) then
#endif
         if (icgns .eq. 1) then
#if defined CGNS
            call readeqn(iccg,ibase,igrid,iv1,iv2,iv3,gammaw,
     +           prw,prtw,cbarw)
            if(gammaw .ne. gamma) then
               write(901,'('' WARNING:  changing gamma!'')')
               write(901,'(''    old='',f12.5,'' new='',f12.5)')
     +              gammaw,gamma
            end if
            if(prw .ne. pr) then
               write(901,'('' WARNING:  changing pr!'')')
               write(901,'(''    old='',f12.5,'' new='',f12.5)')
     +              prw,pr
            end if
            if(prtw .ne. prt) then
               write(901,'('' WARNING:  changing prt!'')')
               write(901,'(''    old='',f12.5,'' new='',f12.5)')
     +              prtw,prt
            end if
            if(cbarw .ne. cbar) then
               write(901,'('' WARNING:  changing cbar!'')')
               write(901,'(''    old='',f12.5,'' new='',f12.5)')
     +              cbarw,cbar
            end if
#endif
         else
            read(2) iv1,iv2,iv3
         end if
#if defined DIST_MPI
      end if
c     
      mytag = itag_iv1 + nbl
      if (myid .eq. myhost) then
         call MPI_Send(iv1,1,MPI_INTEGER,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(iv1,1,MPI_INTEGER,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
      mytag = itag_iv2 + nbl
      if (myid .eq. myhost) then
         call MPI_Send(iv2,1,MPI_INTEGER,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(iv2,1,MPI_INTEGER,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
      mytag = itag_iv3 + nbl
      if (myid .eq. myhost) then
         call MPI_Send(iv3,1,MPI_INTEGER,
     .        nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(iv3,1,MPI_INTEGER,
     .        myhost,mytag,mycomm,istat,ierr)
      end if
c     
#endif
      if (iv1.ne.ivisc(1).or.iv2.ne.ivisc(2).or.
     .     iv3.ne.ivisc(3)) then
#if defined DIST_MPI
         if (myid.eq.mblk2nd(nbl)) then
#endif
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' WARNING: restarting block'',i6,
     .           '' with DIFFERENT'',
     .           '' ivisc: old='',3i3,'', new='',3i3)') nbl,iv1,iv2,
     .           iv3,ivisc(1),ivisc(2),ivisc(3)
#if defined DIST_MPI
         end if
#endif
      end if
      ivmxold=max(iv1,iv2)
      ivmxold=max(ivmxold,iv3)
      iread = 1
      if (ivmxold.ne.ivmx) iread=0
      if ((ivmxold.ge.6.and.ivmxold.le.15) .and.
     .     (ivmx.ge.6.and.ivmx.le.15)) iread = 1
      if ((ivmxold.eq.4 .and. ivmx.eq.25) .or.
     .    (ivmxold.eq.25.and. ivmx.eq. 4)) iread=1
c
c   Note: if you try to allow restarts that change to/from type 30,
c   it may screw up the restart for time accurate: see rrestg.F
      if (real(dt) .gt. 0.) then
      if (ivmxold.eq.30 .and. ivmx.ne.30) then
        write(11,*) ' old ivisc=30; new must also be... stopping'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      if (ivmxold.ne.30 .and. ivmx.eq.30) then
        write(11,*) ' new ivisc=30; old must also be... stopping'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      end if
c   Note: if you try to allow restarts that change to/from type 40,
c   it may screw up the restart for time accurate: see rrestg.F
      if (real(dt) .gt. 0.) then
      if (ivmxold.eq.40 .and. ivmx.ne.40) then
        write(11,*) ' old ivisc=40; new must also be... stopping'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      if (ivmxold.ne.40 .and. ivmx.eq.40) then
        write(11,*) ' new ivisc=40; old must also be... stopping'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      end if
c   Need special checks when ivisc=72 is involved:
      if (iread .eq. 0) then
        if(ivmx==72) then
c         Note: this will only work if not time-accurate: see rrestg.F
          if (ivmxold/=6 .and. ivmxold/=7 .and. ivmxold/=8 .and.
     +        ivmxold/=12.and. ivmxold/=14) then
            write(11,*)"cannot restart ivisc=72 from turb model ",
     +        ivmxold
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' ERROR: cannot restart ivisc=72'',
     +        '' from turb model'',i6)') ivmxold
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          end if
        end if
        if(ivmxold==72 .and. ivmx/=72) then
          write(11,*)"cannot restart from turb model 72"
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' ERROR: cannot restart'',
     +      '' from turb model 72'')')
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
      end if
c     
c     Turbulence model convergence data
c     
      if (iskip.gt.0) then
         if (irest2 .eq. 1) then
#if defined DIST_MPI
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) (dum,n=1,ntr),(dum,n=1,ntr),
     .                 (idum,n=1,ntr),(idum,n=1,ntr)
                  if (iv1.ge.30 .or. iv2.ge.30 .or. iv3.ge.30) then
                    read(2) nummem_read
                    do l=3,nummem_read
                      read(2) (dum,n=1,ntt),(idum,n=1,ntt)
                    enddo
                  end if
               end if
#if defined DIST_MPI
            end if
#endif
         else
#if defined DIST_MPI
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) (rmstr(n,1),n=1,ntr),(rmstr(n,2),n=1,ntr),
     .                 (nneg(n,1), n=1,ntr),(nneg(n,2), n=1,ntr)
                  if (iv1.ge.30 .or. iv2.ge.30 .or. iv3.ge.30) then
                    read(2) nummem_read
                    do l=3,nummem_read
                      read(2) (rmstr(n,l),n=1,ntt),(nneg(n,l),n=1,ntt)
                    enddo
                  end if
               end if
#if defined DIST_MPI
            end if
c     
            mytag = itag_rmstr + nbl
            if (myid .eq. myhost) then
               call MPI_Send(rmstr,ntt*nummem,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(rmstr,ntt*nummem,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
            mytag = itag_nneg + nbl
            if (myid .eq. myhost) then
               call MPI_Send(nneg,ntt*nummem,MPI_INTEGER,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(nneg,ntt*nummem,MPI_INTEGER,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
#endif
         end if
      end if
c     
      if (icgns .eq. 1 .and. myid.eq.myhost) then
#if defined CGNS
         call readturb(iccg,ibase,igrid,nsoluse,idima,jdima,kdima,
     +        idim,jdim,kdim,wk,iv1,iv2,iv3,iread,irind,jrind,krind,
     +        i2d,ivmx,xmachw,reuew,vist3d,tursav,smin,xjb,xkb,
     +        tursav(1,1,1,2),blnum,cmuv,nummem)
#endif
      end if
      if (iread .eq. 0) then
c     
c     If diff turb model, DO NOT restart turb soln from quantities on 
c     binary restart file
c     (isminc=0 do not compute smin, do not call initvist
c     1 compute smin, call initvist
c     2 compute smin, do not call initvist
c     3 do not compute smin, call initvist)
c     Need to call initvist if start from scratch or restart from a
c     different turb model.  Need to call smin if start from scratch
c     or restart from a non-field-eqn turb model or if restart but 
c     mesh sequencing up to a finer level.
c     
         if (iv1.ge.2 .or. iv2.ge.2 .or. iv3.ge.2) then
#if defined DIST_MPI
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                  if (irghost .ne. 0)
     .                 read(2) ((((dum,j=1,jdim),k=1,kdim),l=1,1),
     $                 m=1,4),
     .                 ((((dum,k=1,kdim),i=1,idim1),l=1,1),m=1,4),
     .                 ((((dum,j=1,jdim),i=1,idim1),l=1,1),m=1,4)
               end if
#if defined DIST_MPI
            end if
#endif
         end if
c     
         if (iv1.ge.4 .or. iv2.ge.4 .or. iv3.ge.4) then
#if defined DIST_MPI
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     .                 m=1,2)
                  if (irghost .ne. 0)
     .                 read(2) ((((dum,j=1,jdim),k=1,kdim),l=1,2),
     $                 m=1,4),
     .                 ((((dum,k=1,kdim),i=1,idim1),l=1,2),m=1,4),
     .                 ((((dum,j=1,jdim),i=1,idim1),l=1,2),m=1,4)
                  read(2) (((smin(j,k,i),j=1,jdim1),k=1,kdim1),
     $                 i=1,idim1)
               end if
#if defined DIST_MPI
            end if
c     
            mytag = itag_smin + nbl
            if (myid .eq. myhost) then
               call MPI_Send(smin,jkim,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(smin,jkim,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
#endif
            if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4 .or.
     .          iv1.eq.25.or. iv2.eq.25.or. iv3.eq.25) then
#if defined DIST_MPI
               if (myid.eq.myhost) then
#endif
                  if (icgns .ne. 1) then
                     read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                     read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                     read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                  end if
#if defined DIST_MPI
               end if
#endif
            end if
            if (iv1.eq.8 .or. iv2.eq.8 .or. iv3.eq.8  .or.
     .           iv1.eq.9 .or. iv2.eq.9 .or. iv3.eq.9  .or.
     .           iv1.eq.13.or. iv2.eq.13.or. iv3.eq.13 .or.
     .           iv1.eq.14.or. iv2.eq.14.or. iv3.eq.14) then
#if defined DIST_MPI
               if (myid.eq.myhost) then
#endif
                  if (icgns .ne. 1) then
                     read(2) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
                  end if
#if defined DIST_MPI
               end if
#endif
            end if
            if ((ivisc(1).ge.5.or.ivisc(2).ge.5.or.ivisc(3).ge.5)
     .           .and. isminc.ne.2 .and. (ivmxold.ne.4 .and.
     .           ivmxold.ne.25)) isminc=3
            if ((ivisc(1).eq.25.or.ivisc(2).eq.25.or.ivisc(3).eq.25)
     .           .and. (ivmxold.ne.4 .and. ivmxold.ne.25)) isminc=1
            if ((ivisc(1).eq.25.or.ivisc(2).eq.25.or.ivisc(3).eq.25)
     .           .and. ivmxold.eq.4) isminc=3
            if ((ivisc(1).eq. 4.or.ivisc(2).eq. 4.or.ivisc(3).eq. 4)
     .           .and. ivmxold.eq.25) isminc=3
         end if
         if (isminc .eq. 2) then
            isminc=1
         end if
         if (ismincforce .ne. -1) isminc=ismincforce
      else
c     
         if (isminc .ne. 2) then
            isminc=0
         end if
         if (ismincforce .ne. -1) isminc=ismincforce
         if (iv1.ge.2 .or. iv2.ge.2 .or. iv3.ge.2) then
#if defined DIST_MPI
            if (myid.eq.mblk2nd(nbl)) then
#endif
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'('' reading vist3d data from '',
     .              '' restart file, block'',i5)') nbl
#if defined DIST_MPI
            end if
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) (((vist3d(j,k,i),j=1,jdim1),k=1,kdim1),
     $                 i=1,idim1)
                  if (irghost .ne. 0)
     .                 read(2) ((((vi0(j,k,l,m),j=1,jdim),k=1,kdim),
     $                 l=1,1),m=1,4),
     .                 ((((vj0(k,i,l,m),k=1,kdim),i=1,idim1),l=1,1),
     $                 m=1,4),
     .                 ((((vk0(j,i,l,m),j=1,jdim),i=1,idim1),l=1,1),
     $                 m=1,4)
               end if
#if defined DIST_MPI
            end if
c     
            mytag = itag_qv + nbl
            if (myid .eq. myhost) then
               call MPI_Send(vist3d,jki,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(vist3d,jki,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
            if (irghost .ne.0) then
               mytag = itag_vj0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(vj0,ki4,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(vj0,ki4,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
c     
               mytag = itag_vk0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(vk0,ji4,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(vk0,ji4,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
c     
               mytag = itag_vi0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(vi0,jk4,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(vi0,jk4,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
            end if
c     
#endif
         end if
         if (iv1.ge.4 .or. iv2.ge.4 .or. iv3.ge.4) then
#if defined DIST_MPI
            if (myid.eq.mblk2nd(nbl)) then
#endif
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'('' reading field eqn turb '',
     .              '' quantities from restart file'',
     .              '', block'',i6)') nbl
#if defined DIST_MPI
            end if
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) ((((tursav(j,k,i,m),j=1,jdim1),k=1,kdim1),
     $                 i=1,idim1),
     .                 m=1,2)
                  if (iv1.ge.30.or.iv2.ge.30.or.iv3.ge.30) then
                  read(2) nummem_read
                  read(2) ((((tursav(j,k,i,m),j=1,jdim1),k=1,kdim1),
     $                 i=1,idim1),
     .                 m=3,nummem_read)
                  end if
                  if (irghost .ne. 0) then
                       read(2) ((((ti0(j,k,l,m),j=1,jdim),k=1,kdim),
     $                 l=1,2),m=1,4),
     .                 ((((tj0(k,i,l,m),k=1,kdim),i=1,idim1),l=1,2),
     $                 m=1,4),
     .                 ((((tk0(j,i,l,m),j=1,jdim),i=1,idim1),l=1,2),
     $                 m=1,4)
                  if (iv1.ge.30.or.iv2.ge.30.or.iv3.ge.30) then
                       read(2) ((((ti0(j,k,l,m),j=1,jdim),k=1,kdim),
     $                 l=3,nummem_read),m=1,4),
     .                 ((((tj0(k,i,l,m),k=1,kdim),i=1,idim1),
     .                 l=3,nummem_read),m=1,4),
     .                 ((((tk0(j,i,l,m),j=1,jdim),i=1,idim1),
     .                 l=3,nummem_read),m=1,4)
                  end if
                  end if
               end if
               if ((iv1.eq.6.or.iv2.eq.6.or.iv3.eq.6.or.iv1.eq.7.or.
     .              iv2.eq.7.or.iv3.eq.7) .and. (ivisc(1).eq.8.or.
     .              ivisc(2).eq.8.or.ivisc(3).eq.8.or.ivisc(1).eq.12
     .              .or.ivisc(2).eq.12.or.ivisc(3).eq.12.or.
     .              ivisc(1).eq.14.or.ivisc(2).eq.14.or.
     $              ivisc(3).eq.14)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)*0.09
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.6.or.ivisc(2).eq.6.or.ivisc(3).eq.6
     .              .or.ivisc(1).eq.7.or.ivisc(2).eq.7.or.
     $              ivisc(3).eq.7) .and.
     .              (iv1.eq.8.or.iv2.eq.8.or.iv3.eq.8.or.
     .              iv1.eq.12.or.iv2.eq.12.or.iv3.eq.12.or.
     .              iv1.eq.14.or.iv2.eq.14.or.iv3.eq.14)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)/0.09
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.8.or.ivisc(2).eq.8.or.ivisc(3).eq.8
     .              .or.ivisc(1).eq.12.or.ivisc(2).eq.12.or.ivisc(3)
     $              .eq.12.or.ivisc(1).eq.14.or.ivisc(2).eq.14.or.
     $              ivisc(3).eq.14)
     .              .and.(iv1.eq.9.or.iv2.eq.9.or.iv3.eq.9.or.
     .              iv1.eq.10.or.iv2.eq.10.or.iv3.eq.10.or.
     .              iv1.eq.11.or.iv2.eq.11.or.iv3.eq.11.or.
     .              iv1.eq.13.or.iv2.eq.13.or.iv3.eq.13.or.
     .              iv1.eq.15.or.iv2.eq.15.or.iv3.eq.15)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)/
     $                          tursav(j,k,i,2)
c     Limit new omega to handle extremes
                           tursav(j,k,i,1)=
     $                          ccmincr(tursav(j,k,i,1),1.)
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.6.or.ivisc(2).eq.6.or.ivisc(3).eq.6
     .              .or.ivisc(1).eq.7.or.ivisc(2).eq.7.or.
     $              ivisc(3).eq.7) .and.
     .              (iv1.eq.9.or.iv2.eq.9.or.iv3.eq.9.or.
     .              iv1.eq.10.or.iv2.eq.10.or.iv3.eq.10.or.
     .              iv1.eq.11.or.iv2.eq.11.or.iv3.eq.11.or.
     .              iv1.eq.13.or.iv2.eq.13.or.iv3.eq.13.or.
     .              iv1.eq.15.or.iv2.eq.15.or.iv3.eq.15)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)/
     +                          (0.09*tursav(j,k,i,2))
c     Limit new omega to handle extremes
                           tursav(j,k,i,1)=
     $                          ccmincr(tursav(j,k,i,1),1.)
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.9.or.ivisc(2).eq.9.or.ivisc(3).eq.9
     .              .or.ivisc(1).eq.10.or.ivisc(2).eq.10
     .              .or.ivisc(3).eq.10
     .              .or.ivisc(1).eq.11.or.ivisc(2).eq.11.or.
     $              ivisc(3).eq.11.or.
     .              ivisc(1).eq.13.or.ivisc(2).eq.13.or.
     $              ivisc(3).eq.13.or.
     .              ivisc(1).eq.15.or.ivisc(2).eq.15.or.
     $              ivisc(3).eq.15)
     .              .and.(iv1.eq.8.or.iv2.eq.8.or.iv3.eq.8.or.
     .              iv1.eq.12.or.iv2.eq.12.or.iv3.eq.12.or.
     .              iv1.eq.14.or.iv2.eq.14.or.iv3.eq.14)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)*
     $                          tursav(j,k,i,2)
c     Limit new epsilon to handle extremes
                           tursav(j,k,i,1)=
     $                          ccmincr(tursav(j,k,i,1),1.)
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.9.or.ivisc(2).eq.9.or.
     $              ivisc(3).eq.9.or.
     .              ivisc(1).eq.10.or.ivisc(2).eq.10.or.
     $              ivisc(3).eq.10.or.
     .              ivisc(1).eq.11.or.ivisc(2).eq.11.or.
     $              ivisc(3).eq.11.or.
     .              ivisc(1).eq.13.or.ivisc(2).eq.13.or.
     $              ivisc(3).eq.13.or.
     .              ivisc(1).eq.15.or.ivisc(2).eq.15.or.
     $              ivisc(3).eq.15)
     .              .and.(iv1.eq.6.or.iv2.eq.6.or.iv3.eq.6.or.
     .              iv1.eq.7.or.iv2.eq.7.or.iv3.eq.7)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,1)=tursav(j,k,i,1)*
     +                          0.09*tursav(j,k,i,2)
c     Limit new epsilon to handle extremes
                           tursav(j,k,i,1)=
     $                          ccmincr(tursav(j,k,i,1),1.)
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.72.or.ivisc(2).eq.72.or.
     $              ivisc(3).eq.72)
     .              .and.(iv1.eq.6.or.iv2.eq.6.or.iv3.eq.6.or.
     .              iv1.eq.7.or.iv2.eq.7.or.iv3.eq.7)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,7)=tursav(j,k,i,1)
                           tursav(j,k,i,1)=-tursav(j,k,i,2)*2./3.
                        enddo
                     enddo
                  enddo
               end if
               if ((ivisc(1).eq.72.or.ivisc(2).eq.72.or.
     $              ivisc(3).eq.72)
     .              .and.(iv1.eq.8.or.iv2.eq.8.or.iv3.eq.8.or.
     .              iv1.eq.12.or.iv2.eq.12.or.iv3.eq.12.or.
     .              iv1.eq.14.or.iv2.eq.14.or.iv3.eq.14)) then
                  do i=1,idim1
                     do k=1,kdim1
                        do j=1,jdim1
                           tursav(j,k,i,7)=tursav(j,k,i,1)/0.09
                           tursav(j,k,i,1)=-tursav(j,k,i,2)*2./3.
                        enddo
                     enddo
                  enddo
               end if
#if defined DIST_MPI
c     
c     end host section for tursav
c     
            end if
c     
            if (iv1.eq.30.or.iv2.eq.30.or.iv3.eq.30) then
              jki2  = 3*jki
              jk8   = jdim*kdim*12
              ki8   = kdim*idim1*12
              ji8   = jdim*idim1*12
            end if
            if (iv1.eq.40.or.iv2.eq.40.or.iv3.eq.40) then
              jki2  = 4*jki
              jk8   = jdim*kdim*16
              ki8   = kdim*idim1*16
              ji8   = jdim*idim1*16
            end if
            if (iv1.eq.72.or.iv2.eq.72.or.iv3.eq.72) then
              jki2  = 7*jki
              jk8   = jdim*kdim*28
              ki8   = kdim*idim1*28
              ji8   = jdim*idim1*28
            end if
            mytag = itag_qt + nbl
            if (myid .eq. myhost) then
               call MPI_Send(tursav,jki2,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(tursav,jki2,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
            if (irghost .ne.0) then
               mytag = itag_tj0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(tj0,ki8,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(tj0,ki8,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
c     
               mytag = itag_tk0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(tk0,ji8,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(tk0,ji8,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
c     
               mytag = itag_ti0 + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(ti0,jk8,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(ti0,jk8,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
            end if
c     
            if (myid.eq.myhost) then
#endif
               if (icgns .ne. 1) then
                  read(2) (((smin(j,k,i),j=1,jdim1),
     $                 k=1,kdim1),i=1,idim1)
               end if
#if defined DIST_MPI
            end if
c     
            mytag = itag_smin + nbl
            if (myid .eq. myhost) then
               call MPI_Send(smin,jkim,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(smin,jkim,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
#endif
            if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4 .or.
     .          iv1.eq.25.or. iv2.eq.25.or. iv3.eq.25) then
#if defined DIST_MPI
               if (myid.eq.myhost) then
#endif
                  if (icgns .ne. 1) then
                     read(2) (((xjb(j,k,i),j=1,jdim1),
     $                    k=1,kdim1),i=1,idim1)
                     read(2) (((xkb(j,k,i),j=1,jdim1),
     $                    k=1,kdim1),i=1,idim1)
                     read(2) (((blnum(j,k,i),j=1,jdim1),
     $                    k=1,kdim1),i=1,idim1)
                  end if
#if defined DIST_MPI
               end if
c     
               mytag = itag_xjb + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(xjb,jkim,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(xjb,jkim,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
               mytag = itag_xkb + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(xkb,jkim,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(xkb,jkim,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
               mytag = itag_blnum + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(blnum,jkim,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(blnum,jkim,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
#endif
            end if
            if (iv1.eq.8 .or. iv2.eq.8 .or. iv3.eq.8  .or.
     .           iv1.eq.9 .or. iv2.eq.9 .or. iv3.eq.9  .or.
     .           iv1.eq.13.or. iv2.eq.13.or. iv3.eq.13 .or.
     .           iv1.eq.14.or. iv2.eq.14.or. iv3.eq.14) then
#if defined DIST_MPI
               if (myid.eq.myhost) then
#endif
                  if (icgns .ne. 1) then
                     if(ivmx .eq. 8 .or. ivmx .eq. 9
     $                    .or. ivmx .eq. 13 .or. ivmx .eq. 14) then
                        read(2) (((cmuv(j,k,i),j=1,jdim1),k=1,kdim1),
     $                       i=1,idim1)
                     else
                        read(2) (((dum,j=1,jdim1),k=1,kdim1),
     $                       i=1,idim1)
                     end if
                  end if
#if defined DIST_MPI
               end if
c     
               mytag = itag_cmuv + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(cmuv,jkim,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. mblk2nd(nbl)) then
                  call MPI_Recv(cmuv,jkim,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
               end if
#endif
            end if
         end if
      end if
c     
c     If irest2=1 (results when irest<0 on input), do not save
c     previous histories; otherwise, do.
c     
      if (irest2 .eq. 1) then
         ntt = 0
         ntr = 0
      else
c     
c     this was defined earlier for message passing
c     
c     ntt = ntr
c     
      end if

      if (iteravg .eq. 2) then
         if (ipertavg .eq. 0) then
            if (myid.eq.myhost) then
               read(97,end=1010,err=1010) dum,dum,dum,xnumavg
               read(97,end=1010,err=1010) ((((qavg(j,k,i,m),
     +              i=1,idim-1),j=1,jdim-1),k=1,kdim-1),m=1,5)
c     need to switch y and z if ialph .ne. 0
               if (ialph .ne. 0) then
                  do i=1,idim-1
                     do j=1,jdim-1
                        do k=1,kdim-1
                           temp=qavg(j,k,i,4)
                           qavg(j,k,i,4)=qavg(j,k,i,3)
                           qavg(j,k,i,3)=-temp
                        enddo
                     enddo
                  enddo
               end if
            end if
#if defined DIST_MPI
c     
            mytag = itag_qavg + nbl
            if (myid .eq. myhost) then
               call MPI_Send(qavg,jki5,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(qavg,jki5,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
            mytag = itag_xnum + nbl
            if (myid .eq. myhost) then
               call MPI_Send(xnumavg,1,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(xnumavg,1,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
#endif
         else
            if (myid.eq.myhost) then
               read(97,end=1010,err=1010) dum,dum,dum,xnumavg
               read(97,end=1010,err=1010) ((((qavg(j,k,i,m),
     +              i=1,idim),j=1,jdim),k=1,kdim),m=1,5) 
c     
               read(98,end=1010,err=1010) dum,dum,dum,xnumavg2
               read(98,end=1010,err=1010) ((((q2avg(j,k,i,m), 
     +              i=1,idim),j=1,jdim),k=1,kdim),m=1,5)
            
c     need to switch y and z if ialph .ne. 0
               if (ialph .ne. 0) then
                  do i=1,idim
                     do j=1,jdim
                        do k=1,kdim
c     
                           temp=qavg(j,k,i,4)
                           qavg(j,k,i,4)=qavg(j,k,i,3)
                           qavg(j,k,i,3)=-temp
c     
                           temp=q2avg(j,k,i,4) 
                           q2avg(j,k,i,4)=q2avg(j,k,i,3) 
                           q2avg(j,k,i,3)=temp 
c     
                        enddo
                     enddo
                  enddo
               end if
            end if
c     
#if defined DIST_MPI
c     
            mytag = itag_qavg + nbl
            if (myid .eq. myhost) then
               call MPI_Send(qavg,jki5,MY_MPI_REAL, 
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(qavg,jki5,MY_MPI_REAL, 
     .              myhost,mytag,mycomm,istat,ierr)
            end if
            mytag = itag_xnum + nbl
            if (myid .eq. myhost) then
               call MPI_Send(xnumavg,1,MY_MPI_REAL,
     .              nd_dest,mytag,mycomm,ierr)
            else if (myid .eq. mblk2nd(nbl)) then
               call MPI_Recv(xnumavg,1,MY_MPI_REAL,
     .              myhost,mytag,mycomm,istat,ierr)
            end if
c     
            mytag = itag_q2avg + nbl 
            if (myid .eq. myhost) then 
               call MPI_Send(q2avg,jki5,MY_MPI_REAL, 
     .              nd_dest,mytag,mycomm,ierr) 
            else if (myid .eq. mblk2nd(nbl)) then 
               call MPI_Recv(q2avg,jki5,MY_MPI_REAL, 
     .              myhost,mytag,mycomm,istat,ierr) 
            end if                 
            mytag = itag_xnum2 + nbl 
            if (myid .eq. myhost) then 
               call MPI_Send(xnumavg2,1,MY_MPI_REAL, 
     .              nd_dest,mytag,mycomm,ierr) 
            else if (myid .eq. mblk2nd(nbl)) then 
               call MPI_Recv(xnumavg2,1,MY_MPI_REAL, 
     .              myhost,mytag,mycomm,istat,ierr) 
            end if                 
c     
#endif
         end if ! ipertavg
      end if ! iteravg
c
c     clcd information
c
              
      if (iskip.gt.0) then
         if (myid.eq.myhost) then
            if( iclcd .eq. 1 .or. iclcd .eq. 2 ) then
               clcd(1,:,:) = 1.e21
               clcd(2,:,:) = 0.e0
               i_ntt = 0
            end if
            if (iclcd .eq. 2) then
               read(102,end=1012,err=1012) i_clcd,iblocks_clcd,i_ntt
c     
               if (n_clcd .ne. i_clcd ) then
                  write(*,*) 'Inconsistent number of ClCd ',
     $                 'calcs in restart file ', n_clcd, i_clcd
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
            
               read(102,end=1012,err=1012) clcd(1:2,1:i_clcd,1:i_ntt)
            end if
         end if                 
      end if

#if defined DIST_MPI
      mytag = itag_clcd + nbl 
      if (myid.eq.myhost) then
         call MPI_Send(clcd,ntt*2*n_clcd,MY_MPI_REAL,
     .        nd_dest,mytag,mycomm,ierr) 
      else if (myid .eq. mblk2nd(nbl)) then 
         call MPI_Recv(clcd,ntt*2*n_clcd,MY_MPI_REAL, 
     .        myhost,mytag,mycomm,istat,ierr) 
      end if
#endif

      return
 1010 continue
      write(11,'(/,'' stopping... flag iteravg=2, but cannot read'',
     +     '' running-average Q file, nbl='',i6)') nbl
      write(11,'('' ... it either does not exist or it is an'',
     +     '' incorrect file'')')
      call termn8(myid,-1,ibufdim,nbuf,bou,nou)
 1012 continue
      write(11,'(/,'' stopping... flag iclcd=2, but cannot read'',
     +     '' clcd.bin file, nbl='',i6)') nbl
      write(11,'('' ... it either does not exist or it is an'',
     +     '' incorrect file'')')
      call termn8(myid,-1,ibufdim,nbuf,bou,nou)

c     
      return
      end
