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
      program cfl3d_to_pegbc
c
c     $Id$
c
c***********************************************************************
c     Purpose: Creates a peg.bc.raw file for use with PEGSUS 5.x
c     Note: this program does NOT use list-directed input. Rather,
c     command-line arguments are used to determine which cfl3d input
c     file to process. Usage:
c
c        cfl3d_to_pegbc [cfl3d input file name]
c
c        if no file name specified, cfl3d.inp is used by default
c***********************************************************************
c
      parameter (ibufdim=2000,nbuf=4,mxbcfil=100)
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,inpfile
      character*80 inpstring
      character*120 bou(ibufdim,nbuf)
c
      dimension itrnsfr(27),titlw(20),nou(nbuf)
c
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /unit5/ iunit5
      common /key/ nkey
      common /zero/ iexp
c
c     set some parameters not really need by this program but used by
c     some subroutines borrowed in from other cfl3d modules
c
      myid   = 0
      myhost = 0
      nnodes = 1
      mycomm = 1
      xmach  = 0.1
      iexp   = -7
c
      do ii=1,nbuf
         nou(ii) = 0
         do mm=1,ibufdim
            bou(mm,ii) = ' '
         end do
      end do 
c
c***********************************************************************
c     open files
c***********************************************************************
c
c     get the input file name from the command line
c
      call get_cmd_args(inpfile,ierr)
c
      if (ierr .gt. 0) then
         write(6,*)
         write(6,'(''usage: cfl3d_to_pegbc [cfl3d input file name]'')')
         write(6,'(''if no file name specified, cfl3d.inp is used'')')
         stop
      end if
c
      iunit5 = 32
      open(iunit5,file=inpfile,form='formatted',status='unknown')
c
      rewind(iunit5)
c
      read(iunit5,*)
      read(iunit5,'(a60)')grid
      read(iunit5,'(a60)')plt3dg
      read(iunit5,'(a60)')plt3dq
      read(iunit5,'(a60)')output
      read(iunit5,'(a60)')residual
      read(iunit5,'(a60)')turbres
      read(iunit5,'(a60)')blomx
      read(iunit5,'(a60)')output2
      read(iunit5,'(a60)')printout
      read(iunit5,'(a60)')pplunge
      read(iunit5,'(a60)')ovrlap
      read(iunit5,'(a60)')patch
      read(iunit5,'(a60)')restrt
c
      open(unit=2,file=restrt,form='unformatted',
     .     status='unknown')
      open(unit=7,file='peg.bc.raw',form='formatted',
     .     status='unknown')
      open(unit=11,file='cfl3d_to_pegbc.out',form='formatted',
     .     status='unknown')
      open(unit=22,file=patch,form='unformatted',
     .     status='unknown')
      open(unit=99,file='cfl3d_to_pegbc.error',form='formatted',
     .     status='unknown')
c
      rewind(2)
      rewind(7)
      rewind(11)
      rewind(66)
      rewind(99)
c
c***********************************************************************
c     determine array size requirements
c***********************************************************************
c
c     read input file to get the array dimensions needed by precfl3d
c
      rewind(iunit5)
c
      ibufdim0 = ibufdim
      nbuf0    = nbuf
      mxbcfil0 = mxbcfil
c
c     global0 sets the parameters needed for precfl3d (sizer)
c
      iunit11 = 99
      call global0(nplots0,maxnode0,mxbli0,lbcprd0,lbcemb0,
     .             lbcrad0,maxbl0,maxgr0,maxseg0,maxcs0,ncycmax0,
     .             intmax0,nsub10,intmx0,mxxe0,mptch0,msub10,
     .             ibufdim0,nbuf0,mxbcfil0,nmds0,maxaes0,
     .             maxsegdg0,ntr,nnodes,nou,bou,iunit11,myid,
     .             idm0,jdm0,kdm0)
c
c     rewind the input file so it can be read again
c
      rewind(iunit5)
c
c***********************************************************************
c     write peg.bc file
c***********************************************************************
c
      iunit11 = 11
      call pegbcout(nplots0,maxnode0,mxbli0,lbcprd0,lbcemb0,
     .              lbcrad0,maxbl0,maxgr0,maxseg0,maxcs0,ncycmax0,
     .              intmax0,nsub10,intmx0,mxxe0,mptch0,msub10,
     .              ibufdim0,nbuf0,mxbcfil0,nmds0,maxaes0,
     .              maxsegdg0,ntr,nnodes,nou,bou,iunit11)
c
c***********************************************************************
c     normal program termnination
c***********************************************************************
c
      call termn8(myid,0,ibufdim,nbuf,bou,nou)
c
      stop
      end
c
      subroutine pegbcout(nplots,maxnode,mxbli,lbcprd,lbcemb,
     .                 lbcrad,maxbl,maxgr,maxseg,maxcs,ncycmax,
     .                 intmax,nsub1,intmx,mxxe,mptch,msub1,
     .                 ibufdim,nbuf,mxbcfil,nmds,maxaes,
     .                 maxsegdg,ntr,nnodes,nou,bou,iunit11)
c
      character*120 bou(ibufdim,nbuf)
      character*40 string
c
      integer stats
c
      dimension nou(nbuf)
c
      allocatable :: ibcinfo(:,:,:,:)
      allocatable :: idir(:)
      allocatable :: igridg(:)
      allocatable :: iovrlp(:)
      allocatable :: isva(:,:,:)
      allocatable :: isym(:)
      allocatable :: itype(:)
      allocatable :: jbce(:)
      allocatable :: jbcinfo(:,:,:,:)
      allocatable :: jbcs(:)
      allocatable :: jdimg(:)
      allocatable :: jxie_e(:)
      allocatable :: jxie_s(:)
      allocatable :: kbce(:)
      allocatable :: kbcinfo(:,:,:,:)
      allocatable :: kbcs(:)
      allocatable :: kdimg(:)
      allocatable :: keta_e(:)
      allocatable :: keta_s(:)
      allocatable :: lbce(:)
      allocatable :: lbcs(:)
      allocatable :: ldimg(:)
      allocatable :: levelg(:)
      allocatable :: limblk(:,:,:)
      allocatable :: nbci0(:)
      allocatable :: nbcidim(:)
      allocatable :: nbcj0(:)
      allocatable :: nbcjdim(:)
      allocatable :: nbck0(:)
      allocatable :: nbckdim(:)
      allocatable :: nbl_to(:)
      allocatable :: nblk(:,:)
      allocatable :: nblon(:)
      allocatable :: nm_to(:)
c
      common /unit5/ iunit5
c
c     allocate memory
c
      memuse = 0
      allocate( ibcinfo(maxbl,maxseg,7,2), stat=stats )
      call umalloc_r(maxbl*maxseg*7*2,1,'ibcinfo',memuse,stats)
      allocate( idir(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'idir',memuse,stats)
      allocate( igridg(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'igridg',memuse,stats)
      allocate( iovrlp(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'iovrlp',memuse,stats)
      allocate( isva(2,2,mxbli), stat=stats )
      call umalloc_r(2*2*mxbli,1,'isva',memuse,stats)
      allocate( isym(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'isym',memuse,stats)
      allocate( itype(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'itype',memuse,stats)
      allocate( jbce(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'jbce',memuse,stats)
      allocate( jbcinfo(maxbl,maxseg,7,2), stat=stats )
      call umalloc_r(maxbl*maxseg*7*2,1,'jbcinfo',memuse,stats)
      allocate( jbcs(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'jbcs',memuse,stats)
      allocate( jdimg(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'jdimg',memuse,stats)
      allocate( jxie_e(intmax), stat=stats )
      call umalloc_r(intmax,1,'jxie_e',memuse,stats)
      allocate( jxie_s(intmax), stat=stats )
      call umalloc_r(intmax,1,'jxie_s',memuse,stats)
      allocate( kbce(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'kbce',memuse,stats)
      allocate( kbcinfo(maxbl,maxseg,7,2), stat=stats )
      call umalloc_r(maxbl*maxseg*7*2,1,'kbcinfo',memuse,stats)
      allocate( kbcs(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'kbcs',memuse,stats)
      allocate( kdimg(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'kdimg',memuse,stats)
      allocate( keta_e(intmax), stat=stats )
      call umalloc_r(intmax,1,'keta_e',memuse,stats)
      allocate( keta_s(intmax), stat=stats )
      call umalloc_r(intmax,1,'keta_s',memuse,stats)
      allocate( lbce(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'lbce',memuse,stats)
      allocate( lbcs(maxseg*6+2*mxbli+intmax), stat=stats )
      call umalloc_r((maxseg*6+2*mxbli+intmax),1,'lbcs',memuse,stats)
      allocate( ldimg(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'ldimg',memuse,stats)
      allocate( levelg(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'levelg',memuse,stats)
      allocate( limblk(2,6,mxbli), stat=stats )
      call umalloc_r(2*6*mxbli,1,'limblk',memuse,stats)
      allocate( nbci0(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbci0',memuse,stats)
      allocate( nbcidim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbcidim',memuse,stats)
      allocate( nbcj0(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbcj0',memuse,stats)
      allocate( nbcjdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbcjdim',memuse,stats)
      allocate( nbck0(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbck0',memuse,stats)
      allocate( nbckdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'nbckdim',memuse,stats)
      allocate( nbl_to(intmax), stat=stats )
      call umalloc_r(intmax,1,'nbl_to',memuse,stats)
      allocate( nblk(2,mxbli), stat=stats )
      call umalloc_r(2*mxbli,1,'nblk',memuse,stats)
      allocate( nblon(mxbli), stat=stats )
      call umalloc_r(mxbli,1,'nblon',memuse,stats)
      allocate( nm_to(intmax), stat=stats )
      call umalloc_r(intmax,1,'nm_to',memuse,stats)
c
      do nn=1,nbuf
         nou(nn) = 0
      end do
c
c     read i/o file names
c
      nread = 14
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c
c     read keyword-driven input, if any
c
      call readkey(ititr,myid,ibufdim,nbuf,bou,nou,0,-99)
c
c     read title
c
      if (ititr.eq.0) then
         call echoinp(iunit5,iunit11,1)
      end if
c
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) dum1,dum1,dum1,dum1,dum1,ialph,idum
      write(iunit11,'(3f10.5,e10.3,f10.5,i10,2i10)') dum1,dum1,dum1,
     .      dum1,dum1,ialph,idum
      call echoinp(iunit5,iunit11,0)
      nread = 3
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c
      read(iunit5,*) dt,irest,iflagts,fmax,iunst,cfltau
      write(iunit11,'(f10.5,2i10,f10.5,i10,f10.5)') dt,irest,
     .     iflagts,fmax,iunst,cfltau
      call echoinp(iunit5,iunit11,0)
      if (real(dt).lt.0.) iunst = 0
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) ngrid,nplot3d,nprint,nwrest,ichkd,i2d,ntstep,ita
      write(iunit11,'(8i10)') ngrid,nplot3d,nprint,nwrest,ichkd,
     .     i2d,ntstep,ita
      call echoinp(iunit5,iunit11,0)
c
      ngrid = abs(ngrid)
c
      call echoinp(iunit5,iunit11,1)
c
c     read ncg data
      ncgmax = 0
      nbl    = 0
      do n=1,ngrid
         read(iunit5,*) ncg,idum1,idum2,idum3,idum4,idum5,idum6
         write(iunit11,'(7i10)') ncg,idum1,idum2,idum3,idum4,idum5,idum6
         call echoinp(iunit5,iunit11,0)
         ncgmax = max(ncg,ncgmax)
         nbl = nbl + 1
         igridg(nbl) = n
         levelg(nbl) = 1
         if (ncg .gt. 0) then
            do nn=1,ncg
               nbl = nbl + 1
               igridg(nbl) = n
               levelg(nbl) = nn + 1
            end do
         end if
      end do
c
      nread = ngrid
c
c     read grid dims section
      call echoinp(iunit5,iunit11,1)
      do n=1,nread
         read(iunit5,*) ldimg(n),jdimg(n),kdimg(n)
         write(iunit11,'(3i10)') ldimg(n),jdimg(n),kdimg(n)
         call echoinp(iunit5,iunit11,0)
      end do
c
      nread = ngrid + 1
c     read laminar regions section
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c     read embeded grid section
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c     read idiag/iflim section
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c     read ifds/rkap0 section
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c
c     read no. of bc segments section
c
      call echoinp(iunit5,iunit11,1)
      nread = ngrid
      do n=1,nread
         read(iunit5,*) mdum,nbci0(n),nbcidim(n),nbcj0(n),nbcjdim(n),
     .                  nbck0(n),nbckdim(n)
         write(iunit11,'(8i10)') mdum,nbci0(n),nbcidim(n),nbcj0(n),
     .        nbcjdim(n),nbck0(n),nbckdim(n)
         call echoinp(iunit5,iunit11,0)
      end do
c
c     read through bc section
c
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbci0(n)
            read(iunit5,*) ig,nseg,ibcinfo(n,nn,1,1),
     .                     ibcinfo(n,nn,2,1),
     .                     ibcinfo(n,nn,3,1),
     .                     ibcinfo(n,nn,4,1),
     .                     ibcinfo(n,nn,5,1),ndata
            write(iunit11,'(8i10)') ig,nseg,ibcinfo(n,nn,1,1),
     .                     ibcinfo(n,nn,2,1),
     .                     ibcinfo(n,nn,3,1),
     .                     ibcinfo(n,nn,4,1),
     .                     ibcinfo(n,nn,5,1),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbcidim(n)
            read(iunit5,*) ig,nseg,ibcinfo(n,nn,1,2),
     .                     ibcinfo(n,nn,2,2),
     .                     ibcinfo(n,nn,3,2),
     .                     ibcinfo(n,nn,4,2),
     .                     ibcinfo(n,nn,5,2),ndata
            write(iunit11,'(8i10)') ig,nseg,ibcinfo(n,nn,1,2),
     .                     ibcinfo(n,nn,2,2),
     .                     ibcinfo(n,nn,3,2),
     .                     ibcinfo(n,nn,4,2),
     .                     ibcinfo(n,nn,5,2),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbcj0(n)
            read(iunit5,*) ig,nseg,jbcinfo(n,nn,1,1),
     .                     jbcinfo(n,nn,2,1),
     .                     jbcinfo(n,nn,3,1),
     .                     jbcinfo(n,nn,4,1),
     .                     jbcinfo(n,nn,5,1),ndata
            write(iunit11,'(8i10)') ig,nseg,jbcinfo(n,nn,1,1),
     .                     jbcinfo(n,nn,2,1),
     .                     jbcinfo(n,nn,3,1),
     .                     jbcinfo(n,nn,4,1),
     .                     jbcinfo(n,nn,5,1),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbcjdim(n)
            read(iunit5,*) ig,nseg,jbcinfo(n,nn,1,2),
     .                     jbcinfo(n,nn,2,2),
     .                     jbcinfo(n,nn,3,2),
     .                     jbcinfo(n,nn,4,2),
     .                     jbcinfo(n,nn,5,2),ndata
            write(iunit11,'(8i10)') ig,nseg,jbcinfo(n,nn,1,2),
     .                     jbcinfo(n,nn,2,2),
     .                     jbcinfo(n,nn,3,2),
     .                     jbcinfo(n,nn,4,2),
     .                     jbcinfo(n,nn,5,2),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbck0(n)
            read(iunit5,*) ig,nseg,kbcinfo(n,nn,1,1),
     .                     kbcinfo(n,nn,2,1),
     .                     kbcinfo(n,nn,3,1),
     .                     kbcinfo(n,nn,4,1),
     .                     kbcinfo(n,nn,5,1),ndata
            write(iunit11,'(8i10)') ig,nseg,kbcinfo(n,nn,1,1),
     .                     kbcinfo(n,nn,2,1),
     .                     kbcinfo(n,nn,3,1),
     .                     kbcinfo(n,nn,4,1),
     .                     kbcinfo(n,nn,5,1),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,ngrid
         do nn=1,nbckdim(n)
            read(iunit5,*) ig,nseg,kbcinfo(n,nn,1,2),
     .                     kbcinfo(n,nn,2,2),
     .                     kbcinfo(n,nn,3,2),
     .                     kbcinfo(n,nn,4,2),
     .                     kbcinfo(n,nn,5,2),ndata
            write(iunit11,'(8i10)') ig,nseg,kbcinfo(n,nn,1,2),
     .                     kbcinfo(n,nn,2,2),
     .                     kbcinfo(n,nn,3,2),
     .                     kbcinfo(n,nn,4,2),
     .                     kbcinfo(n,nn,5,2),ndata
            call echoinp(iunit5,iunit11,0)
            if (abs(ndata) .gt.0) then
               call echoinp(iunit5,iunit11,1)
               call echoinp(iunit5,iunit11,1)
            end if
         end do
      end do
c
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) mseq,idum1,idum2,idum3,idum4
      write(iunit11,'(5i10)') mseq,idum1,idum2,idum3,idum4
      call echoinp(iunit5,iunit11,0)
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      ncyctot = 0 
      do n=1,mseq
         read(iunit5,*) ncyc,idum1,idum2,idum3
         write(iunit11,'(4i10)') ncyc,idum1,idum2,idum3
         call echoinp(iunit5,iunit11,0)
         if (real(dt).lt.0.e0) then
            ncyctot = ncyctot+ncyc
         else
            if (ncyc.gt.0) then
               ncyctot = ncyctot+ntstep
            end if
         end if
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,mseq
         call echoinp(iunit5,iunit11,1)
      end do
c
c     read through 1-1 interface data
c
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) nbli0
      write(iunit11,'(i10)') nbli0
      call echoinp(iunit5,iunit11,0)
      call echoinp(iunit5,iunit11,1)
      do n=1,abs(nbli0)
         read(iunit5,*) mdum,nblk(1,n),
     .     (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
         write(iunit11,'(10i8)') mdum,nblk(1,n),
     .     (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
         call echoinp(iunit5,iunit11,0)
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,abs(nbli0)
         read(iunit5,*) mdum,nblk(2,n),
     .     (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
         write(iunit11,'(10i8)') mdum,nblk(2,n),
     .     (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
         call echoinp(iunit5,iunit11,0)
      end do
c
c     read through patch data file
c
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) nint0
      write(iunit11,'(i10)') nint0
      call echoinp(iunit5,iunit11,0)
      if (nint0.eq.0) then
         intmax0 = 1
         nsub10  = 1
         ninter  = 0
      else 
         read(22) ninter
         do n=1,abs(ninter)
            read(22) nfb
            read(22) (idum,ll=1,nfb)
            read(22) nbl_to(n)
            read(22) (idum,ll=1,nfb)
            read(22) nm_to(n)
            read(22) idum2
            read(22) lst
            read(22) jxie_s(n)
            read(22) jxie_e(n)
            read(22) keta_s(n)
            read(22) keta_e(n)
            len = lst + idum2 - 1
            read(22) (ndum,nnn=lst,len)
            read(22) ((dum,nnn=lst,len),ll=1,2)
         end do
      end if
c
c     that should be all that's needed to generate a peg.bc file,
c     except for some mesh names
c
      write(6,'(''enter 0 to specify a name for each mesh'')')
      write(6,'(''enter 1 to  use default names (grid.n)'')')
      read(5,*) iname
c
      izeero = 0
      ione   = 1
c
      do n=1,ngrid
c
         nbc     = 0
         isym(n) = 0
c
c        convert cfl3d bc's to pegsus bc's. pegsus bc's are limited to:
c
c        1...all types not listed below (i.e not 2 or 3)
c        2...solid wall
c        3...blanked out region (no corresponding cfl3d bc)
c
c        any cfl3d zonal bc's that are NOT chimera bc's  (i.e. 1-1
c        and patched bc's) will get flagged as pegsus bc type 1
c        to prevent them from being set as an outer boundary in pegsus.
c
c        1st convert all physical bc's
c
         do nn=1,nbci0(n)
            ibctype = abs(ibcinfo(n,nn,1,1))
            if (ibctype .eq. 1001) isym(n) = 1
            if (ibctype .gt. 0) then
               nbc = nbc + 1
               if (ibctype.eq.1005 .or. abs(ibctype).eq.2004 .or.
     .             ibctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = 3
               jbcs(nbc) = ibcinfo(n,nn,2,1)
               jbce(nbc) = ibcinfo(n,nn,3,1)
               kbcs(nbc) = ibcinfo(n,nn,4,1)
               kbce(nbc) = ibcinfo(n,nn,5,1)
               lbcs(nbc) = 1
               lbce(nbc) = 1
               call shortinp(jbcs(nbc),jbce(nbc),kbcs(nbc),kbce(nbc))
            end if
         end do 
c
         do nn=1,nbcidim(n)
            ibctype = abs(ibcinfo(n,nn,1,2))
            if (ibctype .eq. 1001) isym(n) = 1
            if (ibctype .gt. 0) then
               nbc = nbc + 1
               if (ibctype.eq.1005 .or. abs(ibctype).eq.2004 .or.
     .             ibctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = -3
               jbcs(nbc) = ibcinfo(n,nn,2,2)
               jbce(nbc) = ibcinfo(n,nn,3,2)
               kbcs(nbc) = ibcinfo(n,nn,4,2)
               kbce(nbc) = ibcinfo(n,nn,5,2)
               lbcs(nbc) = -1
               lbce(nbc) = -1
               call shortinp(jbcs(nbc),jbce(nbc),kbcs(nbc),kbce(nbc))
            end if
         end do
c
         do nn=1,nbcj0(n)
            jbctype = abs(jbcinfo(n,nn,1,1))
            if (jbctype .eq. 1001) isym(n) = 1
            if (jbctype .gt. 0) then
               nbc = nbc + 1
               if (jbctype.eq.1005 .or. abs(jbctype).eq.2004 .or.
     .             jbctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = 1
               jbcs(nbc) = 1
               jbce(nbc) = 1
               kbcs(nbc) = jbcinfo(n,nn,4,1)
               kbce(nbc) = jbcinfo(n,nn,5,1)
               lbcs(nbc) = jbcinfo(n,nn,2,1)
               lbce(nbc) = jbcinfo(n,nn,3,1)
               call shortinp(lbcs(nbc),lbce(nbc),kbcs(nbc),kbce(nbc))
            end if
         end do
c
         do nn=1,nbcjdim(n)
            jbctype = abs(jbcinfo(n,nn,1,2))
            if (jbctype .eq. 1001) isym(n) = 1
            if (jbctype .gt. 0) then
               nbc = nbc + 1
               if (jbctype.eq.1005 .or. abs(jbctype).eq.2004 .or.
     .             jbctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = -1
               jbcs(nbc) = -1
               jbce(nbc) = -1
               kbcs(nbc) = jbcinfo(n,nn,4,2)
               kbce(nbc) = jbcinfo(n,nn,5,2)
               lbcs(nbc) = jbcinfo(n,nn,2,2)
               lbce(nbc) = jbcinfo(n,nn,3,2)
               call shortinp(lbcs(nbc),lbce(nbc),kbcs(nbc),kbce(nbc))
            end if
         end do
c
         do nn=1,nbck0(n)
            kbctype = abs(kbcinfo(n,nn,1,1))
            if (kbctype .eq. 1001) isym(n) = 1
            if (kbctype .gt. 0) then
               nbc = nbc + 1
               if (kbctype.eq.1005 .or. abs(kbctype).eq.2004 .or.
     .             kbctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = 2
               jbcs(nbc) = kbcinfo(n,nn,4,1)
               jbce(nbc) = kbcinfo(n,nn,5,1)
               kbcs(nbc) = 1
               kbce(nbc) = 1
               lbcs(nbc) = kbcinfo(n,nn,2,1)
               lbce(nbc) = kbcinfo(n,nn,3,1)
               call shortinp(lbcs(nbc),lbce(nbc),jbcs(nbc),jbce(nbc))
            end if
         end do
c
         do nn=1,nbckdim(n)
            kbctype = abs(kbcinfo(n,nn,1,2))
            if (kbctype .eq. 1001) isym(n) = 1
            if (kbctype .gt. 0) then
               nbc = nbc + 1
               if (kbctype.eq.1005 .or. abs(kbctype).eq.2004 .or.
     .             kbctype.eq.2014) then
                  itype(nbc) = 2
               else
                  itype(nbc) = 1
               end if
               idir(nbc) = -2
               jbcs(nbc) = kbcinfo(n,nn,4,2)
               jbce(nbc) = kbcinfo(n,nn,5,2)
               kbcs(nbc) = -1
               kbce(nbc) = -1
               lbcs(nbc) = kbcinfo(n,nn,2,2)
               lbce(nbc) = kbcinfo(n,nn,3,2)
               call shortinp(lbcs(nbc),lbce(nbc),jbcs(nbc),jbce(nbc))
            end if
         end do
c
c        next, convert all 1-1 data to peg bc type 1
c
         do int=1,nbli0
            do ib = 1,2
               nbl = nblk(ib,int)
               if (n.eq.nbl) then
                  nbc        = nbc + 1
                  itype(nbc) = 1
                  if (limblk(ib,1,int) .eq. limblk(ib,4,int)) then
                     idir(nbc) = 3
                     if (limblk(ib,1,int) .eq. ldimg(nbl)) 
     .                  idir(nbc) = -3
                     jbcs(nbc) = min(limblk(ib,2,int),limblk(ib,5,int))
                     jbce(nbc) = max(limblk(ib,2,int),limblk(ib,5,int))
                     lbcs(nbc) = limblk(ib,1,int)
                     lbce(nbc) = limblk(ib,1,int)
                     kbcs(nbc) = min(limblk(ib,3,int),limblk(ib,6,int))
                     kbce(nbc) = max(limblk(ib,3,int),limblk(ib,6,int))
                  else if (limblk(ib,2,int) .eq. limblk(ib,5,int)) then
                     idir(nbc) = 1
                     if (limblk(ib,1,int) .eq. jdimg(nbl))
     .                   idir(nbc) = -1
                     jbcs(nbc) = limblk(ib,2,int)
                     jbce(nbc) = limblk(ib,2,int)
                     kbcs(nbc) = min(limblk(ib,3,int),limblk(ib,6,int))
                     kbce(nbc) = max(limblk(ib,3,int),limblk(ib,6,int))
                     lbcs(nbc) = min(limblk(ib,1,int),limblk(ib,4,int))
                     lbce(nbc) = max(limblk(ib,1,int),limblk(ib,4,int))
                  else if (limblk(ib,3,int) .eq. limblk(ib,6,int)) then
                     idir(nbc) = 2
                     if (limblk(ib,1,int) .eq. kdimg(nbl))
     .                   idir(nbc) = -2
                     jbcs(nbc) = min(limblk(ib,2,int),limblk(ib,5,int))
                     jbce(nbc) = max(limblk(ib,2,int),limblk(ib,5,int))
                     kbcs(nbc) = limblk(ib,3,int)
                     kbce(nbc) = limblk(ib,3,int)
                     lbcs(nbc) = min(limblk(ib,1,int),limblk(ib,4,int))
                     lbce(nbc) = max(limblk(ib,1,int),limblk(ib,4,int))
                  end if
               end if
            end do
         end do
c
c        finally, convert all patch data to peg bc type 1
c
         do int=1,abs(ninter)
            if (igridg(nbl_to(int)).eq.n .and. 
     .         levelg(nbl_to(int)).eq.1) then
               nbc        = nbc + 1
               itype(nbc) = 1
               iface = nm_to(int)/10
               iend  = nm_to(int) - iface*10
               if (iface .eq. 1) then
                  jbcs(nbc) = jxie_s(int)
                  jbce(nbc) = jxie_e(int)
                  kbcs(nbc) = keta_s(int)
                  kbce(nbc) = keta_e(int)
                  if (iend .eq. 1) then
                     lbcs(nbc) = 1
                     lbce(nbc) = 1
                     idir(nbc) = 3
                  else
                     lbcs(nbc) = -1
                     lbce(nbc) = -1
                     idir(nbc) = -3
                  end if
               else if (iface .eq. 2) then
                  kbcs(nbc) = jxie_s(int)
                  kbce(nbc) = jxie_e(int)
                  lbcs(nbc) = keta_s(int)
                  lbce(nbc) = keta_e(int)
                  if (iend .eq. 1) then
                     jbcs(nbc) = 1
                     jbce(nbc) = 1
                     idir(nbc) = 1
                  else
                     jbcs(nbc) = -1
                     jbce(nbc) = -1
                     idir(nbc) = -1
                  end if
               else if (iface .eq. 3) then
                  lbcs(nbc) = jxie_s(int)
                  lbce(nbc) = jxie_e(int)
                  jbcs(nbc) = keta_s(int)
                  jbce(nbc) = keta_e(int)
                  if (iend .eq. 1) then
                     kbcs(nbc) = 1
                     kbce(nbc) = 1
                     idir(nbc) = 2
                  else
                     kbcs(nbc) = -1
                     kbce(nbc) = -1
                     idir(nbc) = -2
                  end if
               end if 
            end if
         end do
c
c        write out peg bc's to peg.bc file
c
         if (iname .eq. 0) then
c
            write(6,'(''enter name for mesh '',i4,
     .                '' (up to 40 char.)'')') n
            read(5,'(a40)') string
c
         else
            if (n.gt.99) then
               len1 = 8
               write(string,'("grid.",i3)') n
            else if (n.gt.9) then
               len1 = 7
               write(string,'("grid.",i2)') n
            else
               len1 = 6
               write(string,'("grid.",i1)') n
            endif
            do i = len1+1, 40
               string(i:i) = ' '
            end do
         end if
         write(7,'(a40)') string
         if (isym(n) .eq. 0) then
            write(7,'(7i5)') nbc,izeero,izeero,izeero,
     .                       izeero,i2d,izeero
         else
            if (ialph .eq. 0) then 
                write(7,'(7i5)') nbc,izeero,ione,izeero,
     .                           izeero,i2d,izeero 
            else
                write(7,'(7i5)') nbc,izeero,izeero,ione,
     .                           izeero,i2d,izeero
            end if
         end if
         write(7,'(8i5)') (itype(nb),idir(nb),jbcs(nb),jbce(nb),
     .                     kbcs(nb),kbce(nb),lbcs(nb),lbce(nb),nb=1,nbc)
      end do
c
c     free memory
c
      deallocate(jdimg)
      deallocate(kdimg)
      deallocate(ldimg)
      deallocate(nbci0)
      deallocate(nbcj0)
      deallocate(nbck0)
      deallocate(nbcidim)
      deallocate(nbcjdim)
      deallocate(nbckdim)
      deallocate(ibcinfo)
      deallocate(jbcinfo)
      deallocate(kbcinfo)
      deallocate(nblk)
      deallocate(limblk)
      deallocate(isva)
      deallocate(nblon)
      deallocate(iovrlp)
      deallocate(isym)
      deallocate(itype)
      deallocate(idir)
      deallocate(jbcs)
      deallocate(jbce)
      deallocate(kbcs)
      deallocate(kbce)
      deallocate(lbcs)
      deallocate(lbce)
      deallocate(nbl_to)
      deallocate(nm_to)
      deallocate(jxie_s)
      deallocate(jxie_e)
      deallocate(keta_s)
      deallocate(keta_e)
      deallocate(igridg)
      deallocate(levelg)
c
      return
      end
c
      subroutine shortinp(nb1,ne1,nb2,ne2)
c
c     converts cfl3d "shortcut" input of full range of index values
c     with pegsus/overflow "shortcut" input where -1 means the
c     max index
c
      if (nb1.eq.0) nb1 = 1
      if (nb2.eq.0) nb2 = 1
      if (ne1.eq.0) ne1 = -1
      if (ne2.eq.0) ne2 = -1
c
      return
      end
c
      subroutine get_cmd_args(inpfile,ierr)
c
      character*80 inpfile,arg
c
c     default file: cfl3d.inp
c
      inpfile = 'cfl3d.inp'
c
      ierr = 0
      num_args = iargc()
      if (num_args .eq. 1) then
#if defined CRAY_TIME
         call pxfgetarg (1, inpfile, 3, ierror)
#else
         call getarg (1, inpfile)
#endif
            read (inpfile, 100, err = 110) arg
 100        format (a80)
            
            go to 111
 110        continue
            ierr = 1
 111        continue
      else if (num_args .gt. 1) then
         ierr = 1
      end if
c
      return
      end
