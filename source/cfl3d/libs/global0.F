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
      subroutine global0(nplots0,maxnode0,mxbli0,lbcprd0,lbcemb0,
     .                   lbcrad0,maxbl0,maxgr0,maxseg0,maxcs0,ncycmax0,
     .                   intmax0,nsub10,intmx0,mxxe0,mptch0,msub10,
     .                   ibufdim0,nbuf0,mxbcfil0,nmds0,maxaes0,
     .                   maxsegdg0,ntr,nnodes,nou,bou,iunit11,myid,
     .                   idm0,jdm0,kdm0)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Scan/Read in the case input data in order to set some
c     sizing parameters for precfl3d (sizer). As each line is read,
c     it is output to the error file (iunit11), and the output buffer
c     is flushed. This helps insure that in the event of an omission
c     in the input file, the end of error file will fairly closely
c     indicate the location of the omission
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim0,nbuf0)
c
      real realval(20)
c
      dimension titlw(20),nou(nbuf0)
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /unit5/ iunit5
      common /elastic_ss/ idef_ss
c
      do nn=1,nbuf0
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
      call readkey(ititr,myid,ibufdim0,nbuf0,bou,nou,0,-99)
c
c     read title
c
      if (ititr.eq.0) then
         call echoinp(iunit5,iunit11,1)
      end if
c
c     read input file up to irest
c
      nread = 5
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c
      read(iunit5,*) realval(1),irest,iflagts,realval(2),iunst,
     .               realval(3)
      write(iunit11,'(f10.5,2i10,f10.5,i10,f10.5)') realval(1),irest,
     .     iflagts,realval(2),iunst,realval(3) 
      call echoinp(iunit5,iunit11,0)
      dt     = realval(1)
      if (real(dt).lt.0.) iunst = 0
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) maxgr0,nplot3d,nprint,nwrest,ichkd,i2d,ntstep,ita
      write(iunit11,'(8i10)') maxgr0,nplot3d,nprint,nwrest,ichkd,
     .     i2d,ntstep,ita
      call echoinp(iunit5,iunit11,0)
c
      maxgr0 = abs(maxgr0)
c
      call echoinp(iunit5,iunit11,1)
c
c     read ncg data
      ncgmax = 0
      do n=1,maxgr0
         read(iunit5,*) ncg,idum1,idum2,idum3,idum4,idum5,idum6
         write(iunit11,'(7i10)') ncg,idum1,idum2,idum3,idum4,idum5,idum6
         call echoinp(iunit5,iunit11,0)
         ncgmax = max(ncg,ncgmax)
      end do
c
      maxbl0 = (ncgmax+1)*maxgr0
      maxnode0 = max(maxgr0,nnodes)
c
      idm0 = 0
      jdm0 = 0
      kdm0 = 0
      read(iunit5,*)
      nread = maxgr0
      do n=1,nread
         read(iunit5,*) idm,jdm,kdm
         if (idm.gt.idm0) then
            idm0 = idm
         end if
         if (jdm.gt.jdm0) then
            jdm0 = jdm
         end if
         if (kdm.gt.kdm0) then
            kdm0 = kdm
         end if
      end do
      nread = maxgr0+1
c
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
      maxseg0 = 1
      nreadi0 = 0
      nreadid = 0
      nreadj0 = 0
      nreadjd = 0
      nreadk0 = 0
      nreadkd = 0
c
      call echoinp(iunit5,iunit11,1)
      nread = maxgr0
      do n=1,nread
         read(iunit5,*) mdum,mbci0,mbcidim,mbcj0,mbcjdim,
     .                  mbck0,mbckdim
         write(iunit11,'(8i10)') mdum,mbci0,mbcidim,mbcj0,mbcjdim,
     .        mbck0,mbckdim
         call echoinp(iunit5,iunit11,0)
         maxseg0 = max(maxseg0,mbci0,mbcidim,mbcj0,mbcjdim,
     .                 mbck0,mbckdim)
         nreadi0 = nreadi0+mbci0
         nreadid = nreadid+mbcidim
         nreadj0 = nreadj0+mbcj0
         nreadjd = nreadjd+mbcjdim
         nreadk0 = nreadk0+mbck0
         nreadkd = nreadkd+mbckdim
      end do
c
c     read through bc section, counting number of solid surface,
c     periodic, and radial equilibrium bcs
c
      nprd    = 0
      nrad    = 0
      nsol    = 0
      nsolgd  = 0
      iglast  = 0
      lbcrad0 = 1
      lbcprd0 = 1
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadi0
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadid
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadj0
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadjd
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadk0
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
      call echoinp(iunit5,iunit11,1)
      do n=1,nreadkd
         read(iunit5,*) igrid,idum2,ibctyp,idum4,idum5,idum6,idum7,ndat
         write(iunit11,'(8i10)') igrid,idum2,ibctyp,idum4,idum5,
     .        idum6,idum7,ndat
         call echoinp(iunit5,iunit11,0)
         if (igrid.ne.iglast) nsolgd = 0
         if (abs(ibctyp).eq.2005) nprd = nprd + 1
         if (abs(ibctyp).eq.2006) nrad = nrad + 1
         if (abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .       abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .       abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .       abs(ibctyp).eq.2016) then
            nsol   = nsol + 1
            nsolgd = nsolgd + 1
         end if
         if (abs(ndat) .gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         iglast = igrid
      end do
c
      lbcprd0 = max(lbcprd0,nprd*(ncgmax+1))
      lbcrad0 = max(lbcrad0,nrad*(ncgmax+1))
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
      nread = 2*nbli0 + 2
      do n=1,nread
         call echoinp(iunit5,iunit11,1)
      end do
c
      mxbli0 = max(1,nbli0*(ncgmax+1))
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
      else 
         nsub10 = 0
         read(22) nint0
         call echoinp(iunit5,iunit11,0)
         do  n=1,abs(nint0)
            read(22) nsb10
            call echoinp(iunit5,iunit11,0)
            nsub10 = max(nsub10,nsb10)
            do nn=1,12
               read(22)
            end do
         end do
         rewind(22)
c
         intmax0 = nint0*(ncgmax+1)
c
      end if
c
c     plot3d, printout, and control surface data
c
      nplots0 = 1
      if (nplot3d.gt.0) then
         nplots0 = nplot3d
      else if (nplot3d.lt.0) then
         nplots0 = max(nsol,abs(nplot3d))
         if (abs(i2d).gt.0) nplots0 = max(maxgr0,abs(nplot3d))
      end if
      if (nprint.gt.0) then
         nplots0 = max(nplots0,nprint)
      else if (nprint.lt.0) then
         nplots0 = max(nplots0,nsol)
      end if
c
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      do n=1,abs(nplot3d)
         call echoinp(iunit5,iunit11,1)
      end do
c
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) movie
      write(iunit11,'(i10)') movie
      call echoinp(iunit5,iunit11,0)
c
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      do n=1,abs(nprint)
         call echoinp(iunit5,iunit11,1)
      end do
c
      call echoinp(iunit5,iunit11,1)
      call echoinp(iunit5,iunit11,1)
      read(iunit5,*) ncs
      write(iunit11,'(i10)') ncs
      call echoinp(iunit5,iunit11,0)
      call echoinp(iunit5,iunit11,1)
      do n=1,ncs
         call echoinp(iunit5,iunit11,1)
      end do
c
      maxcs0 = max(1,ncs)
c
c     rigid grid motion data 
c
      if (iunst.eq.1 .or. iunst.eq.3) then
c        translation data
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) ntrans
         write(iunit11,'(i10)') ntrans
         call echoinp(iunit5,iunit11,0)
         call echoinp(iunit5,iunit11,1)
         if (ntrans.gt.0) then
            call echoinp(iunit5,iunit11,1)
         end if
         call echoinp(iunit5,iunit11,1)
         do n=1,ntrans
             call echoinp(iunit5,iunit11,1)
         end do
         call echoinp(iunit5,iunit11,1)
         do n=1,ntrans
             call echoinp(iunit5,iunit11,1)
         end do
c        rotation data
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) nrotat
         write(iunit11,'(i10)') nrotat
         call echoinp(iunit5,iunit11,0)
         call echoinp(iunit5,iunit11,1)
         if (nrotat.gt.0) then
            call echoinp(iunit5,iunit11,1)
         end if
         call echoinp(iunit5,iunit11,1)
         do n=1,nrotat
             call echoinp(iunit5,iunit11,1)
         end do
         call echoinp(iunit5,iunit11,1)
         do n=1,nrotat
             call echoinp(iunit5,iunit11,1)
         end do
      end if
c
c    deforming grid motion data
c
      nmds0     = 1
      maxaes0   = 1
      maxsegdg0 = 1
c
      if (iunst.gt.1 .or. idef_ss.gt.0) then
c        deforming mesh data
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) ndefrm
         write(iunit11,'(i10)') ndefrm
         call echoinp(iunit5,iunit11,0)
         if (abs(ndefrm).gt.0) then
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
            do ndef=1,abs(ndefrm)
               call echoinp(iunit5,iunit11,1)
            end do
            call echoinp(iunit5,iunit11,1)
            do ndef=1,abs(ndefrm)
               call echoinp(iunit5,iunit11,1)
            end do
         else
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
         if (ndefrm.gt.0) then
            maxsegdg0 = nsol + ndefrm
         else
            maxsegdg0 = nsol
         end if
c        aeroelastic mesh data
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) naesrf
         write(iunit11,'(i10)') naesrf
         call echoinp(iunit5,iunit11,0)
         maxaes0 = max(maxaes0,naesrf)
         if (naesrf.gt.0) then
            do naes=1,naesrf
               call echoinp(iunit5,iunit11,1)
               read(iunit5,*) iaes,ngd,(realval(i),i=1,3),nmodes,iskyhk
               write(iunit11,'(i8,i9,3f9.4,2i9)') iaes,ngd,
     .              (realval(i),i=1,3),nmodes,iskyhk
               call echoinp(iunit5,iunit11,0)
               nmds0 = max(nmds0,nmodes)
               call echoinp(iunit5,iunit11,1)
               do nm = 1,nmodes
                  call echoinp(iunit5,iunit11,1)
               end do
               call echoinp(iunit5,iunit11,1)
               do nm = 1,nmodes
                  call echoinp(iunit5,iunit11,1)
               end do
               call echoinp(iunit5,iunit11,1)
               do ng = 1,abs(ngd)
                  call echoinp(iunit5,iunit11,1)
               end do
               if (ngd.ge.0) then
                  maxsegdg0 = maxsegdg0 + ngd
               end if
            end do
         else
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
            call echoinp(iunit5,iunit11,1)
         end if
c
c        offbody/multiblock mesh deformation data
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) nskip,idum1,realval(1),realval(2),realval(3)
     .                  ,realval(4),idum2
         write(iunit11,'(2i8,4f9.6,i9)') nskip,idum1,realval(1),
     .                  realval(2),realval(3),realval(4),idum2
         if(abs(idum1).eq.1) then
           call echoinp(iunit5,iunit11,1)
           if (abs(nskip).gt.0) then
              do ng = 1,abs(nskip)
                 call echoinp(iunit5,iunit11,1)
              end do
           end if
         else
           call echoinp(iunit5,iunit11,1)
           if(abs(nskip).gt.0) then
             do ng = 1,abs(nskip)
               call echoinp(iunit5,iunit11,1)
               read(iunit5,*) igr,nskpi,nskpj,nskpk
               call echoinp(iunit5,iunit11,1)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpi) iie = nskpi
                 call echoinp(iunit5,iunit11,1)
                 if(iie.eq.nskpi) goto 325
               enddo
325            continue
               call echoinp(iunit5,iunit11,1)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpj) iie = nskpj
                 call echoinp(iunit5,iunit11,1)
                 if(iie.eq.nskpj) goto 350
               enddo
350            continue
               call echoinp(iunit5,iunit11,1)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpk) iie = nskpk
                 call echoinp(iunit5,iunit11,1)
                 if(iie.eq.nskpk) goto 375
               enddo
375            continue
             enddo
           end if
         end if
c
c        multi-motion coupling
         call echoinp(iunit5,iunit11,1)
         call echoinp(iunit5,iunit11,1)
         read(iunit5,*) ncoupl
         write(iunit11,'(i8,f9.6,2i9)') ncoupl
         call echoinp(iunit5,iunit11,0)
         call echoinp(iunit5,iunit11,1)
         if (abs(ncoupl).gt.0) then
            do ng = 1,abs(ncoupl)
               call echoinp(iunit5,iunit11,1)
            end do
         end if
c
         maxsegdg0 = max(1,maxsegdg0)
c
      end if
c
c     dynamic patch data
c
      intmx0 = 1
      msub10 = 1
      nint1  = 0
c
      if (iunst.gt.0) then
         read(iunit5,'(a80)',end=888) string
         write(iunit11,'(a80)') string
         call echoinp(iunit5,iunit11,0)
         read(iunit5,'(a80)',end=888) string
         write(iunit11,'(a80)') string
         call echoinp(iunit5,iunit11,0)
         read(iunit5,*,end=888) nint1
         write(iunit11,'(i10)') nint1
         call echoinp(iunit5,iunit11,0)
  888    continue
c
         if (nint1.gt.0) then
c
             intmx0 = nint1*(ncgmax+1)
c
             call echoinp(iunit5,iunit11,1)
             do n=1,nint1
                call echoinp(iunit5,iunit11,1)
             end do
             do n=1,nint1
                call echoinp(iunit5,iunit11,1)
                read(iunit5,*) int,ito,i1,i2,j1,j2,nfb
                write(iunit11,'(8i10)') int,ito,i1,i2,j1,j2,nfb
                call echoinp(iunit5,iunit11,0)
                do nn=1,nfb
                   call echoinp(iunit5,iunit11,1)
                   call echoinp(iunit5,iunit11,1)
                   call echoinp(iunit5,iunit11,1)
                   call echoinp(iunit5,iunit11,1)
                end do
                msub10 = max(msub10,nfb)
             end do
         end if
c
      end if
c
      nsub10 = max(nsub10,msub10)
c
c     if a restart, read restart file to determine how many cycles
c     or time steps have already been done
c
      ntr = 0
      if (irest.ne.0) then
         if (icgns .ne. 1) then
         read(2,end=999) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
 999     continue
         rewind(2)
         close(2)
         else
#if defined CGNS
c          the following opens cgns file, gets ntr, and closes it
c          it assumes only one base
           call getntr(grid,ntr,istp)
           if (istp .eq. 1) then
             call termn8(myid,-99,ibufdim0,nbuf0,bou,nou)
           end if
#endif
         end if
      end if
      ncycmax0 = max(ncyctot+ntr,1)
c
c     set the remaining parameters neede by sizer (precfl3d).
c     values of 1 for mxxe and mptch are sufficient for sizer;
c     don't have a good way of setting lbcemb0 based on the minimal
c     data read by this routine. The following value should be big
c     enough for any possible embedded case...sizer will set the
c     exact values for cfl3d
c 
      mxxe0    =          1
      mptch0   =          1
      lbcemb0  =   6*maxbl0
c
c     print out the current parameters, if desired, or
c     if one or more turn out to be zero (an error)
c
      iprint = 0
      ichek = 1
c
      if (nplots0   .le. 0 .or.
     .    maxnode0  .le. 0 .or.
     .    mxbli0    .le. 0 .or.
     .    lbcprd0   .le. 0 .or.
     .    lbcemb0   .le. 0 .or.
     .    lbcrad0   .le. 0 .or.
     .    maxbl0    .le. 0 .or.
     .    maxgr0    .le. 0 .or.
     .    maxseg0   .le. 0 .or.
     .    maxcs0    .le. 0 .or.
     .    ncycmax0  .le. 0 .or.
     .    intmax0   .le. 0 .or.
     .    nsub10    .le. 0 .or.
     .    intmx0    .le. 0 .or.
     .    mxxe0     .le. 0 .or.
     .    mptch0    .le. 0 .or.
     .    msub10    .le. 0 .or.
     .    ibufdim0  .le. 0 .or.
     .    nbuf0     .le. 0 .or.
     .    maxsegdg0 .le. 0 .or.
     .    mxbcfil0  .le. 0) ichek = 0
c
      if (ichek.eq.0) iprint = 1
c
      if (iprint.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'nplots0   = ',nplots0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxnode0  = ',maxnode0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'mxbli0    = ',mxbli0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'lbcprd0   = ',lbcprd0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'lbcemb0   = ',lbcemb0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'lbcrad0   = ',lbcrad0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxbl0    = ',maxbl0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxgr0    = ',maxgr0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxseg0   = ',maxseg0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxcs0    = ',maxcs0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'ncycmax0  = ',ncycmax0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'intmax0   = ',intmax0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'nsub10    = ', nsub10
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'intmx0    = ',intmx0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'mxxe0     = ',mxxe0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'mptch0    = ',mptch0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'msub10    = ',msub10
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'ibufdim0  = ',ibufdim0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'nbuf0     = ',nbuf0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'mxbcfil0  = ',mxbcfil0
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),*)'maxsegdg0 = ',maxsegdg0
      end if
c
      if (ichek.eq.0) then
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),'(''error in routine global0 - one of '',
     .                     ''the parameters listed above is zero'')')
         ierrflg = -99
         call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
      end if
c
      rewind(iunit11)
c
      return
      end

      subroutine echoinp(iunit5,iunit11,irw)
c***********************************************************************
c     Purpose: if irw > 0, read and write one line from the input file,
c     and flush the output buffer of unit iunit11; if irw = 0, simply
c     flush the buffer. This helps insure that if there is an error in
c     the input file, the end of the file (iunit11) will fairly
c     closely reflect the location of the error.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 string
c
      if (irw .gt. 0) then
         read(iunit5,'(a80)') string
         write(iunit11,'(a80)') string
      end if
      call my_flush(iunit11)
c
      return
      end
