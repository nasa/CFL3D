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
      program main
c
c     $Id$
c
c***********************************************************************
c
c     preprocessor for the patched-grid utility for cfl3d (ronnie) 
c
c     If you find a coding error or have suggestions for improved 
c     implementation, contact Bob Biedron, 804-864-2156
c                                          r.t.biedron@larc.nasa.gov
c
c***********************************************************************
c
c  input file:
c
c        unit 5       input file which contains setup information 
c                     (ascii) read from standard input
c
c  output files:
c
c       "preronnie.out"  name specified by user in main input file
c                        primary information file (ascii)
c
c        unit 99         contains error flag (= 0, normal termination,
c                                             < 0, abnormal termination)
c
c     maxbl  - maximum number of blocks
c     maxgr  - maximum number of grids
c     mwork  - size of 1-d array used to store xyz
c     nwork  - size of work array
c     intmax - maximum number of block interpolations
c     mpatch - maximum dimension of any block face involved in a patch
c     maxxe  - size of 1-d array used to store interpolation 
c              coefficients for all patch interfaces, including
c              those on coarser blocks
c     nsub1  - maximum number of blocks a single patch face may be
c              interpolated from
c
c     see main routine for cfl3d for definition/usage of the internal
c     buffer bou
c
c**********************************************************************
c
      parameter (ibufdim = 2000,nbuf = 4)
c
      external usrint
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist
      character*80 inpstring
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
      common /unit5/ iunit5
c
      iunit5 = 55
c
c***********************************************************************
c     set some signals for receiving a SIGTERM
c***********************************************************************
c
c      8...arithmetic exception
c      9...killed
c     10...bus error
c     11...segmentation fault
c
#if defined IBM
      call signal(8,usrint)
      call signal(9,usrint)
      call signal(10,usrint)
      call signal(11,usrint)
#else
      call signal(8,usrint,-1)
      call signal(9,usrint,-1)
#  if defined LINUX
      call signal(7,usrint,-1)
#  else
      call signal(10,usrint,-1)
#  endif
      call signal(11,usrint,-1)
#endif
c
c***********************************************************************
c     initialization
c***********************************************************************
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
      open(iunit5,file='tempz3y2x1.inp',form='formatted',
     .status='unknown')
c
      do n=1,9999
         read(5,'(a80)',end=999) inpstring
         write(iunit5,'(a80)') inpstring
      end do
  999 continue
      rewind(iunit5)
c
      read(iunit5,*)
c
      read(iunit5,'(a60)')grid
      read(iunit5,'(a60)')output
      open(unit=25,file=output,form='formatted',status='unknown')
c
      grdmov = output
c
      read(iunit5,'(a60)')patch
      open(unit=2,file=patch,form='unformatted',status='unknown')
c
      errfile='preronnie.error'
      open(unit=99,file=errfile,form='formatted',status='unknown')
c
c***********************************************************************
c     determine array size requirements
c***********************************************************************
c
c     read part of the the input file to get the primary array
c     dimensions needed by sizer for this case
c
      read(iunit5,*)
      read(iunit5,*) ioflag,itrace
      read(iunit5,*)
      read(iunit5,*)
      read(iunit5,*) ngrid
      ngrid = abs(ngrid)
c
      maxgr0   = ngrid
c
      read(iunit5,*)
      ncgmax = 0
      do ll=1,ngrid
         read(iunit5,*) ncg,iem,idim,jdim,kdim
         ncgmax = max(ncgmax,ncg)
      end do
c
      maxbl0   = ngrid*(ncgmax+1)
c
      read(iunit5,*)
      read(iunit5,*) nint0
c
      intmax0 = nint0*(ncgmax+1)
c
      read(iunit5,*) 
      do n=1,nint0
         read(iunit5,*)int,ifit,limit,itmax,mcxie,mceta,ic0,iorph
      end do
c
      nsub10 = 1
      if (nint0.gt.0) then
          if (ioflag .eq. 0) then
             read(iunit5,*)
             do n=1,nint0
                read(iunit5,*) int,ito,i1,i2,j1,j2,nfb
                nsub10 = max(nsub10,nfb)
             end do
          end if
          if (ioflag .eq. 1) then
             read(iunit5,*)
             read(iunit5,*)
             do n=1,nint0
                read(iunit5,*) int,ito,i1,i2,j1,j2,nfb
                do nn=1,nfb
                   read(iunit5,*)
                end do
                nsub10 = max(nsub10,nfb)
             end do
          end if
          if (ioflag .eq. 2) then
             do n=1,nint0
                read(iunit5,*)
                read(iunit5,*) int,ito,i1,i2,j1,j2,nfb
                do nn=1,nfb
                   read(iunit5,*)
                   read(iunit5,*)
                   read(iunit5,*)
                   read(iunit5,*)
                end do
                nsub10 = max(nsub10,nfb)
             end do
          end if
      end if
c
      rewind(iunit5)
c
c     the remaining *precfl3d* parameters *should* be large
c     enough for just about any case
c
      nwork0   = 1
      mpatch0  = 1
      maxxe0   = 1
      ibufdim0 = ibufdim
      nbuf0    = nbuf
c
      imode = 0
      call sizer(maxgr,maxbl,maxxe,nsub1,intmax,mpatch,nwork,
     .           nbuf0,ibufdim0,maxgr0,maxbl0,maxxe0,nsub10,
     .           intmax0,mpatch0,imode)
c
c***********************************************************************
c     normal program termnination
c***********************************************************************
c
c     remove the temporary input file
c
      close(iunit5,status='delete')
c
      call termn8(myid,0,ibufdim,nbuf,bou,nou)
c
      stop
      end
