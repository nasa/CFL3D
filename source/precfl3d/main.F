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
c
c***********************************************************************
c
c     PRECFL3D: preprocessor to cfl3d for determining memory requirements
c
c     VERSION 6.7 :  Computational Fluids Laboratory, Mail Stop 128,
c                    NASA Langley Research Center, Hampton, VA
c
c***********************************************************************
c     this program reads in a cfl3d input file (version 5/6 format) and
c     evaluates the parameter sizes required for version 6.0 of cfl3d,
c     sized for the particular problem in the cfl3d input file that is
c     specified as standard input to precfl3d.
c
c     usage: precfl3d -np NUM_PROCS <input_file
c
c     where: NUM_PROCS is the total number of processors on which
c            the case will be run (only NUM_PROCS-1 processors are
c            used for solution updates; the additional processor
c            serves as a host, performing I/O). For building only the
c            sequential version of the code, specify NUM_PROCS=1
c
c     note: In version 6, dynamic memory allocation has been implemented
c           in cfl3d, so this preprocessor is needed only to determine
c           apriori the memory that cfl3d will request at execution time.
c
c     Description of the files read/written.
c
c     unit  2 - Input/output file containing binary restart data.
c     unit  5 - Input case data file.
c     unit 21 - Input file for overlapped grid information,
c     unit 22 - Input file for patched grid information.
c     unit 32 - Temporary (tempz3y2x1.inp, copy of input file, deleted at
c               normal termination)
c     unit 66 - Primary output file - contains detailed info on how the
c               permanent and temporary work array dimensions were
c               determined
c     unit 67 - Secondary output file - contains estimation of the 
c               speedup possible for the case at hand using multiple
c               processors
c
c     see main routine for cfl3d for definition/usage of the internal
c     buffer bou
c
c***********************************************************************
c
      parameter (ibufdim = 2000,nbuf = 4, mxbcfil=100)
c
      external usrint
c
      character*120 bou(ibufdim,nbuf)
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist
      character*80  inpstring
      character*80 bcfiles(mxbcfil)
c
      dimension nou(nbuf)
c
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /unit5/ iunit5
      common /proces/ numprocs
      common /zero/ iexp
c
      call get_cmd_args(ierr)
      if (ierr .eq. 1) then
         write(6,599)
 599     format('error in command line input for precfl3d,',
     .   ' version 6.0')
         write(6,*)
         write(6,600)
 600     format('correct usage: precfl3d -np num_nodes < input_file')
         write(6,601)
 601     format('   where num_nodes is the number of nodes',
     .   ' (processors), INCLUDING host')
         write(6,*)
         write(6,602)
 602     format('to size only the sequential code, use -np 1')
         call exit (1)
      endif
c
      myhost     = 0
      myid       = 0
      mycomm     = 0
      if (numprocs.gt.1) then
         nnodes = numprocs - 1
      else
         nnodes = 1
      end if
c
      do ii=1,nbuf
         nou(ii) = 0
         do mm=1,ibufdim
            bou(mm,ii) = ' '
         end do
      end do
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
      call signal(10,usrint,-1)
      call signal(11,usrint,-1)
#endif
c
c***********************************************************************
c     open files
c***********************************************************************
c
      iunit5 = 32
c
#if defined NOREDIRECT
      open(iunit5,file='cfl3d.inp',form='formatted',status='old')
# else
      open(iunit5,file='tempz3y2x1.inp',form='formatted',
     .status='unknown')
c
c     first copy the command-line input file to a temporary input file
c     this is necessary because the input file must be read more
c     more than once to accomodate dynamic memory allocation, and
c     some operating systems do not permit a rewind of a file
c     that has not been explicitly opened
c
c     NOTE: this assumes the cfl3d input file contains lines of no
c           more than 80 characters!!!
c
      do n=1,9999
         read(5,'(a80)',end=999) inpstring
         write(iunit5,'(a80)') inpstring
      end do
  999 continue
      rewind(iunit5)
# endif
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
c     the following file names are not currently set from the input
c
      preout  = 'precfl3d.out'
      errfile = 'precfl3d.error'
      grdmov  = 'cfl3d.dynamic_patch'
c
      open(unit=21,file=ovrlap,form='unformatted',status='unknown')
      open(unit=22,file=patch,form='unformatted',status='unknown')
      open(unit=25,file=grdmov,form='formatted',status='unknown')
      open(unit=2,file=restrt,form='unformatted',status='unknown')
      open(unit=66,file=preout,form='formatted',status='unknown')
      open(unit=99,file=errfile,form='formatted',status='unknown')
c
c***********************************************************************
c     determine machine zero for use in setting tolerances
c     (10.**(-iexp) is machine zero)
c***********************************************************************
c
      icount = 0
      compare = 1.0
      do i = 1,20
         icount = icount + 1
         add = 1.
         do n=1,i
            add = add*.1
         enddo
         x11 = compare + add
         if (x11.eq.compare)then
            iexp = i-1
            goto 4010
         end if
      end do
 4010 continue
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
c     sizer sets the parameters needed for cfl3d
c
      imode = 0
      call sizer(mwork,mworki,nplots,minnode,iitot,intmax,
     .           maxxe,mxbli,nsub1,lbcprd,lbcemb,lbcrad,
     .           maxbl,maxgr,maxseg,maxcs,ncycmax,intmx,mxxe,
     .           mptch,msub1,nmds,maxaes,maxsegdg,nnodes,nslave,
     .           nmaster,myhost,myid,mycomm,nplots0,maxnode0,mxbli0,
     .           lbcprd0,lbcemb0,lbcrad0,maxbl0,maxgr0,maxseg0,
     .           maxcs0,ncycmax0,intmax0,nsub10,intmx0,mxxe0,
     .           mptch0,msub10,ibufdim0,nbuf0,mxbcfil0,nmds0,
     .           maxaes0,maxsegdg0,imode,ntr,bcfiles,bou,nou)
c
c***********************************************************************
c     normal program termnination
c***********************************************************************
c
c     remove the temporary input file, if used
c
#if defined NOREDIRECT
#else
      if (myid.eq.myhost) then
         close(iunit5,status='delete')
      end if
#endif
c
      call termn8(myid,0,ibufdim,nbuf,bou,nou)
c 
      end
