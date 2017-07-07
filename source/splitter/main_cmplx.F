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
c     Given a grid and a CFL3D input file (and corresponding ronnie
c     input file, if applicable) and grid, this program will split the
c     grid along user-specified grid lines to create a new grid file
c     with more blocks, and will split and create a new CFL3D input
c     file (and new ronnie input file) to correspond to the new, split
c     grid. For the versions of CFL3D that compute sensitivity derivatives,
c     this splitter can also split the "seed matrix" or grid sensitivity
c     file.
c
c     File formats: 1) grid file is PLOT3D GRID, multiblock style
c                      (formated or unformatted). NOTE: cfl3d style
c                      grid files are currently not supported!
c                   2) grid sensitivity file (if applicable)
c                      is PLOT3D FUNCTION, multiblock style
c
c***********************************************************************
c
c     Original coding by E. B. Parlette for TLNS3D (Summer 1995)
c
c     Modifications for CFL3D by R. T. Biedron (Winter 1998)
c
c     This is a modified version of the splitter originally developed
c     by Ed Parlette of Vigyan for TLNS3D. The original splitter has
c     been coupled to file conversion programs that convert CFL3D input
c     (and patch) files to TLNS3D map files, and vice versa. The result
c     is a block and input file splitter for CFL3D.
c
c     Caveats: The conversions from the CFL3D input file to a TLNS3D
c              map file are not perfect! The user is urged check the
c              resulting split CFL3D input (and patch) files. A useful
c              check before actually splitting the files is to run this
c              splitter with the number of splittings = 0, and the
c              output grid file as null. This will cause to code to go
c              through the translations, but the "split" files will
c              have the same numbers of blocks, and the split grid
c              will not be output. A "diff" or "gdiff" will point
c              to translation-induced differences that should be
c              easier to sort out than when coupled with true splitting.
c              Note that the 2-step procedure almost always results in a
c              *reordering* of some boundary condition segments.
c
c     Files and I/O unit numbers:
c
c              unit  5...main input file controling the splitting
c              unit 51...copy of unit 5 data (temporary)
c              unit  6...information output file
c
c              unit 30...input (unsplit) grid file
c              unit 10...input (unsplit) cfl3d input file
c              unit 15...input (unsplit) ronnie input file
c              unit 50...input (unsplit) grid sensitivity file
c
c              unit  7...output (split) cfl3d input file
c              unit  9...output (split) ronnie input file
c              unit 40...output (split) grid file
c              unit 60...output (split) grid sensitivity file
c
c              unit 20...input (unplit) tlns3d map file (temporary)
c              unit  8...output (split) tlns3d map file (temporary)
c
c        (the temporary files above are deleted upon program completion)
c
c***********************************************************************
c
c     the following parameters are needed by the sizer routine and 
c     should be quite sufficient for any conceivable case. note:
c     some additional parameters are set in the sizing routine; 
c     these should also be more than sufficient for any case, but
c     in case of difficulty, see subroutine sizer
c
c     mbloc0...max number of input blocks
c     msplt0...max number of splits that can be done on any given block
c
      parameter (mbloc0=2000,msplt0=200)
      parameter (ibufdim=2000,nbuf=4)
c
      character*80 inpstring
      character*120 bou(ibufdim,nbuf)
      dimension nou(nbuf)
c
      common /unit5/ iunit5
      common /zero/ iexp
c
      do ii=1,nbuf
         nou(ii) = 0
         do mm=1,ibufdim
            bou(mm,ii) = ' '
         end do
      end do
c
c     determine machine zero for use in setting tolerances
c     (10.**(-iexp) is machine zero)
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
c     copy command-line splitter inpu file to a temp file
c     that is explicitly opened so that this temp file can
c     be rewound
c
      iunit5 = 51
      open(unit=iunit5,file='tempz3y2x1.inp',form='formatted',
     .status='unknown')
c
      do n=1,9999
         read(5,'(a80)',end=999) inpstring
         write(iunit5,'(a80)') inpstring
      end do
  999 continue
      rewind(iunit5)
c
      open(unit=99,file='splitter.error',
     .form='formatted',status='unknown')
c
      call sizer(mbloc,maxseg,npmax,mxbli,intmax,nsub1,msegt,
     .           msegn,msplt,mtot,mbloc0,msplt0,nou,bou,nbuf,ibufdim)
c
      call splitter_cmplx(mbloc,maxseg,npmax,mxbli,intmax,nsub1,
     .                    msegt,msegn,msplt,mtot,nou,bou,nbuf,ibufdim)
c
      close(iunit5,status='delete')
c
      call termn8(0,0,ibufdim,nbuf,bou,nou)
      end
