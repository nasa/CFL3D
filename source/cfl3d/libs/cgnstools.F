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
c
c     $Id$
c
#if defined CGNS
      subroutine getsolinfo(iccg,ibase,igrid,nsoluse,irinddata)
c
c**********************************************************************
c     Purpose: Gets nsoluse and rind data from CGNS data base 
c     for zone number "igrid". The CGNS file must already be 
c     opened. If no solution node exists, the routine stops 
c     execution. If the solutions are not stored at CellCenter, 
c     the routine stops execution.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c   OUTPUTS
c      nsoluse..........CGNS "FlowSolution" index number (integer)
c      irinddata(6).....CGNS rind information (integer)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension irinddata(6)
c
      character*32 solname
c
c   Find out how many Flow solutions exist
      call cg_nsols_f(iccg, ibase, igrid, nsols, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nsols .gt. 1) then
        write(901,'('' nsols='',i5,''.  Expecting 1.'')') nsols
        write(901,'('' Will read LAST one.'')')
        nsoluse=nsols
      else if (nsols .eq. 1) then
        write(901,'('' solution node exists.  Reading it...'')')
        nsoluse=nsols
      else
        write(901,'('' Error.  Solution node does not exist!'')')
        stop
      end if
c
c   Get solution info
      call cg_sol_info_f(iccg, ibase, igrid, nsoluse, solname,
     +  location, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (location .ne. CellCenter) then
        write(901,'('' GridLocation must be CellCenter!  Currently:'',
     +   a32)') GridLocationName(location+1)
        stop
      end if
      write(901,'('' cgns read solution of name:  '',a32)') solname
c
c   Check if Rind data exists
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'FlowSolution_t',
     + nsoluse,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_rind_read_f(irinddata,ier)
      if (ier .ne. 0) then
        write(901,'('' rind info does not exist.  setting to zero'')')
        do n=1,6
          irinddata(n)=0
        enddo
      end if
c
      do n=1,6
        if (irinddata(n) .gt. 1) then
          write(901,'('' error... currently cannot handle rind'',
     +     '' data > 1'')')
          stop
        end if
      enddo
      if (irinddata(1).ne.irinddata(2) .or.
     +    irinddata(3).ne.irinddata(4) .or.
     +    irinddata(5).ne.irinddata(6)) then
        write(901,'('' error... currently must have same rind'',
     +    '' info on both sides of given index direction'')')
        write(901,'('' irinddata='',6i3)') irinddata(1),
     +    irinddata(2),irinddata(3),irinddata(4),irinddata(5),
     +    irinddata(6)
        stop
      end if
c
      return
      end
c
      subroutine getgrid(iccg,ibase,igrid,idim,jdim,kdim,wk,
     +                   ialph,x,y,z)
c***********************************************************************
c     Purpose: Gets x,y, and z from CGNS data base for zone number
c     "igrid" The CGNS file must already be opened. If there is a
c     grid-size inconsistency with what is expected (idim,jdim,kdim)
c     or if the 3 coordinates CoordinateX, CoordinateY, and CoordinateZ 
c     do not exist, the routine stops execution.
c     Gets results in the following order:
c
c          (jdim,kdim,idim) <- necessary for CFL3D!,
c
c     and assumes they are stored in the CGNS database in the order:
c     (idim,jdim,kdim)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idim,jdim,kdim...expected dimensions of this zone (zone "igrid")
c                       (integers)
c      wk...............working space needed, of dimension
c                       (idim*jdim*kdim) or larger (real)
c      ialph............parameter in CFL3D for determining whether 
c                       y or z is "up" (integer)
c   OUTPUTS:
c      x,y,z............coordinates, output in (j,k,i) order (real)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),
     +          z(jdim,kdim,idim)
      dimension wk(idim,jdim,kdim)
      dimension isize(3*3),irmin(3),irmax(3)
c
      character*32 zonename,coordname,testname(3)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Coordinate names that we are looking for:
        testname(1)='CoordinateX'
        testname(2)='CoordinateY'
        testname(3)='CoordinateZ'
c
c   Read general zone information
      call cg_zone_type_f(iccg, ibase, igrid, izonetype, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (izonetype .ne. Structured) then
        write(901,'('' Error, zone must be Structured type.'')')
        write(901,'(''   Currently, it is '',a32)')
     +    ZoneTypeName(izonetype)
        stop
      end if
      call cg_zone_read_f(iccg, ibase, igrid, zonename, isize, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if(isize(1) .ne. idim .or. isize(2) .ne. jdim .or.
     +   isize(3) .ne. kdim) then
        write(901,'('' Grid index inconsistencies:  isize='',3i5,
     +   ''idim,jdim,kdim='',3i5)') isize(1),isize(2),isize(3),
     +   idim,jdim,kdim
        write(901,'('' Be sure to order the zones alphabetically'',
     +   '' in the input file!'')')
        stop
      end if
c   Find out how many grid coordinates exist
      call cg_ncoords_f(iccg, ibase, igrid, ncoords, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (ncoords .ne. 3) then
        write(901,'('' ncoords='',i5,''.  Expecting 3.'')') ncoords
        stop
      end if
c   Check coordinate names in data base
      do icoord=1,3
        call cg_coord_info_f(iccg, ibase, igrid, icoord, itype,
     +    coordname,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(coordname .eq. testname(1) .or.
     +     coordname .eq. testname(2) .or.
     +     coordname .eq. testname(3)) then
          continue
        else
          write(901,'('' coordname of '',a32,'' unrecognized.'')') 
     +     coordname
          write(901,'('' Looking for CoordinateX, CoordinateY, and'',
     +     '' CoordinateZ'')')
          stop
        end if
      enddo
c   Set up array bounds:
      irmin(1)=1
      irmin(2)=1
      irmin(3)=1
      irmax(1)=idim
      irmax(2)=jdim
      irmax(3)=kdim
c   Read x,y,z
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateX', 
     &                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateX', 
     &                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Reorder i,j,k to j,k,i:
      call reorderg(1,1,idim,jdim,kdim,idim,jdim,kdim,wk,x)
c
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateY', 
     &                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateY', 
     &                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Reorder i,j,k to j,k,i:
      if (ialph .eq. 0) then
        call reorderg(1,1,idim,jdim,kdim,idim,jdim,kdim,wk,y)
      else
        call reorderg(1,1,idim,jdim,kdim,idim,jdim,kdim,wk,z)
      end if
c
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateZ', 
     &                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateZ', 
     &                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Reorder i,j,k to j,k,i:
      if (ialph .eq. 0) then
        call reorderg(1,1,idim,jdim,kdim,idim,jdim,kdim,wk,z)
      else
        call reorderg(1,1,idim,jdim,kdim,idim,jdim,kdim,-wk,y)
      end if
c
      return
      end
c
      subroutine getntr(cgnsdesired,ntr,istp)
c***********************************************************************
c     Purpose: Opens CGNS data file for READ ONLY, gets ntr,
c     then closes it. NOTE: ibase is assumed to be 1
c     (This routine is called during the precfl3d phase of
c     execution, so error messages are written to unit 66)
c
c   INPUT:
c     cgnsdesired.....name of CGNS file to open (character*80)
c
c   OUTPUT:
c     ntr.............number of iterations run so far
c     istp............=0 no error, =1 error (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      character*80 cgnsdesired
      character*32 name2
c
      istp=0
c   Open CGNS file
      write(901,'('' Opening CGNS database from cgnstools/getntr'')')
      call cg_open_f(cgnsdesired,CG_MODE_READ,iccg,ier)
      if (ier .ne. 0) then
        write(66,'('' Error opening cgns file'')')
        write(66,'('' Make sure the file on 1st line of input IS'',
     +   '' a cgns file!'')')
        call cg_error_print_f
        istp=1
        return
      end if
c   Goto base index
      ibase=1
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) then
        write(66,'('' Error attempting to find Base node'')')
        write(66,'('' Make sure the cgns file is written correctly'',
     +   '' (according to SIDS standards)'')')
        call cg_error_print_f
        istp=1
        return
      end if
c
c   Read number of iterations of data in ConvergenceHistory node
      call cg_convergence_read_f(ntr,name2,ier)
      if (ier .ne. 0) then
c        ConvergenceHistory node does not exist...
c        setting ntr=0 and continuing:
         ntr=0
      end if
c
c   Close CGNS file
      write(901,'('' Closing CGNS data base from cgnstools/getntr'')')
      call cg_close_f(iccg,ier)
      if (ier .ne. 0) then
        call cg_error_print_f
        istp=1
        return
      end if
c
      return
      end
c
      subroutine ropencgns(cgnsdesired,basedesired,idimdesired,iccg,
     .                     ibase,nzones)
c***********************************************************************
c     Purpose: Opens CGNS data file for READ ONLY
c     If cannot find base of name "basedesired" or if the dimensionality
c     is not "idimdesired", the routine stops execution.
c     Returns CGNS file index, base index, and number of zones
c
c   INPUT:
c     cgnsdesired.....name of CGNS file to open (character*80)
c     basedesired.....name of base desired (character*32)
c     idimdesired.....dimensionality desired (integer)
c     
c   OUTPUT:
c     iccg............index number of CGNS file (integer)
c     ibase...........index number of base desired (integer)
c     nzones..........number of zones in base desired (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      character*80 cgnsdesired
      character*32 basedesired,basename
c
c   Open CGNS file
      write(901,'(/,'' Accessing cgns database... (must'',
     +  '' be linked to CGNS Version 2.5.2 or later)'')')
      write(901,'('' Opening CGNS database from cgnstools/ropencgns'')')
      call cg_open_f(cgnsdesired,CG_MODE_READ,iccg,ier)
      if (ier .ne. 0) then
        write(901,'('' Error opening cgns file'')')
        write(901,'('' Make sure the file on 1st line of input IS'',
     +   '' a cgns file!'')')
        call cg_error_exit_f
      end if
      write(901,'(/,'' Opening cgns database '',a32)') cgnsdesired
      write(901,'(''   Note: if restarting from a cgns file created'',
     +  '' prior to early 2002, then'')')
      write(901,'(''   will not work if the old case was (1) second'',
     +  '' order time accurate,'')')
      write(901,'(''   (2) moving or deforming grid, or (3) ialph>0.'',
     +  ''  In all of these cases,'')')
      write(901,'(''   changes have been made in how variables are'',
     +  '' stored in cgns.'',/)')
c
c   Find out how many "bases" there are
      call cg_nbases_f(iccg, nbases, ier)
      if (ier .ne. 0) call cg_error_exit_f
      write(901,'('' Found ''i5,'' base(s) in cgns file: '',a32)') 
     +  nbases,cgnsdesired
c
c   Get base name and index dimension (dimensionality)
      ifind=0
      do ib=1,nbases
        call cg_base_read_f(iccg, ib, basename, indx1, indx2, ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (indx1 .ne. indx2) then
          write(901,'('' Error, celldim MUST be same as physdim'',
     +     '' for CFL3D to use this file'')')
          stop
        end if
        if (basename .eq. basedesired) then
          ibase=ib
          index_dim=indx2
          ifind=1
        end if
      end do
      if (ifind .eq. 0) then
        write(901,'('' Desired base name '',a32,'' not found'')')
     +   basedesired
        write(901,'('' Available bases found are:'')')
        do ib=1,nbases
          call cg_base_read_f(iccg, ib, basename, indx1, index_dim, ier)
          if (ier .ne. 0) call cg_error_exit_f
          write(901,'(a32)') basename
        enddo
        stop
      end if
c
      if (index_dim .ne. idimdesired) then
        write(901,'('' Index dimension wrong in base'',a32)')
     +   basedesired
        write(901,'('' index_dim='',i5,''.  Wanted = '',i5)')
     +   index_dim,idimdesired
        stop
      end if
c
c   Get number of zones
      call cg_nzones_f(iccg, ibase, nzones, ier)
      if (ier .ne. 0) call cg_error_exit_f
      write(901,'('' Found ''i5,'' zone(s) in basename: '',a32)')nzones,
     +  basedesired
c
      return
      end
c
      subroutine wopencgns(cgnsdesired,basedesired,idimdesired,iccg,
     .                     ibase,nzones)
c***********************************************************************
c     Purpose: Opens CGNS data file for MODIFY (WRITE as well as READ)
c     If cannot find base of name "basedesired" or if the dimensionality
c     is not "idimdesired", the routine stops execution.
c     Returns CGNS file index, base index, and number of zones
c
c   INPUT:
c     cgnsdesired.....name of CGNS file to open (character*80)
c     basedesired.....name of base desired (character*32)
c     idimdesired.....dimensionality desired (integer)
c     
c   OUTPUT:
c     iccg............index number of CGNS file (integer)
c     ibase...........index number of base desired (integer)
c     nzones..........number of zones in base desired (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      character*80 cgnsdesired
      character*32 basedesired,basename
c
c   Open CGNS file
      write(901,'('' Opening CGNS database from cgnstools/wopencgns'')')
      call cg_open_f(cgnsdesired,CG_MODE_MODIFY,iccg,ier)
      if (ier .ne. 0) then
        write(901,'('' Error opening cgns file'')')
        write(901,'('' Make sure the file on 1st line of input IS'',
     +   '' a cgns file!'')')
        call cg_error_exit_f
      end if
c
c   Find out how many "bases" there are
      call cg_nbases_f(iccg, nbases, ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get base name and index dimension (dimensionality)
      ifind=0
      do ib=1,nbases
        call cg_base_read_f(iccg, ib, basename, indx1, indx2, ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (indx1 .ne. indx2) then
          write(901,'('' Error, celldim MUST be same as physdim'',
     +     '' for CFL3D to use this file'')')
          stop
        end if
        if (basename .eq. basedesired) then
          ibase=ib
          index_dim=indx2
          ifind=1
        end if
      end do
      if (ifind .eq. 0) then
        write(901,'('' Desired base name '',a32,'' not found'')')
     +   basedesired
        write(901,'('' Available bases found are:'')')
        do ib=1,nbases
          call cg_base_read_f(iccg, ib, basename, indx1, index_dim, ier)
          if (ier .ne. 0) call cg_error_exit_f
          write(901,'(a32)') basename
        enddo
        stop
      end if
c
      if (index_dim .ne. idimdesired) then
        write(901,'('' Index dimension wrong in base'',a32)')
     +   basedesired
        write(901,'('' index_dim='',i5,''.  Wanted = '',i5)')
     +   index_dim,idimdesired
        stop
      end if
c
c   Get number of zones
      call cg_nzones_f(iccg, ibase, nzones, ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
c
      subroutine writesoln(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +  jdim,kdim,wk,q,qj0,qk0,qi0,bcj,bck,bci,i2d,ialph,nsoluse)
c***********************************************************************
c     Purpose: Writes solution to CGNS data base for zone number "igrid"
c     The CGNS file must be already opened in "MODE_MODIFY" mode.
c     Writes results in the following order:
c
c     (idim,jdim,kdim) <- necessary for CGNS!,
c
c     even though they are stored in CFL3D in the order: jdim,kdim,idim
c     Also writes rind cell information:  all FACES of an IxJxK grid 
c     are accurate, but all edges and corners have no value... 
c     currently we use nearby-cell information (somewhat arbitrary).
c     For 2-D solutions, rind cells are not written in the i-direction.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid"), at finest level (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) - if actual dimensions are a
c                       subset idima,jdima,kdima, results will be
c                       prolonged to the fine grid and written at that
c                       level (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      q................cell-centered q-values, in (j,k,i) order,
c                       CONSERVED quantities (real)
c      qi0,qj0,qk0......B.C. q-values from CFL3D, 
c                       CONSERVED quantities (real)
c      bci,bcj,bck......flag from CFL3D; 0=ghost-cell B.C., 1=B.C. at
c                       interface (real)
c      i2d..............0 if 3-D (rind data to be output in all 3 
c                          directions),
c                       1 if 2-D (rind data to be output only in j- 
c                          and k- directions) (integer)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c   OUTPUTS:
c      nsoluse..........CGNS "FlowSolution" index number (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension q(jdim,kdim,idim,5)
      dimension qi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4),
     +          qk0(jdim,idim-1,5,4),bci(jdim,kdim,2),
     +          bcj(kdim,idim-1,2),bck(jdim,idim-1,2),
     +          irinddata(6)
      dimension wk(i2d:idima-i2d,0:jdima,0:kdima)
c
      character*32 solname
      character*80 string1
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Write FlowSolution node
      write(901,'('' cgns writing to <<FlowSolution>> node'')')
      write(901,'(''  ...any flow solution node of any other name'',
     + '' is left alone'')')
      write(901,'(''  ...if a FlowSolution node already exists, it'',
     + '' is overwritten'')')
      solname='FlowSolution'
      call cg_sol_write_f(iccg,ibase,igrid,solname,CellCenter,
     +                    nsoluse,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Write q
c   First, get "factor" in case want to write every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
        string1='Fine grid solution (on same level as GridCoordinates)'
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
        string1='Every-2nd-cell solution (1 level down)'
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
        string1='Every-4th-cell solution (2 levels down)'
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
        string1='Every-8th-cell solution (3 levels down)'
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
        string1='Every-16th-cell solution (4 levels down)'
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'FlowSolution_t',
     + nsoluse,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_descriptor_write_f('Information',string1,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Write rind information
      do n=1,6
        irinddata(n)=1
      enddo
      if (i2d .eq. 1) then
        irinddata(1)=0
        irinddata(2)=0
      end if
      call cg_rind_write_f(irinddata,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Density:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,1,0,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble, 
     +                      'Density', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'Density', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   X-momentum:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,2,0,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'MomentumX', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'MomentumX', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Y-momentum:
      if (ialph .eq. 0) then
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,3,0,wk)
      else
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,4,0,wk)
      end if
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'MomentumY', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'MomentumY', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Z-momentum:
      if (ialph .eq. 0) then
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,4,0,wk)
      else
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,3,0,wk)
      do i=i2d,idima-i2d
        do j=0,jdima
          do k=0,kdima
            wk(i,j,k)=-wk(i,j,k)
          enddo
        enddo
      enddo
      end if
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'MomentumZ', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'MomentumZ', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Stagnation energy per unit volume:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,5,5,0,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'EnergyStagnationDensity', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'EnergyStagnationDensity', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
c
      subroutine readsoln(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +  jdim,kdim,wk,nsoluse,irind,jrind,krind,i2d,ialph,q,iprim)
c***********************************************************************
c     Purpose: Reads solution from CGNS data base for zone number
c     "igrid". The CGNS file must be already opened. If the desired 
c     solution names are not in the CGNS file, the routine stops
c     execution. Must know (through irind,jrind,krind) what Rind 
c     Cell info exists. Gets results in the following order:
c
c      (jdim,kdim,idim) <- necessary for CFL3D!,
c
c     and assumes they are stored in the CGNS database in the order:
c     (idim,jdim,kdim)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idima,jdima,kdima...expected (GRIDPOINT) dimensions of this zone
c                          (zone "igrid"), at finest level (integers)
c      idim,jdim,kdim...desired (GRIDPOINT) dimensions of this zone
c                       (e.g., could be every other point of fine
c                       grid) (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      nsoluse..........CGNS "FlowSolution" index number (determined
c                       outside this routine) (integer)
c      irind,jrind,krind...= 0 if rind data exists for the direction,
c                          = 1 if it does not (integers)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c   OUTPUTS:
c      q................cell-centered q-values, in (j,k,i) order,
c                       CONSERVED or PRIMITIVE quantities (real)
c      iprim............0 if conserved read, 1 if primitive (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter(numfield=50)
c
      dimension q(jdim,kdim,idim,5)
      dimension wk(idima-2*irind+1,jdima-2*jrind+1,kdima-2*krind+1)
      dimension irmin(3),irmax(3)
c
      character*32 solname,name
      character*32 fieldname(numfield),fndesired(numfield)
      character*80 text
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Find out how many fields there are and their names:
      call cg_nfields_f(iccg, ibase, igrid, nsoluse, nfields, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nfields .gt. numfield) then
        write(901,'('' Need to increase numfield parameter'',
     +    '' in readsoln'')')
        stop
      end if
      do n=1,nfields
      call cg_field_info_f(iccg, ibase, igrid, nsoluse, n, itype,
     +  fieldname(n), ier)
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Check to make sure the conserved field names needed exist:
      numwanted=5
      fndesired(1)='Density'
      fndesired(2)='MomentumX'
      fndesired(3)='MomentumY'
      fndesired(4)='MomentumZ'
      fndesired(5)='EnergyStagnationDensity'
      iprim=0
      do n=1,numwanted
        ido=0
        do nn=1,nfields
          if(fndesired(n) .eq. fieldname(nn)) then
            ido=1
          end if
        enddo
        if (ido .eq. 0) then
          write(901,'('' Cannot find solution of name '',a32)') 
     +      fndesired(n)
          write(901,'(''  for conservative... try primitive instead'')')
          iprim=1
        end if
      enddo
      if (iprim .eq. 1) then
c       if conserved do not exist, try reading primitive variables
      fndesired(1)='Density'
      fndesired(2)='VelocityX'
      fndesired(3)='VelocityY'
      fndesired(4)='VelocityZ'
      fndesired(5)='Pressure'
      do n=1,numwanted
        ido=0
        do nn=1,nfields
          if(fndesired(n) .eq. fieldname(nn)) then
            ido=1
          end if
        enddo
        if (ido .eq. 0) then
          write(901,'('' Cannot find solution of name '',a32)')
     +      fndesired(n)
          write(901,'(''   for primitive either...'')')
          iprim=2
        end if
      enddo
      end if
      if (iprim .eq. 2) then
c       if cannot find conserved or primitive names, stop execution
        stop
      end if
c
c   Read Q's
c   First, get "factor" in case want to read every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
c   Check that correct grid level is stored
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'FlowSolution_t',
     + nsoluse,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_ndescriptors_f(ndesc,ier)
      if (ier .ne. 0) call cg_error_exit_f
      istop=0
      do n=1,ndesc
        call cg_descriptor_read_f(n,name,text,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (name .eq. 'Information') then
          if (text .eq. 
     +     'Fine grid solution (on same level as GridCoordinates)' .and.
     +     nfac .ne. 1) istop=1
          if (text .eq.
     +     'Every-2nd-cell solution (1 level down)' .and.
     +     nfac .ne. 2) istop=2
          if (text .eq.
     +     'Every-4th-cell solution (2 levels down)' .and.
     +     nfac .ne. 4) istop=3
          if (text .eq.
     +     'Every-8th-cell solution (3 levels down)' .and.
     +     nfac .ne. 8) istop=4
          if (text .eq.
     +     'Every-16th-cell solution (3 levels down)' .and.
     +     nfac .ne. 8) istop=5
          goto 333
        end if
      enddo
      if (nfac .eq. 1) then
      write(901,'('' Information node not found.  Defaulting to read'',
     + '' of fine level flowfield'')')
      else if (nfac .eq. 2) then
      write(901,'('' Information node not found.  Defaulting to read'',
     + '' of every 2nd pt of fine level flowfield'')')
      else if (nfac .eq. 4) then
      write(901,'('' Information node not found.  Defaulting to read'',
     + '' of every 4th pt of fine level flowfield'')')
      else
      write(901,'('' Information node not found.  Defaulting to read'',
     + '' of every 8th pt of fine level flowfield'')')
      end if
      goto 334
 333  continue
      if (istop .gt. 0) then
        if (nfac .eq. 1) then
        write(901,'('' Expecting read on level '',i5,'' of grid...'',
     +   '' you are asking for level 1'')') istop
        else if (nfac .eq. 2) then
        write(901,'('' Expecting read on level '',i5,'' of grid...'',
     +   '' you are asking for level 2'')') istop
        else if (nfac .eq. 4) then
        write(901,'('' Expecting read on level '',i5,'' of grid...'',
     +   '' you are asking for level 3'')') istop
        else
        write(901,'('' Expecting read on level '',i5,'' of grid...'',
     +   '' you are asking for level 4'')') istop
        end if
        stop
      end if
 334  continue
c   Set up array bounds:
      irmin(1)=1
      irmin(2)=1
      irmin(3)=1
      irmax(1)=idima-2*irind+1
      irmax(2)=jdima-2*jrind+1
      irmax(3)=kdima-2*krind+1
c   read Q's:
      do m=1,5
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, fndesired(m),
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, fndesired(m),
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (m .eq. 3 .and. ialph .ne. 0) then
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,q(1,1,1,4))
      else if (m .eq. 4 .and. ialph .ne. 0) then
      do i=1,idima-2*irind+1
        do j=1,jdima-2*jrind+1
          do k=1,kdima-2*krind+1
            wk(i,j,k)=-wk(i,j,k)
          enddo
        enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,q(1,1,1,3))
      else
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,q(1,1,1,m))
      end if
      enddo
c
      return
      end
c
      subroutine writeblnk(iccg,ibase,igrid,idim,jdim,kdim,npnts,ipnt,
     .                     blank)
c***********************************************************************
c     Purpose: Writes OversetHoles info (at cell centers) to CGNS data 
c     base for zone number "igrid" (does this only if the zone in 
c     question has overset information)
c     The CGNS file must be already opened in "MODE_MODIFY" mode.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (integers)
c      npnts............number of points where blank=0 in this zone
c                       (integer)
c      ipnt.............input as work array, used internal to this 
c                       routine for storing OversetHoles info,
c                       must be at least of size (3*npnts) in the
c                       call to this routine (integer)
c      blank............cell-centered blank-values, in (j,k,i) order 
c                       0.0=blanked, 1.0=not-blanked (real)
c   OUTPUTS:
c      none
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension blank(jdim,kdim,idim)
      dimension ipnt(3,npnts)
c
c   Create integer array where there are holes:
      n=0
      do k=1,kdim-1
      do j=1,jdim-1
      do i=1,idim-1
        if(blank(j,k,i) .eq. 0.) then
          n=n+1
          ipnt(1,n)=i
          ipnt(2,n)=j
          ipnt(3,n)=k
        end if
      enddo
      enddo
      enddo
      npnts=n
c   Write OversetHoles
      if (npnts .gt. 0) then
        call cg_hole_write_f(iccg, ibase, igrid, 'OversetHoles',
     +                       CellCenter, PointList, 1, npnts, ipnt,
     +                       iindex, ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
c
      return
      end
c
      subroutine writeturb(iccg,ibase,igrid,nsoluse,idima,jdima,kdima,
     +   idim,jdim,kdim,ivisc,wk,vist3d,tursav,smin,xjb,xkb,
     +   xib,blnum,cmuv,vj0,vk0,vi0,tj0,tk0,ti0,i2d,
     +   xmach,reue,nummem)
c***********************************************************************
c     Purpose: Writes turb info to CGNS data base for zone number 
c     "igrid". The CGNS file must be already opened in "MODE_MODIFY"
c     mode. Writes results in the following order:
c
c       (idim,jdim,kdim) <- necessary for CGNS!,
c
c     even though they are stored in CFL3D in the order: jdim,kdim,idim
c     Also writes rind cell information for vist3d and tursav:  all 
c     FACES of an IxJxK grid are accurate, but all edges and corners
c     have no value...currently we use nearby-cell information (somewhat
c     arbitrary). No B.C info exists for other output quantities (like
c     smin, xjb, etc). Currently we output nearby values for these.
c     For 2-D solutions, rind cells are not written in the i-direction.
c     Note: turb quantities need to be renormalized to be put into SIDS-
c     standard "consistent" normalization
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      nsoluse..........CGNS "FlowSolution" index number (determined
c                       outside this routine) (integer)
c      idima,jdima,kdima...expected dimensions of this zone (zone
c                          "igrid"), at finest level (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) - if actual dimensions are a subset
c                       of idima,jdima,kdima, results will be prolonged 
c                       to the fine grid and written at that level
c                       (integers)
c      ivisc............ivisc parameters, of dimension 3 (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      vist3d...........cell-centered eddy viscosity-values, in (j,k,i) 
c                       order (real) (CFL3D-normalization)
c      tursav...........set of nummem cell-centered turbulence quantities, in 
c                       (j,k,i,nummem) order (real) (CFL3D-normalization)
c      smin.............minimum distance function, in (j,k,i) order 
c                       (real)
c      xjb..............nearest-wall j-location for B-B/LES, in (j,k,i) 
c                       order (real)
c      xkb..............nearest-wall k-location for B-B/LES, in (j,k,i) 
c                       order (real)
c      xib..............nearest-wall i-location for B-B/LES, in (j,k,i) 
c                       order (real)
c      blnum............nearest-wall block number for B-B/LES, in (j,k,i) 
c                       order (real)
c      cmuv.............variable cmu coefficient for EASM, var g, in
c                       (j,k,i) order (real)
c      vi0,vj0,vk0......B.C. vist3d-values (eddy viscosity) from CFL3D 
c                       (real) (CFL3D-normalization)
c      ti0,tj0,tk0......B.C. turbulence-values from CFL3D (real) 
c                       (CFL3D-normalization)
c      i2d..............0 if 3-D (rind data to be output in all 3 
c                       directions), 1 if 2-D (rind data to be output 
c                       only in j- and k- directions) (integer)
c      xmach............Mach number (real)
c      reue.............Reynolds number (real)
c      nummem...........one of dimensions for tursav, ti0,tj0,tk0 (integer)
c   OUTPUTS:
c      none
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension vist3d(jdim,kdim,idim),tursav(jdim,kdim,idim,nummem),
     .          smin(jdim-1,kdim-1,idim-1),xjb(jdim-1,kdim-1,idim-1),
     .          xkb(jdim-1,kdim-1,idim-1),xib(jdim,kdim,idim),
     .          blnum(jdim-1,kdim-1,idim-1),cmuv(jdim-1,kdim-1,idim-1)
      dimension wk(i2d:idima-i2d,0:jdima,0:kdima)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
      dimension ivisc(3)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Write turb quantities
c   First, get "factor" in case want to write every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
c
      if (ivisc(1).ge.2 .or. ivisc(2).ge.2 .or. ivisc(3).ge.2) then
         write(901,'(''   writing vist3d data to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Eddy viscosity:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     vist3d,wk,wk,wk,vi0,vj0,vk0,1,1,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'Mu_tOverMu_inf', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'Mu_tOverMu_inf', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*xmach/reue
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble, 
     +                      'ViscosityEddy', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'ViscosityEddy', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
         write(901,'(''   writing field eqn turb quantity'',
     .              '' min dist to cgns file'')')
c        Min distance:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     smin,wk,wk,wk,wk,wk,wk,1,1,2,wk)
c        Need to be sure to write out absolute value of smin
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=abs(wk(i,j,k))
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentDistance', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentDistance', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.4 .or. ivisc(2).eq.4 .or. ivisc(3).eq.4) then
         write(901,'(''   writing B-B field eqn turb info'',
     .              '' to cgns file'')')
c        BBReynolds:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*xmach/reue
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentBBReynolds', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentBBReynolds', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        BBNearestWallJ:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     xjb,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentBBNearestWallJ', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentBBNearestWallJ', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        BBNearestWallK:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     xkb,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentBBNearestWallK', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentBBNearestWallK', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        BBNearestWallI:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     xib,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentBBNearestWallI', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentBBNearestWallI', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        BBNearestWallBlk:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     blnum,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentBBNearestWallBlk', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentBBNearestWallBlk', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.5 .or. ivisc(2).eq.5 .or. ivisc(3).eq.5) then
         write(901,'(''   writing S-A field eqn turb info'',
     .              '' to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        SANuTilde:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*xmach/reue
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentSANuTilde', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentSANuTilde', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.6 .or. ivisc(2).eq.6 .or. ivisc(3).eq.6 .or.
     .    ivisc(1).eq.7 .or. ivisc(2).eq.7 .or. ivisc(3).eq.7 .or.
     .    ivisc(1).eq.8 .or. ivisc(2).eq.8 .or. ivisc(3).eq.8 .or.
     .    ivisc(1).eq.12.or. ivisc(2).eq.12.or. ivisc(3).eq.12.or.
     .    ivisc(1).eq.14.or. ivisc(2).eq.14.or. ivisc(3).eq.14) then
         write(901,'(''   writing k-omega field eqn turb info'',
     .              '' to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Omega:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*reue/xmach
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        K:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,2,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.9 .or. ivisc(2).eq.9 .or. ivisc(3).eq.9 .or.
     .    ivisc(1).eq.10.or. ivisc(2).eq.10.or. ivisc(3).eq.10.or.
     .    ivisc(1).eq.11.or. ivisc(2).eq.11.or. ivisc(3).eq.11.or.
     .    ivisc(1).eq.13.or. ivisc(2).eq.13.or. ivisc(3).eq.13) then
         write(901,'(''   writing k-epsilon field eqn turb info'',
     .              '' to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Epsilon:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*reue/xmach
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentDissipation', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentDissipation', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        K:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,2,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.15.or. ivisc(2).eq.15.or. ivisc(3).eq.15) then
         write(901,'(''   writing k-enstrophy field eqn turb info'',
     .              '' to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Enstrophy:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*reue/xmach
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnstrophy', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnstrophy', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        K:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,2,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.25.or. ivisc(2).eq.25.or. ivisc(3).eq.25) then
         write(901,'(''   writing LES-related distance info'',
     .              '' to cgns file'')')
c        LESNearestWallJ:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     xjb,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentLESNearestWallJ', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentLESNearestWallJ', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        LESNearestWallK:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     xkb,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentLESNearestWallK', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentLESNearestWallK', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        LESNearestWallI:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     xib,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentLESNearestWallI', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentLESNearestWallI', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        LESNearestWallBlk:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     blnum,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentLESNearestWallBlk', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentLESNearestWallBlk', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.30.or. ivisc(2).eq.30.or. ivisc(3).eq.30) then
         write(901,'(''   writing k-omega field eqn turb info'',
     .              '' plus transition to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Omega:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*reue/xmach
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        K:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,2,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         end if
c        Intermittency:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,3,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'Intermittency', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'Intermittency', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.40.or. ivisc(2).eq.40.or. ivisc(3).eq.40) then
         write(901,'(''   writing k-omega field eqn turb info'',
     .              '' plus intermittency and Retheta_t to cgns file'',
     +     '' (SIDS-consistent-normalization)'')')
c        Omega:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,1,1,wk)
c        put in SIDS-like-"consistent normalization"
         do i=i2d,idima-i2d
         do j=0,jdima
         do k=0,kdima
           wk(i,j,k)=wk(i,j,k)*reue/xmach
         enddo
         enddo
         enddo
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentDissipationRate', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
c        K:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,2,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentEnergyKinetic', wk, ifindex, ier)
         end if
c        Intermittency:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,3,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'Intermittency', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'Intermittency', wk, ifindex, ier)
         end if
c        Retheta_t:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +     tursav,wk,wk,wk,ti0,tj0,tk0,nummem,4,1,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'Retheta_t', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'Retheta_t', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      if (ivisc(1).eq.8 .or. ivisc(2).eq.8 .or. ivisc(3).eq.8 .or.
     .    ivisc(1).eq.9 .or. ivisc(2).eq.9 .or. ivisc(3).eq.9 .or.
     .    ivisc(1).eq.13.or. ivisc(2).eq.13.or. ivisc(3).eq.13.or.
     .    ivisc(1).eq.14.or. ivisc(2).eq.14.or. ivisc(3).eq.14) then
         write(901,'(''   writing EASM variable cmu info'',
     .              '' to cgns file'')')
c        Cmu:
         call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim-1,kdim-1,
     +     cmuv,wk,wk,wk,wk,wk,wk,1,1,2,wk)
         if (idouble .eq. 1) then
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                   'TurbulentCmu', wk, ifindex, ier)
         else
         call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                   'TurbulentCmu', wk, ifindex, ier)
         end if
         if (ier .ne. 0) call cg_error_exit_f
      end if
c
      return
      end
c
      subroutine readturb(iccg,ibase,igrid,nsoluse,idima,jdima,kdima,
     .  idim,jdim,kdim,wk,iv1,iv2,iv3,iread,irind,jrind,krind,i2d,
     .  ivmx,xmach,reue,vist3d,tursav,smin,xjb,xkb,xib,
     .  blnum,cmuv,nummem)
c***********************************************************************
c     Purpose: Reads turb quantities from CGNS data base for zone number
c     "igrid". The CGNS file must be already opened. If the desired
c     solution names are not in the CGNS file, the routine stops 
c     execution. Must know (through irind,jrind,krind) what Rind 
c     Cell info exists. Gets results in the following order:
c
c       (jdim,kdim,idim) <- necessary for CFL3D!,
c
c     and assumes they are stored in the CGNS database in the order:
c     (idim,jdim,kdim)
c     Note: turb quantities need to be renormalized to be taken out of 
c     SIDS-standard "consistent" normalization, and back to 
c     CFL3D-type-normalization
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      nsoluse..........CGNS "FlowSolution" index number (determined
c                       outside this routine) (integer)
c      idima,jdima,kdima...expected (GRIDPOINT) dimensions of this
c                          zone (zone "igrid"), at finest level
c                          (integers)
c      idim,jdim,kdim...desired (GRIDPOINT) dimensions of this zone
c                       (e.g., could be every other point of fine
c                       grid) (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      iv1,iv2,iv3......ivisc(n) parameters from previous run (integers)
c      iread............flag that indicates whether to read all info
c                       or not:
c                       =0 turb model has changed, do NOT read all info
c                       =1 turb model same (or has remained a 2-eqn
c                          model), read info
c                       This is important, because some memory may not
c                       be set up
c      irind,jrind,krind...= 0 if rind data exists for the direction,
c                          = 1 if it does not (integers)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      ivmx.............turb model number for CFL3D (integer)
c      xmach............Mach number given in CGNS file (real)
c      reue.............Reynolds number given in CGNS file (real)
c      nummem...........one of dimensions of tursav (integer)
c   OUTPUTS:
c      vist3d...........cell-centered eddy viscosity-values, in (j,k,i)
c                       order (real)
c      tursav...........set of nummem cell-centered turbulence quantities,
c                       in (j,k,i,nummem) order (real)
c      smin.............minimum distance function, in (j,k,i)
c                       order (real)
c      xjb..............nearest-wall j-location for B-B/LES, in (j,k,i)
c                       order (real)
c      xkb..............nearest-wall k-location for B-B/LES, in (j,k,i)
c                       order (real)
c      xib..............nearest-wall i-location for B-B/LES, in (j,k,i)
c                       order (real)
c      blnum............nearest-wall block number for B-B/LES, in (j,k,i)
c                       order (real)
c      cmuv.............variable cmu coefficient for EASM, var g, in
c                       (j,k,i) order (real)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension vist3d(jdim,kdim,idim),tursav(jdim,kdim,idim,nummem),
     .          smin(jdim-1,kdim-1,idim-1),xjb(jdim-1,kdim-1,idim-1),
     .          xkb(jdim-1,kdim-1,idim-1),xib(jdim,kdim,idim),
     .          blnum(jdim-1,kdim-1,idim-1),cmuv(jdim-1,kdim-1,idim-1)
      dimension wk(idima-2*irind+1,jdima-2*jrind+1,kdima-2*krind+1)
      dimension irmin(3),irmax(3)
c
      character*32 solname
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
      ivmxold = max(iv1,iv2)
      ivmxold = max(ivmxold,iv3)
c
c   First, get "factor" in case want to read every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
c   Set up array bounds:
      irmin(1)=1
      irmin(2)=1
      irmin(3)=1
      irmax(1)=idima-2*irind+1
      irmax(2)=jdima-2*jrind+1
      irmax(3)=kdima-2*krind+1
c
c   Read vist3d
      if (ivmxold .ge. 2 .and. iread .eq. 1) then
      solname='ViscosityEddy'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname, 
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname, 
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*reue/xmach
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,vist3d)
      end if
c
c   Read smin
      if (ivmxold .ge. 4) then
      solname='TurbulentDistance'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,smin)
      end if
c
c   Read B-B info
      if (ivmxold .eq. 4 .and. iread .eq. 1) then
      solname='TurbulentBBReynolds'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*reue/xmach
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentBBNearestWallJ'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,         
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xjb)
c
      solname='TurbulentBBNearestWallK'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xkb)
c
      solname='TurbulentBBNearestWallI'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xib)
c
      solname='TurbulentBBNearestWallBlk'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,blnum)
      end if
c
c   Read S-A info
      if (ivmxold .eq. 5 .and. iread .eq. 1) then
      solname='TurbulentSANuTilde'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*reue/xmach
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
      end if
c
c   Read k-omega info
      if ((ivmxold .eq. 6 .or. ivmxold .eq. 7 .or.
     +     ivmxold .eq. 8 .or. ivmxold .eq. 12.or.
     +     ivmxold .eq. 14) .and. iread .eq. 1) then
      solname='TurbulentDissipationRate'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*xmach/reue
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentEnergyKinetic'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,2))
      end if
c
c   Read k-epsilon info
      if ((ivmxold .eq. 9 .or. ivmxold .eq. 10 .or.
     +     ivmxold .eq.11 .or. ivmxold .eq. 13) .and.
     +     iread .eq. 1) then
      solname='TurbulentDissipation'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*xmach/reue
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentEnergyKinetic'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,2))
      end if
c
c   Read k-enstrophy info
      if (ivmxold .eq.15 .and. iread .eq. 1) then
      solname='TurbulentEnstrophy'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*xmach/reue
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentEnergyKinetic'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,2))
      end if
c
c   Read LES-related info
      if (ivmxold .eq. 25 .and. iread .eq. 1) then
      solname='TurbulentLESNearestWallJ'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xjb)
c
      solname='TurbulentLESNearestWallK'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xkb)
c
      solname='TurbulentLESNearestWallI'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,xib)
c
      solname='TurbulentLESNearestWallBlk'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,blnum)
      end if
c
c   Read k-omega+transition info
      if (ivmxold .eq. 30 .and. iread .eq. 1) then
      solname='TurbulentDissipationRate'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*xmach/reue
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentEnergyKinetic'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,2))
c
      solname='Intermittency'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,3))
      end if
c
c   Read k-omega+intermittency+Retheta_t info
      if (ivmxold .eq. 40 .and. iread .eq. 1) then
      solname='TurbulentDissipationRate'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Change from SIDS-"consistent-normalization" back to CFL3D-type
      do i=1,idima-2*irind+1
      do j=1,jdima-2*jrind+1
      do k=1,kdima-2*krind+1
        wk(i,j,k)=wk(i,j,k)*xmach/reue
      enddo
      enddo
      enddo
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,1))
c
      solname='TurbulentEnergyKinetic'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,2))
c
      solname='Intermittency'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,3))
c
      solname='Retheta_t'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav(1,1,1,4))
      end if
c
c   Read variable cmuv
      if ((ivmxold .eq. 8 .or. ivmxold .eq. 9 .or.
     +     ivmxold .eq.13 .or. ivmxold .eq.14) .and.
     +    (ivmx    .eq. 8 .or. ivmx    .eq. 9 .or.
     +     ivmxold .eq.13 .or. ivmxold .eq.14) .and.
     +     iread .eq. 1) then
      solname='TurbulentCmu'
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse, solname,
     +                     RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderq(idima,jdima,kdima,idim-1,jdim-1,kdim-1,nfac,
     +    nfaci,wk,irind,jrind,krind,cmuv)
      end if
c
      return
      end
c
      subroutine read1to1(iccg,ibase,nzones,n1to1_global,iblk,
     +   irange,index_vary)
c***********************************************************************
c     Purpose: Gets 1-to-1 connectivity information.  The CGNS file must
c     already be opened.
c     Note: in current release of CFL3D, this routine is not used
c     (1-to-1 info is still obtained from standard input file)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      nzones...........number of zones in the CGNS file (integer)
c      n1to1_global.....number of 1-to-1 interfaces in base ibase
c                       (integer)
c                       NOTE:  this is 1/2 the number of 1-to-1 NODES
c   OUTPUTS:
c      iblk.............block number of receiver patch,
c                       dimensioned (2,n1to1_global)
c                       1st index is receiver, 2nd is donor (integer)
c      irange...........min-max range defining receiver patch, 
c                       dimensioned (2,6,n1to1_global) 
c                       1st index is receiver, 2nd is donor (integer)
c      index_vary.......which indices vary on receiver patch, 
c                       dimensioned (2,2,n1to1_global) 
c                       1st index is receiver, 2nd is donor (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension iblk(2,n1to1_global),idrang(6,n1to1_global),
     +          irange(2,6,n1to1_global),irang(6,n1to1_global),
     +          index_vary(2,2,n1to1_global),
     +          itransform(3,n1to1_global),isize(9)
c
      character*32 connectname(n1to1_global),zonename(n1to1_global),
     +             donorname(n1to1_global)
      character*32 znam(nzones)
c
        ier=0
        call cg_1to1_read_global_f(iccg,ibase,connectname,zonename,
     +   donorname,irang,idrang,itransform,ier)
        if (ier .ne. 0) call cg_error_exit_f
c   find out all zone names
        do n=1,nzones
          call cg_zone_read_f(iccg,ibase,n,znam(n),isize,ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
c   assign 1-to-1 information appropriately
        nbli=n1to1_global
        do n=1,nbli
          iflg1=0
          iflg2=0
          do nm=1,nzones
            if (zonename(n) .eq. znam(nm)) then
              iblk(1,n)=nm
              iflg1=1
            end if
            if (donorname(n) .eq. znam(nm)) then
              iblk(2,n)=nm
              iflg2=1
            end if
          enddo
          if (iflg1 .eq. 0 .or. iflg2 .eq. 0) then
            write(901,'('' Error, cgns 1-to-1 zonename not'',
     +        '' matched'')')
            stop
          end if
          do nm=1,6
            irange(1,nm,n)=irang(nm,n)
            irange(2,nm,n)=idrang(nm,n)
          enddo
c   find out which index not varying in 1st zone
          if (irang(1,n) .ne. irang(4,n)) then
            index_vary(1,1,n)=1
            if (irang(2,n) .ne. irang(5,n)) then
              index_vary(1,2,n)=2
            else
              index_vary(1,2,n)=3
            end if
          else
            index_vary(1,1,n)=2
            index_vary(1,2,n)=3
          end if
c   use itransform to get appropriate index of 2nd zone
          index_vary(2,1,n)=abs(itransform(index_vary(1,1,n),n))
          index_vary(2,2,n)=abs(itransform(index_vary(1,2,n),n))
c
c   check to make sure ranges match, and direction of indices (lo to hi vs. hi to lo)
c   is as expected, according to value of itransform
c   
c   first, get values of hi-lo
          ival1rec=irange(1,index_vary(1,1,n)+3,n)-
     +            irange(1,index_vary(1,1,n),n)
          ival2rec=irange(1,index_vary(1,2,n)+3,n)-
     +            irange(1,index_vary(1,2,n),n)
          ival1don=irange(2,index_vary(2,1,n)+3,n)-
     +            irange(2,index_vary(2,1,n),n)
          ival2don=irange(2,index_vary(2,2,n)+3,n)-
     +            irange(2,index_vary(2,2,n),n)
c   the absolute value of the differences should match
          if ((abs(ival1rec) .ne. abs(ival1don)) .or. 
     +        (abs(ival2rec) .ne. abs(ival2don))) then
            write(901,'('' 1-to-1 donor & reciever ranges do not'',
     +       '' match, interface #'',i5)') n
            write(901,'(''    irange1='',6i6)') irange(1,1,n),
     +       irange(1,2,n),irange(1,3,n),irange(1,4,n),
     +       irange(1,5,n),irange(1,6,n)
            write(901,'(''    irange2='',6i6)') irange(2,1,n),
     +       irange(2,2,n),irange(2,3,n),irange(2,4,n),
     +       irange(2,5,n),irange(2,6,n)
            stop
          end if
c   check the directionality
          icheck1=ival1rec*ival1don*itransform(index_vary(1,1,n),n)
          icheck2=ival2rec*ival2don*itransform(index_vary(1,2,n),n)
          if (icheck1 .lt. 0 .or. icheck2 .lt. 0) then
            write(901,'('' directionality of transform wrong'',
     +       '', interface #'',i5)') n
            write(901,'(''    irange1='',6i6)') irange(1,1,n),
     +       irange(1,2,n),irange(1,3,n),irange(1,4,n),
     +       irange(1,5,n),irange(1,6,n)
            write(901,'(''    irange2='',6i6)') irange(2,1,n),
     +       irange(2,2,n),irange(2,3,n),irange(2,4,n),
     +       irange(2,5,n),irange(2,6,n)
            write(901,'(''    itransform='',3i5)') itransform(1,n),
     +       itransform(2,n),itransform(3,n)
            stop
          end if
        enddo
        return
        end
c
      subroutine getnumseg(iccg,ibase,izone,numseg)
c***********************************************************************
c     Purpose: Determine the number of BC segments for each of 6 faces
c     of a particular zone.
c     For now, this routine only counts BC segments and 1-to-1 segments.
c     Note: in current release of CFL3D, this routine is not used
c     (BC info is still obtained from standard input file)
c
c   INPUT:
c     iccg.............CGNS file index number (determined outside this
c                      routine) (integer)
c     ibase............CGNS base index number (determined outside this
c                      routine) (integer)
c     izone............zone (or grid) number (integer)
c   OUTPUT:
c     numseg....number of segments on each face (1-6 represent 
c               ilo,jlo,klo,ihi,jhi,khi, respectively)
c               dimensioned (6)  (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      character*32 connectname,donorname,boconame
      dimension irange(6),idrange(6),itransform(3),ipts(3,2)
      dimension numseg(6),normalindex(3)
c
c   First, need to get all connectivity info so can specify those segments
      call cg_n1to1_f(iccg,ibase,izone,n1to1,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   numseg(1-6) represents ilo,jlo,klo,ihi,jhi,khi, respectively
      do i=1,6
        numseg(i)=0
      enddo
      do m=1,n1to1
        call cg_1to1_read_f(iccg,ibase,izone,m,connectname,donorname,
     +   irange,idrange,itransform,ier)
        if (ier .ne. 0) call cg_error_exit_f
c       Increment if correct surface only
        if (    irange(1) .eq. 1 .and. irange(4) .eq. 1) then
          numseg(1)=numseg(1)+1
        else if(irange(1) .eq. irange(4)) then
          numseg(4)=numseg(4)+1
        else if(irange(2) .eq. 1 .and. irange(5) .eq. 1) then
          numseg(2)=numseg(2)+1
        else if(irange(2) .eq. irange(5)) then
          numseg(5)=numseg(5)+1
        else if(irange(3) .eq. 1 .and. irange(6) .eq. 1) then
          numseg(3)=numseg(3)+1
        else if(irange(3) .eq. irange(6)) then
          numseg(6)=numseg(6)+1
        end if
      enddo
c   Next get all BC info
      call cg_nbocos_f(iccg,ibase,izone,nbocos,ier)
      if (ier .ne. 0) call cg_error_exit_f
      do m=1,nbocos
        call cg_boco_info_f(iccg,ibase,izone,m,boconame,ibocotype,
     +    iptset,npts,normalindex,normallistflag,normaldatatype,
     +    ndataset,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (iptset .ne. PointRange) then
          write(901,'('' Error.  Currently BCs must be set up as'',
     +     '' PointRange type'',a32)') PointSetTypeName(iptset)
          stop
        end if
        call cg_boco_read_f(iccg,ibase,izone,m,ipts,normallist,ier)
        if (ier .ne. 0) call cg_error_exit_f
c       Increment if ilo surface only:
        if (    ipts(1,1) .eq. 1 .and. ipts(1,2) .eq. 1) then
          numseg(1)=numseg(1)+1
        else if(ipts(1,1) .eq. ipts(1,2)) then
          numseg(4)=numseg(4)+1
        else if(ipts(2,1) .eq. 1 .and. ipts(2,2) .eq. 1) then
          numseg(2)=numseg(2)+1
        else if(ipts(2,1) .eq. ipts(2,2)) then
          numseg(5)=numseg(5)+1
        else if(ipts(3,1) .eq. 1 .and. ipts(3,2) .eq. 1) then
          numseg(3)=numseg(3)+1
        else if(ipts(3,1) .eq. ipts(3,2)) then
          numseg(6)=numseg(6)+1
        end if
      enddo
c   
      return
      end
c
      subroutine getbc(iccg,ibase,izone,maxseg,ifacedesired,ibctyp,
     +                 ibcinf,iforc,numdata,bcval)
c***********************************************************************
c     Purpose: Determine the BCs for each segment for each of 6 faces
c     of a particular zone.
c     For now, this routine only counts BC segments and 1-to-1 segments.
c     Note: in current release of CFL3D, this routine is not used
c     (BC info is still obtained from standard input file)
c
c   INPUT:
c     iccg.............CGNS file index number (determined outside this
c                      routine) (integer)
c     ibase............CGNS base index number (determined outside this
c                      routine) (integer)
c     izone............zone (or grid) number (integer)
c     maxseg.........max number of segments per face allowed (integer)
c     ifacedesired...face desired (1-6 = ilo,jlo,klo,ihi,jhi,khi, 
c                    respectively) (integer)
c   OUTPUT:
c     ibctyp.........CFL3D-BC-type; dimensioned (maxseg) (integer)
c     ibcinf.........BC info (1-4 = ilo1,ihi1,ilo2,ihi2, where ilo1
c                    and ihi1 are the first non-face index (i or j),
c                    and ilo2 and ihi2 are the second non-face index
c                    (j or k); dimensioned (maxseg,4) (integer)
c     iforc..........indicates whether a force computation is to be
c                    ignored (1=computes force if appropriate, 
c                    0=override and do NOT compute);
c                    dimensioned (maxseg) (integer)
c     numdata........number of additional pieces of data needed;
c                    dimensioned (maxseg) (integer)
c     bcval..........additional BC info (if numdata .ne. 0)
c                    dimensioned (maxseg,5) (real)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      character*32 connectname,donorname,boconame
      character*7 nameseg
      dimension irange(6),idrange(6),itransform(3),ipts(3,2)
      dimension ibcinf(maxseg,4),bcval(maxseg,5),iforc(maxseg)
      dimension ibctyp(maxseg),numdata(maxseg),normalindex(3)
c
c   Determine appropriate indices for faces:
      if (    ifacedesired .eq. 1 .or. ifacedesired .eq. 4) then
        i1=1
        i2=4
        n1=2
        n2=5
        n3=3
        n4=6
      else if(ifacedesired .eq. 2 .or. ifacedesired .eq. 5) then
        i1=2
        i2=5
        n1=1
        n2=4
        n3=3
        n4=6
      else if(ifacedesired .eq. 3 .or. ifacedesired .eq. 6) then
        i1=3
        i2=6
        n1=1
        n2=4
        n3=2
        n4=5
      else
        write(901,'('' given value of ifacedesired not allowed='',i5)')
     +   ifacedesired
        stop
      end if
c   First, need to get all connectivity info so can specify those segments
      call cg_n1to1_f(iccg,ibase,izone,n1to1,ier)
      if (ier .ne. 0) call cg_error_exit_f
      numseg=0
      do m=1,n1to1
        call cg_1to1_read_f(iccg,ibase,izone,m,connectname,donorname,
     +   irange,idrange,itransform,ier)
        if (ier .ne. 0) call cg_error_exit_f
c       Increment if correct surface only
        if (ifacedesired .le. 3) then
          if (irange(i1) .eq. 1 .and. irange(i2) .eq. 1) then
            numseg=numseg+1
            ibctyp(numseg)=0
            ibcinf(numseg,1)=irange(n1)
            ibcinf(numseg,2)=irange(n2)
            ibcinf(numseg,3)=irange(n3)
            ibcinf(numseg,4)=irange(n4)
            iforc(numseg)=1
            numdata(numseg)=0
          end if
        else
          if ((irange(i1) .eq. irange(i2)) .and. irange(i1) .ne. 1) then
            numseg=numseg+1
            ibctyp(numseg)=0
            ibcinf(numseg,1)=irange(n1)
            ibcinf(numseg,2)=irange(n2)
            ibcinf(numseg,3)=irange(n3)
            ibcinf(numseg,4)=irange(n4)
            iforc(numseg)=1
            numdata(numseg)=0
          end if
        end if
      enddo
c   Next get all BC info
      call cg_nbocos_f(iccg,ibase,izone,nbocos,ier)
      if (ier .ne. 0) call cg_error_exit_f
      do m=1,nbocos
        call cg_boco_info_f(iccg,ibase,izone,m,boconame,ibocotype,
     +    iptset,npts,normalindex,normallistflag,normaldatatype,
     +    ndataset,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (iptset .ne. PointRange) then
          write(901,'('' Error.  Currently BCs must be set up as'',
     +     '' PointRange type'',a32)') PointSetTypeName(iptset)
          stop
        end if
        call cg_boco_read_f(iccg,ibase,izone,m,ipts,normallist,ier)
        if (ier .ne. 0) call cg_error_exit_f
c       Increment if correct surface only
        if (ifacedesired .le. 3) then
          if (ipts(i1,1) .eq. 1 .and. ipts(i1,2) .eq. 1) then
            numseg=numseg+1
            numdata(numseg)=0
            if (     ibocotype .eq. BCInflowSupersonic) then
              ibctyp(numseg)=1000
            else if (ibocotype .eq. BCSymmetryPlane) then
              ibctyp(numseg)=1001
            else if (ibocotype .eq. BCExtrapolate) then
              ibctyp(numseg)=1002
            else if (ibocotype .eq. BCFarfield) then
              ibctyp(numseg)=1003
            else if (ibocotype .eq. BCWallInviscid) then
              ibctyp(numseg)=1005
            else if (ibocotype .eq. BCTunnelInflow) then
              ibctyp(numseg)=1008
            else if (ibocotype .eq. BCSymmetryPolar) then
              ibctyp(numseg)=1011
            else if (ibocotype .eq. BCDegenerateLine) then
              ibctyp(numseg)=1013
            else if (ibocotype .eq. BCWallViscousHeatFlux) then
              ibctyp(numseg)=2004
              numdata(numseg)=2
              bcval(numseg,1)=0.
              bcval(numseg,2)=0.
            else
              write(901,'('' BCtype '',i5,'' not implemented yet'')')
     +         ibocotype
            end if
            ibcinf(numseg,1)=ipts(n1,1)
            ibcinf(numseg,2)=ipts(n1,2)
            ibcinf(numseg,3)=ipts(n3,1)
            ibcinf(numseg,4)=ipts(n3,2)
            nameseg=boconame
            if(nameseg .eq. 'NoForce') then
              iforc(numseg)=0
            else
              iforc(numseg)=1
            end if
          end if
        else
          if ((ipts(i1,1) .eq. ipts(i1,2)) .and. ipts(i1,1) .ne. 1) then
            numseg=numseg+1
            numdata(numseg)=0
            if (     ibocotype .eq. BCInflowSupersonic) then
              ibctyp(numseg)=1000
            else if (ibocotype .eq. BCSymmetryPlane) then
              ibctyp(numseg)=1001
            else if (ibocotype .eq. BCExtrapolate) then
              ibctyp(numseg)=1002
            else if (ibocotype .eq. BCFarfield) then
              ibctyp(numseg)=1003
            else if (ibocotype .eq. BCWallInviscid) then
              ibctyp(numseg)=1005
            else if (ibocotype .eq. BCTunnelInflow) then
              ibctyp(numseg)=1008
            else if (ibocotype .eq. BCSymmetryPolar) then
              ibctyp(numseg)=1011
            else if (ibocotype .eq. BCDegenerateLine) then
              ibctyp(numseg)=1013
            else if (ibocotype .eq. BCWallViscousHeatFlux) then
              ibctyp(numseg)=2004
              numdata(numseg)=2
              bcval(numseg,1)=0.
              bcval(numseg,2)=0.
            else
              write(901,'('' BCtype '',i5,'' not implemented yet'')')
     +         ibocotype
            end if
            ibcinf(numseg,1)=ipts(n1,1)
            ibcinf(numseg,2)=ipts(n1,2)
            ibcinf(numseg,3)=ipts(n3,1)
            ibcinf(numseg,4)=ipts(n3,2)
            nameseg=boconame
            if(nameseg .eq. 'NoForce') then
              iforc(numseg)=0
            else
              iforc(numseg)=1
            end if
          end if
        end if
      enddo
c   
      return
      end
c
      subroutine writehist(iccg,ibase,ntt,rms,clw,cdw,cdpw,cdvw,cxw,
     +   cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,cftpw,cftvw,
     +   cfttotw,rmstr1,rmstr2,rmstr3,rmstr4,rmstr5,rmstr6,rmstr7,
     +   nneg1,nneg2,nneg3,nneg4,nneg5,nneg6,nneg7,nummem)
c***********************************************************************
c     Purpose: Writes history info to CGNS file.  
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      ntt..............no of iterations (integer)
c      rms..............rms of density residual (real)
c      clw,cdw,...,rmstr......other history variables (real)
c      nneg.............neg pts in turb solution (integer)
c      nummem...........one of dimensions for rmstr, nneg (integer)
c   OUTPUTS:
c      none
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension rms(ntt),clw(ntt),cdw(ntt),cdpw(ntt),cdvw(ntt),
     +  cxw(ntt),cyw(ntt),czw(ntt),cmxw(ntt),cmyw(ntt),cmzw(ntt),
     +  fmdotw(ntt),cftmomw(ntt),cftpw(ntt),cftvw(ntt),
     +  cfttotw(ntt)
      dimension rmstr1(ntt),rmstr2(ntt),rmstr3(ntt),rmstr4(ntt),
     +  rmstr5(ntt),rmstr6(ntt),rmstr7(ntt)
      dimension nneg1(ntt),nneg2(ntt),nneg3(ntt),nneg4(ntt),
     +  nneg5(ntt),nneg6(ntt),nneg7(ntt)
      character*32 name
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Goto base index
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Create History node
      call cg_convergence_write_f(ntt,'',ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Goto new History node
      call cg_goto_f(iccg,ibase,ier,'ConvergenceHistory_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create Residual node
      write(901,'('' Creating RSDMassRMS history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('RSDMassRMS',RealDouble,1,ntt,rms,ier)
      else
        call cg_array_write_f('RSDMassRMS',RealSingle,1,ntt,rms,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create Lift Coefficient node
      write(901,'('' Creating CoefLift history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefLift',RealDouble,1,ntt,clw,ier)
      else
        call cg_array_write_f('CoefLift',RealSingle,1,ntt,clw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create Drag Coefficient node
      write(901,'('' Creating CoefDrag history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefDrag',RealDouble,1,ntt,cdw,ier)
      else
        call cg_array_write_f('CoefDrag',RealSingle,1,ntt,cdw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create Pressure Drag Coefficient node
      write(901,'('' Creating CoefPressureDrag history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefPressureDrag',RealDouble,1,ntt,
     +       cdpw,ier)
      else
        call cg_array_write_f('CoefPressureDrag',RealSingle,1,ntt,
     +       cdpw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create Viscous Drag Coefficient node
      write(901,'('' Creating CoefViscousDrag history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefViscousDrag',RealDouble,1,ntt,
     +       cdvw,ier)
      else
        call cg_array_write_f('CoefViscousDrag',RealSingle,1,ntt,
     +       cdvw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create x-direction force coefficient node
      write(901,'('' Creating CoefForceX history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefForceX',RealDouble,1,ntt,
     +       cxw,ier)
      else
        call cg_array_write_f('CoefForceX',RealSingle,1,ntt,
     +       cxw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create y-direction force coefficient node
      write(901,'('' Creating CoefForceY history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefForceY',RealDouble,1,ntt,
     +       cyw,ier)
      else
        call cg_array_write_f('CoefForceY',RealSingle,1,ntt,
     +       cyw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create z-direction force coefficient node
      write(901,'('' Creating CoefForceZ history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefForceZ',RealDouble,1,ntt,
     +       czw,ier)
      else
        call cg_array_write_f('CoefForceZ',RealSingle,1,ntt,
     +       czw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create moment coefficient about x node
      write(901,'('' Creating CoefMomentX history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefMomentX',RealDouble,1,ntt,
     +       cmxw,ier)
      else
        call cg_array_write_f('CoefMomentX',RealSingle,1,ntt,
     +       cmxw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create moment coefficient about y node
      write(901,'('' Creating CoefMomentY history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefMomentY',RealDouble,1,ntt,
     +       cmyw,ier)
      else
        call cg_array_write_f('CoefMomentY',RealSingle,1,ntt,
     +       cmyw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create moment coefficient about z node
      write(901,'('' Creating CoefMomentZ history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('CoefMomentZ',RealDouble,1,ntt,
     +       cmzw,ier)
      else
        call cg_array_write_f('CoefMomentZ',RealSingle,1,ntt,
     +       cmzw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create mass-flow-thru-control-surface node
      write(901,'('' Creating ControlSurfaceMassFlow history node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('ControlSurfaceMassFlow',
     +       RealDouble,1,ntt,fmdotw,ier)
      else
        call cg_array_write_f('ControlSurfaceMassFlow',
     +       RealSingle,1,ntt,fmdotw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create momentum-force-on-control-surface node
      write(901,'('' Creating ControlSurfaceMomentumForce history'',
     +  '' node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('ControlSurfaceMomentumForce',
     +       RealDouble,1,ntt,cftmomw,ier)
      else
        call cg_array_write_f('ControlSurfaceMomentumForce',
     +       RealSingle,1,ntt,cftmomw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create pressure-force-on-control-surface node
      write(901,'('' Creating ControlSurfacePressureForce history'',
     +  '' node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('ControlSurfacePressureForce',
     +       RealDouble,1,ntt,cftpw,ier)
      else
        call cg_array_write_f('ControlSurfacePressureForce',
     +       RealSingle,1,ntt,cftpw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create viscous-force-on-control-surface node
      write(901,'('' Creating ControlSurfaceViscousForce history'',
     +  '' node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('ControlSurfaceViscousForce',
     +       RealDouble,1,ntt,cftvw,ier)
      else
        call cg_array_write_f('ControlSurfaceViscousForce',
     +       RealSingle,1,ntt,cftvw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create total-force-on-control-surface node
      write(901,'('' Creating ControlSurfaceForce history'',
     +  '' node.'')')
      if (idouble .eq. 1) then
        call cg_array_write_f('ControlSurfaceForce',
     +       RealDouble,1,ntt,cfttotw,ier)
      else
        call cg_array_write_f('ControlSurfaceForce',
     +       RealSingle,1,ntt,cfttotw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create turb-equation-residual-history nodes
      do l=1,nummem
        if (l .eq. 1) name='RSDTurbEqn1RMS'
        if (l .eq. 2) name='RSDTurbEqn2RMS'
        if (l .eq. 3) name='RSDTurbEqn3RMS'
        if (l .eq. 4) name='RSDTurbEqn4RMS'
        if (l .eq. 5) name='RSDTurbEqn5RMS'
        if (l .eq. 6) name='RSDTurbEqn6RMS'
        if (l .eq. 7) name='RSDTurbEqn7RMS'
      write(901,'('' Creating RSDTurbEqnxRMS history'',
     +  '' node'')')
      write(901,'('' (where x='',i5,'')'')') l
      if (idouble .eq. 1) then
        if (l .eq. 1) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr1,ier)
        if (l .eq. 2) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr2,ier)
        if (l .eq. 3) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr3,ier)
        if (l .eq. 4) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr4,ier)
        if (l .eq. 5) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr5,ier)
        if (l .eq. 6) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr6,ier)
        if (l .eq. 7) call cg_array_write_f(name,
     +       RealDouble,1,ntt,rmstr7,ier)
      else
        if (l .eq. 1) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr1,ier)
        if (l .eq. 2) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr2,ier)
        if (l .eq. 3) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr3,ier)
        if (l .eq. 4) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr4,ier)
        if (l .eq. 5) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr5,ier)
        if (l .eq. 6) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr6,ier)
        if (l .eq. 7) call cg_array_write_f(name,
     +       RealSingle,1,ntt,rmstr7,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Create turb-equation-1-bad-points-history node
      do l=1,nummem
        if (l .eq. 1) name='BadPointsTurbEqn1'
        if (l .eq. 2) name='BadPointsTurbEqn2'
        if (l .eq. 3) name='BadPointsTurbEqn3'
        if (l .eq. 4) name='BadPointsTurbEqn4'
        if (l .eq. 5) name='BadPointsTurbEqn5'
        if (l .eq. 6) name='BadPointsTurbEqn6'
        if (l .eq. 7) name='BadPointsTurbEqn7'
      write(901,'('' Creating BadPointsTurbEqnx history'',
     +  '' node'')')
      write(901,'('' (where x='',i5,'')'')') l
      if (l .eq. 1) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg1,ier)
      if (l .eq. 2) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg2,ier)
      if (l .eq. 3) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg3,ier)
      if (l .eq. 4) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg4,ier)
      if (l .eq. 5) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg5,ier)
      if (l .eq. 6) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg6,ier)
      if (l .eq. 7) call cg_array_write_f(name,
     +       Integer,1,ntt,nneg7,ier)
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
      return
      end
c
      subroutine readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     +   cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     +   cftpw,cftvw,cfttotw,rmstr1,rmstr2,rmstr3,rmstr4,rmstr5,
     +   rmstr6,rmstr7,nneg1,nneg2,nneg3,nneg4,nneg5,nneg6,nneg7,
     +   nummem)
c***********************************************************************
c     Purpose: Reads history info from CGNS file.  
c      (NOTE:  cg_array_read_as_f may eventually be converted
c      to cg_array_read_f)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      ncycmax..........max number of iterations that can be saved
c                       (integer)
c      nummem...........one of dimensions for rmstr, nneg (integer)
c   OUTPUTS:
c      ntr..............no of iterations (integer)
c      rms..............rms of density residual (real)
c      clw,cdw,...,rmstr......other history variables (real)
c      nneg.............neg pts in turb solution (integer)
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=30)
c
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),
     +  cdpw(ncycmax),cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),
     +  czw(ncycmax),cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     +  fmdotw(ncycmax),cftmomw(ncycmax),cftpw(ncycmax),
     +  cftvw(ncycmax),cfttotw(ncycmax)
      dimension rmstr1(ncycmax),rmstr2(ncycmax),rmstr3(ncycmax),
     +  rmstr4(ncycmax),rmstr5(ncycmax),rmstr6(ncycmax),
     +  rmstr7(ncycmax)
      dimension nneg1(ncycmax),nneg2(ncycmax),nneg3(ncycmax),
     +  nneg4(ncycmax),nneg5(ncycmax),nneg6(ncycmax),
     +  nneg7(ncycmax)
      dimension idimvec(4)
      character*32 name(numnames),name2,searchname
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Goto base index
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Read number of iterations of data in ConvergenceHistory node
      call cg_convergence_read_f(ntr,name2,ier)
      if (ier .ne. 0) then
         write(901,'('' ConvergenceHistory node does not exist.'')')
         write(901,'(''   setting ntr=0 and continuing'')')
         ntr=0
         return
      end if
      write(901,'('' reading convergence hist data... ntr='',i5)') ntr
      if (ntr .ge. ncycmax) then
      write(901,1239)
 1239 format(/,1x,11hstopping...,
     .       40hprevious number of iterations computed >,
     .       1x,18h dimension ncycmax)
      write(901,*)' ntr,ncycmax = ',ntr,ncycmax
      write(901,*)' increase value of ncycmax to at LEAST ',
     .ntr+ncycmax
      stop
      end if
c
c   Get Information about what is under History node
      call cg_goto_f(iccg,ibase,ier,'ConvergenceHistory_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' WARNING.  Too many history arrays.'')')
         write(901,'(''   setting ntr=0 and continuing'')')
         ntr=0
         return
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Get Residual
      do n=1,narrays
        if (name(n) .eq. 'RSDMassRMS') goto 101
      enddo
      write(901,'('' WARNING. No RSDMassRMS node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 101  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rms,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rms,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Lift Coefficient
      do n=1,narrays
        if (name(n) .eq. 'CoefLift') goto 102
      enddo
      write(901,'('' WARNING. No CoefLift node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 102  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,clw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,clw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Drag Coefficient
      do n=1,narrays
        if (name(n) .eq. 'CoefDrag') goto 103
      enddo
      write(901,'('' WARNING. No CoefDrag node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 103  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cdw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cdw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Pressure Drag Coefficient
      do n=1,narrays
        if (name(n) .eq. 'CoefPressureDrag') goto 104
      enddo
      write(901,'('' WARNING. No CoefPressureDrag node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 104  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cdpw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cdpw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Viscous Drag Coefficient
      do n=1,narrays
        if (name(n) .eq. 'CoefViscousDrag') goto 105
      enddo
      write(901,'('' WARNING. No CoefViscousDrag node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cdvw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cdvw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Force coeff in x-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefForceX') goto 106
      enddo
      write(901,'('' WARNING. No CoefForceX node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cxw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cxw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Force coeff in y-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefForceY') goto 107
      enddo
      write(901,'('' WARNING. No CoefForceY node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 107  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cyw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cyw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Force coeff in z-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefForceZ') goto 108
      enddo
      write(901,'('' WARNING. No CoefForceZ node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,czw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,czw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Moment coeff in x-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefMomentX') goto 109
      enddo
      write(901,'('' WARNING. No CoefMomentX node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 109  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cmxw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cmxw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Moment coeff in y-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefMomentY') goto 110
      enddo
      write(901,'('' WARNING. No CoefMomentY node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 110  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cmyw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cmyw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Moment coeff in z-direction
      do n=1,narrays
        if (name(n) .eq. 'CoefMomentZ') goto 111
      enddo
      write(901,'('' WARNING. No CoefMomentZ node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 111  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cmzw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cmzw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Control surface mass flow
      do n=1,narrays
        if (name(n) .eq. 'ControlSurfaceMassFlow') goto 112
      enddo
      write(901,'('' WARNING. No ControlSurfaceMassFlow node exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 112  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,fmdotw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,fmdotw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Control surface momentum force
      do n=1,narrays
        if (name(n) .eq. 'ControlSurfaceMomentumForce') goto 113
      enddo
      write(901,'('' WARNING. No ControlSurfaceMomentumForce node'',
     +  '' exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 113  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cftmomw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cftmomw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Control surface pressure force
      do n=1,narrays
        if (name(n) .eq. 'ControlSurfacePressureForce') goto 114
      enddo
      write(901,'('' WARNING. No ControlSurfacePressureForce node'',
     +  '' exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 114  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cftpw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cftpw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Control surface viscous force
      do n=1,narrays
        if (name(n) .eq. 'ControlSurfaceViscousForce') goto 115
      enddo
      write(901,'('' WARNING. No ControlSurfaceViscousForce node'',
     +  '' exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 115  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cftvw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cftvw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Control surface total force
      do n=1,narrays
        if (name(n) .eq. 'ControlSurfaceForce') goto 116
      enddo
      write(901,'('' WARNING. No ControlSurfaceForce node'',
     +  '' exists'')')
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 116  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cfttotw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cfttotw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Residual of turb eqns
      do l=1,nummem
        if (l .eq. 1) searchname='RSDTurbEqn1RMS'
        if (l .eq. 2) searchname='RSDTurbEqn2RMS'
        if (l .eq. 3) searchname='RSDTurbEqn3RMS'
        if (l .eq. 4) searchname='RSDTurbEqn4RMS'
        if (l .eq. 5) searchname='RSDTurbEqn5RMS'
        if (l .eq. 6) searchname='RSDTurbEqn6RMS'
        if (l .eq. 7) searchname='RSDTurbEqn7RMS'
      do n=1,narrays
        if (name(n) .eq. searchname) goto 117
      enddo
      write(901,'('' WARNING. No RSDTurbEqnxRMS node exists'')')
      write(901,'(''   (where x='',i5,'')'')') l
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 117  continue
      if (idouble .eq. 1) then
        if (l.eq.1) call cg_array_read_as_f(n,RealDouble,rmstr1,ier)
        if (l.eq.2) call cg_array_read_as_f(n,RealDouble,rmstr2,ier)
        if (l.eq.3) call cg_array_read_as_f(n,RealDouble,rmstr3,ier)
        if (l.eq.4) call cg_array_read_as_f(n,RealDouble,rmstr4,ier)
        if (l.eq.5) call cg_array_read_as_f(n,RealDouble,rmstr5,ier)
        if (l.eq.6) call cg_array_read_as_f(n,RealDouble,rmstr6,ier)
        if (l.eq.7) call cg_array_read_as_f(n,RealDouble,rmstr7,ier)
      else
        if (l.eq.1) call cg_array_read_as_f(n,RealSingle,rmstr1,ier)
        if (l.eq.2) call cg_array_read_as_f(n,RealSingle,rmstr2,ier)
        if (l.eq.3) call cg_array_read_as_f(n,RealSingle,rmstr3,ier)
        if (l.eq.4) call cg_array_read_as_f(n,RealSingle,rmstr4,ier)
        if (l.eq.5) call cg_array_read_as_f(n,RealSingle,rmstr5,ier)
        if (l.eq.6) call cg_array_read_as_f(n,RealSingle,rmstr6,ier)
        if (l.eq.7) call cg_array_read_as_f(n,RealSingle,rmstr7,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get turb eqn bad points
      do l=1,nummem
        if (l .eq. 1) searchname='BadPointsTurbEqn1'
        if (l .eq. 2) searchname='BadPointsTurbEqn2'
        if (l .eq. 3) searchname='BadPointsTurbEqn3'
        if (l .eq. 4) searchname='BadPointsTurbEqn4'
        if (l .eq. 5) searchname='BadPointsTurbEqn5'
        if (l .eq. 6) searchname='BadPointsTurbEqn6'
        if (l .eq. 7) searchname='BadPointsTurbEqn7'
      do n=1,narrays
        if (name(n) .eq. searchname) goto 119
      enddo
      write(901,'('' WARNING. No BadPointsTurbEqnx node exists'')')
      write(901,'(''   (where x='',i5,'')'')') l
      write(901,'(''   setting ntr=0 and continuing'')')
      ntr=0
      return
 119  continue
      if (l.eq.1) call cg_array_read_as_f(n,Integer,nneg1,ier)
      if (l.eq.2) call cg_array_read_as_f(n,Integer,nneg2,ier)
      if (l.eq.3) call cg_array_read_as_f(n,Integer,nneg3,ier)
      if (l.eq.4) call cg_array_read_as_f(n,Integer,nneg4,ier)
      if (l.eq.5) call cg_array_read_as_f(n,Integer,nneg5,ier)
      if (l.eq.6) call cg_array_read_as_f(n,Integer,nneg6,ier)
      if (l.eq.7) call cg_array_read_as_f(n,Integer,nneg7,ier)
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
      return
      end
c
      subroutine writeeqn(iccg,ibase,igrid,iv1,iv2,iv3,ideg,gamma,
     +   prandtl,prandtlt,suth,i2d,ifullns,ides,isarcxd,ieasmcc2d,
     +   isstrc,isar)
c***********************************************************************
c     Purpose: Writes equation info to CGNS data base for zone number 
c     "igrid". The CGNS file must be already opened in "MODE_MODIFY"
c      mode.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      iv1,iv2,iv3......ivisc(n) parameters from CFL3D (integers)
c      ideg(3)..........degani-schiff turb model parameter from CFL3D 
c                       (integers)
c      gamma............specific heat ratio (real)
c      prandtl..........Prandtl number (real)
c      prandtlt.........turbulent Prandtl number (real)
c      suth.............Sutherland law constant (real)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      ifullns..........0 if thin layer, 1 if full (integer)
c      ides.............0 if DES off, 1 if DES, 2+ if DDES (integer)
c      isarcxd..........0 if SARC curv correction off, 1 if on (integer)
c      ieasmcc2d........0 if EASMCC curv correction off, 1 if on (integer)
c      isstrc...........0 if SSTRC off, 1 if on (integer)
c      isar.............0 if SAR curv correction off, 1 if on (integer)
c   OUTPUTS:
c      none
c***********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension ideg(3),idata(6),units(5)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
      ivmx=max(iv1,iv2)
      ivmx=max(ivmx,iv3)
c
c   Go to appropriate node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Write 'FlowEquationSet' node
      call cg_equationset_write_f(3,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create 'GoverningEquations' node under 'FlowEquationSet'
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if(ivmx .eq. 0) then
        call cg_governing_write_f(Euler,ier)
      else if (ivmx .eq. 1) then
        call cg_governing_write_f(NSLaminar,ier)
      else
        call cg_governing_write_f(NSTurbulent,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create 'DiffusionModel' node under 'GoverningEquations'
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'GoverningEquations_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      idata(1)=0
      idata(2)=0
      idata(3)=0
      idata(4)=0
      idata(5)=0
      idata(6)=0
      if(iv1 .gt. 0) idata(1)=1
      if(iv2 .gt. 0) idata(2)=1
      if(iv3 .gt. 0) idata(3)=1
      if(ifullns .eq. 1) then
        if(i2d .eq. 1) then
          idata(5)=1
        else
          idata(4)=1
          idata(5)=1
          idata(6)=1
        end if
      end if
      call cg_diffusion_write_f(idata,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create other nodes under 'FlowEquationSet'
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_model_write_f('GasModel_t',CaloricallyPerfect,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_model_write_f('ViscosityModel_t',
     + SutherlandLaw,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_model_write_f('ThermalConductivityModel_t',
     + ConstantPrandtl,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create DataArray containing Specific Heat Ratio under GasModel
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'GasModel_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('SpecificHeatRatio',RealDouble,
     +   1,1,gamma,ier)
      else
        call cg_array_write_f('SpecificHeatRatio',RealSingle,
     +   1,1,gamma,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'GasModel_t',1,'DataArray_t',
     +  1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NondimensionalParameter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create DataArray containing Prandtl number under ThermalConductivityModel
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'ThermalConductivityModel_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Prandtl',RealDouble,
     +   1,1,prandtl,ier)
      else
        call cg_array_write_f('Prandtl',RealSingle,
     +   1,1,prandtl,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'ThermalConductivityModel_t',1,
     +  'DataArray_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NondimensionalParameter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create DataArray containing Sutherland law constant under ViscosityModel
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'ViscosityModel_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('SutherlandLawConstant',RealDouble,
     +   1,1,suth,ier)
      else
        call cg_array_write_f('SutherlandLawConstant',RealSingle,
     +   1,1,suth,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'ViscosityModel_t',1,
     +  'DataArray_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(Dimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_units_write_f(Null,Null,Null,Rankine,Null,ier)
      if (ier .ne. 0) call cg_error_exit_f
      units(1)=0.
      units(2)=0.
      units(3)=0.
      units(4)=1.
      units(5)=0.
      if (idouble .eq. 1) then
        call cg_exponents_write_f(RealDouble,units,ier)
      else
        call cg_exponents_write_f(RealSingle,units,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c  Only write turbulence-type nodes if turbulent
      if (ivmx .le. 1) return
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (ivmx .eq. 11 .or. ivmx .eq. 12 .or.
     +    ivmx .eq. 13 .or. ivmx .eq. 14) then
        call cg_model_write_f('TurbulenceClosure_t',
     +    ReynoldsStressAlgebraic,ier)
      else
        call cg_model_write_f('TurbulenceClosure_t',
     +    EddyViscosity,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      idegg = max(ideg(1),ideg(2))
      idegg = max(ideg(3),idegg)
      if(ivmx .eq. 2 .and. idegg .eq. 0) then
c   (standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   Algebraic_BaldwinLomax,ier)
      else if(ivmx .eq. 2 .and. idegg .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'Algebraic_BaldwinLomax_DS',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 4) then
c   (standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   OneEquation_BaldwinBarth,ier)
      else if(ivmx .eq. 5 .and. ides .ne. 1 .and. isarcxd .ne. 1
     +  .and. isar .ne. 1) then
c   (standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   OneEquation_SpalartAllmaras,ier)
      else if(ivmx .eq. 5 .and. ides .eq. 1 .and. isarcxd .ne. 1
     +  .and. isar .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DetachedEddySimulation',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .ge. 2 .and. isarcxd .ne. 1
     +  .and. isar .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DDES',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .ne. 1 .and. isarcxd .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'SpalartAllmarasSARC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .eq. 1 .and. isarcxd .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DetachedEddySimulationSARC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .ge. 2 .and. isarcxd .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DDES-SARC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .ne. 1 .and. isar .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'SpalartAllmarasSAR',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .eq. 1 .and. isar .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DetachedEddySimulationSAR',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 5 .and. ides .ge. 2 .and. isar .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'DDES-SAR',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 6 .and. ides .ne. 1 .and. isstrc .ne. 1) then
c   (standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   TwoEquation_Wilcox,ier)
      else if(ivmx .eq. 6 .and. ides .ne. 1 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_Wilcox-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 6 .and. ides .eq. 1 .and. isstrc .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_WilcoxDES',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 6 .and. ides .eq. 1 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_WilcoxDES-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 6 .and. ides .ge. 2 .and. isstrc .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_WilcoxDDES',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 6 .and. ides .ge. 2 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_WilcoxDDES-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 7 .and. ides .ne. 1 .and. isstrc .ne. 1) then
c   (standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   TwoEquation_MenterSST,ier)
      else if(ivmx .eq. 7 .and. ides .ne. 1 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_MenterSST-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 7 .and. ides .eq. 1 .and. isstrc .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_MenterSSTDES',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 7 .and. ides .eq. 1 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_MenterSSTDES-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 7 .and. ides .ge. 2 .and. isstrc .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_MenterSSTDDES',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 7 .and. ides .ge. 2 .and. isstrc .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_MenterSSTDDES-RC',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 8 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASM_KO_Lin',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 8 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASMCC_KO_Lin',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 9 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASM_KE_Lin',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 9 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASMCC_KE_Lin',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 10) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_Abid_KE',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 11 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_GS_EASM_KE',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 11 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_GS_EASMCC_KE',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 12 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_GS_EASM_KO',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 12 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_GS_EASMCC_KO',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 13 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASM_KE',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 13 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASMCC_KE',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 14 .and. ieasmcc2d .ne. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASM_KO',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 14 .and. ieasmcc2d .eq. 1) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_EASMCC_KO',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 15) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'TwoEquation_K_Enstrophy',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 25) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'LES_Smagorinsky_Cs_const',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 30) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'ThreeEquation_SST+transition',ier)
        if (ier .ne. 0) call cg_error_exit_f
      else if(ivmx .eq. 40) then
c   (non-standard)
        call cg_model_write_f('TurbulenceModel_t',
     +   UserDefined,ier)
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_descriptor_write_f('TurbulenceModelName',
     +   'FourEquation_SST+transition',ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
c
c   Create DataArray containing turb Prandtl number under TurbulenceClosure
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'TurbulenceClosure_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('PrandtlTurbulent',RealDouble,
     +   1,1,prandtlt,ier)
      else
        call cg_array_write_f('PrandtlTurbulent',RealSingle,
     +   1,1,prandtlt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'TurbulenceClosure_t',1,
     +  'DataArray_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NondimensionalParameter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
c
      subroutine readeqn(iccg,ibase,igrid,iv1,iv2,iv3,gamma,
     +   prandtl,prandtlt,suth)
c**********************************************************************
c     Purpose: Reads equation info from CGNS data base for zone number
c     "igrid".
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c   OUTPUTS:
c      iv1,iv2,iv3......ivisc(n) parameters (integers)
c      gamma............specific heat ratio (real)
c      prandtl..........Prandtl number (real)
c      prandtlt.........turbulent Prandtl number (real)
c      suth.............Sutherland law constant (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension idata(6),idimvec(4)
      character name*32,text*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Go to appropriate node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Read 'FlowEquationSet' node
      call cg_equationset_read_f(id,ige,igm,ivm,itcm,itc,itm,ier)
      if (ier .gt. 0) then
        write(901,'('' FlowEquationSet node does not exist.'')')
        write(901,'(''   assuming gamma=1.4, pr=0.72, prt=0.9,'',
     +   '' suth=198.6, and restart from Euler'')')
        iv1=0
        iv2=0
        iv3=0
        gamma=1.4
        prandtl=0.72
        prandtlt=0.9
        suth=198.6
        return
      end if
c   Read 'GoverningEquations' node
      if (ige .eq. 1) then
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_governing_read_f(itype,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (itype .eq. Euler .or.
     +      itype .eq. NSLaminar .or.
     +      itype .eq. NSTurbulent) then
          write(901,'('' cgns read solution of:  '',a32)') 
     +     GoverningEquationsTypeName(itype)
        else
          write(901,'('' WARNING! GoverningEquations type not'',
     +     '' recognized:  '',a32)') GoverningEquationsTypeName(itype)
          write(901,'(''   assuming gamma=1.4, pr=0.72, prt=0.9,'',
     +     '' suth=198.6, and restart from Euler'')')
          iv1=0
          iv2=0
          iv3=0
          gamma=1.4
          prandtl=0.72
          prandtlt=0.9
          suth=198.6
          return
        end if
      else
        write(901,'('' GoverningEquations node does not exist. '')')
        write(901,'(''   assuming gamma=1.4, pr=0.72, prt=0.9,'',
     +   '' suth=198.6, and restart from Euler'')')
        iv1=0
        iv2=0
        iv3=0
        gamma=1.4
        prandtl=0.72
        prandtlt=0.9
        suth=198.6
        return
      end if
c   Get gamma
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'GasModel_t',1,'end')
      if (ier .ne. 0) then
        write(901,'('' No GasModel exists'')')
        write(901,'(''   Setting gamma=1.4 & continuing'')')
        gamma=1.4
        goto 205
      end if
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 1) then
        write(901,'('' Expecting only one array under GasModel'')')
        write(901,'(''   Read '',i5)') narrays
        write(901,'(''   Setting gamma=1.4 & continuing'')')
        gamma=1.4
        goto 205
      end if
      call cg_array_info_f(1,name,ltype,idatadim,idimvec,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (name .eq. 'SpecificHeatRatio') goto 105
      write(901,'('' WARNING. No SpecificHeatRatio node exists'')')
      write(901,'(''   setting to 1.4 & continuing'')')
      gamma=1.4
      goto 205
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(1,RealDouble,gamma,ier)
      else
        call cg_array_read_as_f(1,RealSingle,gamma,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 205  continue
c   Get prandtl
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'ThermalConductivityModel_t',1,'end')
      if (ier .ne. 0) then
        write(901,'('' No ThermalConductivityModel exists'')')
        write(901,'(''   Setting prandtl=0.72 & continuing'')')
        prandtl=0.72
        goto 206
      end if
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 1) then
        write(901,'('' Expecting only one array under '',
     +    ''ThermalConductivityModel'')')
        write(901,'(''   Read '',i5)') narrays
        write(901,'(''   Setting prandtl=0.72 & continuing'')')
        prandtl=0.72
        goto 206
      end if
      call cg_array_info_f(1,name,ltype,idatadim,idimvec,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (name .eq. 'Prandtl') goto 106
      write(901,'('' WARNING. No Prandtl node exists'')')
      write(901,'(''   setting to 0.72 & continuing'')')
      prandtl=0.72
      goto 206
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(1,RealDouble,prandtl,ier)
      else
        call cg_array_read_as_f(1,RealSingle,prandtl,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 206  continue
c   Get prandtlt
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceClosure_t',1,'end')
      if (ier .ne. 0) then
        write(901,'('' No TurbulenceClosure exists'')')
        write(901,'(''   Setting prandtlt=0.9 & continuing'')')
        prandtlt=0.9
        goto 207
      end if
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 1) then
        write(901,'('' Expecting only one array under '',
     +    ''TurbulenceClosure'')')
        write(901,'(''   Read '',i5)') narrays
        write(901,'(''   Setting prandtlt=0.9 & continuing'')')
        prandtlt=0.9
        goto 207
      end if
      call cg_array_info_f(1,name,ltype,idatadim,idimvec,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (name .eq. 'PrandtlTurbulent') goto 107
      write(901,'('' WARNING. No PrandtlTurbulent node exists'')')
      write(901,'(''   setting to 0.9 & continuing'')')
      prandtlt=0.9
      goto 207
 107  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(1,RealDouble,prandtlt,ier)
      else
        call cg_array_read_as_f(1,RealSingle,prandtlt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 207  continue
c   Get suth
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'ViscosityModel_t',1,'end')
      if (ier .ne. 0) then
        write(901,'('' No ViscosityModel exists'')')
        write(901,'(''   Setting suth=198.6 & continuing'')')
        suth=198.6
        goto 208
      end if
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 1) then
        write(901,'('' Expecting only one array under '',
     +    ''ViscosityModel'')')
        write(901,'(''   Read '',i5)') narrays
        write(901,'(''   Setting suth=198.6 & continuing'')')
        suth=198.6
        goto 208
      end if
      call cg_array_info_f(1,name,ltype,idatadim,idimvec,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (name .eq. 'SutherlandLawConstant') goto 108
      write(901,'('' WARNING. No SutherlandLawConstant node exists'')')
      write(901,'(''   setting to 198.6 & continuing'')')
      suth=198.6
      goto 208
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(1,RealDouble,suth,ier)
      else
        call cg_array_read_as_f(1,RealSingle,suth,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 208  continue
c
c   Read 'DiffusionModel' node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'FlowEquationSet_t',1,'GoverningEquations_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_diffusion_read_f(idata,ier)
      if (ier .gt. 0) then
        write(901,'('' DiffusionModel node does not exist.'')')
        write(901,'(''   assuming restart from Euler'')')
        iv1=0
        iv2=0
        iv3=0
        return
      end if
      if (idata(4) .ne. 0 .or. idata(5) .ne. 0 .or.
     +    idata(6) .ne. 0) then
        write(901,'('' Restarting from FULL N-S'',
     +  '' (X-derivatives included)'')')
      end if
      if (itype .eq. Euler .and. 
     + (idata(1).eq.1 .or. idata(2).eq.1 .or. idata(3).eq.1)) then
        write(901,'('' Inconsistency in ADF file between'',
     +  '' GoverningEquations and DiffusionModel'')')
        write(901,'(''    using DiffusionModel...'')')
      end if
c   Read 'TurbulenceModel' node, if appropriate
      if (itype .eq. Euler .or.
     +    itype .eq. NSLaminar) then
        iv1=idata(1)
        iv2=idata(2)
        iv3=idata(3)
      else
      if (itm .eq. 1) then
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_model_read_f('TurbulenceModel_t',itype,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (itype .eq. UserDefined) then
          call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'FlowEquationSet_t',1,'TurbulenceModel_t',1,'end')
          if (ier .ne. 0) call cg_error_exit_f
          call cg_descriptor_read_f(1,name,text,ier)
          if (ier .ne. 0) call cg_error_exit_f
        end if
        if (itype .eq. Algebraic_BaldwinLomax) then
          iv1=idata(1)*2
          iv2=idata(2)*2
          iv3=idata(3)*2
        else if(itype .eq. UserDefined .and. 
     +   text .eq. 'Algebraic_BaldwinLomax_DS') then
          iv1=idata(1)*2
          iv2=idata(2)*2
          iv3=idata(3)*2
        else if(itype .eq. OneEquation_BaldwinBarth) then
          iv1=idata(1)*4
          iv2=idata(2)*4
          iv3=idata(3)*4
        else if(itype .eq. OneEquation_SpalartAllmaras) then
          iv1=idata(1)*5
          iv2=idata(2)*5
          iv3=idata(3)*5
        else if(itype .eq. TwoEquation_Wilcox) then
          iv1=idata(1)*6
          iv2=idata(2)*6
          iv3=idata(3)*6
        else if(itype .eq. TwoEquation_MenterSST) then
          iv1=idata(1)*7
          iv2=idata(2)*7
          iv3=idata(3)*7
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASM_KO_Lin') then
          iv1=idata(1)*8
          iv2=idata(2)*8
          iv3=idata(3)*8
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASM_KE_Lin') then
          iv1=idata(1)*9
          iv2=idata(2)*9
          iv3=idata(3)*9
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_Abid_KE') then
          iv1=idata(1)*10
          iv2=idata(2)*10
          iv3=idata(3)*10
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_GS_EASM_KE') then
          iv1=idata(1)*11
          iv2=idata(2)*11
          iv3=idata(3)*11
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_GS_EASM_KO') then
          iv1=idata(1)*12
          iv2=idata(2)*12
          iv3=idata(3)*12
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASM_KE') then
          iv1=idata(1)*13
          iv2=idata(2)*13
          iv3=idata(3)*13
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASM_KO') then
          iv1=idata(1)*14
          iv2=idata(2)*14
          iv3=idata(3)*14
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_K_Enstrophy') then
          iv1=idata(1)*15
          iv2=idata(2)*15
          iv3=idata(3)*15
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'DetachedEddySimulation') then
          iv1=idata(1)*5
          iv2=idata(2)*5
          iv3=idata(3)*5
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'SpalartAllmarasSARC') then
          iv1=idata(1)*5
          iv2=idata(2)*5
          iv3=idata(3)*5
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'DetachedEddySimulationSARC') then
          iv1=idata(1)*5
          iv2=idata(2)*5
          iv3=idata(3)*5
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASMCC_KO_Lin') then
          iv1=idata(1)*8
          iv2=idata(2)*8
          iv3=idata(3)*8
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASMCC_KE_Lin') then
          iv1=idata(1)*9
          iv2=idata(2)*9
          iv3=idata(3)*9
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_GS_EASMCC_KE') then
          iv1=idata(1)*11
          iv2=idata(2)*11
          iv3=idata(3)*11
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_GS_EASMCC_KO') then
          iv1=idata(1)*12
          iv2=idata(2)*12
          iv3=idata(3)*12
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASMCC_KE') then
          iv1=idata(1)*13
          iv2=idata(2)*13
          iv3=idata(3)*13
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'TwoEquation_EASMCC_KO') then
          iv1=idata(1)*14
          iv2=idata(2)*14
          iv3=idata(3)*14
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'LES_Smagorinsky_Cs_const') then
          iv1=idata(1)*25
          iv2=idata(2)*25
          iv3=idata(3)*25
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'ThreeEquation_SST+transition') then
          iv1=idata(1)*30
          iv2=idata(2)*30
          iv3=idata(3)*30
        else if(itype .eq. UserDefined .and.
     +   text .eq. 'FourEquation_SST+transition') then
          iv1=idata(1)*40
          iv2=idata(2)*40
          iv3=idata(3)*40
        else
          write(901,'('' Model name not recognized:'',a32)')
     +     ModelTypeName(itype)
          write(901,'(''   assuming NSLaminar'')')
          iv1=idata(1)
          iv2=idata(2)
          iv3=idata(3)
        end if
      else
        write(901,'('' TurbulenceModel node does not exist. '')')
        if (itype .eq. Euler) then
          write(901,'(''   (restart is Euler)'')')
        else
          write(901,'(''   (restart is NSLaminar)'')')
        end if
        iv1=idata(1)
        iv2=idata(2)
        iv3=idata(3)
      end if
      end if
c
      return
      end
c
      subroutine writetime(iccg,ibase,time,iter,dt)
c**********************************************************************
c     Purpose: Writes latest time info to CGNS file.  This includes
c     the time step dt, which will only be read in (by rsecord) if
c     the run is 2nd-order time-accurate.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      time.............last time (for time-accurate) (real)
c      iter.............last iteration (integer)
c      dt...............time step (real)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=5)
c
      character name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   For CFL3D, nsteps written to CGNS file is ALWAYS 1:
      nsteps=1
      call cg_biter_write_f(iccg,ibase,'TimeIterValues',nsteps,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Go to BaseIterativeData level and write time value
      call cg_goto_f(iccg,ibase,ier,'BaseIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('TimeValues',RealDouble,1,1,time,ier)
      else
        call cg_array_write_f('TimeValues',RealSingle,1,1,time,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'BaseIterativeData_t',1,
     +  'DataArray_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Go to BaseIterativeData level and write iteration value
      call cg_goto_f(iccg,ibase,ier,'BaseIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('IterationValues',Integer,1,1,iter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to Base node, write CFL3DTimeStep array as UserDefinedData
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_user_data_write_f('CFL3DTimeStep', ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nuser_data_f(nuserdata, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .gt. numnames) then
        write(901,'('' Error... numnames too small in writetime'')')
        stop
      end if
      do n=1,nuserdata
        call cg_user_data_read_f(n, name(n), ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(n) .eq. 'CFL3DTimeStep') then
          nset=n
        end if
      enddo
      call cg_goto_f(iccg,ibase,ier,'UserDefinedData_t',nset,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('LatestDTUsed',RealDouble,1,1,dt,ier)
      else
        call cg_array_write_f('LatestDTUsed',RealSingle,1,1,dt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'UserDefinedData_t',nset,
     +  'DataArray_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine readtime(iccg,ibase,time)
c**********************************************************************
c     Purpose: Reads latest time info from CGNS file (latest 
c     iteration is not read from here).
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c   OUTPUTS:
c      time.............last time (for time-accurate) (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=3)
c
      dimension idimvec(4)
      character name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Read BaseIterativeData node
      call cg_biter_read_f(iccg,ibase,name,nsteps,ier)
      if (ier .gt. 0) then
         write(901,'('' WARNING: BaseIterativeData node does not '',
     +    ''exist.'')')
         write(901,'(''   Assuming time=0 and continuing'')')
         time=0.
         return
      end if
c
c   Go to BaseIterativeData node and get info
      call cg_goto_f(iccg,ibase,ier,'BaseIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
         write(901,'('' WARNING.  Too many BaseIterativeData'',
     +    '' arrays.'')')
         write(901,'(''   Setting time=0 and continuing'')')
         time=0.
         return
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get TimeValues
      do n=1,narrays
        if (name(n) .eq. 'TimeValues') goto 101
      enddo
      write(901,'('' WARNING. No TimeValues node exists'')')
      write(901,'(''   setting time=0 & continuing'')')
      time=0.
      goto 201
 101  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,time,ier)
      else
        call cg_array_read_as_f(n,RealSingle,time,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 201  continue
      return
      end
c
      subroutine writeref(iccg,ibase,xmach,reue,
     +                    rho0,c0,p0,vk0,xlength0,tinf,alphw,betaw,
     +                    u0,v0,w0,ialph)
c**********************************************************************
c     Purpose: Writes reference (nondimensionalization) info to
c     CGNS file.  
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      xmach............Mach number (real)
c      reue.............Reynolds number (real)
c      rho0.............reference nondim density (real)
c      c0...............reference nondim speed of sound (real)
c      p0...............reference nondim pressure (real)
c      vk0..............reference nondim molecular viscosity (real) 
c                       (CFL3D-normalization)
c      xlength0.........reference nondim length (real)
c      tinf.............reference temperature in Rankine (real)
c      alphw............angle of attack in deg (real)
c      betaw............angle of yaw in deg (real)
c      u0...............reference nondim x-velocity (real)
c      v0...............reference nondim y-velocity (real)
c      w0...............reference nondim z-velocity (real)
c      ialph............parameter in CFL3D for determining whether y 
c                       or z is "up" (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension units(5)
      character chardata*80
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Define reinverse
      reinverse=xmach/reue
c   Need to alter vk0 in order to have "consistent normalization" in the
c   CGNS file (see SIDS)
      vk0n=vk0*xmach/reue
c   Go to Base level
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create/get DataClass node
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create/get ReferenceState node
      chardata='Reference state variables for nondimensional data'
      call cg_state_write_f(chardata,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write Mach array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Mach',RealDouble,1,1,xmach,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Mach_Velocity',RealDouble,1,1,xmach,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Mach_VelocitySound',RealDouble,1,1,
     +   c0,ier)
        if (ier .ne. 0) call cg_error_exit_f
      else
        call cg_array_write_f('Mach',RealSingle,1,1,xmach,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Mach_Velocity',RealSingle,1,1,xmach,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Mach_VelocitySound',RealSingle,1,1,
     +   c0,ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NondimensionalParameter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write Reynolds array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Reynolds',RealDouble,1,1,reue,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_Velocity',RealDouble,1,1,xmach,
     +   ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_Length',RealDouble,1,1,xlength0,
     +   ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_ViscosityKinematic',RealDouble,
     +   1,1,reinverse,ier)
        if (ier .ne. 0) call cg_error_exit_f
      else
        call cg_array_write_f('Reynolds',RealSingle,1,1,reue,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_Velocity',RealSingle,1,1,xmach,
     +   ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_Length',RealSingle,1,1,xlength0,
     +   ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('Reynolds_ViscosityKinematic',RealSingle,
     +   1,1,reinverse,ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  4,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NondimensionalParameter,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write Temperature array and dataclass & units
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Temperature',RealDouble,1,1,tinf,ier)
      else
        call cg_array_write_f('Temperature',RealSingle,1,1,tinf,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  8,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(Dimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_units_write_f(Null,Null,Null,Rankine,Null,ier)
      if (ier .ne. 0) call cg_error_exit_f
      units(1)=0.
      units(2)=0.
      units(3)=0.
      units(4)=1.
      units(5)=0.
      if (idouble .eq. 1) then
        call cg_exponents_write_f(RealDouble,units,ier)
      else
        call cg_exponents_write_f(RealSingle,units,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write AngleofAttack array and dataclass & units
c   (this is a CFL3D-specific node:  it is not part of the SIDS standard)
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('AngleofAttack',RealDouble,1,1,alphw,ier)
      else
        call cg_array_write_f('AngleofAttack',RealSingle,1,1,alphw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  9,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(Dimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_units_write_f(Null,Null,Null,Null,Degree,ier)
      if (ier .ne. 0) call cg_error_exit_f
      units(1)=0.
      units(2)=0.
      units(3)=0.
      units(4)=0.
      units(5)=1.
      if (idouble .eq. 1) then
        call cg_exponents_write_f(RealDouble,units,ier)
      else
        call cg_exponents_write_f(RealSingle,units,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write AngleofYaw array and dataclass & units
c   (this is a CFL3D-specific node:  it is not part of the SIDS standard)
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('AngleofYaw',RealDouble,1,1,betaw,ier)
      else
        call cg_array_write_f('AngleofYaw',RealSingle,1,1,betaw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  10,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(Dimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_units_write_f(Null,Null,Null,Null,Degree,ier)
      if (ier .ne. 0) call cg_error_exit_f
      units(1)=0.
      units(2)=0.
      units(3)=0.
      units(4)=0.
      units(5)=1.
      if (idouble .eq. 1) then
        call cg_exponents_write_f(RealDouble,units,ier)
      else
        call cg_exponents_write_f(RealSingle,units,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write Density array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Density',RealDouble,1,1,rho0,ier)
      else
        call cg_array_write_f('Density',RealSingle,1,1,rho0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  11,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write VelocitySound array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('VelocitySound',RealDouble,1,1,c0,ier)
      else
        call cg_array_write_f('VelocitySound',RealSingle,1,1,c0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  12,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write ViscosityMolecular array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('ViscosityMolecular',RealDouble,
     +   1,1,vk0n,ier)
      else
        call cg_array_write_f('ViscosityMolecular',RealSingle,
     +   1,1,vk0n,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  13,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write ViscosityMolecularCFL3D array and dataclass
c   (this is a CFL3D-specific node:  it is not part of the SIDS standard)
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('ViscosityMolecularCFL3D',RealDouble,
     +   1,1,vk0,ier)
      else
        call cg_array_write_f('ViscosityMolecularCFL3D',RealSingle,
     +   1,1,vk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  14,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_descriptor_write_f('Explanation',
     + 'not consistently-normalized (as defined in SIDS)',ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write LengthReference array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('LengthReference',RealDouble,1,1,
     +   xlength0,ier)
      else
        call cg_array_write_f('LengthReference',RealSingle,1,1,
     +   xlength0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  15,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write Pressure array and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Pressure',RealDouble,1,1,p0,ier)
      else
        call cg_array_write_f('Pressure',RealSingle,1,1,p0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  16,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Go to ReferenceState node, write VelocityX, Y, and Z arrays and dataclass
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('VelocityX',RealDouble,1,1,u0,ier)
      else
        call cg_array_write_f('VelocityX',RealSingle,1,1,u0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  17,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
        if (idouble .eq. 1) then
          call cg_array_write_f('VelocityY',RealDouble,1,1,v0,ier)
        else
          call cg_array_write_f('VelocityY',RealSingle,1,1,v0,ier)
        end if
      else
        if (idouble .eq. 1) then
          call cg_array_write_f('VelocityY',RealDouble,1,1,w0,ier)
        else
          call cg_array_write_f('VelocityY',RealSingle,1,1,w0,ier)
        end if
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  18,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
        if (idouble .eq. 1) then
          call cg_array_write_f('VelocityZ',RealDouble,1,1,w0,ier)
        else
          call cg_array_write_f('VelocityZ',RealSingle,1,1,w0,ier)
        end if
      else
        if (idouble .eq. 1) then
          call cg_array_write_f('VelocityZ',RealDouble,1,1,-v0,ier)
        else
          call cg_array_write_f('VelocityZ',RealSingle,1,1,-v0,ier)
        end if
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  19,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_write_f(NormalizedByUnknownDimensional,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
c
      subroutine readref(iccg,ibase,xmach,reue,
     +                   rho0,c0,vk0,xlength0,tinf,alphw)
c**********************************************************************
c     Purpose: Reads reference (nondimensionalization) info from
c     CGNS file. If the DataClass is not right 
c     (NormalizedByUnknownDimensional), the routine stops. Doesn't 
c     bother to read Pressure, AngleofYaw, VelocityX, VelocityY, or
c     VelocityZ because they are not needed. Tries to read 
c     ViscosityMolecularCFL3D (CFL3D-type-normalization) first - if
c     this does not exist, then reads from ViscosityMolecular (SIDS 
c     standard)
c      (NOTE:  cg_array_read_as_f may eventually be converted
c      to cg_array_read_f)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c   OUTPUTS:
c      xmach............Mach number (real)
c      reue.............Reynolds number (real)
c      rho0.............reference nondim density (real)
c      c0...............reference nondim speed of sound (real)
c      vk0..............reference nondim molecular viscosity (real) 
c                       (CFL3D-normalization)
c      xlength0.........reference nondim length (real)
c      tinf.............reference temperature in Rankine (real)
c      alphw............angle of attack in deg (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=19)
c
      dimension idimvec(4)
      character namea*80,name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Go to Base level
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Read DataClass node
      call cg_dataclass_read_f(itype,ier)
      if (ier .gt. 0) then
         write(901,'('' WARNING: DataClass node does not exist.'')')
         write(901,'(''   Assuming NormalizedByUnknownDimensional'')')
         write(901,'(''   and continuing'')')
         ntr=0
      else
         if(itype .ne. NormalizedByUnknownDimensional) then
           write(901,'('' Error:  Restart file must have data that'',
     +       '' is NormalizedByUnknownDimensional'')')
           write(901,'('' Actual classname is: '',a32)') 
     +      DataClassName(itype)
           stop
         end if
      end if
c
c   Read ReferenceState node
      call cg_state_read_f(namea,ier)
      if (ier .gt. 0) then
         write(901,'('' WARNING: ReferenceState node does not exist'')')
         write(901,'(''   Assuming CFL3D <<standard>> for'')')
         write(901,'(''   nondim quantities, and assuming M=0.5,'')')
         write(901,'(''   Re=1.e6, Tinf=520R, alp=0 and continuing'')')
         rho0=1.
         c0=1.
         vk0=1.
         xlength0=1.
         xmach=0.5
         reue=1.e6
         tinf=520.
         alphw=0.
         return
      end if
c
c   Go to ReferenceState node and get info
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
         write(901,'('' WARNING.  To many ReferenceState arrays.'')')
         write(901,'(''   Assuming CFL3D <<standard>> for'')')
         write(901,'(''   nondim quantities, and assuming M=0.5,'')')
         write(901,'(''   Re=1.e6, Tinf=520R, alph=0 and continuing'')')
         rho0=1.
         c0=1.
         vk0=1.
         xlength0=1.
         xmach=0.5
         reue=1.e6
         tinf=520.
         alphw=0.
         return
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get Mach
      do n=1,narrays
        if (name(n) .eq. 'Mach') goto 101
      enddo
      write(901,'('' WARNING. No Mach node exists'')')
      write(901,'(''   setting arbitrarily to 0.5 & continuing'')')
      xmach=0.5
      goto 201
 101  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,xmach,ier)
      else
        call cg_array_read_as_f(n,RealSingle,xmach,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 201  continue
c
c   Get Reynolds number
      do n=1,narrays
        if (name(n) .eq. 'Reynolds') goto 102
      enddo
      write(901,'('' WARNING. No Reynolds node exists'')')
      write(901,'(''   setting arbitrarily to 1.e6 & continuing'')')
      reue=1.e6
      goto 202
 102  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,reue,ier)
      else
        call cg_array_read_as_f(n,RealSingle,reue,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 202  continue
c
c   Get reference density
      do n=1,narrays
        if (name(n) .eq. 'Density') goto 105
      enddo
      write(901,'('' WARNING. No Density node exists'')')
      write(901,'(''   setting to 1.0 & continuing'')')
      rho0=1.0
      goto 205
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rho0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rho0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 205  continue
c
c   Get reference speed of sound
      do n=1,narrays
        if (name(n) .eq. 'VelocitySound') goto 106
      enddo
      write(901,'('' WARNING. No VelocitySound node exists'')')
      write(901,'(''   setting to 1.0 & continuing'')')
      c0=1.0
      goto 206
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,c0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,c0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 206  continue
c
c   Get reference molecular viscosity
c   try to read CFL3D-normalization first (not SIDS standard)
      iflag1=0
      iflag2=0
      do n=1,narrays
        if (name(n) .eq. 'ViscosityMolecularCFL3D') then
          iflag1=1
          nsav1=n
        end if
        if (name(n) .eq. 'ViscosityMolecular') then
          iflag2=1
          nsav2=n
        end if
      enddo
      if (iflag1.eq.1) then
        write(901,'('' Reading ViscosityMolecular'',
     +   '' (CFL3D-normalization)'')')
        if (idouble .eq. 1) then
          call cg_array_read_as_f(nsav1,RealDouble,vk0,ier)
        else
          call cg_array_read_as_f(nsav1,RealSingle,vk0,ier)
        end if
        if (ier .ne. 0) call cg_error_exit_f
      else if (iflag2.eq.1) then
        write(901,'('' Reading ViscosityMolecular'')')
        if (idouble .eq. 1) then
          call cg_array_read_as_f(nsav2,RealDouble,vk0n,ier)
        else
          call cg_array_read_as_f(nsav2,RealSingle,vk0n,ier)
        end if
        if (ier .ne. 0) call cg_error_exit_f
c   put vk0 back into CFL3D normalization
        vk0=vk0n*reue/xmach
      else
        write(901,'('' WARNING. No ViscosityMolecular nodes exist'')')
        write(901,'(''   setting to 1.0 & continuing'')')
        vk0=1.0
      end if
c
c   Get reference length
      do n=1,narrays
        if (name(n) .eq. 'LengthReference') goto 108
      enddo
      write(901,'('' WARNING. No LengthReference node exists'')')
      write(901,'(''   setting to 1.0 & continuing'')')
      xlength0=1.0
      goto 208
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,xlength0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,xlength0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 208  continue
c
c   Get reference temperature
      do n=1,narrays
        if (name(n) .eq. 'Temperature') goto 110
      enddo
      write(901,'('' WARNING. No Temperature node exists'')')
      write(901,'(''   setting arbitrarily to 520R & continuing'')')
      tinf=520.0
      goto 210
 110  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,tinf,ier)
      else
        call cg_array_read_as_f(n,RealSingle,tinf,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  n,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_read_f(itype,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (itype .ne. Dimensional) then
         write(901,'('' WARNING: Restart file must have Temp data '',
     +    ''that is Dimensional'')')
         write(901,'(''   Actual classname is: '',a32)') 
     +    DataClassName(itype)
         write(901,'(''   Ignoring this, and continuing.'')')
      end if
      call cg_units_read_f(im,il,it,ite,ia,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (ite .ne. Rankine) then
         write(901,'('' WARNING: Temp DimensionalUnits should be'',
     +    '' Rankine.'')')
         write(901,'(''   Instead, they are:'',a32)') 
     +    TemperatureUnitsName(ite)
         write(901,'(''   Resetting to 520R and continuing.'')')
         tinf=520.
      end if
 210  continue
c
c   Get angle of attack
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      do n=1,narrays
        if (name(n) .eq. 'AngleofAttack') goto 410
      enddo
      write(901,'('' WARNING. No AngleofAttack node exists'')')
      write(901,'(''   setting arbitrarily to 0 & continuing'')')
      alphw=0.
      goto 510
 410  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,alphw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,alphw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'ReferenceState_t',1,'DataArray_t',
     +  n,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_dataclass_read_f(itype,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (itype .ne. Dimensional) then
         write(901,'('' WARNING: Restart file must have alpha data '',
     +    ''that is Dimensional'')')
         write(901,'(''   Actual classname is: '',a32)')
     +    DataClassName(itype)
         write(901,'(''   Ignoring this, and continuing.'')')
      end if
      call cg_units_read_f(im,il,it,ite,ia,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (ia .ne. Degree) then
         write(901,'('' WARNING: alpha DimensionalUnits should be'',
     +    '' Degree.'')')
         write(901,'(''   Instead, they are:'',a32)')
     +    AngleUnitsName(ia)
         write(901,'(''   Resetting to 0 and continuing.'')')
         alphw=0.
      end if
 510  continue
      return
      end
c
      subroutine writeinput(iccg,ibase,iunit5,irest)
c**********************************************************************
c     Purpose: Writes a "copy" of the CFL3D input file to CGNS file,
c     as a descriptor node called "InputFileUsed" under the Base.
c     (Currently limited to a maximum of 2969 lines)
c     Also writes the case title to a separate descriptor node called
c     "CaseTitle"
c   
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      iunit5...........unit number of the input file to copy (integer)
c      irest............flag is zero if not a restart run (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (inummax=296)
c   stx must be dimensioned 81*(inummax+1), stxsav 81*(inummax+1)*10
      character s(inummax)*80, stx(10)*24057, dum*80, dum1*1, extra*18, 
     .          title*98, stxsav*240570
      character name*32, filenum*2, title2*32
c
c   Use ASCII code 10 (line feed)
      i=10
c
c   Go to Base level
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Determine existing names of descriptor nodes
c   (we want to keep adding new InputFileUsedXX nodes)
c   Currently Max of 50 allowed - beyond this, the last one keeps getting overwritten
      call cg_ndescriptors_f(ndesc,ier)
      if (ier .ne. 0) call cg_error_exit_f
      inum=0
      do n=ndesc,1,-1
        call cg_descriptor_read_f(n,name,stxsav,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (name.eq.'InputFileUsed 1'.or.name.eq.'InputFileUsed 2'.or.
     .      name.eq.'InputFileUsed 3'.or.name.eq.'InputFileUsed 4'.or.
     .      name.eq.'InputFileUsed 5'.or.name.eq.'InputFileUsed 6'.or.
     .      name.eq.'InputFileUsed 7'.or.name.eq.'InputFileUsed 8'.or.
     .      name.eq.'InputFileUsed 9'.or.name.eq.'InputFileUsed10'.or.
     .      name.eq.'InputFileUsed11'.or.name.eq.'InputFileUsed12'.or.
     .      name.eq.'InputFileUsed13'.or.name.eq.'InputFileUsed14'.or.
     .      name.eq.'InputFileUsed15'.or.name.eq.'InputFileUsed16'.or.
     .      name.eq.'InputFileUsed17'.or.name.eq.'InputFileUsed18'.or.
     .      name.eq.'InputFileUsed19'.or.name.eq.'InputFileUsed20'.or.
     .      name.eq.'InputFileUsed21'.or.name.eq.'InputFileUsed22'.or.
     .      name.eq.'InputFileUsed23'.or.name.eq.'InputFileUsed24'.or.
     .      name.eq.'InputFileUsed25'.or.name.eq.'InputFileUsed26'.or.
     .      name.eq.'InputFileUsed27'.or.name.eq.'InputFileUsed28'.or.
     .      name.eq.'InputFileUsed29'.or.name.eq.'InputFileUsed30'.or.
     .      name.eq.'InputFileUsed31'.or.name.eq.'InputFileUsed32'.or.
     .      name.eq.'InputFileUsed33'.or.name.eq.'InputFileUsed34'.or.
     .      name.eq.'InputFileUsed35'.or.name.eq.'InputFileUsed36'.or.
     .      name.eq.'InputFileUsed37'.or.name.eq.'InputFileUsed38'.or.
     .      name.eq.'InputFileUsed39'.or.name.eq.'InputFileUsed40'.or.
     .      name.eq.'InputFileUsed41'.or.name.eq.'InputFileUsed42'.or.
     .      name.eq.'InputFileUsed43'.or.name.eq.'InputFileUsed44'.or.
     .      name.eq.'InputFileUsed45'.or.name.eq.'InputFileUsed46'.or.
     .      name.eq.'InputFileUsed47'.or.name.eq.'InputFileUsed48'.or.
     .      name.eq.'InputFileUsed49'.or.name.eq.'InputFileUsed50') then
          if (irest .eq. 0) then
            call cg_delete_node_f(name,ier)
            if (ier .ne. 0) call cg_error_exit_f
          else
            inum=inum+1
          end if
        end if
      enddo
      write(filenum,102) inum+1
 102  format(i2)
      title2='InputFileUsed' // filenum
c
c   Re-read input file
      rewind(iunit5)
      idvalssav=0
      do m=1,10
        stx(m)=''
      enddo
      do 100 m=1,10
      iend=0
      do n=1,inummax
        s(n)=''
      enddo
      do n=1,inummax
      read(iunit5,'(a80)',end=1999) s(n)
      enddo
      n=n+1
      read(iunit5,'(a80)',end=1999) dum
      iend=-1
 1999 continue
      if(iend .eq. 0) then
        dum=''
      else if(iend .eq. -1 .and. m .eq. 10) then
        dum='   ... input file was cut off here ...'
      end if
      nfinal=n-1
      idvals=idvalssav+81*(nfinal+1)
      stx(m)=s(  1)//char(i)//s(  2)//char(i)//s(  3)//char(i)
     + // s(  4)//char(i)//s(  5)//char(i)//s(  6)//char(i)
     + // s(  7)//char(i)//s(  8)//char(i)//s(  9)//char(i)
     + // s( 10)//char(i)//s( 11)//char(i)//s( 12)//char(i)
     + // s( 13)//char(i)//s( 14)//char(i)//s( 15)//char(i)
     + // s( 16)//char(i)//s( 17)//char(i)//s( 18)//char(i)
     + // s( 19)//char(i)//s( 20)//char(i)//s( 21)//char(i)
     + // s( 22)//char(i)//s( 23)//char(i)//s( 24)//char(i)
     + // s( 25)//char(i)//s( 26)//char(i)//s( 27)//char(i)
     + // s( 28)//char(i)//s( 29)//char(i)//s( 30)//char(i)
     + // s( 31)//char(i)//s( 32)//char(i)//s( 33)//char(i)
     + // s( 34)//char(i)//s( 35)//char(i)//s( 36)//char(i)
     + // s( 37)//char(i)//s( 38)//char(i)//s( 39)//char(i)
     + // s( 40)//char(i)//s( 41)//char(i)//s( 42)//char(i)
     + // s( 43)//char(i)//s( 44)//char(i)//s( 45)//char(i)
     + // s( 46)//char(i)//s( 47)//char(i)//s( 48)//char(i)
     + // s( 49)//char(i)//s( 50)//char(i)//s( 51)//char(i)
     + // s( 52)//char(i)//s( 53)//char(i)//s( 54)//char(i)
     + // s( 55)//char(i)//s( 56)//char(i)//s( 57)//char(i)
     + // s( 58)//char(i)//s( 59)//char(i)//s( 60)//char(i)
     + // s( 61)//char(i)//s( 62)//char(i)//s( 63)//char(i)
     + // s( 64)//char(i)//s( 65)//char(i)//s( 66)//char(i)
     + // s( 67)//char(i)//s( 68)//char(i)//s( 69)//char(i)
     + // s( 70)//char(i)//s( 71)//char(i)//s( 72)//char(i)
     + // s( 73)//char(i)//s( 74)//char(i)//s( 75)//char(i)
     + // s( 76)//char(i)//s( 77)//char(i)//s( 78)//char(i)
     + // s( 79)//char(i)//s( 80)//char(i)//s( 81)//char(i)
     + // s( 82)//char(i)//s( 83)//char(i)//s( 84)//char(i)
     + // s( 85)//char(i)//s( 86)//char(i)//s( 87)//char(i)
     + // s( 88)//char(i)//s( 89)//char(i)//s( 90)//char(i)
     + // s( 91)//char(i)//s( 92)//char(i)//s( 93)//char(i)
     + // s( 94)//char(i)//s( 95)//char(i)//s( 96)//char(i)
     + // s( 97)//char(i)//s( 98)//char(i)//s( 99)//char(i)
     + // s(100)//char(i)
     + // s(101)//char(i)//s(102)//char(i)//s(103)//char(i)
     + // s(104)//char(i)//s(105)//char(i)//s(106)//char(i)
     + // s(107)//char(i)//s(108)//char(i)//s(109)//char(i)
     + // s(110)//char(i)//s(111)//char(i)//s(112)//char(i)
     + // s(113)//char(i)//s(114)//char(i)//s(115)//char(i)
     + // s(116)//char(i)//s(117)//char(i)//s(118)//char(i)
     + // s(119)//char(i)//s(120)//char(i)//s(121)//char(i)
     + // s(122)//char(i)//s(123)//char(i)//s(124)//char(i)
     + // s(125)//char(i)//s(126)//char(i)//s(127)//char(i)
     + // s(128)//char(i)//s(129)//char(i)//s(130)//char(i)
     + // s(131)//char(i)//s(132)//char(i)//s(133)//char(i)
     + // s(134)//char(i)//s(135)//char(i)//s(136)//char(i)
     + // s(137)//char(i)//s(138)//char(i)//s(139)//char(i)
     + // s(140)//char(i)//s(141)//char(i)//s(142)//char(i)
     + // s(143)//char(i)//s(144)//char(i)//s(145)//char(i)
     + // s(146)//char(i)//s(147)//char(i)//s(148)//char(i)
     + // s(149)//char(i)//s(150)//char(i)//s(151)//char(i)
     + // s(152)//char(i)//s(153)//char(i)//s(154)//char(i)
     + // s(155)//char(i)//s(156)//char(i)//s(157)//char(i)
     + // s(158)//char(i)//s(159)//char(i)//s(160)//char(i)
     + // s(161)//char(i)//s(162)//char(i)//s(163)//char(i)
     + // s(164)//char(i)//s(165)//char(i)//s(166)//char(i)
     + // s(167)//char(i)//s(168)//char(i)//s(169)//char(i)
     + // s(170)//char(i)//s(171)//char(i)//s(172)//char(i)
     + // s(173)//char(i)//s(174)//char(i)//s(175)//char(i)
     + // s(176)//char(i)//s(177)//char(i)//s(178)//char(i)
     + // s(179)//char(i)//s(180)//char(i)//s(181)//char(i)
     + // s(182)//char(i)//s(183)//char(i)//s(184)//char(i)
     + // s(185)//char(i)//s(186)//char(i)//s(187)//char(i)
     + // s(188)//char(i)//s(189)//char(i)//s(190)//char(i)
     + // s(191)//char(i)//s(192)//char(i)//s(193)//char(i)
     + // s(194)//char(i)//s(195)//char(i)//s(196)//char(i)
     + // s(197)//char(i)//s(198)//char(i)//s(199)//char(i)
     + // s(200)//char(i)
     + // s(201)//char(i)//s(202)//char(i)//s(203)//char(i)
     + // s(204)//char(i)//s(205)//char(i)//s(206)//char(i)
     + // s(207)//char(i)//s(208)//char(i)//s(209)//char(i)
     + // s(210)//char(i)//s(211)//char(i)//s(212)//char(i)
     + // s(213)//char(i)//s(214)//char(i)//s(215)//char(i)
     + // s(216)//char(i)//s(217)//char(i)//s(218)//char(i)
     + // s(219)//char(i)//s(220)//char(i)//s(221)//char(i)
     + // s(222)//char(i)//s(223)//char(i)//s(224)//char(i)
     + // s(225)//char(i)//s(226)//char(i)//s(227)//char(i)
     + // s(228)//char(i)//s(229)//char(i)//s(230)//char(i)
     + // s(231)//char(i)//s(232)//char(i)//s(233)//char(i)
     + // s(234)//char(i)//s(235)//char(i)//s(236)//char(i)
     + // s(237)//char(i)//s(238)//char(i)//s(239)//char(i)
     + // s(240)//char(i)//s(241)//char(i)//s(242)//char(i)
     + // s(243)//char(i)//s(244)//char(i)//s(245)//char(i)
     + // s(246)//char(i)//s(247)//char(i)//s(248)//char(i)
     + // s(249)//char(i)//s(250)//char(i)//s(251)//char(i)
     + // s(252)//char(i)//s(253)//char(i)//s(254)//char(i)
     + // s(255)//char(i)//s(256)//char(i)//s(257)//char(i)
     + // s(258)//char(i)//s(259)//char(i)//s(260)//char(i)
     + // s(261)//char(i)//s(262)//char(i)//s(263)//char(i)
     + // s(264)//char(i)//s(265)//char(i)//s(266)//char(i)
     + // s(267)//char(i)//s(268)//char(i)//s(269)//char(i)
     + // s(270)//char(i)//s(271)//char(i)//s(272)//char(i)
     + // s(273)//char(i)//s(274)//char(i)//s(275)//char(i)
     + // s(276)//char(i)//s(277)//char(i)//s(278)//char(i)
     + // s(279)//char(i)//s(280)//char(i)//s(281)//char(i)
     + // s(282)//char(i)//s(283)//char(i)//s(284)//char(i)
     + // s(285)//char(i)//s(286)//char(i)//s(287)//char(i)
     + // s(288)//char(i)//s(289)//char(i)//s(290)//char(i)
     + // s(291)//char(i)//s(292)//char(i)//s(293)//char(i)
     + // s(294)//char(i)//s(295)//char(i)//s(296)//char(i)//dum
c
      if (iend .eq. 0) goto 101
      idvalssav=idvals
 100  continue
 101  continue
      stxsav=stx(1)//char(i)//stx(2)//char(i)//stx(3)//char(i)//
     .       stx(4)//char(i)//stx(5)//char(i)//stx(6)//char(i)//
     .       stx(7)//char(i)//stx(8)//char(i)//stx(9)//char(i)//
     .       stx(10)
c   now write data
      if (idvals .le. 10000) then
        call wout1(title2,stxsav)
      else if (idvals .gt.  10000 .and. idvals .le.  20000) then
        call wout2(title2,stxsav)
      else if (idvals .gt.  20000 .and. idvals .le.  30000) then
        call wout3(title2,stxsav)
      else if (idvals .gt.  30000 .and. idvals .le.  40000) then
        call wout4(title2,stxsav)
      else if (idvals .gt.  40000 .and. idvals .le.  50000) then
        call wout5(title2,stxsav)
      else if (idvals .gt.  50000 .and. idvals .le.  60000) then
        call wout6(title2,stxsav)
      else if (idvals .gt.  60000 .and. idvals .le.  80000) then
        call wout7(title2,stxsav)
      else if (idvals .gt.  80000 .and. idvals .le. 100000) then
        call wout8(title2,stxsav)
      else if (idvals .gt. 100000 .and. idvals .le. 120000) then
        call wout9(title2,stxsav)
      else if (idvals .gt. 120000 .and. idvals .le. 140000) then
        call wout10(title2,stxsav)
      else if (idvals .gt. 140000 .and. idvals .le. 160000) then
        call wout11(title2,stxsav)
      else if (idvals .gt. 160000 .and. idvals .le. 200000) then
        call wout12(title2,stxsav)
      else if (idvals .gt. 200000 .and. idvals .le. 240000) then
        call wout13(title2,stxsav)
      else
        call wout14(title2,stxsav)
      end if
c
c   Create/get CaseTitle node
c
c   Re-read input file
      rewind(iunit5)
      do n=1,14
        read(iunit5,'(a80)') dum
      enddo
      read(iunit5,'(a1)') dum1
      if (dum1 .eq. '>') then
        do n=1,500
          read(iunit5,'(a1)') dum1
          if (dum1 .eq. '<') goto 1102
        enddo
        write(901,'('' CaseTitle in cgns file incorrect because'',
     .    '' more than 500 lines of keyword input'')')
1102    continue
      else
        backspace(iunit5)
      end if
      read(iunit5,'(a80)') dum
c   Tack on this to beginning of title (18 characters):
      extra='CFL3D V6.7 code:  '
      title = extra // dum
      call cg_descriptor_write_f('CaseTitle',title,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
c**********************************************************************
c     The following subroutines are used to write a copy of the input
c     file to a descriptor node.  There are many (rather than one) so
c     that the length of the string will be only somewhat bigger than
c     the actual length of the input file (rather than the max
c     length), making "adfedit" reading easier
c**********************************************************************
c
      subroutine wout1(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*10000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout2(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*20000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout3(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*30000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout4(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*40000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout5(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*50000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout6(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*60000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout7(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*80000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout8(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*100000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout9(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*120000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout10(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*140000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout11(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*160000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout12(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*200000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout13(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*240000
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wout14(title,stx)
#     include "cgnslib_f.h"
      character title*32,stx*240570
      call cg_descriptor_write_f(title,stx,ier)
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine winfot(iccg,ibase,dt)
c**********************************************************************
c     Purpose: Writes out Descriptor node at base level, describing
c     SimulationType 
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      dt...............delta t parameter in CFL3D 
c                       (+ = time accurate, - = SS) (real)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      if (dt .lt. 0.) then
        call cg_simulation_type_write_f(iccg,ibase,NonTimeAccurate,ier)
      else
        call cg_simulation_type_write_f(iccg,ibase,TimeAccurate,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      return
      end
c
      subroutine wsecord(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +                   jdim,kdim,wk,qc0,tursav2,nsoluse,i2d,ialph,
     +                   ivmx,nummem)
c**********************************************************************
c     Purpose: Writes qc0 & tursav2 (Q's and turb quantities
c     from previous time step) for 2nd order
c     time-accurate solns, for zone number "igrid".
c     The CGNS file must be already opened in "MODE_MODIFY" mode.
c     Writes results in the following order:
c
c        (idim,jdim,kdim) <- necessary for CGNS!,
c
c     even though they are stored in CFL3D in the order: jdim,kdim,idim
c     Also writes rind cell information: No B.C info exists for these.
c     Currently we output nearby values for these.
c     For 2-D solutions, rind cells are not written in the i-direction.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idima,jdima,kdima...expected dimensions of this zone (zone 
c                          "igrid"), at finest level (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) - if actual dimensions are a  
c                       subset of idima,jdima,kdima, results will be
c                       prolonged to the fine grid and written at that
c                       level (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      qc0..............cell-centered q-values, in (j,k,i) order,
c                       PRIMITIVE quantities (real)
c      tursav2..........cell-centered turb values from last time step, 
c                       in (j,k,i) order; dimensioned by 2*nummem, but only
c                       need to store nummem of them (real)
c      nsoluse..........FlowSolution index number where data to be
c                       stored (integer)
c      i2d..............0 if 3-D (rind data to be output in all 3
c                          directions),
c                       1 if 2-D (rind data to be output only in j-
c                          and k- directions) (integer)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c      ivmx............ turb model number for CFL3D (integer)
c      nummem...........one of dimensions for tursav2 (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension qc0(jdim,kdim,idim-1,5),
     .          tursav2(jdim,kdim,idim,2*nummem),
     .          wk(i2d:idima-i2d,0:jdima,0:kdima)
      character*32 name
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   First, get "factor" in case want to write every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
c
c   Density:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,1,2,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'DensityLastDT', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'DensityLastDT', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   U-velocity:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,2,2,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'VelocityXLastDT', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'VelocityXLastDT', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   V-velocity:
      if (ialph .eq. 0) then
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,3,2,wk)
      else
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,4,2,wk)
      end if
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'VelocityYLastDT', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'VelocityYLastDT', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   W-velocity:
      if (ialph .eq. 0) then
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,4,2,wk)
      else
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,3,2,wk)
      do i=i2d,idima-i2d
        do j=0,jdima
          do k=0,kdima
            wk(i,j,k)=-wk(i,j,k)
          enddo
        enddo
      enddo
      end if
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'VelocityZLastDT', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'VelocityZLastDT', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Pressure:
      call reordsln(nfac,i2d,idima,jdima,kdima,idim-1,jdim,kdim,
     +  qc0,wk,wk,wk,wk,wk,wk,5,5,2,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      'PressureLastDT', wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      'PressureLastDT', wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Turb quantities:
      if (ivmx .ge. 4) then
      do l=1,nummem
        if (l .eq. 1) name='TurbVal1LastDT'
        if (l .eq. 2) name='TurbVal2LastDT'
        if (l .eq. 3) name='TurbVal3LastDT'
        if (l .eq. 4) name='TurbVal4LastDT'
        if (l .eq. 5) name='TurbVal5LastDT'
        if (l .eq. 6) name='TurbVal6LastDT'
        if (l .eq. 7) name='TurbVal7LastDT'
      call reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  tursav2,wk,wk,wk,wk,wk,wk,4,l,2,wk)
      if (idouble .eq. 1) then
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealDouble,
     +                      name, wk, ifindex, ier)
      else
      call cg_field_write_f(iccg, ibase, igrid, nsoluse, RealSingle,
     +                      name, wk, ifindex, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      enddo
      end if
c
      return
      end
c
      subroutine rsecord(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +                   jdim,kdim,wk,nsoluse,irind,jrind,krind,i2d,dt,
     +                   ialph,qc0,tursav2,dtold,ivmx,nummem)
c**********************************************************************
c     Purpose: Reads qc0 & tursav2 (Q's and turb quantities
c     from previous time step) & dt for 2nd order
c     time-accurate solns, for zone number "igrid". Must know (through 
c     irind,jrind,krind) what Rind Cell info exists.
c     Gets results in the following order:
c
c       (jdim,kdim,idim) <- necessary for CFL3D!,
c
c     and assumes they are stored in the CGNS database in the order:
c     (idim,jdim,kdim)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idima,jdima,kdima...expected (GRIDPOINT) dimensions of this zone
c                          (zone "igrid"), at finest level (integers)
c      idim,jdim,kdim...desired (GRIDPOINT) dimensions of this zone
c                       (e.g., could be every other point of fine grid)
c                       (integers)
c      wk...............working space needed, of dimension
c                       (idima+1)*(jdima+1)*(kdima+1) or larger (real)
c      nsoluse..........CGNS "FlowSolution" index number (integer)
c      irind,jrind,krind...= 0 if rind data exists for the direction,
c                          = 1 if it does not (integers)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      dt...............current time step (real)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c      ivmx............ turb model number for CFL3D (integer)
c      nummem...........one of dimensions for tursav2 (integer)
c   OUTPUTS:
c      qc0..............cell-centered q-values, in (j,k,i) order,
c                       PRIMITIVE quantities (real)
c      tursav2..........cell-centered turb values from last time step,
c                       in (j,k,i) order; dimensioned by 2*nummem, but only
c                       need to store numem of them (real)
c      dtold............previous (saved) time step (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter(numfield=50)
c
      dimension qc0(jdim,kdim,idim-1,5),
     +   tursav2(jdim,kdim,idim,2*nummem),
     +   wk(idima-2*irind+1,jdima-2*jrind+1,kdima-2*krind+1)
      dimension irmin(3),irmax(3)
      character*32 fieldname(numfield),fndesired(numfield)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c   Find out how many fields there are and their names:
      call cg_nfields_f(iccg, ibase, igrid, nsoluse, nfields, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nfields .gt. numfield) then
        write(901,'('' Need to increase numfield parameter'',
     +    '' in rsecord'')')
        stop
      end if
      do n=1,nfields
      call cg_field_info_f(iccg, ibase, igrid, nsoluse, n, itype,
     +  fieldname(n), ier)
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Check to make sure the field names needed exist:
      numwanted=5
      fndesired(1)='DensityLastDT'
      fndesired(2)='VelocityXLastDT'
      fndesired(3)='VelocityYLastDT'
      fndesired(4)='VelocityZLastDT'
      fndesired(5)='PressureLastDT'
      ibad=0
      do n=1,numwanted
        ido=0
        do nn=1,nfields
          if(fndesired(n) .eq. fieldname(nn)) then
            ido=1
          end if
        enddo
        if (ido .eq. 0) then
          write(901,'('' Cannot find solution of name '',a32)')
     +      fndesired(n)
          ibad=1
        end if
      enddo
      if (ibad .eq. 1) then
        stop
      end if
c   Read Q's
c   First, get "factor" in case want to read every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level of soln not'',
     +   '' supported'')')
        stop
      end if
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
c   Set up array bounds:
      irmin(1)=1
      irmin(2)=1
      irmin(3)=1
      irmax(1)=idima-2*irind+1
      irmax(2)=jdima-2*jrind+1
      irmax(3)=kdima-2*krind+1
c   Read Q's:
      do m=1,5
      if (idouble .eq. 1) then
      call cg_field_read_f(iccg, ibase, igrid, nsoluse,
     +     fndesired(m), RealDouble, irmin, irmax, wk, ier)
      else
      call cg_field_read_f(iccg, ibase, igrid, nsoluse,
     +     fndesired(m), RealSingle, irmin, irmax, wk, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (m .eq. 3 .and. ialph .ne. 0) then
      call reorderq(idima,jdima,kdima,idim-1,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,qc0(1,1,1,4))
      else if (m .eq. 4 .and. ialph .ne. 0) then
      do i=1,idima-2*irind+1
        do j=1,jdima-2*jrind+1
          do k=1,kdima-2*krind+1
            wk(i,j,k)=-wk(i,j,k)
          enddo
        enddo
      enddo
      call reorderq(idima,jdima,kdima,idim-1,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,qc0(1,1,1,3))
      else
      call reorderq(idima,jdima,kdima,idim-1,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,qc0(1,1,1,m))
      end if
      enddo
c
c   Now get 2nd order turb quantities (they may not exist if this
c   is reading an older file, or if last run was not turbulent)
      if (ivmx .ge. 4) then
      fndesired(1)='TurbVal1LastDT'
      fndesired(2)='TurbVal2LastDT'
      fndesired(3)='TurbVal3LastDT'
      fndesired(4)='TurbVal4LastDT'
      fndesired(5)='TurbVal5LastDT'
      fndesired(6)='TurbVal6LastDT'
      fndesired(7)='TurbVal7LastDT'
      ibad=0
      do n=1,nummem
        ido=0
        do nn=1,nfields
          if(fndesired(n) .eq. fieldname(nn)) then
            ido=1
          end if
        enddo
        if (ido .eq. 0) then
          write(901,'('' Cannot find solution of name '',a32)')
     +      fndesired(n)
          write(901,'('' ...leaving old turbs as 0 and continuing'')')
          ibad=1
        end if
      enddo
      if (ibad .eq. 1) then
c       Do not read - values will be left as zero (as initialized)
        continue
      else
c   Read old turb quantities:
        do m=1,nummem
        if (idouble .eq. 1) then
        call cg_field_read_f(iccg, ibase, igrid, nsoluse,
     +     fndesired(m), RealDouble, irmin, irmax, wk, ier)
        else
        call cg_field_read_f(iccg, ibase, igrid, nsoluse,
     +     fndesired(m), RealSingle, irmin, irmax, wk, ier)
        end if
        if (ier .ne. 0) call cg_error_exit_f
        call reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,tursav2(1,1,1,m))
        enddo
      end if
      end if
c
c   Go to Base node, read CFL3DTimeStep array from UserDefinedData
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .eq. 0) then
        write(901,'('' no user data exists in base, so'',
     +    '' no old time step is available: set to current'')')
        dtold=dt
      else if (nuserdata .gt. numfield) then
        write(901,'('' Error... numfield too small in rsecord'')')
        stop
      else
        do n=1,nuserdata
          call cg_user_data_read_f(n,fieldname(n),ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
        do n=1,nuserdata
          if(fieldname(n) .eq. 'CFL3DTimeStep') goto 101
        enddo
        write(901,'('' WARNING. No CFL3DTimeStep node exists'')')
        write(901,'(''   setting dtold=dt and continuing'')')
        dtold=dt
        goto 201
 101    continue
        call cg_goto_f(iccg,ibase,ier,'UserDefinedData_t',n,'end')
        if (ier .ne. 0) call cg_error_exit_f
        if (idouble .eq. 1) then
          call cg_array_read_as_f(1,RealDouble,dtold,ier)
        else
          call cg_array_read_as_f(1,RealSingle,dtold,ier)
        end if
        if (ier .ne. 0) call cg_error_exit_f
      end if
 201  continue
      return
      end
c
      subroutine reorderg(nfac,nfaci,idima,jdima,kdima,
     .   idim,jdim,kdim,wk,xyz)
c**********************************************************************
c     Purpose: Given wk array, reorder i,j,k to j,k,i for grid
c     (this routine contains no CGNS calls - it merely manipulates
c     the data).
c
c   INPUTS:
c      nfac.............factor denoting "level" of grid (e.g.,
c                       1=grid on finest level, 2=1-level down,
c                       4=2-levels down, 8=3-levels down,
c                       16=4-levels down (integer)
c      nfaci............same as nfac, except for i-direction only
c                       (when 2-D, the i-direction does not coarsen)
c                       (integer)
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid") at finest level (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) (integers)
c      wk...............grid coordinates array (usually x,
c                       y, or z) in i,j,k-ordering, as stored
c                       in CGNS file (real)
c   OUTPUTS:
c      xyz..............grid coordinates array (usually x,
c                       y, or z) in j,k,i-ordering (real)
c**********************************************************************
c
      dimension xyz(jdim,kdim,idim),wk(idima,jdima,kdima)
c
      kk=0
      do k=1,kdima,nfac
      kk=kk+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
        xyz(jj,kk,ii)=wk(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
c
      subroutine reorderq(idima,jdima,kdima,idim,jdim,kdim,nfac,
     +    nfaci,wk,irind,jrind,krind,q)
c**********************************************************************
c     Purpose: given wk array, reorder i,j,k to j,k,i for q
c     (this routine contains no CGNS calls - it merely manipulates
c     the data).
c
c   INPUTS:
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid") at finest level, 
c                          modified by rind storage (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) - if actual dimensions are a
c                       subset idima,jdima,kdima, results will be
c                       prolonged to the fine grid and written at that
c                       level (integers)
c      nfac.............factor denoting "level" of solution (e.g.,
c                       1=soln on finest level, 2=1-level down,
c                       4=2-levels down, 8=3-levels down,
c                       16=4-levels down (integer)
c      nfaci............same as nfac, except for i-direction only
c                       (when 2-D, the i-direction does not coarsen)
c                       (integer)
c      wk...............solution array in i,j,k-ordering, as stored 
c                       in CGNS file (real)
c        
c      irind,jrind,krind...= 0 if rind data exists for the direction,
c                          = 1 if it does not (integers)
c   OUTPUTS:
c      q................solution array in j,k,i ordering (real)
c**********************************************************************
c
      dimension q(jdim,kdim,idim)
      dimension wk(idima-2*irind+1,jdima-2*jrind+1,kdima-2*krind+1)
c
      kk=0
      do k=1,kdima-1,nfac
      kk=kk+1
      jj=0
      do j=1,jdima-1,nfac
      jj=jj+1
      ii=0
      do i=1,idima-1,nfaci
      ii=ii+1
        q(jj,kk,ii)=wk(i+1-irind,j+1-jrind,k+1-krind)
      enddo
      enddo
      enddo
c
      return
      end
c
      subroutine reordsln(nfac,i2d,idima,jdima,kdima,idim,jdim,kdim,
     +  q,bci,bcj,bck,qi0,qj0,qk0,ldim,luse,iflag,wk)
c**********************************************************************
c     Purpose: Given q array, reorder solution j,k,i to i,j,k (include
c     factor if on coarser level, and include rind (BC) data).
c     Note:  when iflag=1, bci,bcj,bck are not accessed, and when
c     iflag=2, bci,bcj,bck,qi0,qj0,qk0 are not accessed.
c
c   INPUTS:
c      nfac.............factor denoting "level" of solution (e.g.,
c                       1=soln on finest level, 2=1-level down,
c                       4=2-levels down, 8=3-levels down,
c                       16=4-levels down (integer)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid") at finest level,
c                          modified by i2d (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) - if actual dimensions are a
c                       subset idima,jdima,kdima, results will be
c                       prolonged to the fine grid and written at that
c                       level (integers)
c      q................solution array in j,k,i ordering (real)
c      bci,bcj,bck......arrays denoting where BCs stored (0=at cell
c                       center, 1=at face center) (real)
c      qi0,qj0,qk0......arrays of BC data (real)
c      ldim.............last dimension of q-array (integer)
c      luse.............which last-dimension of q-array to use 
c                       (integer)
c      iflag............=0 means std BCs used (keys off of bci,bcj,bck)
c                       =1 means BC is at ghost cell
c                       =2 means no BC exists for rind cell (fills with 
c                        nearby values)
c   OUTPUTS:
c      wk...............solution array in i,j,k-ordering, as stored
c                       in CGNS file (real)
c**********************************************************************
c
      dimension wk(i2d:idima-i2d,0:jdima,0:kdima)
      dimension q(jdim,kdim,idim,ldim)
      dimension qi0(jdim,kdim,ldim,4),qj0(kdim,idim-1,ldim,4),
     +          qk0(jdim,idim-1,ldim,4),bci(jdim,kdim,2),
     +          bcj(kdim,idim-1,2),bck(jdim,idim-1,2)
c
c   Set up nfaci parameter
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
c   Reorder j,k,i to i,j,k:
c   if writing every other pt, etc, actually write it out for finest grid
c   (data is repeated as if it is on a coarser grid)
      do n=1,nfac
        kk=0
        do k=n,kdima-1,nfac
          kk=kk+1
          do m=1,nfac
            jj=0
            do j=m,jdima-1,nfac
              jj=jj+1
              do l=1,nfaci
                ii=0
                do i=l,idima-1,nfaci
                  ii=ii+1
                  wk(i,j,k)=q(jj,kk,ii,luse)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c   Load nearby values for rind cells if iflag=2
      if(iflag .eq. 2) then
        if (i2d .ne. 1) then
        i=0
        do j=1,jdima-1
        do k=1,kdima-1
          wk(i,j,k)=wk(i+1,j,k)
        enddo
        enddo
        do j=1,jdima-1
          wk(i,j,0)=wk(i+1,j,1)
          wk(i,j,kdima)=wk(i+1,j,kdima-1)
        enddo
        do k=1,kdima-1
          wk(i,0,k)=wk(i+1,1,k)
          wk(i,jdima,k)=wk(i+1,jdima-1,k)
        enddo
        wk(i,0,0)=wk(i+1,1,1)
        wk(i,jdima,0)=wk(i+1,jdima-1,1)
        wk(i,0,kdima)=wk(i+1,1,kdima-1)
        wk(i,jdima,kdima)=wk(i+1,jdima-1,kdima-1)
        i=idima
        do j=1,jdima-1
        do k=1,kdima-1
          wk(i,j,k)=wk(i-1,j,k)
        enddo
        enddo
        do k=1,kdima-1
          wk(i,0,k)=wk(i-1,1,k)
          wk(i,jdima,k)=wk(i-1,jdima-1,k)
        enddo
        do j=1,jdima-1
          wk(i,j,0)=wk(i-1,j,1)
          wk(i,j,kdima)=wk(i-1,j,kdima-1)
        enddo
        wk(i,0,0)=wk(i-1,1,1)
        wk(i,jdima,0)=wk(i-1,jdima-1,1)
        wk(i,0,kdima)=wk(i-1,1,kdima-1)
        wk(i,jdima,kdima)=wk(i-1,jdima-1,kdima-1)
        end if
c
        j=0
        do i=1,idima-1
        do k=1,kdima-1
          wk(i,j,k)=wk(i,j+1,k)
        enddo
        enddo
        do i=1,idima-1
          wk(i,j,0)=wk(i,j+1,1)
          wk(i,j,kdima)=wk(i,j+1,kdima-1)
        enddo
        if (i2d .ne. 1) then
          do k=1,kdima-1
            wk(0,j,k)=wk(1,j+1,k)
            wk(idima,j,k)=wk(idima-1,j+1,k)
          enddo
          wk(0,j,0)=wk(1,j+1,1)
          wk(idima,j,0)=wk(idima-1,j+1,1)
          wk(0,j,kdima)=wk(1,j+1,kdima-1)
          wk(idima,j,kdima)=wk(idima-1,j+1,kdima-1)
        end if
        j=jdima
        do i=1,idima-1
        do k=1,kdima-1
          wk(i,j,k)=wk(i,j-1,k)
        enddo
        enddo
        do i=1,idima-1
          wk(i,j,0)=wk(i,j-1,1)
          wk(i,j,kdima)=wk(i,j-1,kdima-1)
        enddo
        if (i2d .ne. 1) then
          do k=1,kdima-1
            wk(0,j,k)=wk(1,j-1,k)
            wk(idima,j,k)=wk(idima-1,j-1,k)
          enddo
          wk(0,j,0)=wk(1,j-1,1)
          wk(idima,j,0)=wk(idima-1,j-1,1)
          wk(0,j,kdima)=wk(1,j-1,kdima-1)
          wk(idima,j,kdima)=wk(idima-1,j-1,kdima-1)
        end if
c
        k=0
        do i=1,idima-1
        do j=1,jdima-1
          wk(i,j,k)=wk(i,j,k+1)
        enddo
        enddo
        do i=1,idima-1
          wk(i,0,k)=wk(i,1,k+1)
          wk(i,jdima,k)=wk(i,jdima-1,k+1)
        enddo
        if (i2d .ne. 1) then
          do j=1,jdima-1
            wk(0,j,k)=wk(1,j,k+1)
            wk(idima,j,k)=wk(idima-1,j,k+1)
          enddo
          wk(0,0,k)=wk(1,1,k+1)
          wk(idima,0,k)=wk(idima-1,1,k+1)
          wk(0,jdima,k)=wk(1,jdima-1,k+1)
          wk(idima,jdima,k)=wk(idima-1,jdima-1,k+1)
        end if
        k=kdima
        do i=1,idima-1
        do j=1,jdima-1
          wk(i,j,k)=wk(i,j,k-1)
        enddo
        enddo
        do i=1,idima-1
          wk(i,0,k)=wk(i,1,k-1)
          wk(i,jdima,k)=wk(i,jdima-1,k-1)
        enddo
        if (i2d .ne. 1) then
          do j=1,jdima-1
            wk(0,j,k)=wk(1,j,k-1)
            wk(idima,j,k)=wk(idima-1,j,k-1)
          enddo
          wk(0,0,k)=wk(1,1,k-1)
          wk(idima,0,k)=wk(idima-1,1,k-1)
          wk(0,jdima,k)=wk(1,jdima-1,k-1)
          wk(idima,jdima,k)=wk(idima-1,jdima-1,k-1)
        end if
        return
      end if
c
c   Load BCs into "Rind" regions
c    i0 and idim faces:
      if (i2d .ne. 1) then
      do n=1,nfac
        kk=0
        do k=n,kdima-1,nfac
          kk=kk+1
          do m=1,nfac
            jj=0
            do j=m,jdima-1,nfac
              jj=jj+1
              if(iflag .eq. 1) then
              wk(    0,j,k)=qi0(jj,kk,luse,1)
              wk(idima,j,k)=qi0(jj,kk,luse,3)
              else
              wk(    0,j,k)=(1.-bci(jj,kk,1))*qi0(jj,kk,luse,1)+
     +         bci(jj,kk,1)*(2.*qi0(jj,kk,luse,1)-q(jj,kk,1,luse))
              wk(idima,j,k)=(1.-bci(jj,kk,2))*qi0(jj,kk,luse,3)+
     +         bci(jj,kk,2)*(2.*qi0(jj,kk,luse,3)-q(jj,kk,idima-1,luse))
              end if
            enddo
          enddo
        enddo
      enddo
      end if
c    j0 and jdim faces:
      do n=1,nfac
        kk=0
        do k=n,kdima-1,nfac
          kk=kk+1
          do l=1,nfaci
            ii=0
            do i=l,idima-1,nfaci
              ii=ii+1
              if(iflag .eq. 1) then
              wk(i,    0,k)=qj0(kk,ii,luse,1)
              wk(i,jdima,k)=qj0(kk,ii,luse,3)
              else
              wk(i,    0,k)=(1.-bcj(kk,ii,1))*qj0(kk,ii,luse,1)+
     +         bcj(kk,ii,1)*(2.*qj0(kk,ii,luse,1)-q(1,kk,ii,luse))
              wk(i,jdima,k)=(1.-bcj(kk,ii,2))*qj0(kk,ii,luse,3)+
     +         bcj(kk,ii,2)*(2.*qj0(kk,ii,luse,3)-q(jdima-1,kk,ii,luse))
              end if
            enddo
          enddo
        enddo
      enddo
c    k0 and kdim faces:
      do m=1,nfac
        jj=0
        do j=m,jdima-1,nfac
          jj=jj+1
          do l=1,nfaci
            ii=0
            do i=l,idima-1,nfaci
              ii=ii+1
              if(iflag .eq. 1) then
              wk(i,j,    0)=qk0(jj,ii,luse,1)
              wk(i,j,kdima)=qk0(jj,ii,luse,3)
              else
              wk(i,j,    0)=(1.-bck(jj,ii,1))*qk0(jj,ii,luse,1)+
     +         bck(jj,ii,1)*(2.*qk0(jj,ii,luse,1)-q(jj,1,ii,luse))
              wk(i,j,kdima)=(1.-bck(jj,ii,2))*qk0(jj,ii,luse,3)+
     +         bck(jj,ii,2)*(2.*qk0(jj,ii,luse,3)-q(jj,kdima-1,ii,luse))
              end if
            enddo
          enddo
        enddo
      enddo
c    no unique way to define edge and corner Rind cells:
      do i=1,idima-1
        wk(i,    0,    0)=wk(i,      1,    0)
        wk(i,jdima,    0)=wk(i,jdima-1,    0)
        wk(i,    0,kdima)=wk(i,      1,kdima)
        wk(i,jdima,kdima)=wk(i,jdima-1,kdima)
      enddo
      if (i2d .ne. 1) then
      do k=1,kdima-1
        wk(    0,    0,k)=wk(      1,    0,k)
        wk(    0,jdima,k)=wk(      1,jdima,k)
        wk(idima,    0,k)=wk(idima-1,    0,k)
        wk(idima,jdima,k)=wk(idima-1,jdima,k)
      enddo
      do j=1,jdima-1
        wk(    0,j,    0)=wk(      1,j,    0)
        wk(    0,j,kdima)=wk(      1,j,kdima)
        wk(idima,j,    0)=wk(idima-1,j,    0)
        wk(idima,j,kdima)=wk(idima-1,j,kdima)
      enddo
      wk(    0,    0,    0)=wk(      1,      1,      1)
      wk(idima,    0,    0)=wk(idima-1,      1,      1)
      wk(    0,jdima,    0)=wk(      1,jdima-1,      1)
      wk(    0,    0,kdima)=wk(      1,      1,kdima-1)
      wk(idima,jdima,    0)=wk(idima-1,jdima-1,      1)
      wk(idima,    0,kdima)=wk(idima-1,      1,kdima-1)
      wk(    0,jdima,kdima)=wk(      1,jdima-1,kdima-1)
      wk(idima,jdima,kdima)=wk(idima-1,jdima-1,kdima-1)
      end if
c
      return
      end

      subroutine getiflagg(iccg,ibase,nsoluse,iflagg)
c**********************************************************************
c     Purpose: Checks CGNS file to determine if restart solution is 
c     2nd order in time and/or has a moving grid.
c     Specific to CFL3D:
c     If DensityLastDT exists in any zone, then it means that it is
c     a second order solution.
c     If ANY RigidGridMotion_t nodes exist in CFL3D file, then it 
c     means that it is an unsteady run (there is unsteady data that 
c     needs to be read).
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this
c                       routine) (integer)
c      ibase............CGNS base index number (determined outside this
c                       routine) (integer)
c      nsoluse..........FlowSolution index number where data is stored 
c                       (integer)
c   OUTPUTS:
c      iflagg...........parameter in CFL3D for determining whether 1st
c                       or 2nd order in time and whether moving grid
c                       or not:
c                       0=1st order, non-moving;
c                       1=2nd order, non-moving;
c                       2=1st order, moving;
c                       3=2nd order, moving)
c                       (used only if dt > 0) (integer)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter(numfield=50)
c
      dimension idimvec(4)
      character*32 fieldname(numfield)
c
      iflagg=0
c   Find out how many zones there are
      call cg_nzones_f(iccg, ibase, nzones, ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Check over all zones
      do igrid=1,nzones
      call cg_nfields_f(iccg, ibase, igrid, nsoluse, nfields, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nfields .gt. numfield) then
        write(901,'('' Need to increase numfield parameter'',
     +    '' in getiflagg'')')
        stop
      end if
      do n=1,nfields
      call cg_field_info_f(iccg, ibase, igrid, nsoluse, n, itype,
     +  fieldname(n), ier)
      if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Check to see if any of the fields have the name 'DensityLastDT'
      do n=1,nfields
        if (fieldname(n) .eq. 'DensityLastDT') then
          iflagg=max(iflagg,1)
        end if
      enddo
c   Now check for moving grid info
c   If ANY RigidGridMotion_t nodes exist in CFL3D file -and- if
c   ZoneIterativeData contains RigidGridMotionPointers, then it 
c   it is an unsteady moving grid run (there is unsteady data that  
c   needs to be read)
      call cg_n_rigid_motions_f(iccg, ibase, igrid, nmotions, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nmotions .eq. 0) then
        continue
      else if (nmotions .gt. 1) then
        write(901,'('' Error, more than 1 RigidGridMotion_t nodes!'')')
        stop
      else
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +   'ZoneIterativeData_t',1,'end')
        if (ier .eq. 0) then
c         Find out how many arrays
          call cg_narrays_f(narrays,ier)
          if (ier .ne. 0) call cg_error_exit_f
          if (narrays .gt. numfield) then
            write(901,'('' Error.  Too many RigidGridMotion arrays.'')')
            write(901,'(''   increase numfield in getiflagg'')')
            stop
          end if
          do n=1,narrays
            call cg_array_info_f(n,fieldname(n),itype,idatadim,
     +        idimvec,ier)
            if (ier .ne. 0) call cg_error_exit_f
          enddo
          imov=1
          do n=1,narrays
            if (fieldname(n) .eq. 'RigidGridMotionPointers') goto 101
          enddo
          write(901,'('' No RigidGridMotionPointers exists... assume'',
     +      '' restart from non-moving grid case'')')
          imov=0
 101      continue
          if (imov .eq. 1) then
            if (iflagg .eq. 0) then
              iflagg = 2
            else if (iflagg .eq. 1) then
              iflagg = 3
            end if
          end if
        end if
      end if
      enddo
      return
      end

      subroutine writeziter(iccg,ibase,igrid)
c**********************************************************************
c     Purpose:  Creates 'ZoneIterativeData' node, and write pointer
c     to FlowSolution.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension idata(2)
c
c   Create ZoneIterativeData node
      call cg_ziter_write_f(iccg,ibase,igrid,'ZoneIterativeData',ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ZoneIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      idata(1)=32
      idata(2)=1
      call cg_array_write_f('FlowSolutionPointers',Character,2,idata,
     +    'FlowSolution',ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
 
      subroutine writebcs(iccg,ibase,igrid,idim,jdim,kdim,qj0,qk0,qi0,
     +  vj0,vk0,vi0,tj0,tk0,ti0,i2d,nummem)
c**********************************************************************
c     Purpose:  Store BC info necessary for CFL3D V6 in a
c     'UserDefinedData'.  If i2d=1, don't write qi0, vi0, or ti0 
c     (not ever needed, and writing them triples the size of the
c     CGNS file unnecessarily!).
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (integers)
c      qi0,qj0,qk0......arrays of BC data (PRIMITIVE quantities) (real)
c      vi0,vj0,vk0......arrays of eddy viscosity BC data (real)
c      ti0,tj0,tk0......arrays of turb quantity BC data (real)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      nummem...........one of dimensions for ti0,tj0,tk0 (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter(numnames=20)
c
      dimension qi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4), 
     .          qk0(jdim,idim-1,5,4)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
      dimension jdata(4)
      character*32 name(numnames)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Go to Zone node and create a user-defined node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_user_data_write_f('CFL3DBoundaryValues', ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nuser_data_f(nuserdata, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .gt. numnames) then
        write(901,'('' Error... numnames too small in writebcs'')')
        stop
      end if
      do n=1,nuserdata
        call cg_user_data_read_f(n, name(n), ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(n) .eq. 'CFL3DBoundaryValues') then
          nset=n
        end if
      enddo
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +  'UserDefinedData_t',nset,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Write BC qi0 data
      if (i2d .ne. 1) then
      jdata(1)=jdim
      jdata(2)=kdim
      jdata(3)=5
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Qi0_CFL3D',RealDouble,4,jdata,
     +    qi0,ier)
      else
        call cg_array_write_f('Qi0_CFL3D',RealSingle,4,jdata,
     +    qi0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      end if
c   Write BC qj0 data
      jdata(1)=kdim
      jdata(2)=idim-1
      jdata(3)=5
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Qj0_CFL3D',RealDouble,4,jdata,
     +    qj0,ier)
      else
        call cg_array_write_f('Qj0_CFL3D',RealSingle,4,jdata,
     +    qj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write BC qk0 data
      jdata(1)=jdim
      jdata(2)=idim-1
      jdata(3)=5
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Qk0_CFL3D',RealDouble,4,jdata,
     +    qk0,ier)
      else
        call cg_array_write_f('Qk0_CFL3D',RealSingle,4,jdata,
     +    qk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write BC vi0 data
      if (i2d .ne. 1) then
      jdata(1)=jdim
      jdata(2)=kdim
      jdata(3)=1
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Vi0_CFL3D',RealDouble,4,jdata,
     +    vi0,ier)
      else
        call cg_array_write_f('Vi0_CFL3D',RealSingle,4,jdata,
     +    vi0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      end if
c   Write BC vj0 data
      jdata(1)=kdim
      jdata(2)=idim-1
      jdata(3)=1
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Vj0_CFL3D',RealDouble,4,jdata,
     +    vj0,ier)
      else
        call cg_array_write_f('Vj0_CFL3D',RealSingle,4,jdata,
     +    vj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write BC vk0 data
      jdata(1)=jdim
      jdata(2)=idim-1
      jdata(3)=1
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Vk0_CFL3D',RealDouble,4,jdata,
     +    vk0,ier)
      else
        call cg_array_write_f('Vk0_CFL3D',RealSingle,4,jdata,
     +    vk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write BC ti0 data
      if (i2d .ne. 1) then
      jdata(1)=jdim
      jdata(2)=kdim
      jdata(3)=nummem
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Ti0_CFL3D',RealDouble,4,jdata,
     +    ti0,ier)
      else
        call cg_array_write_f('Ti0_CFL3D',RealSingle,4,jdata,
     +    ti0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      end if
c   Write BC tj0 data
      jdata(1)=kdim
      jdata(2)=idim-1
      jdata(3)=nummem
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Tj0_CFL3D',RealDouble,4,jdata,
     +    tj0,ier)
      else
        call cg_array_write_f('Tj0_CFL3D',RealSingle,4,jdata,
     +    tj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write BC tk0 data
      jdata(1)=jdim
      jdata(2)=idim-1
      jdata(3)=nummem
      jdata(4)=4
      if (idouble .eq. 1) then
        call cg_array_write_f('Tk0_CFL3D',RealDouble,4,jdata,
     +    tk0,ier)
      else
        call cg_array_write_f('Tk0_CFL3D',RealSingle,4,jdata,
     +    tk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine readbcs(iccg,ibase,igrid,idim,jdim,kdim,i2d,
     +  qj0,qk0,qi0,vj0,vk0,vi0,tj0,tk0,ti0,nummem)
c**********************************************************************
c     Purpose:  Reads CFL3D's BC info.
c     If i2d=1, don't need to read qi0, vi0, or ti0.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (integers)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      nummem...........one of dimensions for ti0,tj0,tk0 (integer)
c   OUTPUTS:
c      qi0,qj0,qk0......arrays of BC data (PRIMITIVE quantities
c                       (real)
c      vi0,vj0,vk0......arrays of eddy viscosity BC data (real)
c      ti0,tj0,tk0......arrays of turb quantity BC data (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=20)
c
      dimension qi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4), 
     .          qk0(jdim,idim-1,5,4)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
      dimension idata(2),jdata(4),idimvec(4)
      character name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   We must locate the BCs:
c   In new method, they are stored under Zone in a UserDefinedData node
c   called CFL3DBoundaryValues.  In old method, they were put under
c   ZoneIterativeData (because there was no better place to put them
c   at the time).  Check both places, just in case.
c
c   Go to Zone and look for UserDefinedData
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .gt. numnames) then
        write(901,'('' WARNING... numnames too small in readbcs'')')
        write(901,'(''   not changing default BCs... continuing'')')
        return
      end if
      if (nuserdata .ne. 0) then
        do n=1,nuserdata
          call cg_user_data_read_f(n,name(n),ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
        do n=1,nuserdata
          if(name(n) .eq. 'CFL3DBoundaryValues') goto 1001
        enddo
c   No CFL3DBoundaryValues... Try reading BCs from old location
        write(901,'('' no BC data exists in UserDefinedData nodes'')')
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     .    'ZoneIterativeData_t',1,'end')
        if (ier .ne. 0) then
          write(901,'('' WARNING... BCs not in ZoneIterativeData or'',
     +      '' a UserDefinedData array'')')
          write(901,'(''   not changing default BCs... continuing'')')
          return
        end if
        call cg_narrays_f(narrays,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (narrays .gt. numnames) then
          write(901,'('' WARNING.  Too many ZoneIterativeData'',
     +     '' arrays.'')')
          write(901,'(''   not changing default BCs... continuing'')')
          return
        end if
        write(901,'('' Reading BC info from ZoneIterativeData (old'',
     +    '' method)'')')
        write(901,'(''  ... it will be moved to new location of'',
     +    '' UserDefinedData node CFL3DBoundaryValues when written'')')
        do n=1,narrays
          call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
        goto 1002
c   Read BCs from new location
 1001   continue
        write(901,'('' Reading BC info from UserDefinedData'')')
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     +    'UserDefinedData_t',n,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_narrays_f(narrays,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (narrays .gt. numnames) then
          write(901,'('' WARNING.  Too many arrays in'',
     +      '' CFL3DBoundaryValues...'')')
          write(901,'(''   not changing default BCs... continuing'')')
          return
        end if
        do n=1,narrays
          call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
 1002   continue
      else
c   Try reading BCs from old location
        write(901,'('' no BC data exists in UserDefinedData nodes'')')
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     .    'ZoneIterativeData_t',1,'end')
        if (ier .ne. 0) then
          write(901,'('' WARNING... BCs not in ZoneIterativeData or'',
     +      '' a UserDefinedData array'')')
          write(901,'(''   not changing default BCs... continuing'')')
          return
        end if
        call cg_narrays_f(narrays,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if (narrays .gt. numnames) then
          write(901,'('' WARNING.  Too many ZoneIterativeData'',
     +      '' arrays.'')')
          write(901,'(''   not changing default BCs... continuing'')')
          return
        end if
        write(901,'('' Reading BC info from ZoneIterativeData (old'',
     +    '' method)'')')
        write(901,'(''  ... it will be moved to new location of'',
     +    '' UserDefinedData node CFL3DBoundaryValues when written'')')
        do n=1,narrays
          call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
      end if
c
c   Get BC qi0
      if (i2d .ne. 1) then
      do n=1,narrays
        if (name(n) .eq. 'Qi0_CFL3D') goto 101
      enddo
      write(901,'('' WARNING. No Qi0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 201
 101  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,qi0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,qi0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 201  continue
      end if
c
c   Get BC qj0
      do n=1,narrays
        if (name(n) .eq. 'Qj0_CFL3D') goto 102
      enddo
      write(901,'('' WARNING. No Qj0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 202
 102  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,qj0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,qj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 202  continue
c
c   Get BC qk0
      do n=1,narrays
        if (name(n) .eq. 'Qk0_CFL3D') goto 103
      enddo
      write(901,'('' WARNING. No Qk0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 203
 103  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,qk0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,qk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 203  continue
c
c   Get BC vi0
      if (i2d .ne. 1) then
      do n=1,narrays
        if (name(n) .eq. 'Vi0_CFL3D') goto 104
      enddo
      write(901,'('' WARNING. No Vi0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 204
 104  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,vi0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,vi0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 204  continue
      end if
c
c   Get BC vj0
      do n=1,narrays
        if (name(n) .eq. 'Vj0_CFL3D') goto 105
      enddo
      write(901,'('' WARNING. No Vj0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 205
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,vj0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,vj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 205  continue
c
c   Get BC vk0
      do n=1,narrays
        if (name(n) .eq. 'Vk0_CFL3D') goto 106
      enddo
      write(901,'('' WARNING. No Vk0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 206
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,vk0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,vk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 206  continue
c
c   Get BC ti0
      if (i2d .ne. 1) then
      do n=1,narrays
        if (name(n) .eq. 'Ti0_CFL3D') goto 107
      enddo
      write(901,'('' WARNING. No Ti0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 207
 107  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,ti0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,ti0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 207  continue
      end if
c
c   Get BC tj0
      do n=1,narrays
        if (name(n) .eq. 'Tj0_CFL3D') goto 108
      enddo
      write(901,'('' WARNING. No Tj0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 208
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,tj0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,tj0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 208  continue
c
c   Get BC tk0
      do n=1,narrays
        if (name(n) .eq. 'Tk0_CFL3D') goto 109
      enddo
      write(901,'('' WARNING. No Tk0_CFL3D node exists'')')
      write(901,'(''   keeping default and continuing'')')
      goto 209
 109  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,tk0,ier)
      else
        call cg_array_read_as_f(n,RealSingle,tk0,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
 209  continue
c
      return
      end

      subroutine wgrdmov(iccg,ibase,igrid,iuns,
     .                  itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt,ialph)
c**********************************************************************
c     Purpose: Creates 'RigidGridMotion' node and writes moving grid 
c     info under it (including CFL3D-specific info).
c     Also adds RigidGridMotionPointers to ZoneIterativeData.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      iuns.............flag to indicate whether unsteady or not 
c                       (integer)
c      itrans...........flag to indicate whether translating or not 
c                       (integer)
c      etc...
c      irotat...........flag to indicate whether rotating or not 
c                       (integer)
c      etc...
c      dt...............time step (real)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension data(3),data1(3,2),idata(2)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Add RigidGridMotionPointers to ZoneIterativeData node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ZoneIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      idata(1)=32
      idata(2)=1
      call cg_array_write_f('RigidGridMotionPointers',Character,2,idata,
     +    'RigidGridMotion',ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create SIDS-standard node for rigid motion
      if (itrans .eq. 0 .and. irotat .eq. 0) then
        call cg_rigid_motion_write_f(iccg,ibase,igrid,
     +    'RigidGridMotion',Null,ir,ier)
      else if (itrans .le. 1 .and. irotat .le. 1) then
        call cg_rigid_motion_write_f(iccg,ibase,igrid,
     +    'RigidGridMotion',ConstantRate,ir,ier)
      else
        call cg_rigid_motion_write_f(iccg,ibase,igrid,
     +    'RigidGridMotion',VariableRate,ir,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write Origin, Theta info
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      idata(1)=3
      idata(2)=2
      if (ialph .eq. 0) then
        data1(1,1)=xorig0
        data1(2,1)=yorig0
        data1(3,1)=zorig0
        data1(1,2)=xorig
        data1(2,2)=yorig
        data1(3,2)=zorig
        data(1)=thetax
        data(2)=thetay
        data(3)=thetaz
      else
        data1(1,1)=xorig0
        data1(2,1)=zorig0
        data1(3,1)=-yorig0
        data1(1,2)=xorig
        data1(2,2)=zorig
        data1(3,2)=-yorig
        data(1)=thetax
        data(2)=thetaz
        data(3)=-thetay
      end if
      if (idouble .eq. 1) then
        call cg_array_write_f('OriginLocation',RealDouble,2,idata,
     +    data1,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('RigidRotationAngle',RealDouble,1,3,
     +    data,ier)
        if (ier .ne. 0) call cg_error_exit_f
      else
        call cg_array_write_f('OriginLocation',RealSingle,2,idata,
     +    data1,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_array_write_f('RigidRotationAngle',RealSingle,1,3,
     +    data,ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
c
c   Next write other (cfl3d-specific) data:
      call cg_user_data_write_f('CFL3DParameters', ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'UserDefinedData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfreqt',RealDouble,1,1,rfreqt,ier)
      else
        call cg_array_write_f('Rfreqt',RealSingle,1,1,rfreqt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=utrans
      data(2)=vtrans
      data(3)=wtrans
      if (idouble .eq. 1) then
        call cg_array_write_f('Veltrans',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Veltrans',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=dxmx
      data(2)=dymx
      data(3)=dzmx
      if (idouble .eq. 1) then
        call cg_array_write_f('Dmx',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Dmx',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfreqtmc',RealDouble,1,1,rfreqtmc,ier)
      else
        call cg_array_write_f('Rfreqtmc',RealSingle,1,1,rfreqtmc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=xorigmc
      data(2)=yorigmc
      data(3)=zorigmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Origmc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Origmc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=xorig0mc
      data(2)=yorig0mc
      data(3)=zorig0mc
      if (idouble .eq. 1) then
        call cg_array_write_f('Orig0mc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Orig0mc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=utransmc
      data(2)=vtransmc
      data(3)=wtransmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Veltransmc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Veltransmc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=xmc
      data(2)=ymc
      data(3)=zmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Mc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Mc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=dxmxmc
      data(2)=dymxmc
      data(3)=dzmxmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Dmxmc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Dmxmc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfreqr',RealDouble,1,1,rfreqr,ier)
      else
        call cg_array_write_f('Rfreqr',RealSingle,1,1,rfreqr,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=omegax
      data(2)=omegay
      data(3)=omegaz
      if (idouble .eq. 1) then
        call cg_array_write_f('Omega',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Omega',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=dthxmx
      data(2)=dthymx
      data(3)=dthzmx
      if (idouble .eq. 1) then
        call cg_array_write_f('Dthmx',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Dthmx',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfreqrmc',RealDouble,1,1,rfreqrmc,ier)
      else
        call cg_array_write_f('Rfreqrmc',RealSingle,1,1,rfreqrmc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=thetaxmc
      data(2)=thetaymc
      data(3)=thetazmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Thetamc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Thetamc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=omegaxmc
      data(2)=omegaymc
      data(3)=omegazmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Omegamc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Omegamc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      data(1)=dthxmxmc
      data(2)=dthymxmc
      data(3)=dthzmxmc
      if (idouble .eq. 1) then
        call cg_array_write_f('Dthmxmc',RealDouble,1,3,data,ier)
      else
        call cg_array_write_f('Dthmxmc',RealSingle,1,3,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Time2',RealDouble,1,1,time2,ier)
      else
        call cg_array_write_f('Time2',RealSingle,1,1,time2,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Time2mc',RealDouble,1,1,time2mc,ier)
      else
        call cg_array_write_f('Time2mc',RealSingle,1,1,time2mc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Dtvar',RealDouble,1,1,dt,ier)
      else
        call cg_array_write_f('Dtvar',RealSingle,1,1,dt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      call cg_array_write_f('Itrans',Integer,1,1,itrans,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Itransmc',Integer,1,1,itransmc,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Irotat',Integer,1,1,irotat,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Irotatmc',Integer,1,1,irotatmc,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Iuns',Integer,1,1,iuns,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine rgrdmov(iccg,ibase,igrid,ialph,
     .                  itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt)
c**********************************************************************
c     Purpose:  Reads 'GridMotionDefinition' node moving grid info.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c   OUTPUTS:
c      itrans...........flag to indicate whether translating or not 
c                       (integer)
c      etc...
c      irotat...........flag to indicate whether rotating or not 
c                       (integer)
c      etc...
c      dt...............time step (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=50)
c
      dimension data(3),data1(3,2),idimvec(4)
      character*32 name(numnames)
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Initial values of parameters
      itrans=0
      rfreqt=0.
      xorig=0.
      yorig=0.
      zorig=0.
      xorig0=0.
      yorig0=0.
      zorig0=0.
      utrans=0.
      vtrans=0.
      wtrans=0.
      dxmx=0.
      dymx=0.
      dzmx=0.
      itransmc=0
      rfreqtmc=0.
      xorigmc=0.
      yorigmc=0.
      zorigmc=0.
      xorig0mc=0.
      yorig0mc=0.
      zorig0mc=0.
      utransmc=0.
      vtransmc=0.
      wtransmc=0.
      xmc=0.
      ymc=0.
      zmc=0.
      dxmxmc=0.
      dymxmc=0.
      dzmxmc=0.
      irotat=0
      rfreqr=0.
      thetax=0.
      thetay=0.
      thetaz=0.
      omegax=0.
      omegay=0.
      omegaz=0.
      dthxmx=0.
      dthymx=0.
      dthzmx=0.
      irotatmc=0
      rfreqrmc=0.
      thetaxmc=0.
      thetaymc=0.
      thetazmc=0.
      omegaxmc=0.
      omegaymc=0.
      omegazmc=0.
      dthxmxmc=0.
      dthymxmc=0.
      dthzmxmc=0.
      time2=0.
      time2mc=0.
      dt=0.
c   Find out how many RigidGridMotion_t nodes there are
      call cg_n_rigid_motions_f(iccg,ibase,igrid,nmotions,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nmotions .eq. 0) return
      if (nmotions .gt. 1) then
        write(901,'('' Error, more than 1 RigidGridMotion_t nodes!'')')
        stop
      end if
c   Go to RigidGridMotion_t node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many RigidGridMotion arrays.'')')
        write(901,'(''   increase numnames in rgrdmov'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Get Origin
      do n=1,narrays
        if (name(n) .eq. 'OriginLocation') goto 101
      enddo
      write(901,'('' Error. No OriginLocation node exists'')')
      stop
 101  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data1,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data1,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
        xorig0=data1(1,1)
        yorig0=data1(2,1)
        zorig0=data1(3,1)
        xorig=data1(1,2)
        yorig=data1(2,2)
        zorig=data1(3,2)
      else
        xorig0=data1(1,1)
        zorig0=data1(2,1)
        yorig0=-data1(3,1)
        xorig=data1(1,2)
        zorig=data1(2,2)
        yorig=-data1(3,2)
      end if
c
c   Get Thetas
      do n=1,narrays
        if (name(n) .eq. 'RigidRotationAngle') goto 102
      enddo
      write(901,'('' Error. No RigidRotationAngle node exists'')')
      stop
 102  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
        thetax=data(1)
        thetay=data(2)
        thetaz=data(3)
      else
        thetax=data(1)
        thetaz=data(2)
        thetay=-data(3)
      end if
c
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .ne. 1) then
        write(901,'('' Error... should be 1 user data node in'', 
     +    '' RigidGridMotion'')')
        stop
      else
        call cg_user_data_read_f(1,name(1),ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(1) .eq. 'CFL3DParameters') goto 1004
        write(901,'('' Error. No CFL3DParameters node exists'')')
        stop
 1004   continue
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'UserDefinedData_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
      end if
c
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many RigidGridMotion arrays.'')')
        write(901,'(''   increase numnames in rgrdmov'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get Rfreqt
      do n=1,narrays
        if (name(n) .eq. 'Rfreqt') goto 103
      enddo
      write(901,'('' Error. No Rfreqt node exists'')')
      stop
 103  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfreqt,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfreqt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Veltrans
      do n=1,narrays
        if (name(n) .eq. 'Veltrans') goto 104
      enddo
      write(901,'('' Error. No Veltrans node exists'')')
      stop
 104  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      utrans=data(1)
      vtrans=data(2)
      wtrans=data(3)
c
c   Get Dmx
      do n=1,narrays
        if (name(n) .eq. 'Dmx') goto 105
      enddo
      write(901,'('' Error. No Dmx node exists'')')
      stop
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      dxmx=data(1)
      dymx=data(2)
      dzmx=data(3)
c
c   Get Rfreqtmc
      do n=1,narrays
        if (name(n) .eq. 'Rfreqtmc') goto 106
      enddo
      write(901,'('' Error. No Rfreqtmc node exists'')')
      stop
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfreqtmc,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfreqtmc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Origmc
      do n=1,narrays
        if (name(n) .eq. 'Origmc') goto 107
      enddo
      write(901,'('' Error. No Origmc node exists'')')
      stop
 107  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      xorigmc=data(1)
      yorigmc=data(2)
      zorigmc=data(3)
c
c   Get Orig0mc
      do n=1,narrays
        if (name(n) .eq. 'Orig0mc') goto 108
      enddo
      write(901,'('' Error. No Orig0mc node exists'')')
      stop
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      xorig0mc=data(1)
      yorig0mc=data(2)
      zorig0mc=data(3)
c
c   Get Veltransmc
      do n=1,narrays
        if (name(n) .eq. 'Veltransmc') goto 109
      enddo
      write(901,'('' Error. No Veltransmc node exists'')')
      stop
 109  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      utransmc=data(1)
      vtransmc=data(2)
      wtransmc=data(3)
c
c   Get Mc
      do n=1,narrays
        if (name(n) .eq. 'Mc') goto 110
      enddo
      write(901,'('' Error. No Mc node exists'')')
      stop
 110  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      xmc=data(1)
      ymc=data(2)
      zmc=data(3)
c
c   Get Dmxmc
      do n=1,narrays
        if (name(n) .eq. 'Dmxmc') goto 111
      enddo
      write(901,'('' Error. No Dmxmc node exists'')')
      stop
 111  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      dxmxmc=data(1)
      dymxmc=data(2)
      dzmxmc=data(3)
c
c   Get Rfreqr
      do n=1,narrays
        if (name(n) .eq. 'Rfreqr') goto 112
      enddo
      write(901,'('' Error. No Rfreqr node exists'')')
      stop
 112  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfreqr,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfreqr,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Omega
      do n=1,narrays
        if (name(n) .eq. 'Omega') goto 113
      enddo
      write(901,'('' Error. No Omega node exists'')')
      stop
 113  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      omegax=data(1)
      omegay=data(2)
      omegaz=data(3)
c
c   Get Dthmx
      do n=1,narrays
        if (name(n) .eq. 'Dthmx') goto 114
      enddo
      write(901,'('' Error. No Dthmx node exists'')')
      stop
 114  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      dthxmx=data(1)
      dthymx=data(2)
      dthzmx=data(3)
c
c   Get Rfreqrmc
      do n=1,narrays
        if (name(n) .eq. 'Rfreqrmc') goto 115
      enddo
      write(901,'('' Error. No Rfreqrmc node exists'')')
      stop
 115  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfreqrmc,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfreqrmc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Thetamc
      do n=1,narrays
        if (name(n) .eq. 'Thetamc') goto 116
      enddo
      write(901,'('' Error. No Thetamc node exists'')')
      stop
 116  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      thetaxmc=data(1)
      thetaymc=data(2)
      thetazmc=data(3)
c
c   Get Omegamc
      do n=1,narrays
        if (name(n) .eq. 'Omegamc') goto 117
      enddo
      write(901,'('' Error. No Omegamc node exists'')')
      stop
 117  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      omegaxmc=data(1)
      omegaymc=data(2)
      omegazmc=data(3)
c
c   Get Dthmxmc
      do n=1,narrays
        if (name(n) .eq. 'Dthmxmc') goto 118
      enddo
      write(901,'('' Error. No Dthmxmc node exists'')')
      stop
 118  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,data,ier)
      else
        call cg_array_read_as_f(n,RealSingle,data,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      dthxmxmc=data(1)
      dthymxmc=data(2)
      dthzmxmc=data(3)
c
c   Get Time2
      do n=1,narrays
        if (name(n) .eq. 'Time2') goto 119
      enddo
      write(901,'('' Error. No Time2 node exists'')')
      stop
 119  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,time2,ier)
      else
        call cg_array_read_as_f(n,RealSingle,time2,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Time2mc
      do n=1,narrays
        if (name(n) .eq. 'Time2mc') goto 120
      enddo
      write(901,'('' Error. No Time2mc node exists'')')
      stop
 120  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,time2mc,ier)
      else
        call cg_array_read_as_f(n,RealSingle,time2mc,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Dtvar
      do n=1,narrays
        if (name(n) .eq. 'Dtvar') goto 121
      enddo
      write(901,'('' Error. No Dtvar node exists'')')
      stop
 121  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,dt,ier)
      else
        call cg_array_read_as_f(n,RealSingle,dt,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Itrans
      do n=1,narrays
        if (name(n) .eq. 'Itrans') goto 122
      enddo
      write(901,'('' Error. No Itrans node exists'')')
      stop
 122  continue
      call cg_array_read_as_f(n,Integer,itrans,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Itransmc
      do n=1,narrays
        if (name(n) .eq. 'Itransmc') goto 123
      enddo
      write(901,'('' Error. No Itransmc node exists'')')
      stop
 123  continue
      call cg_array_read_as_f(n,Integer,itransmc,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Irotat
      do n=1,narrays
        if (name(n) .eq. 'Irotat') goto 124
      enddo
      write(901,'('' Error. No Irotat node exists'')')
      stop
 124  continue
      call cg_array_read_as_f(n,Integer,irotat,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Irotatmc
      do n=1,narrays
        if (name(n) .eq. 'Irotatmc') goto 125
      enddo
      write(901,'('' Error. No Irotatmc node exists'')')
      stop
 125  continue
      call cg_array_read_as_f(n,Integer,irotatmc,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine readiuns(iccg,ibase,igrid,iuns)
c**********************************************************************
c     Purpose: Reads 'GridMotionDefinition' node and gets iuns 
c     parameter.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c   OUTPUTS:
c      iuns.............flag that tells what kind of unsteady 
c                       motion it is (integer)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=50)
c
      dimension idimvec(4)
      character*32 name(numnames)
c
c   Find out how many RigidGridMotion_t nodes there are
      call cg_n_rigid_motions_f(iccg,ibase,igrid,nmotions,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nmotions .eq. 0) then
        write(901,'('' No RigidGridMotion nodes, igrid='',i5)') igrid
        write(901,'(''   setting iuns=0 and continuing'')')
        iuns=0
        return
      end if
      if (nmotions .gt. 1) then
        write(901,'('' Error, more than 1 RigidGridMotion_t nodes!'')')
        stop
      end if
c   Go to RigidGridMotion_t node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Find out how many UserDefinedData_t nodes there are
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .ne. 1) then
        write(901,'('' Error... should be 1 user data node in'',
     +    '' RigidGridMotion'')')
        stop
      else
        call cg_user_data_read_f(1,name(1),ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(1) .eq. 'CFL3DParameters') goto 1004
        write(901,'('' Error. No CFL3DParameters node exists'')')
        stop
 1004   continue
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'RigidGridMotion_t',1,'UserDefinedData_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
      end if
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many RigidGridMotion arrays.'')')
        write(901,'(''   increase numnames in readiuns'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Get Iuns
      do n=1,narrays
        if (name(n) .eq. 'Iuns') goto 101
      enddo
      write(901,'('' Error. No Iuns node exists'')')
      stop
 101  continue
      call cg_array_read_as_f(n,Integer,iuns,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine wdeform(iccg,ibase,igrid,idefrm,
     .            nsegdfrm,idfrmseg,utrnsae,
     .            vtrnsae,wtrnsae,omgxae,
     .            omgyae,omgzae,xorgae,
     .            yorgae,zorgae,thtxae,
     .            thtyae,thtzae,rfrqtae,
     .            rfrqrae,icsi,icsf,jcsi,
     .            jcsf,kcsi,kcsf,jdima,kdima,idima,
     .            jdim,kdim,idim,x,y,z,xnm2,ynm2,znm2,wk,ialph,i2d)
c**********************************************************************
c     Purpose: Creates 'ArbitraryGridMotion' node and writes 
c     arbitrary moving grid info (including CFL3D-specific info).
c     Also adds ArbitraryGridMotionPointers and GridCoordinatesPointers
c     to ZoneIterativeData.
c     Writes grids in the following order:
c
c     (idim,jdim,kdim) <- necessary for CGNS!,
c
c     even though they are stored in CFL3D in the order: jdim,kdim,idim
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idefrm...........CFL3D-specific parameter (integer)
c      nsegdfrm.........CFL3D-specific parameter (integer)
c      idfrmseg.........CFL3D-specific parameter (integer)
c      utrnsae..........CFL3D-specific parameter (real)
c      vtrnsae..........CFL3D-specific parameter (real)
c      wtrnsae..........CFL3D-specific parameter (real)
c      omgxae...........CFL3D-specific parameter (real)
c      omgyae...........CFL3D-specific parameter (real)
c      omgzae...........CFL3D-specific parameter (real)
c      xorgae...........CFL3D-specific parameter (real)
c      yorgae...........CFL3D-specific parameter (real)
c      zorgae...........CFL3D-specific parameter (real)
c      thtxae...........CFL3D-specific parameter (real)
c      thtyae...........CFL3D-specific parameter (real)
c      thtzae...........CFL3D-specific parameter (real)
c      rfrqtae..........CFL3D-specific parameter (real)
c      rfrqrae..........CFL3D-specific parameter (real)
c      icsi.............CFL3D-specific parameter (integer)
c      icsf.............CFL3D-specific parameter (integer)
c      jcsi.............CFL3D-specific parameter (integer)
c      jcsf.............CFL3D-specific parameter (integer)
c      kcsi.............CFL3D-specific parameter (integer)
c      kcsf.............CFL3D-specific parameter (integer)
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid"), at finest level (integers)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) (integers)
c      x................new grid x-value (real)
c      y................new grid y-value (real)
c      z................new grid z-value (real)
c      xnm2.............last time-step grid x-value (real)
c      ynm2.............last time-step grid y-value (real)
c      znm2.............last time-step grid z-value (real)
c      wk...............working space needed, of dimension
c                       (idima*jdima*kdima) or larger (real)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      dimension idfrmseg(nsegdfrm),utrnsae(nsegdfrm),
     .          vtrnsae(nsegdfrm),wtrnsae(nsegdfrm),omgxae(nsegdfrm),
     .          omgyae(nsegdfrm),omgzae(nsegdfrm),xorgae(nsegdfrm),
     .          yorgae(nsegdfrm),zorgae(nsegdfrm),thtxae(nsegdfrm),
     .          thtyae(nsegdfrm),thtzae(nsegdfrm),rfrqtae(nsegdfrm),
     .          rfrqrae(nsegdfrm),icsi(nsegdfrm),icsf(nsegdfrm),
     .          jcsi(nsegdfrm),jcsf(nsegdfrm),kcsi(nsegdfrm),
     .          kcsf(nsegdfrm)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     .          xnm2(jdim,kdim,idim),ynm2(jdim,kdim,idim),
     .          znm2(jdim,kdim,idim)
      dimension wk(idima,jdima,kdima)
      dimension idata(2),jdata(3)
      character*52 string1
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Add ArbitraryGridMotionPointers and GridCoordinatesPointers
c   to ZoneIterativeData node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ZoneIterativeData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      idata(1)=32
      idata(2)=1
      call cg_array_write_f('ArbitraryGridMotionPointers',Character,2,
     +    idata,'ArbitraryGridMotion',ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('GridCoordinatesPointers',Character,2,
     +    idata,'MovedGrid',ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Create SIDS-standard node for arbitrary motion
      call cg_arbitrary_motion_write_f(iccg,ibase,igrid,
     +  'ArbitraryGridMotion',DeformingGrid,ir,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Set dimensions
      jdata(1)=idima
      jdata(2)=jdima
      jdata(3)=kdima
c   Initialize wk array to zero
      do i=1,idima
      do j=1,jdima
      do k=1,kdima
        wk(i,j,k)=0.
      enddo
      enddo
      enddo
c   Write CFL3D-specific info
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'end')
      call cg_user_data_write_f('CFL3DParameters', ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'UserDefinedData_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
      call cg_array_write_f('Idefrm',Integer,1,1,idefrm,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Nsegdfrm',Integer,1,1,nsegdfrm,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Idfrmseg',Integer,1,nsegdfrm,idfrmseg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Utrnsae',RealDouble,1,nsegdfrm,
     +   utrnsae,ier)
      else
        call cg_array_write_f('Utrnsae',RealSingle,1,nsegdfrm,
     +   utrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Vtrnsae',RealDouble,1,nsegdfrm,
     +   vtrnsae,ier)
      else
        call cg_array_write_f('Vtrnsae',RealSingle,1,nsegdfrm,
     +   vtrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Wtrnsae',RealDouble,1,nsegdfrm,
     +   wtrnsae,ier)
      else
        call cg_array_write_f('Wtrnsae',RealSingle,1,nsegdfrm,
     +   wtrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Omgxae',RealDouble,1,nsegdfrm,
     +   omgxae,ier)
      else
        call cg_array_write_f('Omgxae',RealSingle,1,nsegdfrm,
     +   omgxae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Omgyae',RealDouble,1,nsegdfrm,
     +   omgyae,ier)
      else
        call cg_array_write_f('Omgyae',RealSingle,1,nsegdfrm,
     +   omgyae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Omgzae',RealDouble,1,nsegdfrm,
     +   omgzae,ier)
      else
        call cg_array_write_f('Omgzae',RealSingle,1,nsegdfrm,
     +   omgzae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Xorgae',RealDouble,1,nsegdfrm,
     +   xorgae,ier)
      else
        call cg_array_write_f('Xorgae',RealSingle,1,nsegdfrm,
     +   xorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Yorgae',RealDouble,1,nsegdfrm,
     +   yorgae,ier)
      else
        call cg_array_write_f('Yorgae',RealSingle,1,nsegdfrm,
     +   yorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Zorgae',RealDouble,1,nsegdfrm,
     +   zorgae,ier)
      else
        call cg_array_write_f('Zorgae',RealSingle,1,nsegdfrm,
     +   zorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Thtxae',RealDouble,1,nsegdfrm,
     +   thtxae,ier)
      else
        call cg_array_write_f('Thtxae',RealSingle,1,nsegdfrm,
     +   thtxae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Thtyae',RealDouble,1,nsegdfrm,
     +   thtyae,ier)
      else
        call cg_array_write_f('Thtyae',RealSingle,1,nsegdfrm,
     +   thtyae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Thtzae',RealDouble,1,nsegdfrm,
     +   thtzae,ier)
      else
        call cg_array_write_f('Thtzae',RealSingle,1,nsegdfrm,
     +   thtzae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfrqtae',RealDouble,1,nsegdfrm,
     +   rfrqtae,ier)
      else
        call cg_array_write_f('Rfrqtae',RealSingle,1,nsegdfrm,
     +   rfrqtae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
        call cg_array_write_f('Rfrqrae',RealDouble,1,nsegdfrm,
     +   rfrqrae,ier)
      else
        call cg_array_write_f('Rfrqrae',RealSingle,1,nsegdfrm,
     +   rfrqrae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Icsi',Integer,1,nsegdfrm,icsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Icsf',Integer,1,nsegdfrm,icsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Jcsi',Integer,1,nsegdfrm,jcsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Jcsf',Integer,1,nsegdfrm,jcsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Kcsi',Integer,1,nsegdfrm,kcsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_array_write_f('Kcsf',Integer,1,nsegdfrm,kcsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Write moved grid
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_grid_write_f(iccg,ibase,igrid,'MovedGrid',ig,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'GridCoordinates_t',ig,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   First, get "factor" in case want to write every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
        string1='Fine grid level'
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
        string1='Every-2nd-gridpoint written (1 level down)'
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
        string1='Every-4th-gridpoint written (2 levels down)'
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
        string1='Every-8th-gridpoint written (3 levels down)'
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
        string1='Every-16th-gridpoint written (4 levels down)'
      else
        write(901,'('' Error.  Desired grid level not supported'')')
        stop
      end if
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
      call cg_descriptor_write_f('Information',string1,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Write X
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=x(jj,kk,ii)
      enddo
      enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateX',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateX',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write Y
      if (ialph .eq. 0) then
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=y(jj,kk,ii)
      enddo
      enddo
      enddo
      else
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=z(jj,kk,ii)
      enddo
      enddo
      enddo
      end if
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateY',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateY',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write Z
      if (ialph .eq. 0) then
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=z(jj,kk,ii)
      enddo
      enddo
      enddo
      else
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=-y(jj,kk,ii)
      enddo
      enddo
      enddo
      end if
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateZ',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateZ',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Write moved grid at earlier time
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_grid_write_f(iccg,ibase,igrid,'MovedGridLastDT',ig,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'GridCoordinates_t',ig,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_descriptor_write_f('Information',string1,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   Write X
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=xnm2(jj,kk,ii)
      enddo
      enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateX',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateX',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write Y
      if (ialph .eq. 0) then
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=ynm2(jj,kk,ii)
      enddo
      enddo
      enddo
      else
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=znm2(jj,kk,ii)
      enddo
      enddo
      enddo
      end if
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateY',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateY',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c   Write Z
      if (ialph .eq. 0) then
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=znm2(jj,kk,ii)
      enddo
      enddo
      enddo
      else
      ii=0
      do i=1,idima,nfaci
      ii=ii+1
      jj=0
      do j=1,jdima,nfac
      jj=jj+1
      kk=0
      do k=1,kdima,nfac
        kk=kk+1
        wk(i,j,k)=-ynm2(jj,kk,ii)
      enddo
      enddo
      enddo
      end if
      if (idouble .eq. 1) then
        call cg_array_write_f('CoordinateZ',RealDouble,3,jdata,
     +   wk,ier)
      else
        call cg_array_write_f('CoordinateZ',RealSingle,3,jdata,
     +   wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine getnsegdfrm(iccg,ibase,igrid,nsegdfrm)
c**********************************************************************
c     Purpose:  Get parameter nsegdfrm from deforming grid info.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c   OUTPUTS:
c      nsegdfrm.........CFL3D-specific parameter (integer)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=50)
c
      dimension idimvec(4)
      character*32 name(numnames)
c
c   Find out how many ArbitraryGridMotion_t nodes there are
      call cg_n_arbitrary_motions_f(iccg,ibase,igrid,nmotions,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nmotions .eq. 0) then
        write(901,'('' Error.  No ArbitraryGridMotion nodes,'',
     +   '' igrid='',i5)') igrid
        write(901,'(''   cannot get nsegdfrm'')')
        stop
      end if
      if (nmotions .gt. 1) then
        write(901,'('' Error, more than 1 ArbitraryGridMotion_t'',
     +   '' nodes!'')')
        stop
      end if
c   Go to ArbitraryGridMotion_t node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Find out how many UserDefinedData_t nodes there are
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .ne. 1) then
        write(901,'('' Error... should be 1 user data node in'',
     +    '' ArbitraryGridMotion'')')
        stop
      else
        call cg_user_data_read_f(1,name(1),ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(1) .eq. 'CFL3DParameters') goto 1004
        write(901,'('' Error. No CFL3DParameters node exists'')')
        stop
 1004   continue
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'UserDefinedData_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
      end if
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many ArbitraryGridMotion arrays.'')')
        write(901,'(''   increase numnames in readiuns'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get Nsegdfrm
      do n=1,narrays
        if (name(n) .eq. 'Nsegdfrm') goto 101
      enddo
      write(901,'('' Error. No Nsegdfrm node exists'')')
      stop
 101  continue
      call cg_array_read_as_f(n,Integer,nsegdfrm,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine rdeform(iccg,ibase,igrid,nsegdfrm,
     .            wk,jdima,kdima,idima,i2d,
     .            jdim,kdim,idim,ialph,idefrm,utrnsae,
     .            vtrnsae,wtrnsae,omgxae,
     .            omgyae,omgzae,xorgae,
     .            yorgae,zorgae,thtxae,
     .            thtyae,thtzae,rfrqtae,
     .            rfrqrae,icsi,icsf,jcsi,
     .            jcsf,kcsi,kcsf,
     .            x,y,z,xnm2,ynm2,znm2)
c**********************************************************************
c     Purpose:  Reads deforming grid info (note:  don't need
c     to get idfrmseg).
c     Gets results in the following order:
c
c          (jdim,kdim,idim) <- necessary for CFL3D!,
c
c     and assumes they are stored in the CGNS database in the order:
c     (idim,jdim,kdim)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside 
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside 
c                       this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      nsegdfrm.........CFL3D-specific parameter that defines the
c                       array sizes (integer)
c      wk...............working space needed, of dimension
c                       (idim*jdim*kdim) or larger (real)
c      idima,jdima,kdima...expected dimensions of this zone
c                          (zone "igrid"), at finest level (integers)
c      i2d..............0 if 3-D, 1 if 2-D (integer)
c      idim,jdim,kdim...actual (GRIDPOINT) dimensions of existing data
c                       for this zone (e.g., could be every other point
c                       of fine grid) (integers)
c      ialph............parameter in CFL3D for determining whether y
c                       or z is "up" (integer)
c   OUTPUTS:
c      idefrm...........CFL3D-specific parameter (integer)
c      utrnsae..........CFL3D-specific parameter (real)
c      vtrnsae..........CFL3D-specific parameter (real)
c      wtrnsae..........CFL3D-specific parameter (real)
c      omgxae...........CFL3D-specific parameter (real)
c      omgyae...........CFL3D-specific parameter (real)
c      omgzae...........CFL3D-specific parameter (real)
c      xorgae...........CFL3D-specific parameter (real)
c      yorgae...........CFL3D-specific parameter (real)
c      zorgae...........CFL3D-specific parameter (real)
c      thtxae...........CFL3D-specific parameter (real)
c      thtyae...........CFL3D-specific parameter (real)
c      thtzae...........CFL3D-specific parameter (real)
c      rfrqtae..........CFL3D-specific parameter (real)
c      rfrqrae..........CFL3D-specific parameter (real)
c      icsi.............CFL3D-specific parameter (integer)
c      icsf.............CFL3D-specific parameter (integer)
c      jcsi.............CFL3D-specific parameter (integer)
c      jcsf.............CFL3D-specific parameter (integer)
c      kcsi.............CFL3D-specific parameter (integer)
c      kcsf.............CFL3D-specific parameter (integer)
c      x................new grid x-value (real)
c      y................new grid y-value (real)
c      z................new grid z-value (real)
c      xnm2.............last time-step grid x-value (real)
c      ynm2.............last time-step grid y-value (real)
c      znm2.............last time-step grid z-value (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=50)
c
      dimension utrnsae(nsegdfrm),
     .          vtrnsae(nsegdfrm),wtrnsae(nsegdfrm),omgxae(nsegdfrm),
     .          omgyae(nsegdfrm),omgzae(nsegdfrm),xorgae(nsegdfrm),
     .          yorgae(nsegdfrm),zorgae(nsegdfrm),thtxae(nsegdfrm),
     .          thtyae(nsegdfrm),thtzae(nsegdfrm),rfrqtae(nsegdfrm),
     .          rfrqrae(nsegdfrm),icsi(nsegdfrm),icsf(nsegdfrm),
     .          jcsi(nsegdfrm),jcsf(nsegdfrm),kcsi(nsegdfrm),
     .          kcsf(nsegdfrm)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     .          xnm2(jdim,kdim,idim),ynm2(jdim,kdim,idim),
     .          znm2(jdim,kdim,idim)
      dimension wk(idima,jdima,kdima)
      dimension idimvec(4)
      character*32 name(numnames),nameg(3),namei
      character*52 text
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Find out how many ArbitraryGridMotion_t nodes there are
      call cg_n_arbitrary_motions_f(iccg,ibase,igrid,nmotions,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nmotions .eq. 0) then
        write(901,'('' Error.  No ArbitraryGridMotion nodes,'',
     +   '' igrid='',i5)') igrid
        write(901,'(''   cannot get nsegdfrm'')')
        stop
      end if
      if (nmotions .gt. 1) then
        write(901,'('' Error, more than 1 ArbitraryGridMotion_t'',
     +   '' nodes!'')')
        stop
      end if
c   Go to ArbitraryGridMotion_t node
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Find out how many UserDefinedData_t nodes there are
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .ne. 1) then
        write(901,'('' Error... should be 1 user data node in'',
     +    '' ArbitraryGridMotion'')')
        stop
      else
        call cg_user_data_read_f(1,name(1),ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(1) .eq. 'CFL3DParameters') goto 1004
        write(901,'('' Error. No CFL3DParameters node exists'')')
        stop
 1004   continue
        call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'ArbitraryGridMotion_t',1,'UserDefinedData_t',1,'end')
        if (ier .ne. 0) call cg_error_exit_f
      end if
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many ArbitraryGridMotion arrays.'')')
        write(901,'(''   increase numnames in rdeform'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get Idefrm
      do n=1,narrays
        if (name(n) .eq. 'Idefrm') goto 101
      enddo
      write(901,'('' Error. No Idefrm node exists'')')
      stop
 101  continue
      call cg_array_read_as_f(n,Integer,idefrm,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Utrnsae
      do n=1,narrays
        if (name(n) .eq. 'Utrnsae') goto 102
      enddo
      write(901,'('' Error. No Utrnsae node exists'')')
      stop
 102  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,utrnsae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,utrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Vtrnsae
      do n=1,narrays
        if (name(n) .eq. 'Vtrnsae') goto 103
      enddo
      write(901,'('' Error. No Vtrnsae node exists'')')
      stop
 103  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,vtrnsae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,vtrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Wtrnsae
      do n=1,narrays
        if (name(n) .eq. 'Wtrnsae') goto 104
      enddo
      write(901,'('' Error. No Wtrnsae node exists'')')
      stop
 104  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wtrnsae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wtrnsae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Omgxae
      do n=1,narrays
        if (name(n) .eq. 'Omgxae') goto 105
      enddo
      write(901,'('' Error. No Omgxae node exists'')')
      stop
 105  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,omgxae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,omgxae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Omgyae
      do n=1,narrays
        if (name(n) .eq. 'Omgyae') goto 106
      enddo
      write(901,'('' Error. No Omgyae node exists'')')
      stop
 106  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,omgyae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,omgyae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Omgzae
      do n=1,narrays
        if (name(n) .eq. 'Omgzae') goto 107
      enddo
      write(901,'('' Error. No Omgzae node exists'')')
      stop
 107  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,omgzae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,omgzae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Xorgae
      do n=1,narrays
        if (name(n) .eq. 'Xorgae') goto 108
      enddo
      write(901,'('' Error. No Xorgae node exists'')')
      stop
 108  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,xorgae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,xorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Yorgae
      do n=1,narrays
        if (name(n) .eq. 'Yorgae') goto 109
      enddo
      write(901,'('' Error. No Yorgae node exists'')')
      stop
 109  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,yorgae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,yorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Zorgae
      do n=1,narrays
        if (name(n) .eq. 'Zorgae') goto 110
      enddo
      write(901,'('' Error. No Zorgae node exists'')')
      stop
 110  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,zorgae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,zorgae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Thtxae
      do n=1,narrays
        if (name(n) .eq. 'Thtxae') goto 111
      enddo
      write(901,'('' Error. No Thtxae node exists'')')
      stop
 111  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,thtxae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,thtxae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Thtyae
      do n=1,narrays
        if (name(n) .eq. 'Thtyae') goto 112
      enddo
      write(901,'('' Error. No Thtyae node exists'')')
      stop
 112  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,thtyae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,thtyae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Thtzae
      do n=1,narrays
        if (name(n) .eq. 'Thtzae') goto 113
      enddo
      write(901,'('' Error. No Thtzae node exists'')')
      stop
 113  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,thtzae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,thtzae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Rfrqtae
      do n=1,narrays
        if (name(n) .eq. 'Rfrqtae') goto 114
      enddo
      write(901,'('' Error. No Rfrqtae node exists'')')
      stop
 114  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfrqtae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfrqtae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Rfrqrae
      do n=1,narrays
        if (name(n) .eq. 'Rfrqrae') goto 115
      enddo
      write(901,'('' Error. No Rfrqrae node exists'')')
      stop
 115  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rfrqrae,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rfrqrae,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Icsi
      do n=1,narrays
        if (name(n) .eq. 'Icsi') goto 116
      enddo
      write(901,'('' Error. No Icsi node exists'')')
      stop
 116  continue
      call cg_array_read_as_f(n,Integer,icsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Icsf
      do n=1,narrays
        if (name(n) .eq. 'Icsf') goto 117
      enddo
      write(901,'('' Error. No Icsf node exists'')')
      stop
 117  continue
      call cg_array_read_as_f(n,Integer,icsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Jcsi
      do n=1,narrays
        if (name(n) .eq. 'Jcsi') goto 118
      enddo
      write(901,'('' Error. No Jcsi node exists'')')
      stop
 118  continue
      call cg_array_read_as_f(n,Integer,jcsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Jcsf
      do n=1,narrays
        if (name(n) .eq. 'Jcsf') goto 119
      enddo
      write(901,'('' Error. No Jcsf node exists'')')
      stop
 119  continue
      call cg_array_read_as_f(n,Integer,jcsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Kcsi
      do n=1,narrays
        if (name(n) .eq. 'Kcsi') goto 120
      enddo
      write(901,'('' Error. No Kcsi node exists'')')
      stop
 120  continue
      call cg_array_read_as_f(n,Integer,kcsi,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Kcsf
      do n=1,narrays
        if (name(n) .eq. 'Kcsf') goto 121
      enddo
      write(901,'('' Error. No Kcsf node exists'')')
      stop
 121  continue
      call cg_array_read_as_f(n,Integer,kcsf,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   First, get "factor" in case want to read every other point, etc
      if(idim.eq.idima .and. jdim.eq.jdima .and. kdim.eq.kdima) then
        nfac=1
      else if(jdim.eq.(jdima+1)/2 .and. kdim.eq.(kdima+1)/2) then
        nfac=2
      else if(jdim.eq.(jdima+3)/4 .and. kdim.eq.(kdima+3)/4) then
        nfac=4
      else if(jdim.eq.(jdima+7)/8 .and. kdim.eq.(kdima+7)/8) then
        nfac=8
      else if(jdim.eq.(jdima+15)/16 .and. kdim.eq.(kdima+15)/16) then
        nfac=16
      else
        write(901,'('' Error.  Desired grid level not supported'')')
        stop
      end if
      if(i2d .eq. 1) then
        nfaci=1
      else
        nfaci=nfac
      end if
c
c   Find out how many grids there are
      call cg_ngrids_f(iccg,ibase,igrid,ngrids,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (ngrids .ne. 3) then
        write(901,'('' Error.  Expecting 3 grids for deforming'',
     +   '' info'')')
        stop
      end if
      do n=1,ngrids
        call cg_grid_read_f(iccg,ibase,igrid,n,nameg(n),ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get MovedGrid info
      do n=1,ngrids
        if (nameg(n) .eq. 'MovedGrid') goto 122
      enddo
      write(901,'('' Error. No MovedGrid node exists'')')
      stop
 122  continue
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'GridCoordinates_t',n,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 3) then
        write(901,'('' Error.  Expecting exactly 3 arrays under'',
     +   '' MovedGrid.'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Get X
      do n=1,narrays
        if (name(n) .eq. 'CoordinateX') goto 123
      enddo
      write(901,'('' Error. No CoordinateX node exists'',
     +  '' under MovedGrid'')')
      stop
 123  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,x)
c   Get Y
      do n=1,narrays
        if (name(n) .eq. 'CoordinateY') goto 124
      enddo
      write(901,'('' Error. No CoordinateY node exists'',
     +  '' under MovedGrid'')')
      stop
 124  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,y)
      else
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,z)
      end if
c   Get Z
      do n=1,narrays
        if (name(n) .eq. 'CoordinateZ') goto 125
      enddo
      write(901,'('' Error. No CoordinateZ node exists'',
     +  '' under MovedGrid'')')
      stop
 125  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,z)
      else
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,-wk,y)
      end if
c
c   Get MovedGridLastDT info
      do n=1,ngrids
        if (nameg(n) .eq. 'MovedGridLastDT') goto 126
      enddo
      write(901,'('' Error. No MovedGridLastDT node exists'')')
      stop
 126  continue
      call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,
     + 'GridCoordinates_t',n,'end')
      if (ier .ne. 0) call cg_error_exit_f
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. 3) then
        write(901,'('' Error.  Expecting exactly 3 arrays under'',
     +   '' MovedGridLastDT.'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec,ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c   Get X
      do n=1,narrays
        if (name(n) .eq. 'CoordinateX') goto 127
      enddo
      write(901,'('' Error. No CoordinateX node exists'',
     +  '' under MovedGrid'')')
      stop
 127  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,xnm2)
c   Get Y
      do n=1,narrays
        if (name(n) .eq. 'CoordinateY') goto 128
      enddo
      write(901,'('' Error. No CoordinateY node exists'',
     +  '' under MovedGrid'')')
      stop
 128  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,ynm2)
      else
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,znm2)
      end if
c   Get Z
      do n=1,narrays
        if (name(n) .eq. 'CoordinateZ') goto 129
      enddo
      write(901,'('' Error. No CoordinateZ node exists'',
     +  '' under MovedGrid'')')
      stop
 129  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (ialph .eq. 0) then
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,wk,znm2)
      else
      call reorderg(nfac,nfaci,idima,jdima,kdima,
     +   idim,jdim,kdim,-wk,ynm2)
      end if
c
      return
      end

      subroutine waeromode(iccg,ibase,timekeep,xxn,gforcn,
     .         gforcnm,aehist,ncycmax,nmds,maxaes,
     .         ntt,naesrf,wk,maxnum,wk2,maxnum2)
c**********************************************************************
c     Purpose:  Write aeroelastic mode data.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      timekeep, xnn, gforcn, gforcnm, aehist
c           ............aeroelastic mode data (real)
c      ncycmax, nmds, maxaes, ntt
c           ............dimension arguments (integer)
c      naesrf...........number of aeroelastic surfaces (integer)
c      wk...............working space needed, of dimension
c                       (maxnum) or larger (real)
c      maxnum...........dimension of wk array (integer)
c      wk2..............working space needed, of dimension
c                       (maxnum2) or larger (real)
c      maxnum2..........dimension of wk2 array (integer)
c   OUTPUTS:
c      none
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=5)
c
      dimension timekeep(ntt),xxn(2*nmds,maxaes),
     + gforcn(2*nmds,maxaes),gforcnm(2*nmds,maxaes),
     + aehist(ncycmax,3,nmds,maxaes)
      dimension wk(maxnum),wk2(maxnum2)
      character name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Go to Base node, write CFL3DAeroModeData array as UserDefinedData
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_user_data_write_f('CFL3DAeroModeData', ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nuser_data_f(nuserdata, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .gt. numnames) then
        write(901,'('' Error... numnames too small in waeromode'')')
        stop
      end if
      do n=1,nuserdata
        call cg_user_data_read_f(n, name(n), ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(name(n) .eq. 'CFL3DAeroModeData') then
          nset=n
        end if
      enddo
      call cg_goto_f(iccg,ibase,ier,'UserDefinedData_t',nset,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Write data
c
      if (idouble .eq. 1) then
        call cg_array_write_f('Timekeep',RealDouble,1,ntt,timekeep,ier)
      else
        call cg_array_write_f('Timekeep',RealSingle,1,ntt,timekeep,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          wk(j)=xxn(n,i)
        enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('Xxn',RealDouble,1,j,wk,ier)
      else
        call cg_array_write_f('Xxn',RealSingle,1,j,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          wk(j)=gforcn(n,i)
        enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('Gforcn',RealDouble,1,j,wk,ier)
      else
        call cg_array_write_f('Gforcn',RealSingle,1,j,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          wk(j)=gforcnm(n,i)
        enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('Gforcnm',RealDouble,1,j,wk,ier)
      else
        call cg_array_write_f('Gforcnm',RealSingle,1,j,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      j=0
      do i=1,naesrf
        do n=1,nmds
          do m=1,3
            do l=1,ntt
              j=j+1
              wk2(j)=aehist(l,m,n,i)
            enddo
          enddo
        enddo
      enddo
      if (idouble .eq. 1) then
        call cg_array_write_f('Aehist',RealDouble,1,j,wk2,ier)
      else
        call cg_array_write_f('Aehist',RealSingle,1,j,wk2,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end

      subroutine raeromode(iccg,ibase,ncycmax,nmds,maxaes,
     .         ntt,naesrf,wk,maxnum,wk2,maxnum2,
     .         timekeep,xxn,gforcn,gforcnm,aehist)
c**********************************************************************
c     Purpose:  Read aeroelastic mode data.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside
c                       this routine) (integer)
c      ibase............CGNS base index number (determined outside
c                       this routine) (integer)
c      ncycmax, nmds, maxaes, ntt
c           ............dimension arguments (integer)
c      naesrf...........number of aeroelastic surfaces (integer)
c      wk...............working space needed, of dimension
c                       (maxnum) or larger (real)
c      maxnum...........dimension of wk array (integer)
c      wk2..............working space needed, of dimension
c                       (maxnum2) or larger (real)
c      maxnum2..........dimension of wk2 array (integer)
c   OUTPUTS:
c      timekeep, xnn, gforcn, gforcnm, aehist
c           ............aeroelastic mode data (real)
c**********************************************************************
c
#     include "cgnslib_f.h"
c
      parameter (numnames=5)
c
      dimension timekeep(ntt),xxn(2*nmds,maxaes),
     + gforcn(2*nmds,maxaes),gforcnm(2*nmds,maxaes),
     + aehist(ncycmax,3,nmds,maxaes)
      dimension wk(maxnum),wk2(maxnum2),idimvec(numnames)
      character name(numnames)*32
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c   Go to Base node:
      call cg_goto_f(iccg,ibase,ier,'end')
c   Find out how many UserDefinedData_t nodes there are
      call cg_nuser_data_f(nuserdata,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nuserdata .gt. numnames) then
        write(901,'('' Error... numnames not big enough'',
     +    '' in raeromode'')')
        write(901,'(''   increase numnames in raeromode'')')
        stop
      else
        do n=1,nuserdata
          call cg_user_data_read_f(n, name(n), ier)
          if (ier .ne. 0) call cg_error_exit_f
          if(name(n) .eq. 'CFL3DAeroModeData') goto 100
        enddo
        write(901,'('' Error. No CFL3DAeroModeData node exists'')')
        stop
 100    continue
        call cg_goto_f(iccg,ibase,ier,'UserDefinedData_t',n,'end')
        if (ier .ne. 0) call cg_error_exit_f
      end if
c   Find out how many arrays
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .gt. numnames) then
        write(901,'('' Error.  Too many arrays under'',
     +    '' CFL3DAeroModeData'')')
        write(901,'(''   increase numnames in raeromode'')')
        stop
      end if
      do n=1,narrays
        call cg_array_info_f(n,name(n),itype,idatadim,idimvec(n),ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c
c   Get Timekeep
      do n=1,narrays
        if (name(n) .eq. 'Timekeep') goto 101
      enddo
      write(901,'('' Error. No Timekeep node exists'')')
      stop
 101  continue
      if (idimvec(n) .ne. ntt) then
        write(901,'('' Error, Timekeep is wrong size'')')
        stop
      end if
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,timekeep,ier)
      else
        call cg_array_read_as_f(n,RealSingle,timekeep,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Xxn
      do n=1,narrays
        if (name(n) .eq. 'Xxn') goto 201
      enddo
      write(901,'('' Error. No Xxn node exists'')')
      stop
 201  continue
      if (idimvec(n) .ne. maxnum) then
        write(901,'('' Error, Xxn is wrong size'')')
        stop
      end if
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          xxn(n,i)=wk(j)
        enddo
      enddo
c
c   Get Gforcn
      do n=1,narrays
        if (name(n) .eq. 'Gforcn') goto 301
      enddo
      write(901,'('' Error. No Gforcn node exists'')')
      stop
 301  continue
      if (idimvec(n) .ne. maxnum) then
        write(901,'('' Error, Gforcn is wrong size'')')
        stop
      end if
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          gforcn(n,i)=wk(j)
        enddo
      enddo
c
c   Get Gforcnm
      do n=1,narrays
        if (name(n) .eq. 'Gforcnm') goto 401
      enddo
      write(901,'('' Error. No Gforcnm node exists'')')
      stop
 401  continue
      if (idimvec(n) .ne. maxnum) then
        write(901,'('' Error, Gforcnm is wrong size'')')
        stop
      end if
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      j=0
      do i=1,naesrf
        do n=1,2*nmds
          j=j+1
          gforcnm(n,i)=wk(j)
        enddo
      enddo
c
c   Get Aehist
      do n=1,narrays
        if (name(n) .eq. 'Aehist') goto 501
      enddo
      write(901,'('' Error. No Aehist node exists'')')
      stop
 501  continue
      if (idimvec(n) .ne. maxnum2) then
        write(901,'('' Error, Aehist is wrong size'')')
        stop
      end if
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,wk2,ier)
      else
        call cg_array_read_as_f(n,RealSingle,wk2,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      j=0
      do i=1,naesrf
        do n=1,nmds
          do m=1,3
            do l=1,ntt
              j=j+1
              aehist(l,m,n,i)=wk2(j)
            enddo
          enddo
        enddo
      enddo
c
      return
      end

#else

      subroutine dummycgns
c**********************************************************************
c     Purpose: Provide a dummy routine to compile if CGNS option is
c     not installed
c**********************************************************************
      return
      end
#endif
