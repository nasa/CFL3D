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
      program cgns_readhist
c
c     $Id$
c
c***********************************************************************
c   Purpose: reads in CGNS file and writes out a formatted file with history
c   information (as a function of iteration)
c***********************************************************************
c
#ifdef CGNS
#     include "cgnslib_f.h"
      character*80 file2
c
c  Read CGNS file
      write(6,'('' Note:  currently, you must be linked to CGNS'',
     +  '' V2.5.2 or later'',/)')
      write(6,'('' input name of CGNS file to read'')')
      read(5,'(a80)') file2
      call cg_open_f(file2,CG_MODE_READ,iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nbases_f(iccg,nbases,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nbases .ne. 1) then
        write(6,'('' Error nbases ='',i5,'' (should be 1)'')') nbases
        write(6,'('' stopping'')')
        stop
      end if
      ibase=1
c
c   Goto base index
      call cg_goto_f(iccg,ibase,ier,'end')
      if (ier .ne. 0) call cg_error_exit_f
c
c   Read number of iterations of data in ConvergenceHistory node
      call cg_convergence_read_f(ncycmax,'',ier)
      write(6,'('' Database contains '',i7,'' iterations'')') ncycmax
      if (ier .ne. 0) then
         write(6,'('' ConvergenceHistory node does not exist.'')')
         write(6,'(''   Stopping.'')')
         stop
      end if
c
c   Call subroutine to allocate memory, read data, and write file
      call readdata(iccg,ibase,ncycmax)
c
c   Close data base
      call cg_close_f(iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      stop
      end

c**************************************************************************

      subroutine readdata(iccg,ibase,ncycmax)
c
      character*80 file1
c
      integer stats
c
      allocatable :: cdpw(:)
      allocatable :: cdvw(:)
      allocatable :: cdw(:)
      allocatable :: cftmomw(:)
      allocatable :: cftpw(:)
      allocatable :: cfttotw(:)
      allocatable :: cftvw(:)
      allocatable :: clw(:)
      allocatable :: cmxw(:)
      allocatable :: cmyw(:)
      allocatable :: cmzw(:)
      allocatable :: cxw(:)
      allocatable :: cyw(:)
      allocatable :: czw(:)
      allocatable :: fmdotw(:)
      allocatable :: nneg1(:)
      allocatable :: nneg2(:)
      allocatable :: rms(:)
      allocatable :: rmstr1(:)
      allocatable :: rmstr2(:)
c
c     allocate memory
c
      memuse = 0
      allocate( cdpw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdpw',memuse,stats)
      allocate( cdvw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdvw',memuse,stats)
      allocate( cdw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdw',memuse,stats)
      allocate( cftmomw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftmomw',memuse,stats)
      allocate( cftpw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftpw',memuse,stats)
      allocate( cfttotw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cfttotw',memuse,stats)
      allocate( cftvw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftvw',memuse,stats)
      allocate( clw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'clw',memuse,stats)
      allocate( cmxw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmxw',memuse,stats)
      allocate( cmyw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmyw',memuse,stats)
      allocate( cmzw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmzw',memuse,stats)
      allocate( cxw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cxw',memuse,stats)
      allocate( cyw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cyw',memuse,stats)
      allocate( czw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'czw',memuse,stats)
      allocate( fmdotw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'fmdotw',memuse,stats)
      allocate( nneg1(ncycmax), stat=stats )
      call umalloc_r(ncycmax,1,'nneg1',memuse,stats)
      allocate( nneg2(ncycmax), stat=stats )
      call umalloc_r(ncycmax,1,'nneg2',memuse,stats)
      allocate( rms(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rms',memuse,stats)
      allocate( rmstr1(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rmstr1',memuse,stats)
      allocate( rmstr2(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rmstr2',memuse,stats)
c
c   Get history information
      call readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     + cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     + cftpw,cftvw,cfttotw,rmstr1,rmstr2,nneg1,nneg2)
c
c   Write out data to a formatted file
c   NOTE:  modify here if you want anything written other than
c          log(residual), cl, and cd
      write(6,'('' what formatted filename do you want to write to?'')')
      read(5,'(a80)') file1
      write(6,'('' input 0 for short output (it,res,cl,cd,cdp,cdv),'',
     . '' 1 for long output:'')')
      read(5,*) ilong
      open(2,file=file1,form='formatted',status='unknown')
      if (ilong .eq. 0) then
        write(2,'(''variables="iter","log(res)","cl","cd","cdp",'',
     .   ''"cdv"'')')
      else
        write(2,'(''variables="iter","log(res)","cl","cd","cdp",'',
     .   ''"cdv","cx","cy","cz","cmx","cmy","cmz","fmdot","cftmom",'',
     .   ''"cftp","cftv","cfttot","rmstr1","rmstr2","nneg1",'',
     .   ''"nneg2"'')')
      end if
      write(6,'('' ntr='',i7)') ntr
      do n=1,ntr
        if(rms(n) .le. 0.) rms(n)=1.
        rms(n)=log10(rms(n))
        if (ilong .eq. 0) then
        write(2,'(i7,5e15.5)') n,rms(n),clw(n),cdw(n),cdpw(n),
     .   cdvw(n)
        else
        write(2,'(i7,18e15.5,2i7)') n,rms(n),clw(n),cdw(n),cdpw(n),
     .   cdvw(n),cxw(n),cyw(n),czw(n),cmxw(n),cmyw(n),cmzw(n),
     .   fmdotw(n),cftmomw(n),cftpw(n),cftvw(n),cfttotw(n),
     .   rmstr1(n),rmstr2(n),nneg1(n),nneg2(n)
        end if
      enddo
      write(6,'('' History data written to file '',a80)') file1
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(rms)
         deallocate(clw)
         deallocate(cdw)
         deallocate(cdpw)
         deallocate(cdvw)
         deallocate(cxw)
         deallocate(cyw)
         deallocate(czw)
         deallocate(cmxw)
         deallocate(cmyw)
         deallocate(cmzw)
         deallocate(fmdotw)
         deallocate(cftmomw)
         deallocate(cftpw)
         deallocate(cftvw)
         deallocate(cfttotw)
         deallocate(rmstr1)
         deallocate(rmstr2)
         deallocate(nneg1)
         deallocate(nneg2)
      end if
c
      return
      end

c**************************************************************************

      subroutine readhist(iccg,ibase,ncycmax,ntr,rms,clw,cdw,
     +   cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,fmdotw,cftmomw,
     +   cftpw,cftvw,cfttotw,rmstr1,rmstr2,nneg1,nneg2)
c
c   Reads history info from CGNS file.  
c      (NOTE:  cg_array_read_as_f may eventually be converted
c      to cg_array_read_f)
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this routine) (integer)
c      ibase............CGNS base index number (determined outside this routine) (integer)
c      ncycmax..........max number of iterations that can be saved (integer)
c   OUTPUTS:
c      ntr..............no of iterations (integer)
c      rms..............rms of density residual (real)
c      clw,cdw,...,rmstr2......other history variables (real)
c      nneg1,nneg2.............neg pts in turb solution (integer)
c
#     include "cgnslib_f.h"
c
      parameter (numnames=20)
c
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),
     +  cdpw(ncycmax),cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),
     +  czw(ncycmax),cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     +  fmdotw(ncycmax),cftmomw(ncycmax),cftpw(ncycmax),
     +  cftvw(ncycmax),cfttotw(ncycmax),rmstr1(ncycmax),
     +  rmstr2(ncycmax),nneg1(ncycmax),nneg2(ncycmax)
      character*32 name(numnames)
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
      call cg_convergence_read_f(ntr,'',ier)
      if (ier .ne. 0) then
         write(6,'('' ConvergenceHistory node does not exist.'')')
         write(6,'(''   setting ntr=0 and continuing'')')
         ntr=0
         return
      end if
      write(6,'('' reading convergence history data... ntr='',i5)') ntr
      if (ntr .gt. ncycmax) then
      write(6,1239)
 1239 format(/,1x,11hstopping...,
     .       40hprevious number of iterations computed >,
     .       1x,18h dimension ncycmax)
      write(6,*)' ntr,ncycmax = ',ntr,ncycmax
      write(6,*)' increase value of ncycmax to at LEAST ',
     .ntr+ncycmax
      stop
      end if
c
c   Get Information about what is under History node
      call cg_goto_f(iccg,ibase,ier,'ConvergenceHistory_t',1,'end')
      if (ier .ne. 0) call cg_error_exit_f
      call cg_narrays_f(narrays,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (narrays .ne. numnames) then
        write(6,'('' Incorrect number of history arrays.'')')
         write(6,'(''   setting ntr=0 and continuing'')')
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
      write(6,'('' Error. No RSDMassRMS node exists'')')
      stop
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
      write(6,'('' Error. No CoefLift node exists'')')
      stop
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
      write(6,'('' Error. No CoefDrag node exists'')')
      stop
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
      write(6,'('' Error. No CoefPressureDrag node exists'')')
      stop
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
      write(6,'('' Error. No CoefViscousDrag node exists'')')
      stop
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
      write(6,'('' Error. No CoefForceX node exists'')')
      stop
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
      write(6,'('' Error. No CoefForceY node exists'')')
      stop
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
      write(6,'('' Error. No CoefForceZ node exists'')')
      stop
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
      write(6,'('' Error. No CoefMomentX node exists'')')
      stop
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
      write(6,'('' Error. No CoefMomentY node exists'')')
      stop
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
      write(6,'('' Error. No CoefMomentZ node exists'')')
      stop
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
      write(6,'('' Error. No ControlSurfaceMassFlow node exists'')')
      stop
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
      write(6,'('' Error. No ControlSurfaceMomentumForce node'',
     +  '' exists'')')
      stop
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
      write(6,'('' Error. No ControlSurfacePressureForce node'',
     +  '' exists'')')
      stop
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
      write(6,'('' Error. No ControlSurfaceViscousForce node'',
     +  '' exists'')')
      stop
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
      write(6,'('' Error. No ControlSurfaceForce node'',
     +  '' exists'')')
      stop
 116  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,cfttotw,ier)
      else
        call cg_array_read_as_f(n,RealSingle,cfttotw,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Residual of turb eqn 1
      do n=1,narrays
        if (name(n) .eq. 'RSDTurbEqn1RMS') goto 117
      enddo
      write(6,'('' Error. No RSDTurbEqn1RMS node exists'')')
      stop
 117  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rmstr1,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rmstr1,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get Residual of turb eqn 2
      do n=1,narrays
        if (name(n) .eq. 'RSDTurbEqn2RMS') goto 118
      enddo
      write(6,'('' Error. No RSDTurbEqn2RMS node exists'')')
      stop
 118  continue
      if (idouble .eq. 1) then
        call cg_array_read_as_f(n,RealDouble,rmstr2,ier)
      else
        call cg_array_read_as_f(n,RealSingle,rmstr2,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get turb eqn 1 bad points
      do n=1,narrays
        if (name(n) .eq. 'BadPointsTurbEqn1') goto 119
      enddo
      write(6,'('' Error. No BadPointsTurbEqn1 node exists'')')
      stop
 119  continue
      call cg_array_read_as_f(n,Integer,nneg1,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
c   Get turb eqn 2 bad points
      do n=1,narrays
        if (name(n) .eq. 'BadPointsTurbEqn2') goto 120
      enddo
      write(6,'('' Error. No BadPointsTurbEqn2 node exists'')')
      stop
 120  continue
      call cg_array_read_as_f(n,Integer,nneg2,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
#else
c     this is now just a dummy routine since CGNS libs have not been
c     used
c
      write(6,*)'This code is non-functional since the installation'
      write(6,*)'of cfl3d was done without cgns libraries'
c
      stop
      end
#endif
