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
      subroutine get_cmd_args(ierr)
c
c     $Id$
c
c**********************************************************************
c     Purpose: Interrogate the command line for specific options.
c
c     Usage:
c
c     call get_cmd_args (ierr)
c
c     ierr     : (out) 0 if no errors, 1 if some error occurred.
c
c     Original coding: Thomas R. Faulkner, MRJ Technology Solutions
c**********************************************************************

      character*8 arg
c
      common /proces/ numprocs
c
c     Begin get_cmd_args.
c
c     Initialize common data.
c
      numprocs = -1
c
c     Get the total number of command line arguments.
c
      num_args = iargc()
      if (num_args .eq. 2) then
#if defined CRAY_TIME
         call pxfgetarg (1, arg, 3, ierror)
#else
         call getarg (1, arg)
#endif
         if (arg .eq. "-np" .or. arg .eq. "-nt" .or. arg .eq. "-n") then
#if defined CRAY_TIME
            call pxfgetarg (2, arg, 4, ierror)
#else
            call getarg (2, arg)
#endif
            read (arg, 100, err = 110) n
 100        format (i3)
            numprocs = n
            ierr = 0
            go to 111
 110        continue
            ierr = 1
 111        continue
         else
            ierr = 1
         end if
      else
         ierr = 1
      end if
      return
      end
