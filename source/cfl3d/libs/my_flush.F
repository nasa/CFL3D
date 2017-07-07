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
      subroutine my_flush(iunit)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Generic routine to flush an output buffer
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      integer iunit, istat
c
#   ifdef CRAY
c     no flush on the crays....
#   endif
c
#   ifdef IBM
      call flush_(iunit)
#   endif
c
#   ifdef SGI
      call flush(iunit,istat)
#   endif
c
#   ifdef SUN
      call flush(iunit)
#   endif
c
#   ifdef HP
      call flush(iunit)
#   endif
c
#   ifdef ALPHA
      call flush(iunit)
#   endif
c
#   ifdef INTEL
c     no flush on the intel...
#   endif
c
#   ifdef PG
      call flush(iunit)
#   endif
c
#   ifdef LAHEY
      call flush(iunit)
#   endif
c
#   ifdef SUN
      call flush(iunit)
#   endif
c
#   ifdef GENERIC
      call flush(iunit)
#   endif

c
      return
      end

