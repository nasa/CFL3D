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
      subroutine blkmax(jdim,kdim,idim,res,resmax,jm,km,im)
c
c     $Id$
c
c*************************************************************************
c     Purpose:  Determine the location of the maximum residual.
c*************************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5)
      resmax = res(1,1,1,1)
      jm     = 1
      km     = 1
      im     = 1
      do 10 j=1,jdim
      do 10 k=1,kdim
      do 10 i=1,idim-1
      if (real(res(j,k,i,1)).gt.real(resmax)) then
         resmax = res(j,k,i,1)
         jm = j
         km = k
         im = i
      end if
   10 continue
      return
      end
