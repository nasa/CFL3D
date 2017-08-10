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
      subroutine dthole(jdim,kdim,idim,dtj,vol,blank,dtmin,
     .                  nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Update the dt values for the hole and fringe cells;
c     the values will be replaced by dtmin.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension blank(jdim,kdim,idim)
      dimension dtj(jdim,kdim,idim-1),vol(jdim,kdim,idim-1) 
c
c     nou(1) = min(nou(1)+1,ibufdim)
c     if (isklton.gt.0) then
c        nou(1) = min(nou(1)+1,ibufdim)
c        write(bou(nou(1),1),*)' in dthole, dtmin= ',real(dtmin)
c     end if
c
      kn = jdim*kdim*(idim-1)
cdir$ ivdep
      do 12 n=1,kn
      dtj(n,1,1) = ccvmgt(vol(n,1,1)/dtmin,dtj(n,1,1),
     .                  (real(blank(n,1,1)).eq.0.e0))
   12 continue
      return
      end 
