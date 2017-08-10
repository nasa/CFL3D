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
      subroutine fill(jdim,kdim,idim,q,ll)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Fills the edges of the q array for safety using
c     multi-plane vectorization technique.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,ll)
c
c     fill edges of q array for safety
c
      idim1 = idim-1
      kdim1 = kdim-1
      do 3227 l=1,ll
      do 3226 i=1,idim1
cdir$ ivdep
      do 3229 k=1,kdim1
      q(jdim,k,i,l) = q(jdim-1,k,i,l)
 3229 continue
      do 1000 izz=1,jdim
      q(izz,kdim,i,l) = q(izz,kdim-1,i,l)
 1000 continue
 3226 continue
      n = jdim*kdim
cdir$ ivdep
      do 1001 izz=1,n
      q(izz,1,idim,l) = q(izz,1,idim-1,l)
 1001 continue
 3227 continue
      return
      end
