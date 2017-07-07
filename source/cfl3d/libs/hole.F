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
      subroutine hole(i,npl,jdim,kdim,idim,res,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Zero out the rhs residuals for the blanked points.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension blank(jdim,kdim,idim)
      dimension res(jdim,kdim,idim-1,5) 
c
      nn = jdim*kdim
      do 130 ipl=1,npl
      ii = i+ipl-1
      do 100 l=1,5
cdir$ ivdep
      do 101 n=1,nn 
      res(n,1,ii,l) = res(n,1,ii,l)*blank(n,1,ii)
  101 continue
  100 continue
  130 continue
      return
      end 
