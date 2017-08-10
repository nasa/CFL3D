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
      subroutine collx(x,y,z,xx,yy,zz,jdim,kdim,idim,jj2,kk2,ii2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict x, y, and z values to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xx(jj2,kk2,ii2),yy(jj2,kk2,ii2),zz(jj2,kk2,ii2)
c
      ii   = 0
      iinc = 2
      if (idim.eq.2) iinc = 1
      do 10 i=1,idim,iinc
      ii   = ii+1
      kk   = 0
      do 10 k=1,kdim,2
      kk   = kk+1
      jj   = 0
      do 10 j=1,jdim,2
      jj   = jj+1
      xx(jj,kk,ii) = x(j,k,i)
      yy(jj,kk,ii) = y(j,k,i)
      zz(jj,kk,ii) = z(j,k,i)
   10 continue
      return
      end
