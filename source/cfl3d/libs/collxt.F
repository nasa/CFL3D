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
      subroutine collxt(xt,xtt,jdim,kdim,idim,jj2,kk2,ii2,nbl,
     .                  nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict xt (array containing grid speeds)
c     to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension xt(jdim,kdim,idim,3),xtt(jj2,kk2,ii2,3)
c
      common /sklton/ isklton
c
c      restrict xt  to coarser mesh
c
c      jdim,kdim,idim  finer mesh
c      jj2,kk2,ii2     coarser mesh
c
      nbl1  = nbl+1
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl,nbl1
      end if
    7 format(1x,38h    restricting grid speeds from finer,
     .       6h block,i4,1x,16hto coarser block,i4)
c
      do 10 ll=1,3
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
      xtt(jj,kk,ii,ll) = xt(j,k,i,ll)
   10 continue
      return
      end
