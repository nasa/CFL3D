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
      subroutine collqc0(qc0,qc0c,vol,volc,jdim,kdim,idim,jj2,kk2,ii2,
     .                  dqc0,dqc0c,qw,nbl,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict conservative variables Q(n) and Q(n)-Q(n-1)
c     to coarser meshes via summation over fine-grid cells for use
c     in time-accurate multigrid.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension qc0(jdim,kdim,idim-1,5),qc0c(jj2,kk2,ii2-1,5),
     .          qw(jdim,kdim,idim,5)
      dimension dqc0(jdim,kdim,idim-1,5),dqc0c(jj2,kk2,ii2-1,5)
      dimension vol(jdim,kdim,idim-1),volc(jj2,kk2,ii2-1)
c
      common /sklton/ isklton
c
c     qc0.....Q(n)
c     dqc0....Q(n)-Q(n-1)
c
c     qc0,dqc0,jdim,kdim,idim  finer mesh
c     qc0c,dqc0c,jj2,kk2,ii2   coarser mesh
c
      nbl1  = nbl+1
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl,nbl1
      end if
    7 format(1x,45hrestricting qc0       and     dqc0 from finer,
     .       6h block,i4,1x,16hto coarser block,i4)
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
      jjl   = jj2-1
      kkl   = kk2-1
      iil   = ii2-1
      n     = jj2*kk2
c
c     restrict Q(n)
c
      do 60 n=1,5
      if (idim.gt.2) then
      kk = 0
      do 40 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 40 j=1,jdim1,2
      jj = jj+1
      ii = 0
      do 40 i=1,idim1,2
      ii = ii+1
      qc0c(jj,kk,ii,n)  = qc0(j,k,i,n)     + qc0(j,k,i+1,n)
     .                  + qc0(j+1,k,i,n)   + qc0(j+1,k,i+1,n)
     .                  + qc0(j,k+1,i,n)   + qc0(j,k+1,i+1,n)
     .                  + qc0(j+1,k+1,i,n) + qc0(j+1,k+1,i+1,n)
   40 continue
c
      else
c
      ii = 1
      i  = 1
      kk = 0
      do 404 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 404 j=1,jdim1,2
      jj = jj+1
      qc0c(jj,kk,ii,n)  = qc0(j,k,i,n)   + qc0(j+1,k,i,n)    
     .                  + qc0(j,k+1,i,n) + qc0(j+1,k+1,i,n)
  404 continue
      end if
   60 continue
c
c     restrict Q(n)-Q(n-1)
c
      do 61 n=1,5
      if (idim.gt.2) then
      kk = 0
      do 41 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 41 j=1,jdim1,2
      jj = jj+1
      ii = 0
      do 41 i=1,idim1,2
      ii = ii+1
      dqc0c(jj,kk,ii,n)  = dqc0(j,k,i,n)     + dqc0(j,k,i+1,n)
     .                   + dqc0(j+1,k,i,n)   + dqc0(j+1,k,i+1,n)
     .                   + dqc0(j,k+1,i,n)   + dqc0(j,k+1,i+1,n)
     .                   + dqc0(j+1,k+1,i,n) + dqc0(j+1,k+1,i+1,n)
   41 continue
c
      else
c
      ii = 1
      i  = 1
      kk = 0
      do 414 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 414 j=1,jdim1,2
      jj = jj+1
      dqc0c(jj,kk,ii,n)  = dqc0(j,k,i,n)   + dqc0(j+1,k,i,n)
     .                   + dqc0(j,k+1,i,n) + dqc0(j+1,k+1,i,n)
  414 continue
      end if
   61 continue
      return
      end
