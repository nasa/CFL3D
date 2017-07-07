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
      subroutine coll2q(q,qc,vol,volc,jdim,kdim,idim,jc,kc,ic,
     .                  res,qr,qw,js,ks,is,je,ke,ie,nbl,nblc,
     .                  vistf,vistc,tursavf,tursavc,nou,bou,nbuf,
     .                  ibufdim,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict q (the primative variables) with a volume-
c     weighted interpolation and residuals from finer embedded meshes
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
      dimension q(jdim,kdim,idim,5),qc(jc,kc,ic,5)
      dimension vol(jdim,kdim,idim-1),volc(jc,kc,ic-1)
      dimension qr(je-js+1,ke-ks+1,ie-is,5),qw(jdim,kdim,idim,5)
      dimension res(jdim,kdim,idim-1,5)
      dimension vistf(jdim,kdim,idim),vistc(jc,kc,ic),
     .          tursavf(jdim,kdim,idim,nummem),tursavc(jc,kc,ic,nummem)
c
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
c
c      restrict q(volume-weighted) from finer embedded mesh
c       to corresponding location in coarser mesh
c
c      restrict r from finer embedded mesh to coarser mesh
c
c      jdim,kdim,idim  finer mesh
c      jc,kc,ic        coarser mesh
c      js,ks,is        coarser mesh starting indices
c      je,ke,ie        coarser mesh ending indices
c
      if (isklton.gt.0)  then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl,nblc
      end if
    7 format(1x,39hrestricting variables and residual from,
     .       15h embedded block,i3,1x,16hto coarser block,i3)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 31 n=1,5
      nv    = jdim*kdim
      do 30 i=1,idim1
cdir$ ivdep
      do 30 izz=1,nv
      qw(izz,1,i,n) = q(izz,1,i,n)*vol(izz,1,i)
   30 continue
   31 continue
c
c      semi-coarsening / directional refinement
c
      nsi = (idim-1)/(ie-is)
      if (nsi.eq.2) then
c
      do 41 n=1,5
      kk = 0
      kq = ks-1
      do 40 k=1,kdim1,2
      kk = kk+1
      kq = kq+1
      ii = 0
      iq = is-1
      do 40 i=1,idim1,2
      ii = ii+1
      iq = iq+1
      jj = 0
      jq = js-1
      do 40 j=1,jdim1,2
      jj = jj+1
      jq = jq+1
      qc(jq,kq,iq,n)  = (qw(j,k,i,n)      +qw(j,k,i+1,n)
     .                  +qw(j+1,k,i,n)    +qw(j+1,k,i+1,n)
     .                  +qw(j,k+1,i,n)    +qw(j,k+1,i+1,n)
     .                  +qw(j+1,k+1,i,n)  +qw(j+1,k+1,i+1,n))/
     .                  (vol(j,k,i)      +vol(j,k,i+1)
     .                  +vol(j+1,k,i)    +vol(j+1,k,i+1)
     .                  +vol(j,k+1,i)    +vol(j,k+1,i+1)
     .                  +vol(j+1,k+1,i)  +vol(j+1,k+1,i+1))
      qr(jj,kk,ii,n)  = (res(j,k,i,n)      +res(j,k,i+1,n)
     .                  +res(j+1,k,i,n)    +res(j+1,k,i+1,n)
     .                  +res(j,k+1,i,n)    +res(j,k+1,i+1,n)
     .                  +res(j+1,k+1,i,n)  +res(j+1,k+1,i+1,n))
      if ((ivisc(1).gt.1 .and. n.eq.1) .or. (ivisc(2).gt.1 .and. n.eq.1)
     .   .or. (ivisc(3).gt.1 .and. n.eq.1)) then
         vistc(jq,kq,iq) = .125e0*(vistf(j,k,i)     +vistf(j,k,i+1)
     .                            +vistf(j+1,k,i)   +vistf(j+1,k,i+1)
     .                            +vistf(j,k+1,i)   +vistf(j,k+1,i+1)
     .                            +vistf(j+1,k+1,i) +vistf(j+1,k+1,i+1))
      end if
      if ((ivisc(1).eq.4 .and. n.eq.1) .or. (ivisc(2).eq.4 .and. n.eq.1)
     .   .or. (ivisc(3).eq.4 .and. n.eq.1) .or.
     .    (ivisc(1).eq.5 .and. n.eq.1) .or. (ivisc(2).eq.5 .and. n.eq.1)
     .   .or. (ivisc(3).eq.5 .and. n.eq.1)) then
         tursavc(jq,kq,iq,n)=.125e0*(tursavf(j,k,i,n)+tursavf(j,k,i+1,n)
     .                   +tursavf(j+1,k,i,n)   +tursavf(j+1,k,i+1,n)
     .                   +tursavf(j,k+1,i,n)   +tursavf(j,k+1,i+1,n)
     .                   +tursavf(j+1,k+1,i,n) +tursavf(j+1,k+1,i+1,n))
      end if
   40 continue
   41 continue
      if (ivisc(1).gt.5 .or. ivisc(2).gt.5 .or. ivisc(3).gt.5) then
      do 141 n=1,nummem
      kk = 0
      kq = ks-1
      do 140 k=1,kdim1,2
      kk = kk+1
      kq = kq+1
      ii = 0
      iq = is-1
      do 140 i=1,idim1,2
      ii = ii+1
      iq = iq+1
      jj = 0
      jq = js-1
      do 140 j=1,jdim1,2
      jj = jj+1
      jq = jq+1
      tursavc(jq,kq,iq,n)=.125e0*(tursavf(j,k,i,n)+tursavf(j,k,i+1,n)
     .                +tursavf(j+1,k,i,n)   +tursavf(j+1,k,i+1,n)
     .                +tursavf(j,k+1,i,n)   +tursavf(j,k+1,i+1,n)
     .                +tursavf(j+1,k+1,i,n) +tursavf(j+1,k+1,i+1,n))
  140 continue
  141 continue
      end if
c
      else
c
      do 46 n=1,5
      kk = 0
      kq = ks-1
      do 45 k=1,kdim1,2
      kk = kk+1
      kq = kq+1
      ii = 0
      iq = is-1
      do 45 i=1,idim1,1
      ii = ii+1
      iq = iq+1
      jj = 0
      jq = js-1
      do 45 j=1,jdim1,2
      jj = jj+1
      jq = jq+1
      qc(jq,kq,iq,n)  =( qw(j,k,i,n)    + qw(j+1,k,i,n) 
     .                 + qw(j,k+1,i,n)  + qw(j+1,k+1,i,n) )/
     .                 ( vol(j,k,i)     + vol(j+1,k,i)
     .                 + vol(j,k+1,i)   + vol(j+1,k+1,i)  )
      qr(jj,kk,ii,n)  =  res(j,k,i,n)   + res(j+1,k,i,n)
     .                 + res(j,k+1,i,n) + res(j+1,k+1,i,n)
      if ((ivisc(1).gt.1 .and. n.eq.1) .or. (ivisc(2).gt.1 .and. n.eq.1)
     .   .or. (ivisc(3).gt.1 .and. n.eq.1)) then
         vistc(jq,kq,iq) = .25e0*(vistf(j,k,i)      +vistf(j+1,k,i)
     .                           +vistf(j,k+1,i)    +vistf(j+1,k+1,i))
      end if
      if ((ivisc(1).eq.4 .and. n.eq.1) .or. (ivisc(2).eq.4 .and. n.eq.1)
     .   .or. (ivisc(3).eq.4 .and. n.eq.1) .or.
     .    (ivisc(1).eq.5 .and. n.eq.1) .or. (ivisc(2).eq.5 .and. n.eq.1)
     .   .or. (ivisc(3).eq.5 .and. n.eq.1)) then
         tursavc(jq,kq,iq,n)=.25e0*(tursavf(j,k,i,n)+tursavf(j+1,k,i,n)
     .                       +tursavf(j,k+1,i,n)+tursavf(j+1,k+1,i,n))
      end if
   45 continue
   46 continue
      if (ivisc(1).gt.5 .or. ivisc(2).gt.5 .or. ivisc(3).gt.5) then
      do 146 n=1,nummem
      kk = 0
      kq = ks-1
      do 145 k=1,kdim1,2
      kk = kk+1
      kq = kq+1
      ii = 0
      iq = is-1
      do 145 i=1,idim1,1
      ii = ii+1
      iq = iq+1
      jj = 0
      jq = js-1
      do 145 j=1,jdim1,2
      jj = jj+1
      jq = jq+1
      tursavc(jq,kq,iq,n)=.25e0*(tursavf(j,k,i,n)+tursavf(j+1,k,i,n)
     .                       +tursavf(j,k+1,i,n)+tursavf(j+1,k+1,i,n))
  145 continue
  146 continue
      end if
c
      end if
      return
      end
