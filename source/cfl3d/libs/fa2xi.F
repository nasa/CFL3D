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
      subroutine fa2xi(jf,kf,if,qif,js,ks,is,je,ke,ie,
     .                 ibctyp,w,iw,nw,niw,neta,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c      Purpose:  Accumulate fluxes in I-direction for use on twice
c     coarser mesh to ensure conservation.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension qif(jf,kf,5,4)
      dimension w(1)
      dimension iw(1)
      dimension ibctyp(2)
      common /sklton/ isklton
c
c     js,ks,is  - starting indices of coarser grid 
c                 defining twice finer embedded grid
c     je,ke,ie  - ending indices of coarser grid
c                 defining twice finer embedded grid
c
c     ibctyp(1) - mesh indicator for i=0 boundary
c     ibctyp(2) - mesh indicator for i=idim boundary
c
c     jf,kf,if  - grid dimension of twice finer mesh
c     qif       - array containing fine grid fluxes
c
c     w         - storage array for coarse grid flux
c     iw        - storage array for coarse grid
c                 starting and ending indices
c
c     nw        - counter for w
c     niw       - counter for iw (number of fine
c                 grid flux accumulations)
c     neta      - number of edges on which to accumulate
c
      jfm1 = jf-1
      kfm1 = kf-1
      ifm1 = if-1
c
      if (ibctyp(1).eq.21) then
c
c     left boundary
c
         mm        = 2
         iw(niw+1) = js
         iw(niw+2) = ks
                         iw(niw+3) = is
         iw(niw+4) = je
         iw(niw+5) = ke
                         iw(niw+6) = is
         iw(niw+7) = nw+1
         niw       = niw+7
         neta      = neta+1
c
         do 100 l=1,5
         do 100 k=1,kfm1,2
         do 100 j=1,jfm1,2
         nw    = nw+1
         w(nw) = qif(j,k,l,mm)
     .          +qif(j+1,k,l,mm)
     .          +qif(j,k+1,l,mm)
     .          +qif(j+1,k+1,l,mm)
  100    continue
         if (isklton.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)  '  summing flux at i=is  niw=',niw,
     .      ' nw=',nw
         end if
      end if
c
c
      if (ibctyp(2).eq.21) then
c
c     right boundary
c
         mm        = 4
         iw(niw+1) = js
         iw(niw+2) = ks
                         iw(niw+3)=ie
         iw(niw+4) = je
         iw(niw+5) = ke
                         iw(niw+6) = ie
         iw(niw+7) = nw+1
         niw       = niw+7
         neta      = neta+1
c
         do 200 l=1,5
         do 200 k=1,kfm1,2
         do 200 j=1,jfm1,2
         nw    = nw+1
         w(nw) = qif(j,k,l,mm)
     .          +qif(j+1,k,l,mm)
     .          +qif(j,k+1,l,mm)
     .          +qif(j+1,k+1,l,mm)
  200    continue
         if (isklton.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)  '  summing flux at i=ie  niw=',niw,
     .      ' nw=',nw
         end if
c
      end if
      return
      end
