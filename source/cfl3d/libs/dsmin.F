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
      subroutine dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
     .                 j1,k1,l1,lout,lflag,xif1,xif2,etf1,etf2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Find closest point in grid to point (xc,yc,zc).
c     lflag > 0 : search over all "from" blocks not yet searched.
c                (lout(l)>0 indicates block l has already been searched)
c     lflag < 0 : search only block l1
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,nsub),y(jdim,kdim,nsub),z(jdim,kdim,nsub)
      dimension jjmax(nsub),kkmax(nsub)
      integer   lout(nsub),xif1(nsub),xif2(nsub),etf1(nsub),
     .          etf2(nsub)
c
      dmin = 1.0e+20
c
      ls = 1
      le = lmax
      if(lflag.lt.0) then
        ls = l1
        le = l1
      end if
c
      do 1235 l=ls,le
c     
c     skip over blocks already tried
c
      if(ls.ne.le) then
        if (lout(l).gt.0) go to 1235
      end if
c
c     search only over specified range, and only inside non-expanded limits
c
      js = xif1(l)
      je = xif2(l)
      ks = etf1(l)
      ke = etf2(l)
      je = min(je-2,jjmax(l)-2)
      ke = min(ke-2,kkmax(l)-2)
      js = max(js+1,2)
      ks = max(ks+1,2)
c
      do 1234 k=ks,ke
      do 1234 j=js,je
      d1 = (xc-x(j,k,l))**2 + (yc-y(j,k,l))**2 + (zc-z(j,k,l))**2
      if (real(d1).lt.real(dmin)) then
         j1   = j
         k1   = k
         l1   = l
         dmin = d1
      end if
 1234 continue
 1235 continue
      return
      end 
