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
      subroutine expand(mdim,ndim,msub1,jmax1,kmax1,l,x1,y1,z1,xte,yte,
     .                  zte,factjlo,factjhi,factklo,factkhi,
     .                  jmax2,kmax2,x2,y2,z2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Expand grid at boundaries.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c    
c      expand "from"  grid(s) at boundaries to insure that the 
c      "to" grid is completely covered
c
      dimension x1(mdim,ndim),y1(mdim,ndim),z1(mdim,ndim)
      dimension x2(mdim,ndim),y2(mdim,ndim),z2(mdim,ndim)
      dimension xte(mdim,ndim,msub1),yte(mdim,ndim,msub1),
     .          zte(mdim,ndim,msub1)
c
      sk1 = 0.
      sk2 = 0.
      do 10 j=2,jmax1
      sk1 = sk1+sqrt((x1(j,1)-x1(j-1,1))**2+(y1(j,1)-y1(j-1,1))**2
     .              +(z1(j,1)-z1(j-1,1))**2)
      sk2 = sk2+sqrt((x1(j,kmax1)-x1(j-1,kmax1))**2
     .              +(y1(j,kmax1)-y1(j-1,kmax1))**2
     .              +(z1(j,kmax1)-z1(j-1,kmax1))**2)
   10 continue
      sj1 = 0.
      sj2 = 0.
      do 20 k=2,kmax1
      sj1 = sj1+sqrt((x1(1,k)-x1(1,k-1))**2+(y1(1,k)-y1(1,k-1))**2
     .              +(z1(1,k)-z1(1,k-1))**2)
      sj2 = sj2+sqrt((x1(jmax1,k)-x1(jmax1,k-1))**2
     .              +(y1(jmax1,k)-y1(jmax1,k-1))**2
     .              +(z1(jmax1,k)-z1(jmax1,k-1))**2)
   20 continue
      factk1 = factklo*sk1
      factk2 = factkhi*sk2
      factj1 = factjlo*sj1
      factj2 = factjhi*sj2
c
      do 210 j=1,jmax1
c      dx = .5*(-x1(j,3)-3.*x1(j,1)+4.*x1(j,2))
c      dy = .5*(-y1(j,3)-3.*y1(j,1)+4.*y1(j,2))
c      dz = .5*(-z1(j,3)-3.*z1(j,1)+4.*z1(j,2))
      dx = x1(j,2)-x1(j,1)
      dy = y1(j,2)-y1(j,1)
      dz = z1(j,2)-z1(j,1)
      ds = sqrt(dx*dx+dy*dy+dz*dz)
      if(real(ds).le.0)ds = 1.0
      xte(j+1,1,l) = x1(j,1)-factk1*dx/ds
      yte(j+1,1,l) = y1(j,1)-factk1*dy/ds
      zte(j+1,1,l) = z1(j,1)-factk1*dz/ds
  210 continue
c
      do 220 j=1,jmax1
c      dx = .5*(x1(j,kmax1-2)-4.*x1(j,kmax1-1)+3.*x1(j,kmax1))
c      dy = .5*(y1(j,kmax1-2)-4.*y1(j,kmax1-1)+3.*y1(j,kmax1))
c      dz = .5*(z1(j,kmax1-2)-4.*z1(j,kmax1-1)+3.*z1(j,kmax1))
      dx = x1(j,kmax1)-x1(j,kmax1-1)
      dy = y1(j,kmax1)-y1(j,kmax1-1)
      dz = z1(j,kmax1)-z1(j,kmax1-1)
      ds = sqrt(dx*dx+dy*dy+dz*dz)
      if(real(ds).le.0)ds = 1.0
      xte(j+1,kmax1+2,l) = x1(j,kmax1)+factk2*dx/ds
      yte(j+1,kmax1+2,l) = y1(j,kmax1)+factk2*dy/ds
      zte(j+1,kmax1+2,l) = z1(j,kmax1)+factk2*dz/ds
  220 continue
c 
      do 250 k=2,kmax1+1
      do 250 j=1,jmax1
      xte(j+1,k,l) = x1(j,k-1)
      yte(j+1,k,l) = y1(j,k-1)
      zte(j+1,k,l) = z1(j,k-1)
  250 continue
c
      do 310 k=1,kmax1+2
c      dx = .5*(4.*xte(3,k,l)-3.*xte(2,k,l)-xte(4,k,l))
c      dy = .5*(4.*yte(3,k,l)-3.*yte(2,k,l)-yte(4,k,l))
c      dz = .5*(4.*zte(3,k,l)-3.*zte(2,k,l)-zte(4,k,l))
      dx = xte(3,k,l)-xte(2,k,l)
      dy = yte(3,k,l)-yte(2,k,l)
      dz = zte(3,k,l)-zte(2,k,l)
      ds = sqrt(dx*dx+dy*dy+dz*dz)
      if(real(ds).le.0)ds = 1.0
      xte(1,k,l) = xte(2,k,l)-factj1*dx/ds
      yte(1,k,l) = yte(2,k,l)-factj1*dy/ds
      zte(1,k,l) = zte(2,k,l)-factj1*dz/ds
  310 continue
c
      do 320 k=1,kmax1+2
c      dx = .5*(3.*xte(jmax1+1,k,l)-4.*xte(jmax1,k,l)+xte(jmax1-1,k,l))
c      dy = .5*(3.*yte(jmax1+1,k,l)-4.*yte(jmax1,k,l)+yte(jmax1-1,k,l))
c      dz = .5*(3.*zte(jmax1+1,k,l)-4.*zte(jmax1,k,l)+zte(jmax1-1,k,l))
      dx = xte(jmax1+1,k,l)-xte(jmax1,k,l)
      dy = yte(jmax1+1,k,l)-yte(jmax1,k,l)
      dz = zte(jmax1+1,k,l)-zte(jmax1,k,l)
      ds = sqrt(dx*dx+dy*dy+dz*dz)
      if(real(ds).le.0)ds = 1.0
      xte(jmax1+2,k,l) = xte(jmax1+1,k,l)+factj2*dx/ds
      yte(jmax1+2,k,l) = yte(jmax1+1,k,l)+factj2*dy/ds
      zte(jmax1+2,k,l) = zte(jmax1+1,k,l)+factj2*dz/ds
  320 continue
      return
      end
