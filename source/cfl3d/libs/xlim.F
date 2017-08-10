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
      subroutine xlim(xkap,n,x1,x2,xc,iswi,npts,leq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform monotone interpolations to the interfaces
c     of the cells.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(n),x2(n),xc(n) 
c
      common /chk/ ichk
c
c     on input: 
c        x1 is left interface gradient del(-)
c        x2 is right interface gradient del(+) of cell
c        xc is cell center value
c     on output:  
c        x2 is gradient for right interface construction
c        x1 is gradient for left interface construction
c     for unlimited fully-upwind (kappa=-1;iswi=0)
c        x1(output)=x2(input)
c        x2(output)=x1(input)
c
      ibias = 1
      eps   = 1.0e-06
      phi   = (3.-xkap)/(1.-xkap+eps)
      eps2  = 0.5e-06
      zeroc = 0.
c
      if (real(xkap).lt.-1.0) then
c
c      first order
c
cdir$ ivdep
      do 8200 izz=1,n
      x2(izz) = 0.
      x1(izz) = 0.
 8200 continue
      return
      end if
c
c      bias gradients by average values (leq = 1 or 5)
c
      if (ibias.eq.1) then
      if (leq.eq.1 .or. leq.eq.5) then
cdir$ ivdep
      do 8000 izz=1,n
      x2avg   = xc(izz) + x2(izz)*0.5
      x1avg   = xc(izz) - x1(izz)*0.5
      x2(izz) = x2(izz)/x2avg
      x1(izz) = x1(izz)/x1avg
 8000 continue
      end if
      end if
c
      if (iswi.eq.-1) then
c
c     od - second order
c
cdir$ ivdep
      do 1000 izz=1,n
      x3      = ccabs(x1(izz))
      x4      = ccabs(x2(izz))
      x2(izz) = ccvmgt(x1(izz)*0.5e0,x2(izz),(real(x3).lt.real(x4)))
      x1(izz) = ccvmgt(x2(izz),x1(izz),(real(x3).lt.real(x4)))
      x1(izz) = x2(izz)
 1000 continue
      else if (iswi.eq.1) then
c
c      smooth (van albeda) limiter
c
      term    = xkap*4.e0
cdir$ ivdep
      do 1001 izz=1,n
      x3      = x1(izz)*x1(izz)+x2(izz)*x2(izz)+eps
      x4      = (eps2+x1(izz)*x2(izz))/x3
      x4      = x4*0.5e0
      x3      = (x1(izz)+x2(izz))*x4
      x4      = (x2(izz)-x1(izz))*x4*term*x4
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1001 continue
      else if (iswi.eq.0) then
c
c      unlimited  (kappa scheme)
c
      term    = xkap*0.25e0
cdir$ ivdep
      do 1002 izz=1,n
      x3      = .25e0*(x1(izz)+x2(izz))
      x4      = term*(x2(izz)-x1(izz))
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1002 continue
      else if (iswi.eq.2) then
c
c      limited  (min-mod scheme)
c
      term    = xkap*0.25e0
cdir$ ivdep
      do 1003 izz=1,n
      x4      = x1(izz)*x2(izz)
      x1(izz) = ccvmgt(zeroc,x1(izz),(real(x4).lt.0.e0))
      x2(izz) = ccvmgt(zeroc,x2(izz),(real(x4).lt.0.e0))
      x4      = phi*x1(izz)
      x4      = ccabs(x4)
      x3      = phi*x2(izz)
      x3      = ccabs(x3)
      x5      = ccabs(x2(izz))
      x2(izz) = ccvmgt(x1(izz)*phi,x2(izz),(real(x4).lt.real(x5)))
      x5      = ccabs(x1(izz))
      x1(izz) = ccvmgt(x2(izz)*phi,x1(izz),(real(x3).lt.real(x5)))
      x3      = .25e0*(x1(izz)+x2(izz))
      x4      = term*(x2(izz)-x1(izz))
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1003 continue
      else if (iswi.eq.3 .or. iswi.eq.4) then
c
c     tuned k=1/3 limiter - Spekreijse - Venkat                   AIAA-90-0429
c
      delx    = 10./float(npts)
      eps2    = delx**3
cdir$ ivdep
      do 7500 izz=1,n
      t3      = x1(izz)*x1(izz)
      t4      = x2(izz)*x2(izz)
      t5      = x1(izz)*x2(izz)
      t6      = x1(izz)+x2(izz)
      term    = 0.5*(t5+eps2)/(2.*(t3 + t4) - t5 + 3.*eps2)
      x2(izz) = (x2(izz)+t6)*term 
      x1(izz) = (x1(izz)+t6)*term 
 7500 continue
      end if
c
c      bias gradients by average values (leq = 1 or 5)
c
      if (ibias.eq.1) then
      if (leq.eq.1 .or. leq.eq.5) then
cdir$ ivdep
      do 9000 izz=1,n
      x2(izz) = x2(izz)*xc(izz)
      x1(izz) = x1(izz)*xc(izz)
 9000 continue
      end if
      end if
c
c     cap density and pressure
c     - ensures that they stay positive
c     - limits their maxima to their stagnation values
c
      if (ichk.eq.2) then
         call prolim(n,x1,x2,xc,leq)
      end if
      if (ichk.eq.3) then
         call prolim2(n,x1,x2,xc,leq)
      end if
c
      return
      end
