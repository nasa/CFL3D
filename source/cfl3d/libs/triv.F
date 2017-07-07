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
      subroutine triv(jdim,kdim,jl,ju,kl,ku,x,a,b,c,f)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve a scalar tridiagonal system of equations
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(jdim,kdim),b(jdim,kdim),c(jdim,kdim)
      dimension x(jdim,kdim),f(jdim,kdim)
c
      do 10 j=jl,ju
      x(j,kl)=c(j,kl)/b(j,kl)
      f(j,kl)=f(j,kl)/b(j,kl)
 10   continue
      klp1 = kl +1
      do 1 i=klp1,ku
      do 20 j=jl,ju
         z=1./(b(j,i)-a(j,i)*x(j,i-1))
         x(j,i)=c(j,i)*z
         f(j,i)=(f(j,i)-a(j,i)*f(j,i-1))*z
 20   continue
1     continue
c
      kupkl=ku+kl
      do 2 i1=klp1,ku
         i=kupkl-i1
         do 30 j=jl,ju
         f(j,i)=f(j,i)-x(j,i)*f(j,i+1)
 30      continue
2     continue
c
      return
      end
