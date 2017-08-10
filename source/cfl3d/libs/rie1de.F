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
      subroutine rie1de(jvdim,t,jv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine far-field boundary data using quasi 1-d 
c     characteristic relations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension t(jvdim,25)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
c
c     unsteady quasi-1-d characteristic boundary conditions
c     for inflow/outflow
c
cdir$ ivdep
      do 3117 izz=1,jv
      t(izz,10) = gamma*t(izz,5)/t(izz,1)
 3117 continue
c
      x2gm1=2.e0/gm1
cdir$ ivdep
      do 3118 izz=1,jv
      t(izz,11) = t(izz,6)*t(izz,2)+t(izz,7)*t(izz,3)+t(izz,8)*t(izz,4)
      t(izz,12) = sqrt(t(izz,10))
      t(izz,12) = t(izz,11)+x2gm1*t(izz,12)
      t(izz,13) = t(izz,6)*t(izz,21)+t(izz,7)*t(izz,22)
     .          + t(izz,8)*t(izz,23)
      t(izz,14) = t(izz,13)-x2gm1*t(izz,24)
 3118 continue
c
      x4gm1 = 0.25e0*gm1
cdir$ ivdep
      do 3119 izz=1,jv
      t(izz,15) = .5e0*(t(izz,12)+t(izz,14))
      t(izz,16) = x4gm1*(t(izz,12)-t(izz,14))
      t(izz,17) = t(izz,15)-t(izz,13)
c  put unsteady metrics into t(5)
      t(izz,5)  = t(izz,20)
      t(izz,18) = t(izz,21)+t(izz,6)*t(izz,17)
      t(izz,19) = t(izz,22)+t(izz,7)*t(izz,17)
      t(izz,20) = t(izz,23)+t(izz,8)*t(izz,17)
 3119 continue
c
cdir$ ivdep
      do 3120 izz=1,jv
      t(izz,13) = t(izz,15)+t(izz,5)
      t(izz,17) = ccvmgt(t(izz,15)-t(izz,11),t(izz,17),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,18) = ccvmgt(t(izz,2)+t(izz,6)*t(izz,17),t(izz,18),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,19) = ccvmgt(t(izz,3)+t(izz,7)*t(izz,17),t(izz,19),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,20) = ccvmgt(t(izz,4)+t(izz,8)*t(izz,17),t(izz,20),
     .                  (real(t(izz,13)).ge.0.e0))
 3120 continue
c
      xgm1 = 1.e0/gm1
cdir$ ivdep
      do 3121 izz=1,jv
      rmloc2=(t(izz,18)**2+t(izz,19)**2+t(izz,20)**2)/
     .       (t(izz,16)**2)
      t(izz,5)  = t(izz,25)/(gamma*(1.0+0.5*gm1*rmloc2)**(gamma/gm1))
      t(izz,1)  = gamma*t(izz,5)/(t(izz,16)**2)
      t(izz,2)  = t(izz,18)
      t(izz,3)  = t(izz,19)
      t(izz,4)  = t(izz,20)
 3121 continue
      return
      end
