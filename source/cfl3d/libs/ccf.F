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
      subroutine ccf(x,y,ca,sa,cl,xm,uf,wf,cf,pi)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Modify uf, wf, and cf (velocities and speed of sound) at
c     the farfield, based on point vortex correction (used when i2d=-1 and
c     farfield bctype=1003).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
c   far field boundary condition (point vortex)
c
      xe= ca*x+sa*y
      ye=-sa*x+ca*y
      beta=sqrt(1.e0-xm*xm)
      re=xe**2+(ye*beta)**2
      fact=cl*beta*xm/(re*4.e0*pi)
      ue=xm+fact*ye
      ve=  -fact*xe
      uf= ca*ue - sa*ve
      wf= sa*ue + ca*ve
      cf=sqrt(1.e0+0.2e0*(xm*xm-uf*uf-wf*wf))
      return
      end
