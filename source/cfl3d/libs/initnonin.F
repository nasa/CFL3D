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
      subroutine initnonin(nbl,jdim,kdim,idim,q,qj0,qk0,qi0,
     .                vol,volj0,volk0,voli0,
     .                x,y,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Increment the initial conditions on a mesh to include
c     the rotational component for noninertial calculations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension volj0(kdim,idim-1,4),volk0(jdim,idim-1,4),
     .          voli0(jdim,kdim  ,4),vol(jdim,kdim,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
c
      wx = xrotrate
      wy = yrotrate
      wz = zrotrate
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 1000 j=1,jdim1
        do 1000 k=1,kdim1
          do 1000 i=1,idim1
 
c    compute cell centers 
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
c     increase velocity with rigid body rotaion component:       
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      q(j,k,i,2) = q(j,k,i,2) + ( cy * wz - cz * wy )
      q(j,k,i,3) = q(j,k,i,3) + ( cz * wx - cx * wz )
      q(j,k,i,4) = q(j,k,i,4) + ( cx * wy - cy * wx )
 
 1000 continue
c
      return
      end
