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
      subroutine direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .                   a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute (normalized) directed area components, 
c     or equivalently, components of unit normal to cell face
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      xxie = x8-x7
      yxie = y8-y7
      zxie = z8-z7
      xeta = x6-x5
      yeta = y6-y5
      zeta = z6-z5
      a1 = yxie*zeta-zxie*yeta
      a2 = zxie*xeta-xxie*zeta
      a3 = xxie*yeta-yxie*xeta
      d  = sqrt(a1*a1+a2*a2+a3*a3)
      if(real(d) .le. 0.) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) ' WARNINING: a cell with zero area has ',
     .             'been detected in subroutine direct'
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) '  - severe problem'
        d=1.
      end if
      a1 = a1/d
      a2 = a2/d
      a3 = a3/d
c
c     find coordinate direction with maximum area component
c
      amax  = ccabs(a1)
      imaxa = 1
      if (abs(real(a2)).gt.real(amax)) then
         amax  = ccabs(a2)
         imaxa = 2
      end if
      if (abs(real(a3)).gt.real(amax)) then
         amax  = ccabs(a3)
         imaxa = 3
      end if
      return
      end      
