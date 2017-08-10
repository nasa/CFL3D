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
      subroutine transmc(itran,rfreqt,utran,vtran,wtran,xorg,yorg,
     .                   zorg,xorg0,yorg0,zorg0,xmc,ymc,zmc,iupdat,
     .                   time2)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determines increment to moment center due to translation
c
c     itran....modulation for translational motion
c              = 0 no translation
c              = 1 constant velocity
c              = 2 sinusoidal variation of displacement
c              = 3 smooth increase in displacement, asypmtotically
c                 reaching a fixed displacement
c
c     iupdat..flag to update moment center position
c             = 0 don't update position
c             > 0 update position
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
c     ft modulates the displacement
c
      if (itran.eq. 0)  then
         return
      else if (itran .eq. 1)  then
         ft     = time2
      else if (itran .eq. 2)  then
         ft     = sin(rfreqt*time2)
      else if (itran .eq. 3)  then
         expt   = exp(-rfreqt*time2)
         ft     = 1.-expt
      end if
c
      xold = xorg
      yold = yorg
      zold = zorg
c
      xnew = utran*ft + xorg0
      ynew = vtran*ft + yorg0
      znew = wtran*ft + zorg0
c
      dx     = xnew - xold
      dy     = ynew - yold
      dz     = znew - zold
c
c*************************************************
c     update moment center to new position;
c     update rotation point for moment center
c*************************************************
c
      if (iupdat .gt. 0) then
c
         xmc = xmc+dx
         ymc = ymc+dy
         zmc = zmc+dz
c
         xorg = xorg+dx
         yorg = yorg+dy
         zorg = zorg+dz
c
      end if
c
      return
      end
