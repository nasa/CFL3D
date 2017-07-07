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
      subroutine rotateq(jd,kd,id,q,qrot,ista,iend,jsta,jend,ksta,kend,
     .                   dthtx,dthty,dthtz)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Rotate solution contained in array q through a specified
c     angle; rotated solution stored in qrot
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jd,kd,id,5),qrot(jd,kd,id,5)
c
      if (abs(real(dthtx)) .gt. 0.) then
c
c        rotate q about an axis parallel to the x-axis
c
         ca = cos(dthtx)
         sa = sin(dthtx)
c
         do 110 i=ista,iend
         do 120 j=jsta,jend
         do 130 k=ksta,kend
         qrot(j,k,i,1) =  q(j,k,i,1)
         qrot(j,k,i,2) =  q(j,k,i,2)
         temp          =  q(j,k,i,3)
         qrot(j,k,i,3) =  q(j,k,i,3)*ca - q(j,k,i,4)*sa
         qrot(j,k,i,4) =  temp      *sa + q(j,k,i,4)*ca
         qrot(j,k,i,5) =  q(j,k,i,5)
 130     continue
 120     continue
 110     continue
c
      else if (abs(real(dthty)) .gt. 0.) then
c
c        rotate q about an axis parallel to the y-axis
c
         ca = cos(dthty)
         sa = sin(dthty)
c
         do 210 i=ista,iend
         do 220 j=jsta,jend
         do 230 k=ksta,kend
         qrot(j,k,i,1) =  q(j,k,i,1)
         temp          =  q(j,k,i,2)
         qrot(j,k,i,2) =  q(j,k,i,2)*ca + q(j,k,i,4)*sa
         qrot(j,k,i,3) =  q(j,k,i,3)
         qrot(j,k,i,4) = -temp      *sa + q(j,k,i,4)*ca
         qrot(j,k,i,5) =  q(j,k,i,5)
 230     continue
 220     continue
 210     continue
c
      else if (abs(real(dthtz)) .gt. 0.) then
c
c        rotate q about an axis parallel to the z-axis
c
         ca = cos(dthtz)
         sa = sin(dthtz)
c
         do 310 i=ista,iend
         do 320 j=jsta,jend
         do 330 k=ksta,kend
         qrot(j,k,i,1) =  q(j,k,i,1)
         temp          =  q(j,k,i,2)
         qrot(j,k,i,2) =  q(j,k,i,2)*ca - q(j,k,i,3)*sa
         qrot(j,k,i,3) =  temp      *sa + q(j,k,i,3)*ca
         qrot(j,k,i,4) =  q(j,k,i,4)
         qrot(j,k,i,5) =  q(j,k,i,5)
 330     continue
 320     continue
 310     continue
c
      else
c
c     no rotation
c
         do 410 i=ista,iend
         do 420 j=jsta,jend
         do 430 k=ksta,kend
         qrot(j,k,i,1) =  q(j,k,i,1)
         qrot(j,k,i,2) =  q(j,k,i,2)
         qrot(j,k,i,3) =  q(j,k,i,3)
         qrot(j,k,i,4) =  q(j,k,i,4)
         qrot(j,k,i,5) =  q(j,k,i,5)
 430     continue
 420     continue
 410     continue
c
      end if
c
      return
      end
