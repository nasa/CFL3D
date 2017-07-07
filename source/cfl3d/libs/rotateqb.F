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
      subroutine rotateqb(nbl,dthtx,dthty,dthtz,maxbl,iitot,ibcg,
     .                    lig,lbg,ibpntsg,iipntsg,qb)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Rotate solution in qb array through angle dthtx/y/z for 
c     chimera scheme with rotating grids
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ibcg(iitot),lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),
     .          iipntsg(maxbl),qb(iitot,5,3)
c
      iset=1
c
c     lsta = lbg(nbl)
c     lend = lsta-1
      lsta = lig(nbl)
      lend = lsta+iipntsg(nbl)-1
c
      if (abs(real(dthtx)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the x-axis
c
         ca = cos(dthtx)
         sa = sin(dthtx)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 10 l=lsta,lend
c           qb3 = qb(ibcg(l),3,iset)
c           qb(ibcg(l),3,iset) = qb3*ca - qb(ibcg(l),4,iset)*sa
c           qb(ibcg(l),4,iset) = qb3*sa + qb(ibcg(l),4,iset)*ca
            qb3 = qb(l,3,iset)
            qb(l,3,iset) = qb3*ca - qb(l,4,iset)*sa
            qb(l,4,iset) = qb3*sa + qb(l,4,iset)*ca
   10       continue
c        end if
c
      end if
c
      if (abs(real(dthty)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the y-axis
c
         ca = cos(dthty)
         sa = sin(dthty)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 20 l=lsta,lend
c           qb2 = qb(ibcg(l),2,iset)
c           qb(ibcg(l),2,iset) =  qb2*ca + qb(ibcg(l),4,iset)*sa
c           qb(ibcg(l),4,iset) = -qb2*sa + qb(ibcg(l),4,iset)*ca
            qb2 = qb(l,2,iset)
            qb(l,2,iset) =  qb2*ca + qb(l,4,iset)*sa
            qb(l,4,iset) = -qb2*sa + qb(l,4,iset)*ca
   20       continue
c        end if
c
      end if
c
      if (abs(real(dthtz)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the z-axis
c
         ca = cos(dthtz)
         sa = sin(dthtz)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 30 l=lsta,lend
c           qb2 = qb(ibcg(l),2,iset)
c           qb(ibcg(l),2,iset) = qb2*ca - qb(ibcg(l),3,iset)*sa
c           qb(ibcg(l),3,iset) = qb2*sa + qb(ibcg(l),3,iset)*ca
            qb2 = qb(l,2,iset)
            qb(l,2,iset) = qb2*ca - qb(l,3,iset)*sa
            qb(l,3,iset) = qb2*sa + qb(l,3,iset)*ca
   30       continue
c        end if
c
      end if
c
      return
      end
