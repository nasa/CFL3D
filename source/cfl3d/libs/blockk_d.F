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
      subroutine blockk_d (q,qk0,idimr,jdimr,kdimr,idimt,jdimt,limblk,  
     .                     isva,it,ir,nvals,ldim,bck,iedge)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Transfer information from block (ir) to qk0 array of 
c      block (it).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(nvals,ldim,2),qk0(jdimt,idimt-1,ldim,2)
      dimension limblk(2,6),isva(2,2)
      dimension bck(jdimt,idimt-1,2)
c
      ist = limblk(it,1)
      iet = limblk(it,4)
      if (ist .eq. iet) then
         iinct = 1
      else
         iinct = (iet-ist)/abs(iet-ist)
      end if
c
      jst = limblk(it,2)
      jet = limblk(it,5)
      if (jst .eq. jet) then
         jinct = 1
      else
         jinct = (jet-jst)/abs(jet-jst)
      end if
c
c     determine the side of the q array to transfer from
c
c     k = constant side
c
      if (isva(ir,1)+isva(ir,2) .eq. 3) then
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i     and     j varies with j
c
            ij = 0
            do 200 i=ist,iet,iinct
               do 100 j=jst,jet,jinct
                  ij = ij + 1
                  do 50 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
   50             continue
  100          continue
  200       continue
         else
c
c     j varies with i     and     i varies with j
c
            ij = 0
            do 500 i=ist,iet,iinct
               do 400 j=jst,jet,jinct
                  ij = ij + 1
                  do 350 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
  350             continue
  400          continue
  500       continue
         end if
c  
c     j = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 4) then
c
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i    and    k varies with j
c
            ij = 0
            do 800 i=ist,iet,iinct
               do 700 j=jst,jet,jinct
                  ij = ij + 1
                  do 650 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
  650             continue
  700          continue
  800       continue
         else
c
c     k varies with i    and    i varies with j
c
            ij = 0
            do 1100 i=ist,iet,iinct
               do 1000 j=jst,jet,jinct
                  ij = ij + 1
                  do 950 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
  950             continue
 1000          continue
 1100       continue
         end if
c 
c     i = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 5) then
c
      if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .    (isva(ir,2) .eq. isva(it,2))) then
c
c     k varies with i    and    j varies with j
c
            ij = 0
            do 1400 i=ist,iet,iinct
               do 1300 j=jst,jet,jinct
                  ij = ij + 1
                  do 1250 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
 1250             continue
 1300          continue
 1400       continue
         else
c
c     j varies with i    and    k varies with j
c
            ij = 0
            do 1700 i=ist,iet,iinct
               do 1600 j=jst,jet,jinct
                  ij = ij + 1
                  do 1550 l = 1,ldim
                     qk0(j,i,l,1) = q(ij,l,1)
                     qk0(j,i,l,2) = q(ij,l,2)
                     bck(j,i,iedge) = 0.0
 1550             continue
 1600          continue
 1700       continue
         end if
c
      end if
c
      return
      end
