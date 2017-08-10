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
      subroutine pre_blockj (idimr,jdimr,kdimr,limblk,isva,it,ir,
     .                       q,qlocal,maxdims,ldim,ij,ivolflg,
     .                       myid,mblk2nd,maxbl,bou,ibufdim,nbuf,nou)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Check information transferred from block (ir) to 
c      qj0 array of block (it).
c***********************************************************************
c
c     ivolflg...a flag to indicate if the "q" array being passed to
c               this routine is the cell volume array (ivolflg=1)
c               or q/vist3d/turres (ivolflg=0). this is needed
c               because the volume array contains one less i-plane
c               than the other arrays
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension mblk2nd(maxbl)
      dimension limblk(2,6),isva(2,2)
      dimension q(jdimr,kdimr,idimr,ldim)
      dimension qlocal(maxdims,ldim,2)
c
      ist = limblk(it,1)
      iet = limblk(it,4)
      if (ist .eq. iet) then
         iinct = 1
      else
         iinct = (iet-ist)/abs(iet-ist)
      end if
c
      kst = limblk(it,3)
      ket = limblk(it,6)
      if (kst .eq. ket) then
         kinct = 1
      else
         kinct = (ket-kst)/abs(ket-kst)
      end if
c
      eps = 0.
c
      isr = limblk(ir,1)
      ier = limblk(ir,4)
      jsr = limblk(ir,2)
      jer = limblk(ir,5)
      ksr = limblk(ir,3)
      ker = limblk(ir,6)
c 
c     determine the side of the q array to transfer from
c
c     k = constant side
c
      if (isva(ir,1)+isva(ir,2) .eq. 3) then
         if (ksr.eq.1) then
            kloc1r = 1
            kloc2r = 2
         else
            kloc1r = kdimr-1
            kloc2r = kdimr-2
         end if
c
         if (kdimr.eq.2) then
            kloc1r = 1
            kloc2r = 1
         end if
c
         if (jer .eq. jsr) then
            jincr = 1
         else
            jincr = (jer-jsr)/abs(jer-jsr)
         end if
c
         if (ier .eq. isr) then
            iincr = 1
         else
            iincr = (ier-isr)/abs(ier-isr)
         end if
c
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i     and     j varies with k
c
            ij = 0
            icount = -1
            do 200 i=ist,iet,iinct
               icount = icount + 1
               jcount = -1
               do 100 k=kst,ket,kinct
                  jcount = jcount + 1
                  ilocr  = isr + iincr*icount
                  jlocr  = jsr + jincr*jcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 50 l = 1, ldim
                     qlocal(ij,l,1) = q(jlocr,kloc1r,ilocr,l)
                     qlocal(ij,l,2) = q(jlocr,kloc2r,ilocr,l)
   50             continue
  100          continue
  200       continue
         else
c
            ij = 0
            jcount = -1
            do 500 i=ist,iet,iinct
               jcount = jcount + 1
               icount = -1
               do 400 k=kst,ket,kinct
                  icount = icount + 1
                  ilocr  = isr + iincr*icount
                  jlocr  = jsr + jincr*jcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 350 l = 1, ldim
                     qlocal(ij,l,1) = q(jlocr,kloc1r,ilocr,l)
                     qlocal(ij,l,2) = q(jlocr,kloc2r,ilocr,l)
  350             continue
  400          continue
  500       continue
         end if
c  
c     j = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 4) then
         if (jsr.eq.1) then
            jloc1r = 1
            jloc2r = 2
         else
            jloc1r = jdimr-1
            jloc2r = jdimr-2
         end if
c
         if (jdimr.eq.2) then
            jloc1r = 1
            jloc2r = 1
         end if
c
         if (ier .eq. isr) then
            iincr = 1
         else
            iincr = (ier-isr)/abs(ier-isr)
         end if
c
         if (ker .eq. ksr) then
            kincr = 1
         else
            kincr = (ker-ksr)/abs(ker-ksr)
         end if
c
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i    and    k varies with k
c
            ij = 0
            icount = -1
            do 800 i=ist,iet,iinct
               icount = icount + 1
               kcount = -1
               do 700 k=kst,ket,kinct
                  kcount = kcount + 1
                  ilocr  = isr + iincr*icount
                  klocr  = ksr + kincr*kcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 650 l = 1, ldim
                     qlocal(ij,l,1) = q(jloc1r,klocr,ilocr,l)
                     qlocal(ij,l,2) = q(jloc2r,klocr,ilocr,l)
  650             continue
  700          continue
  800       continue
         else
c
c     k varies with i    and    i varies with k
c
            ij = 0
            kcount = -1
            do 1100 i=ist,iet,iinct
               kcount = kcount + 1
               icount = -1
               do 1000 k=kst,ket,kinct
                  icount = icount + 1
                  ilocr  = isr + iincr*icount
                  klocr  = ksr + kincr*kcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 950 l = 1, ldim
                     qlocal(ij,l,1) = q(jloc1r,klocr,ilocr,l)
                     qlocal(ij,l,2) = q(jloc2r,klocr,ilocr,l)
  950             continue
 1000          continue
 1100       continue
         end if
c 
c     i = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 5) then
         if (isr.eq.1) then
            iloc1r = 1
            iloc2r = 2
         else
            if (ivolflg.eq.0) then
               iloc1r = idimr-1
               iloc2r = idimr-2
            else
               iloc1r = idimr
               iloc2r = idimr-1
            end if
         end if
c
         if (idimr.eq.2) then
            iloc1r = 1
            iloc2r = 1
         end if
c
         if (jer .eq. jsr) then
            jincr = 1
         else
            jincr = (jer-jsr)/abs(jer-jsr)
         end if
c
         if (ker .eq. ksr) then
            kincr = 1
         else
            kincr = (ker-ksr)/abs(ker-ksr)
         end if
c
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     k varies with k    and    j varies with i
c
            ij = 0
            jcount = -1
            do 1400 i=ist,iet,iinct
               jcount = jcount + 1
               kcount = -1
               do 1300 k=kst,ket,kinct
                  kcount = kcount + 1
                  jlocr  = jsr + jincr*jcount
                  klocr  = ksr + kincr*kcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 1250 l = 1, ldim
                     qlocal(ij,l,1) = q(jlocr,klocr,iloc1r,l)
                     qlocal(ij,l,2) = q(jlocr,klocr,iloc2r,l)
 1250             continue
 1300          continue
 1400       continue
         else
c
c     k varies with i    and    j varies with k
c
            ij = 0
            kcount = -1
            do 1700 i=ist,iet,iinct
               kcount = kcount + 1
               jcount = -1
               do 1600 k=kst,ket,kinct
                  jcount = jcount + 1
                  jlocr  = jsr + jincr*jcount
                  klocr  = ksr + kincr*kcount
                  ij = ij + 1
                  if (ij.gt.maxdims) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in pre_blockj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do 1550 l = 1, ldim
                     qlocal(ij,l,1) = q(jlocr,klocr,iloc1r,l)
                     qlocal(ij,l,2) = q(jlocr,klocr,iloc2r,l)
 1550             continue
 1600          continue
 1700       continue
         end if
      end if
c
      return
      end
