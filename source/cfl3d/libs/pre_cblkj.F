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
      subroutine pre_cblkj (idimr,jdimr,kdimr,limblk,isva,it,ir,
     .                      xr,yr,zr,xyzn,maxdims,ij,myid,ibufdim,
     .                      nbuf,bou,nou)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Check information transferred from block (ir) to 
c      qj0 array of block (it).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension limblk(2,6),isva(2,2)
      dimension xr(jdimr,kdimr,idimr),
     .          yr(jdimr,kdimr,idimr),
     .          zr(jdimr,kdimr,idimr)
      dimension xyzn(maxdims,3)
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
         kr = kloc1r
         if (kr.gt.1) kr = kdimr
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jlocr  ,kr,ilocr) +
     .                               xr(jlocr  ,kr,ilocr+1) +
     .                               xr(jlocr+1,kr,ilocr) +
     .                               xr(jlocr+1,kr,ilocr+1) )
                  xyzn(ij,2) = .25*( yr(jlocr  ,kr,ilocr) +
     .                               yr(jlocr  ,kr,ilocr+1) +
     .                               yr(jlocr+1,kr,ilocr) +
     .                               yr(jlocr+1,kr,ilocr+1) )
                  xyzn(ij,3) = .25*( zr(jlocr  ,kr,ilocr) +
     .                               zr(jlocr  ,kr,ilocr+1) +
     .                               zr(jlocr+1,kr,ilocr) +
     .                               zr(jlocr+1,kr,ilocr+1) )
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jlocr  ,kr,ilocr) +
     .                               xr(jlocr  ,kr,ilocr+1) +
     .                               xr(jlocr+1,kr,ilocr) +
     .                               xr(jlocr+1,kr,ilocr+1) )
                  xyzn(ij,2) = .25*( yr(jlocr  ,kr,ilocr) +
     .                               yr(jlocr  ,kr,ilocr+1) +
     .                               yr(jlocr+1,kr,ilocr) +
     .                               yr(jlocr+1,kr,ilocr+1) )
                  xyzn(ij,3) = .25*( zr(jlocr  ,kr,ilocr) +
     .                               zr(jlocr  ,kr,ilocr+1) +
     .                               zr(jlocr+1,kr,ilocr) +
     .                               zr(jlocr+1,kr,ilocr+1) )
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
         jr = jloc1r
         if (jr.gt.1) jr = jdimr
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jr,klocr  ,ilocr) +
     .                               xr(jr,klocr  ,ilocr+1) +
     .                               xr(jr,klocr+1,ilocr) +
     .                               xr(jr,klocr+1,ilocr+1) )
                  xyzn(ij,2) = .25*( yr(jr,klocr  ,ilocr) +
     .                               yr(jr,klocr  ,ilocr+1) +
     .                               yr(jr,klocr+1,ilocr) +
     .                               yr(jr,klocr+1,ilocr+1) )
                  xyzn(ij,3) = .25*( zr(jr,klocr  ,ilocr) +
     .                               zr(jr,klocr  ,ilocr+1) +
     .                               zr(jr,klocr+1,ilocr) +
     .                               zr(jr,klocr+1,ilocr+1) )
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jr,klocr  ,ilocr) +
     .                               xr(jr,klocr  ,ilocr+1) +
     .                               xr(jr,klocr+1,ilocr) +
     .                               xr(jr,klocr+1,ilocr+1) )
                  xyzn(ij,2) = .25*( yr(jr,klocr  ,ilocr) +
     .                               yr(jr,klocr  ,ilocr+1) +
     .                               yr(jr,klocr+1,ilocr) +
     .                               yr(jr,klocr+1,ilocr+1) )
                  xyzn(ij,3) = .25*( zr(jr,klocr  ,ilocr) +
     .                               zr(jr,klocr  ,ilocr+1) +
     .                               zr(jr,klocr+1,ilocr) +
     .                               zr(jr,klocr+1,ilocr+1) )
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
            iloc1r = idimr-1
            iloc2r = idimr-2
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
         irn = iloc1r
         if (irn.gt.1) irn = idimr
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jlocr  ,klocr,irn) +
     .                               xr(jlocr  ,klocr+1,irn) +
     .                               xr(jlocr+1,klocr,irn) +
     .                               xr(jlocr+1,klocr+1,irn) )
                  xyzn(ij,2) = .25*( yr(jlocr  ,klocr,irn) +
     .                               yr(jlocr  ,klocr+1,irn) +
     .                               yr(jlocr+1,klocr,irn) +
     .                               yr(jlocr+1,klocr+1,irn) )
                  xyzn(ij,3) = .25*( zr(jlocr  ,klocr,irn) +
     .                               zr(jlocr  ,klocr+1,irn) +
     .                               zr(jlocr+1,klocr,irn) +
     .                               zr(jlocr+1,klocr+1,irn) )
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
                     write(bou(nou(1),1),*)' stopping in pre_cblkj',
     .               ' ij.gt.maxdims'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  xyzn(ij,1) = .25*( xr(jlocr  ,klocr,irn) +
     .                               xr(jlocr  ,klocr+1,irn) +
     .                               xr(jlocr+1,klocr,irn) +
     .                               xr(jlocr+1,klocr+1,irn) )
                  xyzn(ij,2) = .25*( yr(jlocr  ,klocr,irn) +
     .                               yr(jlocr  ,klocr+1,irn) +
     .                               yr(jlocr+1,klocr,irn) +
     .                               yr(jlocr+1,klocr+1,irn) )
                  xyzn(ij,3) = .25*( zr(jlocr  ,klocr,irn) +
     .                               zr(jlocr  ,klocr+1,irn) +
     .                               zr(jlocr+1,klocr,irn) +
     .                               zr(jlocr+1,klocr+1,irn) )
 1600          continue
 1700       continue
         end if
      end if
c
      return
      end
