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
      subroutine cblkj_d(nbli,idimt,jdimt,kdimt,
     .                   limblk,isva,it,ir,iedge,xyzr,nvals,xt,yt,zt,
     .                   ntime,lcnt,geom_miss,mxbli)
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
      dimension limblk(2,6),isva(2,2)
      dimension xyzr(nvals,3),xt(jdimt,kdimt,idimt),
     .          yt(jdimt,kdimt,idimt),
     .          zt(jdimt,kdimt,idimt)
      dimension geom_miss(2*mxbli)
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
      jt  = 1
      if (iedge.eq.2) jt = jdimt
c
c     determine the side of the q array to transfer from
c
c
c     k = constant side
c
      if (isva(ir,1)+isva(ir,2) .eq. 3) then
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i     and     j varies with k
c
            ij = 0
            do 200 i=ist,iet,iinct
               do 100 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k   ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i)  + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k   ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i)  + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k   ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i)  + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
  100          continue
  200       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         else
c
            ij = 0
            do 500 i=ist,iet,iinct
               do 400 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k  ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i) + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k  ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i) + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k  ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i) + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
  400          continue
  500       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         end if
c  
c     j = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 4) then
          if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .        (isva(ir,2) .eq. isva(it,2))) then
c
c     i varies with i    and    k varies with k
c
            ij = 0
            do 800 i=ist,iet,iinct
               do 700 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k  ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i) + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k  ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i) + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k  ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i) + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
  700          continue
  800       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         else
c
c     k varies with i    and    i varies with k
c
            ij = 0
            do 1100 i=ist,iet,iinct
               do 1000 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k  ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i) + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k  ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i) + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k  ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i) + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
 1000          continue
 1100       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         end if
c 
c     i = constant side
c
      else if (isva(ir,1)+isva(ir,2) .eq. 5) then
         if ((isva(ir,1) .eq. isva(it,1)) .or. 
     .       (isva(ir,2) .eq. isva(it,2))) then
c
c     k varies with k    and    j varies with i
c
            ij = 0
            do 1400 i=ist,iet,iinct
               do 1300 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k  ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i) + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k  ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i) + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k  ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i) + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
 1300          continue
 1400       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         else
c
c     k varies with i    and    j varies with k
c
            ij = 0
      do 1700 i=ist,iet,iinct
      do 1600 k=kst,ket,kinct
                  ij = ij + 1
                  xr1 = xyzr(ij,1)
                  yr1 = xyzr(ij,2)
                  zr1 = xyzr(ij,3)
                  xt1    = .25*( xt(jt,k  ,i) + xt(jt,k  ,i+1)
     .                        +  xt(jt,k+1,i) + xt(jt,k+1,i+1) )
                  yt1    = .25*( yt(jt,k  ,i) + yt(jt,k  ,i+1)
     .                        +  yt(jt,k+1,i) + yt(jt,k+1,i+1) )
                  zt1    = .25*( zt(jt,k  ,i) + zt(jt,k  ,i+1)
     .                        +  zt(jt,k+1,i) + zt(jt,k+1,i+1) )
                  eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+
     &                      (zr1-zt1)**2))
 1600          continue
 1700       continue
            if (ntime.eq.1) then
               geom_miss(lcnt) = eps
            else
               geom_miss(lcnt) = 0.
            end if
c
         end if
      end if
      return
      end
