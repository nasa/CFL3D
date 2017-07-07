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
      subroutine pre_patch (nbl,lw,icount,ninter,iindex,intmax,nsub1,
     .                      isav_pat,isav_pat_b,jjmax1,kkmax1,
     .                      iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                      ierrflg)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set up data arrays for patch bc message passing
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension lw(65,maxbl)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl)
      dimension isav_pat(intmax,17), isav_pat_b(intmax,nsub1,6)
      dimension jjmax1(nsub1),kkmax1(nsub1)
      dimension iiint1(nsub1),iiint2(nsub1)
      dimension iindex(intmax,6*nsub1+9)
c
      if(abs(ninter).gt.0) then
c
c     loop over all patch interfaces
c
      do 8001 icheck=1,abs(ninter)
      lmax1 =  iindex(icheck,1)
      nblcc =  iindex(icheck,lmax1+2)
      if (nblcc.ne.nbl) go to 8001
      lst   =  iindex(icheck,2*lmax1+5)
      npt   =  iindex(icheck,2*lmax1+4)
c
c     set range of points requiring interpolation on "to" side of
c     patch interface
c
      j21       = iindex(icheck,2*lmax1+6)
      j22       = iindex(icheck,2*lmax1+7)
      k21       = iindex(icheck,2*lmax1+8)
      k22       = iindex(icheck,2*lmax1+9)
c
c     patch surface to be interpolated to is an i=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.1) then
          lqedge    = lw(4,nbl)
          lqedgt    = lw(25,nbl)
          lqedgv    = lw(30,nbl)
          lqedgb    = lw(33,nbl)
          jmax2     = jdimg(nbl)
          kmax2     = kdimg(nbl)
          if (iindex(icheck,2*lmax1+3).eq.11) then
             mint1 = 1
             mint2 = 2
          else
             mint1 = 3
             mint2 = 4
          end if
      end if
c
c     patch surface to be interpolated to is a j=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.2) then
         lqedge    = lw(2,nbl)
         lqedgt    = lw(23,nbl)
         lqedgv    = lw(28,nbl)
         lqedgb    = lw(31,nbl)
         jmax2     = kdimg(nbl)
         kmax2     = idimg(nbl)-1
         if (iindex(icheck,2*lmax1+3).eq.21) then
            mint1 = 1
            mint2 = 2
         else
            mint1 = 3
            mint2 = 4
         end if
      end if
c
c     patch surface to be interpolated to is a k=constant surface
c
      if (iindex(icheck,2*lmax1+3)/10.eq.3) then
         lqedge    = lw(3,nbl)
         lqedgt    = lw(24,nbl)
         lqedgv    = lw(29,nbl)
         lqedgb    = lw(32,nbl)
         jmax2     = jdimg(nbl)
         kmax2     = idimg(nbl)-1
         if (iindex(icheck,2*lmax1+3).eq.31) then
            mint1 = 1
            mint2 = 2
         else
            mint1 = 3
            mint2 = 4
         end if
      end if
c
c     put some relevant values into isav_pat array
c
      icount = icount + 1
      isav_pat(icount,1)  = nbl
      isav_pat(icount,2)  = lmax1
      isav_pat(icount,3)  = j21
      isav_pat(icount,4)  = j22  
      isav_pat(icount,5)  = k21  
      isav_pat(icount,6)  = k22  
      isav_pat(icount,7)  = lqedge
      isav_pat(icount,8)  = lqedgt
      isav_pat(icount,9)  = lqedgv
      isav_pat(icount,10) = lqedgb
      isav_pat(icount,11) = jmax2 
      isav_pat(icount,12) = kmax2 
      isav_pat(icount,13) = mint1 
      isav_pat(icount,14) = mint2 
      isav_pat(icount,15) = lst
      isav_pat(icount,16) = npt
      isav_pat(icount,17) = icheck
c
      do 1705 l=1,lmax1
      mbl = iindex(icheck,l+1)
      mtype  = iindex(icheck,l+lmax1+2)
c
c     patch surface to be interpolated from is an i=constant surface
c
      if (mtype/10.eq.1) then
         jjmax1(l) = jdimg(mbl)
         kkmax1(l) = kdimg(mbl)
         if (mtype.eq.11) then
            iiint1(l) = 1
            iiint2(l) = min(2,idimg(mbl)-1)
         else
            iiint1(l) = idimg(mbl)-1
            iiint2(l) = max(1,idimg(mbl)-2)
         end if
      end if
c
c     patch surface to be interpolated from is a j=constant surface
c
      if (mtype/10.eq.2) then
         jjmax1(l) = kdimg(mbl)
         kkmax1(l) = idimg(mbl)
         if (mtype.eq.21) then
            iiint1(l) = 1
            iiint2(l) = min(2,jdimg(mbl)-1)
         else
            iiint1(l) = jdimg(mbl)-1
            iiint2(l) = max(1,jdimg(mbl)-2)
         end if
      end if

c
c     patch surface to be interpolated from is a k=constant surface
c
      if (mtype/10.eq.3) then
         jjmax1(l) = jdimg(mbl)
         kkmax1(l) = idimg(mbl)
         if (mtype.eq.31) then
            iiint1(l) = 1
            iiint2(l) = min(2,kdimg(mbl)-1)
         else
            iiint1(l) = kdimg(mbl)-1
            iiint2(l) = max(1,kdimg(mbl)-2)
         end if
      end if
c
      isav_pat_b(icount,l,1) = mbl
      isav_pat_b(icount,l,2) = mtype/10
      isav_pat_b(icount,l,3) = jjmax1(l)
      isav_pat_b(icount,l,4) = kkmax1(l)
      isav_pat_b(icount,l,5) = iiint1(l)
      isav_pat_b(icount,l,6) = iiint2(l)
c
 1705 continue
c
 8001 continue
c
      end if
c
      return
      end
