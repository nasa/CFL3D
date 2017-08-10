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
      subroutine delv(jdim,kdim,idim,q,sj,sk,si,vol,ux,wt,blank,iover,
     .                qj0,qk0,qi0,bcj,bck,bci,nbl,volj0,volk0,voli0,
     .                maxbl,vormax,ivmax,jvmax,kvmax)
c
c     $Id$
c
c**********************************************************************
c     Purpose:  Evaluate velocity derivatives at cell centers
c**********************************************************************
c
c      input arrays   : 
c       q             : primitive variables at cell centers
c                       (rho,u,v,w,p)
c       vol           : cell volumes
c       sj,sk,si      : metrics
c                       (direction cosines,areas,speeds of cell faces)
c       ifj,ifk,ifi   : =0 don't compute contribution from that direction
c                       >0 do      "
c      output arrays  :
c       ux            : velocity derivatives
c                           1 = ux    4 = vx    7 = wx
c                           2 = uy    5 = vy    8 = wy
c                           3 = uz    6 = vz    9 = wz
c      scratch arrays :
c       wt(1-3)       : cross-flow temporaries
c
c**********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),ux(jdim-1,kdim-1,idim-1,9)
      dimension sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension vol(jdim,kdim,idim-1),blank(jdim,kdim,idim)
      dimension wt(jdim,kdim,9)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension volj0(kdim,idim-1,4),
     .          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
c
      common /twod/ i2d
c
      ifi=1
      ifj=1
      ifk=1
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
      do 999 l=1,9
      do 999 i=1,idim1
      do 999 k=1,kdim1
      do 999 j=1,jdim1
        ux(j,k,i,l) = 0.
 999  continue
      onec = 1.
c
c      J-direction contributions
c
      if(ifj .gt. 0) then
      do 1000 i=1,idim1
      do 1000 k=1,kdim1
c
c      cycle through interfaces
c
      do 100 j=2,jdim1
      term      = sj(j,k,i,4)/(vol(j,k,i)+vol(j-1,k,i))
      wt(j,1,1) = term*( (q(j,k,i,2)-q(j-1,k,i,2))*sj(j,k,i,1) )
      wt(j,1,2) = term*( (q(j,k,i,2)-q(j-1,k,i,2))*sj(j,k,i,2) )
      wt(j,1,3) = term*( (q(j,k,i,2)-q(j-1,k,i,2))*sj(j,k,i,3) )
      wt(j,1,4) = term*( (q(j,k,i,3)-q(j-1,k,i,3))*sj(j,k,i,1) )
      wt(j,1,5) = term*( (q(j,k,i,3)-q(j-1,k,i,3))*sj(j,k,i,2) )
      wt(j,1,6) = term*( (q(j,k,i,3)-q(j-1,k,i,3))*sj(j,k,i,3) )
      wt(j,1,7) = term*( (q(j,k,i,4)-q(j-1,k,i,4))*sj(j,k,i,1) )
      wt(j,1,8) = term*( (q(j,k,i,4)-q(j-1,k,i,4))*sj(j,k,i,2) )
      wt(j,1,9) = term*( (q(j,k,i,4)-q(j-1,k,i,4))*sj(j,k,i,3) )
  100 continue
      j=1
      jp1=min(jdim-1,j+1)
      term      = sj(j,k,i,4)/(volj0(k,i,1)+vol(j,k,i))
      factor    = bcj(k,i,1) + 1.0
      term = term*factor
      wt(j,1,1) = term*( (q(j,k,i,2)-qj0(k,i,2,1))*sj(j,k,i,1) )
      wt(j,1,2) = term*( (q(j,k,i,2)-qj0(k,i,2,1))*sj(j,k,i,2) )
      wt(j,1,3) = term*( (q(j,k,i,2)-qj0(k,i,2,1))*sj(j,k,i,3) )
      wt(j,1,4) = term*( (q(j,k,i,3)-qj0(k,i,3,1))*sj(j,k,i,1) )
      wt(j,1,5) = term*( (q(j,k,i,3)-qj0(k,i,3,1))*sj(j,k,i,2) )
      wt(j,1,6) = term*( (q(j,k,i,3)-qj0(k,i,3,1))*sj(j,k,i,3) )
      wt(j,1,7) = term*( (q(j,k,i,4)-qj0(k,i,4,1))*sj(j,k,i,1) )
      wt(j,1,8) = term*( (q(j,k,i,4)-qj0(k,i,4,1))*sj(j,k,i,2) )
      wt(j,1,9) = term*( (q(j,k,i,4)-qj0(k,i,4,1))*sj(j,k,i,3) )
      j=jdim
      jm2=max(1,j-2)
      term      = sj(j,k,i,4)/(volj0(k,i,3)+vol(j-1,k,i))
      factor    = bcj(k,i,2) + 1.0
      term = term*factor
      wt(j,1,1) = term*( (qj0(k,i,2,3)-q(j-1,k,i,2))*sj(j,k,i,1) )
      wt(j,1,2) = term*( (qj0(k,i,2,3)-q(j-1,k,i,2))*sj(j,k,i,2) )
      wt(j,1,3) = term*( (qj0(k,i,2,3)-q(j-1,k,i,2))*sj(j,k,i,3) )
      wt(j,1,4) = term*( (qj0(k,i,3,3)-q(j-1,k,i,3))*sj(j,k,i,1) )
      wt(j,1,5) = term*( (qj0(k,i,3,3)-q(j-1,k,i,3))*sj(j,k,i,2) )
      wt(j,1,6) = term*( (qj0(k,i,3,3)-q(j-1,k,i,3))*sj(j,k,i,3) )
      wt(j,1,7) = term*( (qj0(k,i,4,3)-q(j-1,k,i,4))*sj(j,k,i,1) )
      wt(j,1,8) = term*( (qj0(k,i,4,3)-q(j-1,k,i,4))*sj(j,k,i,2) )
      wt(j,1,9) = term*( (qj0(k,i,4,3)-q(j-1,k,i,4))*sj(j,k,i,3) )
c
c      cycle through cell centers
c
      do 200 l=1,9
      do 200 j=1,jdim1
      ux(j,k,i,l)   = ux(j,k,i,l) + wt(j,1,l) + wt(j+1,1,l)
  200 continue
 1000 continue
      end if
c
c      K-direction contributions
c
      if(ifk .gt. 0) then
      do 2000 i=1,idim1
c
c      cycle through interfaces
c
      do 1200 k=2,kdim-1
      do 1200 j=1,jdim1
      term      = sk(j,k,i,4)/(vol(j,k,i)+vol(j,k-1,i))
      wt(j,k,1) = term*( (q(j,k,i,2)-q(j,k-1,i,2))*sk(j,k,i,1) )
      wt(j,k,2) = term*( (q(j,k,i,2)-q(j,k-1,i,2))*sk(j,k,i,2) )
      wt(j,k,3) = term*( (q(j,k,i,2)-q(j,k-1,i,2))*sk(j,k,i,3) )
      wt(j,k,4) = term*( (q(j,k,i,3)-q(j,k-1,i,3))*sk(j,k,i,1) )
      wt(j,k,5) = term*( (q(j,k,i,3)-q(j,k-1,i,3))*sk(j,k,i,2) )
      wt(j,k,6) = term*( (q(j,k,i,3)-q(j,k-1,i,3))*sk(j,k,i,3) )
      wt(j,k,7) = term*( (q(j,k,i,4)-q(j,k-1,i,4))*sk(j,k,i,1) )
      wt(j,k,8) = term*( (q(j,k,i,4)-q(j,k-1,i,4))*sk(j,k,i,2) )
      wt(j,k,9) = term*( (q(j,k,i,4)-q(j,k-1,i,4))*sk(j,k,i,3) )
 1200 continue
      k=1
      kp1=min(kdim-1,k+1)
      do 1201 j=1,jdim1
      term      = sk(j,k,i,4)/(volk0(j,i,1)+vol(j,k,i))
      factor    = bck(j,i,1) + 1.0
      term = term*factor
        wt(j,k,1) = term*( (q(j,k,i,2)-qk0(j,i,2,1))*sk(j,k,i,1) )
        wt(j,k,2) = term*( (q(j,k,i,2)-qk0(j,i,2,1))*sk(j,k,i,2) )
        wt(j,k,3) = term*( (q(j,k,i,2)-qk0(j,i,2,1))*sk(j,k,i,3) )
        wt(j,k,4) = term*( (q(j,k,i,3)-qk0(j,i,3,1))*sk(j,k,i,1) )
        wt(j,k,5) = term*( (q(j,k,i,3)-qk0(j,i,3,1))*sk(j,k,i,2) )
        wt(j,k,6) = term*( (q(j,k,i,3)-qk0(j,i,3,1))*sk(j,k,i,3) )
        wt(j,k,7) = term*( (q(j,k,i,4)-qk0(j,i,4,1))*sk(j,k,i,1) )
        wt(j,k,8) = term*( (q(j,k,i,4)-qk0(j,i,4,1))*sk(j,k,i,2) )
        wt(j,k,9) = term*( (q(j,k,i,4)-qk0(j,i,4,1))*sk(j,k,i,3) )
 1201 continue
      k=kdim
      km2=max(1,k-2)
      do 1202 j=1,jdim1
      term      = sk(j,k,i,4)/(volk0(j,i,3)+vol(j,k-1,i))
      factor    = bck(j,i,2) + 1.0
      term = term*factor
        wt(j,k,1) = term*( (qk0(j,i,2,3)-q(j,k-1,i,2))*sk(j,k,i,1) )
        wt(j,k,2) = term*( (qk0(j,i,2,3)-q(j,k-1,i,2))*sk(j,k,i,2) )
        wt(j,k,3) = term*( (qk0(j,i,2,3)-q(j,k-1,i,2))*sk(j,k,i,3) )
        wt(j,k,4) = term*( (qk0(j,i,3,3)-q(j,k-1,i,3))*sk(j,k,i,1) )
        wt(j,k,5) = term*( (qk0(j,i,3,3)-q(j,k-1,i,3))*sk(j,k,i,2) )
        wt(j,k,6) = term*( (qk0(j,i,3,3)-q(j,k-1,i,3))*sk(j,k,i,3) )
        wt(j,k,7) = term*( (qk0(j,i,4,3)-q(j,k-1,i,4))*sk(j,k,i,1) )
        wt(j,k,8) = term*( (qk0(j,i,4,3)-q(j,k-1,i,4))*sk(j,k,i,2) )
        wt(j,k,9) = term*( (qk0(j,i,4,3)-q(j,k-1,i,4))*sk(j,k,i,3) )
 1202 continue
c
c      cycle through cell centers
c
      do 1400 l=1,9
      do 1400 k=1,kdim1
      do 1400 j=1,jdim1
      ux(j,k,i,l) = ux(j,k,i,l) + wt(j,k,l) + wt(j,k+1,l)
 1400 continue
 2000 continue
      end if
c
c      I-direction contributions
c
      if(ifi .gt. 0) then
c      cycle through interfaces
c
      if(i2d .eq. 0) then
      do 1700 i=2,idim1
      do 1700 k=1,kdim1
      do 1700 j=1,jdim1
      term      = si(j,k,i,4)/(vol(j,k,i)+vol(j,k,i-1))
      wt(j,k,1) = term*( (q(j,k,i,2)-q(j,k,i-1,2))*si(j,k,i,1) )
      wt(j,k,2) = term*( (q(j,k,i,2)-q(j,k,i-1,2))*si(j,k,i,2) )
      wt(j,k,3) = term*( (q(j,k,i,2)-q(j,k,i-1,2))*si(j,k,i,3) )
      wt(j,k,4) = term*( (q(j,k,i,3)-q(j,k,i-1,3))*si(j,k,i,1) )
      wt(j,k,5) = term*( (q(j,k,i,3)-q(j,k,i-1,3))*si(j,k,i,2) )
      wt(j,k,6) = term*( (q(j,k,i,3)-q(j,k,i-1,3))*si(j,k,i,3) )
      wt(j,k,7) = term*( (q(j,k,i,4)-q(j,k,i-1,4))*si(j,k,i,1) )
      wt(j,k,8) = term*( (q(j,k,i,4)-q(j,k,i-1,4))*si(j,k,i,2) )
      wt(j,k,9) = term*( (q(j,k,i,4)-q(j,k,i-1,4))*si(j,k,i,3) )
c
      ux(j,k,i-1,1) = ux(j,k,i-1,1) + wt(j,k,1)
      ux(j,k,i-1,2) = ux(j,k,i-1,2) + wt(j,k,2)
      ux(j,k,i-1,3) = ux(j,k,i-1,3) + wt(j,k,3)
      ux(j,k,i-1,4) = ux(j,k,i-1,4) + wt(j,k,4)
      ux(j,k,i-1,5) = ux(j,k,i-1,5) + wt(j,k,5)
      ux(j,k,i-1,6) = ux(j,k,i-1,6) + wt(j,k,6)
      ux(j,k,i-1,7) = ux(j,k,i-1,7) + wt(j,k,7)
      ux(j,k,i-1,8) = ux(j,k,i-1,8) + wt(j,k,8)
      ux(j,k,i-1,9) = ux(j,k,i-1,9) + wt(j,k,9)
c
      ux(j,k,i,1) = ux(j,k,i,1) + wt(j,k,1)
      ux(j,k,i,2) = ux(j,k,i,2) + wt(j,k,2)
      ux(j,k,i,3) = ux(j,k,i,3) + wt(j,k,3)
      ux(j,k,i,4) = ux(j,k,i,4) + wt(j,k,4)
      ux(j,k,i,5) = ux(j,k,i,5) + wt(j,k,5)
      ux(j,k,i,6) = ux(j,k,i,6) + wt(j,k,6)
      ux(j,k,i,7) = ux(j,k,i,7) + wt(j,k,7)
      ux(j,k,i,8) = ux(j,k,i,8) + wt(j,k,8)
      ux(j,k,i,9) = ux(j,k,i,9) + wt(j,k,9)
 1700 continue
c
      ii=1
      iip1=min(idim-1,ii+1)
      do 1900 k=1,kdim1
      do 1900 j=1,jdim1
         term       = si(j,k,ii,4)/(voli0(j,k,1)+vol(j,k,ii))
         factor     = bci(j,k,1) + 1.0
         term = term*factor
         wt(j,k,1) = term*( (q(j,k,ii,2)-qi0(j,k,2,1))*si(j,k,ii,1) )
         wt(j,k,2) = term*( (q(j,k,ii,2)-qi0(j,k,2,1))*si(j,k,ii,2) )
         wt(j,k,3) = term*( (q(j,k,ii,2)-qi0(j,k,2,1))*si(j,k,ii,3) )
         wt(j,k,4) = term*( (q(j,k,ii,3)-qi0(j,k,3,1))*si(j,k,ii,1) )
         wt(j,k,5) = term*( (q(j,k,ii,3)-qi0(j,k,3,1))*si(j,k,ii,2) )
         wt(j,k,6) = term*( (q(j,k,ii,3)-qi0(j,k,3,1))*si(j,k,ii,3) )
         wt(j,k,7) = term*( (q(j,k,ii,4)-qi0(j,k,4,1))*si(j,k,ii,1) )
         wt(j,k,8) = term*( (q(j,k,ii,4)-qi0(j,k,4,1))*si(j,k,ii,2) )
         wt(j,k,9) = term*( (q(j,k,ii,4)-qi0(j,k,4,1))*si(j,k,ii,3) )
         ux(j,k,1,1) = ux(j,k,1,1) + wt(j,k,1)
         ux(j,k,1,2) = ux(j,k,1,2) + wt(j,k,2)
         ux(j,k,1,3) = ux(j,k,1,3) + wt(j,k,3)
         ux(j,k,1,4) = ux(j,k,1,4) + wt(j,k,4)
         ux(j,k,1,5) = ux(j,k,1,5) + wt(j,k,5)
         ux(j,k,1,6) = ux(j,k,1,6) + wt(j,k,6)
         ux(j,k,1,7) = ux(j,k,1,7) + wt(j,k,7)
         ux(j,k,1,8) = ux(j,k,1,8) + wt(j,k,8)
         ux(j,k,1,9) = ux(j,k,1,9) + wt(j,k,9)
 1900 continue
c
      ii=idim
      iim2=max(1,ii-2)
      do 2200 k=1,kdim1
      do 2200 j=1,jdim1
         term      = si(j,k,ii,4)/(voli0(j,k,3)+vol(j,k,ii-1))
         factor    = bci(j,k,2) + 1.0
         term = term*factor
         wt(j,k,1) = term*( (qi0(j,k,2,3)-q(j,k,ii-1,2))*si(j,k,ii,1) )
         wt(j,k,2) = term*( (qi0(j,k,2,3)-q(j,k,ii-1,2))*si(j,k,ii,2) )
         wt(j,k,3) = term*( (qi0(j,k,2,3)-q(j,k,ii-1,2))*si(j,k,ii,3) )
         wt(j,k,4) = term*( (qi0(j,k,3,3)-q(j,k,ii-1,3))*si(j,k,ii,1) )
         wt(j,k,5) = term*( (qi0(j,k,3,3)-q(j,k,ii-1,3))*si(j,k,ii,2) )
         wt(j,k,6) = term*( (qi0(j,k,3,3)-q(j,k,ii-1,3))*si(j,k,ii,3) )
         wt(j,k,7) = term*( (qi0(j,k,4,3)-q(j,k,ii-1,4))*si(j,k,ii,1) )
         wt(j,k,8) = term*( (qi0(j,k,4,3)-q(j,k,ii-1,4))*si(j,k,ii,2) )
         wt(j,k,9) = term*( (qi0(j,k,4,3)-q(j,k,ii-1,4))*si(j,k,ii,3) )
         ux(j,k,idim1,1) = ux(j,k,idim1,1) + wt(j,k,1)
         ux(j,k,idim1,2) = ux(j,k,idim1,2) + wt(j,k,2)
         ux(j,k,idim1,3) = ux(j,k,idim1,3) + wt(j,k,3)
         ux(j,k,idim1,4) = ux(j,k,idim1,4) + wt(j,k,4)
         ux(j,k,idim1,5) = ux(j,k,idim1,5) + wt(j,k,5)
         ux(j,k,idim1,6) = ux(j,k,idim1,6) + wt(j,k,6)
         ux(j,k,idim1,7) = ux(j,k,idim1,7) + wt(j,k,7)
         ux(j,k,idim1,8) = ux(j,k,idim1,8) + wt(j,k,8)
         ux(j,k,idim1,9) = ux(j,k,idim1,9) + wt(j,k,9)
 2200 continue
      end if
      end if
c
c     set velocity derivatives for hole cells to one
c
      if (iover.eq.1) then
         do 2403 nk=1,9
         do 2402 i=1,idim1
         do 2402 k=1,kdim1
         do 2402 j=1,jdim1
         ux(j,k,i,nk) = ccvmgt(onec,ux(j,k,i,nk),
     .                  (real(blank(j,k,i)).eq.0.e0))
 2402    continue
 2403    continue
      end if
c
      return
      end
