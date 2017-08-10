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
      subroutine aesurf(nbl,jdim,kdim,idim,deltj,deltk,delti,xmdj,
     .                  xmdk,xmdi,wkj,wkk,wki,maxbl,maxseg,nmds,
     .                  maxaes,aesrfdat,xs,xxn,icsi,icsf,jcsi,
     .                  jcsf,kcsi,kcsf,nsegdfrm,idfrmseg,iaesurf,
     .                  maxsegdg)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Update surface deflections due to aeroelastic motion
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension aesrfdat(5,maxaes)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension xmdj(kdim,idim,6,nmds,maxaes),
     .          xmdk(jdim,idim,6,nmds,maxaes),
     .          xmdi(jdim,kdim,6,nmds,maxaes),
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes)
      dimension wkk(jdim,idim,2),wki(jdim,kdim,2),wkj(kdim,idim,2)
c
      common /elastic/ ndefrm,naesrf
c
      do iaes=1,naesrf
c
         iskyhk = aesrfdat(1,iaes)
         grefl  = aesrfdat(2,iaes)
         uinf   = aesrfdat(3,iaes)
         qinf   = aesrfdat(4,iaes)
         nmodes = aesrfdat(5,iaes)
         rgrefl = 1./grefl
c
         do is=1,nsegdfrm(nbl)
c
         if (iaesurf(nbl,is).eq.iaes .and. idfrmseg(nbl,is).eq.99) then
c
c           convert modal deflections to x,y,z - j=const surfaces
c
            if (jcsi(nbl,is) .eq. jcsf(nbl,is)) then
c
               j  = jcsi(nbl,is)
               ll = 0
               m  = 1
               if (jcsi(nbl,is) .eq.jdim) then
                  ll = 3
                  m  = 2
               end if
               ist = icsi(nbl,is)
               ifn = icsf(nbl,is)
               kst = kcsi(nbl,is)
               kfn = kcsf(nbl,is)
c
               do 3020 i = ist,ifn
               do 3020 k = kst,kfn
c
c                save any existing displacements
                 deltj1 = deltj(k,i,1,m)
                 deltj2 = deltj(k,i,2,m)
                 deltj3 = deltj(k,i,3,m)
c                initialize modal displacements
                 deltj(k,i,1,m) = 0.
                 deltj(k,i,2,m) = 0.
                 deltj(k,i,3,m) = 0.
c                sum delta contributions from all modes
                 do 3010 n = 1,nmodes
                 deltj(k,i,1,m) = deltj(k,i,1,m) +
     .                            wkj(k,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdj(k,i,ll+1,n,iaes)*rgrefl
                 deltj(k,i,2,m) = deltj(k,i,2,m) +
     .                            wkj(k,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdj(k,i,ll+2,n,iaes)*rgrefl
                 deltj(k,i,3,m) = deltj(k,i,3,m) +
     .                            wkj(k,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdj(k,i,ll+3,n,iaes)*rgrefl
3010             continue
c                add modal displacements to existing displacements
                 deltj(k,i,1,m) = deltj(k,i,1,m) + deltj1
                 deltj(k,i,2,m) = deltj(k,i,2,m) + deltj2
                 deltj(k,i,3,m) = deltj(k,i,3,m) + deltj3
                 wkj(k,i,m)    = 0.
3020           continue
c
            end if
c
c           convert modal deflections to x,y,z - k=const surfaces
c
            if (kcsi(nbl,is) .eq. kcsf(nbl,is)) then
c
               k  = kcsi(nbl,is)
               ll = 0
               m  = 1
               if (kcsi(nbl,is) .eq.kdim) then
                  ll = 3
                  m  = 2
               end if
               ist = icsi(nbl,is)
               ifn = icsf(nbl,is)
               jst = jcsi(nbl,is)
               jfn = jcsf(nbl,is)
c
               do 3050 i = ist,ifn
               do 3050 j = jst,jfn
c
c                save any existing displacements
                 deltk1 = deltk(j,i,1,m)
                 deltk2 = deltk(j,i,2,m)
                 deltk3 = deltk(j,i,3,m)
c                initialize modal displacements
                 deltk(j,i,1,m) = 0.
                 deltk(j,i,2,m) = 0.
                 deltk(j,i,3,m) = 0.
c                sum delta contributions from all modes
                 do 3040 n = 1,nmodes
                 deltk(j,i,1,m) = deltk(j,i,1,m) +
     .                            wkk(j,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdk(j,i,ll+1,n,iaes)*rgrefl
                 deltk(j,i,2,m) = deltk(j,i,2,m) +
     .                            wkk(j,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdk(j,i,ll+2,n,iaes)*rgrefl
                 deltk(j,i,3,m) = deltk(j,i,3,m) +
     .                            wkk(j,i,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdk(j,i,ll+3,n,iaes)*rgrefl
3040             continue
c                add modal displacements to existing displacements
                 deltk(j,i,1,m)=deltk(j,i,1,m) + deltk1
                 deltk(j,i,2,m)=deltk(j,i,2,m) + deltk2
                 deltk(j,i,3,m)=deltk(j,i,3,m) + deltk3
                 wkk(j,i,m)   = 0.
3050           continue
c
            end if
c
c           convert modal deflections to x,y,z - i=const surfaces
c
            if (icsi(nbl,is) .eq. icsf(nbl,is)) then
c
               i  = icsi(nbl,is)
               ll = 0
               m  = 1
               if (icsi(nbl,is) .eq.idim) then
                  ll = 3
                  m  = 2
               end if
               kst = kcsi(nbl,is)
               kfn = kcsf(nbl,is)
               jst = jcsi(nbl,is)
               jfn = jcsf(nbl,is)
c
               do 3080 k = kst,kfn
               do 3080 j = jst,jfn
c
c                save any existing displacements
                 delti1 = delti(j,k,1,m)
                 delti2 = delti(j,k,2,m)
                 delti3 = delti(j,k,3,m)
c                initialize modal displacements
                 delti(j,k,1,m) = 0.
                 delti(j,k,2,m) = 0.
                 delti(j,k,3,m) = 0.
c                sum delta contributions from all modes
                 do 3070 n = 1,nmodes
                 delti(j,k,1,m) = delti(j,k,1,m) +
     .                            wki(j,k,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdi(j,k,ll+1,n,iaes)*rgrefl
                 delti(j,k,2,m) = delti(j,k,2,m) +
     .                            wki(j,k,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdi(j,k,ll+2,n,iaes)*rgrefl
                 delti(j,k,3,m) = delti(j,k,3,m) +
     .                            wki(j,k,m)*(xs(2*n-1,iaes)
     .                  - xxn(2*n-1,iaes))*xmdi(j,k,ll+3,n,iaes)*rgrefl
3070             continue
c                add modal displacements to existing displacements
                 delti(j,k,1,m)=delti(j,k,1,m) + delti1
                 delti(j,k,2,m)=delti(j,k,2,m) + delti2
                 delti(j,k,3,m)=delti(j,k,3,m) + delti3
                 wki(j,k,m)   = 0.
3080           continue
c
            end if
c
         end if
c
         end do
c
      end do
c
      return
      end
