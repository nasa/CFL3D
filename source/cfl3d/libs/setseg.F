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
      subroutine setseg(maxgr,maxbl,maxseg,nblg,ncgg,idimg,jdimg,
     .                  kdimg,ibcinfo,jbcinfo,kbcinfo,nbci0,
     .                  nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                  ndefrm0,idefrm,nsegdfrm,icsi,icsf,
     .                  jcsi,jcsf,kcsi,kcsf,maxsegdg,ngrid)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set deforming segment counter and index ranges to
c     correspond to all solid surfaces.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ncgg(maxgr),nblg(maxgr),idimg(maxbl),jdimg(maxbl),
     .          kdimg(maxbl)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          jbcinfo(maxbl,maxseg,7,2),ibcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idefrm(maxbl)

      do ng=1,ngrid
         nbl  = nblg(ng)
         ncg  = ncgg(ng)
         idim = idimg(nbl)
         jdim = jdimg(nbl)
         kdim = kdimg(nbl)
         do nseg=1,nbci0(nbl)
            if (abs(ibcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2016 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = 1
               icsf(nbl,iseg)     = 1
               jcsi(nbl,iseg)     = ibcinfo(nbl,nseg,2,1)
               jcsf(nbl,iseg)     = ibcinfo(nbl,nseg,3,1)
               kcsi(nbl,iseg)     = ibcinfo(nbl,nseg,4,1)
               kcsf(nbl,iseg)     = ibcinfo(nbl,nseg,5,1)
            end if
         end do
         do nseg=1,nbcidim(nbl)
            if (abs(ibcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2016 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = idimg(nbl)
               icsf(nbl,iseg)     = idimg(nbl)
               jcsi(nbl,iseg)     = ibcinfo(nbl,nseg,2,2)
               jcsf(nbl,iseg)     = ibcinfo(nbl,nseg,3,2)
               kcsi(nbl,iseg)     = ibcinfo(nbl,nseg,4,2)
               kcsf(nbl,iseg)     = ibcinfo(nbl,nseg,5,2)
            end if
         end do
         do nseg=1,nbcj0(nbl)
            if (abs(jbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2016 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = jbcinfo(nbl,nseg,2,1)
               icsf(nbl,iseg)     = jbcinfo(nbl,nseg,3,1)
               jcsi(nbl,iseg)     = 1
               jcsf(nbl,iseg)     = 1
               kcsi(nbl,iseg)     = jbcinfo(nbl,nseg,4,1)
               kcsf(nbl,iseg)     = jbcinfo(nbl,nseg,5,1)
            end if
         end do
         do nseg=1,nbcjdim(nbl)
            if (abs(jbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2016 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = jbcinfo(nbl,nseg,2,2)
               icsf(nbl,iseg)     = jbcinfo(nbl,nseg,3,2)
               jcsi(nbl,iseg)     = jdimg(nbl)
               jcsf(nbl,iseg)     = jdimg(nbl)
               kcsi(nbl,iseg)     = jbcinfo(nbl,nseg,4,2)
               kcsf(nbl,iseg)     = jbcinfo(nbl,nseg,5,2)
            end if
         end do
         do nseg=1,nbck0(nbl)
            if (abs(kbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2016 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = kbcinfo(nbl,nseg,2,1)
               icsf(nbl,iseg)     = kbcinfo(nbl,nseg,3,1)
               jcsi(nbl,iseg)     = kbcinfo(nbl,nseg,4,1)
               jcsf(nbl,iseg)     = kbcinfo(nbl,nseg,5,1)
               kcsi(nbl,iseg)     = 1
               kcsf(nbl,iseg)     = 1
            end if
         end do
         do nseg=1,nbckdim(nbl)
            if (abs(kbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2016 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.1006) then
               ndefrm0           = ndefrm0 + 1
               idefrm(nbl)        = 1
               iseg               = nsegdfrm(nbl) + 1
               nsegdfrm(nbl)      = iseg
               icsi(nbl,iseg)     = kbcinfo(nbl,nseg,2,2)
               icsf(nbl,iseg)     = kbcinfo(nbl,nseg,3,2)
               jcsi(nbl,iseg)     = kbcinfo(nbl,nseg,4,2)
               jcsf(nbl,iseg)     = kbcinfo(nbl,nseg,5,2)
               kcsi(nbl,iseg)     = kdimg(nbl)
               kcsf(nbl,iseg)     = kdimg(nbl)
            end if
         end do
         if (ncg.gt.0) then
            do n=1,ncg
               nbl = nbl + 1
               idefrm(nbl)   = idefrm(nbl-1)
               nsegdfrm(nbl) = nsegdfrm(nbl-1)
               do iseg=1,nsegdfrm(nbl)
                  icsi(nbl,iseg) = icsi(nbl-1,iseg)/2 + 1
                  icsf(nbl,iseg) = icsf(nbl-1,iseg)/2 + 1
                  jcsi(nbl,iseg) = jcsi(nbl-1,iseg)/2 + 1
                  jcsf(nbl,iseg) = jcsf(nbl-1,iseg)/2 + 1
                  kcsi(nbl,iseg) = kcsi(nbl-1,iseg)/2 + 1
                  kcsf(nbl,iseg) = kcsf(nbl-1,iseg)/2 + 1
               end do
            end do
         end if
      end do
c
      return
      end
