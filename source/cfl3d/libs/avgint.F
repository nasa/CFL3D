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
      subroutine avgint(xiec,etac,nblkc,nptc,xief,etaf,nblkf,nptf,
     .               j1c,j2c,k1c,k2c,j1f,j2f,k1f,k2f)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Create interpolation data on coarser levels by
c               averaging finer-level interpolation data
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xiec(nptc),etac(nptc),nblkc(nptc),xief(nptf),
     .          etaf(nptf),nblkf(nptf)
c
      do kc = k1c,k2c-1
         kf = (kc-1)*2 + 1
         kfp = kf + 1
         kfp = min(kfp,k2f-1)
         do jc = j1c,j2c-1
            jf  = (jc-1)*2 + 1
            jfp = jf + 1
            jfp = min(jfp,j2f-1)
            llc = (j2c-j1c)*(kc-k1c) + (jc-j1c+1)
            llf = (j2f-j1f)*(kf-k1f) + (jf-j1f+1)
            llfjp = (j2f-j1f)*(kf-k1f) + (jfp-j1f+1)
            llfkp = (j2f-j1f)*(kfp-k1f) + (jf-j1f+1)
            llfjpkp = (j2f-j1f)*(kfp-k1f) + (jfp-j1f+1)
c
c           all finer points must lie in same block for averaging;
c           otherwise, use only data from lower corner point
c
            if (nblkf(llfjp)   .ne. nblkf(llf)) llfjp   = llf
            if (nblkf(llfkp)   .ne. nblkf(llf)) llfkp   = llf
            if (nblkf(llfjpkp) .ne. nblkf(llf)) llfjpkp = llf
c
            xiefavg = 0.25*(xief(llf)   + xief(llfjp) 
     .                    + xief(llfkp) + xief(llfjpkp))
            etafavg = 0.25*(etaf(llf)   + etaf(llfjp)
     .                    + etaf(llfkp) + etaf(llfjpkp))
            jfhat = int(xiefavg)
            xiefhat = xiefavg - jfhat
            kfhat = int(etafavg)
            etafhat = etafavg - kfhat
            jchat = (jfhat+1)/2
            kchat = (kfhat+1)/2
            if (j1f.eq.1 .and. j2f.eq.2) then
c              2d case
               xiechat = xiefhat
            else
               if (jfhat/2*2 .eq. jfhat) then
                  xiechat = 0.5*(1.+xiefhat)
               else
                  xiechat = 0.5*xiefhat
               end if
            end if
            if (k1f.eq.1 .and. k2f.eq.2) then
c              2d case
               etachat = etafhat
            else
               if (kfhat/2*2 .eq. kfhat) then
                  etachat = 0.5*(1.+etafhat)
               else
                  etachat = 0.5*etafhat
               end if
            end if
            xiec(llc) = jchat + xiechat
            etac(llc) = kchat + etachat
            nblkc(llc) = nblkf(llf)
         end do
      end do
c
      return
      end
