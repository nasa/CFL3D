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
      subroutine ae_pred(aesrfdat,stm,stmi,gforcn,gforcnm,xs,xxn,
     .                x0,perturb,cmyt,cnwt,xorig,yorig,
     .                zorig,nmds,maxaes,irbtrim,maxbl,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose: either update the modal displacements and velocities via
c              a predictor step of the aeroelastic equations of motion,
c              or specify them (via a call to moddfl) if the modal
c              shapes are subject to forced motion.
c
c        Reference: Cunningham, H.J., Batina, J.T., and Bennett, R.M,
c                  "Modern Wing Flutter Analysis by Computational Fluid
c                   Dynamics Methods," J. Aircraft, Vol. 25, No. 10,
c                   October 1988, pp. 962-968.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension gforcn(2*nmds,maxaes),gforcnm(2*nmds,maxaes),
     .          stm(2*nmds,2*nmds,maxaes),stmi(2*nmds,2*nmds,maxaes),
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes),x0(2*nmds,maxaes),
     .          aesrfdat(5,maxaes),perturb(nmds,maxaes,4)
c
      common /elastic/ ndefrm,naesrf
c
c     modal displacement and velocity prediction via aeroelastic
c     equations of motion
c
      do iaes=1,naesrf
         nmodes = aesrfdat(5,iaes)
         iskyhk = aesrfdat(1,iaes)
         do n=1,2*nmodes
c           don't update if the modal time variation is specified
            moddfl = perturb((n+1)/2,iaes,1)
            if (moddfl .eq. 0) then
               xs(n,iaes) = 0.
               do j=1,2*nmodes
                  xs(n,iaes) = xs(n,iaes) + stm(n,j,iaes)*xxn(j,iaes)
     .                       + .5*stmi(n,j,iaes)*(3.*gforcn(j,iaes)
     .                       - gforcnm(j,iaes))
               end do
            end if
         end do
      end do
c
c     modal displacement and velocity as specified functions of time
c
      call moddefl(xs,xxn,aesrfdat,perturb,cmyt,cnwt,xorig,yorig,
     .             zorig,maxaes,nmds,irbtrim,maxbl,myid)
c
      return
      end
