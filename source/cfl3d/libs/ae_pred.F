      subroutine ae_pred(aesrfdat,stm,stmi,gforcn,gforcnm,xs,xxn,
     .                x0,perturb,cmyt,cnwt,xorig,yorig,
     .                zorig,nmds,maxaes,irbtrim,maxbl)
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
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
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
      call moddefl(xs,aesrfdat,perturb,cmyt,cnwt,xorig,yorig,
     .             zorig,maxaes,nmds,irbtrim,maxbl)
c
      return
      end
