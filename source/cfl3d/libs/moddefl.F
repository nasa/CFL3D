      subroutine moddefl(xs,aesrfdat,perturb,cmyt,cnwt,xorig,yorig,
     .                   zorig,maxaes,nmds,irbtrim,maxbl)
c***********************************************************************
c     Purpose: Specify temporal variation of modal deflections for
c              either flutter initiation or determination of generalized
c              force response
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension aesrfdat(5,maxaes),perturb(nmds,maxaes,4),
     .          xs(2*nmds,maxaes)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
c
      common /elastic/ ndefrm,naesrf
      common /zero/ iexp
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4)
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
c
      dtimemx  = -log(10.**(-iexp))
c
      do iaes=1,naesrf
c
         nmodes = aesrfdat(5,iaes)
         grefl  = aesrfdat(2,iaes)
         uinf   = aesrfdat(3,iaes)
c
         do nm=1,nmodes
c
            moddfl = perturb(nm,iaes,1)
c
c           fixed mode
c
            if (moddfl.lt.0) then
               xs(2*nm-1,iaes) = 0.
               xs(2*nm,iaes)   = 0.
            end if
c
c           harmonic modal oscillation
c
            if (moddfl.eq.1) then
               amp             = perturb(nm,iaes,2)
               freqp           = perturb(nm,iaes,3)
               xs(2*nm-1,iaes) = amp*sin(freqp*time)
               xs(2*nm,iaes)   = amp*freqp*cos(freqp*time)
            end if
c
c           Gaussian pulse modal deflection (half is the pulse half-life)
c
            if (moddfl.eq.2) then
               amp             = perturb(nm,iaes,2)
               half            = perturb(nm,iaes,3)
               t0              = perturb(nm,iaes,4)
               const           = log(2.)/half**2
               dtime           = const*(time-t0)**2
               dtime           = ccmin(dtime,dtimemx)
               expterm         = exp(-dtime)
               xs(2*nm-1,iaes) = amp*expterm
               xs(2*nm,iaes)   = -2.*const*(time-t0)*xs(2*nm-1,iaes)
            end if
c
c           step pulse modal deflection
c
            if (moddfl.eq.3) then
               amp         = perturb(nm,iaes,2)
               t0          = perturb(nm,iaes,4)
               if (real(time).lt.real(t0-dt/2.)) then
                  xs(2*nm-1,iaes) = 0.
                  xs(2*nm,iaes)   = 0.
               else if (real(time).gt.real(t0-dt/2.) .and. 
     .                  real(time).lt.real(t0+dt/2.)) then 
                  xs(2*nm-1,iaes) = amp
                  xs(2*nm,iaes)   = amp/dt
               else
                  xs(2*nm-1,iaes) = amp
                  xs(2*nm,iaes)   = 0.
               end if
            end if
c
c           rigid body mode
c
            if (moddfl.eq.4) then
c
               if (ntt/itrminc*itrminc .eq. ntt) then
c
c                 Update the angle of attack and control surface trim
c                 angle to reach the trim lift coefficient
c                 (tmass*g/(qinf*areat*gfrefl*gfrefl))
c                 and zero moment coefficient.
c
c                 Level 1 g flight
c
                  dmtrmnm = dmtrmn
c
c                 assumes ialph = 0
c                 crefrb was set equal to cref in init_rb
                  cmy     = cmyt + cnwt*(xorig(1)-xmc0)/crefrb
                  dmtrmn  = cmy
                  cnw     = cnwt
                  dlclnm  = dlcln
                  dlcln   = cnw - ((tmass*gaccel)/
     .                        (2.*qinfrb*areat*greflrb*greflrb))
                  if ((abs(real(dlcln)).lt.real(trtol)) .and.
     .                    (abs(real(dmtrmn)).lt.real(trtol))) return

                  ddtrm0 = xs(2*nm-1,iaes)
                  alf0   = alf1
                  relax  = .20
                  alf1   =alf0  +relax*(-fp(1,1)*dlcln-fp(1,2)*dmtrmn)
                  ddtrm1 =ddtrm0+relax*(-fp(2,1)*dlcln-fp(2,2)*dmtrmn)
                  write(79,21928) itrmt,ntt,cmy,ddtrm0,ddtrm1,dlcln
     .                            ,alf0,alf1
                  if(real(ddtrm1).gt.real(dtrmsmx)) ddtrm1 = dtrmsmx
                  if(real(ddtrm1).lt.real(dtrmsmn)) ddtrm1 = dtrmsmn
                  if(abs(real(ddtrm1-ddtrm0)).gt.real(ddtmx)) then
                     ddtrm1 = ddtrm0 + ddtmx*(ddtrm1-ddtrm0)
     .                                    /ccabs(ddtrm1-ddtrm0)
                  end if
                  xs(2*nm-1,iaes) = ddtrm1
                  if (abs(real(alf1-alf0)).gt.real(dalfmx)) then
                        alf1 = alf0+dalfmx*(alf1-alf0)/ccabs(alf1-alf0)
                  end if
21928             format(2i8,6(1x,f16.9))
                  itrmt = itrmt + 1
                  irbtrim = 1
               end if
            end if
c
         end do
      end do
c
      return
      end
