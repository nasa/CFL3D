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
      subroutine moddefl(xs,xxn,aesrfdat,perturb,cmyt,cnwt,xorig,yorig,
     .                   zorig,maxaes,nmds,irbtrim,maxbl,myid)
c
c     $Id$
c
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
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes) 
c
      common /elastic/ ndefrm,naesrf
      common /zero/ iexp
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax 
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst 
      common /trim1/ dcl(5000),ddclda(5000),ddcmda(5000),r33,r44,ittr
     .              ,dcm(5000),dd(5000),da(5000),a11,a12,a22,r11,r22
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
         qinf   = aesrfdat(4,iaes)
         ainf   = uinf/xmach
         times  = time*grefl/ainf
         dts    = dt*grefl/ainf
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
               t0              = perturb(nm,iaes,4)
               xs(2*nm-1,iaes) = amp*sin(freqp*(times-t0))
               xs(2*nm,iaes)   = amp*freqp*cos(freqp*(times-t0))
            end if
c
c           Gaussian pulse modal deflection (half is the pulse half-life)
c
            if (moddfl.eq.2) then
               amp             = perturb(nm,iaes,2)
               half            = perturb(nm,iaes,3)
               t0              = perturb(nm,iaes,4)
               const           = log(2.)/half**2
               dtime           = const*(times-t0)**2
               dtime           = ccmin(dtime,dtimemx)
               expterm         = exp(-dtime)
               xs(2*nm-1,iaes) = amp*expterm
               xs(2*nm,iaes)   = -2.*const*(times-t0)*
     .                            xs(2*nm-1,iaes)
            end if
c
c           step pulse modal deflection
c
            if (moddfl.eq.3) then
               amp         = perturb(nm,iaes,2)
               t0          = perturb(nm,iaes,4)
               if (real(times).lt.real(t0-dts/2.))
     .             then
                  xs(2*nm-1,iaes) = 0.
                  xs(2*nm,iaes)   = 0.
               else if (real(times).gt.real(t0-dts/2.) .and. 
     .                  real(times).lt.real(t0+dts/2.)) then 
                  xs(2*nm-1,iaes) = amp
                  xs(2*nm,iaes)   = amp/dts 
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
                if ((abs(real(dlcln)) .lt.real(trtol)) .and.
     .              (abs(real(dmtrmn)).lt.real(trtol)) .and.
     .              (itrmt .gt. 0)) return
                  itrmt = itrmt + 1
                  alf0  = alf1
                  ddtrm0 = xs(2*nm-1,iaes)
c                 call sqrtcumm(itrmt,dcl,dcm,ddclda,ddcmda,da,dd
c    .                         ,dlcln,dmtrmn,tp,fp,a11,a12,a22
c    .                         ,r11,r22,r33,r44,alf1,alf0,ddtrm1
c    .                         ,ddtrm0)
c
c                 Level 1 g flight
c
                  dmtrmnm = dmtrmn
                  qinfrb  = qinf
c
c                 assumes ialph = 0
c                 crefrb was set equal to cref in init_rb
c
                  cmy     = cmyt + cnwt*(xorig-xmc0)/crefrb
                  dmtrmn  = cmy
                  cnw     = cnwt
                  dlclnm  = dlcln
                  dlcln   = cnw - ((tmass*gaccel)/
     .                        (qinfrb*areat*grefl*grefl))
                  if ((abs(real(dlcln)).lt.real(trtol)) .and.
     .                    (abs(real(dmtrmn)).lt.real(trtol))) return
 
                  alf1   =alf0  +relax*(-fp(1,1)*dlcln-fp(1,2)*dmtrmn)
                  ddtrm1 =ddtrm0+relax*(-fp(2,1)*dlcln-fp(2,2)*dmtrmn)
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
                  if(myid.eq.0) then
                    write(79+myid,21928) itrmt,ntt,cmyt,cmy,cnwt,dlcln
     .                            ,ddtrm1*180./3.14159,alf1*180./3.14159
     .                          ,xorig,xmc0,crefrb
21928               format(2i8,9(1x,f16.9))
                  end if
                else
                  xs(2*nm-1,iaes) = xxn(2*nm-1,iaes)
               end if
            end if
c
         end do
      end do
c
      return
      end
      subroutine sqrtcumm(ittr,dcl,dcm,ddclda,ddcmda,da,dd
     .                    ,dlcln,dmtrmn,tp,fp,a11,a12,a22
     .                    ,r11,r22,r33,r44,alf1,alf0,ddtrm1
     .                    ,ddtrm0) 
      dimension dcl(5000),ddclda(5000),ddcmda(5000)
     .         ,dcm(5000),dd(5000),da(5000),tp(4,4),fp(4,4)
      dcl(ittr) = dlcln 
      dcm(ittr) = dmtrmn
      if(ittr.eq.1) then
        da(ittr) = 0.
        dd(ittr) = 0.
      else
        da(ittr) = alf1-alf0
        dd(ittr) = ddtrm1-ddtrm0
        ddclda(ittr) = (dcl(ittr)-dcl(ittr-1))
        ddcmda(ittr) = (dcm(ittr)-dcm(ittr-1))
        a1  = da(ittr-1)*da(ittr-1) 
        a2  = da(ittr-1)*dd(ittr-1)
        a3  = dd(ittr-1)*dd(ittr-1)
        r1  = ddclda(ittr)*da(ittr-1)
        r2  = ddclda(ittr)*dd(ittr-1)
        r3  = ddcmda(ittr)*da(ittr-1)
        r4  = ddcmda(ittr)*dd(ittr-1)
        a11 =a11+a1
        a12 =a12+a2
        a22 =a22+a3
        r11 = r11 + r1
        r22 = r22 + r2
        r33 = r33 + r3
        r44 = r44 + r4
        dtrm = a11*a22-a12*a12
        if(ittr.gt.2) then
          if(dtrm.ne.0.) then
            tp(1,1) =( a22*r11 - a12*r22)/dtrm
            tp(1,2) =(-a12*r11 + a11*r22)/dtrm
            tp(2,1) =( a22*r33 - a12*r44)/dtrm
            tp(2,2) =(-a12*r33 + a11*r44)/dtrm
          end if
          dtr     = tp(1,1)*tp(2,2)-tp(1,2)*tp(2,1)
          if(dtr.ne.0.) then
            fp(1,1) = tp(2,2)/dtr
            fp(1,2) =-tp(1,2)/dtr
            fp(2,1) =-tp(2,1)/dtr
            fp(2,2) = tp(1,1)/dtr
          end if
        end if
      end if
      return
      end
