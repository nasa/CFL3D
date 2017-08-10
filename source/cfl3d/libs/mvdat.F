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
      subroutine mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                 omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                 dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                 thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,time2,
     .                 thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                 aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                 maxsegdg,wk,nwork)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Move dynamic mesh data to/from common blocks
c     into an array for sending data as a single message
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),omegax(maxbl),
     .          omegay(maxbl),omegaz(maxbl),xorig(maxbl),yorig(maxbl),
     .          zorig(maxbl),dxmx(maxbl),dymx(maxbl),dzmx(maxbl),
     .          dthxmx(maxbl),dthymx(maxbl),dthzmx(maxbl),thetax(maxbl),
     .          thetay(maxbl),thetaz(maxbl),rfreqt(maxbl),rfreqr(maxbl),
     .          xorig0(maxbl),yorig0(maxbl),zorig0(maxbl),time2(maxbl),
     .          thetaxl(maxbl),thetayl(maxbl),thetazl(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4)
      dimension aesrfdat(5,maxaes),wk(nwork)
c
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
c
      if (idir .eq. 0) then
         wk( 1) = float(itrans(nbl))
         wk( 2) = rfreqt(nbl)
         wk( 3) = xorig(nbl)
         wk( 4) = yorig(nbl)
         wk( 5) = zorig(nbl)
         wk( 6) = xorig0(nbl)
         wk( 7) = yorig0(nbl)
         wk( 8) = zorig0(nbl)
         wk( 9) = utrans(nbl)
         wk(10) = vtrans(nbl)
         wk(11) = wtrans(nbl)
         wk(12) = dxmx(nbl)
         wk(13) = dymx(nbl)
         wk(14) = dzmx(nbl)
         wk(15) = float(itransmc)
         wk(16) = rfreqtmc
         wk(17) = xorigmc
         wk(18) = yorigmc
         wk(19) = zorigmc
         wk(20) = xorig0mc
         wk(21) = yorig0mc
         wk(22) = zorig0mc
         wk(23) = utransmc
         wk(24) = vtransmc
         wk(25) = wtransmc
         wk(26) = xmc
         wk(27) = ymc
         wk(28) = zmc
         wk(29) = dxmxmc
         wk(30) = dymxmc
         wk(31) = dzmxmc
         wk(32) = float(irotat(nbl))
         wk(33) = rfreqr(nbl)
         wk(34) = thetax(nbl)
         wk(35) = thetay(nbl)
         wk(36) = thetaz(nbl)
         wk(37) = omegax(nbl)
         wk(38) = omegay(nbl)
         wk(39) = omegaz(nbl)
         wk(40) = dthxmx(nbl)
         wk(41) = dthymx(nbl)
         wk(42) = dthzmx(nbl)
         wk(43) = float(irotatmc)
         wk(44) = rfreqrmc
         wk(45) = thetaxmc
         wk(46) = thetaymc
         wk(47) = thetazmc
         wk(48) = omegaxmc
         wk(49) = omegaymc
         wk(50) = omegazmc
         wk(51) = dthxmxmc
         wk(52) = dthymxmc
         wk(53) = dthzmxmc
         wk(54) = time2(nbl)
         wk(55) = time2mc
         wk(56) = dt
         wk(57) = float(idefrm(nbl))
         wk(58) = float(nsegdfrm(nbl))
         nlast  = 58
         nseg   = nsegdfrm(nbl)
         do is=1,nseg
            wk(is+nlast)         = utrnsae(nbl,is)
            wk(is+nlast+   nseg) = vtrnsae(nbl,is)
            wk(is+nlast+ 2*nseg) = wtrnsae(nbl,is)
            wk(is+nlast+ 3*nseg) = omgxae(nbl,is)
            wk(is+nlast+ 4*nseg) = omgyae(nbl,is)
            wk(is+nlast+ 5*nseg) = omgzae(nbl,is)
            wk(is+nlast+ 6*nseg) = xorgae(nbl,is)
            wk(is+nlast+ 7*nseg) = yorgae(nbl,is)
            wk(is+nlast+ 8*nseg) = zorgae(nbl,is)
            wk(is+nlast+ 9*nseg) = thtxae(nbl,is)
            wk(is+nlast+10*nseg) = thtyae(nbl,is)
            wk(is+nlast+11*nseg) = thtzae(nbl,is)
            wk(is+nlast+12*nseg) = rfrqtae(nbl,is)
            wk(is+nlast+13*nseg) = rfrqrae(nbl,is)
            wk(is+nlast+14*nseg) = float(icsi(nbl,is))
            wk(is+nlast+15*nseg) = float(icsf(nbl,is))
            wk(is+nlast+16*nseg) = float(jcsi(nbl,is))
            wk(is+nlast+17*nseg) = float(jcsf(nbl,is))
            wk(is+nlast+18*nseg) = float(kcsi(nbl,is))
            wk(is+nlast+19*nseg) = float(kcsf(nbl,is))
         end do
      else
         itrans(nbl)   = int(wk( 1))
         rfreqt(nbl)   = wk( 2)
         xorig(nbl)    = wk( 3)
         yorig(nbl)    = wk( 4)
         zorig(nbl)    = wk( 5)
         xorig0(nbl)   = wk( 6)
         yorig0(nbl)   = wk( 7)
         zorig0(nbl)   = wk( 8)
         utrans(nbl)   = wk( 9)
         vtrans(nbl)   = wk(10)
         wtrans(nbl)   = wk(11)
         dxmx(nbl)     = wk(12)
         dymx(nbl)     = wk(13)
         dzmx(nbl)     = wk(14)
         itransmc      = int(wk(15))
         rfreqtmc      = wk(16)
         xorigmc       = wk(17)
         yorigmc       = wk(18)
         zorigmc       = wk(19)
         xorig0mc      = wk(20)
         yorig0mc      = wk(21)
         zorig0mc      = wk(22)
         utransmc      = wk(23)
         vtransmc      = wk(24)
         wtransmc      = wk(25)
         xmc           = wk(26)
         ymc           = wk(27)
         zmc           = wk(28)
         dxmxmc        = wk(29)
         dymxmc        = wk(30)
         dzmxmc        = wk(31)
         irotat(nbl)   = int(wk(32))
         rfreqr(nbl)   = wk(33)
         thetax(nbl)   = wk(34)
         thetay(nbl)   = wk(35)
         thetaz(nbl)   = wk(36)
         omegax(nbl)   = wk(37)
         omegay(nbl)   = wk(38)
         omegaz(nbl)   = wk(39)
         dthxmx(nbl)   = wk(40)
         dthymx(nbl)   = wk(41)
         dthzmx(nbl)   = wk(42)
         irotatmc      = int(wk(43))
         rfreqrmc      = wk(44)
         thetaxmc      = wk(45)
         thetaymc      = wk(46)
         thetazmc      = wk(47)
         omegaxmc      = wk(48)
         omegaymc      = wk(49)
         omegazmc      = wk(50)
         dthxmxmc      = wk(51)
         dthymxmc      = wk(52)
         dthzmxmc      = wk(53)
         time2(nbl)    = wk(54)
         time2mc       = wk(55)
         dt            = wk(56)
         idefrm(nbl)   = int(wk(57))
         nsegdfrm(nbl) = int(wk(58))
         nlast  = 58
         nseg   = nsegdfrm(nbl)
         do is=1,nseg
            utrnsae(nbl,is) = wk(is+nlast)
            vtrnsae(nbl,is) = wk(is+nlast+   nseg)
            wtrnsae(nbl,is) = wk(is+nlast+ 2*nseg)
            omgxae(nbl,is)  = wk(is+nlast+ 3*nseg)
            omgyae(nbl,is)  = wk(is+nlast+ 4*nseg)
            omgzae(nbl,is)  = wk(is+nlast+ 5*nseg)
            xorgae(nbl,is)  = wk(is+nlast+ 6*nseg)
            yorgae(nbl,is)  = wk(is+nlast+ 7*nseg) 
            zorgae(nbl,is)  = wk(is+nlast+ 8*nseg)
            thtxae(nbl,is)  = wk(is+nlast+ 9*nseg)
            thtyae(nbl,is)  = wk(is+nlast+10*nseg)
            thtzae(nbl,is)  = wk(is+nlast+11*nseg)
            rfrqtae(nbl,is) = wk(is+nlast+12*nseg)
            rfrqrae(nbl,is) = wk(is+nlast+13*nseg)
            icsi(nbl,is)    = int(wk(is+nlast+14*nseg))
            icsf(nbl,is)    = int(wk(is+nlast+15*nseg))
            jcsi(nbl,is)    = int(wk(is+nlast+16*nseg))
            jcsf(nbl,is)    = int(wk(is+nlast+17*nseg))
            kcsi(nbl,is)    = int(wk(is+nlast+18*nseg))
            kcsf(nbl,is)    = int(wk(is+nlast+19*nseg))
         end do
      end if
c
      return
      end
