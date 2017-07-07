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
      subroutine int2_d(iq,jdim1,kdim1,idim1,jmax2,kmax2,nsub1,
     .               jmax1,kmax1,l1,i1,q1,bc2,i2,q2,xie2,eta2,
     .               q1g,q2g,temp,nblkpt,intmax,icheck,mtype,
     .               iindex,ifo,ldim,npt,j21,j22,k21,k22,
     .               q1wk,dthtx,dthty,dthtz,lim_ptch)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Linearly interpolate q from one grid to ghost cells of  
c     another grid using generalized coordinates.
c***********************************************************************
c
c     lim_ptch = 1...limiter employed for interpolation
c                0...no limiter employed
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension iindex(intmax,2*nsub1+9)
      dimension q1(jmax1,kmax1,ldim)
      dimension q2(jmax2,kmax2,ldim,4),bc2(jmax2,kmax2,2)
      dimension xie2(npt),eta2(npt)
      dimension q1g(jmax1,kmax1),q2g(jmax1,kmax1)
      dimension nblkpt(npt)
      dimension temp(jmax1+kmax1)
      dimension q1wk(jmax1,kmax1,ldim)
c
      eps = 1.0e-06
c
c     if non-zero, ensure lim_ptch is 1 and no other value
c
      if (lim_ptch.ne.0) lim_ptch = 1
c
c     load q into work array q1wk
c
      nv = jmax1*kmax1
      do 1000 l = 1,ldim
      do 1010 izz = 1,nv
      q1wk(izz,1,l) = q1(izz,1,l)
 1010 continue
 1000 continue
c
c     use injection (e.g. piecewise constant interpolation) if 
c     solution is first order (ifo=0) 
c
      factfo = 0.
      if (ifo.gt.0) factfo = 1.
c
      if(i2.eq.1 .or. i2.eq.2) i22 = 1
      if(i2.eq.3 .or. i2.eq.4) i22 = 2
c
c     interpolation on an i=constant surface
c
      if(mtype.eq.1) then
c
c     rotate "from" q array, if needed - rotated values stored in q1wk
c     (never need to rotate scalar turbulence quantities - iq <=  0)
c
      if (iq .gt. 0) then
        adthtx = ccabs(dthtx)
        adthty = ccabs(dthty)
        adthtz = ccabs(dthtz)
        if (real(adthtx).gt.0. .or. real(adthty).gt.0. .or.
     .      real(adthtz).gt.0.)then
            jd   = jdim1
            kd   = kdim1
            jsta = 1
            jend = jdim1
            ksta = 1
            kend = kdim1
            call rotateq2_d(jd,kd,q1,q1wk,jsta,jend,
     .                   ksta,kend,dthtx,dthty,dthtz)
        end if
      end if
c
      do 3800 l=1,ldim
c
c     determine gradient values in xie
c
      do 1800 k=1,kmax1-1
      do 1750 j=1,jmax1-2
      q1jk    = q1wk(j,k,l)
      q1j1k   = q1wk(j+1,k,l)
      temp(j) = q1j1k - q1jk
 1750 continue
      do 1751 j=2,jmax1-2
      phi      = temp(j)*temp(j-1)
     .         / (temp(j)**2 + temp(j-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q1g(j,k) = (temp(j) + temp(j-1))*phi
 1751 continue
      q1g(1,k)       = temp(1)
      q1g(jmax1-1,k) = temp(jmax1-2)
 1800 continue
c
c     determine gradient values in eta
c
      do 1900 j=1,jmax1-1
      do 1850 k=1,kmax1-2
      q1jk    = q1wk(j,k,l)
      q1jk1   = q1wk(j,k+1,l)
      temp(k) = q1jk1 - q1jk
 1850 continue
      do 1851 k=2,kmax1-2
      phi      = temp(k)*temp(k-1)
     .         / (temp(k)**2 + temp(k-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q2g(j,k) = (temp(k) + temp(k-1))*phi
 1851 continue
      q2g(j,1)       = temp(1)
      q2g(j,kmax1-1) = temp(kmax1-2)
 1900 continue
 4800 continue
c
c      determine interpolated values
c
      do 2000 k=k21,k22-1
      do 2000 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      lc = nblkpt(ll)
      if (lc.ne.l1) go to 2001
      jc   = int( xie2(ll) )
      kc   = int( eta2(ll) )
c
c     keep within bounds of "from" grid
c
      jc       = max( 1 , jc)
      kc       = max( 1,  kc)
      jc       = min( jc , jmax1-1)
      kc       = min( kc,  kmax1-1)
      xie2(ll) = ccmaxcr(xie2(ll),1.0)
      eta2(ll) = ccmaxcr(eta2(ll),1.0)
      xie2(ll) = ccmincr(xie2(ll),float(jmax1))
      eta2(ll) = ccmincr(eta2(ll),float(kmax1))
c
      xiec = float(jc) + 0.5
      etac = float(kc) + 0.5
c
      bc2(j,k,i22)  = 0.0
      q2(j,k,l,i2) = q1wk(jc,kc,l) + (q1g(jc,kc)*(xie2(ll) - xiec)
     .              +q2g(jc,kc)*(eta2(ll) - etac))*factfo
 2001 continue
 2000 continue
 3800 continue
c
c     interpolation on a j=constant surface
c
      else if(mtype.eq.2) then
c
c     rotate "from" q array, if needed - rotated values stored in q1wk
c     (never need to rotate scalar turbulence quantities - iq <=  0)
c
      if (iq .gt. 0) then
        adthtx = ccabs(dthtx)
        adthty = ccabs(dthty)
        adthtz = ccabs(dthtz)
        if (real(adthtx).gt.0. .or. real(adthty).gt.0. .or.
     .      real(adthtz).gt.0.)then
            jd   = kdim1
            kd   = idim1
            jsta = 1
            jend = kdim1
            ksta = 1
            kend = idim1
            call rotateq2_d(jd,kd,q1,q1wk,jsta,jend,
     .                   ksta,kend,dthtx,dthty,dthtz)
        end if
      end if
c
      do 3810 l=1,ldim
c
c     determine gradient values in xie
c
      do 1810 k=1,kmax1-1
      do 1760 j=1,jmax1-2
      q1jk    = q1wk(j,k,l)
      q1j1k   = q1wk(j+1,k,l)
      temp(j) = q1j1k - q1jk
 1760 continue
      do 1761 j=2,jmax1-2
      phi      = temp(j)*temp(j-1)
     .         / (temp(j)**2 + temp(j-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q1g(j,k) = (temp(j) + temp(j-1))*phi
 1761 continue
      q1g(1,k)       = temp(1)
      q1g(jmax1-1,k) = temp(jmax1-2)
 1810 continue
c
c     determine gradient values in eta
c
      do 1910 j=1,jmax1-1
      do 1860 k=1,kmax1-2
      q1jk    = q1wk(j,k,l)
      q1jk1   = q1wk(j,k+1,l)
      temp(k) = q1jk1 - q1jk
 1860 continue
      do 1861 k=2,kmax1-2
      phi      = temp(k)*temp(k-1)
     .         / (temp(k)**2 + temp(k-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q2g(j,k) = (temp(k) + temp(k-1))*phi
 1861 continue
      q2g(j,1)       = temp(1)
      q2g(j,kmax1-1) = temp(kmax1-2)
 1910 continue
 4810 continue
c
c      determine interpolated values
c
      do 2010 k=k21,k22-1
      do 2010 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      lc = nblkpt(ll)
      if (lc.ne.l1) go to 2011
      jc   = int( xie2(ll) )
      kc   = int( eta2(ll) )
c
c     keep within bounds of "from" grid
c
      jc       = max( 1 , jc)
      kc       = max( 1,  kc)
      jc       = min( jc , jmax1-1)
      kc       = min( kc,  kmax1-1)
      xie2(ll) = ccmaxcr(xie2(ll),1.0)
      eta2(ll) = ccmaxcr(eta2(ll),1.0)
      xie2(ll) = ccmincr(xie2(ll),float(jmax1))
      eta2(ll) = ccmincr(eta2(ll),float(kmax1))
c
      xiec = float(jc) + 0.5
      etac = float(kc) + 0.5
c
      bc2(j,k,i22)  = 0.0
      q2(j,k,l,i2) = q1wk(jc,kc,l) + (q1g(jc,kc)*(xie2(ll) - xiec)
     .              +q2g(jc,kc)*(eta2(ll) - etac))*factfo
 2011 continue
 2010 continue
 3810 continue
c
c     interpolation on a k=constant surface
c
      else if(mtype.eq.3) then
c
c     rotate "from" q array, if needed - rotated values stored in q1wk
c     (never need to rotate scalar turbulence quantities - iq <=  0)
c
      if (iq .gt. 0) then
        adthtx = ccabs(dthtx)
        adthty = ccabs(dthty)
        adthtz = ccabs(dthtz)
        if (real(adthtx).gt.0. .or. real(adthty).gt.0. .or.
     .      real(adthtz).gt.0.)then
            jd   = jdim1
            kd   = idim1
            jsta = 1
            jend = jdim1
            ksta = 1
            kend = idim1
            call rotateq2_d(jd,kd,q1,q1wk,jsta,jend,
     .                   ksta,kend,dthtx,dthty,dthtz)
        end if
      end if
c
      do 3820 l=1,ldim
c
c     determine gradient values in xie
c
      do 1820 k=1,kmax1-1
      do 1770 j=1,jmax1-2
      q1jk    = q1wk(j,k,l)
      q1j1k   = q1wk(j+1,k,l)
      temp(j) = q1j1k - q1jk
 1770 continue
      do 1771 j=2,jmax1-2
      phi      = temp(j)*temp(j-1)
     .         / (temp(j)**2 + temp(j-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q1g(j,k) = (temp(j) + temp(j-1))*phi
 1771 continue
      q1g(1,k)       = temp(1)
      q1g(jmax1-1,k) = temp(jmax1-2)
 1820 continue
c
c     determine gradient values in eta
c
      do 1920 j=1,jmax1-1
      do 1870 k=1,kmax1-2
      q1jk    = q1wk(j,k,l)
      q1jk1   = q1wk(j,k+1,l)
      temp(k) = q1jk1 - q1jk
 1870 continue
      do 1871 k=2,kmax1-2
      phi      = temp(k)*temp(k-1)
     .         / (temp(k)**2 + temp(k-1)**2 + eps)
      phi      = 0.5 + float(lim_ptch)*(phi - 0.5)
      q2g(j,k) = (temp(k) + temp(k-1))*phi
 1871 continue
      q2g(j,1)       = temp(1)
      q2g(j,kmax1-1) = temp(kmax1-2)
 1920 continue
 4820 continue
c
c      determine interpolated values
c
      do 2020 k=k21,k22-1
      do 2020 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      lc = nblkpt(ll)
      if (lc.ne.l1) go to 2021
      jc   = int( xie2(ll) )
      kc   = int( eta2(ll) )
c
c     keep within bounds of "from" grid
c
      jc       = max( 1 , jc)
      kc       = max( 1,  kc)
      jc       = min( jc , jmax1-1)
      kc       = min( kc,  kmax1-1)
      xie2(ll) = ccmaxcr(xie2(ll),1.0)
      eta2(ll) = ccmaxcr(eta2(ll),1.0)
      xie2(ll) = ccmincr(xie2(ll),float(jmax1))
      eta2(ll) = ccmincr(eta2(ll),float(kmax1))
c
      xiec = float(jc) + 0.5
      etac = float(kc) + 0.5
c
      bc2(j,k,i22)  = 0.0
      q2(j,k,l,i2) = q1wk(jc,kc,l) + (q1g(jc,kc)*(xie2(ll) - xiec)
     .              +q2g(jc,kc)*(eta2(ll) - etac))*factfo
 2021 continue
 2020 continue
 3820 continue
      end if     
      return
      end
