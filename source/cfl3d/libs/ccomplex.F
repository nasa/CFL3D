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
c
c A collection of routines to substitue for doing complex arithmetic.
c (revert to the standard routines for the real case)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To return a value of a or b, depending on whether the
c     logical expression c is true or false
c***********************************************************************
#   ifdef CMPLX
      complex function ccvmgt(a,b,c)
      implicit complex(a-h,o-z)
#   else
      function ccvmgt(a,b,c)
#   endif
      logical c
c
#   ifdef CRAY
      ccvmgt = merge(a,b,c)
#   else
      if (c) then
         ccvmgt=a
      else
         ccvmgt=b
      end if
#   endif
c
      return
      end
c
c***********************************************************************
c     Purpose: complex/real max
c***********************************************************************
#ifdef CMPLX
      complex function ccmax(a1,a2)
      implicit complex(a-h,o-z)
      if (real(a1).gt.real(a2)) then
        ccmax=a1
      else
         ccmax=a2
      endif
#else
      function ccmax(a1,a2)
      ccmax = max(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real max with 4 complex arguments
c***********************************************************************
#ifdef CMPLX
      complex function ccmax4(a1,a2,a3,a4)
      implicit complex(a-h,o-z)
      dimension a(4)
      a(1)  = a1
      a(2)  = a2
      a(3)  = a3
      a(4)  = a4
      ccmax4 = a(1)
      do n=2,4
         if (real(a(n)).gt.real(ccmax4)) then
            ccmax4  = a(n)
         end if
      end do
#else
      function ccmax4(a1,a2,a3,a4)
      ccmax4 = max(a1,a2,a3,a4)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real max with 1st argument complex, and
c              2nd argument real
c***********************************************************************
#ifdef CMPLX
      complex function ccmaxcr(a1,a2)
      implicit complex(a-h,o-z)
      real a2
      if (real(a1).gt.a2) then
         ccmaxcr=a1
      else
         ccmaxcr=a2
      endif
#else
      function ccmaxcr(a1,a2)
      ccmaxcr = max(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real max with 1st argument real, and
c              2nd argument complex
c***********************************************************************
#ifdef CMPLX
      complex function ccmaxrc(a1,a2)
      implicit complex(a-h,o-z)
      real a1
      if (a1.gt.real(a2)) then
         ccmaxrc=a1
      else
         ccmaxrc=a2
      endif
#else
      function ccmaxrc(a1,a2)
      ccmaxrc = max(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real min
c***********************************************************************
#ifdef CMPLX
      complex function ccmin(a1,a2)
      implicit complex(a-h,o-z)
      if (real(a1).lt.real(a2)) then
         ccmin=a1
      else
         ccmin=a2
      endif
#else
      function ccmin(a1,a2)
      ccmin = min(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real min with 4 complex arguments
c***********************************************************************
#ifdef CMPLX
      complex function ccmin4(a1,a2,a3,a4)
      implicit complex(a-h,o-z)
      dimension a(4)
      a(1)  = a1
      a(2)  = a2
      a(3)  = a3
      a(4)  = a4
      ccmin4 = a(1)
      do n=2,4 
         if (real(a(n)).lt.real(ccmin4)) then
            ccmin4  = a(n)
         end if
      end do
#else
      function ccmin4(a1,a2,a3,a4)
      ccmin4 = min(a1,a2,a3,a4)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real min with 8 complex arguments
c***********************************************************************
#ifdef CMPLX
      complex function ccmin8(a1,a2,a3,a4,a5,a6,a7,a8)
      implicit complex(a-h,o-z)
      dimension a(8)
      a(1)  = a1
      a(2)  = a2
      a(3)  = a3
      a(4)  = a4
      a(5)  = a5
      a(6)  = a6
      a(7)  = a7
      a(8)  = a8
      ccmin8  = a(1)
      do n=2,8
         if (real(a(n)).lt.real(ccmin8)) then
            ccmin8  = a(n)
         end if
      end do
#else
      function ccmin8(a1,a2,a3,a4,a5,a6,a7,a8)
      ccmin8 = min(a1,a2,a3,a4,a5,a6,a7,a8)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real min with 1st argument complex, and
c              2nd argument real
c***********************************************************************
#ifdef CMPLX
      complex function ccmincr(a1,a2)
      implicit complex(a-h,o-z)
      real a2
      if (real(a1).lt.a2) then
         ccmincr=a1
      else
         ccmincr=a2
      endif
#else
      function ccmincr(a1,a2)
      ccmincr = min(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real min with 1st argument real, and
c              2nd argument complex
c***********************************************************************
#ifdef CMPLX
      complex function ccminrc(a1,a2)
      implicit complex(a-h,o-z)
      real a1
      if (a1.lt.real(a2)) then
         ccminrc=a1
      else
         ccminrc=a2
      endif
#else
      function ccminrc(a1,a2)
      ccminrc = min(a1,a2)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real absolute value
c***********************************************************************
#ifdef CMPLX
      complex function ccabs(a)
      implicit complex(a-h,o-z)
      ccabs = a
      if (real(a).lt.0.) ccabs = -a
c     real a1,a2,aa1,aa2
c     a1=real(a)
c     a2=imag(a)
c     aa1=abs(a1)
c     aa2=abs(a2)
c     ccabs=cmplx(aa1,aa2)
c
#else
      function ccabs(a)
      ccabs = abs(a)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real sign
c***********************************************************************
#ifdef CMPLX
      complex function ccsign(a,b)
      implicit complex(a-h,o-z)
      real sgnb
      sgnb = 1.
      if (real(b).lt.0.) sgnb = -1.
      ccsign = ccabs(a)*sgnb
c     ccsign = abs(a)*sgnb
#else
      function ccsign(a,b)
      ccsign = sign(a,b)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real sign with 1st argument real, and 
c     2nd argument complex
c***********************************************************************
#ifdef CMPLX
      complex function ccsignrc(a,b)
      implicit complex(a-h,o-z)
      real a
      real sgnb
      sgnb = 1.
      if (real(b).lt.0.) sgnb = -1.
      ccsignrc = abs(a)*sgnb
#else
      function ccsignrc(a,b)
      ccsignrc = sign(a,b)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex/real hyperbolic tangent
c***********************************************************************
#ifdef CMPLX
      complex function cctanh(a)
      implicit complex(a-h,o-z)
      if (real(a).gt.50) then
         cctanh = 1.
      else if (real(a).lt.-50) then
         cctanh = -1.
      else
         eplus = exp(a)
         eminus = exp(-a)
         cctanh = (eplus - eminus)/(eplus + eminus)
      end if
#else
      function cctanh(a)
      cctanh = tanh(a)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex arc tangent in correct quadrant
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function ccatan2(csn,ccs)
      complex csn,ccs,ccatan
      real pi
      pi = 4.0*atan(1.0)
      if (abs(ccs).eq.0.) then
         if (abs(csn).eq.0.) then
            ccatan2 = cmplx(0.0)
         else
            ccatan2 = cmplx(sign(0.5*pi,real(csn)),0.0)
         end if
      else
         ccatan2 = ccatan(csn/ccs)
         if (real(ccs).lt.0.) ccatan2 = ccatan2 + pi
         if (real(ccatan2).gt.pi) ccatan2 = ccatan2 - 2.0*pi
      end if
#else
      function ccatan2(csn,ccs)
      ccatan2 = atan2(csn,ccs)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex arc tangent (1st and 4th quadrants)
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function ccatan(z)
      complex z,z2
      real pi2
      pi2    = 2.0*atan(1.0)
      r      = abs(z)
      x      = real(z)
      y      = imag(z)
      r2     = r*r
      xans   = 0.5*atan2 (2.0*x, 1.0-r2)
      yans   = 0.25*log((r2+2.0*y+1.0)/(r2-2.0*y+1.0))
      ccatan = cmplx(xans,yans)
#else
      function ccatan(z)
      ccatan = atan(z)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex arccosine
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function ccacos(z)
      complex z,i,one
      i      = (0.0,1.0)
      one    = (1.0,0.0)
      ccacos = -i*log(z+i*sqrt(one-z*z))
#else
      function ccacos(z)
      ccacos = acos(z)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex arcsine
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function ccasin(z)
      complex z,i,one
      i      = (0.0,1.0)
      one    = (1.0,0.0)
      ccasin = -i*log(i*z+sqrt(one-z*z))
#else
      function ccasin(z)
      ccasin = asin(z)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex tangent
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function cctan(z)
      implicit complex(a-h,o-z)
      complex z
      cctan = sin(z)/cos(z)
#else
      function cctan(z)
      cctan = tan(z)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex log base 10
c     David Rodriguez (April 1999)
c***********************************************************************
#ifdef CMPLX
      complex function cclog10(z)
      implicit complex(a-h,o-z)
      cclog10 = log(z)/log((10.0,0.0))
#else
      function cclog10(z)
      cclog10 = log10(z)
#endif
      return
      end
c
c***********************************************************************
c     Purpose: complex positive difference
c***********************************************************************
#ifdef CMPLX
      complex function ccdim(a1,a2)
      implicit complex(a-h,o-z)
      real a1r,a2r
      a1r  = real(a1)
      a2r  = real(a2)
      ccdim = 0
      if(a1r.gt.a2r) ccdim = a1 - a2
#else
      function ccdim(a1,a2)
      ccdim = dim(a1,a2)
#endif
      return
      end
