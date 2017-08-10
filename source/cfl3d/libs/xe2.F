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
      subroutine xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,y6,y7,y8,
     .               yc,xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,
     .               myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set up coefficients for (locally) fitting a polynomial
c     variation in the xie and eta directions .
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /tacos/ iretry
      common /tol/ epsc,epsc0,epsreen,epscoll
      common /zero/iexp
c
      dx2 = x2 - x1
      dy2 = y2 - y1
      dx3 = x3 - x1
      dy3 = y3 - y1
      dx4 = x4 - x1
      dy4 = y4 - y1
      dx5 = x5 - x1
      dy5 = y5 - y1
      dx6 = x6 - x1
      dy6 = y6 - y1
      dx7 = x7 - x1
      dy7 = y7 - y1
      dx8 = x8 - x1
      dy8 = y8 - y1
c
      if (ifit.eq.1) then
c
c     coefficients for bi-linear fit
c
      a2 = dx2
      a3 = dx4
      a4 = dx3 - a2 - a3
      a5 = 0.
      a6 = 0.
      a7 = 0.
      a8 = 0.
      b2 = dy2
      b3 = dy4
      b4 = dy3 - b2 - b3
      b5 = 0.
      b6 = 0.
      b7 = 0.
      b8 = 0.
      end if
c
      if (ifit.eq.2) then
c
c     coefficients for (degenerate) bi-quadradratic fit
c     (quadratic in xie and eta, but without the (xie**2)*(eta**2) term)
c
      a2  = -dx2 + 4.*dx5
      a3  = -dx4 + 4.*dx7
      a5  = 2.*dx2 - 4.*dx5
      a7  = 2.*dx4 - 4.*dx7
      df1 = dx3 - a2 - a3 - a5 - a7
      df2 = dx6 - .5*a2 - a3 - .25*a5 - a7
      df3 = dx8 - a2 - .5*a3 - a5 - .25*a7
      a4  = -3.*df1 + 4.*df2 + 4.*df3
      a6  = 2.*df1 - 4.*df2
      a8  = 2.*df1 - 4.*df3
      b2  = -dy2 + 4.*dy5
      b3  = -dy4 + 4.*dy7
      b5  = 2.*dy2 - 4.*dy5
      b7  = 2.*dy4 - 4.*dy7
      df1 = dy3 - b2 - b3 - b5 - b7
      df2 = dy6 - .5*b2 - b3 - .25*b5 - b7
      df3 = dy8 - b2 - .5*b3 - b5 - .25*b7
      b4  = -3.*df1 + 4.*df2 + 4.*df3
      b6  = 2.*df1 - 4.*df2
      b8  = 2.*df1 - 4.*df3
      end if
c
c     coefficients for quadratic fit in xie, linear fit in eta
c
      if (ifit.eq.3) then
      a3  = dx4
      b3  = dy4
      a2  = -dx2 + 4.*dx5
      b2  = -dy2 + 4.*dy5
      a5  = 2.*dx2  - 4.*dx5
      b5  = 2.*dy2  - 4.*dy5
      df1 = x3 - x2 - a3
      df2 = x6 - x5 - a3
      a4  =   -df1 + 4.*df2
      a6  = 2.*df1 - 4.*df2
      df1 = y3 - y2 - b3
      df2 = y6 - y5 - b3
      b4  =   -df1 + 4.*df2
      b6  = 2.*df1 - 4.*df2
      a7  = 0.
      a8  = 0.
      b7  = 0.
      b8  = 0.
      end if
c
c     coefficients for linear fit in xie,quadratic fit in eta
c
      if (ifit.eq.4) then
      a2  = dx2
      a3  = -dx4 + 4.*dx7
      a7  = 2.*dx4 -4.*dx7
      df1 = dx3 - a2 - a3 - a7
      df2 = dx8 - a2 -.5*a3 - .25*a7
      a4  = -df1 + 4.*df2
      a8  = 2.*df1 - 4.*df2
      a5  = 0.
      a6  = 0.
      b2  = dy2
      b3  = -dy4 + 4.*dy7
      b7  = 2.*dy4 -4.*dy7
      df1 = dy3 - b2 - b3 - b7
      df2 = dy8 - b2 -.5*b3 - .25*b7
      b4  = -df1 + 4.*df2
      b8  = 2.*df1 - 4.*df2
      b5  = 0.
      b6  = 0.
      end if
c
c     newton iteration to invert  for xie and eta (using two 
c     of the three equations defining the surface patch)
c
      if(iretry.eq.0) then
        xie   = 0.5
        eta   = 0.5
      else
        jj  = xie
        kk  = eta
        xie = xie - jj
        eta = eta - kk
      end if
      limit = 50
      iter  = 0
      itry  = 0
c
c     convergence criterion for Newton iteration...require 4 orders of 
c     reduction to starting error=abs(f1)+abs(f2) 
c
      f1 = x1 + a3*eta + eta*( a7*eta + a8*xie*eta )
     .   + xie*( a2 + a4*eta + a5*xie + a6*xie*eta )- xc  
      f2 = y1 + b3*eta + eta*( b7*eta + b8*xie*eta )
     .   + xie*( b2 + b4*eta + b5*xie + b6*xie*eta )- yc  
      error0 = ccabs(f1)+ccabs(f2)
      epsf   = 1.0e-4*error0
c     10.**(-iexp) is machine zero
      epsf1  = max(1.0e-09,10.**(-iexp+1))
      if(real(epsf) .lt. real(epsf1)) epsf = epsf1
c
      errm1 = error0
      error = error0
    2 continue 
c
      errm2 = errm1
      errm1 = error
c
      f1 = x1 + a3*eta + eta*( a7*eta + a8*xie*eta )
     .   + xie*( a2 + a4*eta + a5*xie + a6*xie*eta )- xc  
      f2 = y1 + b3*eta + eta*( b7*eta + b8*xie*eta )
     .   + xie*( b2 + b4*eta + b5*xie + b6*xie*eta )- yc  
c
c d(f1)/d(xie)
c
      a2b = a2 + 2.*xie*( a5 + a6*eta ) + eta*( a4 + a8*eta )
c
c d(f1)/d(eta)
c
      a3b = a3 + xie*( a4 + a6*xie ) + 2.*eta*( a7 + a8*xie ) 
c
c d(f2)/d(xie)
c
      b2b = b2 + 2.*xie*( b5 + b6*eta ) + eta*( b4 + b8*eta )
c
c d(f2)/d(eta)
c
      b3b = b3 + xie*( b4 + b6*xie )
     .     + 2.*eta*( b7 + b8*xie )
c
      iter  = iter + 1
      det   = 1./( a2b*b3b - b2b*a3b )
      xie   = xie - det*( b3b*f1 - a3b*f2 )
      eta   = eta - det*( a2b*f2 - b2b*f1 )
c
c     exit newton iteration if huge xie or eta encountered
c
      huge = 1.e6
      if(abs(real(xie)).ge.real(huge).or.abs(real(eta)).ge.real(huge))
     .  then
        if(abs(real(xie)) .ge. real(huge)) xie = huge
        if(abs(real(eta)) .ge. real(huge)) eta = huge
        imiss = 1
        return
      end if
c
      error = ccabs(f1) + ccabs(f2)
      if (real(error).gt.real(epsf) .and. iter.le.limit) go to 2
c
c     check to ensure point is inside cell
c
      imiss = 0
      if (real(xie).lt.-real(epsc) .or. real(xie).gt.1.+real(epsc) .or.
     .    real(eta).lt.-real(epsc) .or. real(eta).gt.1.+real(epsc)) then
         imiss = 1
      end if
      if (iter.gt.limit .and. imiss.eq.0) then
c
c        newton iteration may be hung, but this may be the correct
c        cell. this behavior is sometimes seen when the code is run in
c        single precision. typically, 4 orders reduction in starting
c        error is required; if iteration is hung and at least 2 orders
c        have been obtained, then just accept it and print a warning. If
c        less than 2 orders reduction in error have been obtained, stop.
c        errm1 = error at iter-1; errm2 = error at iter-2
c
         aerrm1 = ccabs(error-errm1)
         aerrm2 = ccabs(error-errm2)
         if (real(aerrm1).le.real(epsf) .or. real(aerrm2).le.real(epsf))
     .      then
c           iteration is hung
            if (real(error).lt.real(epsf)*100.) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),101)
 101           format('WARNING: Newton iter. hung, but',
     .         ' >= 2 orders reduction in error was obtained')
               return
            else
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),102)
 102           format('stopping...Newton iter. hung with',
     .         ' < 2 orders reduction in error')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         else
c           not converging...must be wrong cell
            xie = xie + 1.
            eta = eta + 1.
            imiss = 1
         end if
      end if
      return
      end
