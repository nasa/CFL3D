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
      subroutine extrae(jdim,kdim,nsub,l,x,y,z,jcell,kcell,kcl,kcr,
     .                  x7,y7,z7,icase,ifit)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute extra midcell points (x7,y7,z7) at 
c     (xie,eta) = (0.,.5)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,nsub),y(jdim,kdim,nsub),z(jdim,kdim,nsub)
c
c      check if only bilinear fit reasonable 
c
      icase = 999
      x1    = x(jcell,kcell,l)
      y1    = y(jcell,kcell,l)
      z1    = z(jcell,kcell,l)
      x2    = x(jcell,kcell+1,l)
      y2    = y(jcell,kcell+1,l)
      z2    = z(jcell,kcell+1,l)
      x7    = 0.5*( x1 + x2 )
      y7    = 0.5*( y1 + y2 )
      z7    = 0.5*( z1 + z2 )
c     if(ifit.eq.1 .or. ifit.eq.4) return
      if (kcell.lt.kcl .or. kcell.gt.kcr) go to 1500
      if (kcl.eq.kcr) go to 1500
      dxp   = x2 - x1
      dyp   = y2 - y1
      dzp   = z2 - z1
      xlen  = sqrt( dxp*dxp + dyp*dyp +dzp*dzp )
      if(real(xlen).le.0) go to 1500
      if (kcell.eq.kcl) go to 1000
      if (kcell.eq.kcr) go to 2000
c
c        interior points
c
      dxq   = x1 - x(jcell,kcell-1,l)
      dyq   = y1 - y(jcell,kcell-1,l)
      dzq   = z1 - z(jcell,kcell-1,l)
      xlenq = sqrt( dxq*dxq + dyq*dyq +dzq*dzq )
      if(real(xlenq).le.0) go to 1500
      check = ( dxq*dxp + dyq*dyp +dzq*dzp )/(xlenq*xlen)
      if (real(check).lt.0.707107) go to 2000
      xlenq = xlenq / xlen
      if (real(xlenq).lt.0.333 .or. real(xlenq).gt.3.) go to 2000
      xlens = - 1.
      b2    = x(jcell,kcell-1,l) - x1 - (x2-x1)*xlens
      c2    = y(jcell,kcell-1,l) - y1 - (y2-y1)*xlens
      d2    = z(jcell,kcell-1,l) - z1 - (z2-z1)*xlens
      a2    = xlens*( 1. - xlens )
      x7    = x7 +.25*b2/a2
      y7    = y7 +.25*c2/a2
      z7    = z7 +.25*d2/a2
      dxr   = x(jcell,kcell+2,l) - x2
      dyr   = y(jcell,kcell+2,l) - y2
      dzr   = z(jcell,kcell+2,l) - z2
      xlenr = sqrt( dxr*dxr + dyr*dyr +dzr*dzr )
      if(real(xlenr).le.0) go to 1500
      check = ( dxr*dxp + dyr*dyp +dzr*dzp )/(xlenr*xlen)
      if (real(check).lt.0.707107) go to 1500
      xlenr = xlenr / xlen
      if (real(xlenr).lt.0.333 .or. real(xlenr).gt.3.) go to 1500
      xlens = 2.
      b1    = x(jcell,kcell+2,l) - x1 - (x2-x1)*xlens
      c1    = y(jcell,kcell+2,l) - y1 - (y2-y1)*xlens
      d1    = z(jcell,kcell+2,l) - z1 - (z2-z1)*xlens
      a1    = xlens*( 1. - xlens ) 
      trat  = xlenq*xlenq/(xlenr*xlenr)
      trat  = 1.
      term  = .25/( a1*a1 + trat*a2*a2 )
      x7    = 0.5*( x1 + x2 ) + term*( b1*a1 + b2*a2*trat )
      y7    = 0.5*( y1 + y2 ) + term*( c1*a1 + c2*a2*trat )
      z7    = 0.5*( z1 + z2 ) + term*( d1*a1 + d2*a2*trat )
      icase = 0
      go to 1500
 1000 continue
c
c     left edge
c
      icase = -1
      dxr   = x(jcell,kcell+2,l) - x2
      dyr   = y(jcell,kcell+2,l) - y2
      dzr   = z(jcell,kcell+2,l) - z2
      xlenr = sqrt( dxr*dxr + dyr*dyr +dzr*dzr )
      if(real(xlenr).le.0) go to 1500
      check = ( dxr*dxp + dyr*dyp +dzr*dzp )/(xlenr*xlen)
      if (real(check).lt.0.707107) go to 1500
      xlenr = xlenr / xlen
      if (real(xlenr).lt.0.333 .or. real(xlenr).gt.3.) go to 1500
      xlens = 2.
      b1    = x(jcell,kcell+2,l) - x1 - (x2-x1)*xlens
      c1    = y(jcell,kcell+2,l) - y1 - (y2-y1)*xlens
      d1    = z(jcell,kcell+2,l) - z1 - (z2-z1)*xlens
      a1    = 1./( xlens*(1. - xlens ) )
      x7    = x7 +.25*( b1 * a1 )
      y7    = y7 +.25*( c1 * a1 )
      z7    = z7 +.25*( d1 * a1 )
      icase = 1
      go to 1500
 2000 continue
c
c     right edge
c
      icase = -2
      dxq   = x1 - x(jcell,kcell-1,l)
      dyq   = y1 - y(jcell,kcell-1,l)
      dzq   = z1 - z(jcell,kcell-1,l)
      xlenq = sqrt( dxq*dxq + dyq*dyq +dzq*dzq )
      if(real(xlenq).le.0) go to 1500
      check = ( dxq*dxp + dyq*dyp +dzq*dzp )/(xlenq*xlen)
      if (real(check).lt.0.707107) go to 1500
      xlenq = xlenq / xlen
      if (real(xlenq).lt.0.333 .or. real(xlenq).gt.3.) go to 1500
      xlens = -1.
      b2    = x(jcell,kcell-1,l) - x1 - (x2-x1)*xlens
      c2    = y(jcell,kcell-1,l) - y1 - (y2-y1)*xlens
      d2    = z(jcell,kcell-1,l) - z1 - (z2-z1)*xlens
      a2    = 1./( xlens*( 1. - xlens ) )
      x7    = x7 +.25*( b2 * a2 )
      y7    = y7 +.25*( c2 * a2 )
      z7    = z7 +.25*( d2 * a2 )
      icase = 2
 1500 continue
      return
      end
