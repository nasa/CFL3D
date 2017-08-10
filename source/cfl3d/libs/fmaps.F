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
      subroutine fmaps(nx,ny,nz,area,at,f,qr,ql,n,nvtq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute flux at the interface using Cord Rossow's
c     MAPS+ scheme, given the left and right  states at the interface.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#   ifdef CMPLX
      complex m_r,m_l,m_av,m_up,m_lo,m_max1,m_min1,m_dif,m_upd,m_max,
     .     m_0,m_m1
      complex nx(n),ny(n),nz(n)
#   else
      real m_r,m_l,m_av,m_up,m_lo,m_max1,m_min1,m_dif,m_upd,m_max,
     .     m_0,m_m1
      real nx(n),ny(n),nz(n)
#   endif
c
      dimension area(n),at(n),qr(nvtq,5),ql(nvtq,5),f(nvtq,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /precond/ cprec,uref,avn
c
      x1   = gamma/gm1
      c1   = 1.e0/gm1
      eps4 = 0.e0
c
cdir$ ivdep
      do 1000 i=1,n 
c
c     ------------------------------------------------------------------
c     left/right density, velocity, pressure, enthalpy and sound speed
c     ------------------------------------------------------------------
c
      r_l  = ql(i,1)
      vx_l = ql(i,2)
      vy_l = ql(i,3)
      vz_l = ql(i,4)
      p_l  = ql(i,5)
      h_l  = x1*ql(i,5)/ql(i,1)+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                              +ql(i,4)*ql(i,4))
      c_l  = sqrt(gamma*ql(i,5)/ql(i,1))
c
      r_r  = qr(i,1)
      vx_r = qr(i,2)
      vy_r = qr(i,3)
      vz_r = qr(i,4)
      p_r  = qr(i,5)
      h_r  = x1*qr(i,5)/qr(i,1)+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                              +qr(i,4)*qr(i,4))
      c_r  = sqrt(gamma*qr(i,5)/qr(i,1))
c
      c_av = 0.5*(c_l+c_r)
c
c     nx, ny, nz are the components of the normal vector at a cell face
c     anx, any, anz are the face vector components { (nx, ny, nz) * area }
c
      anx = nx(i)*area(i)
      any = ny(i)*area(i)
      anz = nz(i)*area(i)
c
      q_l   = (vx_l * nx(i) + vy_l * ny(i) + vz_l * nz(i))
      q_r   = (vx_r * nx(i) + vy_r * ny(i) + vz_r * nz(i))
c      
      m_l    = q_l / c_av
      m_r    = q_r / c_av
c
c     ------------------------------------------------------------------
c     compute conditions at cell face:
c     ------------------------------------------------------------------
c
      am_l   = ccabs(m_l)
      am_r   = ccabs(m_r)
c
      sm_l   = ccsignrc( 1.0, m_l)
      sm_r   = ccsignrc( 1.0, m_r)
c
      m_av   = 0.5*(m_l + m_r)
      m_up   = ccmax(am_l, am_r)
      m_lo   = ccmin(am_l, am_r)
      m_max1 = ccmincr(m_up, 1.0)
      m_min1 = ccmincr(m_lo, 1.0)
      m_dif  = 0.5*(am_r - am_l)
c
      b_m    = ccmaxrc(0.0, 2.0*m_max1 - 1.0)
      b_p    = ccmaxrc(0.0, 2.0*m_min1 - 1.0)
c
c     no lower bound on dissipation
c
      m_upd  = m_up
      dfak   = 1.0
c
c     take aspect ratio into account
c
c     fasp   = scale**zeta
c     dfak   = ccmaxcr(fasp, 1.0)
c     dcut   = ccmincr(fasp, 1.0)
c     dcut   = ccmax(hfl1*dcut, hfl2)
c     m_upd  = ccmax(m_up, dcut):
c
      q_cd   = 0.5 * area(i) * c_av * ( m_av - b_m*m_dif )
      q_up   = 0.5 * area(i) * c_av *   m_upd * dfak
c
      p_av  = 0.5*(p_r+p_l)
      p_ds  = 0.5*(p_r*sm_r - p_l*sm_l)
      p_12  = p_av - b_p*p_ds
c                                    
c     ------------------------------------------------------------------
c     compute flux through cell face (maps scheme)
c     ------------------------------------------------------------------
c
      f(i,1)       = q_cd * (r_r      + r_l) -
     .               q_up * (r_r      - r_l)
      f(i,2)       = q_cd * (r_r*vx_r + r_l*vx_l) -
     .               q_up * (r_r*vx_r - r_l*vx_l) + anx * p_12
      f(i,3)       = q_cd * (r_r*vy_r + r_l*vy_l) -
     .               q_up * (r_r*vy_r - r_l*vy_l) + any * p_12
      f(i,4)       = q_cd * (r_r*vz_r + r_l*vz_l) -
     .               q_up * (r_r*vz_r - r_l*vz_l) + anz * p_12
      f(i,5)       = q_cd * (r_r*h_r  + r_l*h_l)  -
     .               q_up * (r_r*h_r  - r_l*h_l)
c
c     need to check for moving grids (gs related to at(i))
c
c     f(i,5)       = q_cd * (r_r*h_r  + r_l*h_l)  -
c    .               q_up * (r_r*h_r  - r_l*h_l)  + gs * area(i) * p_12
c
c     ------------------------------------------------------------------
c     additional terms (maps+ scheme)
c     ------------------------------------------------------------------
c
c     rho   = ccmin(r_l, r_r)
c     h     = ccmin(h_l, h_r)
c
      p_d   = p_r - p_l
      q_d   = q_r - q_l
      rho   = 0.5*( r_l+ r_r)
      u     = 0.5*(vx_l+vx_r)
      v     = 0.5*(vy_l+vy_r)
      w     = 0.5*(vz_l+vz_r)
      h     = 0.5*( h_l+ h_r)
      c_min = ccmin(c_l, c_r)
      c_max = ccmax(c_l, c_r)
      q_min = m_lo*c_min
      q_max = m_up*c_max
      qq    = u*u + v*v + w*w
      cc    = c_max*c_max
c
c     ------------------------------------------------------------------
c     rescale speed of sound for preconditioning
c     ------------------------------------------------------------------
c
c     qs     = q_max
c     qref2  = ccmin( ccmax( upc*qq, epslocm2*cc), cc)
c     alp    = (1.-qref2/cc)/2.0
c     q_maxp = q_max*(1.-alp)
c     q_minp = q_maxp
c     c_maxp = sqrt(qs*qs * alp*alp + qref2)
c     c_minp = c_maxp
c 
c     if (no preconditioning) then (var, ...) else (varp, ...) */
c     (globals.precond.preconditioning == 0) ? 
c     (q_max = q_max , q_min = q_min , c_max = c_max , c_min = c_min ) :
c     (q_max = q_maxp, q_min = q_minp, c_max = c_maxp, c_min = c_minp) 
c
      rc    = 1.0/c_max
      m_max = q_max/c_min
      m_0   = ccmincr(m_max, 1.0)
c     m_m1  = 2.0*ccmaxcr(0.5-m_0, 0.0) 
      m_m1  = ccmaxcr(1.0-m_0, 0.0)
c
      p_scal = 0.5 * rc  * m_m1         * p_d
      q_scal = 0.5 * rho * m_m1 * c_min * q_d
c
      f(i,1) = f(i,1) - area(i) *     p_scal
      f(i,2) = f(i,2) - area(i) * u * p_scal - anx * q_scal
      f(i,3) = f(i,3) - area(i) * v * p_scal - any * q_scal
      f(i,4) = f(i,4) - area(i) * w * p_scal - anz * q_scal
      f(i,5) = f(i,5) - area(i) * h * p_scal
c
 1000 continue
c
      return
      end
