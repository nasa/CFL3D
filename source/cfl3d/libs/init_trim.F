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
      subroutine init_trim
c
c     $Id$
c
c***********************************************************************
c     Purpose: initialize trim data.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst 
      common /trim1/ dcl(5000),ddclda(5000),ddcmda(5000),r33,r44,ittr
     .              ,dcm(5000),dd(5000),da(5000),a11,a12,a22,r11,r22
c
c     The following relaxation coefficients and limits are hardwired:
c
      dtrmsmx =  0.50
      dtrmsmn = -0.50
      dmtrmn  =  0.
      dmtrmnm =  0.
      dlcln   =  0.
      dlclnm  =  0.
      ddtmx   =  0.0005
      dalfmx  =  0.0005
      trtol   =  0.0000002
      itrmt   =  0
      epstr   =  1.e-7
      a11     = 0.
      a12     = 0.
      a22     = 0.
      r11     = 0.
      r22     = 0.
      r33     = 0.
      r44     = 0.
c
      dtr     = tp(1,1)*tp(2,2)-tp(1,2)*tp(2,1)
      dtr     = dtr + ccsignrc(1e-30,dtr)
      fp(1,1) = tp(2,2)/dtr
      fp(1,2) =-tp(1,2)/dtr
      fp(2,1) =-tp(2,1)/dtr
      fp(2,2) = tp(1,1)/dtr
c
      return
      end
