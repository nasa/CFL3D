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
#if defined WKSTN_OFF
      function ismin_wkstn(n,x,incx)
#else
      function ismin(n,x,incx)
#endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To return the index corresponding to the minimum value 
c     in an arrary; a replacement for the standard Cray function
c***********************************************************************
c
      dimension x(1+(n-1)*incx)
c
      ismin = 1
      xmin = x(1)
c
      do 100 nn=2,n,incx
      if(x(nn).lt.xmin) then
        ismin = nn
        xmin  = x(nn)
      end if
100   continue
c
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      ismin_wkstn=-999
#endif
c
      return
      end 
#if defined WKSTN_OFF
      function ismax_wkstn(n,x,incx)
#else
      function ismax(n,x,incx)
#endif
c***********************************************************************
c     Purpose:  To return the index corresponding to the maximum value
c     in an arrary; a replacement for the standard Cray function
c***********************************************************************
c
c
      dimension x(1+(n-1)*incx)
c
      ismax = 1
      xmax  = x(1)
c
      do 100 nn=2,n,incx
      if(x(nn).gt.xmax) then
        ismax= nn
        xmax  = x(nn)
      end if
100   continue
c
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      ismax_wkstn=-999
#endif
c
      return
      end
