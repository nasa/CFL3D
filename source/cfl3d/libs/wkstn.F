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
      integer function isrcheq_wkstn(n,x,incx,target)
#else
      integer function isrcheq(n,x,incx,target)
#endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To find the first occurance in the array x that is equal
c     to target; a replacement for the standard Cray function
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(*)
      j=1
      isrcheq=0
      if(n.le.0) return
      if(incx.lt.0) j=1-(n-1)*incx
      do 100 i=1,n
      if(real(x(j)).eq.real(target)) go to 200
      j=j+incx
  100 continue
  200 isrcheq=i
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      isrcheq_wkstn=-999
#endif
      return
      end
