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
      subroutine blnkfr(nbl,iibg,kkbg,jjbg,ibpntsg,lbg,iitot,blank,
     .                  jdim,kdim,idim,maxbl,blnkval)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Undo/redo blanking at fringe points (but not hole points)
c     blnkval = 1...set blank=1 at fringe points for plotting
c               0...set blank=0 (normal setting for fringe points)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension iibg(iitot),kkbg(iitot),jjbg(iitot),ibpntsg(maxbl,4),
     .          lbg(maxbl),blank(jdim,kdim,idim)
      real blnkval
c
      lsta = lbg(nbl)
      lend = lsta-1
      if (ibpntsg(nbl,1).gt.0) then
         lend = lsta+ibpntsg(nbl,1)-1
         do l=lsta,lend
            blank(jjbg(l),kkbg(l),iibg(l)) = blnkval
         end do
      end if
c
      return
      end
