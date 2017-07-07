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
      subroutine shortinp(nb1,ne1,nb2,ne2,idm,jdm,kdm,iface)
c
c     $Id$
c
c     converts cfl3d "shortcut" input to full range of index values
c
      if (nb1.eq.0) nb1 = 1
      if (nb2.eq.0) nb2 = 1
      if (ne1.eq.0) then
         if (iface.eq.1 .or. iface.eq.2) ne1 = jdm
         if (iface.eq.3 .or. iface.eq.4) ne1 = kdm
         if (iface.eq.5 .or. iface.eq.6) ne1 = idm
      end if
      if (ne2.eq.0) then
         if (iface.eq.1 .or. iface.eq.2) ne2 = kdm
         if (iface.eq.3 .or. iface.eq.4) ne2 = idm
         if (iface.eq.5 .or. iface.eq.6) ne2 = jdm
      end if
c 
      return
      end 
