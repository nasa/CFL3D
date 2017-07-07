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
      subroutine ingrd(nf,ib,x,y,z,ni,nj,nk,ip3dgrd)
c
c     $Id$
c
      dimension x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
c
      if (ip3dgrd.gt.0) then      
         if (ib.eq.0) then
            read(nf,*)           (((x(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((y(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((z(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk)
         else
            read(nf)             (((x(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((y(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((z(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk)
         end if
      else
         if (ib.eq.0) then
            read(nf,*)           (((x(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((y(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((z(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni)
         else
            read(nf)             (((x(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((y(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((z(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni)
         end if
      end if
c
      return
      end
