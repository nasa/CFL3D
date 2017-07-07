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
      subroutine unld_qc(qc,qtemp,jc,kc,ic,is,ie,js,je,ks,ke,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Install the qtemp array used for message passing into
c     the appropriate section of the qc array. Only the cell-center
c     locations are unloaded, not the fill in locations on the idim,
c     jdim and kdim planes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,kc,ic,ldim),qtemp(je-js+1,ke-ks+1,ie-is+1,ldim)
c
      do l =1,ldim
         jj = js - 1
         do j=1,je-js
            jj = jj + 1
            kk = ks - 1
            do k=1,ke-ks
               kk = kk + 1
               ii = is - 1
               do i=1,ie-is
                  ii = ii + 1
                  qc(jj,kk,ii,l) = qtemp(j,k,i,l)
               end do
            end do
         end do
      end do
      return
      end
