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
      subroutine ld_qc(qc,qtemp,jc,kc,ic,is,ie,js,je,ks,ke,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Install the embedded subset of the qc array into qtemp
c     for massage passing.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,kc,ic,ldim),qtemp(je-js+1,ke-ks+1,ie-is+1,ldim)
c
      do l =1,ldim
         jj = 0
         do j=js,je
            jj = jj + 1
            kk = 0
            do k=ks,ke
               kk = kk + 1
               ii = 0
               do i=is,ie
                  ii = ii + 1
                  qtemp(jj,kk,ii,l) = qc(j,k,i,l)
               end do
            end do
         end do
      end do
      return
      end
