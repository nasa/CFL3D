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
      subroutine delintr(jdimf,kdimf,idimf,jdimc,kdimc,idimc,deltjf,
     .                   deltkf,deltif,deltjc,deltkc,deltic)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Interpolate coarse grid deltas to fine grid for mesh
c              sequencing in deforming mesh cases
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension deltjf(kdimf,idimf,3,2),deltkf(jdimf,idimf,3,2),
     .          deltif(jdimf,kdimf,3,2)
      dimension deltjc(kdimc,idimc,3,2),deltkc(jdimc,idimc,3,2),
     .          deltic(jdimc,kdimc,3,2)
c
      common /twod/ i2d
c
      if (i2d.eq.0) then
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do kf=1,kdimf,2
               kc = (kf+1)/2
               do ll=1,2
                  deltif(jf,kf,1,ll) = deltic(jc,kc,1,ll)
                  deltif(jf,kf,2,ll) = deltic(jc,kc,2,ll)
                  deltif(jf,kf,3,ll) = deltic(jc,kc,3,ll)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do ll=1,2
                  deltif(jf,kf,1,ll) = 0.5*(deltif(jf-1,kf,1,ll)
     .                                     +deltif(jf+1,kf,1,ll))
                  deltif(jf,kf,2,ll) = 0.5*(deltif(jf-1,kf,2,ll)
     .                                     +deltif(jf+1,kf,2,ll))
                  deltif(jf,kf,3,ll) = 0.5*(deltif(jf-1,kf,3,ll)
     .                                     +deltif(jf+1,kf,3,ll))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do ll=1,2
                  deltif(jf,kf,1,ll) = 0.5*(deltif(jf,kf-1,1,ll)
     .                                     +deltif(jf,kf+1,1,ll))
                  deltif(jf,kf,2,ll) = 0.5*(deltif(jf,kf-1,2,ll)
     .                                     +deltif(jf,kf+1,2,ll))
                  deltif(jf,kf,3,ll) = 0.5*(deltif(jf,kf-1,3,ll)
     .                                     +deltif(jf,kf+1,3,ll))
               end do
            end do
         end do
c
         do kf=1,kdimf,2
            kc = (kf+1)/2
            do if=1,idimf,2
               ic = (if+1)/2
               do ll=1,2
                  deltjf(kf,if,1,ll) = deltjc(kc,ic,1,ll)
                  deltjf(kf,if,2,ll) = deltjc(kc,ic,2,ll)
                  deltjf(kf,if,3,ll) = deltjc(kc,ic,3,ll)
               end do
            end do
         end do
         do kf=2,kdimf,2
            do if=1,idimf,2
               do ll=1,2
                  deltjf(kf,if,1,ll) = 0.5*(deltjf(kf-1,if,1,ll)
     .                                     +deltjf(kf+1,if,1,ll))
                  deltjf(kf,if,2,ll) = 0.5*(deltjf(kf-1,if,2,ll)
     .                                     +deltjf(kf+1,if,2,ll))
                  deltjf(kf,if,3,ll) = 0.5*(deltjf(kf-1,if,3,ll)
     .                                     +deltjf(kf+1,if,3,ll))
               end do
            end do
         end do
         do kf=1,kdimf
            do if=2,idimf,2
               do ll=1,2
                  deltjf(kf,if,1,ll) = 0.5*(deltjf(kf,if-1,1,ll)
     .                                     +deltjf(kf,if+1,1,ll))
                  deltjf(kf,if,2,ll) = 0.5*(deltjf(kf,if-1,2,ll)
     .                                     +deltjf(kf,if+1,2,ll))
                  deltjf(kf,if,3,ll) = 0.5*(deltjf(kf,if-1,3,ll)
     .                                     +deltjf(kf,if+1,3,ll))
               end do
            end do
         end do
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do if=1,idimf,2
               ic = (if+1)/2
               do ll=1,2
                  deltkf(jf,if,1,ll) = deltkc(jc,ic,1,ll)
                  deltkf(jf,if,2,ll) = deltkc(jc,ic,2,ll)
                  deltkf(jf,if,3,ll) = deltkc(jc,ic,3,ll)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do if=1,idimf,2
               do ll=1,2
                  deltkf(jf,if,1,ll) = 0.5*(deltkf(jf-1,if,1,ll)
     .                                     +deltkf(jf+1,if,1,ll))
                  deltkf(jf,if,2,ll) = 0.5*(deltkf(jf-1,if,2,ll)
     .                                     +deltkf(jf+1,if,2,ll))
                  deltkf(jf,if,3,ll) = 0.5*(deltkf(jf-1,if,3,ll)
     .                                     +deltkf(jf+1,if,3,ll))
               end do
            end do
         end do
         do jf=1,jdimf
            do if=2,idimf,2
               do ll=1,2
                  deltkf(jf,if,1,ll) = 0.5*(deltkf(jf,if-1,1,ll)
     .                                     +deltkf(jf,if+1,1,ll))
                  deltkf(jf,if,2,ll) = 0.5*(deltkf(jf,if-1,2,ll)
     .                                     +deltkf(jf,if+1,2,ll))
                  deltkf(jf,if,3,ll) = 0.5*(deltkf(jf,if-1,3,ll)
     .                                     +deltkf(jf,if+1,3,ll))
               end do
            end do
         end do
c
      else
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do kf=1,kdimf,2
               kc = (kf+1)/2
               do ll=1,2
                  deltif(jf,kf,1,ll) = deltic(jc,kc,1,ll)
                  deltif(jf,kf,2,ll) = deltic(jc,kc,2,ll)
                  deltif(jf,kf,3,ll) = deltic(jc,kc,3,ll)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do ll=1,2
                  deltif(jf,kf,1,ll) = 0.5*(deltif(jf-1,kf,1,ll)
     .                                     +deltif(jf+1,kf,1,ll))
                  deltif(jf,kf,2,ll) = 0.5*(deltif(jf-1,kf,2,ll)
     .                                     +deltif(jf+1,kf,2,ll))
                  deltif(jf,kf,3,ll) = 0.5*(deltif(jf-1,kf,3,ll)
     .                                     +deltif(jf+1,kf,3,ll))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do ll=1,2
                  deltif(jf,kf,1,ll) = 0.5*(deltif(jf,kf-1,1,ll)
     .                                     +deltif(jf,kf+1,1,ll))
                  deltif(jf,kf,2,ll) = 0.5*(deltif(jf,kf-1,2,ll)
     .                                     +deltif(jf,kf+1,2,ll))
                  deltif(jf,kf,3,ll) = 0.5*(deltif(jf,kf-1,3,ll)
     .                                     +deltif(jf,kf+1,3,ll))
               end do
            end do
         end do
c
         do kf=1,kdimf,2
            kc = (kf+1)/2
            do i=1,2
               do ll=1,2
                  deltjf(kf,i,1,ll) = deltjc(kc,i,1,ll)
                  deltjf(kf,i,2,ll) = deltjc(kc,i,2,ll)
                  deltjf(kf,i,3,ll) = deltjc(kc,i,3,ll)
               end do
            end do
         end do
         do kf=2,kdimf,2
            do i=1,2
               do ll=1,2
                  deltjf(kf,i,1,ll) = 0.5*(deltjf(kf-1,i,1,ll)
     .                                    +deltjf(kf+1,i,1,ll))
                  deltjf(kf,i,2,ll) = 0.5*(deltjf(kf-1,i,2,ll)
     .                                    +deltjf(kf+1,i,2,ll))
                  deltjf(kf,i,3,ll) = 0.5*(deltjf(kf-1,i,3,ll)
     .                                    +deltjf(kf+1,i,3,ll))
               end do
            end do
         end do
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do i=1,2
               do ll=1,2
                  deltkf(jf,i,1,ll) = deltkc(jc,i,1,ll)
                  deltkf(jf,i,2,ll) = deltkc(jc,i,2,ll)
                  deltkf(jf,i,3,ll) = deltkc(jc,i,3,ll)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do i=1,2
               do ll=1,2
                  deltkf(jf,i,1,ll) = 0.5*(deltkf(jf-1,i,1,ll)
     .                                    +deltkf(jf+1,i,1,ll))
                  deltkf(jf,i,2,ll) = 0.5*(deltkf(jf-1,i,2,ll)
     .                                    +deltkf(jf+1,i,2,ll))
                  deltkf(jf,i,3,ll) = 0.5*(deltkf(jf-1,i,3,ll)
     .                                    +deltkf(jf+1,i,3,ll))
               end do
            end do
         end do
c
      end if
c
      return
      end
