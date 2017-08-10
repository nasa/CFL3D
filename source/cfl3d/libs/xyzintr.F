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
      subroutine xyzintr(jdimf,kdimf,idimf,jdimc,kdimc,idimc,xf,yf,
     .                   zf,xc,yc,zc)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Interpolate coarse grid x,y,z to fine grid for mesh
c              sequencing in deforming mesh cases
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xf(jdimf,kdimf,idimf),yf(jdimf,kdimf,idimf),
     .          zf(jdimf,kdimf,idimf)
      dimension xc(jdimc,kdimc,idimc),yc(jdimc,kdimc,idimc),
     .          zc(jdimc,kdimc,idimc)
c
      common /twod/ i2d
c
      if (i2d.eq.0) then
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do kf=1,kdimf,2
               kc = (kf+1)/2
               do if=1,idimf,2
                  ic = (if+1)/2
                  xf(jf,kf,if) = xc(jc,kc,ic)
                  yf(jf,kf,if) = yc(jc,kc,ic)
                  zf(jf,kf,if) = zc(jc,kc,ic)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do if=1,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf-1,kf,if)+xf(jf+1,kf,if))
                  yf(jf,kf,if) = 0.5*(yf(jf-1,kf,if)+yf(jf+1,kf,if))
                  zf(jf,kf,if) = 0.5*(zf(jf-1,kf,if)+zf(jf+1,kf,if))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do if=1,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf,kf-1,if)+xf(jf,kf+1,if))
                  yf(jf,kf,if) = 0.5*(yf(jf,kf-1,if)+yf(jf,kf+1,if))
                  zf(jf,kf,if) = 0.5*(zf(jf,kf-1,if)+zf(jf,kf+1,if))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=1,kdimf
               do if=2,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf,kf,if-1)+xf(jf,kf,if+1))
                  yf(jf,kf,if) = 0.5*(yf(jf,kf,if-1)+yf(jf,kf,if+1))
                  zf(jf,kf,if) = 0.5*(zf(jf,kf,if-1)+zf(jf,kf,if+1))
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
               do i=1,2
                  xf(jf,kf,i) = xc(jc,kc,i)
                  yf(jf,kf,i) = yc(jc,kc,i)
                  zf(jf,kf,i) = zc(jc,kc,i)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do i=1,2
                  xf(jf,kf,i) = 0.5*(xf(jf-1,kf,i)+xf(jf+1,kf,i))
                  yf(jf,kf,i) = 0.5*(yf(jf-1,kf,i)+yf(jf+1,kf,i))
                  zf(jf,kf,i) = 0.5*(zf(jf-1,kf,i)+zf(jf+1,kf,i))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do i=1,2
                  xf(jf,kf,i) = 0.5*(xf(jf,kf-1,i)+xf(jf,kf+1,i))
                  yf(jf,kf,i) = 0.5*(yf(jf,kf-1,i)+yf(jf,kf+1,i))
                  zf(jf,kf,i) = 0.5*(zf(jf,kf-1,i)+zf(jf,kf+1,i))
               end do
            end do
         end do
c
      end if
c
      return
      end
