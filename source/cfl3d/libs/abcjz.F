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
      subroutine abcjz(i,npl,jdim,kdim,idim,aj,bj,cj,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Zero out the LHS matrix element with the help of the 
c     blank array.  For a point with blank=0, all elements of matrices
c     a and c become zero. Only diagonal elements matrix b is changed
c     to 1.0 for J-implicit ; K-sweep.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension aj(kdim-1,npl,jdim,5,5),bj(kdim-1,npl,jdim,5,5),
     .          cj(kdim-1,npl,jdim,5,5),blank(jdim,kdim,idim)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      do 2750 m=1,5 
      do 2750 n=1,5 
      do 2750 j=1,jdim1
      do 2750 lpl=1,npl 
      ii = i+lpl-1
cdir$ ivdep
      do 2765 k=1,kdim1 
      aj(k,lpl,j,m,n) = aj(k,lpl,j,m,n)*blank(j,k,ii)
      cj(k,lpl,j,m,n) = cj(k,lpl,j,m,n)*blank(j,k,ii)
 2765 continue
 2750 continue
c
      do 2850 m=1,5 
      do 2850 j=1,jdim1
      do 2850 lpl=1,npl 
      ii = i+lpl-1
      do 2865 k=1,kdim1 
      bj(k,lpl,j,m,m) =  (bj(k,lpl,j,m,m)*blank(j,k,ii))
     .                  +(1.0-blank(j,k,ii))
 2865 continue
 2850 continue
      return
      end 
