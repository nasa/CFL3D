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
      subroutine collapse(jdim,kdim,jmaxo,kmaxo,xo,yo,zo,
     .                    nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Check for collapsed points in grid, and epand any
c     collapsed lines detected. The threshold for adjacent points
c     being considered as collapsed is eps; 1)for lines on which all
c     segments are not collapsed, collapsed segments between two 
c     adjacent points are expanded to fact*eps in size, 2) for lines
c     on which all segments are collapsed, a polar singularity is
c     assumed, and the expansion factor is epss2*r, where r is the
c     local distance to the next (uncollapsed) grid line.
c
c     specifically excluded is the case where the entire block face
c     is collapsed - in that case there are no local uncollapsed
c     grid lines to indicate in which direction the segments should
c     be uncollapsed.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /tol/ epsc,epsc0,epsreen,eps
      common /sklt1/isklt1
c
      dimension xo(jdim,kdim),yo(jdim,kdim),zo(jdim,kdim)
c
      fact  = 1.e5
      epss1 = fact*eps
c
      epss2 = 1.e-01
c
c******************************
c     check along lines k=const
c******************************
c
      do 111 k=1,kmaxo
      kc     = k
      j1     = 1
      ncount = 0
  110 j2     = j1 + 1 + ncount
      if (abs(real(xo(j2,k)-xo(j1,k))).lt.real(eps) .and.
     .    abs(real(yo(j2,k)-yo(j1,k))).lt.real(eps) .and.
     .    abs(real(zo(j2,k)-zo(j1,k))).lt.real(eps)) then
          ncount = ncount+1
          if (j1.eq.1.and.j2.eq.jmaxo) go to 109
          go to 110
      end if
      if (ncount.gt.0) then
         j2 = j2-1
         if (isklt1 .gt. 0) then
            nou(4) = min(nou(4)+1,ibufdim)
            write(bou(nou(4),4),101) kc,j1,j2
         end if
  101    format(' ',5x,'collapsed boundary on k =',i3,2x,
     .          'between j =',i3,2x,'and j =',i3)
c
c        1) collapsed segments only over part of the k=const line
c
         if (j1.gt.1) then
            sx0 = xo(j1,k)-xo(j1-1,k)
            sy0 = yo(j1,k)-yo(j1-1,k)
            sz0 = zo(j1,k)-zo(j1-1,k)
         else
            sx0 = xo(j2+1,k)-xo(j2,k)
            sy0 = yo(j2+1,k)-yo(j2,k)
            sz0 = zo(j2+1,k)-zo(j2,k)
         end if
c
         xoo = xo(j1,k)
         yoo = yo(j1,k)
         zoo = zo(j1,k)
c
         if (j1.gt.1 .and. j2 .lt.jmaxo) then
            jc  = j1 + ncount/2
         else if (j1.eq.1) then
            jc  = j1
         else
            jc  = j2
         end if
c
         do 88 jj=j1,j2     
         fact1  = float(jj-jc)*epss1  
         xo(jj,k) = xoo+fact1*sx0
         yo(jj,k) = yoo+fact1*sy0
         zo(jj,k) = zoo+fact1*sz0
   88    continue
c
      end if 
c
      j1     = j1 + ncount + 1
      ncount = 0 
      if (j1.lt.jmaxo-1) go to 110
c
      go to 111
c
c     2) collapsed segments over the entire k=const line, but k+1 
c        line is not collapsed (or k-1 for k=kmaxo) - polar singularity
c
  109 continue
      if (isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),101)kc,j1,j2
      end if
      m     = 1
      if (k.eq.kmaxo) m = -1
      do 50 j=1,jmaxo
      sx      = xo(j,k+m)-xo(j,k)
      sy      = yo(j,k+m)-yo(j,k)
      sz      = zo(j,k+m)-zo(j,k)
      xo(j,k) = xo(j,k)+sx*epss2
      yo(j,k) = yo(j,k)+sy*epss2
      zo(j,k) = zo(j,k)+sz*epss2
   50 continue
c
  111 continue  
c
c******************************
c     check along lines j=const
c******************************
c
      do 1111 j=1,jmaxo
      jc     =   j
      k1     =   1
      ncount =   0
 1101 k2     =   k1 + 1 + ncount
      if (abs(real(xo(j,k2)-xo(j,k1))).lt.real(eps) .and.
     .    abs(real(yo(j,k2)-yo(j,k1))).lt.real(eps) .and.
     .    abs(real(zo(j,k2)-zo(j,k1))).lt.real(eps)) then
          ncount = ncount+1
          if (k1.eq.1.and.k2.eq.kmaxo) go to 1009
          go to 1101
      end if
      if (ncount.gt.0) then
         k2 = k2-1
         if (isklt1.gt.0) then
            nou(4) = min(nou(4)+1,ibufdim)
            write(bou(nou(4),4),201) jc,k1,k2
         end if
  201    format(' ',5x,'collapsed boundary on j =',i3,2x,
     .          'between k =',i3,2x,'and k =',i3)
c
c        1) collapsed segments only over part of the j=const line
c
         if (k1.gt.1) then
            sx0 = xo(j,k1)-xo(j,k1-1)
            sy0 = yo(j,k1)-yo(j,k1-1)
            sz0 = zo(j,k1)-zo(j,k1-1)
         else 
            sx0 = xo(j,k2+1)-xo(j,k2)
            sy0 = yo(j,k2+1)-yo(j,k2)
            sz0 = zo(j,k2+1)-zo(j,k2)
         end if
c
         xoo = xo(j,k1)
         yoo = yo(j,k1)
         zoo = zo(j,k1)
c
         if (k1.gt.1 .and. k2 .lt.kmaxo) then
            kc  = k1 + ncount/2
         else if (k1.eq.1) then
            kc  = k1
         else
            kc  = k2
         end if

         do 881 kk=k1,k2
         fact1  = float(kk-kc)*epss1
         xo(j,kk) = xoo+fact1*sx0
         yo(j,kk) = yoo+fact1*sy0
         zo(j,kk) = zoo+fact1*sz0
  881    continue
c
      end if 
c
      k1     = k1 + ncount + 1
      ncount = 0 
      if (k1.lt.kmaxo-1) go to 1101
c
      go to 1111
c
c     2) collapsed segments over the entire j=const line, but j+1
c        line is not collapsed (or j-1 for j=jmaxo) - polar singularity
c
 1009 continue
      if (isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),201)jc,k1,k2
      end if
      m     = 1
      if (j.eq.jmaxo) m = -1
      do 60 k=1,kmaxo
      sx      = xo(j+m,k)-xo(j,k)
      sy      = yo(j+m,k)-yo(j,k)
      sz      = zo(j+m,k)-zo(j,k)
      xo(j,k) = xo(j,k)+sx*epss2
      yo(j,k) = yo(j,k)+sy*epss2
      zo(j,k) = zo(j,k)+sz*epss2
   60 continue
c
 1111 continue
c
      return
      end
