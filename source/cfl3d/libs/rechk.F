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
      subroutine rechk(mdim1,ndim1,jimage,kimage,msub1,jmax1,kmax1,l,
     .                 x1,y1,z1,xie1,xie2,eta1,eta2,nou,bou,nbuf,
     .                 ibufdim,myid,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Check for branch cuts 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      integer xie1,xie2,eta1,eta2
c
      dimension nou(nbuf)
      dimension mblk2nd(maxbl)
      common /tol/ epsc,epsc0,eps,epscoll
      common /sklt1/isklt1
c
      dimension x1(mdim1,ndim1),y1(mdim1,ndim1),z1(mdim1,ndim1)
      dimension jimage(msub1,mdim1,ndim1),kimage(msub1,mdim1,ndim1)
c      
c jimage(l,j,k)...array for "from" block l (local block numbering)
c                   which contains the j-index for the image of
c                   the point j,k
c kimage(l,j,k)...array for "from" block l which contains the k-index for 
c                   the image of the point j,k
c
c   The indicies of the image points found in this subroutine anticipate
c   that 1) the (sub)grid(s) is(are) eventually to be expanded and 
c        2) the indicies need to point to the lower left hand corner 
c           of the image cell 
c
c     default to no branch cuts 
c
      do 5 j=1,jmax1+1
      do 5 k=1,kmax1+1
      jimage(l,j,k) = j
    5 kimage(l,j,k) = k  
c
c     check for C-grid type branch cut on k=1 and k=kmax1 boundaries 
c
c     branch cuts occur only on block edges, so if search range does not
c     include block edges, don't bother to check for cuts
c
      if (eta1.gt.1 .and. eta2.lt.kmax1) go to 11
c
      do 10 k=1,kmax1,kmax1-1
      kcheck = 0
      do 1711 js=xie1,xie2-1
      do 1712 jq=js+1,xie2
      delxyz = ccabs( x1(js,k)-x1(jq,k) ) 
     .       + ccabs( y1(js,k)-y1(jq,k) )
     .       + ccabs( z1(js,k)-z1(jq,k) )
      if (real(delxyz).lt.real(eps)) then
         kcheck = kcheck+1
         if (k.eq.1) then
            jimage(l,js+1,1) = jq
            kimage(l,js+1,1) = 2
            jimage(l,jq,1)   = js+1
            kimage(l,jq,1)   = 2
         end if
         if (k.eq.kmax1) then
            jimage(l,js+1,kmax1+1) = jq
            kimage(l,js+1,kmax1+1) = kmax1
            jimage(l,jq,kmax1+1)   = js+1
            kimage(l,jq,kmax1+1)   = kmax1
         end if
      end if
 1712 continue
 1711 continue
      if (kcheck.eq.0) go to 10 
      if (k.eq.1) then
c
c     correction for border points between start/end 
c     of branch cuts on the boundary
c
      do 1811 j=xie1+1,xie2-1
      if (jimage(l,j+1,1).eq.j+1) then
         jimage(l,jimage(l,j,1),1) = jimage(l,j,1)
         kimage(l,jimage(l,j,1),1) = 1
         jimage(l,j,1)             = j
         kimage(l,j,1)             = 1
      go to 1809
      end if
 1811 continue
 1809 continue
c
c     image point at ends of boundary of extended grid
c
      if (jimage(l,xie1+1,1).eq.xie2) then
         jimage(l,xie1,1)   = xie2+1
         kimage(l,xie1,1)   = 2
         jimage(l,xie2+1,1) = xie1
         kimage(l,xie2+1,1) = 2
      end if
      end if
c
      if (k.eq.kmax1) then
c
c     correction for border points between start/end 
c     of branch cuts on the boundary
c
      do 1812 j=xie1+1,xie2-1
      if (jimage(l,j+1,kmax1+1).eq.j+1) then
         jimage(l,jimage(l,j,kmax1+1),kmax1+1) = jimage(l,j,kmax1+1)
         kimage(l,jimage(l,j,kmax1+1),kmax1+1) = kmax1+1
         jimage(l,j,kmax1+1)                   = j
         kimage(l,j,kmax1+1)                   = kmax1+1
         go to 1808
      end if
 1812 continue
 1808 continue
c
c     image point at left boundary of extended grid
c
      if (jimage(l,xie1+1,kmax1+1).eq.xie2) then
         jimage(l,xie1,kmax1+1)   = xie2+1
         kimage(l,xie1,kmax1+1)   = kmax1
         jimage(l,xie2+1,kmax1+1) = xie1
         kimage(l,xie2+1,kmax1+1) = kmax1
      end if
      end if
   10 continue
c
   11 continue
c
c     check for C-grid type branch cut on j=1 and j=jmax1 boundaries
c
c     branch cuts occur only on block edges, so if search range does not
c     include block edges, don't bother to check for cuts
c
      if (xie1.gt.1 .and. xie2.lt.jmax1) go to 21
c
      do 20 j=1,jmax1,jmax1-1
      jcheck = 0
      do 2712 ks=eta1,eta2-1
      do 2711 kq=ks+1,eta2
      delxyz = ccabs( x1(j,ks)-x1(j,kq) ) 
     .       + ccabs( y1(j,ks)-y1(j,kq) )
     .       + ccabs( z1(j,ks)-z1(j,kq) )
      if (real(delxyz).lt.real(eps)) then
      jcheck = jcheck+1
        if (j.eq.1) then
           jimage(l,1,ks+1) = 2
           kimage(l,1,ks+1) = kq
           jimage(l,1,kq)   = 2
           kimage(l,1,kq)   = ks+1
        end if
        if (j.eq.jmax1) then
           jimage(l,jmax1+1,ks+1) = jmax1
           kimage(l,jmax1+1,ks+1) = kq
           jimage(l,jmax1+1,kq)   = jmax1
           kimage(l,jmax1+1,kq)   = ks+1
        end if
      go to 2712
      end if
 2711 continue
 2712 continue
      if (jcheck.eq.0) go to 20
      if (j.eq.1) then
c
c     correction for border points between start/end 
c     of branch cuts on the boundary
c
      do 2811 k=eta1+1,eta2-1
      if  (kimage(l,1,k+1).eq.k+1) then
          jimage(l,1,kimage(l,1,k)) = 1
          kimage(l,1,kimage(l,1,k)) = kimage(l,1,k)
          jimage(l,1,k)             = 1
          kimage(l,1,k)             = k
          go to 2809
      end if
 2811 continue
 2809 continue
c
c     image points at ends of boundary of extended grid
c
      if (kimage(l,1,eta1+1).eq.eta2) then
         jimage(l,1,eta1)   = 2
         kimage(l,1,eta1)   = eta2+1
         jimage(l,1,eta2+1) = 2
         kimage(l,1,eta2+1) = eta1
      end if
      end if
c
      if (j.eq.jmax1) then
c
c     correction for border points between start/end 
c     of branch cuts on the boundary
c
      do 2812 k=eta1+1,eta2-1
      if (kimage(l,jmax1+1,k+1).eq.k+1) then
         jimage(l,jmax1+1,kimage(l,jmax1+1,k)) = jmax1+1
         kimage(l,jmax1+1,kimage(l,jmax1+1,k)) = kimage(l,jmax1+1,k)
         jimage(l,jmax1+1,k)                   = jmax1+1
         kimage(l,jmax1+1,k)                   = k
         go to 2808
      end if
 2812 continue
 2808 continue
c
c     image point at bottom boundary of extended grid
c
      if (kimage(l,jmax1+1,eta1+1).eq.eta2) then
         jimage(l,jmax1+1,eta1)   = jmax1
         kimage(l,jmax1+1,eta1)   = eta2+1
         jimage(l,jmax1+1,eta2+1) = jmax1
         kimage(l,jmax1+1,eta2+1) = eta1
      end if
      end if
   20 continue
c
   21 continue
c
c     check for O-grid type branch cut on k=1 and k=kmax1 boundaries
c
      jcount = 0
c
c     branch cuts occur only on block edges, so if search range does not
c     include block edges, don't bother to check for cuts
c
      if (eta1.gt.1 .and. eta2.lt.kmax1) go to 31
c
      do 30 j=xie1,xie2
      delxyz = ccabs( x1(j,1)-x1(j,kmax1) ) 
     .       + ccabs( y1(j,1)-y1(j,kmax1) )
     .       + ccabs( z1(j,1)-z1(j,kmax1) )
      if (real(delxyz).lt.real(eps)) jcount = jcount+1
   30 continue
      if (jcount.gt.1.and.jcount.ne.xie2-xie1+1) then
      if (isklt1.gt.1) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),'('' WARNING...O-type branch cut does'',
     .           '' not extend over entire k=constant  boundary, as'',
     .           '' assumed'')')
      end if
      end if
      if (jcount.gt.1) then
         do 35 j=xie1,xie2+1
         jimage(l,j,1)       = j
         kimage(l,j,1)       = kmax1
         jimage(l,j,kmax1+1) = j
         kimage(l,j,kmax1+1) = 2
   35    continue
      end if  
c
   31 continue    
c     
c     check for O-grid type branch cut on j=1 and j=jmax1 boundaries
c
      kcount = 0
c
c     branch cuts occur only on block edges, so if search range does not
c     include block edges, don't bother to check for cuts
c
      if (xie1.gt.1 .and. xie2.lt.jmax1) go to 41
c
      do 40 k=eta1,eta2
      delxyz = ccabs( x1(1,k)-x1(jmax1,k) ) 
     .       + ccabs( y1(1,k)-y1(jmax1,k) )
     .       + ccabs( z1(1,k)-z1(jmax1,k) )
      if (real(delxyz).lt.real(eps)) kcount = kcount+1
   40 continue
      if (kcount.gt.1.and.kcount.ne.eta2-eta1+1) then
      if (isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),'('' WARNING...O-type branch cut does'',
     .          '' not extend over entire k=constant  boundary, as'',
     .          '' assumed'')')
      end if
      end if
      if (kcount.gt.1) then
         do 45 k=eta1,eta2+1
         jimage(l,1,k)       = jmax1
         kimage(l,1,k)       = k
         jimage(l,jmax1+1,k) = 2
         kimage(l,jmax1+1,k) = k
   45    continue
      end if      
c
   41 continue
c
c     print out branch cut topology
c
      do 50 j=xie1,xie2+1,xie2+1-xie1
      iflag = 0
      jpr   = j
      if (j.eq.jmax1+1) jpr = jmax1
      do 51 k=eta1+1,eta2
      if (jimage(l,j,k).ne.j .or. kimage(l,j,k).ne.k) iflag = 1
   51 continue
      iflago = 0
      do 52 k=eta1+1,eta2
      if (jimage(l,j,k).eq.j .or. kimage(l,j,k).ne.k) iflago = 1
   52 continue
      if (iflag.ne.0 .and. iflago.eq.0 .and. isklt1.gt.0) then 
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),101)jpr
      end if
  101 format(' ','     j=',i3,' is an O-type (periodic)',
     .       ' branch cut boundary')
      if (iflag.ne.0 .and. iflago.ne.0 .and. isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),102) jpr
      end if
  102 format(' ','     j=',i3,' is a C-type branch cut boundary')
   50 continue
      do 60 k=eta1,eta2+1,eta2+1-eta1
      kpr = k
      if (k.eq.kmax1+1) kpr = kmax1
      iflag = 0
      do 61 j=xie1+1,xie2
      if (jimage(l,j,k).ne.j .or. kimage(l,j,k).ne.k) iflag = 1
   61 continue
      iflago = 0
      do 62 j=xie1+1,xie2
      if (jimage(l,j,k).ne.j .or. kimage(l,j,k).eq.k) iflago = 1
   62 continue
      if (iflag.ne.0 .and. iflago.eq.0 .and. isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),201) kpr
      end if
  201 format(' ','     k=',i3,' is an O-type (periodic)',
     .       ' branch cut boundary')
      if (iflag.ne.0 .and. iflago.ne.0 .and. isklt1.gt.0) then 
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),202) kpr
      end if
  202 format(' ','     k=',i3,' is a C-type branch cut boundary')
   60 continue
      return
      end 
