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
      subroutine modread(idim,jdim,kdim,n,nbl,iunit,jbcinfo,kbcinfo,
     .                   ibcinfo,nbcj0,nbcjdim,nbck0,nbckdim,
     .                   nbci0,nbcidim,maxbl,maxseg,nmds,xmdj,
     .                   xmdk,xmdi,iaes,maxaes)
c
c     $Id$
c
c***********************************************************************
c     Purpose: read modal deflections from a file
c***********************************************************************
c
#ifdef CMPLX
      implicit complex(a-h,o-z)
c
      integer stats
c
      real, dimension(:,:,:), allocatable :: xmdir
      real, dimension(:,:,:), allocatable :: xmdjr
      real, dimension(:,:,:), allocatable :: xmdkr
c
#else
c
      integer stats
c
      real, dimension(:,:,:), allocatable :: xmdir1 
      real, dimension(:,:,:), allocatable :: xmdjr1 
      real, dimension(:,:,:), allocatable :: xmdkr1 
#endif
c
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension xmdj(kdim,idim,6,nmds,maxaes),
     .          xmdk(jdim,idim,6,nmds,maxaes),
     .          xmdi(jdim,kdim,6,nmds,maxaes)
      common /elastic/ ndefrm,naesrf
      common /igrdtyp/ ip3dgrd,ialph

c

#ifdef CMPLX
c     allocate real arrays in which to read modal data
c
      memuse = 0
      allocate( xmdir(jdim,kdim,6), stat=stats )
      call umalloc(jdim*kdim*6,0,'xmdir',memuse,stats)
      allocate( xmdjr(kdim,idim,6), stat=stats )
      call umalloc(kdim*idim*6,0,'xmdjr',memuse,stats)
      allocate( xmdkr(jdim,idim,6), stat=stats )
      call umalloc(jdim*idim*6,0,'xmdkr',memuse,stats)
#else
c     allocate real arrays in which to read modal data
c
      memuse = 0
      allocate( xmdir1(jdim,kdim,6), stat=stats )
      call umalloc(jdim*kdim*6,0,'xmdir1',memuse,stats)
      allocate( xmdjr1(kdim,idim,6), stat=stats )
      call umalloc(kdim*idim*6,0,'xmdjr1',memuse,stats)
      allocate( xmdkr1(jdim,idim,6), stat=stats )
      call umalloc(jdim*idim*6,0,'xmdkr1',memuse,stats)
#endif
c
c     read in modal deflections for j=const surfaces
c
      do 110 m = 1,2
c
      if (m.eq.1) then
         j    = 1
         l    = 0
         nseg = nbcj0(nbl)
      else
         j    = jdim
         l    = 3
         nseg = nbcjdim(nbl)
      end if
c
      do 100 ns = 1,nseg
c
      lbc = iabs(jbcinfo(nbl,ns,1,m))
c
      iflag = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)*
     .        (lbc-2024)*(lbc-2034)*(lbc-2016)
c
      if (iflag.ne.0) goto 100
c
      ist = jbcinfo(nbl,ns,2,m)
      ifn = jbcinfo(nbl,ns,3,m)
      kst = jbcinfo(nbl,ns,4,m)
      kfn = jbcinfo(nbl,ns,5,m)
c
      if(ialph.eq.0) then
#ifdef CMPLX
        read(iunit,*)((xmdjr(k,i,l+1),xmdjr(k,i,l+2),
     .                 xmdjr(k,i,l+3),k=kst,kfn),i=ist,ifn)
        do i=ist,ifn
           do k=kst,kfn
              xmdj(k,i,l+1,n,iaes) = xmdjr(k,i,l+1)
              xmdj(k,i,l+2,n,iaes) = xmdjr(k,i,l+2)
              xmdj(k,i,l+3,n,iaes) = xmdjr(k,i,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdj(k,i,l+1,n,iaes),xmdj(k,i,l+2,n,iaes),
     .                 xmdj(k,i,l+3,n,iaes),k=kst,kfn),i=ist,ifn)
#endif
      else
#ifdef CMPLX
        read(iunit,*)((xmdjr(k,i,l+1),xmdjr(k,i,l+3),
     .                 xmdjr(k,i,l+2),k=kst,kfn),i=ist,ifn)
        do i=ist,ifn
           do k=kst,kfn
              xmdj(k,i,l+1,n,iaes) = xmdjr(k,i,l+1)
              xmdj(k,i,l+2,n,iaes) =-xmdjr(k,i,l+2)
              xmdj(k,i,l+3,n,iaes) = xmdjr(k,i,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdjr1(k,i,l+1),xmdjr1(k,i,l+3),
     .                 xmdjr1(k,i,l+2),k=kst,kfn),i=ist,ifn)
        do i=ist,ifn
           do k=kst,kfn
              xmdj(k,i,l+1,n,iaes) = xmdjr1(k,i,l+1)
              xmdj(k,i,l+2,n,iaes) =-xmdjr1(k,i,l+2)
              xmdj(k,i,l+3,n,iaes) = xmdjr1(k,i,l+3)
           end do
        end do
#endif
      end if
100   continue
110   continue
c
c     read in modal deflections for k=const surfaces
c
      do 210 m = 1,2
c
      if (m.eq.1) then
         k    = 1
         l    = 0
         nseg = nbck0(nbl)
      else
         k    = kdim
         l    = 3
         nseg = nbckdim(nbl)
      end if
c
      do 200 ns = 1,nseg
c
      lbc = iabs(kbcinfo(nbl,ns,1,m))
c
      iflag = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)*
     .        (lbc-2024)*(lbc-2034)*(lbc-2016)
c
      if (iflag.ne.0) goto 200
c
      ist = kbcinfo(nbl,ns,2,m)
      ifn = kbcinfo(nbl,ns,3,m)
      jst = kbcinfo(nbl,ns,4,m)
      jfn = kbcinfo(nbl,ns,5,m)
c
      if(ialph.eq.0) then
#ifdef CMPLX
        read(iunit,*)((xmdkr(j,i,l+1),xmdkr(j,i,l+2),
     .                 xmdkr(j,i,l+3),j=jst,jfn),i=ist,ifn)
        do i=ist,ifn
           do j=jst,jfn
              xmdk(j,i,l+1,n,iaes) = xmdkr(j,i,l+1)
              xmdk(j,i,l+2,n,iaes) = xmdkr(j,i,l+2)
              xmdk(j,i,l+3,n,iaes) = xmdkr(j,i,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdk(j,i,l+1,n,iaes),xmdk(j,i,l+2,n,iaes),
     .                 xmdk(j,i,l+3,n,iaes),j=jst,jfn),i=ist,ifn)
#endif
      else
#ifdef CMPLX
        read(iunit,*)((xmdkr(j,i,l+1),xmdkr(j,i,l+3),
     .                 xmdkr(j,i,l+2),j=jst,jfn),i=ist,ifn)
        do i=ist,ifn
           do j=jst,jfn
              xmdk(j,i,l+1,n,iaes) = xmdkr(j,i,l+1)
              xmdk(j,i,l+2,n,iaes) =-xmdkr(j,i,l+2)
              xmdk(j,i,l+3,n,iaes) = xmdkr(j,i,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdkr1(j,i,l+1),xmdkr1(j,i,l+3),
     .                 xmdkr1(j,i,l+2),j=jst,jfn),i=ist,ifn)
        do i=ist,ifn
           do j=jst,jfn
              xmdk(j,i,l+1,n,iaes) = xmdkr1(j,i,l+1)
              xmdk(j,i,l+2,n,iaes) =-xmdkr1(j,i,l+2)
              xmdk(j,i,l+3,n,iaes) = xmdkr1(j,i,l+3)
           end do
        end do
#endif
      end if
200   continue
210   continue
c
c     read in modal deflections for i=const surfaces
c
      do 310 m = 1,2
c
      if (m.eq.1) then
         i    = 1
         l    = 0
         nseg = nbci0(nbl)
      else
         i    = idim
         l    = 3
         nseg = nbcidim(nbl)
      end if
c
      do 300 ns = 1,nseg
c
      lbc = iabs(ibcinfo(nbl,ns,1,m))
c
      iflag = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)*
     .        (lbc-2024)*(lbc-2034)*(lbc-2016)
c
      if (iflag.ne.0) goto 300
c
      jst = ibcinfo(nbl,ns,2,m)
      jfn = ibcinfo(nbl,ns,3,m)
      kst = ibcinfo(nbl,ns,4,m)
      kfn = ibcinfo(nbl,ns,5,m)
c
      if(ialph.eq.0) then
#ifdef CMPLX
        read(iunit,*)((xmdir(j,k,l+1),xmdir(j,k,l+2),
     .                 xmdir(j,k,l+3),j=jst,jfn),k=kst,kfn)
        do k=kst,kfn
           do j=jst,jfn
              xmdi(j,k,l+1,n,iaes) = xmdir(j,k,l+1)
              xmdi(j,k,l+2,n,iaes) = xmdir(j,k,l+2)
              xmdi(j,k,l+3,n,iaes) = xmdir(j,k,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdi(j,k,l+1,n,iaes),xmdi(j,k,l+2,n,iaes),
     .                 xmdi(j,k,l+3,n,iaes),j=jst,jfn),k=kst,kfn)
#endif
      else
#ifdef CMPLX
        read(iunit,*)((xmdir(j,k,l+1),xmdir(j,k,l+3),
     .                 xmdir(j,k,l+2),j=jst,jfn),k=kst,kfn)
        do k=kst,kfn
           do j=jst,jfn
              xmdi(j,k,l+1,n,iaes) = xmdir(j,k,l+1)
              xmdi(j,k,l+2,n,iaes) =-xmdir(j,k,l+2)
              xmdi(j,k,l+3,n,iaes) = xmdir(j,k,l+3)
           end do
        end do
#else
        read(iunit,*)((xmdir1(j,k,l+1),xmdir1(j,k,l+3),
     .                 xmdir1(j,k,l+2),j=jst,jfn),k=kst,kfn)
        do k=kst,kfn
           do j=jst,jfn
              xmdi(j,k,l+1,n,iaes) = xmdir1(j,k,l+1)
              xmdi(j,k,l+2,n,iaes) =-xmdir1(j,k,l+2)
              xmdi(j,k,l+3,n,iaes) = xmdir1(j,k,l+3)
           end do
        end do
#endif
      end if
300   continue
310   continue
#ifdef CMPLX
c
c     deallocate real arrays in which modal data was read
c
      deallocate(xmdjr)
      deallocate(xmdkr)
      deallocate(xmdir)
#else
c
c     deallocate real arrays in which modal data was read
c
      deallocate(xmdjr1)
      deallocate(xmdkr1)
      deallocate(xmdir1)
#endif
c
      return
      end
