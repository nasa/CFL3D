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
      program moovmaker
c
c     $Id$
c
c***********************************************************************
c     Purpose: postprocesses 2D/3D SURFACE cfl3d multi-frame plot3d
c     files (generated when cfl3d is run with dt > 0 and abs(movie)>0) 
c     so as to have time as the third (k) direction in the output
c     plot3d files  g.bin and q.bin - the output grid has iblanking
c     and is 3d, multiblock.
c***********************************************************************
c
c     deliberately single precision
      implicit real*4(a-h,o-z)
c
      character*80 g_in,q_in
c
      write(6,*)'name of plot3d grid file to read'
      read(5,'(a80)') g_in
c
      write(6,*)'name of plot3d q file to read'
      read(5,'(a80)') q_in
c
#if defined ASN_P3D
      call asnfile(g_in, '-F f77 -N ieee', IER)
      call asnfile(q_in, '-F f77 -N ieee', IER)
#endif
      open(unit=1,file=g_in,form='unformatted',status='old')      
      open(unit=2,file=q_in,form='unformatted',status='old')      
c
c     NOTE: ibflag indicates if input  plot3d files have iblank array
c     ibflag = 0 no iblank in input grid file
c     ibflag = 1 iblank in input grid file (cfl3d standard)
c
      ibflag = 1
c
c     NOTE: numgrd indicates whether plot3d grid file contains just one
c           copy of the grid or ntime copies of the grid
c     numgrd = 0 just one copy
c     numgrd = 1 ntime copies of the grid
c
      write(6,*)'stationary or moving grid?...0=stationary, 1=moving'
      read(5,*) numgrd
c
      write(6,*)'3D or 2D grid?...0=3D, 1=2d'
      read(5,*) i2d
c            
c     read plot3d grid and solution files (3d/mg) to get sizing info
c
      read(1) maxbl
c
      write(6,*)'how many frames are in the plot3d file?'
      read(5,*) ntime
c
      call sizejk(maxbl,jmax,kmax,i2d)
c
      rewind(1)
c
      call makemoov(maxbl,jmax,kmax,ntime,ibflag,numgrd,i2d)
c
      stop
      end
c
      subroutine sizejk(maxbl,jmax,kmax,i2d)
c***********************************************************************
c     Purpose: determine max j,k dimensions needed for this case
c***********************************************************************
c
c     deliberately single precision
      implicit real*4(a-h,o-z)
c
      integer stats
c
      allocatable :: itemp(:)
      allocatable :: jtemp(:)
      allocatable :: ktemp(:)
c
c     allocate memory
c
      memuse = 0
      allocate( itemp(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'itemp',memuse,stats)
      allocate( jtemp(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'jtemp',memuse,stats)
      allocate( ktemp(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'ktemp',memuse,stats)
c
       if (i2d.eq.0) then
          read(1) (itemp(l),jtemp(l),ktemp(l),l=1,maxbl)
       else
          read(1) (jtemp(l),ktemp(l),l=1,maxbl)
          do l=1,maxbl
             itemp(l) = 1
          end do
       end if
       rewind(1)
c
       jmax = 1
       kmax = 1
       do l=1,maxbl
          if (itemp(l) .eq. 1) then
             jmax = max(jmax,jtemp(l))
             kmax = max(kmax,ktemp(l))
          else if (jtemp(l) .eq. 1) then
             jmax = max(jmax,itemp(l))
             kmax = max(kmax,ktemp(l))
          else if (ktemp(l) .eq. 1) then
             jmax = max(jmax,itemp(l))
             kmax = max(kmax,jtemp(l))
          else
             write(6,*)'error - for 3-D movies, one output dimension',
     +        ' must be of dimension 1'
             write(6,*)'... (kstart=kend=1 or jstart=jend=1 or',
     +        ' istart=iend=1 in input file)'
             write(6,*)'stopping'
             stop
          end if
       end do
       write(6,*)'maxbl,jmax,kmax',maxbl,jmax,kmax
c
c     dealocate memory
c
      deallocate(itemp)
      deallocate(jtemp)
      deallocate(ktemp)
c
      return
      end
c
      subroutine makemoov(maxbl,jmax,kmax,ntime,ibflag,numgrd,i2d)
c***********************************************************************
c     Purpose: read in concatenated plot3d frames generated by cfl3d
c     and ouput as a standard plot3d file with time as the third (k)
c     direction. The new files are g.bin and q.bin, and are iblanked
c     unless ibflag = 0.
c***********************************************************************
c
c     deliberately single precision
      implicit real*4(a-h,o-z)
c
      integer stats

      allocatable :: iblank(:,:,:)
      allocatable :: idim(:)
      allocatable :: jdim(:)
      allocatable :: jjdim(:)
      allocatable :: kdim(:)
      allocatable :: kkdim(:)
      allocatable :: q(:,:,:,:,:)
      allocatable :: time(:)
      allocatable :: x(:,:,:,:)
      allocatable :: y(:,:,:,:)
      allocatable :: z(:,:,:,:)
c
      memuse = 0
c
      allocate( iblank(jmax,kmax,maxbl), stat=stats )
      call umalloc_r(jmax*kmax*maxbl,1,'iblank',memuse,stats)
      allocate( idim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'idim',memuse,stats)
      allocate( jdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'jdim',memuse,stats)
      allocate( jjdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'jjdim',memuse,stats)
      allocate( kdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'kdim',memuse,stats)
      allocate( kkdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'kkdim',memuse,stats)
      allocate( q(jmax,kmax,ntime,maxbl,5), stat=stats )
      call umalloc_r(jmax*kmax*ntime*maxbl*5,0,'q',memuse,stats)
      allocate( time(ntime), stat=stats )
      call umalloc_r(ntime,0,'time',memuse,stats)
      allocate( x(jmax,kmax,ntime,maxbl), stat=stats )
      call umalloc_r(jmax*kmax*ntime*maxbl,0,'x',memuse,stats)
      allocate( y(jmax,kmax,ntime,maxbl), stat=stats )
      call umalloc_r(jmax*kmax*ntime*maxbl,0,'y',memuse,stats)
      allocate( z(jmax,kmax,ntime,maxbl), stat=stats )
      call umalloc_r(jmax*kmax*ntime*maxbl,0,'z',memuse,stats)
c
      do 100 nt = 1,ntime
c
      if(numgrd.eq.0 .and. nt.gt.1) go to 11
      read(1) nbl
      write(6,*)'nbl = ',nbl
      if (i2d.eq.0) then
         read(1) (idim(nn),jdim(nn),kdim(nn),nn=1,nbl)
      else
         read(1) (jdim(nn),kdim(nn),nn=1,nbl)
         do nn=1,nbl
            idim(nn) = 1
         end do
      end if
      do 10 nn=1,nbl
      write(6,*)'block ',nn,' i x j x k ',idim(nn),jdim(nn),kdim(nn)
      call my_flush(6)
10    continue
11    continue
      read(2) idum
      if (i2d.eq.0) then
         read(2)(idum,idum,idum,nn=1,nbl)
      else
         read(2)(idum,idum,nn=1,nbl)
      end if
c
      do 21 n=1,nbl
      if (idim(n).eq.1) then
         jd = jdim(n)
         kd = kdim(n)
         jjdim(n) = jd
         kkdim(n) = kd
      else if (jdim(n).eq.1) then
         jd = idim(n)
         kd = kdim(n)
         jjdim(n) = jd
         kkdim(n) = kd
      else if (kdim(n).eq.1) then
         jd = idim(n)
         kd = jdim(n)
         jjdim(n) = jd
         kkdim(n) = kd
      end if
      if(numgrd.eq.0 .and. nt.gt.1) go to 12
      write(6,'('' reading surface grid of dimensions i x j x k = '',
     .          3i5)') idim(n),jdim(n),kdim(n)
      if (ibflag.eq.0) then
         if (i2d.eq.0) then
            read(1) ((x(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((y(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((z(j,k,nt,n),j=1,jd),k=1,kd)
         else
            read(1) ((x(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((y(j,k,nt,n),j=1,jd),k=1,kd)
            do j=1,jd
               do k=1,kd
                  z(j,k,nt,n) = 0.
               end do
            end do
         end if
      else
         if (i2d.eq.0) then
            read(1) ((x(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((y(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((z(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((iblank(j,k,n),j=1,jd),k=1,kd)
         else
            read(1) ((x(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((y(j,k,nt,n),j=1,jd),k=1,kd),
     .              ((iblank(j,k,n),j=1,jd),k=1,kd)
            do j=1,jd
               do k=1,kd
                  z(j,k,nt,n) = 0.
               end do
            end do
         end if
      end if
      write(6,*)'completed read of grid file for time step ',nt
12    continue
c
      read(2) xmach,alpha,re,time(nt)
      if (i2d.eq.0) then
         read(2) (((q(j,k,nt,n,l),j=1,jd),k=1,kd),l=1,5)
      else
         read(2) (((q(j,k,nt,n,l),j=1,jd),k=1,kd),l=1,4)
         do k=1,kd
            do j=1,jd
               q(j,k,nt,n,5) = q(j,k,nt,n,4)
               q(j,k,nt,n,4) = 0.
            end do
         end do
      end if
      write(6,*)'completed read of soln file for time step ',nt
c
21    continue
c
100   continue
c
c     output new plot3d files with time as the third direction
c
      write(6,*) 'creating file g.bin (plot3d grid file with time',
     .           ' as k-direction)'
      if (ibflag.gt.0) then
         write(6,*) 'grid file is iblanked'
      end if
      write(6,*) 'creating file q.bin (plot3d q file with time',
     .           ' as k-direction)'
c
#if defined ASN_P3D
      call asnfile('g.bin', '-F f77 -N ieee', IER)
      call asnfile('q.bin', '-F f77 -N ieee', IER)
#endif
      open(unit=3,file='g.bin',form='unformatted',status='unknown')
      open(unit=4,file='q.bin',form='unformatted',status='unknown')
c
      write(3) nbl
      write(3) (jjdim(nn),kkdim(nn),ntime,nn=1,nbl)
      write(4) nbl
      write(4) (jjdim(nn),kkdim(nn),ntime,nn=1,nbl)
      do 31 n=1,nbl
      jd = jjdim(n)
      kd = kkdim(n)
      write(6,*) ' writing 3d (time) grid WITH IBLANK of dimensions',
     .' j x k x nt = ',jd,kd,ntime
      if (numgrd.gt.0) then
         if (ibflag.eq.0) then
            write(3) (((x(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((y(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((z(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime)
         else
            write(3) (((x(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((y(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((z(j,k,nt,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((iblank(j,k,n),j=1,jd),k=1,kd),nt=1,ntime)
         end if
      else
         if (ibflag.eq.0) then
            write(3) (((x(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((y(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((z(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime)
         else
            write(3) (((x(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((y(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((z(j,k,1,n),j=1,jd),k=1,kd),nt=1,ntime),
     .               (((iblank(j,k,n),j=1,jd),k=1,kd),nt=1,ntime)
         end if
      end if
      write(4) xmach,alpha,re,time(1)
      write(4) ((((q(j,k,nt,n,l),j=1,jd),k=1,kd),nt=1,ntime),l=1,5)
31    continue
c
      write(6,*)'completed writing new plot3d grid file: g.bin'
      write(6,*)'completed writing new plot3d q    file: q.bin'
c
c     free memory
c
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(q)
      deallocate(iblank)
      deallocate(idim)
      deallocate(jdim)
      deallocate(kdim)
      deallocate(jjdim)
      deallocate(kkdim)
      deallocate(time)
c
      return
      end
