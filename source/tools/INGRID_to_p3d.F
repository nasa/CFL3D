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
      program INGRID_to_p3d
c
c     $Id$
c
c***********************************************************************
c     Purpose: Converts  PEGSUS 4.x INGRID file to a plot3d file that
c     can be used in CFL3D. Note that the INGRID file must correspond
c     to grid points rather than "augmented" cell centers.
c***********************************************************************
c
      character*80 p3dfile,ingfile
c
c     query for input
c
      write(6,*)
      write(6,'(''input name of pegsus INGRID file to read '',
     .          ''(up to 80 characters)'')')
      read(5,'(a80)') ingfile
      write(6,*)
      write(6,'(''enter 0 if this file is unformatted'')')
      write(6,'(''enter 1 if this file is formatted'')')
      read(5,*) ibin_ing
      write(6,*)
      write(6,'(''enter 0 to create a plot3d-type grid'')')
      write(6,'(''enter 1 to create a  cfl3d-type grid'')')
      read(5,*) ip3d
      write(6,*)
      write(6,'(''enter 0 to create a 3d grid'')')
      write(6,'(''enter 1 to create a 2d grid (2 i-planes)'')')
      write(6,'(''(use the 2d option for INGRID files with 3 '',
     .          ''i-planes)'')')
      read(5,*) i2d
      write(6,*)
      if (ip3d .eq. 0) then
         write(6,'(''input name of plot3d grid file to create '',
     .             ''(up to 80 characters)'')')
      else
         write(6,'(''input name of cfl3d grid file to create '',
     .             ''(up to 80 characters)'')')
      end if
      read(5,'(a80)') p3dfile
      write(6,*)
      write(6,'(''enter 0 to create an unformatted grid file'')')
      write(6,'(''enter 1 to create a formatted grid file'')')
      read(5,*) ibin_p3d
c
c     open files
c
      if (ibin_ing .eq. 0) then
         open(unit=1,file=ingfile,form='unformatted',status='old')
      else
         open(unit=1,file=ingfile,form='formatted',status='old')
      end if
      if (ibin_p3d .eq. 0) then
         open(unit=2,file=p3dfile,form='unformatted',status='unknown')
      else
         open(unit=2,file=p3dfile,form='formatted',status='unknown')
      end if
c
c     obtain required array dimensions
c
      maxbl = 1
      lmax  = 1
      jmax  = 1
      kmax  = 1
      write(6,*)
      do ig = 1,9999
         if (ibin_ing .eq. 0) then
            read(1,end=999)
            read(1) jd,kd,ld
         else
           read(1,*,end=999)
            read(1,*) jd,kd,ld
         end if
         write(6,'(''INGRID grid '',i3)') ig
         write(6,'(''  dimensions j x k x l = '',3i4)') jd,kd,ld
         maxbl = max(maxbl,ig)
         lmax  = max(lmax,ld)
         jmax  = max(jmax,jd)
         kmax  = max(kmax,kd)
         if (ibin_ing .eq. 0) then
            read(1) (xx,ll=1,jd*kd*ld),
     .              (yy,ll=1,jd*kd*ld),
     .              (zz,ll=1,jd*kd*ld)
         else
            read(1,*) (xx,ll=1,jd*kd*ld),
     .                (yy,ll=1,jd*kd*ld),
     .                (zz,ll=1,jd*kd*ld)
         end if
      end do
 999  continue
      rewind(1)
c
      ngrid = maxbl
c
      write(6,*)
      write(6,'(''required array sizes: maxbl = '',i6)') maxbl
      write(6,'(''                       lmax = '',i6)') lmax
      write(6,'(''                       jmax = '',i6)') jmax
      write(6,'(''                       kmax = '',i6)') kmax
c
      call convert(maxbl,lmax,jmax,kmax,ngrid,ibin_ing,ibin_p3d,ip3d,
     .             i2d)
c
      write(6,*)
      write(6,'(''conversion completed'')')
c
      stop
      end
c
      subroutine convert(maxbl,lmax,jmax,kmax,ngrid,ibin_ing,ibin_p3d,
     .                   ip3d,i2d)
c
      integer stats
c
      allocatable :: iphntm(:)
      allocatable :: jdim(:)
      allocatable :: kdim(:)
      allocatable :: ldim(:)
      allocatable :: x(:,:,:)
      allocatable :: y(:,:,:)
      allocatable :: z(:,:,:)
c
c     allocate memory
c
      memuse = 0
      allocate( iphntm(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'iphntm',memuse,stats)
      allocate( jdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'jdim',memuse,stats)
      allocate( kdim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'kdim',memuse,stats)
      allocate( ldim(maxbl), stat=stats )
      call umalloc_r(maxbl,1,'ldim',memuse,stats)
      allocate( x(jmax,kmax,lmax), stat=stats )
      call umalloc_r(jmax*kmax*lmax,0,'x',memuse,stats)
      allocate( y(jmax,kmax,lmax), stat=stats )
      call umalloc_r(jmax*kmax*lmax,0,'y',memuse,stats)
      allocate( z(jmax,kmax,lmax), stat=stats )
      call umalloc_r(jmax*kmax*lmax,0,'z',memuse,stats)
c
c     default to no phantom grids
c
      nphntm  = 0
      do ig = 1,ngrid
         iphntm(ig) = 0
      end do
c
      write(6,*)
      write(6,'(''enter number of phantom grids in the'',
     .          '' INGRID file'')')
      read(5,*) nphntm
      if (nphntm .gt. 0) then
         do igph = 1,nphntm
            write(6,'(''enter the global grid number of phantom'',
     .                '' grid '',i2)') igph
            read(5,*) ig
           iphntm(ig) = 1
         end do
      end if
c
      ngout = 0
      do ig = 1,ngrid
         if (ibin_ing .eq. 0) then
            read(1)
            read(1) jd,kd,ld
         else
            read(1,*)
            read(1,*) jd,kd,ld
         end if
         if (iphntm(ig) .eq. 0) then
            ngout = ngout + 1
            jdim(ngout) = jd
            kdim(ngout) = kd
            ldim(ngout) = ld
            if (i2d .gt. 0) then 
               ldim(ngout) = 2
            end if
         end if
         if (ibin_ing .eq. 0) then
            read(1)  (((x(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .               (((y(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .               (((z(j,k,l),j=1,jd),k=1,kd),l=1,ld)
         else
            read(1,*)  (((x(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .                 (((y(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .                 (((z(j,k,l),j=1,jd),k=1,kd),l=1,ld)
         end if
      end do
      rewind(1)
c
      if (ip3d .eq. 0) then
         if (ibin_p3d .eq. 0) then
            write(2) ngout
            write(2) (ldim(n),jdim(n),kdim(n),n=1,ngout)
         else
            write(2,*) ngout
            write(2,*) (ldim(n),jdim(n),kdim(n),n=1,ngout)
         end if
      end if
c
      do ig=1,ngrid
         if (ibin_ing .eq. 0) then
            read(1) 
            read(1) jd,kd,ld
         else
            read(1,*)
            read(1,*) jd,kd,ld
         end if
         ldout = ld
         if (i2d .gt. 0) then
            ldout = 2
         end if
         if (ip3d .ne. 0 .and. iphntm(ig) .eq. 0) then
            if (ibin_p3d .eq. 0) then
               write(2) jd,kd,ldout
            else
               write(2,*) jd,kd,ldout
            end if
         end if
         if (ibin_ing .eq. 0) then
            read(1)  (((x(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .               (((y(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .               (((z(j,k,l),j=1,jd),k=1,kd),l=1,ld)
         else
            read(1,*)  (((x(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .                 (((y(j,k,l),j=1,jd),k=1,kd),l=1,ld),
     .                 (((z(j,k,l),j=1,jd),k=1,kd),l=1,ld)
         end if
         if (iphntm(ig) .eq. 0) then
            if (ip3d .eq. 0) then
               if (ibin_p3d .eq. 0) then
                  write(2) (((x(j,k,l),l=1,ldout),j=1,jd),k=1,kd),
     .                     (((y(j,k,l),l=1,ldout),j=1,jd),k=1,kd),
     .                     (((z(j,k,l),l=1,ldout),j=1,jd),k=1,kd)
                  else
                  write(2,*) (((x(j,k,l),l=1,ldout),j=1,jd),k=1,kd),
     .                       (((y(j,k,l),l=1,ldout),j=1,jd),k=1,kd),
     .                       (((z(j,k,l),l=1,ldout),j=1,jd),k=1,kd)
               end if
            else
               if (ibin_p3d .eq. 0) then
                  write(2) (((x(j,k,l),j=1,jd),k=1,kd),l=1,ldout),
     .                     (((y(j,k,l),j=1,jd),k=1,kd),l=1,ldout),
     .                     (((z(j,k,l),j=1,jd),k=1,kd),l=1,ldout)
               else
                  write(2,*) (((x(j,k,l),j=1,jd),k=1,kd),l=1,ldout),
     .                       (((y(j,k,l),j=1,jd),k=1,kd),l=1,ldout),
     .                       (((z(j,k,l),j=1,jd),k=1,kd),l=1,ldout)
               end if
            end if
         end if
      end do
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(x)
         deallocate(y)
         deallocate(z)
         deallocate(ldim)
         deallocate(jdim)
         deallocate(kdim)
      end if
c
      return
      end
