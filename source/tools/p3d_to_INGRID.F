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
      program p3d_to_INGRID
c
c     $Id$
c
c***********************************************************************
c     Purpose: Converts either PLOT3D or CFL3D type grids into either
c     INGRID type grids that can be used with PEGSUS 4.x, or PLOT3D
c     type grids that can be used with PEGSUS 5.x. The converted grids
c     can contain either the grid points as given in the input grids,
c     or "augmented" cell centers of the input grids.
c***********************************************************************
c
      dimension ldim(10000),jdim(10000),kdim(10000)
c
      character*80 ingfile,p3dfile
      character*1 yesno
c
c     set file defaults:
c
c     input grid file(s) is plot3d type 
      ip3d = 0
c
c     input grid file(s) is unformatted
      ibin_in = 0
c
c     output INGRID/plot3d file is unformatted
      ibin_out = 0
c
      write(6,*)
      write(6,'(''file defaults:'')')
      write(6,'(''   input grid files are plot3d type'')')
      write(6,'(''   input grid files are unformatted'')')
      write(6,'(''   output INGRID/plot3d file is unformatted'')')
      write(6,*)
      write(6,'(''do you wish to use these defaults (y/n)?'')')
      write(6,'(''   (alternate options include formatted'')')
      write(6,'(''   files and cfl3d-type input grid files)'')')
      read(5,'(a1)') yesno 
      if (yesno .eq. 'n' .or. yesno .eq. 'N') then
         write(6,*)
         write(6,'(''enter 0 to convert plot3d-type grids'')')
         write(6,'(''enter 1 to convert  cfl3d-type grids'')')
         read(5,*) ip3d
         write(6,*)
         write(6,'(''enter 0 for unformatted grid files '')')
         write(6,'(''enter 1 for   formatted grid files '')')
         read(5,*) ibin_in
         write(6,*)
         write(6,'(''enter 0 to create an unformatted '',
     .             ''INGRID/plot3d file'')')
         write(6,'(''enter 1 to create an   formatted '',
     .             ''INGRID/plot3d file'')')
         read(5,*) ibin_out
      end if
      write(6,*)
      write(6,'(''choose the type of output file to create'')')
      write(6,'(''enter 0 to create an INGRID file (for PEGSUS 4.x)'')')
      write(6,'(''enter 1 to create a plot3d file  (for PEGSUS 5.x)'')')
      read(5,*) iouttyp
      write(6,*)
      write(6,'(''enter the name of the output file to create '',
     .          ''(up to 80 char.)'')')
      read(5,'(a80)') ingfile
      write(6,*)
      write(6,'(''enter 0 to create an output file with grid points'')')
      write(6,'(''enter 1 to create an output file with augmented '',
     .          ''cell centers'')')
      read(5,*) icc
      write(6,*)
      write(6,'(''enter 0 to preserve input-grid i,j,k index '',
     .          ''definitions in the output grid'')')
      write(6,'(''you may want this option if you have an '',
     .          ''existing PEGSUS input file that'')')
      write(6,'(''was generated for OVERFLOW'')')
      write(6,'(''  Note: this will require the following '',
     .          ''translation of'')')
      write(6,'(''  indicies between CFL3D and PEGSUS input files:'')')
      write(6,'(''            PEGSUS   CFL3D'')')
      write(6,'(''              J    =   I'')')
      write(6,'(''              K    =   J'')')
      write(6,'(''              L    =   K'')')
      write(6,*)
      write(6,'(''enter 1 to swap input-grid i,j,k index '',
     .          ''definitions in the output grid'')')
      write(6,'(''  Note: this will require NO translation of'')')
      write(6,'(''  indicies between CFL3D and PEGSUS input files:'')')
      write(6,'(''            PEGSUS   CFL3D'')')
      write(6,'(''              L    =   I'')')
      write(6,'(''              J    =   J'')')
      write(6,'(''              K    =   K'')')
      read(5,*) iswap
c
      write(6,'(''enter number of separate grid files to convert '',
     .          ''into one output file'')')
      if (iouttyp.eq.1) then
         write(6,'(''NOTE: for plot3d output, only one grid file '',
     .             ''can be converted at a time'')')
      end if
      read(5,*) nfiles
      if (iouttyp.eq.1 .and. nfiles.gt.1) then 
         write(6,'(''sorry...will only process 1 file'')')
         nfiles = 1
      end if
c
      iname = 0
      nname = 0
      if (iouttyp .eq. 0) then
         write(6,'(''enter 0 to specify a name for each mesh'')')
         write(6,'(''enter 1 to  use default names (grid.n)'')')
         read(5,*) iname
      end if
c
      if (ibin_out .eq. 0) then
         open(unit=2,file=ingfile,form='unformatted',status='unknown')
      else
         open(unit=2,file=ingfile,form='formatted',status='unknown')
      end if
c
      do nf = 1,nfiles
c
         write(6,*)
         write(6,'(''begining processing of grid file number '',i3)') nf
         write(6,*)
         if (ip3d .eq. 0) then
            if (ibin_in .eq. 0) then
               write(6,'(''input name of unformatted plot3d grid '',
     .                   ''file to read (up to 80 characters)'')')
            else
               write(6,'(''input name of formatted plot3d grid '',
     .                   ''file to read (up to 80 characters)'')')
            end if
            read(5,'(a80)') p3dfile
            write(6,*)
            write(6,'(''enter 0 if a   single-grid plot3d file'')')
            write(6,'(''enter 1 if a multiple-grid plot3d file'')')
            read(5,*) img
         else
            if (ibin_in .eq. 0) then 
               write(6,'(''input name of unformatted cfl3d grid '',
     .                   ''file to read (up to 80 characters)'')')
            else
               write(6,'(''input name of formatted cfl3d grid '',
     .                   ''file to read (up to 80 characters)'')')
            end if
            read(5,'(a80)') p3dfile
         end if
c
         if (ibin_in .eq. 0) then
            open(unit=1,file=p3dfile,form='unformatted',status='old')
         else
            open(unit=1,file=p3dfile,form='formatted',status='old')
         end if
c
c        obtain required array dimensions
c
         if (ip3d .eq. 0) then
            if (ibin_in .eq. 0) then
               if (img .ne. 0) then
                  read(1) ngrid
               else
                  ngrid = 1
               end if
               read(1) (ldim(n),jdim(n),kdim(n),n=1,ngrid)
            else
               if (img .ne. 0) then
                  read(1,*) ngrid
               else
                  ngrid = 1
               end if
               read(1,*) (ldim(n),jdim(n),kdim(n),n=1,ngrid)
            end if
         else
            ngrid = 0
            if (ibin_in .eq. 0) then
               do n=1,9999
                  read(1,end=9997) jdim(n),kdim(n),ldim(n)
                  ngrid = ngrid + 1
                  read(1) (xx,ll=1,jdim(n)*kdim(n)*ldim(n)),
     .                    (yy,ll=1,jdim(n)*kdim(n)*ldim(n)),
     .                    (zz,ll=1,jdim(n)*kdim(n)*ldim(n))
               end do
 9997          continue
            else
               do n=1,9999
                  read(1,*,end=9998) jdim(n),kdim(n),ldim(n)
                  ngrid = ngrid + 1
                  read(1,*) (xx,ll=1,jdim(n)*kdim(n)*ldim(n)),
     .                      (yy,ll=1,jdim(n)*kdim(n)*ldim(n)),
     .                      (zz,ll=1,jdim(n)*kdim(n)*ldim(n))
               end do
 9998          continue
            end if
         end if
         rewind(1)
c
         maxbl = ngrid
         lmax  = 1
         jmax  = 1
         kmax  = 1
         do n=1,ngrid
            lmax  = max(lmax,ldim(n))
            jmax  = max(jmax,jdim(n))
            kmax  = max(kmax,kdim(n))
         end do
         i2d = 0
         if (lmax .eq. 2) then
            i2d = 1
            lmax = lmax + 1
         end if
c
c        increase for augmented cell-center grid
c
         if (icc .eq. 1) then
            lmax = lmax + 1
            jmax = jmax + 1
            kmax = kmax + 1
         end if
c
         write(6,*)
         write(6,'(''required array sizes: maxbl = '',i6)') maxbl
         write(6,'(''                       lmax = '',i6)') lmax
         write(6,'(''                       jmax = '',i6)') jmax
         write(6,'(''                       kmax = '',i6)') kmax
c
         call convert(maxbl,lmax,jmax,kmax,ngrid,ibin_out,ibin_in,ip3d,
     .                img,i2d,iouttyp,icc,iswap,jdim,kdim,ldim,iname,
     .                nname)
c
         write(6,*)
         write(6,'(''conversion of grid file '',i3,'' complete'')') nf
      end do
c
      write(6,*)
      write(6,'(''conversion of all grid files completed'')')
c
      write(6,*)
      if (iouttyp .eq. 0) then
         write(6,'(''the INGRID output file:'')')
      else
         write(6,'(''the plot3d output file:'')')
      endif
      write(6,'(a80)') ingfile
      if (icc .eq. 0) then
         write(6,'(''contains grid points of the input grid(s)'')')
      else
         write(6,'(''contains (augmented) cell centers of the '',
     .             ''input grid(s)'')')
      end if
      if (iouttyp .eq. 0) then
         write(6,'(''note: the mesh names in the INGRID file contain '',
     .             ''40 characters'')')
         write(6,'(''make sure the PEGSUS parameter ICHAR is set '',
     .             ''to 40'')')
      end if
c
      stop
      end
c
      subroutine convert(maxbl,lmax,jmax,kmax,ngrid,ibin_out,ibin_in,
     .                   ip3d,img,i2d,iouttyp,icc,iswap,jdim,kdim,ldim,
     .                   iname,nname)
c
      character*40 string
c
      integer stats
c
      dimension jdim(10000)
      dimension kdim(10000)
      dimension ldim(10000)
c
      allocatable :: x(:,:,:)
      allocatable :: xcc(:,:,:)
      allocatable :: y(:,:,:)
      allocatable :: ycc(:,:,:)
      allocatable :: z(:,:,:)
      allocatable :: zcc(:,:,:)
c
c     allocate memory
c
      memuse = 0
      allocate( x(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'x',memuse,stats)
      allocate( xcc(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'xcc',memuse,stats)
      allocate( y(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'y',memuse,stats)
      allocate( ycc(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'ycc',memuse,stats)
      allocate( z(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'z',memuse,stats)
      allocate( zcc(lmax,jmax,kmax), stat=stats )
      call umalloc_r(lmax*jmax*kmax,0,'zcc',memuse,stats)
c
      write(6,*)
c
      if (ip3d .eq. 0) then
         if (ibin_in .eq. 0) then
            if (img .ne. 0) then
               read(1) ngrid
            else
               ngrid = 1
            end if
            read(1) (ldim(n),jdim(n),kdim(n),n=1,ngrid)
         else
            if (img .ne. 0) then
               read(1,*) ngrid
            else
               ngrid = 1
            end if
            read(1,*) (ldim(n),jdim(n),kdim(n),n=1,ngrid)
         end if
      end if
c
      icchold = icc
      do n=1,ngrid
         if (ldim(n).eq.1 .or. jdim(n).eq.1 .or. kdim(n).eq.1) icc = 0
      end do
      iadd = 0
      if (icc .gt. 0) iadd = 1
c
      if (iouttyp .eq. 1) then
         if (ibin_out .eq. 0) then
            write(2) ngrid
            if (iswap .eq. 0) then
               if (i2d .eq. 0) then 
                  write(2) (ldim(n)+iadd,jdim(n)+iadd,kdim(n)+iadd,
     .                      n=1,ngrid)
               else
                  write(2) (ldim(n)+iadd+1,jdim(n)+iadd,kdim(n)+iadd,
     .                      n=1,ngrid)
               end if
            else
               if (i2d .eq. 0) then
                  write(2) (jdim(n)+iadd,kdim(n)+iadd,ldim(n)+iadd,
     .                      n=1,ngrid)
               else
                  write(2) (jdim(n)+iadd,kdim(n)+iadd,ldim(n)+iadd+1,
     .                      n=1,ngrid)
               end if
            end if
         else
            write(2,*) ngrid
            if (iswap .eq. 0) then
               if (i2d .eq. 0) then
                  write(2,*) (ldim(n)+iadd,jdim(n)+iadd,kdim(n)+iadd,
     .                      n=1,ngrid)
               else
                  write(2,*) (ldim(n)+iadd+1,jdim(n)+iadd,kdim(n)+iadd,
     .                      n=1,ngrid)
               end if
            else
               if (i2d .eq. 0) then
                  write(2,*) (jdim(n)+iadd,kdim(n)+iadd,ldim(n)+iadd,
     .                      n=1,ngrid)
               else
                  write(2,*) (jdim(n)+iadd,kdim(n)+iadd,ldim(n)+iadd+1,
     .                      n=1,ngrid)
               end if
            end if
         end if
      end if
c
      do n=1,ngrid
         if (ip3d .ne. 0) then
            if (ibin_in  .eq. 0) then
               read(1) jdim(n),kdim(n),ldim(n)
            else
               read(1,*) jdim(n),kdim(n),ldim(n)
            end if
         end if
c
c        query for mesh name if an INGRID file is to be created
c
         nname = nname + 1
         if (iouttyp .eq. 0) then
            if (iname.eq.0) then
               write(6,*)
               write(6,'(''input name (up to 40 char.) for zone '',
     .         i4)') nname
               read(5,'(a40)') string
            else
               if (nname.gt.99) then
                  len1 = 8
                  write(string,'("grid.",i3)') nname
               else if (nname.gt.9) then
                  len1 = 7
                  write(string,'("grid.",i2)') nname
               else
                  len1 = 6
                  write(string,'("grid.",i1)') nname
               endif
               do i = len1+1, 40
                  string(i:i) = ' '
               end do
            end if
         end if
c
         jd =jdim(n)
         kd =kdim(n)
         ld =ldim(n)
c
         write(6,*)
         write(6,'(''reading zone '',i4)') n
         write(6,'(''  input dimensions '',3i4)') ld,jd,kd
c
         if (ip3d .eq. 0) then
            if (ibin_in .eq. 0) then
               read(1) (((x(l,j,k),l=1,ld),j=1,jd),k=1,kd),
     .                 (((y(l,j,k),l=1,ld),j=1,jd),k=1,kd),
     .                 (((z(l,j,k),l=1,ld),j=1,jd),k=1,kd)
            else
               read(1,*) (((x(l,j,k),l=1,ld),j=1,jd),k=1,kd),
     .                   (((y(l,j,k),l=1,ld),j=1,jd),k=1,kd),
     .                   (((z(l,j,k),l=1,ld),j=1,jd),k=1,kd)
            end if
         else
            if (ibin_in .eq. 0) then
               read(1) (((x(l,j,k),j=1,jd),k=1,kd),l=1,ld),
     .                 (((y(l,j,k),j=1,jd),k=1,kd),l=1,ld),
     .                 (((z(l,j,k),j=1,jd),k=1,kd),l=1,ld)
            else
               read(1,*) (((x(l,j,k),j=1,jd),k=1,kd),l=1,ld),
     .                   (((y(l,j,k),j=1,jd),k=1,kd),l=1,ld),
     .                   (((z(l,j,k),j=1,jd),k=1,kd),l=1,ld)
            end if
         end if
c
c        add extra y plane if "2d" case
c
         if (i2d .gt. 0) then
            ygrad = y(ld,1,1)-y(ld-1,1,1)
            do j=1,jd
               do k=1,kd
                  x(ld+1,j,k) = x(ld,j,k)
                  y(ld+1,j,k) = y(ld,j,k) + ygrad
                  z(ld+1,j,k) = z(ld,j,k)
               end do
            end do
            ld = ld + 1
         end if
c
         if (icc .gt. 0) then
c
c           create augmented cell-center grid
c
c           irind = 0...augmented "rind" cells correspond to 
c                       cell-face centers on original block boundaries.
c                       this implies that the physical domain of the
c                       augmented grid is exactly the same as the 
c                       original grid.
c                   1...augmented "rind" cells are extrapolated from 
c                       interior cell centers, and therefore extend 
c                       outside of the original block boundaries.
c                       this is the method used on PEGSUS 4.2, but
c                       Stu Rogers and I decided that this could lead
c                       to problems...PEGSUS 5.1+ is expected to use
c                       the equivalent of irind=0.
c
            irind = 0
            call cellcen(x,y,z,xcc,ycc,zcc,lmax,jmax,kmax,
     .                   jd,kd,ld,irind)
c
         else
c
c           just copy input x,y,z to xcc,ycc,zcc to simplify
c           output coding
c
            do k=1,kd
               do j=1,jd
                  do l=1,ld
                     xcc(l,j,k) = x(l,j,k)
                     ycc(l,j,k) = y(l,j,k)
                     zcc(l,j,k) = z(l,j,k)
                  end do
               end do
            end do
               
         end if
c
         jd1 = jd+iadd
         kd1 = kd+iadd
         ld1 = ld+iadd
c
         if (iouttyp .eq. 0) then
            write(6,'(''writing zone '',i4,'' with name '',a40)') 
     .      n,string
         else
            write(6,'(''writing zone '',i4)') n
         end if
         if (iswap .eq. 0) then
            write(6,'(''  output dimensions '',3i4)') ld1,jd1,kd1
         else
            write(6,'(''  output dimensions '',3i4)') jd1,kd1,ld1
         end if
c
         if (ibin_out .eq. 0) then
            if (iswap .eq. 0) then
               if (iouttyp .eq. 0) then
                  write(2) string
                  write(2) ld1,jd1,kd1
               end if
               write(2) (((xcc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1),
     .                  (((ycc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1),
     .                  (((zcc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1)
            else
               if (iouttyp .eq. 0) then
                  write(2) string
                  write(2) jd1,kd1,ld1
               end if
               write(2) (((xcc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1),
     .                  (((ycc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1),
     .                  (((zcc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1)
            end if
         else
            if (iswap .eq. 0) then
               if (iouttyp .eq. 0) then
                  write(2,*) string
                  write(2,*) ld1,jd1,kd1
               end if
               write(2,*) (((xcc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1),
     .                    (((ycc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1),
     .                    (((zcc(l,j,k),l=1,ld1),j=1,jd1),k=1,kd1)
            else
               if (iouttyp .eq. 0) then
                  write(2,*) string
                  write(2,*) jd1,kd1,ld1
               end if
               write(2,*) (((xcc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1),
     .                    (((ycc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1),
     .                    (((zcc(l,j,k),j=1,jd1),k=1,kd1),l=1,ld1)
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
         deallocate(xcc)
         deallocate(ycc)
         deallocate(zcc)
      end if
      icc = icchold
c
      return
      end
c
      subroutine cellcen(x,y,z,xcc,ycc,zcc,lmax,jmax,kmax,
     .                   jd,kd,ld,irind)
c***************************************************************** 
c     Purpose: find the cell centers (xcc,ycc,zcc) of the grid
c     (x,y,z) and augment those with a "rind" layer of cells that
c     corresponds either to cell-face centers or cell centers
c     extrapolated from interior cell-center data.
c
c     irind = 0...rind layer is cell-face center
c             1...rind layer is extrapolated from interior
c                 cell-centers
c***************************************************************** 
c
      dimension x(lmax,jmax,kmax),y(lmax,jmax,kmax),
     .          z(lmax,jmax,kmax),
     .          xcc(lmax,jmax,kmax),ycc(lmax,jmax,kmax),
     .          zcc(lmax,jmax,kmax)
c
c     interior cell centers
c 
      do k=1,kd-1 
         kk = k+1
         do j=1,jd-1
            jj = j+1
            do l=1,ld-1
               ll = l+1
               xcc(ll,jj,kk) = (x(l,j,k)     + x(l,j+1,k)
     .                       +  x(l,j,k+1)   + x(l,j+1,k+1)
     .                       +  x(l+1,j,k)   + x(l+1,j+1,k)
     .                       +  x(l+1,j,k+1) + x(l+1,j+1,k+1))/8.0
               ycc(ll,jj,kk) = (y(l,j,k)     + y(l,j+1,k)
     .                       +  y(l,j,k+1)   + y(l,j+1,k+1)
     .                       +  y(l+1,j,k)   + y(l+1,j+1,k)
     .                       +  y(l+1,j,k+1) + y(l+1,j+1,k+1))/8.0
               zcc(ll,jj,kk) = (z(l,j,k)     + z(l,j+1,k)
     .                       +  z(l,j,k+1)   + z(l,j+1,k+1)
     .                       +  z(l+1,j,k)   + z(l+1,j+1,k)
     .                       +  z(l+1,j,k+1) + z(l+1,j+1,k+1))/8.0
            end do
         end do
      end do
c
c     rind cells
c
      if (irind .eq. 0) then
c
c        cell-face centers on l=const faces
c
         do l=1,ld,ld-1
            ll = l
            if (l.eq.ld) ll = l+1
            do k=1,kd-1
               kk = k+1
               do j=1,jd-1
                  jj = j+1
                  xcc(ll,jj,kk) = (x(l,j,k)   + x(l,j+1,k)
     .                          +  x(l,j,k+1) + x(l,j+1,k+1))/4.0
                  ycc(ll,jj,kk) = (y(l,j,k)   + y(l,j+1,k)
     .                          +  y(l,j,k+1) + y(l,j+1,k+1))/4.0
                  zcc(ll,jj,kk) = (z(l,j,k)   + z(l,j+1,k)
     .                          +  z(l,j,k+1) + z(l,j+1,k+1))/4.0
               end do
            end do
         end do
c
c        cell-face centers on j=const faces
c
         do j=1,jd,jd-1
            jj = j
            if (j.eq.jd) jj = j+1
            do k=1,kd-1
               kk = k+1
               do l=1,ld-1
                  ll = l+1
                  xcc(ll,jj,kk) = (x(l,j,k)   + x(l+1,j,k)
     .                          +  x(l,j,k+1) + x(l+1,j,k+1))/4.0
                  ycc(ll,jj,kk) = (y(l,j,k)   + y(l+1,j,k)
     .                          +  y(l,j,k+1) + y(l+1,j,k+1))/4.0
                  zcc(ll,jj,kk) = (z(l,j,k)   + z(l+1,j,k)
     .                          +  z(l,j,k+1) + z(l+1,j,k+1))/4.0
               end do
            end do
         end do
c
c        cell-face centers on k=const faces
c
         do k=1,kd,kd-1
            kk = k
            if (k.eq.kd) kk = k+1
            do j=1,jd-1
               jj = j+1
               do l=1,ld-1
                  ll = l+1
                  xcc(ll,jj,kk) = (x(l,j,k)   + x(l,j+1,k)
     .                          +  x(l+1,j,k) + x(l+1,j+1,k))/4.0
                  ycc(ll,jj,kk) = (y(l,j,k)   + y(l,j+1,k)
     .                          +  y(l+1,j,k) + y(l+1,j+1,k))/4.0
                  zcc(ll,jj,kk) = (z(l,j,k)   + z(l,j+1,k)
     .                          +  z(l+1,j,k) + z(l+1,j+1,k))/4.0
               end do
            end do
         end do
c
c        fill in edge values with grid-edge midpoints
c
         do k=1,kd,kd-1
            kk = k
            if (k.eq.kd) kk = k+1
            do j=1,jd,jd-1
               jj = j
               if (j.eq.jd) jj = j+1
               do l=1,ld-1
                  ll = l+1
                  xcc(ll,jj,kk) = (x(l,j,k) + x(l+1,j,k))/2.0
                  ycc(ll,jj,kk) = (y(l,j,k) + y(l+1,j,k))/2.0
                  zcc(ll,jj,kk) = (z(l,j,k) + z(l+1,j,k))/2.0
               end do
            end do
         end do
c
         do k=1,kd-1
            kk = k+1
            do j=1,jd,jd-1
               jj = j
               if (j.eq.jd) jj = j+1
               do l=1,ld,ld-1
                  ll = l
                  if (l.eq.ld) ll = l+1
                  xcc(ll,jj,kk) = (x(l,j,k) + x(l,j,k+1))/2.0
                  ycc(ll,jj,kk) = (y(l,j,k) + y(l,j,k+1))/2.0
                  zcc(ll,jj,kk) = (z(l,j,k) + z(l,j,k+1))/2.0
               end do
            end do
         end do
c
         do k=1,kd,kd-1
            kk = k
            if (k.eq.kd) kk = k+1
            do j=1,jd-1
               jj = j+1
               do l=1,ld,ld-1
                  ll = l
                  if (l.eq.ld) ll = l+1
                  xcc(ll,jj,kk) = (x(l,j,k) + x(l,j+1,k))/2.0
                  ycc(ll,jj,kk) = (y(l,j,k) + y(l,j+1,k))/2.0
                  zcc(ll,jj,kk) = (z(l,j,k) + z(l,j+1,k))/2.0
               end do
            end do
         end do
c
c        fill in corner values with grid corner values
c
         xcc(1,1,1)          = x(1,1,1)
         xcc(1,jd+1,1)       = x(1,jd,1)
         xcc(1,1,kd+1)       = x(1,1,kd)
         xcc(1,jd+1,kd+1)    = x(1,jd,kd)
         xcc(ld+1,1,1)       = x(ld,1,1)
         xcc(ld+1,jd+1,1)    = x(ld,jd,1)
         xcc(ld+1,1,kd+1)    = x(ld,1,kd)
         xcc(ld+1,jd+1,kd+1) = x(ld,jd,kd)
         ycc(1,1,1)          = y(1,1,1)
         ycc(1,jd+1,1)       = y(1,jd,1)
         ycc(1,1,kd+1)       = y(1,1,kd)
         ycc(1,jd+1,kd+1)    = y(1,jd,kd)
         ycc(ld+1,1,1)       = y(ld,1,1)
         ycc(ld+1,jd+1,1)    = y(ld,jd,1)
         ycc(ld+1,1,kd+1)    = y(ld,1,kd)
         ycc(ld+1,jd+1,kd+1) = y(ld,jd,kd)
         zcc(1,1,1)          = z(1,1,1)
         zcc(1,jd+1,1)       = z(1,jd,1)
         zcc(1,1,kd+1)       = z(1,1,kd)
         zcc(1,jd+1,kd+1)    = z(1,jd,kd)
         zcc(ld+1,1,1)       = z(ld,1,1)
         zcc(ld+1,jd+1,1)    = z(ld,jd,1)
         zcc(ld+1,1,kd+1)    = z(ld,1,kd)
         zcc(ld+1,jd+1,kd+1) = z(ld,jd,kd)
c
      else
c
c        exratapolate from interior in l-direction for l=const faces
c
         do l=1,ld+1,ld
            l1 = 1
            l2 = 2
            if (l.eq.ld+1) then
               l1 = -1
               l2 = -2
            end if
            do k=2,kd
               do j=2,jd
                  dx = xcc(l+l1,j,k) - xcc(l+l2,j,k)
                  dy = ycc(l+l1,j,k) - ycc(l+l2,j,k)
                  dz = zcc(l+l1,j,k) - zcc(l+l2,j,k)
                  xcc(l,j,k) = xcc(l+l1,j,k) + dx
                  ycc(l,j,k) = ycc(l+l1,j,k) + dy
                  zcc(l,j,k) = zcc(l+l1,j,k) + dz
               end do
            end do
         end do
c
c        exratapolate from interior in j-direction for j=const faces
c
         do j=1,jd+1,jd
            j1 = 1
            j2 = 2
            if (j.eq.jd+1) then
               j1 = -1
               j2 = -2
            end if
            do k=1,kd
               do l=1,ld
                  dx = xcc(l,j+j1,k) - xcc(l,j+j2,k)
                  dy = ycc(l,j+j1,k) - ycc(l,j+j2,k)
                  dz = zcc(l,j+j1,k) - zcc(l,j+j2,k)
                  xcc(l,j,k) = xcc(l,j+j1,k) + dx
                  ycc(l,j,k) = ycc(l,j+j1,k) + dy
                  zcc(l,j,k) = zcc(l,j+j1,k) + dz
               end do
            end do
         end do
c
c        exratapolate from interior in k-direction for k=const faces
c
         do k=1,kd+1,kd
            k1 = 1
            k2 = 2
            if (k.eq.kd+1) then
               k1 = -1
               k2 = -2
            end if
            do j=1,jd
               do l=1,ld
                  dx = xcc(l,j,k+k1) - xcc(l,j,k+k2)
                  dy = ycc(l,j,k+k1) - ycc(l,j,k+k2)
                  dz = zcc(l,j,k+k1) - zcc(l,j,k+k2)
                  xcc(l,j,k) = xcc(l,j,k+k1) + dx
                  ycc(l,j,k) = ycc(l,j,k+k1) + dy
                  zcc(l,j,k) = zcc(l,j,k+k1) + dz
               end do
            end do
         end do
c
c        fill in j-k edge values via extapolation in both j and k
c        directions and averaging
c
         do k=1,kd+1,kd
            k1 = 1
            k2 = 2
            if (k.eq.kd+1) then
               k1 = -1
               k2 = -2
            end if
            do j=1,jd+1,jd
                  j1 = 1
                  j2 = 2
               if (j.eq.jd+1) then
                  j1 = -1
                  j2 = -2
               end if
               do l=2,ld
                  dx1 = xcc(l,j+j1,k) - xcc(l,j+j2,k)
                  dy1 = ycc(l,j+j1,k) - ycc(l,j+j2,k)
                  dz1 = zcc(l,j+j1,k) - zcc(l,j+j2,k)
                  x1  = xcc(l,j+j1,k) + dx1
                  y1  = ycc(l,j+j1,k) + dy1
                  z1  = zcc(l,j+j1,k) + dz1
                  dx2 = xcc(l,j,k+k1) - xcc(l,j,k+k2)
                  dy2 = ycc(l,j,k+k1) - ycc(l,j,k+k2)
                  dz2 = zcc(l,j,k+k1) - zcc(l,j,k+k2)
                  x2  = xcc(l,j,k+k1) + dx2
                  y2  = ycc(l,j,k+k1) + dy2
                  z2  = zcc(l,j,k+k1) + dz2
                  xcc(l,j,k) = (x1 + x2)/2.0
                  ycc(l,j,k) = (y1 + y2)/2.0
                  zcc(l,j,k) = (z1 + z2)/2.0
               end do
            end do
         end do
c
c        fill in l-j edge values via extapolation in both l and j
c        directions and averaging
c
         do k=2,kd
            do j=1,jd+1,jd
               j1 = 1
               j2 = 2
               if (j.eq.jd+1)  then
                  j1 = -1
                  j2 = -2
               end if
               do l=1,ld+1,ld
                  l1 = 1
                  l2 = 2
                  if (l.eq.ld+1) then
                     l1 = -1
                     l2 = -2 
                  end if
                  dx1 = xcc(l,j+j1,k) - xcc(l,j+j2,k)
                  dy1 = ycc(l,j+j1,k) - ycc(l,j+j2,k)
                  dz1 = zcc(l,j+j1,k) - zcc(l,j+j2,k)
                  x1  = xcc(l,j+j1,k) + dx1
                  y1  = ycc(l,j+j1,k) + dy1
                  z1  = zcc(l,j+j1,k) + dz1
                  dx2 = xcc(l+l1,j,k) - xcc(l+l2,j,k)
                  dy2 = ycc(l+l1,j,k) - ycc(l+l2,j,k)
                  dz2 = zcc(l+l1,j,k) - zcc(l+l2,j,k)
                  x2  = xcc(l+l1,j,k) + dx2
                  y2  = ycc(l+l1,j,k) + dy2
                  z2  = zcc(l+l1,j,k) + dz2
                  xcc(l,j,k) = (x1 + x2)/2.0
                  ycc(l,j,k) = (y1 + y2)/2.0
                  zcc(l,j,k) = (z1 + z2)/2.0
               end do
            end do
         end do
c
c        fill in l-k edge values via extapolation in both l and k
c        directions and averaging
c
         do k=1,kd+1,kd
            k1 = 1
            k2 = 2
            if (k.eq.kd+1) then
               k1 = -1
               k2 = -2
            end if
            do j=2,jd
               do l=1,ld+1,ld
                  l1 = 1
                  l2 = 2
                  if (l.eq.ld+1) then
                     l1 = -1
                     l2 = -2
                  end if
                  dx1 = xcc(l+l1,j,k) - xcc(l+l2,j,k)
                  dy1 = ycc(l+l1,j,k) - ycc(l+l2,j,k)
                  dz1 = zcc(l+l1,j,k) - zcc(l+l2,j,k)
                  x1  = xcc(l+l1,j,k) + dx1
                  y1  = ycc(l+l1,j,k) + dy1
                  z1  = zcc(l+l1,j,k) + dz1
                  dx2 = xcc(l,j,k+k1) - xcc(l,j,k+k2)
                  dy2 = ycc(l,j,k+k1) - ycc(l,j,k+k2)
                  dz2 = zcc(l,j,k+k1) - zcc(l,j,k+k2)
                  x2  = xcc(l,j,k+k1) + dx2
                  y2  = ycc(l,j,k+k1) + dy2
                  z2  = zcc(l,j,k+k1) + dz2
                  xcc(l,j,k) = (x1 + x2)/2.0
                  ycc(l,j,k) = (y1 + y2)/2.0
                  zcc(l,j,k) = (z1 + z2)/2.0
               end do
            end do
         end do
c
c        fill in corner values via extapolation in all 3 directions
c        and averaging
c
         do k=1,kd+1,kd
            k1 = 1
            k2 = 2
            if (k.eq.kd+1) then
               k1 = -1
               k2 = -2
            end if
            do j=1,jd+1,jd
               j1 = 1
               j2 = 2
               if (j.eq.jd+1) then
                  j1 = -1
                  j2 = -2
               end if
               do l=1,ld+1,ld
                  l1 = 1
                  l2 = 2
                  if (l.eq.ld+1) then
                     l1 = -1
                     l2 = -2
                  end if         
                  dx1 = xcc(l+l1,j,k) - xcc(l+l2,j,k)
                  dy1 = ycc(l+l1,j,k) - ycc(l+l2,j,k)
                  dz1 = zcc(l+l1,j,k) - zcc(l+l2,j,k)
                  x1  = xcc(l+l1,j,k) + dx1
                  y1  = ycc(l+l1,j,k) + dy1
                  z1  = zcc(l+l1,j,k) + dz1
                  dx2 = xcc(l,j+j1,k) - xcc(l,j+j2,k)
                  dy2 = ycc(l,j+j1,k) - ycc(l,j+j2,k)
                  dz2 = zcc(l,j+j1,k) - zcc(l,j+j2,k)
                  x2  = xcc(l,j+j1,k) + dx2
                  y2  = ycc(l,j+j1,k) + dy2
                  z2  = zcc(l,j+j1,k) + dz2
                  dx3 = xcc(l,j,k+k1) - xcc(l,j,k+k2)
                  dy3 = ycc(l,j,k+k1) - ycc(l,j,k+k2)
                  dz3 = zcc(l,j,k+k1) - zcc(l,j,k+k2)
                  x3  = xcc(l,j,k+k1) + dx3
                  y3  = ycc(l,j,k+k1) + dy3
                  z3  = zcc(l,j,k+k1) + dz3
                  xcc(l,j,k) = (x1 + x2 + x3)/3.
                  ycc(l,j,k) = (y1 + y2 + y3)/3.
                  zcc(l,j,k) = (z1 + z2 + z3)/3.
               end do
            end do
         end do
c
      end if
c
      return 
      end
