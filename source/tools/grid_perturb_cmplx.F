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
      program grid_perturb_cmplx
c
c     $Id$
c
c***********************************************************************
c     Purpose: Generate a complex-valued grid (plot3d multiblock form)
c     by reading in a real-valued grid (plot3d multiblock form) and a
c     corresponding real-valued matrix of grid-sensitivity derivatives
c     (plot3d multiblock function file form, with 3*ndv variables for
c     the x,y,z components of the ndv design variables). The data in the
c     sensitivity file should correspond to d(x)/d(DV),d(y)/d(DV),
c     d(z)/d(DV) of the grid, for each of the ndv design variables. The
c     user chooses one of the design variables and a step size, eps. The
c     new grid is created by adding an imaginary component to the input 
c     grid via:
c             (x,y,z) = (x,y,z) +i*eps*d(x,y,z)/d(DV)  where i=sqrt(-1)
c
c     The output grid may be read into the *complex version* of cfl3d
c     (cfl3dcmplx_mpi or cfl3dcmplx_seq) to determine the solution
c     derivatives with respect to the chosen design variable.
c*********************************************************************** 
c
      implicit complex(a-h,o-z)
c
      real reps
c
      character*80 gfilein,gfileout,sdfilein
      common /files/gfilein,gfileout,sdfilein
c
      write(6,*)'input name of baseline grid file'
      read(5,'(a80)') gfilein
      write(6,*)'unformatted or formatted (unform = 0; form=1)'
      read(5,*) iginform
      write(6,*)'input name of grid sensitivity file to use'
      read(5,'(a80)') sdfilein
      write(6,*)'unformatted or formatted (unform = 0; form=1)'
      read(5,*) isdform
      write(6,*)'input name of perturbed grid to create'
      read(5,'(a80)') gfileout
      write(6,*)'unformatted or formatted (unform = 0; form=1)'
      read(5,*) igoutform
      write(6,*)'input DV number to base perturbation on'
      read(5,*) idv
      write(6,*)'input step size for design variable ',idv
      read(5,*) reps
      eps = reps
c
c     check for input and output files named the same
c
      if (gfileout .eq. gfilein)then
         write(6,*)'stopping input and perturbed grid files ',
     .             'are the same'
         stop
      end if
      if (gfileout .eq. sdfilein) then
         write(6,*)'stopping grid sensitivity and perturbed ',
     .             'grid files are the same'
         stop
      end if
c
c     read files to determine array sizes
c
      if (iginform.eq.0) then
         open(unit=10,file=gfilein,form='unformatted',status='old')
      else
         open(unit=10,file=gfilein,form='formatted',status='old')
      end if
      if (isdform.eq.0) then
         open(unit=11,file=sdfilein,form='unformatted',status='old')
      else
         open(unit=11,file=sdfilein,form='formatted',status='old')
      end if
c
      if (iginform.eq.0) then
         read(10) ngrd
      else
         read(10,*) ngrd
      end if
      maxbl = ngrd
c
c     determine array i,j,k sizes
c
      call sizeijk(maxbl,imax,jmax,kmax,iginform)
c
      if (isdform.eq.0) then
         read(11) ngrdsd
      else
         read(11,*) ngrdsd
      end if
      if (ngrd.ne.ngrdsd) then
         write(6,*)'stopping...grid and SD files contain ',
     .             'different number of zones'
         write(6,*)'zones in grid file: ',ngrd
         write(6,*)'zones in  SD  file: ',ngrdsd
         stop
      end if
      if (isdform.eq.0) then
         read(11) (idum,jdum,kdum,ndvx3,m=1,ngrdsd)
      else
         read(11,*) (idum,jdum,kdum,ndvx3,m=1,ngrdsd)
      end if
      ndv   = ndvx3/3
      ndvmx = ndv
c
      write(6,*)
      write(6,'(''required array dimensions:'')')
      write(6,'(''  maxbl = '',i3)') maxbl
      write(6,'(''  imax  = '',i3)') imax
      write(6,'(''  jmax  = '',i3)') jmax
      write(6,'(''  kmax  = '',i3)') kmax
      write(6,'(''  ndvmx = '',i3)') ndvmx
c
c     rewind and close and input files
c
      rewind(10)
      rewind(11)
      close(10)
      close(11)
c
      call mkgrd(jmax,kmax,imax,maxbl,ndvmx,eps,iginform,
     .           igoutform,isdform,idv)
c
      stop
      end
c 
      subroutine sizeijk(maxbl,imax,jmax,kmax,iginform)
c***********************************************************************
c     Purpose: determine max i,j,k dimensions needed for this case
c***********************************************************************
c
      implicit complex(a-h,o-z)
c
      integer stats

      allocatable :: itemp(:)
      allocatable :: jtemp(:)
      allocatable :: ktemp(:)
c
c     allocate memory
c
      memuse = 0
      allocate( itemp(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'itemp',memuse,stats)
      allocate( jtemp(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'jtemp',memuse,stats)
      allocate( ktemp(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'ktemp',memuse,stats)
c
      if (iginform.eq.0) then
         read(10) (itemp(n),jtemp(n),ktemp(n),n=1,maxbl)
      else
         read(10,*) (itemp(n),jtemp(n),ktemp(n),n=1,maxbl)
      end if
c
      imax = itemp(1)
      jmax = jtemp(1)
      kmax = ktemp(1)
      do n=1,maxbl
         if (itemp(n).gt.imax) imax=itemp(n)
         if (jtemp(n).gt.jmax) jmax=jtemp(n)
         if (ktemp(n).gt.kmax) kmax=ktemp(n)
      end do
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
      subroutine mkgrd(jmax,kmax,imax,maxbl,ndvmx,eps,iginform,
     .                 igoutform,isdform,idv)
c***********************************************************************
c     Purpose: generate a new, perturbed grid from the baseline grid
c     and the grid sensitivities
c***********************************************************************
c
      implicit complex(a-h,o-z)
c
      integer stats

      real, allocatable :: g_x(:,:,:,:)
      real, allocatable :: g_y(:,:,:,:)
      real, allocatable :: g_z(:,:,:,:)
      allocatable :: idim(:)
      allocatable :: idimsd(:)
      allocatable :: jdim(:)
      allocatable :: jdimsd(:)
      allocatable :: kdim(:)
      allocatable :: kdimsd(:)
      real, allocatable :: x1(:,:,:)
      allocatable :: x2(:,:,:)
      real, allocatable :: y1(:,:,:)
      allocatable :: y2(:,:,:)
      real, allocatable :: z1(:,:,:)
      allocatable :: z2(:,:,:)

      character*80 gfilein,gfileout,sdfilein
      common /files/gfilein,gfileout,sdfilein
c
      sqrtm1 = cmplx(0.,1.)
      memuse = 0
c
c     allocate memory
c
      allocate( g_x(ndvmx,jmax,kmax,imax), stat=stats )
      call umalloc_r(ndvmx*jmax*kmax*imax,0,'g_x',memuse,stats)
      allocate( g_y(ndvmx,jmax,kmax,imax), stat=stats )
      call umalloc_r(ndvmx*jmax*kmax*imax,0,'g_y',memuse,stats)
      allocate( g_z(ndvmx,jmax,kmax,imax), stat=stats )
      call umalloc_r(ndvmx*jmax*kmax*imax,0,'g_z',memuse,stats)
      allocate( idim(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'idim',memuse,stats)
      allocate( idimsd(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'idimsd',memuse,stats)
      allocate( jdim(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'jdim',memuse,stats)
      allocate( jdimsd(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'jdimsd',memuse,stats)
      allocate( kdim(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'kdim',memuse,stats)
      allocate( kdimsd(maxbl), stat=stats )
      call umalloc_c(maxbl,1,'kdimsd',memuse,stats)
      allocate( x1(jmax,kmax,imax), stat=stats )
      call umalloc_r(jmax*kmax*imax,0,'x1',memuse,stats)
      allocate( x2(jmax,kmax,imax), stat=stats )
      call umalloc_c(jmax*kmax*imax,0,'x2',memuse,stats)
      allocate( y1(jmax,kmax,imax), stat=stats )
      call umalloc_r(jmax*kmax*imax,0,'y1',memuse,stats)
      allocate( y2(jmax,kmax,imax), stat=stats )
      call umalloc_c(jmax*kmax*imax,0,'y2',memuse,stats)
      allocate( z1(jmax,kmax,imax), stat=stats )
      call umalloc_r(jmax*kmax*imax,0,'z1',memuse,stats)
      allocate( z2(jmax,kmax,imax), stat=stats )
      call umalloc_c(jmax*kmax*imax,0,'z2',memuse,stats)
c
c     open files
c
      if (iginform.eq.0) then 
         open(unit=10,file=gfilein,form='unformatted',status='old')
      else
         open(unit=10,file=gfilein,form='formatted',status='old')
      end if
      if (isdform.eq.0) then
         open(unit=11,file=sdfilein,form='unformatted',status='old')
      else
         open(unit=11,file=sdfilein,form='formatted',status='old')
      end if
      if (igoutform.eq.0) then
         open(unit=12,file=gfileout,form='unformatted',status='unknown')
      else
         open(unit=12,file=gfileout,form='formatted',status='unknown')
      end if
c
c     read grid file dimensions and check parameters
c
      if (iginform.eq.0) then
         read(10) ngrd
      else
         read(10,*) ngrd
      end if
      if (ngrd.gt.maxbl) then
         write(6,*)'increase parameter maxbl to ',ngrd
         stop
      end if
      if (iginform.eq.0) then
         read(10) (idim(n),jdim(n),kdim(n),n=1,ngrd)
      else
         read(10,*) (idim(n),jdim(n),kdim(n),n=1,ngrd)
      end if
      istop = 0
      jstop = 0
      kstop = 0
      do n=1,ngrd
         if (idim(n).gt.imax) then
            ineed = idim(n)
            istop = 1
         end if
         if (jdim(n).gt.jmax) then
            jneed = jdim(n)
            jstop = 1
         end if
         if (kdim(n).gt.kmax) then
            kneed = kdim(n)
            kstop = 1
         end if
      end do
      if (istop.gt.0) then
         write(6,*)'stopping...increase parameter imax to ',ineed
      end if
      if (jstop.gt.0) then
         write(6,*)'stopping...increase parameter jmax to ',jneed
      end if
      if (kstop.gt.0) then
         write(6,*)'stopping...increase parameter jmax to ',kneed
      end if
      if (istop.gt.0 .or. jstop.gt.0 .or.
     .    kstop.gt.0) then
          stop
      end if
c
c     read sensitivity file dimensions and make sure they match the 
c     grid dimensions
c
      if (isdform.eq.0) then
         read(11) ngrdsd
      else
         read(11,*) ngrdsd
      end if
      if (ngrd.ne.ngrdsd) then 
         write(6,*)'stopping...grid and SD files contain ',
     .             'different number of zones'
         write(6,*)'zones in grid file: ',ngrd
         write(6,*)'zones in  SD  file: ',ngrdsd
         stop
      end if
      if (isdform.eq.0) then
         read(11) (idimsd(m),jdimsd(m),kdimsd(m),ndvx3,m=1,ngrdsd)
      else
         read(11,*) (idimsd(m),jdimsd(m),kdimsd(m),ndvx3,m=1,ngrdsd)
      end if
      ndv = ndvx3/3
      do m=1,ngrdsd
         if (idim(m).ne.idimsd(m) .or. jdim(m).ne.jdimsd(m) .or.
     .       kdim(m).ne.kdimsd(m)) then
            write(6,*)'stopping...grid and SD files have ',
     .                'different dimensions in zone',m
            write(6,*)'grid file i,j,k: ',idim(m),jdim(m),kdim(m)
            write(6,*)'SD   file i,j,k: ',idimsd(m),jdimsd(m),
     .                kdimsd(m)
            stop
         end if
      end do
c
c     write dimensions for perturbed grid file
c 
      if (igoutform.eq.0) then
         write(12) ngrd
         write(12) (idim(n),jdim(n),kdim(n),n=1,ngrd)
      else
         write(12,*) ngrd
         write(12,*) (idim(n),jdim(n),kdim(n),n=1,ngrd)
      end if
c
      do n=1,ngrd
c
         jd = jdim(n)
         kd = kdim(n)
         id = idim(n)
         call rp3d(x1,y1,z1,jmax,kmax,imax,jd,kd,id,iginform)
         call g_rp3d(g_x,g_y,g_z,ndvmx,jmax,kmax,imax,
     .               jd,kd,id,ndv,idv,isdform)
c
         do j=1,jd
         do k=1,kd
         do i=1,id
            x2(j,k,i) = x1(j,k,i) + sqrtm1*eps*g_x(idv,j,k,i)
            y2(j,k,i) = y1(j,k,i) + sqrtm1*eps*g_y(idv,j,k,i)
            z2(j,k,i) = z1(j,k,i) + sqrtm1*eps*g_z(idv,j,k,i)
         end do
         end do
         end do
c
      if (igoutform.eq.0) then
         write(12) (((x2(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .             (((y2(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .             (((z2(j,k,i),i=1,id),j=1,jd),k=1,kd)
      else
         write(12,*) (((x2(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .               (((y2(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .               (((z2(j,k,i),i=1,id),j=1,jd),k=1,kd)
      end if
c
      end do
c
c     free memory
c
      deallocate(idim)
      deallocate(jdim)
      deallocate(kdim)
      deallocate(idimsd)
      deallocate(jdimsd)
      deallocate(kdimsd)
      deallocate(x1)
      deallocate(y1)
      deallocate(z1)
      deallocate(x2)
      deallocate(y2)
      deallocate(z2)
      deallocate(g_x)
      deallocate(g_y)
      deallocate(g_z)
c
      return
      end
c
      subroutine rp3d(x,y,z,jdim,kdim,idim,jd,kd,id,iginform)
c***********************************************************************
c     Purpose: read grids in plot3d format.
c     iginform - flag for grid file type
c              = 0 grid file is unformatted
c              > 0 grid file is formatted
c***********************************************************************
c
      implicit complex(a-h,o-z)
c
      real x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      if (iginform.eq.0) then
         read(10) (((x(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .            (((y(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .            (((z(j,k,i),i=1,id),j=1,jd),k=1,kd)
      else
         read(10,*) (((x(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .              (((y(j,k,i),i=1,id),j=1,jd),k=1,kd),
     .              (((z(j,k,i),i=1,id),j=1,jd),k=1,kd)
      end if
c
      return
      end
c
      subroutine g_rp3d(g_x,g_y,g_z,ndvmx,jdim,kdim,idim,
     .                  jd,kd,id,ndv,idv,isdform)
c***********************************************************************
c     Purpose: read sensitivity derivative file in plot3d function
c     format for ndv design variables (the function file must contain
c     3*ndv functions for the x,y,z components of ndv design variables).
c     iginform - flag for grid file type
c              = 0 grid file is unformatted
c              > 0 grid file is formatted
c***********************************************************************
c
      implicit complex(a-h,o-z)
c
      real g_x(ndvmx,jdim,kdim,idim),g_y(ndvmx,jdim,kdim,idim),
     .     g_z(ndvmx,jdim,kdim,idim)
c
      if (isdform.eq.0) then
         read(11)((((g_x(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),
     .            (((g_y(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),
     .            (((g_z(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),nn=1,ndv)
      else
         read(11,*)((((g_x(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),
     .              (((g_y(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),
     .              (((g_z(nn,j,k,i),i=1,id),j=1,jd),k=1,kd),nn=1,ndv)
      end if
c
      return
      end
