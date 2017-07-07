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
      program everyother_xyz
c
c     $Id$
c
c***********************************************************************
c     Purpose: Reads grid and makes <every other> point grid.
c     User can choose to coarsen only particular index directions,
c     if desired.  Does not keep iblanking information if it is present.
c
c     If using this alone (not in conjunction with cfl3d makefile):
c     f90 -64 -r8 everyother_xyz.f umalloc_r.o -o everyother_xyz
c***********************************************************************
c
      parameter (nbb=10000)
      dimension idim(nbb),jdim(nbb),kdim(nbb)
      character file1*80
c
      write(6,'('' input grid to read from:'')')
      read(5,'(a80)') file1
c
      write(6,'(''enter 0 to convert plot3d-type MG grid'')')
      write(6,'(''enter 1 to convert  cfl3d-type grid'')')
      read(5,*) ip3d
      iblank=0
      if (ip3d .eq. 0) then
        write(6,'('' input 1 if grid is iblanked:'')')
        read(5,*) iblank
      end if
c
      write(6,'(''enter 0 to read formatted file'')')
      write(6,'(''enter 1 to read unformatted file'')')
      read(5,*) iunfi
      if (iunfi .eq. 1) then
        open(2,file=file1,form='unformatted',status='old')
      else
        open(2,file=file1,form='formatted',status='old')
      end if
c
      write(6,'('' input 3 numbers (i,j,k-directions) to'',
     . '' indicate which directions to coarsen'')')
      write(6,'(''   e.g., 1,1,1 means coarsen in all 3'',
     . '' directions,'')')
      write(6,'(''   e.g., 0,1,0 means coarsen in j only:'')')
      read(5,*) ic,jc,kc
      if (ic .lt. 0 .or. ic .gt. 1) ic=0
      if (jc .lt. 0 .or. jc .gt. 1) jc=0
      if (kc .lt. 0 .or. kc .gt. 1) kc=0
      if (ic .eq. 1) then
        write(6,'('' coarsening in i-direction (if idim=2'',
     +   '' it is assumed to be'')')
        write(6,'(''     2-D and is not coarsened)'')')
      end if
      if (jc .eq. 1) then
        write(6,'('' coarsening in j-direction'')')
      end if
      if (kc .eq. 1) then
        write(6,'('' coarsening in k-direction'',/)')
      end if
c
      write(6,'(''enter 0 to write formatted file'')')
      write(6,'(''enter 1 to write unformatted file'')')
      read(5,*) iunfo
      if (iunfo .eq. 1) then
        open(3,file='grid.new',form='unformatted',status='unknown')
      else
        open(3,file='grid.new',form='formatted',status='unknown')
      end if
c
      if (ip3d .eq. 0) then
        if (iunfi .eq. 1) then
          read(2) nbl
        else
          read(2,*) nbl
        end if
        if(nbl .gt. nbb) then
          write(6,'('' need to increase nbb to '',i6)') nbl
          stop
        end if
        write(6,'('' Grid has '',i5,'' zones.'')') nbl
        if (iunfo .eq. 1) then
          write(3) nbl
        else
          write(3,*) nbl
        end if
        if (iunfi .eq. 1) then
          read(2) (idim(n),jdim(n),kdim(n),n=1,nbl)
        else
          read(2,*) (idim(n),jdim(n),kdim(n),n=1,nbl)
        end if
      else
        nbl=0
        do n=1,nbb-1
          if (iunfi .eq. 1) then
            read(2,end=9997) jdim(n),kdim(n),idim(n)
          else
            read(2,*,end=9997) jdim(n),kdim(n),idim(n)
          end if
          nbl=nbl+1
          if (iunfi .eq. 1) then
            if (iblank .eq. 1) then
            read(2) (xx,ll=1,jdim(n)*kdim(n)*idim(n)),
     .              (yy,ll=1,jdim(n)*kdim(n)*idim(n)),
     .              (zz,ll=1,jdim(n)*kdim(n)*idim(n)),
     .              (ibl,ll=1,jdim(n)*kdim(n)*idim(n))
            else
            read(2) (xx,ll=1,jdim(n)*kdim(n)*idim(n)),
     .              (yy,ll=1,jdim(n)*kdim(n)*idim(n)),
     .              (zz,ll=1,jdim(n)*kdim(n)*idim(n))
            end if
          else
            if (iblank .eq. 1) then
            read(2,*) (xx,ll=1,jdim(n)*kdim(n)*idim(n)),
     .                (yy,ll=1,jdim(n)*kdim(n)*idim(n)),
     .                (zz,ll=1,jdim(n)*kdim(n)*idim(n)),
     .                (ibl,ll=1,jdim(n)*kdim(n)*idim(n))
            else
            read(2,*) (xx,ll=1,jdim(n)*kdim(n)*idim(n)),
     .                (yy,ll=1,jdim(n)*kdim(n)*idim(n)),
     .                (zz,ll=1,jdim(n)*kdim(n)*idim(n))
            end if
          end if
        end do
        write(6,'('' Error, nbb not big enough'')')
        stop
 9997   continue
        write(6,'('' Grid has '',i5,'' zones.'')') nbl
        if(nbl .gt. nbb) then
          write(6,'('' need to increase nbb to '',i6)') nbl
          stop
        end if
        rewind(2)
      end if
c
      maxbl = nbl
      imax  = 1
      jmax  = 1
      kmax  = 1
      do n=1,nbl
         write(6,'('' read in zone '',i5,'', i,j,k='',3i6)')
     +     n,idim(n),jdim(n),kdim(n)
         imax  = max(imax,idim(n))
         jmax  = max(jmax,jdim(n))
         kmax  = max(kmax,kdim(n))
      end do
c
      call convert(maxbl,imax,jmax,kmax,nbl,idim,jdim,kdim,
     .             iunfi,iunfo,ip3d,nbb,iblank,ic,jc,kc)
c
      write(6,'('' successful termination'')')
c
      stop
      end
c
c  **************************************************************
      subroutine convert(maxbl,imax,jmax,kmax,nbl,idim,jdim,kdim,
     .                   iunfi,iunfo,ip3d,nbb,iblank,ic,jc,kc)
c
      integer stats
c
      dimension idim(nbb)
      dimension idim2(nbb)
      dimension jdim(nbb)
      dimension jdim2(nbb)
      dimension kdim(nbb)
      dimension kdim2(nbb)
c
      allocatable :: x(:,:,:)
      allocatable :: y(:,:,:)
      allocatable :: z(:,:,:)
c
c     allocate memory
c
      memuse = 0
      allocate( x(imax,jmax,kmax), stat=stats )
      call umalloc_r(imax*jmax*kmax,0,'x',memuse,stats)
      allocate( y(imax,jmax,kmax), stat=stats )
      call umalloc_r(imax*jmax*kmax,0,'y',memuse,stats)
      allocate( z(imax,jmax,kmax), stat=stats )
      call umalloc_r(imax*jmax*kmax,0,'z',memuse,stats)

c
      do n=1,nbl
      if(jc .eq. 1) then
      if(jdim(n)/2 .eq. (jdim(n)+1)/2) then
        write(6,'('' jdim must be odd'')')
        stop
      end if
      jdim2(n)=(jdim(n)+1)/2
      else
      jdim2(n)=jdim(n)
      end if
      if(kc .eq. 1) then
      if(kdim(n)/2 .eq. (kdim(n)+1)/2) then
        write(6,'('' kdim must be odd'')')
        stop
      end if
      kdim2(n)=(kdim(n)+1)/2
      else
      kdim2(n)=kdim(n)
      end if
      if (idim(n) .eq. 2) then
        idim2(n)=2
      else
        if(ic .eq. 1) then
        if(idim(n)/2 .eq. (idim(n)+1)/2) then
          write(6,'('' idim must be odd'')')
          stop
        end if
        idim2(n)=(idim(n)+1)/2
        else
        idim2(n)=idim(n)
        end if
      end if
      write(6,'('' Zone #'',i5,'' output idim,jdim,kdim='',3i6)') 
     +   n,idim2(n),jdim2(n),kdim2(n)
      enddo
      if (ip3d .eq. 0) then
        if (iunfo .eq. 1) then
          write(3) (idim2(n),jdim2(n),kdim2(n),n=1,nbl)
        else
          write(3,*) (idim2(n),jdim2(n),kdim2(n),n=1,nbl)
        end if
      end if
c
c   Do loop over all the zones
      do n=1,nbl
      ni=idim(n)
      nj=jdim(n)
      nk=kdim(n)
      if (ip3d .eq. 0) then
        if (iunfi .eq. 1) then
          if (iblank .eq. 1) then
          read(2) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .            (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .            (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .            (((ibl,i=1,ni),j=1,nj),k=1,nk)
          else
          read(2) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .            (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .            (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
          end if
        else
          if (iblank .eq. 1) then
          read(2,*) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .              (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .              (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .              (((ibl,i=1,ni),j=1,nj),k=1,nk)
          else
          read(2,*) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .              (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     .              (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
          end if
        end if
        if (ni .eq. 2) then
        if (iunfo .eq. 1) then
          write(3) (((x(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1),
     +             (((y(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1),
     +             (((z(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1)
        else
          write(3,*) (((x(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1),
     +               (((y(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1),
     +               (((z(i,j,k),i=1,ni),j=1,nj,jc+1),k=1,nk,kc+1)
        end if
        else
        if (iunfo .eq. 1) then
          write(3) (((x(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1),
     +             (((y(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1),
     +             (((z(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1)
        else
          write(3,*) (((x(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1),
     +               (((y(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1),
     +               (((z(i,j,k),i=1,ni,ic+1),j=1,nj,jc+1),k=1,nk,kc+1)
        end if
        end if
      else
        if (iunfi .eq. 1) then
          read(2) jdum,kdum,idum
          read(2) (((x(i,j,k),j=1,nj),k=1,nk),i=1,ni),
     .            (((y(i,j,k),j=1,nj),k=1,nk),i=1,ni),
     .            (((z(i,j,k),j=1,nj),k=1,nk),i=1,ni)
        else
          read(2,*) jdum,kdum,idum
          read(2,*) (((x(i,j,k),j=1,nj),k=1,nk),i=1,ni),
     .              (((y(i,j,k),j=1,nj),k=1,nk),i=1,ni),
     .              (((z(i,j,k),j=1,nj),k=1,nk),i=1,ni)
        end if
        if (iunfo .eq. 1) then
          write(3) jdim2(n),kdim2(n),idim2(n)
        else
          write(3,*) jdim2(n),kdim2(n),idim2(n)
        end if
        if (ni .eq. 2) then
        if (iunfo .eq. 1) then
          write(3) (((x(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni),
     +             (((y(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni),
     +             (((z(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni)
        else
          write(3,*) (((x(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni),
     +               (((y(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni),
     +               (((z(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni)
        end if
        else
        if (iunfo .eq. 1) then
          write(3) (((x(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1),
     +             (((y(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1),
     +             (((z(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1)
        else
          write(3,*) (((x(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1),
     +               (((y(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1),
     +               (((z(i,j,k),j=1,nj,jc+1),k=1,nk,kc+1),i=1,ni,ic+1)
        end if
        end if
      end if
      enddo
c
      if (ip3d .eq. 0) then
        write(6,'(/,'' New grid to grid.new, plot3d format'')')
        if (iunfo .eq. 1) then
          write(6,'('' ...file is unformatted, MG, no iblanking'')')
        else
          write(6,'('' ...file is formatted, MG, no iblanking'')')
        end if
      else
        write(6,'(/,'' New grid to grid.new, cfl3d format'')')
        if (iunfo .eq. 1) then
          write(6,'('' ...file is unformatted'')')
        else
          write(6,'('' ...file is formatted'')')
        end if
      end if
c
c     free memory
c
      deallocate(x)
      deallocate(y)
      deallocate(z)
c
      return
      end
