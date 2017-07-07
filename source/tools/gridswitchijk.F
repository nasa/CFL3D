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
      program gridswitchijk
c
c***********************************************************************
c     Purpose:
c     Reads a PLOT3D grid and 
c     switches j,k,i to i,j,k (also transposes y and z if desired).
c     This is typically helpful if you have a 2-plane grid with
c     k as the 2-D direction (CFL3D wants it to be i).
c
c     If using this alone (not in conjunction with cfl3d makefile):
c     f90 -64 -r8 gridswitchijk.f umalloc_r.o -o gridswitchijk
c***********************************************************************
c
      allocatable :: jt(:),kt(:),it(:)
      allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
      character*80 file1
c 
      write(6,'('' input existing unformatted plot3d grid name:'')')
      read(5,'(a80)') file1
      open(7,file=file1,form='unformatted',status='old')
      open(8,file='tempijk.p3dbin',form='unformatted',status='unknown')
c
      write(6,'('' reading grid'')')
      read(7) nbin
      write(8) nbin
c   allocate memory for jt, kt, it
      memuse = 0
      allocate( jt(nbin), stat=istats )
      call umalloc_r(nbin,0,'jt',memuse,istats)
      allocate( kt(nbin), stat=istats )
      call umalloc_r(nbin,0,'kt',memuse,istats)
      allocate( it(nbin), stat=istats )
      call umalloc_r(nbin,0,'it',memuse,istats)
c
      read(7) (it(n),jt(n),kt(n),n=1,nbin)
      write(8) (kt(n),it(n),jt(n),n=1,nbin)
      itmax=0
      jtmax=0
      ktmax=0
      do n=1,nbin
        if (it(n) .gt. itmax) itmax=it(n)
        if (jt(n) .gt. jtmax) jtmax=jt(n)
        if (kt(n) .gt. ktmax) ktmax=kt(n)
      enddo
      write(6,'('' itmax,jtmax,ktmax='',3i6)') itmax,jtmax,ktmax
c   allocate memory for x, y, z
      memuse = 0
      allocate( x(ktmax,itmax,jtmax), stat=istats )
      call umalloc_r(ktmax*itmax*jtmax,0,'x',memuse,istats)
      allocate( y(ktmax,itmax,jtmax), stat=istats )
      call umalloc_r(ktmax*itmax*jtmax,0,'y',memuse,istats)
      allocate( z(ktmax,itmax,jtmax), stat=istats )
      call umalloc_r(ktmax*itmax*jtmax,0,'z',memuse,istats)
c
      write(6,'('' transpose y and z? (1=yes)'')')
      write(6,'(''   ... do this, e.g., if y is up and you'',
     +   '' prefer z to be up:'')')
      read(5,*) itrans
      write(6,'('' writing new grid tempijk.p3dbin'')')
c
      do n=1,nbin
      read(7) (((x(k,i,j),i=1,it(n)),j=1,jt(n)),k=1,kt(n)),
     +        (((y(k,i,j),i=1,it(n)),j=1,jt(n)),k=1,kt(n)),
     +        (((z(k,i,j),i=1,it(n)),j=1,jt(n)),k=1,kt(n))
c
      if (itrans .eq. 1) then
      write(8) (((x(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n)),
     +         (((-z(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n)),
     +         (((y(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n))
      else
      write(8) (((x(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n)),
     +         (((y(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n)),
     +         (((z(k,i,j),k=1,kt(n)),i=1,it(n)),j=1,jt(n))
      end if
      enddo
c
      deallocate(it,jt,kt)
      deallocate(x,y,z)
c
      write(6,'('' ijk reversed'')')
      if (itrans .eq. 1) then
        write(6,'('' ... y and z transposed also'')')
      end if
      stop
      end
