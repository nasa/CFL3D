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
      program plot3d_to_cgns
c
c     $Id$
c
c***********************************************************************
c   Purpose: reads in plot3d-type (or cfl3d-type) unformatted grid file 
c   and writes out a CGNS ADF file.  Also reads in as input the 1-to-1
c   section of CFL3D input file, for creating the 1-to-1 CGNS
c   connectivity information.
c   Currently, you must be linked to CGNS V2.5.2 or later
c
c   This does a 3-D grid only, in MG format
c***********************************************************************
c
#ifdef CGNS
      character*80 file2,file3
      parameter (ibufdim0=2000,nbuf0=4,mxbcfil0=100)
c
      character*120 bou(ibufdim0,nbuf0)
c
      dimension nou(nbuf0)
c
      common /unit5/ iunit5
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /zero/ iexp
      common /mydist2/ nnodes,myhost,myid,mycomm
c
c  read the CFL3D input file to determine parameter requirements
c
      write(6,'('' Note:  currently, you must be linked to CGNS'',
     +  '' V2.5.2 or later'',/)')
      write(6,'('' What is the CFL3D input file name?'')')
      write(6,'(''   (be sure that IREST=0 in this file!!!)'')')
      read(5,'(a80)') file2
      open(3,file=file2,form='formatted',status='old')
      write(6,'('' Does this case use a patch file (e.g., '',
     +  ''patch.bin)? (1=yes)'')')
      read(5,*) ipatch
      if (ipatch .eq. 1) then
        write(6,'('' input patch filename:'')')
        read(5,'(a80)') file3
        open(22,file=file3,form='unformatted',status='old')
      end if
c
c     determine machine zero for use in setting tolerances
c     (10.**(-iexp) is machine zero)
c
      icount = 0
      compare = 1.0
      do i = 1,20
         icount = icount + 1
         add = 1.
         do n=1,i
            add = add*.1
         enddo
         x11 = compare + add
         if (x11.eq.compare)then
            iexp = i-1
            goto 4010
         end if
      end do
 4010 continue
c
c     set dummy value of xmach since subroutine readkey uses it
c     (but is not actually needed for this code)
c
      xmach = 0.5
c
      iunit5  = 3
      iunit11 = 99
      myid    = 0
      nnodes  = 1
      call global0(nplots0,maxnode0,n11,lbcprd0,lbcemb0,
     .             lbcrad0,maxbl0,nblock,maxseg0,maxcs0,ncycmax0,
     .             intmax0,nsub10,intmx0,mxxe0,mptch0,msub10,
     .             ibufdim0,nbuf0,mxbcfil0,nmds0,maxaes0,
     .             maxsegdg0,ntr,nnodes,nou,bou,iunit11,myid,
     .             idm,jdm,kdm)
c
      write(6,*)
      write(6,'(''required array dimensions:'')')
      write(6,'(''  nblock = '',i3)') nblock
      write(6,'(''  idm    = '',i3)') idm
      write(6,'(''  jdm    = '',i3)') jdm
      write(6,'(''  kdm    = '',i3)') kdm
      write(6,'(''  n11    = '',i3)') n11  
c
c     rewind and close cfl3d nput file
c
      rewind(3)
      close(3)
c
      call wrt_adf(file2,nblock,idm,jdm,kdm,n11)
c
      stop
      end
c
      subroutine wrt_adf(file2,nblock,idm,jdm,kdm,n11)
c***********************************************************************
c   Purpose: Write out a CGNS ADF file
c***********************************************************************
c
c
      integer stats
c
      dimension data(7)
      dimension ipnts(3,2)
      dimension irange(3,2)
      dimension iranged(3,2)
      dimension isize(9)
      dimension itransform(3)
c
      allocatable :: i0(:)
      allocatable :: iblk1(:)
      allocatable :: iblk2(:)
      allocatable :: idim(:)
      allocatable :: ihi1(:)
      allocatable :: ihi2(:)
      allocatable :: ilo1(:)
      allocatable :: ilo2(:)
      allocatable :: im(:)
      allocatable :: isva1(:)
      allocatable :: isva2(:)
      allocatable :: isvb1(:)
      allocatable :: isvb2(:)
      allocatable :: iv1(:)
      allocatable :: iv2(:)
      allocatable :: iv3(:)
      allocatable :: iz(:)
      allocatable :: j0(:)
      allocatable :: jdim(:)
      allocatable :: jhi1(:)
      allocatable :: jhi2(:)
      allocatable :: jlo1(:)
      allocatable :: jlo2(:)
      allocatable :: jm(:)
      allocatable :: k0(:)
      allocatable :: kdim(:)
      allocatable :: khi1(:)
      allocatable :: khi2(:)
      allocatable :: klo1(:)
      allocatable :: klo2(:)
      allocatable :: km(:)
      allocatable :: x(:,:,:)
      allocatable :: y(:,:,:)
      allocatable :: z(:,:,:)
c
c     real*8 pidz,pidnew
c
      character*80 file1,file2,name,title,textf
      character*32 filename,basename,filenum
      character*8 zonename(nblock)
      character*32 filenamr,connectname,filecon,donorname
      character*32 boconame,boconame2,boconame3
      character*6   texta
      character*21  text
      character*120 newname
      character*1 dum1
c
#     include "cgnslib_f.h"
c
c   Determine if single or double precision is being used:
      idouble=0
#if defined DBLE_PRECSN
      idouble=1
#endif
c
c     allocate memory
c
      memuse = 0
      allocate( i0(nblock), stat=stats )
      call umalloc_r(nblock,1,'i0',memuse,stats)
      allocate( iblk1(n11), stat=stats )
      call umalloc_r(n11,1,'iblk1',memuse,stats)
      allocate( iblk2(n11), stat=stats )
      call umalloc_r(n11,1,'iblk2',memuse,stats)
      allocate( idim(nblock), stat=stats )
      call umalloc_r(nblock,1,'idim',memuse,stats)
      allocate( ihi1(n11), stat=stats )
      call umalloc_r(n11,1,'ihi1',memuse,stats)
      allocate( ihi2(n11), stat=stats )
      call umalloc_r(n11,1,'ihi2',memuse,stats)
      allocate( ilo1(n11), stat=stats )
      call umalloc_r(n11,1,'ilo1',memuse,stats)
      allocate( ilo2(n11), stat=stats )
      call umalloc_r(n11,1,'ilo2',memuse,stats)
      allocate( im(nblock), stat=stats )
      call umalloc_r(nblock,1,'im',memuse,stats)
      allocate( isva1(n11), stat=stats )
      call umalloc_r(n11,1,'isva1',memuse,stats)
      allocate( isva2(n11), stat=stats )
      call umalloc_r(n11,1,'isva2',memuse,stats)
      allocate( isvb1(n11), stat=stats )
      call umalloc_r(n11,1,'isvb1',memuse,stats)
      allocate( isvb2(n11), stat=stats )
      call umalloc_r(n11,1,'isvb2',memuse,stats)
      allocate( iv1(nblock), stat=stats )
      call umalloc_r(nblock,1,'iv1',memuse,stats)
      allocate( iv2(nblock), stat=stats )
      call umalloc_r(nblock,1,'iv2',memuse,stats)
      allocate( iv3(nblock), stat=stats )
      call umalloc_r(nblock,1,'iv3',memuse,stats)
      allocate( iz(nblock), stat=stats )
      call umalloc_r(nblock,1,'iz',memuse,stats)
      allocate( j0(nblock), stat=stats )
      call umalloc_r(nblock,1,'j0',memuse,stats)
      allocate( jdim(nblock), stat=stats )
      call umalloc_r(nblock,1,'jdim',memuse,stats)
      allocate( jhi1(n11), stat=stats )
      call umalloc_r(n11,1,'jhi1',memuse,stats)
      allocate( jhi2(n11), stat=stats )
      call umalloc_r(n11,1,'jhi2',memuse,stats)
      allocate( jlo1(n11), stat=stats )
      call umalloc_r(n11,1,'jlo1',memuse,stats)
      allocate( jlo2(n11), stat=stats )
      call umalloc_r(n11,1,'jlo2',memuse,stats)
      allocate( jm(nblock), stat=stats )
      call umalloc_r(nblock,1,'jm',memuse,stats)
      allocate( k0(nblock), stat=stats )
      call umalloc_r(nblock,1,'k0',memuse,stats)
      allocate( kdim(nblock), stat=stats )
      call umalloc_r(nblock,1,'kdim',memuse,stats)
      allocate( khi1(n11), stat=stats )
      call umalloc_r(n11,1,'khi1',memuse,stats)
      allocate( khi2(n11), stat=stats )
      call umalloc_r(n11,1,'khi2',memuse,stats)
      allocate( klo1(n11), stat=stats )
      call umalloc_r(n11,1,'klo1',memuse,stats)
      allocate( klo2(n11), stat=stats )
      call umalloc_r(n11,1,'klo2',memuse,stats)
      allocate( km(nblock), stat=stats )
      call umalloc_r(nblock,1,'km',memuse,stats)
c
      iflag=0
      write(6,'('' read a 0=formatted or 1=unformatted grid?'')')
      read(5,*) iform
      if (iform .eq. 1) then
      write(6,'('' input unformatted grid filename to read:'')')
      read(5,'(a80)') file1
      open(2,file=file1,form='unformatted',status='old')
      write(6,'('' type 0 = PLOT3D-type, 1 = CFL3D-type grid to'',
     +  '' read:'')')
      read(5,*) itype
      if (itype .eq. 0) then
c  PLOT3D-type
        read(2) nbl
        if (nbl .gt. nblock) then
          write(6,'('' need to increase nblock to '',i6)') nbl
          stop
        end if
        read(2) (idim(n),jdim(n),kdim(n),n=1,nbl)
c       do n=1,nbl
c         id=idim(n)
c         jd=jdim(n)
c         kd=kdim(n)
c         write(6,'('' block#'',i5,'': id,jd,kd='',3i5)') n,id,jd,kd
c         if (id.gt.idm .or. jd.gt.jdm .or. kd.gt.kdm) then
c           write(6,'('' need to increase idm,jdm,kdm to:'',
c    +        3i5)') id,jd,kd
c           stop
c         end if
c         read(2) (((x(i,j,k,n),i=1,id),j=1,jd),k=1,kd),
c    +            (((y(i,j,k,n),i=1,id),j=1,jd),k=1,kd),
c    +            (((z(i,j,k,n),i=1,id),j=1,jd),k=1,kd)
c       enddo
      else
c   CFL3D-type
        nbl = nblock
        do n=1,nblock
          read(2) jdim(n),kdim(n),idim(n)
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          write(6,'('' block#'',i5,'': id,jd,kd='',3i5)') n,id,jd,kd
          if (id.gt.idm .or. jd.gt.jdm .or. kd.gt.kdm) then
            write(6,'('' need to increase idm,jdm,kdm to:'',
     +        3i5)') id,jd,kd
            stop
          end if
c         read(2) (((x(i,j,k,n),j=1,jd),k=1,kd),i=1,id),
c    +            (((y(i,j,k,n),j=1,jd),k=1,kd),i=1,id),
c    +            (((z(i,j,k,n),j=1,jd),k=1,kd),i=1,id)
          read(2) (((dum,j=1,jd),k=1,kd),i=1,id),
     +            (((dum,j=1,jd),k=1,kd),i=1,id),
     +            (((dum,j=1,jd),k=1,kd),i=1,id)
        enddo
      end if
      else
      write(6,'('' input formatted grid filename to read:'')')
      read(5,'(a80)') file1
      open(2,file=file1,form='formatted',status='old')
      write(6,'('' type 0 = PLOT3D-type, 1 = CFL3D-type grid to'',
     +  '' read:'')')
      read(5,*) itype
      if (itype .eq. 0) then
c  PLOT3D-type
        read(2,*) nbl
        if (nbl .gt. nblock) then
          write(6,'('' need to increase nblock to '',i6)') nbl
          stop
        end if
        read(2,*) (idim(n),jdim(n),kdim(n),n=1,nbl)
c       do n=1,nbl
c         id=idim(n)
c         jd=jdim(n)
c         kd=kdim(n)
c         write(6,'('' block#'',i5,'': id,jd,kd='',3i5)') n,id,jd,kd
c         if (id.gt.idm .or. jd.gt.jdm .or. kd.gt.kdm) then
c           write(6,'('' need to increase idm,jdm,kdm to:'',
c    +        3i5)') id,jd,kd
c           stop
c         end if
c         read(2,*) (((x(i,j,k,n),i=1,id),j=1,jd),k=1,kd),
c    +              (((y(i,j,k,n),i=1,id),j=1,jd),k=1,kd),
c    +              (((z(i,j,k,n),i=1,id),j=1,jd),k=1,kd)
c       enddo
      else
c   CFL3D-type
        nbl = nblock
        do n=1,nblock
          read(2,*) jdim(n),kdim(n),idim(n)
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          write(6,'('' block#'',i5,'': id,jd,kd='',3i5)') n,id,jd,kd
          if (id.gt.idm .or. jd.gt.jdm .or. kd.gt.kdm) then
            write(6,'('' need to increase idm,jdm,kdm to:'',
     +        3i5)') id,jd,kd
            stop
          end if
c         read(2,*) (((x(i,j,k,n),j=1,jd),k=1,kd),i=1,id),
c    +              (((y(i,j,k,n),j=1,jd),k=1,kd),i=1,id),
c    +              (((z(i,j,k,n),j=1,jd),k=1,kd),i=1,id)
          read(2,*) (((dum,j=1,jd),k=1,kd),i=1,id),
     +              (((dum,j=1,jd),k=1,kd),i=1,id),
     +              (((dum,j=1,jd),k=1,kd),i=1,id)
        enddo
      end if
      end if
c
      icount_tot_pts=0
      icount_tot_cel=0
      do n=1,nbl
        icount_tot_pts=icount_tot_pts+(idim(n)*jdim(n)*kdim(n))
        icount_tot_cel=icount_tot_cel+((idim(n)-1)*(jdim(n)-1)*
     +   (kdim(n)-1))
      enddo
      write(6,'('' total points='',i12)') icount_tot_pts
      write(6,'('' total cells ='',i12)') icount_tot_cel
      rewind(2)
c
c   open ADF file
      write(6,'(/,'' This program creates a grid file SEPARATE'',
     + '' from the file solutiononly.cgns'')')
      write(6,'(''   (they are LINKED, eliminating redundancy'')')
      write(6,'(''    in having to keep multiple copies of the'')')
      write(6,'(''    grid for multiple restarts)'')')
      write(6,'(/,'' Input desired name for CGNS gridfile (e.g.,'',
     + '' gridname.cgns):'')')
      read(5,'(a32)') filename
      filenamr='solutiononly.cgns'
      write(6,'('' creating adf grid file name:  '',a32)') filename
      call cg_open_f(filename,CG_MODE_WRITE,iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_open_f(filenamr,CG_MODE_WRITE,iccgr,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   create base
      basename='Base'
      index_dim=3
      call cg_base_write_f(iccg,basename,index_dim,index_dim,ibase,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_base_write_f(iccgr,basename,index_dim,index_dim,ibase,ier)
      if (ier .ne. 0) call cg_error_exit_f
c   create zones (name them Zone1, Zone2, ..., ZoneN)
      do n=1,nbl
        if(n .lt. 10000) then
          write(filenum,102) n
        else
          write(6,'('' too many zones - limited to 9999'')')
          write(6,'('' Aborting.  The ADF file is no good.'')')
          stop
        end if
 102    format(i4)
        zonename(n) = 'Zone' // filenum
        isize(1)=idim(n)
        isize(2)=jdim(n)
        isize(3)=kdim(n)
        isize(4)=idim(n)-1
        isize(5)=jdim(n)-1
        isize(6)=kdim(n)-1
        isize(7)=0
        isize(8)=0
        isize(9)=0
        call cg_zone_write_f(iccg,ibase,zonename(n),isize,Structured,
     .    iz(n),ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_zone_write_f(iccgr,ibase,zonename(n),isize,Structured,
     .    iz(n),ier)
        if (ier .ne. 0) call cg_error_exit_f
      enddo
c  create and write the grid coordinates only to grid file
c  will put a link in the soln.cgns file (below)
      if (iform .eq. 1 .and. itype .eq. 0) then
        read(2) idum
        read(2) (idum,idum,idum,n=1,nbl)
        do n=1,nbl
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          write(6,'('' block#'',i5,'': id,jd,kd='',3i5)') n,id,jd,kd
          allocate( x(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'x',memuse,stats)
          allocate( y(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'y',memuse,stats)
          allocate( z(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'z',memuse,stats)
          read(2) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd),
     +            (((y(i,j,k),i=1,id),j=1,jd),k=1,kd),
     +            (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)
          call writexyz(iccg,ibase,idouble,iz(n),id,jd,kd,
     +      id,jd,kd,x,y,z)
          deallocate(x,y,z)
        enddo
      else if (iform .eq. 1 .and. itype .eq. 1) then
        do n=1,nbl
          read(2) idum,idum,idum
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          allocate( x(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'x',memuse,stats)
          allocate( y(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'y',memuse,stats)
          allocate( z(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'z',memuse,stats)
          read(2) (((x(i,j,k),j=1,jd),k=1,kd),i=1,id),
     +            (((y(i,j,k),j=1,jd),k=1,kd),i=1,id),
     +            (((z(i,j,k),j=1,jd),k=1,kd),i=1,id)
          call writexyz(iccg,ibase,idouble,iz(n),id,jd,kd,
     +      id,jd,kd,x,y,z)
          deallocate(x,y,z)
        enddo
      else if (iform .eq. 0 .and. itype .eq. 0) then
        read(2,*) idum
        read(2,*) (idum,idum,idum,n=1,nbl)
        do n=1,nbl
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          allocate( x(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'x',memuse,stats)
          allocate( y(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'y',memuse,stats)
          allocate( z(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'z',memuse,stats)
          read(2,*) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd),
     +              (((y(i,j,k),i=1,id),j=1,jd),k=1,kd),
     +              (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)
          call writexyz(iccg,ibase,idouble,iz(n),id,jd,kd,
     +      id,jd,kd,x,y,z)
          deallocate(x,y,z)
        enddo
      else if (iform .eq. 0 .and. itype .eq. 1) then
        do n=1,nbl
          read(2,*) idum,idum,idum
          id=idim(n)
          jd=jdim(n)
          kd=kdim(n)
          allocate( x(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'x',memuse,stats)
          allocate( y(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'y',memuse,stats)
          allocate( z(id,jd,kd), stat=stats )
          call umalloc_r(id*jd*kd,0,'z',memuse,stats)
          read(2,*) (((x(i,j,k),j=1,jd),k=1,kd),i=1,id),
     +              (((y(i,j,k),j=1,jd),k=1,kd),i=1,id),
     +              (((z(i,j,k),j=1,jd),k=1,kd),i=1,id)
          call writexyz(iccg,ibase,idouble,iz(n),id,jd,kd,
     +      id,jd,kd,x,y,z)
          deallocate(x,y,z)
        enddo
      end if
c  reopen the CFL3D input file
      open(3,file=file2,form='formatted',status='old')
c------------------------------------------------------------------------
c   This section reads in original CFL3D input file
c   and also writes the BC part to the CGNS file,
c   on the fly.  The 1-to-1 info read here is saved and used later to
c   write to the CGNS file.
c------------------------------------------------------------------------
      read(3,'(a80)') name
      read(3,'(a80)') name
      do n=1,11
        read(3,'(a80)') name
      enddo
      read(3,'(a80)') name
c
      read(3,'(a1)') dum1
      if (dum1 .eq. '>') then
        do n=1,500
          read(3,'(a1)') dum1
          if (dum1 .eq. '<') goto 1002
        enddo
        write(6,'('' Error, too many lines (>500) of keyword input'')')
        stop
 1002   continue
      else
        backspace(3)
      end if
c
      read(3,'(a80)') title
      write(6,'('' Reading CFL3D input file with following title:'')')
      write(6,'(a80)') title
      read(3,'(a80)') name
      read(3,*) xm,al,be,re,t,ia,ih
      read(3,'(a80)') name
      read(3,*) sr,cr,br,xm,ym,zm
      read(3,'(a80)') name
      read(3,*) dt,ir,if,fm,iu,cf
c     if(dt .gt. 0. and. iu .ne. 0) then
c       write(6,'('' Error!  CFL3D input file is set up to do'',
c    +   '' moving grid'')')
c       write(6,'('' CGNS not capable of handling this!'')')
c       write(6,'('' Stopping'')')
c       stop
c     end if
      read(3,'(a80)') name
      read(3,*) ngrid,nplot3d,nprint,nw,ic,i2d,nt,it
      read(3,'(a80)') name
      ngr=abs(ngrid)
      do n=1,ngr
        read(3,*) nc,ie,ia,if,iv1(n),iv2(n),iv3(n)
      enddo
      i2d=abs(i2d)
      read(3,'(a80)') name
      do n=1,ngr
        read(3,*) id,jd,kd
        if(id .ne. idim(n) .or. jd .ne. jdim(n) .or.
     +     kd .ne. kdim(n)) then
          write(6,'('' Error!  Input file indices do not agree with'',
     +     '' grid indices'')')
          write(6,'('' Zone '',i5,'' id,jd,kd='',3i5)') n,id,jd,kd
          write(6,'('' Grid indices are: '',3i5)') idim(n),jdim(n),
     +     kdim(n)
          write(6,'('' Aborting.  The ADF file is no good.'')')
          stop
        end if
      enddo
      read(3,'(a80)') name
      do n=1,ngr
        read(3,*) i1,i2,i3,i4,i5,i6
      enddo
      read(3,'(a80)') name
      do n=1,ngr
        read(3,*) in,ig,i1,i2,i3,i4,i5,i6
      enddo
c  Idiag section:
      read(3,'(a80)') name
      do n=1,ngr
        read(3,*) i1,i2,i3,i4,i5,i6
      enddo
c  Ifds section:
      read(3,'(a80)') name
      do n=1,ngr
        read(3,*) if1,if2,if3,rk1,rk2,rk3
      enddo
c  Grid section (DO NOT write this out - it's going into CGNS file!):
      read(3,'(a80)') name
      iovmx=0
      do n=1,ngr
        read(3,*) ig,i0(n),im(n),j0(n),jm(n),k0(n),km(n),iov
        iovmx=max(iov,iovmx)
      enddo
c     write a descriptor if there is overlap data
      if (iovmx .ne. 0) then
        textf='Overset connectivity information is'//
     +   ' NOT included in this file'
        call cg_goto_f(iccg,ibase,ier,'end')
        call cg_descriptor_write_f('OversetInformation',textf,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_goto_f(iccgr,ibase,ier,'end')
        call cg_descriptor_write_f('OversetInformation',textf,ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
c  Grid i0 section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,i0(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=jdim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=kdim(n)
          ipnts(1,1)=1
          ipnts(2,1)=i1
          ipnts(3,1)=i3
          ipnts(1,2)=1
          ipnts(2,2)=i2
          ipnts(3,2)=i4
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Ilo_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Ilo_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv1(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv1(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' i0:   1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv1(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Grid idim section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,im(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=jdim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=kdim(n)
          ipnts(1,1)=idim(n)
          ipnts(2,1)=i1
          ipnts(3,1)=i3
          ipnts(1,2)=idim(n)
          ipnts(2,2)=i2
          ipnts(3,2)=i4
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Ihi_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Ihi_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv1(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv1(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' imax: 1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv1(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Grid j0 section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,j0(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=idim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=kdim(n)
          ipnts(1,1)=i1
          ipnts(2,1)=1
          ipnts(3,1)=i3
          ipnts(1,2)=i2
          ipnts(2,2)=1
          ipnts(3,2)=i4
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Jlo_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Jlo_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv2(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv2(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' j0:   1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv2(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Grid jdim section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,jm(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=idim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=kdim(n)
          ipnts(1,1)=i1
          ipnts(2,1)=jdim(n)
          ipnts(3,1)=i3
          ipnts(1,2)=i2
          ipnts(2,2)=jdim(n)
          ipnts(3,2)=i4
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Jhi_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Jhi_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv2(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv2(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' jmax: 1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv2(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Grid k0 section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,k0(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=idim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=jdim(n)
          ipnts(1,1)=i1
          ipnts(2,1)=i3
          ipnts(3,1)=1
          ipnts(1,2)=i2
          ipnts(2,2)=i4
          ipnts(3,2)=1
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Klo_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Klo_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv3(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv3(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' k0:   1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv3(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Grid kdim section:
      read(3,'(a80)') name
      do n=1,ngr
        do m=1,km(n)
          read(3,*) ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(3,'(a80)') name
            read(3,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(3,'(a80)') name
            read(3,'(a80)') name
          end if
c XXX
          if (i1 .eq. 0) i1=1
          if (i2 .eq. 0) i2=idim(n)
          if (i3 .eq. 0) i3=1
          if (i4 .eq. 0) i4=jdim(n)
          ipnts(1,1)=i1
          ipnts(2,1)=i3
          ipnts(3,1)=kdim(n)
          ipnts(1,2)=i2
          ipnts(2,2)=i4
          ipnts(3,2)=kdim(n)
          if(m .lt. 10000) then
            write(filenum,102) m
          else
            write(6,'('' too many segments - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
          if (is .gt. 0) then
            boconame='Khi_Seg' // filenum
          else
            write(6,'('' Ignore-force-calc on specific BC segs'',
     +       '' (segment<0) not implemented'')')
            write(6,'(''   for CGNS.  Continuing, with CGNS'',
     +       '' file indicating via'')')
            write(6,'(''   a <<NoForce>> at beginning of the BC'',
     +       '' name'')')
            boconame='NoForce_Khi_Seg' // filenum
          end if
          if(     ib .eq. 1000) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCInflowSupersonic,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1001) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPlane,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCExtrapolate,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1003) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCFarfield,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1005 .or. ib .eq. 1006) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallInviscid,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1008 .or. ib .eq. 2008 .or. ib .eq. 2018
     +       .or. ib .eq. 2028 .or. ib .eq. 2009 .or. ib .eq. 2010) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelInflow,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1011) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCSymmetryPolar,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 1013) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCDegenerateLine,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 2002) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCTunnelOutflow,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .eq. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv3(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousHeatFlux,PointRange,2,ipnts,ibc,ier)
          else if((ib .eq. 2004 .or. ib .eq. 2014) .and. 
     +            data(1) .ne. 0 .and.
     +            data(2) .eq. 0 .and. abs(iv3(n)) .gt. 0) then
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        BCWallViscousIsothermal,PointRange,2,ipnts,ibc,ier)
          else if(ib .eq. 0) then
            write(6,'('' kmax: 1-1, patch, or overset BC at '',
     +       ''zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)') n,
     +       ipnts(1,1),ipnts(1,2),ipnts(2,1),ipnts(2,2),
     +       ipnts(3,1),ipnts(3,2)
          else
            call cg_boco_write_f(iccg,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            call cg_boco_write_f(iccgr,ibase,iz(n),boconame,
     +        UserDefined,PointRange,2,ipnts,ibc,ier)
            iflag=1
c           write(6,'('' BCtype '',i6,'' not implemented for'',
c    +       '' CGNS yet'')') ib
          end if
          if (ier .ne. 0) call cg_error_exit_f
          if (ib .ne. 0) then
            boconame2='CFL3DType'
            boconame3='CFL3DAuxData'
            call cg_goto_f(iccg,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            write(texta,'(i6)') ib
            text=texta
            if((abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014) .and. 
     +          iv3(n) .lt. 0) then
              text=texta // ' + WallFunction'
              call cg_bc_wallfunction_write_f(iccg,ibase,iz(n),
     +          ibc,Generic,ier)
              call cg_bc_wallfunction_write_f(iccgr,ibase,iz(n),
     +          ibc,Generic,ier)
              if (ier .ne. 0) call cg_error_exit_f
            end if
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000 .and. nd .gt. 0) then
              backspace 3
              read(3,'(a80)') name
            end if
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',iz(n),
     +        'ZoneBC_t',1,'BC_t',ibc,'end')
            call cg_descriptor_write_f(boconame2,text,ier)
            if (abs(ib) .ge. 2000) then
              call cg_descriptor_write_f(boconame3,name,ier)
            end if
            if (ier .ne. 0) call cg_error_exit_f
          end if
c   XXX
        enddo
      enddo
c  Mseq section
      read(3,'(a80)') name
      read(3,*) mseq,mg,ic,mtt,ng
c  Issc section
      read(3,'(a80)') name
      read(3,*) i1,f2,f3,f4,i5,f6,f7,f8
c  Ncyc section
      read(3,'(a80)') name
      do n=1,mseq
        read(3,*) nc,mg,ne,ni
      enddo
c  Mit section
      read(3,'(a80)') name
      do n=1,mseq
        read(3,'(a80)') name
      enddo
c  1-to-1 data (DO NOT write this out - it's going into CGNS file!):
      read(3,'(a80)') name
      read(3,'(a80)') name
      read(3,*) num1_1
c  Limits1:
      read(3,'(a80)') name
      do n=1,num1_1
        read(3,*) nn,iblk1(n),ilo1(n),jlo1(n),klo1(n),ihi1(n),
     +              jhi1(n),khi1(n),isva1(n),isvb1(n)
      enddo
c  Limits2:
      read(3,'(a80)') name
      do n=1,num1_1
        read(3,*) nn,iblk2(n),ilo2(n),jlo2(n),klo2(n),ihi2(n),
     +              jhi2(n),khi2(n),isva2(n),isvb2(n)
      enddo
c  Patch:
      read(3,'(a80)') name
      read(3,'(a80)') name
      read(3,*) ninter
c     write a descriptor if there is patch data
      if (ninter .lt. 0) then
        textf='Patching connectivity information is'//
     +   ' NOT included in this file'
        call cg_goto_f(iccg,ibase,ier,'end')
        call cg_descriptor_write_f('PatchInformation',textf,ier)
        if (ier .ne. 0) call cg_error_exit_f
        call cg_goto_f(iccgr,ibase,ier,'end')
        call cg_descriptor_write_f('PatchInformation',textf,ier)
        if (ier .ne. 0) call cg_error_exit_f
      end if
c  Plot3d:
      read(3,'(a80)') name
      read(3,'(a80)') name
      if (abs(nplot3d) .gt. 0) then
      do n=1,abs(nplot3d)
        read(3,*) ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
      enddo
      end if
c  Movie:
      read(3,'(a80)') name
      read(3,*) mo
c  Print-out:
      read(3,'(a80)') name
      read(3,'(a80)') name
      if (abs(nprint) .gt. 0) then
      do n=1,abs(nprint)
        read(3,*) ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
      enddo
      end if
c  Control surface:
      read(3,'(a80)') name
      read(3,'(a80)') name
      read(3,*) ncs
      read(3,'(a80)') name
      if (ncs .gt. 0) then
      do n=1,ncs
        read(3,*) ig,i1,i2,i3,i4,i5,i6,iw,inorm
      enddo
      end if
c  Mov grid stuff - trans:
      read(3,'(a80)',end=999) name
      write(6,'('' This cfl3d input file has moving grid'',
     + '' stuff at the bottom'')')
      write(6,'('' CGNS not capable of handling this!'')')
      write(6,'('' Not aborting, but none of the moving grid'',
     + ''stuff is being kept.'')')
 999  continue
c------------------------------------------------------------------------
        write(6,'('' there are '',i5,'' 1-to-1 interfaces being'',
     +   '' read'')') num1_1
        if (num1_1 .gt. n11) then
          write(6,'('' need to increase n11 to '',i6)') num1_1
          write(6,'('' Aborting.  The ADF file is no good.'')')
          stop
        end if
c  now work with the 1-to-1 info (recall there will be 2 nodes written
c  for each interface)
        do n=1,num1_1
c  get interface name
          if(n .lt. 10000) then
            write(filecon,102) n
          else
            write(6,'('' too many zones - limited to 9999'')')
            write(6,'('' Aborting.  The ADF file is no good.'')')
            stop
          end if
c  this section repeats for each side:
c  get connect name
          connectname = '1to1InterfaceA' // filecon
c  get this zone number
          izone=iz(iblk1(n))
c  get donor zone name
          donorname=zonename(iblk2(n))
c  get point ranges (this zone)
          irange(1,1)=ilo1(n)
          irange(2,1)=jlo1(n)
          irange(3,1)=klo1(n)
          irange(1,2)=ihi1(n)
          irange(2,2)=jhi1(n)
          irange(3,2)=khi1(n)
c         write(6,'('' irange='',6i4)') irange(1,1),irange(2,1),
c    +      irange(3,1),irange(1,2),irange(2,2),irange(3,2)
c  get point ranges (donor zone)
          iranged(1,1)=ilo2(n)
          iranged(2,1)=jlo2(n)
          iranged(3,1)=klo2(n)
          iranged(1,2)=ihi2(n)
          iranged(2,2)=jhi2(n)
          iranged(3,2)=khi2(n)
c         write(6,'('' iranged='',6i4)') iranged(1,1),iranged(2,1),
c    +      iranged(3,1),iranged(1,2),iranged(2,2),iranged(3,2)
          write(6,'(''    Pair:'')')
          write(6,'(''    zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)')
     +     iblk1(n),ilo1(n),ihi1(n),jlo1(n),jhi1(n),klo1(n),khi1(n) 
          write(6,'(''    zone,  ilo,ihi,jlo,jhi,klo,khi='',i4,4x,6i4)')
     +     iblk2(n),ilo2(n),ihi2(n),jlo2(n),jhi2(n),klo2(n),khi2(n) 
c   call routine to get transform matrix
          call transget(ilo1(n),jlo1(n),klo1(n),ihi1(n),jhi1(n),
     +      khi1(n),isva1(n),isvb1(n),ilo2(n),jlo2(n),klo2(n),ihi2(n),
     +      jhi2(n),khi2(n),isva2(n),isvb2(n),itransform)
c         write(6,'('' itransform='',3i5)') itransform(1),
c    +      itransform(2),itransform(3)
c   write info
          call cg_1to1_write_f(iccg,ibase,izone,connectname,donorname,
     +       irange,iranged,itransform,iindex,ier)
          if (ier .ne. 0) call cg_error_exit_f
          call cg_1to1_write_f(iccgr,ibase,izone,connectname,donorname,
     +       irange,iranged,itransform,iindex,ier)
          if (ier .ne. 0) call cg_error_exit_f
c  this section repeats for each side:
c  get connect name
          connectname = '1to1InterfaceB' // filecon
c  get this zone number
          izone=iz(iblk2(n))
c  get donor zone name
          donorname=zonename(iblk1(n))
c  get point ranges (this zone)
          irange(1,1)=ilo2(n)
          irange(2,1)=jlo2(n)
          irange(3,1)=klo2(n)
          irange(1,2)=ihi2(n)
          irange(2,2)=jhi2(n)
          irange(3,2)=khi2(n)
c         write(6,'('' irange='',6i4)') irange(1,1),irange(2,1),
c    +      irange(3,1),irange(1,2),irange(2,2),irange(3,2)
c  get point ranges (donor zone)
          iranged(1,1)=ilo1(n)
          iranged(2,1)=jlo1(n)
          iranged(3,1)=klo1(n)
          iranged(1,2)=ihi1(n)
          iranged(2,2)=jhi1(n)
          iranged(3,2)=khi1(n)
c         write(6,'('' iranged='',6i4)') iranged(1,1),iranged(2,1),
c    +      iranged(3,1),iranged(1,2),iranged(2,2),iranged(3,2)
c  call routine to get transform matrix
          call transget(ilo2(n),jlo2(n),klo2(n),ihi2(n),jhi2(n),
     +      khi2(n),isva2(n),isvb2(n),ilo1(n),jlo1(n),klo1(n),ihi1(n),
     +      jhi1(n),khi1(n),isva1(n),isvb1(n),itransform)
c         write(6,'('' itransform='',3i5)') itransform(1),
c    +      itransform(2),itransform(3)
c  write info
          call cg_1to1_write_f(iccg,ibase,izone,connectname,donorname,
     +       irange,iranged,itransform,iindex,ier)
          if (ier .ne. 0) call cg_error_exit_f
          call cg_1to1_write_f(iccgr,ibase,izone,connectname,donorname,
     +       irange,iranged,itransform,iindex,ier)
          if (ier .ne. 0) call cg_error_exit_f
        enddo
c
c  if there are patched or overset connectivities:
c
      if (ninter .lt. 0 .or. iovmx .ne. 0) then
        write(6,'(/,'' There appear to be patched and/or overset'',
     +   '' boundaries'')')
        write(6,'('' There are currently Information Descriptor'',
     +   '' nodes written to the CGNS file'')')
        write(6,'('' If you also want GridConnectivity_t node(s)'',
     +   '' written, you must supply info here:'')')
c  Patched:
        write(6,'(''   Input number of patched faces'',
     +   '' (typically 2 for each interface)'')')
        write(6,'(''   (type 0 if there are none, or if you'',
     +   '' want to leave out this info):'')')
        read(5,*) ipatch
        if (ipatch .gt. 0) then
          do i=1,ipatch
            write(6,'('' Input to-zone index information:'')')
            write(6,'('' Input zone #,ilo,ihi,jlo,jhi,klo,khi:'')')
            read(5,*) izto,iloto,ihito,jloto,jhito,kloto,khito
            irange(1,1)=iloto
            irange(2,1)=jloto
            irange(3,1)=kloto
            irange(1,2)=ihito
            irange(2,2)=jhito
            irange(3,2)=khito
            iok=0
            if (iloto .eq. ihito .or. jloto .eq. jhito .or.
     +          kloto .eq. khito) then
              iok=1
            end if
            if (iok .eq. 0) then
              write(6,'('' Error... one index must be identical.'',
     +         ''  Stopping.'')')
              stop
            end if
c  get donor zone name
            write(6,'('' Input primary number of zone donor for this'',
     +       '' patch:'')')
            read(5,*) izfr
            donorname=zonename(izfr)
c  get interface name
            if(i .lt. 10000) then
              write(filecon,102) i
            else
              write(6,'('' too many zones - limited to 9999'')')
              write(6,'('' Aborting.  The ADF file is no good.'')')
              stop
            end if
c  get connect name
            connectname = 'PatchInterfaceA' // filecon
            call cg_conn_write_short_f(iccg,ibase,iz(izto),connectname,
     +        Vertex,Abutting,PointRange,2,irange,donorname,
     +        iindexto,ier)
            if (ier .ne. 0) call cg_error_exit_f
            call cg_conn_write_short_f(iccgr,ibase,iz(izto),connectname,
     +        Vertex,Abutting,PointRange,2,irange,donorname,
     +        iindexto,ier)
            if (ier .ne. 0) call cg_error_exit_f
c  write descriptors to go along with patching
            call cg_goto_f(iccg,ibase,ier,'Zone_t',izto,
     +       'ZoneGridConnectivity_t',1,'GridConnectivity_t',
     +       iindexto,'end')
            textf='Patching connectivity InterpolantsDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('InterpolantsDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
            textf='Patching connectivity CellListDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('CellListDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
c
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',izto,
     +       'ZoneGridConnectivity_t',1,'GridConnectivity_t',
     +       iindexto,'end')
            textf='Patching connectivity InterpolantsDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('InterpolantsDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
            textf='Patching connectivity CellListDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('CellListDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
          enddo
        end if
c  Overset:
        write(6,'(''   Input number of overset faces'',
     +   '' (typically 2 for each interface)'')')
        write(6,'(''   (type 0 if there are none, or if you'',
     +   '' want to leave out this info):'')')
        read(5,*) iover
        if (iover .gt. 0) then
          do i=1,iover
            write(6,'('' Input to-zone index information:'')')
            write(6,'('' Input zone #,ilo,ihi,jlo,jhi,klo,khi:'')')
            read(5,*) izto,iloto,ihito,jloto,jhito,kloto,khito
            irange(1,1)=iloto
            irange(2,1)=jloto
            irange(3,1)=kloto
            irange(1,2)=ihito
            irange(2,2)=jhito
            irange(3,2)=khito
            iok=0
            if (iloto .eq. ihito .or. jloto .eq. jhito .or.
     +          kloto .eq. khito) then
              iok=1
            end if
            if (iok .eq. 0) then
              write(6,'('' Error... one index must be identical.'',
     +         ''  Stopping.'')')
              stop
            end if
c  get donor zone name
            write(6,'('' Input primary number of zone donor for this'',
     +       '' patch:'')')
            read(5,*) izfr
            donorname=zonename(izfr)
c  get interface name
            if(i .lt. 10000) then
              write(filecon,102) i
            else
              write(6,'('' too many zones - limited to 9999'')')
              write(6,'('' Aborting.  The ADF file is no good.'')')
              stop
            end if
c  get connect name
            connectname = 'OversetInterfaceA' // filecon
            call cg_conn_write_short_f(iccg,ibase,iz(izto),connectname,
     +        Vertex,Overset,PointRange,2,irange,donorname,
     +        iindexto,ier)
            if (ier .ne. 0) call cg_error_exit_f
            call cg_conn_write_short_f(iccgr,ibase,iz(izto),connectname,
     +        Vertex,Overset,PointRange,2,irange,donorname,
     +        iindexto,ier)
            if (ier .ne. 0) call cg_error_exit_f
c  write descriptors to go along with overset
            call cg_goto_f(iccg,ibase,ier,'Zone_t',izto,
     +       'ZoneGridConnectivity_t',1,'GridConnectivity_t',
     +       iindexto,'end')
            textf='Overset connectivity InterpolantsDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('InterpolantsDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
            textf='Overset connectivity CellListDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('CellListDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
c
            call cg_goto_f(iccgr,ibase,ier,'Zone_t',izto,
     +       'ZoneGridConnectivity_t',1,'GridConnectivity_t',
     +       iindexto,'end')
            textf='Overset connectivity InterpolantsDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('InterpolantsDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
            textf='Overset connectivity CellListDonor info is'//
     +       ' NOT included in this file'
            call cg_descriptor_write_f('CellListDonorInfo',
     +       textf,ier)
            if (ier .ne. 0) call cg_error_exit_f
          enddo
        end if
c
      end if
c  close adf file
      call cg_close_f(iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_close_f(iccgr,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      call cg_open_f(filenamr,CG_MODE_MODIFY,iccgr,ier)
      if (ier .ne. 0) call cg_error_exit_f
      do n=1,nbl
        newname='Base/' // zonename(n) // '/GridCoordinates'
c  old way:
c       call cg_zone_id_f(iccgr,ibase,n,pidz,ier)
c       if (ier .ne. 0) call cg_error_exit_f
c       call ADFLINK(pidz,'GridCoordinates',filename,newname,
c    +    pidnew,ier)
c  new way:
        call cg_goto_f(iccgr,ibase,ier,'Zone_t',n,'end')
        if (ier .ne. 0) call cg_error_exit_f
        call cg_link_write_f('GridCoordinates',filename,newname,ier)
        if (ier .ne. 0) call cg_error_exit_f
c   *********
      enddo
      call cg_close_f(iccgr,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (iflag .eq. 1) then
        write(6,'('' One or more BCs not officially supported'',
     +   '' or not fully implemented yet'')')
        write(6,'(''   ...these have been set to <UserDefined>'',
     +   '' at present'')')
      end if
c
      write(6,'(/,'' Successful completion'',/)')
      write(6,'('' CGNS grid file written to:  '',a32)') filename
      write(6,'('' Separate CGNS file with LINK to grid file '',
     +  ''written to:  '',a32)') filenamr
c
      write(6,'(/,'' When you run CFL3D, you have the following'',
     + '' choice:'')')
      write(6,'(''    1. Use LINK (keeping grid file separate from'',
     + '' solution file)'')')
      write(6,'(''       Do this if your grid file is BIG, and you'',
     + '' do not want'')')
      write(6,'(''       it re-written into every restart file.'',
     + ''  For this option,'')')
      write(6,'(''       KEEP the following file available (in'',
     + '' the same directory):'',/,15x,a32)') filename
      write(6,'(''       (this file name MUST NOT BE CHANGED)'')')
      write(6,'(''       and put the following filename on the 2nd'',
     + '' line of your input file:'',/,15x,a32)') filenamr
      write(6,'(''       (this filename can be changed, if desired)'')')
      write(6,'(''    2. Keep grid and solution in the same CGNS'',
     + '' file.  For this option,'')')
      write(6,'(''       DISCARD the following file:'',/,
     +   15x,a32)') filenamr
      write(6,'(''       and put the following filename '',
     + ''on the 2nd line of your input file:'',/,15x,a32)') filename
      write(6,'(''       (this filename can be changed, if desired)'')')
c
      write(6,'(/,'' Remember: the <restart.bin> file is no longer'',
     + '' used (all info is now written'')')
      write(6,'('' to the CGNS file named on the 2nd line of your'',
     + '' input file!)'')')
      write(6,'('' However, your input file should still have a '',
     + ''(dummy) restart name in its place'')')
c
      write(6,'(/,'' Note: although the BCs and 1-to-1 connectivity'',
     + '' have been written to'')')
      write(6,'('' the CGNS file, CFL3D currently still reads this'',
     + '' info from the input file!'')')
      write(6,'('' Also, currently no patched or overset info is put'',
     + '' into the file'')')
c
      write(6,'(/,'' Be sure to set keyword icgns=1 in the input'',
     + '' file'')')
c
c     deallocate memory
c
      deallocate(idim)
      deallocate(jdim)
      deallocate(kdim)
      deallocate(iz)
      deallocate(iv1)
      deallocate(iv2)
      deallocate(iv3)
      deallocate(iblk1)
      deallocate(ilo1)
      deallocate(jlo1)
      deallocate(klo1)
      deallocate(ihi1)
      deallocate(jhi1)
      deallocate(khi1)
      deallocate(isva1)
      deallocate(isvb1)
      deallocate(iblk2)
      deallocate(ilo2)
      deallocate(jlo2)
      deallocate(klo2)
      deallocate(ihi2)
      deallocate(jhi2)
      deallocate(khi2)
      deallocate(isva2)
      deallocate(isvb2)
      deallocate(i0)
      deallocate(im)
      deallocate(j0)
      deallocate(jm)
      deallocate(k0)
      deallocate(km)
c
      return
      end
c
c **********************************************************************
c
      subroutine writexyz(iccg,ib,idouble,iz,idim,jdim,kdim,
     +    idm,jdm,kdm,x,y,z)
c
c   Must use a routine like this, because cg_coord_write_f writes the
c   ENTIRE array.  In the main routine, x,y, and z are dimensioned LARGER
c   than necessary.  Here, a dummy array (wk) is used, which is dimensioned
c   EXACTLY correctly (idim,jdim,kdim).
c
#     include "cgnslib_f.h"
c
      dimension x(idm,jdm,kdm),y(idm,jdm,kdm),z(idm,jdm,kdm)
      dimension wk(idim,jdim,kdim)
c
      do i=1,idim
      do j=1,jdim
      do k=1,kdim
        wk(i,j,k)=x(i,j,k)
      enddo
      enddo
      enddo
      if (idouble .eq. 1) then
        call cg_coord_write_f(iccg,ib,iz,RealDouble,
     +    'CoordinateX',wk,ic,ier)
      else
        call cg_coord_write_f(iccg,ib,iz,RealSingle,
     +    'CoordinateX',wk,ic,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      do i=1,idim
      do j=1,jdim
      do k=1,kdim
        wk(i,j,k)=y(i,j,k)
      enddo
      enddo
      enddo
      if (idouble .eq. 1) then
        call cg_coord_write_f(iccg,ib,iz,RealDouble,
     +    'CoordinateY',wk,ic,ier)
      else
        call cg_coord_write_f(iccg,ib,iz,RealSingle,
     +    'CoordinateY',wk,ic,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      do i=1,idim
      do j=1,jdim
      do k=1,kdim
        wk(i,j,k)=z(i,j,k)
      enddo
      enddo
      enddo
      if (idouble .eq. 1) then
        call cg_coord_write_f(iccg,ib,iz,RealDouble,
     +    'CoordinateZ',wk,ic,ier)
      else
        call cg_coord_write_f(iccg,ib,iz,RealSingle,
     +    'CoordinateZ',wk,ic,ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
c
c **********************************************************************
c
      subroutine transget(ilo1,jlo1,klo1,ihi1,jhi1,
     +      khi1,isva1,isvb1,ilo2,jlo2,klo2,ihi2,
     +      jhi2,khi2,isva2,isvb2,itransform)
c
      dimension itransform(3)
c
          itransform(1)=0
          itransform(2)=0
          itransform(3)=0
c   do 1st varying index
          if      (isva1 .eq. 1 .and. isva2 .eq. 1) then
            itransform(1)=1
            if((ihi1-ilo1)*(ihi2-ilo2) .lt. 0) 
     +        itransform(1)=-itransform(1)
          else if (isva1 .eq. 1 .and. isva2 .eq. 2) then
            itransform(1)=2
            if((ihi1-ilo1)*(jhi2-jlo2) .lt. 0) 
     +        itransform(1)=-itransform(1)
          else if (isva1 .eq. 1 .and. isva2 .eq. 3) then
            itransform(1)=3
            if((ihi1-ilo1)*(khi2-klo2) .lt. 0)
     +        itransform(1)=-itransform(1)
c
          else if (isva1 .eq. 2 .and. isva2 .eq. 1) then
            itransform(2)=1
            if((jhi1-jlo1)*(ihi2-ilo2) .lt. 0)
     +        itransform(2)=-itransform(2)
          else if (isva1 .eq. 2 .and. isva2 .eq. 2) then
            itransform(2)=2
            if((jhi1-jlo1)*(jhi2-jlo2) .lt. 0)
     +        itransform(2)=-itransform(2)
          else if (isva1 .eq. 2 .and. isva2 .eq. 3) then
            itransform(2)=3
            if((jhi1-jlo1)*(khi2-klo2) .lt. 0)
     +        itransform(2)=-itransform(2)
c
          else if (isva1 .eq. 3 .and. isva2 .eq. 1) then
            itransform(3)=1
            if((khi1-klo1)*(ihi2-ilo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          else if (isva1 .eq. 3 .and. isva2 .eq. 2) then
            itransform(3)=2
            if((khi1-klo1)*(jhi2-jlo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          else if (isva1 .eq. 3 .and. isva2 .eq. 3) then
            itransform(3)=3
            if((khi1-klo1)*(khi2-klo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          end if
c   do 2nd varying index
          if      (isvb1 .eq. 1 .and. isvb2 .eq. 1) then
            itransform(1)=1
            if((ihi1-ilo1)*(ihi2-ilo2) .lt. 0)
     +        itransform(1)=-itransform(1)
          else if (isvb1 .eq. 1 .and. isvb2 .eq. 2) then
            itransform(1)=2
            if((ihi1-ilo1)*(jhi2-jlo2) .lt. 0)
     +        itransform(1)=-itransform(1)
          else if (isvb1 .eq. 1 .and. isvb2 .eq. 3) then
            itransform(1)=3
            if((ihi1-ilo1)*(khi2-klo2) .lt. 0)
     +        itransform(1)=-itransform(1)
c
          else if (isvb1 .eq. 2 .and. isvb2 .eq. 1) then
            itransform(2)=1
            if((jhi1-jlo1)*(ihi2-ilo2) .lt. 0)
     +        itransform(2)=-itransform(2)
          else if (isvb1 .eq. 2 .and. isvb2 .eq. 2) then
            itransform(2)=2
            if((jhi1-jlo1)*(jhi2-jlo2) .lt. 0)
     +        itransform(2)=-itransform(2)
          else if (isvb1 .eq. 2 .and. isvb2 .eq. 3) then
            itransform(2)=3
            if((jhi1-jlo1)*(khi2-klo2) .lt. 0)
     +        itransform(2)=-itransform(2)
c
          else if (isvb1 .eq. 3 .and. isvb2 .eq. 1) then
            itransform(3)=1
            if((khi1-klo1)*(ihi2-ilo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          else if (isvb1 .eq. 3 .and. isvb2 .eq. 2) then
            itransform(3)=2
            if((khi1-klo1)*(jhi2-jlo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          else if (isvb1 .eq. 3 .and. isvb2 .eq. 3) then
            itransform(3)=3
            if((khi1-klo1)*(khi2-klo2) .lt. 0)
     +        itransform(3)=-itransform(3)
          end if
c
c   Need to fill the "unused" index with the direction not currently employed
          irev=1
          if ((isva1 .eq. 1 .and. isvb1 .eq. 2) .or.
     +        (isva1 .eq. 2 .and. isvb1 .eq. 1)) then
            if ((isva2 .eq. 1 .and. isvb2 .eq. 2) .or.
     +          (isva2 .eq. 2 .and. isvb2 .eq. 1)) then
              if ((klo1 .eq. 1 .and. klo2 .eq. 1) .or.
     +            (klo1 .ne. 1 .and. klo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 1 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 1)) then
              if ((klo1 .eq. 1 .and. jlo2 .eq. 1) .or.
     +            (klo1 .ne. 1 .and. jlo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 2 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 2)) then
              if ((klo1 .eq. 1 .and. ilo2 .eq. 1) .or.
     +            (klo1 .ne. 1 .and. ilo2 .ne. 1)) then
                irev=-1
              end if
            end if
          else if ((isva1 .eq. 1 .and. isvb1 .eq. 3) .or.
     +             (isva1 .eq. 3 .and. isvb1 .eq. 1)) then
            if ((isva2 .eq. 1 .and. isvb2 .eq. 2) .or.
     +          (isva2 .eq. 2 .and. isvb2 .eq. 1)) then
              if ((jlo1 .eq. 1 .and. klo2 .eq. 1) .or.
     +            (jlo1 .ne. 1 .and. klo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 1 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 1)) then
              if ((jlo1 .eq. 1 .and. jlo2 .eq. 1) .or.
     +            (jlo1 .ne. 1 .and. jlo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 2 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 2)) then
              if ((jlo1 .eq. 1 .and. ilo2 .eq. 1) .or.
     +            (jlo1 .ne. 1 .and. ilo2 .ne. 1)) then
                irev=-1
              end if
            end if
          else if ((isva1 .eq. 2 .and. isvb1 .eq. 3) .or.
     +             (isva1 .eq. 3 .and. isvb1 .eq. 2)) then
            if ((isva2 .eq. 1 .and. isvb2 .eq. 2) .or.
     +          (isva2 .eq. 2 .and. isvb2 .eq. 1)) then
              if ((ilo1 .eq. 1 .and. klo2 .eq. 1) .or.
     +            (ilo1 .ne. 1 .and. klo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 1 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 1)) then
              if ((ilo1 .eq. 1 .and. jlo2 .eq. 1) .or.
     +            (ilo1 .ne. 1 .and. jlo2 .ne. 1)) then
                irev=-1
              end if
            else if ((isva2 .eq. 2 .and. isvb2 .eq. 3) .or.
     +               (isva2 .eq. 3 .and. isvb2 .eq. 2)) then
              if ((ilo1 .eq. 1 .and. ilo2 .eq. 1) .or.
     +            (ilo1 .ne. 1 .and. ilo2 .ne. 1)) then
                irev=-1
              end if
            end if
          end if
c
          if (     itransform(1) .eq. 0) then
            itransform(1)=irev*(6-abs(itransform(2))-abs(itransform(3)))
          else if (itransform(2) .eq. 0) then
            itransform(2)=irev*(6-abs(itransform(1))-abs(itransform(3)))
          else if (itransform(3) .eq. 0) then
            itransform(3)=irev*(6-abs(itransform(1))-abs(itransform(2)))
          end if
          return
          end
#else
c     this is now just a dummy routine since CGNS libs have not been
c     used
c
      write(6,*)'This code is non-functional since the installation'
      write(6,*)'of cfl3d was done without cgns libraries'
c
      stop
      end
#endif
