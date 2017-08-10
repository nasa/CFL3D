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
      program v6_restart_mod
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     Reads in CFL3D V6 restart file as unformatted or formatted,
c     and writes it out either same size, HALF the size, 
c     or as DOUBLE the size.
c     Going to HALF size is useful if one wishes to restart
c     from a fine grid solution and run on a COARSER level.
c     User can choose to coarsen/refine only particular index 
c     directions, if desired.  The program cannot BOTH
c     coarsen and refine different directions simultaneously.
c     For refining, the code uses CFL3D's method for determining
c     fine grid levels from coarse ones for Q's (see addx.F),
c     but uses only a simple-minded approach for BCs (qi0, etc),
c     when written.
c
c     This program CANNOT handle restart files with aeroelastics,
c     mesh deformation, or complex variables
c
c     Note: title (stored in the original restart file) is LOST
c     during conversion
c
c     If using this alone (not in conjunction with cfl3d makefile):
c     f90 -64 -r8 v6_restart_mod.f umalloc_r_r.o -o v6_restart_mod
c     On CRAY:
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 -c umalloc_r_r.F
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 v6_restart_mod.f umalloc_r_r.o -o v6_restart_mod
c*****************************************************************************
c
      parameter (nblmax=1000)
      dimension titlw(20)
      dimension jdd(nblmax),kdd(nblmax),idd(nblmax)
      dimension jddc(nblmax),kddc(nblmax),iddc(nblmax)
      character*80 file1,file2
c
      write(6,'('' input V6 restart file to read from:'')')
      read(5,'(a80)') file1
      write(6,'('' enter 0 to read formatted file'')')
      write(6,'('' enter 1 to read unformatted file'')')
      read(5,*) iunfi
      if (iunfi .eq. 1) then
        open(2,file=file1,form='unformatted',status='old')
      else
        open(2,file=file1,form='formatted',status='old')
      end if
      write(6,'('' input 1 if reading ghost values in restart'',
     + '' file (1=CFL3D default)'')')
      read(5,*) irghost
      write(6,'('' input 0 if restart file is 2nd-order-time-accurate'',
     +  '' and is an older file'')')
      write(6,'(''   (does NOT contain turb info from last time step'',
     +  '' or dt) (1=new CFL3D default)'')')
      read(5,*) itime2read
      if (itime2read .eq. 0) then
        write(6,'(/,'' Note: this program will NOT add this info.'')')
        write(6,'('' You must use keyword itime2read=0 when you run'',
     +    '' CFL3D.'')')
      end if
      write(6,'('' enter 0 to make the restart file every'',
     +  '' other (half-size)'')')
      write(6,'('' enter 1 to make it double'')')
      write(6,'('' enter -1 to keep the same size'',
     +  '' (i.e., change format only)'')')
      read(5,*) ihalf
c
      ic=0
      jc=0
      kc=0
      if (ihalf .ge. 0) then
      write(6,'('' input 3 numbers (i,j,k-directions) to'',
     . '' indicate which directions to modify'')')
      write(6,'(''   e.g., 1,1,1 means modify in all 3'',
     . '' directions,'')')
      write(6,'(''   e.g., 0,1,0 means modify in j only:'')')
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
      end if
c
      write(6,'('' input new V6 '',
     +  ''unformatted restart file name:'')')
      read(5,'(a80)') file2
      write(6,'('' enter 0 to write formatted file'')')
      write(6,'('' enter 1 to write unformatted file'')')
      read(5,*) iunfo
      if (iunfo .eq. 1) then
        open(3,file=file2,form='unformatted',status='unknown')
      else
        open(3,file=file2,form='formatted',status='unknown')
      end if
c
      write(6,'('' How many zones are in the restart file?'')')
      read(5,*) nblk
      if (nblk .gt. nblmax) then
        write(6,'('' need to increase nblmax in parameter stmt'',
     +   '' and recompile'')')
        stop
      end if
c
c     read restart file once through to get dimensions
      imax=0
      jmax=0
      kmax=0
      ntrmax=0
      do 9897 nrty=1,nblk
        if (iunfi .eq. 1) then
        read(2,end=9898) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
        else
        read(2,*,end=9898) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
        end if
        jdim1=jt-1
        kdim1=kt-1
        idim1=it-1
        jdd(nrty)=jdim1
        kdd(nrty)=kdim1
        idd(nrty)=idim1
        if (jdim1 .gt. jmax) jmax=jdim1
        if (kdim1 .gt. kmax) kmax=kdim1
        if (idim1 .gt. imax) imax=idim1
        if (ntr .gt. ntrmax) ntrmax=ntr
        if (nrty .eq. 1) then
        if (iunfi .eq. 1) then
        read(2,end=9898)  (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr)
        else
        read(2,*,end=9898)  (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr),
     .           (dum,n=1,ntr),(dum,n=1,ntr)
        end if
        end if
        if (iunfi .eq. 1) then
        read(2,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
        else
        read(2,*,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     +   l=1,5)
        end if
        if (irghost .eq. 1) then
        if (iunfi .eq. 1) then
        read(2,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
     +          ((((dum,k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     +          ((((dum,j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
        else
        read(2,*,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
     +          ((((dum,k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     +          ((((dum,j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
        end if
        end if
        if (iunfi .eq. 1) then
        read(2,end=9898) iv1,iv2,iv3
        else
        read(2,*,end=9898) iv1,iv2,iv3
        end if
        if (nrty .eq. 1) then
          if (iunfi .eq. 1) then
          read(2,end=9898) (dum,n=1,ntr),(dum,n=1,ntr),
     .            (idum,n=1,ntr),(idum,n=1,ntr)
          else
          read(2,*,end=9898) (dum,n=1,ntr),(dum,n=1,ntr),
     .            (idum,n=1,ntr),(idum,n=1,ntr)
          end if
c   ***rumsey
          if (iv1.ge.30 .or. iv2.ge.30 .or. iv3.ge.30) then
            if (iunfi .eq. 1) then
            read(2,end=9898) nummem
            else
            read(2,*,end=9898) nummem
            end if
            do l=3,nummem
              if (iunfi .eq. 1) then
              read(2,end=9898) (dum,n=1,ntr),(idum,n=1,ntr)
              else
              read(2,*,end=9898) (dum,n=1,ntr),(idum,n=1,ntr)
              end if
            enddo
          end if
c   *********
        end if
        if (iv1.ge.2 .or. iv2.ge.2 .or. iv3.ge.2) then
          if (iunfi .eq. 1) then
          read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          else
          read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          end if
          if (irghost .eq. 1) then
          if (iunfi .eq. 1) then
          read(2,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),l=1,1),
     +      m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=1,1),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=1,1),m=1,4)
          else
          read(2,*,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),l=1,1),
     +      m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=1,1),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=1,1),m=1,4)
          end if
          end if
        end if
        if (iv1.ge.4 .or. iv2.ge.4 .or. iv3.ge.4) then
          if (iunfi .eq. 1) then
          read(2,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     .      m=1,2)
          else
          read(2,*,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     .      m=1,2)
          end if
c   ***rumsey
          if (iv1.ge.30 .or. iv2.ge.30 .or. iv3.ge.30) then
            if (iunfi .eq. 1) then
            read(2,end=9898) nummem
            read(2,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     .      m=3,nummem)
            else
            read(2,*,end=9898) nummem
            read(2,*,end=9898) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),
     .      m=3,nummem)
            end if
          end if
c   *********
          if (irghost .eq. 1) then
          if (iunfi .eq. 1) then
          read(2,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),l=1,2),
     +      m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=1,2),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=1,2),m=1,4)
          else
          read(2,*,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),l=1,2),
     +      m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=1,2),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=1,2),m=1,4)
          end if
c   ***rumsey
          if (iv1.ge.30 .or. iv2.ge.30 .or. iv3.ge.30) then
            if (iunfi .eq. 1) then
            read(2,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),
     +        l=3,nummem),m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=3,nummem),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=3,nummem),m=1,4)
            else
            read(2,*,end=9898) ((((dum,j=1,jdim1+1),k=1,kdim1+1),
     +        l=3,nummem),m=1,4),
     +            ((((dum,k=1,kdim1+1),i=1,idim1),l=3,nummem),m=1,4),
     +            ((((dum,j=1,jdim1+1),i=1,idim1),l=3,nummem),m=1,4)
            end if
          end if
c   *********
          end if
          if (iunfi .eq. 1) then
          read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          else
          read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          end if
c   ***rumsey
c         if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4) then
          if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4 .or.
     +        iv1.eq.25.or. iv2.eq.25.or. iv3.eq.25) then
c   *********
            if (iunfi .eq. 1) then
            read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            else
            read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            end if
          end if
          if (iv1.eq.8 .or. iv2.eq.8 .or. iv3.eq.8  .or.
     .        iv1.eq.9 .or. iv2.eq.9 .or. iv3.eq.9  .or.
     .        iv1.eq.13.or. iv2.eq.13.or. iv3.eq.13 .or.
     .        iv1.eq.14.or. iv2.eq.14.or. iv3.eq.14) then
            if (iunfi .eq. 1) then
            read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            else
            read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
            end if
          end if
        end if
      if (ihalf .eq. 1) then
        if (ic .eq. 1) then
          if (idd(nrty) .eq. 1) then
            iddc(nrty)=idd(nrty)
          else
            iddc(nrty)=idd(nrty)*2
          end if
        else
          iddc(nrty)=idd(nrty)
        end if
        if (jc .eq. 1) then
          jddc(nrty)=jdd(nrty)*2
        else
          jddc(nrty)=jdd(nrty)
        end if
        if (kc .eq. 1) then
          kddc(nrty)=kdd(nrty)*2
        else
          kddc(nrty)=kdd(nrty)
        end if
      else if (ihalf .eq. 0) then
        if (ic .eq. 1) then
          imaxctemp=idd(nrty)/2
          iddc(nrty)=max(1,imaxctemp)
        else
          iddc(nrty)=idd(nrty)
        end if
        if (jc .eq. 1) then
          jddc(nrty)=jdd(nrty)/2
        else
          jddc(nrty)=jdd(nrty)
        end if
        if (kc .eq. 1) then
          kddc(nrty)=kdd(nrty)/2
        else
          kddc(nrty)=kdd(nrty)
        end if
      else
        iddc(nrty)=idd(nrty)
        jddc(nrty)=jdd(nrty)
        kddc(nrty)=kdd(nrty)
      end if
 9897 continue
      rewind(2)
      if (ihalf .eq. 1) then
        if (ic .eq. 1) then
          if (imax .eq. 1) then
            imaxc=imax
          else
            imaxc=imax*2
          end if
        else
          imaxc=imax
        end if
        if (jc .eq. 1) then
          jmaxc=jmax*2
        else
          jmaxc=jmax
        end if
        if (kc .eq. 1) then
          kmaxc=kmax*2
        else
          kmaxc=kmax
        end if
      else if (ihalf .eq. 0) then
        if (ic .eq. 1) then
          imaxctemp=imax/2
          imaxc=max(1,imaxctemp)
        else
          imaxc=imax
        end if
        if (jc .eq. 1) then
          jmaxc=jmax/2
        else
          jmaxc=jmax
        end if
        if (kc .eq. 1) then
          kmaxc=kmax/2
        else
          kmaxc=kmax
        end if
      else
        imaxc=imax
        jmaxc=jmax
        kmaxc=kmax
      end if
      call readwrite(imax,jmax,kmax,ntrmax,iunfi,iunfo,irghost,
     . ihalf,nblk,imaxc,jmaxc,kmaxc,ic,jc,kc,itime2read,titlw,
     . idd,jdd,kdd,iddc,jddc,kddc)
      write(6,'(/,'' successful program completion'')')
      stop
 9898 write(6,'('' Error, end of file reached... are you sure you'')')
      write(6,'('' input the correct number of zones?'')')
      stop
      end
c
c *********************************************************************
      subroutine readwrite(idim,jdim,kdim,ncycmax,iunfi,iunfo,irghost,
     . ihalf,nblk,idimc,jdimc,kdimc,ic,jc,kc,itime2read,titlw,
     . idd,jdd,kdd,iddc,jddc,kddc)
c
      integer stats
c
      dimension titlw(20)
      dimension jdd(nblk),kdd(nblk),idd(nblk)
      dimension jddc(nblk),kddc(nblk),iddc(nblk)
c 
      allocatable :: cdpw(:)
      allocatable :: cdvw(:)
      allocatable :: cdw(:)
      allocatable :: cftmomw(:)
      allocatable :: cftpw(:)
      allocatable :: cfttotw(:)
      allocatable :: cftvw(:)
      allocatable :: clw(:)
      allocatable :: cmxw(:)
      allocatable :: cmyw(:)
      allocatable :: cmzw(:)
      allocatable :: cxw(:)
      allocatable :: cyw(:)
      allocatable :: czw(:)
      allocatable :: fmdotw(:)
      allocatable :: nneg1(:)
      allocatable :: nneg2(:)
      allocatable :: qi0(:,:,:,:)
      allocatable :: qi0c(:,:,:,:)
      allocatable :: qj0(:,:,:,:)
      allocatable :: qj0c(:,:,:,:)
      allocatable :: qk0(:,:,:,:)
      allocatable :: qk0c(:,:,:,:)
      allocatable :: rms(:)
      allocatable :: rmstr1(:)
      allocatable :: rmstr2(:)
c
c     allocate memory
c
      memuse = 0
      allocate( cdpw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdpw',memuse,stats)
      allocate( cdvw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdvw',memuse,stats)
      allocate( cdw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cdw',memuse,stats)
      allocate( cftmomw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftmomw',memuse,stats)
      allocate( cftpw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftpw',memuse,stats)
      allocate( cfttotw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cfttotw',memuse,stats)
      allocate( cftvw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cftvw',memuse,stats)
      allocate( clw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'clw',memuse,stats)
      allocate( cmxw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmxw',memuse,stats)
      allocate( cmyw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmyw',memuse,stats)
      allocate( cmzw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cmzw',memuse,stats)
      allocate( cxw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cxw',memuse,stats)
      allocate( cyw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'cyw',memuse,stats)
      allocate( czw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'czw',memuse,stats)
      allocate( fmdotw(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'fmdotw',memuse,stats)
      allocate( nneg1(ncycmax), stat=stats )
      call umalloc_r(ncycmax,1,'nneg1',memuse,stats)
      allocate( nneg2(ncycmax), stat=stats )
      call umalloc_r(ncycmax,1,'nneg2',memuse,stats)
      allocate( qi0(jdim+1,kdim+1,5,4), stat=stats )
      call umalloc_r((jdim+1)*(kdim+1)*5*4,0,'qi0',memuse,stats)
      allocate( qi0c(jdimc+1,kdimc+1,5,4), stat=stats )
      call umalloc_r((jdimc+1)*(kdimc+1)*5*4,0,'qi0c',memuse,stats)
      allocate( qj0(kdim+1,idim,5,4), stat=stats )
      call umalloc_r((kdim+1)*idim*5*4,0,'qj0',memuse,stats)
      allocate( qj0c(kdimc+1,idimc,5,4), stat=stats )
      call umalloc_r((kdimc+1)*idimc*5*4,0,'qj0c',memuse,stats)
      allocate( qk0(jdim+1,idim,5,4), stat=stats )
      call umalloc_r((jdim+1)*idim*5*4,0,'qk0',memuse,stats)
      allocate( qk0c(jdimc+1,idimc,5,4), stat=stats )
      call umalloc_r((jdimc+1)*idimc*5*4,0,'qk0c',memuse,stats)
      allocate( rms(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rms',memuse,stats)
      allocate( rmstr1(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rmstr1',memuse,stats)
      allocate( rmstr2(ncycmax), stat=stats )
      call umalloc_r(ncycmax,0,'rmstr2',memuse,stats)
c
      iskip = 1
      do 9897 nrty=1,nblk
      if (iunfi .eq. 1) then
        read(2,end=9898) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
      else
        read(2,*,end=9898) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
      end if
      write(6,'('' Reading block '',i5,'':  jt,kt,it='',3i5)')
     . nrty,jt,kt,it
      if (ihalf .eq. 1) then
      if (jc .eq. 1) then
        jtc=jt*2-1
      else
        jtc=jt
      end if
      if (kc .eq. 1) then
        ktc=kt*2-1
      else
        ktc=kt
      end if
      if (ic .eq. 1) then
        if (it .eq. 2) then
          itc=2
          write(6,'('' it=2.  I-direction will not be made finer!'')')
        else
          itc=it*2-1
        end if
      else
        itc=it
      end if
      else if (ihalf .eq. 0) then
      if (jc .eq. 1) then
      if (float(jt/2) .eq. float(jt)/2.) then
        write(6,'('' jt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      jtc=(jt+1)/2
      else 
      jtc=jt
      end if
      if (kc .eq. 1) then
      if (float(kt/2) .eq. float(kt)/2.) then
        write(6,'('' kt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      ktc=(kt+1)/2
      else
      ktc=kt
      end if
      if (ic .eq. 1) then
      if (it .gt. 2 .and. (float(it/2) .eq. float(it)/2.)) then
        write(6,'('' it must be multigriddable!  Cannot proceed'')')
        stop
      end if
      if (it .eq. 2) then
        itc=2
        write(6,'('' it=2.  I-direction will not be coarsened!'')')
      else
        itc=(it+1)/2
      end if
      else
      itc=it
      end if
      else
        jtc=jt
        ktc=kt
        itc=it
      end if
c   title is lost:
      do n=1,20
        titlw(n)=0.
      enddo
c
      if (iunfo .eq. 1) then
        write(3) titlw,xmachw,jtc,ktc,itc,alphw,reuew,ntr,time
      else
        write(3,*) titlw,xmachw,jtc,ktc,itc,alphw,reuew,ntr,time
      end if
      write(6,'('' Writing block '',i5,'':  jtc,ktc,itc='',3i5)')
     . nrty,jtc,ktc,itc
c
c     Convergence data (residual,force coefficients, mass flow, etc.)
c
      jdim1=jt-1
      kdim1=kt-1
      idim1=it-1
      jdim1c=jtc-1
      kdim1c=ktc-1
      idim1c=itc-1
      if (iskip.gt.0) then
        if (iunfi .eq. 1) then
        read(2)  (rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .           (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .           (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .           (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .           (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .           (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .           (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .           (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
        else
        read(2,*)(rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .           (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .           (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .           (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .           (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .           (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .           (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .           (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
        end if
        if (iunfo .eq. 1) then
        write(3) (rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .           (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .           (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .           (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .           (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .           (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .           (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .           (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
        else
        write(3,*)(rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .           (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .           (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .           (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .           (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .           (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .           (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .           (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
        end if
      end if
c
c     Primative variables (rho,u,v,w,p)
c
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,5)
c  ghost values
      if (irghost .eq. 1) then
      if (iunfi .eq. 1) then
      read(2) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
     +        ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     +        ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
      else
      read(2,*) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
     +          ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     +          ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefinbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,
     +    ic,jc,kc,5)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoarbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,5)
      else
c       Don't coarsen
        do m=1,4
        do l=1,5
        do j=1,jdim1
        do k=1,kdim1
          qi0c(j,k,l,m)=qi0(j,k,l,m)
        enddo
        enddo
        do k=1,kdim1
        do i=1,idim1
          qj0c(k,i,l,m)=qj0(k,i,l,m)
        enddo
        enddo
        do j=1,jdim1
        do i=1,idim1
          qk0c(j,i,l,m)=qk0(j,i,l,m)
        enddo
        enddo
        enddo
        enddo
      end if
      if (iunfo .eq. 1) then
      write(3) ((((qi0c(j,k,l,m),j=1,jdim1c),k=1,kdim1c),l=1,5),m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c),i=1,idim1c),l=1,5),m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c),i=1,idim1c),l=1,5),m=1,4)
      else
      write(3,*)((((qi0c(j,k,l,m),j=1,jdim1c),k=1,kdim1c),l=1,5),m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c),i=1,idim1c),l=1,5),m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c),i=1,idim1c),l=1,5),m=1,4)
      end if
      end if
c
c     Turbulence quantities
c
      if (iunfi .eq. 1) then
        read(2) iv1,iv2,iv3
      else
        read(2,*) iv1,iv2,iv3
      end if
      if (iunfo .eq. 1) then
        write(3) iv1,iv2,iv3
      else
        write(3,*) iv1,iv2,iv3
      end if
c
c     Turbulence model convergence data
c
      if (iskip.gt.0) then
         if (iunfi .eq. 1) then
         read(2) (rmstr1(n),n=1,ntr),(rmstr2(n),n=1,ntr),
     .            (nneg1(n), n=1,ntr),(nneg2(n), n=1,ntr)
         else
         read(2,*) (rmstr1(n),n=1,ntr),(rmstr2(n),n=1,ntr),
     .            (nneg1(n), n=1,ntr),(nneg2(n), n=1,ntr)
         end if
         if (iunfo .eq. 1) then
         write(3) (rmstr1(n),n=1,ntr),(rmstr2(n),n=1,ntr),
     .            (nneg1(n), n=1,ntr),(nneg2(n), n=1,ntr)
         else
         write(3,*) (rmstr1(n),n=1,ntr),(rmstr2(n),n=1,ntr),
     .            (nneg1(n), n=1,ntr),(nneg2(n), n=1,ntr)
         end if
c   ***rumsey
         if (iv1.ge.30.or. iv2.ge.30.or. iv3.ge.30) then
           if (iunfi .eq. 1) then
             read(2) nummem
           else
             read(2,*) nummem
           end if
           if (iunfo .eq. 1) then
             write(3) nummem
           else
             write(3,*) nummem
           end if
           do l=3,nummem
             if (iunfi .eq. 1) then
               read(2) (rmstr1(n),n=1,ntr),(nneg1(n),n=1,ntr)
             else
               read(2,*) (rmstr1(n),n=1,ntr),(nneg1(n),n=1,ntr)
             end if
             if (iunfo .eq. 1) then
               write(3) (rmstr1(n),n=1,ntr),(nneg1(n),n=1,ntr)
             else
               write(3,*) (rmstr1(n),n=1,ntr),(nneg1(n),n=1,ntr)
             end if
           enddo
         end if
c   *********
      end if
c
c   vist3d data:
      if (iv1.ge.2 .or. iv2.ge.2 .or. iv3.ge.2) then
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
c  ghost values vi0,vj0,vk0
      if (irghost .eq. 1) then
      if (iunfi .eq. 1) then
      read(2) ((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),l=1,1),m=1,4),
     +        ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,1),m=1,4),
     +        ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,1),m=1,4)
      else
      read(2,*)((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),l=1,1),m=1,4),
     +         ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,1),m=1,4),
     +         ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,1),m=1,4)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefinbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,
     +    ic,jc,kc,1)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoarbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,1)
      else
c       Don't coarsen
        do m=1,4
        do l=1,1
        do j=1,jdim1
        do k=1,kdim1
          qi0c(j,k,l,m)=qi0(j,k,l,m)
        enddo
        enddo
        do k=1,kdim1
        do i=1,idim1
          qj0c(k,i,l,m)=qj0(k,i,l,m)
        enddo
        enddo
        do j=1,jdim1
        do i=1,idim1
          qk0c(j,i,l,m)=qk0(j,i,l,m)
        enddo
        enddo
        enddo
        enddo
      end if
      if (iunfo .eq. 1) then
      write(3) ((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),l=1,1),
     +  m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),l=1,1),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),l=1,1),
     +  m=1,4)
      else
      write(3,*)((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),l=1,1),
     +  m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),l=1,1),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),l=1,1),
     +  m=1,4)
      end if
      end if
      end if
c
c   tursav data:
      if (iv1.ge.4 .or. iv2.ge.4 .or. iv3.ge.4) then
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,2)
c   ***rumsey
      if (iv1.ge.30.or. iv2.ge.30.or. iv3.ge.30) then
        if (iunfi .eq. 1) then
          read(2) nummem
        else
          read(2,*) nummem
        end if
        if (iunfo .eq. 1) then
          write(3) nummem
        else
          write(3,*) nummem
        end if
        call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,nummem-3+1)
      end if
c   *********
c  ghost values ti0,tj0,tk0
      if (irghost .eq. 1) then
      if (iunfi .eq. 1) then
      read(2) ((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),l=1,2),m=1,4),
     +        ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,2),m=1,4),
     +        ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,2),m=1,4)
      else
      read(2,*)((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),l=1,2),m=1,4),
     +         ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,2),m=1,4),
     +         ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,2),m=1,4)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefinbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,
     +    ic,jc,kc,2)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoarbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,2)
      else
c       Don't coarsen
        do m=1,4
        do l=1,2
        do j=1,jdim1
        do k=1,kdim1
          qi0c(j,k,l,m)=qi0(j,k,l,m)
        enddo
        enddo
        do k=1,kdim1
        do i=1,idim1
          qj0c(k,i,l,m)=qj0(k,i,l,m)
        enddo
        enddo
        do j=1,jdim1
        do i=1,idim1
          qk0c(j,i,l,m)=qk0(j,i,l,m)
        enddo
        enddo
        enddo
        enddo
      end if
      if (iunfo .eq. 1) then
      write(3) ((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),l=1,2),
     +  m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),l=1,2),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),l=1,2),
     +  m=1,4)
      else
      write(3,*)((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),l=1,2),
     +  m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),l=1,2),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),l=1,2),
     +  m=1,4)
      end if
c   ***rumsey
      if (iv1.ge.30.or. iv2.ge.30.or. iv3.ge.30) then
      if (iunfi .eq. 1) then
      read(2) ((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),
     +        l=1,nummem-3+1),m=1,4),
     +        ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,nummem-3+1),
     +        m=1,4),
     +        ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,nummem-3+1),
     +        m=1,4)
      else
      read(2,*)((((qi0(j,k,l,m),j=1,jdim1+1),k=1,kdim1+1),
     +         l=1,nummem-3+1),m=1,4),
     +         ((((qj0(k,i,l,m),k=1,kdim1+1),i=1,idim1),l=1,nummem-3+1),
     +         m=1,4),
     +         ((((qk0(j,i,l,m),j=1,jdim1+1),i=1,idim1),l=1,nummem-3+1),
     +         m=1,4)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefinbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,
     +    ic,jc,kc,nummem-3+1)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoarbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,nummem-3+1)
      else
c       Don't coarsen
        do m=1,4
        do l=1,nummem-3+1
        do j=1,jdim1
        do k=1,kdim1
          qi0c(j,k,l,m)=qi0(j,k,l,m)
        enddo
        enddo
        do k=1,kdim1
        do i=1,idim1
          qj0c(k,i,l,m)=qj0(k,i,l,m)
        enddo
        enddo
        do j=1,jdim1
        do i=1,idim1
          qk0c(j,i,l,m)=qk0(j,i,l,m)
        enddo
        enddo
        enddo
        enddo
      end if
      if (iunfo .eq. 1) then
      write(3) ((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),
     +  l=1,nummem-3+1),m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),
     +  l=1,nummem-3+1),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),
     +  l=1,nummem-3+1),
     +  m=1,4)
      else
      write(3,*)((((qi0c(j,k,l,m),j=1,jdim1c+1),k=1,kdim1c+1),
     +  l=1,nummem-3+1),m=1,4),
     +         ((((qj0c(k,i,l,m),k=1,kdim1c+1),i=1,idim1c),
     +  l=1,nummem-3+1),
     +  m=1,4),
     +         ((((qk0c(j,i,l,m),j=1,jdim1c+1),i=1,idim1c),
     +  l=1,nummem-3+1),
     +  m=1,4)
      end if
      end if
c   *********
      end if
c
c   smin data:
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
c   xjb,xkb,blnum data
c   ***rumsey
c        if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4) then
         if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4 .or.
     +       iv1.eq.25.or. iv2.eq.25.or. iv3.eq.25) then
c   *********
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
         end if
c   cmuv data
         if (iv1.eq.8 .or. iv2.eq.8 .or. iv3.eq.8  .or.
     .       iv1.eq.9 .or. iv2.eq.9 .or. iv3.eq.9  .or.
     .       iv1.eq.13.or. iv2.eq.13.or. iv3.eq.13 .or.
     .       iv1.eq.14.or. iv2.eq.14.or. iv3.eq.14) then
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,1)
         end if
      end if
c
      iskip=0
 9897 continue
      goto 9899
 9898 continue
      write(6,'('' Error.  Reached end of file.  Are you sure'',
     . '' about the number of blocks?'')')
      stop
 9899 continue
c
c  2nd order time info in restart file (if dt was < 0, iflagg will not be here)
      if (iunfi .eq. 1) then
        read(2,end=400) iflagg
      else
        read(2,*,end=400) iflagg
      end if
      write(6,'('' iflagg='',i9)') iflagg
      if (iunfo .eq. 1) then
        write(3) iflagg
      else
        write(3,*) iflagg
      end if
      if (iflagg .eq. 1 .or. iflagg .eq. 3) then
        do 8897 nrty=1,nblk
           if (iunfi .eq. 1) then
           read(2) jt,kt,it
           else
           read(2,*) jt,kt,it
           end if
      if (ihalf .eq. 1) then
c     Make finer:
      else if (ihalf .eq. 0) then
      if (float(jt/2) .eq. float(jt)/2.) then
        write(6,'('' jt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      if (float(kt/2) .eq. float(kt)/2.) then
        write(6,'('' kt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      if (it .gt. 2 .and. (float(it/2) .eq. float(it)/2.)) then
        write(6,'('' it must be multigriddable!  Cannot proceed'')')
        stop
      end if
      jtc=(jt+1)/2
      ktc=(kt+1)/2
      if (it .eq. 2) then
        itc=2
        write(6,'('' it=2.  I-direction will not be coarsened!'')')
      else
        itc=(it+1)/2
      end if
      else
        jtc=jt
        ktc=kt
        itc=it
      end if
           if (iunfo .eq. 1) then
             write(3) jtc,ktc,itc
           else
             write(3,*) jtc,ktc,itc
           end if
           jdim1=jt-1
           kdim1=kt-1
           idim1=it-1
           jdim1c=jtc-1
           kdim1c=ktc-1
           idim1c=itc-1
c   qc0 data:
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,5)
      if (itime2read .ne. 0) then
        if (iunfi .eq. 1) then
          read(2) dtold
        else
          read(2,*) dtold
        end if
        if (iunfo .eq. 1) then
          write(3) dtold
        else
          write(3,*) dtold
        end if
c   tursav2 data:
      call getq(idd(nrty),jdd(nrty),kdd(nrty),iddc(nrty),jddc(nrty),
     +  kddc(nrty),iunfi,iunfo,ihalf,ic,jc,kc,it,2)
      end if
 8897   continue
      end if
      if (iflagg .eq. 2 .or. iflagg .eq. 3) then
        do 8898 nrty=1,nblk
           if (iunfi .eq. 1) then
             read(2) iunst
           else
             read(2,*) iunst
           end if
           if (iunfo .eq. 1) then
             write(3) iunst
           else
             write(3,*) iunst
           end if
           if (iunst .ne. 0) then
               if (iunfi .eq. 1) then
               read(2) jt,kt,it
               else
               read(2,*) jt,kt,it
               end if
      if (ihalf .eq. 1) then
c     Make finer:
      else if (ihalf .eq. 0) then
      if (float(jt/2) .eq. float(jt)/2.) then
        write(6,'('' jt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      if (float(kt/2) .eq. float(kt)/2.) then
        write(6,'('' kt must be multigriddable!  Cannot proceed'')')
        stop
      end if
      if (it .gt. 2 .and. (float(it/2) .eq. float(it)/2.)) then
        write(6,'('' it must be multigriddable!  Cannot proceed'')')
        stop
      end if
      jtc=(jt+1)/2
      ktc=(kt+1)/2
      if (it .eq. 2) then
        itc=2
        write(6,'('' it=2.  I-direction will not be coarsened!'')')
      else
        itc=(it+1)/2
      end if
      else
        jtc=jt
        ktc=kt
        itc=it
      end if
               if (iunfo .eq. 1) then
                 write(3) jtc,ktc,itc
               else
                 write(3,*) jtc,ktc,itc
               end if 
               if (iunfi .eq. 1) then
               read(2) itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt
               else
               read(2,*) itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt
               end if
               if (iunfo .eq. 1) then
               write(3) itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt
               else
               write(3,*) itrans,rfreqt,xorig,yorig,zorig,xorig0,yorig0,
     .                  zorig0,utrans,vtrans,wtrans,dxmx,dymx,dzmx,
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat,rfreqr,thetax,thetay,thetaz,
     .                  omegax,omegay,omegaz,dthxmx,dthymx,dthzmx,
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2,time2mc,dt
               end if
           end if
 8898   continue
      end if
      write(6,'(/,'' restart file is from a time-accurate type run'',
     . '' (dt > 0)'',/)')
c
      goto 401
 400  continue
      write(6,'(/,'' restart file is from a steady-state type run'',
     . '' (dt < 0)'',/)')
 401  continue
c
c   Final check to make sure there's no more unread restart data:
      if (iunfi .eq. 1) then
      read(2,end=500) dum
      else
      read(2,*,end=500) dum
      end if
      write(6,'('' WARNING!  There is still unread data at the end'',
     . '' of the file!'')')
      write(6,'('' Are you sure you specified the correct number of'',
     . '' blocks?'',/)')
      stop
c
 500  write(6,'('' New V6 restart file successfully written'',
     . '' for '',i5,'' blocks'',/)') nblk
c
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(qi0)
         deallocate(qj0)
         deallocate(qk0)
         deallocate(qi0c)
         deallocate(qj0c)
         deallocate(qk0c)
         deallocate(rms)
         deallocate(clw)
         deallocate(cdw)
         deallocate(cdpw)
         deallocate(cdvw)
         deallocate(cxw)
         deallocate(cyw)
         deallocate(czw)
         deallocate(cmxw)
         deallocate(cmyw)
         deallocate(cmzw)
         deallocate(fmdotw)
         deallocate(cftmomw)
         deallocate(cftpw)
         deallocate(cftvw)
         deallocate(cfttotw)
         deallocate(rmstr1)
         deallocate(rmstr2)
         deallocate(nneg1)
         deallocate(nneg2)
      end if
c
      return
      end
c   *******************************************************************
      subroutine makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +  idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,ll,
     +  work)
c
      dimension q(jdim,kdim,idim,ll),qc(jdimc,kdimc,idimc,ll)
      dimension work(jdimc,kdimc,idimc,ll)
c
      if (jc .eq. 1) then
        do l=1,ll
        do i=1,idim1
        do k=1,kdim1
          qc(1,k,i,l)=q(1,k,i,l)
          qc(jdim1c,k,i,l)=q(jdim1,k,i,l)
          jj=1
          do j=1,jdim1-1
            jj=jj+1
            qc(jj,k,i,l)=0.75*q(j,k,i,l)+0.25*q(j+1,k,i,l)
            jj=jj+1
            qc(jj,k,i,l)=0.25*q(j,k,i,l)+0.75*q(j+1,k,i,l)
          enddo
        enddo
        enddo
        enddo
      else
        do l=1,ll
        do i=1,idim1
        do k=1,kdim1
        do j=1,jdim1
          qc(j,k,i,l)=q(j,k,i,l)
        enddo
        enddo
        enddo
        enddo
      end if
      if (kc .eq. 1) then
        do l=1,ll
        do i=1,idim1
        do j=1,jdim1c
          work(j,1,i,l)=qc(j,1,i,l)
          work(j,kdim1c,i,l)=qc(j,kdim1,i,l)
          kk=1
          do k=1,kdim1-1
            kk=kk+1
            work(j,kk,i,l)=0.75*qc(j,k,i,l)+0.25*qc(j,k+1,i,l)
            kk=kk+1
            work(j,kk,i,l)=0.25*qc(j,k,i,l)+0.75*qc(j,k+1,i,l)
          enddo
        enddo
        enddo
        enddo
      else
        do l=1,ll
        do i=1,idim1
        do k=1,kdim1
        do j=1,jdim1c
          work(j,k,i,l)=qc(j,k,i,l)
        enddo
        enddo
        enddo
        enddo
      end if
      if (ic .eq. 1 .and. idim1c .ne. 1) then
        do l=1,ll
        do j=1,jdim1c
        do k=1,kdim1c
          qc(j,k,1,l)=work(j,k,1,l)
          qc(j,k,idim1c,l)=work(j,k,idim1,l)
          ii=1
          do i=1,idim1-1
            ii=ii+1
            qc(j,k,ii,l)=0.75*work(j,k,i,l)+0.25*work(j,k,i+1,l)
            ii=ii+1
            qc(j,k,ii,l)=0.25*work(j,k,i,l)+0.75*work(j,k,i+1,l)
          enddo
        enddo
        enddo
        enddo
      else
        do l=1,ll
        do i=1,idim1
        do k=1,kdim1c
        do j=1,jdim1c
          qc(j,k,i,l)=work(j,k,i,l)
        enddo
        enddo
        enddo
        enddo
      end if
      return
      end
c   **********************************************************
      subroutine makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +  idim1,jdim1,kdim1,ic,jc,kc,it,ll)
c
      dimension q(jdim,kdim,idim,ll),qc(jdimc,kdimc,idimc,ll)
c
      if (it .eq. 2) then
        do l=1,ll
        jj=0
        do j=1,jdim1,jc+1
          jj=jj+1
          kk=0
          do k=1,kdim1,kc+1
            kk=kk+1
            qc(jj,kk,1,l)=0.25*(q(j,k   ,1,l)+q(j+jc,k   ,1,l)+
     +                          q(j,k+kc,1,l)+q(j+jc,k+kc,1,l))
          enddo
        enddo
        enddo
      else
        do l=1,ll
        jj=0
        do j=1,jdim1,jc+1
          jj=jj+1
          kk=0
          do k=1,kdim1,kc+1
            kk=kk+1
            ii=0
            do i=1,idim1,ic+1
              ii=ii+1
             qc(jj,kk,ii,l)=0.125*(q(j,k   ,i   ,l)+q(j+jc,k   ,i   ,l)+
     +                             q(j,k+kc,i   ,l)+q(j+jc,k+kc,i   ,l)+
     +                             q(j,k   ,i+ic,l)+q(j+jc,k   ,i+ic,l)+
     +                             q(j,k+kc,i+ic,l)+q(j+jc,k+kc,i+ic,l))
            enddo
          enddo
        enddo
        enddo
      end if
      return
      end
c   ********************************************************************
      subroutine makecoarbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,ll)
c
      dimension qi0(jdim+1,kdim+1,ll,4),qj0(kdim+1,idim,ll,4),
     .          qk0(jdim+1,idim,ll,4)
      dimension qi0c(jdimc+1,kdimc+1,ll,4),qj0c(kdimc+1,idimc,ll,4),
     .          qk0c(jdimc+1,idimc,ll,4)
c
      if (it .eq. 2) then
        do m=1,4
        do l=1,ll
        jj=0
        do j=1,jdim1,jc+1
          jj=jj+1
          kk=0
          do k=1,kdim1,kc+1
            kk=kk+1
            qi0c(jj,kk,l,m)=0.25*(qi0(j,k   ,l,m)+qi0(j+jc,k   ,l,m)+
     +                            qi0(j,k+kc,l,m)+qi0(j+jc,k+kc,l,m))
          enddo
        enddo
        kk=0
        do k=1,kdim1,kc+1
          kk=kk+1
          qj0c(kk,1,l,m)=0.5*(qj0(k,1,l,m)+qj0(k+kc,1,l,m))
        enddo
        jj=0
        do j=1,jdim1,jc+1
          jj=jj+1
          qk0c(jj,1,l,m)=0.5*(qk0(j,1,l,m)+qk0(j+jc,1,l,m))
        enddo
        enddo
        enddo
      else
        do m=1,4
        do l=1,ll
        jj=0
        do j=1,jdim1,jc+1
          jj=jj+1
          kk=0
          do k=1,kdim1,kc+1
            kk=kk+1
            qi0c(jj,kk,l,m)=0.25*(qi0(j,k   ,l,m)+qi0(j+jc,k   ,l,m)+
     +                            qi0(j,k+kc,l,m)+qi0(j+jc,k+kc,l,m))
          enddo
        enddo
        ii=0
        do i=1,idim1,ic+1
          ii=ii+1
          kk=0
          do k=1,kdim1,kc+1
            kk=kk+1
            qj0c(kk,ii,l,m)=0.25*(qj0(k,i   ,l,m)+qj0(k+kc,i   ,l,m)+
     +                            qj0(k,i+ic,l,m)+qj0(k+kc,i+ic,l,m))
          enddo
        enddo
        ii=0
        do i=1,idim1,ic+1
          ii=ii+1
          jj=0
          do j=1,jdim1,jc+1
            jj=jj+1
            qk0c(jj,ii,l,m)=0.25*(qk0(j,i   ,l,m)+qk0(j+jc,i   ,l,m)+
     +                            qk0(j,i+ic,l,m)+qk0(j+jc,i+ic,l,m))
          enddo
        enddo
        enddo
        enddo
      end if
      return
      end
c   *******************************************************************
      subroutine makefinbc(qi0,qj0,qk0,qi0c,qj0c,qk0c,
     +    idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,
     +    ic,jc,kc,ll)
c
      dimension qi0(jdim+1,kdim+1,ll,4),qj0(kdim+1,idim,ll,4),
     .          qk0(jdim+1,idim,ll,4)
      dimension qi0c(jdimc+1,kdimc+1,ll,4),qj0c(kdimc+1,idimc,ll,4),
     .          qk0c(jdimc+1,idimc,ll,4)
c
c   for BCs, simple-minded approach is used to put data at finer level
c
      do m=1,4
      do l=1,ll
        kk=0
        do k=1,kdim1
          kk=kk+1
          jj=0
          do j=1,jdim1
            jj=jj+1
            qi0c(jj,kk,l,m)=qi0(j,k,l,m)
            if (jc .eq. 1) then
              jj=jj+1
              qi0c(jj,kk,l,m)=qi0(j,k,l,m)
            end if
          enddo
          if (kc .eq. 1) then
            kk=kk+1
            jj=0
            do j=1,jdim1
              jj=jj+1
              qi0c(jj,kk,l,m)=qi0(j,k,l,m)
              if (jc .eq. 1) then
                jj=jj+1
                qi0c(jj,kk,l,m)=qi0(j,k,l,m)
              end if
            enddo
          end if
        enddo
c
        ii=0
        do i=1,idim1
          ii=ii+1
          kk=0
          do k=1,kdim1
            kk=kk+1
            qj0c(kk,ii,l,m)=qj0(k,i,l,m)
            if (kc .eq. 1) then
              kk=kk+1
              qj0c(kk,ii,l,m)=qj0(k,i,l,m)
            end if
          enddo
          if (ic .eq. 1 .and. idim1c .ne. 1) then
            ii=ii+1
            kk=0
            do k=1,kdim1
              kk=kk+1
              qj0c(kk,ii,l,m)=qj0(k,i,l,m)
              if (kc .eq. 1) then
                kk=kk+1
                qj0c(kk,ii,l,m)=qj0(k,i,l,m)
              end if
            enddo
          end if
        enddo
c
        ii=0
        do i=1,idim1
          ii=ii+1
          jj=0
          do j=1,jdim1
            jj=jj+1
            qk0c(jj,ii,l,m)=qk0(j,i,l,m)
            if (jc .eq. 1) then
              jj=jj+1
              qk0c(jj,ii,l,m)=qk0(j,i,l,m)
            end if
          enddo
          if (ic .eq. 1 .and. idim1c .ne. 1) then
            ii=ii+1
            jj=0
            do j=1,jdim1
              jj=jj+1
              qk0c(jj,ii,l,m)=qk0(j,i,l,m)
              if (jc .eq. 1) then
                jj=jj+1
                qk0c(jj,ii,l,m)=qk0(j,i,l,m)
              end if
            enddo
          end if
        enddo
c
      enddo
      enddo
      return
      end
c   ***************************************************************
      subroutine getq(idim1,jdim1,kdim1,idim1c,jdim1c,
     +  kdim1c,iunfi,iunfo,ihalf,ic,jc,kc,it,inumb)
c
      integer stats
c
      allocatable :: q(:,:,:,:)
      allocatable :: qc(:,:,:,:)
      allocatable :: work(:,:,:,:)
c
      idim=idim1+1
      jdim=jdim1+1
      kdim=kdim1+1
      idimc=idim1c+1
      jdimc=jdim1c+1
      kdimc=kdim1c+1
      allocate( q(jdim,kdim,idim,inumb), stat=stats )
      call umalloc_r(jdimkdimidim*inumb,0,'q',memuse,stats)
      allocate( qc(jdimc,kdimc,idimc,inumb), stat=stats )
      call umalloc_r(jdimc*kdimc*idimc*inumb,0,'qc',memuse,stats)
      allocate( work(jdimc,kdimc,idimc,inumb), stat=stats )
      call umalloc_r(jdimc*kdimc*idimc*inumb,0,'work',memuse,stats)
c
      if (iunfi .eq. 1) then
      read(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,inumb)
      else
      read(2,*) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     +  l=1,inumb)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,inumb,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,inumb)
      else
c       Don't coarsen
        do l=1,inumb
        do j=1,jdim1
        do k=1,kdim1
        do i=1,idim1
          qc(j,k,i,l)=q(j,k,i,l)
        enddo
        enddo
        enddo
        enddo
      end if
      if (iunfo .eq. 1) then
      write(3) ((((qc(j,k,i,l),j=1,jdim1c),k=1,kdim1c),i=1,idim1c),
     +  l=1,inumb)
      else
      write(3,*) ((((qc(j,k,i,l),j=1,jdim1c),k=1,kdim1c),i=1,idim1c),
     +  l=1,inumb)
      end if
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(q)
         deallocate(qc)
         deallocate(work)
      end if
      return
      end
