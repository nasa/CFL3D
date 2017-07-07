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
      program Get_FD
c
c     $Id$
c
c***********************************************************************
c     Purpose: read 2 cfl3d restart files and calculate finite
c     differences of force and moment coefficients; used to validate
c     complex-variable approach for determining solution derivatives
c***********************************************************************
c
      dimension titlw(20)
c
      character*80 restfile1,restfile2,fdfile
c
      common /files/ restfile1,restfile2,fdfile
c
      write(6,*)'enter first restart file to extract history data from'
      write(6,*)'this should be the "+" step file'
      read(5,'(a80)') restfile1
c
      write(6,*)'enter second restart file to extract history data from'
      write(6,*)'this should be the "-" step file'
      read(5,'(a80)') restfile2
c
      write(6,*)'enter step size'
      read(5,*) stepsize
      write(6,*)'finite diffs to be calculated with central diffs'
c
      write(6,*)'enter file name for output finite differences'
      read(5,'(a80)') fdfile
c
      write(6,*)'enter 0 to output convergence of dcy/ddv,dcmy/ddv'
      write(6,*)'enter 1 to output convergence of dcz/ddv,dcmz/ddv'
      write(6,*)'enter 2 to output convergence of all ',
     .          'force/moment derivatives'
      read(5,*) ialph
c
      open(unit=2,file=restfile1,form='unformatted',status='old')
      open(unit=3,file=restfile2,form='unformatted',status='old')
c
c.....read case title and block info, and check viability
c
      read(2) titlw,xmachw,jt,kt,it,alphw,reuew,ntr1,time
      read(3) titlw,xmachw,jt,kt,it,alphw,reuew,ntr2,time
c
      if (ntr1 .ne. ntr2) then
         write(6,*)'the two restart files have different ',
     .             'number of cycles..cannot do differences'
         stop
      end if
c
      rewind(2)
      rewind(3)
      close(2)
      close(3)
c
      ncycarg = ntr1
c
c     calculate finite differences
c
      call fdiff(ncycarg,ialph,stepsize)
c
      stop
      end
c
      subroutine fdiff(ncycarg,ialph,stepsize)
c***********************************************************************
c     Purpose: read 2 cfl3d restart files and calculate finite
c     differences of force and moment coefficients
c***********************************************************************
c
      integer stats
c
      dimension titlw(20)
c
      allocatable :: cdpw1(:)
      allocatable :: cdpw2(:)
      allocatable :: cdvw1(:)
      allocatable :: cdvw2(:)
      allocatable :: cdw1(:)
      allocatable :: cdw2(:)
      allocatable :: cftmomw1(:)
      allocatable :: cftmomw2(:)
      allocatable :: cftpw1(:)
      allocatable :: cftpw2(:)
      allocatable :: cfttotw1(:)
      allocatable :: cfttotw2(:)
      allocatable :: cftvw1(:)
      allocatable :: cftvw2(:)
      allocatable :: clw1(:)
      allocatable :: clw2(:)
      allocatable :: cmxw1(:)
      allocatable :: cmxw2(:)
      allocatable :: cmyw1(:)
      allocatable :: cmyw2(:)
      allocatable :: cmzw1(:)
      allocatable :: cmzw2(:)
      allocatable :: cxw1(:)
      allocatable :: cxw2(:)
      allocatable :: cyw1(:)
      allocatable :: cyw2(:)
      allocatable :: czw1(:)
      allocatable :: czw2(:)
      allocatable :: fmdotw1(:)
      allocatable :: fmdotw2(:)
      allocatable :: rms1(:)
      allocatable :: rms2(:)
c
      character*80 restfile1,restfile2,fdfile
c
      common /files/ restfile1,restfile2,fdfile
c
c     open files
c
      open(unit=2,file=restfile1,form='unformatted',status='old')
      open(unit=3,file=restfile2,form='unformatted',status='old')
      open(unit=8,file=fdfile,form='formatted',status='unknown')
c
      read(2) titlw,xmachw,jt,kt,it,alphw,reuew,ntr1,time
c
      if (ntr1.gt.ncycarg) then
         write(6,1239)
 1239    format(/,1x,11hstopping...,
     .          40hprevious number of iterations computed >,
     .          1x,18h dimension ncycarg)
         write(6,*)' ntr1,ncycarg = ',ntr1,ncycarg
         write(6,*)' increase value of ncycarg to at LEAST ',
     .   ntr1
      end if
c
c     allocate memory
c
      memuse = 0
      allocate( cdpw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdpw1',memuse,stats)
      allocate( cdpw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdpw2',memuse,stats)
      allocate( cdvw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdvw1',memuse,stats)
      allocate( cdvw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdvw2',memuse,stats)
      allocate( cdw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdw1',memuse,stats)
      allocate( cdw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cdw2',memuse,stats)
      allocate( cftmomw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftmomw1',memuse,stats)
      allocate( cftmomw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftmomw2',memuse,stats)
      allocate( cftpw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftpw1',memuse,stats)
      allocate( cftpw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftpw2',memuse,stats)
      allocate( cfttotw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cfttotw1',memuse,stats)
      allocate( cfttotw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cfttotw2',memuse,stats)
      allocate( cftvw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftvw1',memuse,stats)
      allocate( cftvw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cftvw2',memuse,stats)
      allocate( clw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'clw1',memuse,stats)
      allocate( clw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'clw2',memuse,stats)
      allocate( cmxw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmxw1',memuse,stats)
      allocate( cmxw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmxw2',memuse,stats)
      allocate( cmyw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmyw1',memuse,stats)
      allocate( cmyw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmyw2',memuse,stats)
      allocate( cmzw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmzw1',memuse,stats)
      allocate( cmzw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cmzw2',memuse,stats)
      allocate( cxw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cxw1',memuse,stats)
      allocate( cxw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cxw2',memuse,stats)
      allocate( cyw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cyw1',memuse,stats)
      allocate( cyw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'cyw2',memuse,stats)
      allocate( czw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'czw1',memuse,stats)
      allocate( czw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'czw2',memuse,stats)
      allocate( fmdotw1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'fmdotw1',memuse,stats)
      allocate( fmdotw2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'fmdotw2',memuse,stats)
      allocate( rms1(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'rms1',memuse,stats)
      allocate( rms2(ncycarg), stat=stats )
      call umalloc_r(ncycarg,0,'rms2',memuse,stats)
c
      read(2) (rms1(n),     n=1,ntr1),(clw1(n),     n=1,ntr1),
     .        (cdw1(n),     n=1,ntr1),(cdpw1(n),    n=1,ntr1),
     .        (cdvw1(n),    n=1,ntr1),(cxw1(n),     n=1,ntr1),
     .        (cyw1(n),     n=1,ntr1),(czw1(n),     n=1,ntr1),
     .        (cmxw1(n),    n=1,ntr1),(cmyw1(n),    n=1,ntr1),
     .        (cmzw1(n),    n=1,ntr1),(fmdotw1(n),  n=1,ntr1),
     .        (cftmomw1(n), n=1,ntr1),(cftpw1(n),   n=1,ntr1),
     .        (cftvw1(n),   n=1,ntr1),(cfttotw1(n), n=1,ntr1)
c
      read(3) titlw,xmachw,jt,kt,it,alphw,reuew,ntr2,time
c
      if (ntr2.gt.ncycarg) then
         write(6,1239)
         write(6,*)' ntr2,ncycarg = ',ntr2,ncycarg
         write(6,*)' increase value of ncycarg to at LEAST ',
     .   ntr2
      end if
c
      read(3) (rms2(n),     n=1,ntr2),(clw2(n),     n=1,ntr2),
     .        (cdw2(n),     n=1,ntr2),(cdpw2(n),    n=1,ntr2),
     .        (cdvw2(n),    n=1,ntr2),(cxw2(n),     n=1,ntr2),
     .        (cyw2(n),     n=1,ntr2),(czw2(n),     n=1,ntr2),
     .        (cmxw2(n),    n=1,ntr2),(cmyw2(n),    n=1,ntr2),
     .        (cmzw2(n),    n=1,ntr2),(fmdotw2(n),  n=1,ntr2),
     .        (cftmomw2(n), n=1,ntr2),(cftpw2(n),   n=1,ntr2),
     .        (cftvw2(n),   n=1,ntr2),(cfttotw2(n), n=1,ntr2)
c
      write(8,'('' '')')
      write(8,'(''***********************************************'',
     .  ''********************************'')')
      write(8,'(''                     derivatives via central '',
     . ''differences '')')
      write(8,'(''                        with step size = '',
     . e12.5)') stepsize
      write(8,'(''***********************************************'',
     .  ''********************************'')')
      write(8,'('' '')')
      write(8,2) (titlw(i),i=1,20)
    2 format(20a4)
      write(8,'('' Mach='',e12.4,'', alpha='',e12.4,
     . '', ReUe='',e12.4)') xmachw,alphw,reuew
      write(8,'('' Final d/d() of cl,cd       ='',2e13.5)')
     .          (clw1(ntr1)-clw2(ntr2))/stepsize/2.,
     .          (cdw1(ntr1)-cdw2(ntr2))/stepsize/2.
      write(8,'('' Final d/d() of cx,cy,cz    ='',3e13.5)')
     .          (cxw1(ntr1)-cxw2(ntr2))/stepsize/2.,
     .          (cyw1(ntr1)-cyw2(ntr2))/stepsize/2.,
     .          (czw1(ntr1)-czw2(ntr2))/stepsize/2.
      write(8,'('' Final d/d() of cmx,cmy,cmz ='',3e13.5)')
     .          (cmxw1(ntr1)-cmxw2(ntr2))/stepsize/2.,
     .          (cmyw1(ntr1)-cmyw2(ntr2))/stepsize/2.,
     .          (cmzw1(ntr1)-cmzw2(ntr2))/stepsize/2. 
c
c     output convergence history
c
      write(8,'('' '')')
      if (ialph.eq.0) then
         write(8,'(''  iter     residual      dcl/ddv      dcd/ddv'',
     .             ''      dcy/ddv     dcmy/ddv'')')
      else if (ialph.eq.0) then
         write(8,'(''  iter     residual      dcl/ddv      dcd/ddv'',
     .             ''      dcz/ddv     dcmz/ddv'')')
      else
         write(8,'(''    it     residual    d(cl)/d()'',
     .   ''    d(cd)/d()    d(cx)/d()    d(cy)/d()    d(cz)/d()'',
     .   ''   d(cmx)/d()   d(cmy)/d()   d(cmz)/d()'')')
         write(8,'('' '')')
      end if
c
c     add dummy residual (=1) so that file has same number of
c     variables as regular cfl3d.res file
c
      if (ialph.eq.0) then
         do n=1,ntr1
            write(8,'(i6,5e13.5)') n,1.,
     .                             (clw1(n)-clw2(n))/stepsize/2.,
     .                             (cdw1(n)-cdw2(n))/stepsize/2.,
     .                             (cyw1(n)-cyw2(n))/stepsize/2.,
     .                             (cmyw1(n)-cmyw2(n))/stepsize/2.
        end do
      else if (ialph.eq.1) then
         do n=1,ntr1
            write(8,'(i6,5e13.5)') n,1.,
     .                             (clw1(n)-clw2(n))/stepsize/2.,
     .                             (cdw1(n)-cdw2(n))/stepsize/2.,
     .                             (czw1(n)-czw2(n))/stepsize/2.,
     .                             (cmzw1(n)-cmzw2(n))/stepsize/2.
        end do
      else
         do n=1,ntr1
            write(8,'(i6,9e13.5)') n,1.,
     .                             (clw1(n)-clw2(n))/stepsize/2.,
     .                             (cdw1(n)-cdw2(n))/stepsize/2.,
     .                             (cxw1(n)-cxw2(n))/stepsize/2.,
     .                             (cyw1(n)-cyw2(n))/stepsize/2.,
     .                             (czw1(n)-czw2(n))/stepsize/2.,
     .                             (cmxw1(n)-cmxw2(n))/stepsize/2.,
     .                             (cmyw1(n)-cmyw2(n))/stepsize/2.,
     .                             (cmzw1(n)-cmzw2(n))/stepsize/2.
        end do
      end if
c
c     free memory
c
      deallocate(rms1)
      deallocate(clw1)
      deallocate(cdw1)
      deallocate(cdpw1)
      deallocate(cdvw1)
      deallocate(cxw1)
      deallocate(cyw1)
      deallocate(czw1)
      deallocate(cmxw1)
      deallocate(cmyw1)
      deallocate(cmzw1)
      deallocate(fmdotw1)
      deallocate(cftmomw1)
      deallocate(cftpw1)
      deallocate(cftvw1)
      deallocate(cfttotw1)
      deallocate(rms2)
      deallocate(clw2)
      deallocate(cdw2)
      deallocate(cdpw2)
      deallocate(cdvw2)
      deallocate(cxw2)
      deallocate(cyw2)
      deallocate(czw2)
      deallocate(cmxw2)
      deallocate(cmyw2)
      deallocate(cmzw2)
      deallocate(fmdotw2)
      deallocate(cftmomw2)
      deallocate(cftpw2)
      deallocate(cftvw2)
      deallocate(cfttotw2)
c
      return
      end
