!  ---------------------------------------------------------------------------
!  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
!  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
!  with point-matched, patched, overset, or embedded connectivities. Both
!  multigrid and mesh sequencing are available in time-accurate or
!  steady-state modes.
!
!  Copyright 2001 United States Government as represented by the Administrator
!  of the National Aeronautics and Space Administration. All Rights Reserved.
! 
!  The CFL3D platform is licensed under the Apache License, Version 2.0 
!  (the "License"); you may not use this file except in compliance with the 
!  License. You may obtain a copy of the License at 
!  http://www.apache.org/licenses/LICENSE-2.0. 
! 
!  Unless required by applicable law or agreed to in writing, software 
!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!  License for the specific language governing permissions and limitations 
!  under the License.
!  ---------------------------------------------------------------------------
!
program halfron

  implicit none

  character*80 fname, junk
  character*12 s
  integer i, j, n
  integer n1, n2, n3
  integer ijk(4), ijkc(3)
  integer ngrid, ninter

  write(6,*) 'Enter the name of the ronnie input file to be reduced by half '
  read(5,*) fname

  write(6,'('' input 3 numbers (i,j,k-directions) to'', &
       &'' indicate which directions to coarsen'')')
  write(6,'(''   e.g., 1,1,1 means coarsen in all 3'', &
       &'' directions,'')')
  write(6,'(''   e.g., 0,1,0 means coarsen in j only:'')')
  write(6,'(''   e.g., 0,-1,0 means refine in j only:'')')
  read(5,*) ijkc

  do n = 1, 3
     if (abs(ijkc(n)) .lt. 0 .or. abs(ijkc(n)) .gt. 1) ijkc(n)=0
  end do


  open(unit=10, file=fname, status="old", form="formatted")
  open(unit=11, file="newron.inp", status="unknown", form="formatted")

  do i = 1, 8
     read(10,'(a80)') junk
     write(11,'(a80)') junk
  end do

  read(10,*) ngrid
  write(11,'(I10)') ngrid

  ngrid = abs(ngrid)
  write(*,*) 'ngrid = ', ngrid

  read(10,'(a80)') junk
  write(11,'(a80)') junk


  do i = 1, ngrid

     read(10,*) n1, n2, ijk(1:3)

     call half3(ijk, 3, ijkc)

     if( sum(ijkc) .eq. 3 ) then
        write(11,'(5I10)') n1-1, n2, ijk(1:3)
     else if( sum(ijkc) .eq. -3 ) then
        write(11,'(5I10)') n1+1, n2, ijk(1:3)
     else        
        write(11,'(5I10)') n1, n2, ijk(1:3)
     end if

  end do


  read(10,'(a80)') junk
  write(11,'(a80)') junk

  read(10,*) ninter
  write(11,'(I10)') ninter

  read(10,'(a80)') junk
  write(11,'(a80)') junk

  do i = 1, ninter + 2

  read(10,'(a80)') junk
  write(11,'(a80)') junk

  end do
  
  do i = 1, ninter

     read(10,'(a11)',ADVANCE='NO') s
     read(10,*) ijk(1:4), n1

     read(s,*) n2, n3
     call half(ijk, 4, ijkc, n3)

     write(11,'(a11,2(I6),2(I7),I6)') s, ijk(1:4), n1
     read(s,*) n2, n3

     do j = 1, n1

        read(10,'(a11)',ADVANCE='NO') s
        read(10,*) ijk(1:4)
        read(s,*) n3
        call half(ijk, 4, ijkc, n3)
        write(11,'(a11,2(I6),2(I7))') s, ijk(1:4)

     end do

  end do

  write(*,*) 'Exiting Normally '

  close(10)
  close(11)

end program halfron



subroutine half(ijk, nmax, ijkc, code)
  implicit none

  integer, intent(IN) :: nmax, code
  integer, dimension(3), intent(IN) :: ijkc
  integer, dimension(nmax), intent(INOUT) :: ijk

  integer n, nc, ijkc12, ijkc34

  nc = mod(code / 10, 10)

  if( nc .eq. 1 ) then
     ! i face, indices are jmin, jmax, kmin, kmax

     ijkc12 = ijkc(2)
     ijkc34 = ijkc(3)

  else if( nc .eq. 2 ) then

     ! j face, indices are kmin, kmax, imin, imax

     ijkc12 = ijkc(3)
     ijkc34 = ijkc(1)

  else

     ! k face, indices are jmin, jmax, imin, imax

     ijkc12 = ijkc(2)
     ijkc34 = ijkc(1)

  end if

  if( ijkc12 .eq. 1 ) then
     ijk(1:2) = ijk(1:2) / 2 + 1
  else if (ijkc12 .eq. -1 ) then
     ijk(1:2) = (ijk(1:2)-1) * 2 + 1
  end if

  if( ijkc34 .eq. 1 ) then
     ijk(3:4) = ijk(3:4) / 2 + 1
  else if (ijkc34 .eq. -1 ) then
     ijk(3:4) = (ijk(3:4)-1) * 2 + 1
  end if


end subroutine half


subroutine half3(ijk, nmax, ijkc)
  implicit none

  integer, intent(IN) :: nmax
  integer, dimension(3), intent(IN) :: ijkc
  integer, dimension(nmax), intent(INOUT) :: ijk

  integer n

  do n = 1, 3

     if( ijkc(n) .eq. 1 ) then
        ijk(n) = ijk(n) / 2 + 1
     else if( ijkc(n) .eq. -1 ) then
        ijk(n) = (ijk(n)-1) * 2 + 1
     end if
     
  end do

end subroutine half3
