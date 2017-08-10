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
      program cfl3dinp_to_FVBND
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     This program reads a CFL3D V6 input file and creates a Fieldview
c     *.fvbnd file (useful for finding and plotting viscous/inviscid 
c     walls). The *.fvbnd file takes the root name of the grid by 
c     default.
c***********************************************************************
c
      parameter(nn=10000)
c
      character*80 file1,file2
      character*80 name
      character*60 name2
      character*1 dum1
      dimension i0(nn),im(nn),j0(nn),jm(nn),k0(nn),km(nn)
      dimension idx(nn),jdx(nn),kdx(nn)
      dimension data(5)
c
      write(6,'('' What is input file name to read?'')')
      read(5,'(a80)') file1
c
      open(2,file=file1,form='formatted',status='old')
c
      read(2,'(a80)') name
c   get grid name
      read(2,'(a60)') name2
      file2=trim(name2)//'.fvbnd'
c
      open(3,file=file2,form='formatted',status='unknown')
c
      do n=1,12
        read(2,'(a60)') name2
      enddo
      read(2,'(a1)') dum1
      if (dum1 .eq. '>') then
        backspace(2)
        read(2,'(a80)') name
        do n=1,500
          read(2,'(a1)') dum1
          if (dum1 .eq. '<') goto 1002
          backspace(2)
          read(2,'(a80)') name
        enddo
        write(6,'('' Error, too many lines (>500) of keyword input'')')
        stop
 1002   continue
        backspace(2)
        read(2,'(a80)') name
      else
        backspace(2)
      end if
      read(2,'(a80)') name
      read(2,'(a80)') name
c
      read(2,*) xm,al,be,re,t,ia,ih
      read(2,'(a80)') name
      read(2,*) sr,cr,br,xm,ym,zm
      read(2,'(a80)') name
      read(2,*) dt,ir,if,fm,iu,cf
      read(2,'(a80)') name
      read(2,*) ngrid,nplot3d,nprint,nw,icc,i2d,nt,it
      if (i2d .eq. 1) ic=0
c
      read(2,'(a80)') name
      ngr=abs(ngrid)
      if(ngr .gt. nn) then
        write(6,'('' Need to increase nn to'',i5)') ngr
        stop
      end if
      do n=1,ngr
        read(2,*) nc,ie,ia,if,iv1,iv2,iv3
      enddo
c
      i2d=abs(i2d)
c
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) id,jd,kd
        idx(n)=id
        jdx(n)=jd
        kdx(n)=kd
      enddo
c
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) i1,i2,i3,i4,i5,i6
      enddo
c
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) in,ig,i1,i2,i3,i4,i5,i6
      enddo
c  Idiag section:
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) i1,i2,i3,i4,i5,i6
      enddo
c  Ifds section:
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) if1,if2,if3,rk1,rk2,rk3
      enddo
c  Grid section:
      read(2,'(a80)') name
      do n=1,ngr
        read(2,*) ig,i0(n),im(n),j0(n),jm(n),k0(n),km(n),iov
      enddo
c
c  *** write header for .fvbnd file ***
      write(3,'(''FVBND 1 4'')')
      write(3,'(''viscous wall'')')
      write(3,'(''inviscid wall'')')
      write(3,'(''BOUNDARIES'')')
c  ************************************
c  Grid i0 section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,i0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        indx0,indx0,i1,i2,i3,i4,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        indx0,indx0,i1,i2,i3,i4,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Grid idim section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,im(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=idx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        indx0,indx0,i1,i2,i3,i4,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=idx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        indx0,indx0,i1,i2,i3,i4,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Grid j0 section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,j0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,indx0,indx0,i3,i4,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,indx0,indx0,i3,i4,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Grid jdim section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,jm(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=jdx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,indx0,indx0,i3,i4,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=jdx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,indx0,indx0,i3,i4,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Grid k0 section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,k0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,i3,i4,indx0,indx0,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=1
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,i3,i4,indx0,indx0,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Grid kdim section:
      read(2,'(a80)') name
      do n=1,ngr
        do m=1,km(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
c  *** write to fvbnd file ***
          if (abs(ib) .eq. 2004 .or. abs(ib) .eq. 2014 .or.
     +        abs(ib) .eq. 2024 .or. abs(ib) .eq. 2016) then
            ibd=1
            indx0=kdx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,i3,i4,indx0,indx0,izero
          end if
          if (abs(ib) .eq. 1005 .or. abs(ib) .eq. 1006) then
            ibd=2
            indx0=kdx(n)
            izero=0
            write(3,'(8i5,''  F '',i5)') ibd,n,
     +        i1,i2,i3,i4,indx0,indx0,izero
          end if
c  ***************************
          if (nd .gt. 0) then
            read(2,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            read(2,'(a80)') name
          end if
        enddo
      enddo
c  Mseq section
      read(2,'(a80)') name
      read(2,*) mseq,mg,icc,mtt,ng
c  Issc section
      read(2,'(a80)') name
      read(2,*) i1,f2,f3,f4,i5,f6,f7,f8
c  Ncyc section
      read(2,'(a80)') name
      do n=1,mseq
        read(2,*) nc,mg,ne,ni
      enddo
c  Mit section
      read(2,'(a80)') name
      do n=1,mseq
        read(2,'(a80)') name
      enddo
c  1-to-1 data
      read(2,'(a80)') name
      read(2,'(a80)') name
      read(2,*) nbli
c  Limits1:
      read(2,'(a80)') name
      do n=1,nbli
        read(2,*) num,ng,i1,i2,i3,i4,i5,i6,is1,is2
      enddo
c  Limits2:
      read(2,'(a80)') name
      do n=1,nbli
        read(2,*) num,ng,i1,i2,i3,i4,i5,i6,is1,is2
      enddo
c
      write(6,'('' new file written: '',a80)') file2
      stop
      end
