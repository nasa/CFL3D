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
      program v6inpswitchijk
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     This program reads a CFL3D V6 input file and outputs another input
c     file that transposes i to j, j to k, and k to i.
c     This is typically helpful if you have a 2-plane case with
c     k as the 2-D direction (CFL3D wants it to be i), and you used an
c     automatic tool such as cgns_to_cfl3dinput to create an input file
c     from the k-as-2D-direction initial case.
c
c     This program CANNOT handle input files with aeroelastics or
c     mesh deformation information or moving/sliding grids
c
c     Use the program gridswitchijk to switch the grid.
c***********************************************************************
c
      parameter(nn=1000,maxsec=9)
c
      character*80 file1,file2
      character*80 name
      character*60 name2
      character*1 dum1
      dimension i0(nn),im(nn),j0(nn),jm(nn),k0(nn),km(nn)
      dimension data(5)
c
      dimension i1_ig(maxsec,nn),i1_is(maxsec,nn),i1_ib(maxsec,nn),
     +          i1_i1(maxsec,nn),i1_i2(maxsec,nn),i1_i3(maxsec,nn),
     +          i1_i4(maxsec,nn),i1_nd(maxsec,nn)
      dimension i2_ig(maxsec,nn),i2_is(maxsec,nn),i2_ib(maxsec,nn),
     +          i2_i1(maxsec,nn),i2_i2(maxsec,nn),i2_i3(maxsec,nn),
     +          i2_i4(maxsec,nn),i2_nd(maxsec,nn)
      dimension j1_ig(maxsec,nn),j1_is(maxsec,nn),j1_ib(maxsec,nn),
     +          j1_i1(maxsec,nn),j1_i2(maxsec,nn),j1_i3(maxsec,nn),
     +          j1_i4(maxsec,nn),j1_nd(maxsec,nn)
      dimension j2_ig(maxsec,nn),j2_is(maxsec,nn),j2_ib(maxsec,nn),
     +          j2_i1(maxsec,nn),j2_i2(maxsec,nn),j2_i3(maxsec,nn),
     +          j2_i4(maxsec,nn),j2_nd(maxsec,nn)
      dimension k1_ig(maxsec,nn),k1_is(maxsec,nn),k1_ib(maxsec,nn),
     +          k1_i1(maxsec,nn),k1_i2(maxsec,nn),k1_i3(maxsec,nn),
     +          k1_i4(maxsec,nn),k1_nd(maxsec,nn)
      dimension k2_ig(maxsec,nn),k2_is(maxsec,nn),k2_ib(maxsec,nn),
     +          k2_i1(maxsec,nn),k2_i2(maxsec,nn),k2_i3(maxsec,nn),
     +          k2_i4(maxsec,nn),k2_nd(maxsec,nn)
c
      dimension data_i1(6,maxsec,nn),data_i2(6,maxsec,nn),
     +          data_j1(6,maxsec,nn),data_j2(6,maxsec,nn),
     +          data_k1(6,maxsec,nn),data_k2(6,maxsec,nn)
c
      character*80 i1_namet,i2_namet,j1_namet,j2_namet,
     +             k1_namet,k2_namet
      character*80 i1_name(maxsec,nn),
     +             i2_name(maxsec,nn),
     +             j1_name(maxsec,nn),
     +             j2_name(maxsec,nn),
     +             k1_name(maxsec,nn),
     +             k2_name(maxsec,nn)
      character*80 i1_name1(maxsec,nn),i1_name2(maxsec,nn),
     +             i2_name1(maxsec,nn),i2_name2(maxsec,nn),
     +             j1_name1(maxsec,nn),j1_name2(maxsec,nn),
     +             j2_name1(maxsec,nn),j2_name2(maxsec,nn),
     +             k1_name1(maxsec,nn),k1_name2(maxsec,nn),
     +             k2_name1(maxsec,nn),k2_name2(maxsec,nn)
c
      write(6,'('' What is input file name to read?'')')
      read(5,'(a80)') file1
      write(6,'('' What is name of new input file name to write?'')')
      read(5,'(a80)') file2
c
      open(2,file=file1,form='formatted',status='old')
      open(3,file=file2,form='formatted',status='unknown')
c
      read(2,'(a80)') name
      write(3,'(a80)') name
c   get grid name
      read(2,'(a60)') name2
      write(6,'('' original grid name = '',/,a60)') name2
      write(6,'('' input new grid name:'')')
      read(5,'(a60)') name2
      write(3,'(a60)') name2
      do n=1,12
        read(2,'(a60)') name2
        write(3,'(a60)') name2
      enddo
      read(2,'(a1)') dum1
      if (dum1 .eq. '>') then
        backspace(2)
        read(2,'(a80)') name
        write(3,'(a80)') name
        do n=1,500
          read(2,'(a1)') dum1
          if (dum1 .eq. '<') goto 1002
          backspace(2)
          read(2,'(a80)') name
          write(3,'(a80)') name
        enddo
        write(6,'('' Error, too many lines (>500) of keyword input'')')
        stop
 1002   continue
        backspace(2)
        read(2,'(a80)') name
        write(3,'(a80)') name
      else
        backspace(2)
      end if
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
c
      read(2,*) xm,al,be,re,t,ia,ih
      write(6,'('' is z up (0) or is y up (1) in the new grid?'')')
      read(5,*) ia
      write(3,'(5f10.4,2i10)') xm,al,be,re,t,ia,ih
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) sr,cr,br,xm,ym,zm
      write(3,'(6f10.4)') sr,cr,br,xm,ym,zm
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) dt,ir,if,fm,iu,cf
      write(3,'(e10.4,2i10,f10.4,i10,f10.4)') dt,ir,if,fm,iu,cf
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ngrid,nplot3d,nprint,nw,icc,i2d,nt,it
      write(3,'(8i10)') ngrid,nplot3d,nprint,nw,icc,i2d,nt,it
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      ngr=abs(ngrid)
      if(ngr .gt. nn) then
        write(6,'('' Need to increase nn to'',i5)') ngr
        stop
      end if
      do n=1,ngr
        read(2,*) nc,ie,ia,if,iv2,iv3,iv1
        if (if .eq.   0) ifx=0
        if (if .eq.  10) ifx=1
        if (if .eq.  20) ifx=2
        if (if .eq.  30) ifx=3
        if (if .eq. 100) ifx=10
        if (if .eq. 110) ifx=11
        if (if .eq. 120) ifx=12
        if (if .eq. 130) ifx=13
        if (if .eq. 200) ifx=20
        if (if .eq. 210) ifx=21
        if (if .eq. 220) ifx=22
        if (if .eq. 230) ifx=23
        if (if .eq. 300) ifx=30
        if (if .eq. 310) ifx=31
        if (if .eq. 320) ifx=32
        if (if .eq. 330) ifx=33
c
        if (if .eq.   1) ifx=100
        if (if .eq.  11) ifx=101
        if (if .eq.  21) ifx=102
        if (if .eq.  31) ifx=103
        if (if .eq. 101) ifx=110
        if (if .eq. 111) ifx=111
        if (if .eq. 121) ifx=112
        if (if .eq. 131) ifx=113
        if (if .eq. 201) ifx=120
        if (if .eq. 211) ifx=121
        if (if .eq. 221) ifx=122
        if (if .eq. 231) ifx=123
        if (if .eq. 301) ifx=130
        if (if .eq. 311) ifx=131
        if (if .eq. 321) ifx=132
        if (if .eq. 331) ifx=133
c
        if (if .eq.   2) ifx=200
        if (if .eq.  12) ifx=201
        if (if .eq.  22) ifx=202
        if (if .eq.  32) ifx=203
        if (if .eq. 102) ifx=210
        if (if .eq. 112) ifx=211
        if (if .eq. 122) ifx=212
        if (if .eq. 132) ifx=213
        if (if .eq. 202) ifx=220
        if (if .eq. 212) ifx=221
        if (if .eq. 222) ifx=222
        if (if .eq. 232) ifx=223
        if (if .eq. 302) ifx=230
        if (if .eq. 312) ifx=231
        if (if .eq. 322) ifx=232
        if (if .eq. 332) ifx=233
c
        if (if .eq.   3) ifx=300
        if (if .eq.  13) ifx=301
        if (if .eq.  23) ifx=302
        if (if .eq.  33) ifx=303
        if (if .eq. 103) ifx=310
        if (if .eq. 113) ifx=311
        if (if .eq. 123) ifx=312
        if (if .eq. 133) ifx=313
        if (if .eq. 203) ifx=320
        if (if .eq. 213) ifx=321
        if (if .eq. 223) ifx=322
        if (if .eq. 233) ifx=323
        if (if .eq. 303) ifx=330
        if (if .eq. 313) ifx=331
        if (if .eq. 323) ifx=332
        if (if .eq. 333) ifx=333
        write(3,'(7i10)') nc,ie,ia,ifx,iv1,iv2,iv3
      enddo
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) jd,kd,id
        write(3,'(3i10)') id,jd,kd
      enddo
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) i3,i4,i5,i6,i1,i2
        write(3,'(6i10)') i1,i2,i3,i4,i5,i6
      enddo
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) in,ig,i3,i4,i5,i6,i1,i2
        write(3,'(8i10)') in,ig,i1,i2,i3,i4,i5,i6
      enddo
c  Idiag section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) i2,i3,i1,i5,i6,i4
        write(3,'(6i10)') i1,i2,i3,i4,i5,i6
      enddo
c  Ifds section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) if2,if3,if1,rk2,rk3,rk1
        write(3,'(3i10,3f10.4)') if1,if2,if3,rk1,rk2,rk3
      enddo
c  Grid section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) ig,i0(n),im(n),j0(n),jm(n),k0(n),km(n),iov
        write(3,'(8i10)') ig,k0(n),km(n),i0(n),im(n),
     +                    j0(n),jm(n),iov
        if (i0(n) .gt. maxsec .or. im(n) .gt. maxsec .or.
     +      j0(n) .gt. maxsec .or. jm(n) .gt. maxsec .or.
     +      k0(n) .gt. maxsec .or. km(n) .gt. maxsec) then
          write(6,'('' sections too big... stopping'')')
          stop
        end if
      enddo
c  Grid i0 section:
      read(2,'(a80)') i1_namet
      do n=1,ngr
        do m=1,i0(n)
          read(2,*) i1_ig(m,n),i1_is(m,n),i1_ib(m,n),i1_i1(m,n),
     +              i1_i2(m,n),i1_i3(m,n),i1_i4(m,n),i1_nd(m,n)
          if (i1_nd(m,n) .gt. 0) then
            read(2,'(a80)') i1_name(m,n)
            read(2,*) (data_i1(mh,m,n),mh=1,i1_nd(m,n))
          else if(i1_nd(m,n) .lt. 0) then
            read(2,'(a80)') i1_name1(m,n)
            read(2,'(a80)') i1_name2(m,n)
          end if
        enddo
      enddo
c  Grid idim section:
      read(2,'(a80)') i2_namet
      do n=1,ngr
        do m=1,im(n)
          read(2,*) i2_ig(m,n),i2_is(m,n),i2_ib(m,n),i2_i1(m,n),
     +              i2_i2(m,n),i2_i3(m,n),i2_i4(m,n),i2_nd(m,n)
          if (i2_nd(m,n) .gt. 0) then
            read(2,'(a80)') i2_name(m,n)
            read(2,*) (data_i2(mh,m,n),mh=1,i2_nd(m,n))
          else if(i2_nd(m,n) .lt. 0) then
            read(2,'(a80)') i2_name1(m,n)
            read(2,'(a80)') i2_name2(m,n)
          end if
        enddo
      enddo
c  Grid j0 section:
      read(2,'(a80)') j1_namet
      do n=1,ngr
        do m=1,j0(n)
          read(2,*) j1_ig(m,n),j1_is(m,n),j1_ib(m,n),j1_i1(m,n),
     +              j1_i2(m,n),j1_i3(m,n),j1_i4(m,n),j1_nd(m,n)
          if (j1_nd(m,n) .gt. 0) then
            read(2,'(a80)') j1_name(m,n)
            read(2,*) (data_j1(mh,m,n),mh=1,j1_nd(m,n))
          else if(j1_nd(m,n) .lt. 0) then
            read(2,'(a80)') j1_name1(m,n)
            read(2,'(a80)') j1_name2(m,n)
          end if
        enddo
      enddo
c  Grid jdim section:
      read(2,'(a80)') j2_namet
      do n=1,ngr
        do m=1,jm(n)
          read(2,*) j2_ig(m,n),j2_is(m,n),j2_ib(m,n),j2_i1(m,n),
     +              j2_i2(m,n),j2_i3(m,n),j2_i4(m,n),j2_nd(m,n)
          if (j2_nd(m,n) .gt. 0) then
            read(2,'(a80)') j2_name(m,n)
            read(2,*) (data_j2(mh,m,n),mh=1,j2_nd(m,n))
          else if(j2_nd(m,n) .lt. 0) then
            read(2,'(a80)') j2_name1(m,n)
            read(2,'(a80)') j2_name2(m,n)
          end if
        enddo
      enddo
c  Grid k0 section:
      read(2,'(a80)') k1_namet
      do n=1,ngr
        do m=1,k0(n)
          read(2,*) k1_ig(m,n),k1_is(m,n),k1_ib(m,n),k1_i1(m,n),
     +              k1_i2(m,n),k1_i3(m,n),k1_i4(m,n),k1_nd(m,n)
          if (k1_nd(m,n) .gt. 0) then
            read(2,'(a80)') k1_name(m,n)
            read(2,*) (data_k1(mh,m,n),mh=1,k1_nd(m,n))
          else if(k1_nd(m,n) .lt. 0) then
            read(2,'(a80)') k1_name1(m,n)
            read(2,'(a80)') k1_name2(m,n)
          end if
        enddo
      enddo
c  Grid kdim section:
      read(2,'(a80)') k2_namet
      do n=1,ngr
        do m=1,km(n)
          read(2,*) k2_ig(m,n),k2_is(m,n),k2_ib(m,n),k2_i1(m,n),
     +              k2_i2(m,n),k2_i3(m,n),k2_i4(m,n),k2_nd(m,n)
          if (k2_nd(m,n) .gt. 0) then
            read(2,'(a80)') k2_name(m,n)
            read(2,*) (data_k2(mh,m,n),mh=1,k2_nd(m,n))
          else if(k2_nd(m,n) .lt. 0) then
            read(2,'(a80)') k2_name1(m,n)
            read(2,'(a80)') k2_name2(m,n)
          end if
        enddo
      enddo
c
c  Now write them out switched
c  i0
      write(3,'(a80)') i1_namet
      do n=1,ngr
        do m=1,k0(n)
          write(3,'(8i10)') k1_ig(m,n),k1_is(m,n),k1_ib(m,n),k1_i1(m,n),
     +                      k1_i2(m,n),k1_i3(m,n),k1_i4(m,n),k1_nd(m,n)
          if (k1_nd(m,n) .gt. 0) then
            write(3,'(a80)') k1_name(m,n)
            write(3,'(5f10.4)') (data_k1(mh,m,n),mh=1,k1_nd(m,n))
          else if(k1_nd(m,n) .lt. 0) then
            write(3,'(a80)') k1_name1(m,n)
            write(3,'(a80)') k1_name2(m,n)
          end if
        enddo
      enddo
c  idim
      write(3,'(a80)') i2_namet
      do n=1,ngr
        do m=1,km(n)
          write(3,'(8i10)') k2_ig(m,n),k2_is(m,n),k2_ib(m,n),k2_i1(m,n),
     +                      k2_i2(m,n),k2_i3(m,n),k2_i4(m,n),k2_nd(m,n)
          if (k2_nd(m,n) .gt. 0) then
            write(3,'(a80)') k2_name(m,n)
            write(3,'(5f10.4)') (data_k2(mh,m,n),mh=1,k2_nd(m,n))
          else if(k2_nd(m,n) .lt. 0) then
            write(3,'(a80)') k2_name1(m,n)
            write(3,'(a80)') k2_name2(m,n)
          end if
        enddo
      enddo
c  j0
      write(3,'(a80)') j1_namet
      do n=1,ngr
        do m=1,i0(n)
          write(3,'(8i10)') i1_ig(m,n),i1_is(m,n),i1_ib(m,n),i1_i3(m,n),
     +                      i1_i4(m,n),i1_i1(m,n),i1_i2(m,n),i1_nd(m,n)
          if (i1_nd(m,n) .gt. 0) then
            write(3,'(a80)') i1_name(m,n)
            write(3,'(5f10.4)') (data_i1(mh,m,n),mh=1,i1_nd(m,n))
          else if(i1_nd(m,n) .lt. 0) then
            write(3,'(a80)') i1_name1(m,n)
            write(3,'(a80)') i1_name2(m,n)
          end if
        enddo
      enddo
c  jdim
      write(3,'(a80)') j2_namet
      do n=1,ngr
        do m=1,im(n)
          write(3,'(8i10)') i2_ig(m,n),i2_is(m,n),i2_ib(m,n),i2_i3(m,n),
     +              i2_i4(m,n),i2_i1(m,n),i2_i2(m,n),i2_nd(m,n)
          if (i2_nd(m,n) .gt. 0) then
            write(3,'(a80)') i2_name(m,n)
            write(3,'(5f10.4)') (data_i2(mh,m,n),mh=1,i2_nd(m,n))
          else if(i2_nd(m,n) .lt. 0) then
            write(3,'(a80)') i2_name1(m,n)
            write(3,'(a80)') i2_name2(m,n)
          end if
        enddo
      enddo
c  k0
      write(3,'(a80)') k1_namet
      do n=1,ngr
        do m=1,j0(n)
          write(3,'(8i10)') j1_ig(m,n),j1_is(m,n),j1_ib(m,n),j1_i3(m,n),
     +                      j1_i4(m,n),j1_i1(m,n),j1_i2(m,n),j1_nd(m,n)
          if (j1_nd(m,n) .gt. 0) then
            write(3,'(a80)') j1_name(m,n)
            write(3,'(5f10.4)') (data_j1(mh,m,n),mh=1,j1_nd(m,n))
          else if(j1_nd(m,n) .lt. 0) then
            write(3,'(a80)') j1_name1(m,n)
            write(3,'(a80)') j1_name2(m,n)
          end if
        enddo
      enddo
c  kdim
      write(3,'(a80)') k2_namet
      do n=1,ngr
        do m=1,jm(n)
          write(3,'(8i10)') j2_ig(m,n),j2_is(m,n),j2_ib(m,n),j2_i3(m,n),
     +              j2_i4(m,n),j2_i1(m,n),j2_i2(m,n),j2_nd(m,n)
          if (j2_nd(m,n) .gt. 0) then
            write(3,'(a80)') j2_name(m,n)
            write(3,'(5f10.4)') (data_j2(mh,m,n),mh=1,j2_nd(m,n))
          else if(j2_nd(m,n) .lt. 0) then
            write(3,'(a80)') j2_name1(m,n)
            write(3,'(a80)') j2_name2(m,n)
          end if
        enddo
      enddo

        
c  Mseq section
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) mseq,mg,icc,mtt,ng
      write(3,'(5i10)') mseq,mg,icc,mtt,ng
c  Issc section
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) i1,f2,f3,f4,i5,f6,f7,f8
      write(3,'(i10,3f10.4,i10,3f10.4)') i1,f2,f3,f4,i5,f6,f7,f8
c  Ncyc section
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,mseq
        read(2,*) nc,mg,ne,ni
        write(3,'(4i10)') nc,mg,ne,ni
      enddo
c  Mit section
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,mseq
        read(2,'(a80)') name
        write(3,'(a80)') name
      enddo
c  1-to-1 data
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) nbli
      write(3,'(i10)') nbli
c  Limits1:
      read(2,'(a80)') name
      write(3,'('' NUMBER   GRID     :    ISTA   JSTA   KSTA  '',
     + '' IEND   JEND   KEND  ISVA1  ISVA2'')')
      do n=1,nbli
        read(2,*) num,ng,i2,i3,i1,i5,i6,i4,is1,is2
        if (is1 .eq. 2 .and. is2 .eq. 3) then
          is1=1
          is2=3
        else if (is1 .eq. 1 .and. is2 .eq. 3) then
          is1=1
          is2=2
        else if (is1 .eq. 1 .and. is2 .eq. 2) then
          is1=2
          is2=3
        else
          write(6,'('' cannot handle is1,is2... stopping'')')
          stop
        end if
        write(3,'(2i7,i14,7i7)') num,ng,i1,i2,i3,i4,i5,i6,is1,is2
      enddo
c  Limits2:
      read(2,'(a80)') name
      write(3,'('' NUMBER   GRID     :    ISTA   JSTA   KSTA  '',
     + '' IEND   JEND   KEND  ISVA1  ISVA2'')')
      do n=1,nbli
        read(2,*) num,ng,i2,i3,i1,i5,i6,i4,is1,is2
        if (is1 .eq. 2 .and. is2 .eq. 3) then
          is1=1
          is2=3
        else if (is1 .eq. 1 .and. is2 .eq. 3) then
          is1=1
          is2=2
        else if (is1 .eq. 1 .and. is2 .eq. 2) then
          is1=2
          is2=3
        else
          write(6,'('' cannot handle is1,is2... stopping'')')
          stop
        end if
        write(3,'(2i7,i14,7i7)') num,ng,i1,i2,i3,i4,i5,i6,is1,is2
      enddo
c  Patch:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ninter
      write(3,'(i10)') ninter
c  Plot3d:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      if (abs(nplot3d) .gt. 0) then
      do n=1,abs(nplot3d)
        read(2,*) ig,ip,i4,i5,i6,i7,i8,i9,i1,i2,i3
        write(3,'(11i7)') ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
      enddo
      end if
c  Movie:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) mo
      write(3,'(i7)') mo
c  Print-out:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      if (abs(nprint) .gt. 0) then
      do n=1,abs(nprint)
        read(2,*) ig,ip,i4,i5,i6,i7,i8,i9,i1,i2,i3
        write(3,'(11i7)') ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
      enddo
      end if
c  Control surface:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ncs
      write(3,'(i5)') ncs
      read(2,'(a80)') name
      write(3,'(a80)') name
      if (ncs .gt. 0) then
      do n=1,ncs
        read(2,*) ig,i3,i4,i5,i6,i1,i2,iw,inorm
        write(3,'(3i7,i9,i7,i9,3i7)') ig,i1,i2,i3,i4,i5,i6,iw,inorm
      enddo
      end if
c  Mov grid stuff - trans:
      read(2,'(a80)',end=999) name
      write(6,'('' error... cannot translate input file with'',
     +   '' moving grid info'')')
      stop
c
 999  continue
      write(6,'(/,'' new CFL3D input file successfully output'',
     + '' to '',a80)') file2
      write(6,'(/,'' MOST FILE NAMES IN NEW INPUT DECK ARE UNCHANGED'',
     + '' - BE SURE TO ALTER THEM IF DESIRED'')')
      stop
      end
