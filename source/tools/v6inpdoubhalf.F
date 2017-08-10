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
      program v6inpdoubhalf
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     This program reads a CFL3D V6 input file and outputs another input file
c     that is either for double the grid size or half the grid size.
c
c     This is convenient when, for example, one has a fine grid and would like
c     to run 1/2 of the grid, and avoid working in a big memory queue (if one
c     uses the fine grid with mesh sequencing option where all the iterations 
c     are done on the coarser grid and none on the fine grid, the code runs
c     the coarser grid only, but also requires more memory than if one could
c     run on the coarse grid only).
c     User can also choose to coarsen/double only particular index directions,
c     if desired.
c
c     This program CANNOT handle input files with aeroelastics or
c     mesh deformation information
c
c     Use the program everyother_xyz.f to make a coarser PLOT3D-type grid.
c***********************************************************************
c
      parameter(nn=10000)
c
      character*80 file1,file2
      character*80 name
      character*60 name2
      character*1 dum1
      dimension i0(nn),im(nn),j0(nn),jm(nn),k0(nn),km(nn)
      dimension data(5)
c
      write(6,'('' What is input file name to read?'')')
      read(5,'(a80)') file1
      write(6,'('' What is name of new input file name to write?'')')
      read(5,'(a80)') file2
      write(6,'('' Input 0=halve input file, 1=double:'')')
      read(5,*) idouble
      if (idouble .ne. 0) then
        idouble=1
        write(6,'('' Doubling the input file'')')
      else
        write(6,'('' Halving the input file'')')
      end if
c
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
        write(6,'('' modifying in i-direction (but if it is a'',
     +   '' 2-D case, it is not changed)'')')
      end if
      if (jc .eq. 1) then
        write(6,'('' modifying in j-direction'')')
      end if
      if (kc .eq. 1) then
        write(6,'('' modifying in k-direction'',/)')
      end if
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
      if (i2d .eq. 1) ic=0
      if (ic.eq.0 .and. jc.eq.0 .and. kc.eq.0) then
        write(6,'('' Nothing is changing!  No action taken'')')
        stop
      end if
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      ngr=abs(ngrid)
      if(ngr .gt. nn) then
        write(6,'('' Need to increase nn to'',i5)') ngr
        stop
      end if
      do n=1,ngr
        read(2,*) nc,ie,ia,if,iv1,iv2,iv3
        if (idouble .eq. 1) then
          nc=nc+1
          if (n .eq. 1) then
            write(6,'('' WARNING:  changing ncg to '',i5)') nc
          end if
        else
          if (nc .ne. 0) then
            nc=nc-1
            if (n .eq. 1) then
              write(6,'('' WARNING:  changing ncg to '',i5)') nc
            end if
          end if
        end if
        write(3,'(7i10)') nc,ie,ia,if,iv1,iv2,iv3
      enddo
c
      i2d=abs(i2d)
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) id,jd,kd
        if (idouble .eq. 1) then
          if (ic .eq. 1) id=id*2-1
          if (i2d .eq. 1) id=2
          if (jc .eq. 1) jd=jd*2-1
          if (kc .eq. 1) kd=kd*2-1
        else
          if (i2d .ne. 1) then
            if (ic .eq. 1) then
            if (float(id/2) .eq. (float(id)/2.)) then
              write(6,'('' idim not MGable, zone='',i5)') n
              write(6,'('' new input file only partially written'')')
              stop
            end if
            id=(id+1)/2
            end if
          end if
          if (jc .eq. 1) then
          if (float(jd/2) .eq. (float(jd)/2.)) then
            write(6,'('' jdim not MGable, zone='',i5)') n
            write(6,'('' new input file only partially written'')')
            stop
          end if
          jd=(jd+1)/2
          end if
          if (kc .eq. 1) then
          if (float(kd/2) .eq. (float(kd)/2.)) then
            write(6,'('' kdim not MGable, zone='',i5)') n
            write(6,'('' new input file only partially written'')')
            stop
          end if
          kd=(kd+1)/2
          end if
        end if
        write(3,'(3i10)') id,jd,kd
      enddo
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) i1,i2,i3,i4,i5,i6
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
            if (i5 .ne. 0) then
              if (kc .eq. 1) i5=i5*2-1
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) i6=i6*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ilamlo not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' ilamhi not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' jlamlo not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jlamhi not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
            if (i5 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' klamlo not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i6/2) .eq. (float(i6)/2.)) then
                write(6,'('' klamhi not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i6=(i6+1)/2
              end if
            end if
          end if
        write(3,'(6i10)') i1,i2,i3,i4,i5,i6
      enddo
c
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) in,ig,i1,i2,i3,i4,i5,i6
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) i2=i2*2-1
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i4=i4*2-1
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) i5=i5*2-1
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) i6=i6*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' IS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' JS not MGable, zone='',i5)') n 
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' KS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' IE not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' JE not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i6/2) .eq. (float(i6)/2.)) then
                write(6,'('' KE not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i6=(i6+1)/2
              end if
            end if
          end if
        write(3,'(8i10)') in,ig,i1,i2,i3,i4,i5,i6
      enddo
c  Idiag section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) i1,i2,i3,i4,i5,i6
        write(3,'(6i10)') i1,i2,i3,i4,i5,i6
      enddo
c  Ifds section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) if1,if2,if3,rk1,rk2,rk3
        write(3,'(3i10,3f10.4)') if1,if2,if3,rk1,rk2,rk3
      enddo
c  Grid section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        read(2,*) ig,i0(n),im(n),j0(n),jm(n),k0(n),km(n),iov
        write(3,'(8i10)') ig,i0(n),im(n),j0(n),jm(n),
     +                    k0(n),km(n),iov
      enddo
c  Grid i0 section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,i0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (jc .eq. 1) i1=i1*2-1
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) i2=i2*2-1
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' jsta/I0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' jend/I0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta/I0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' kend/I0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Grid idim section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,im(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (jc .eq. 1) i1=i1*2-1
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) i2=i2*2-1
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' jsta/IDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' jend/IDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta/IDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' kend/IDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Grid j0 section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,j0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/J0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/J0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta/J0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' kend/J0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Grid jdim section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,jm(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/JDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/JDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta/JDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' kend/JDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Grid k0 section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,k0(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/K0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/K0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' jsta/K0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jend/K0 not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Grid kdim section:
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ngr
        do m=1,km(n)
          read(2,*) ig,is,ib,i1,i2,i3,i4,nd
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/KDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/KDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' jsta/KDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jend/KDIM not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
          end if
          write(3,'(8i10)') ig,is,ib,i1,i2,i3,i4,nd
          if (nd .gt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,*) (data(mh),mh=1,nd)
            write(3,'(5f10.4)') (data(mh),mh=1,nd)
          else if(nd .lt. 0) then
            read(2,'(a80)') name
            write(3,'(a80)') name
            read(2,'(a80)') name
            write(3,'(a80)') name
          end if
        enddo
      enddo
c  Mseq section
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) mseq,mg,icc,mtt,ng
      write(3,'(5i10)') mseq,mg,icc,mtt,ng
      if (mseq .gt. 1) then
        write(6,'('' WARNING:  mseq currently > 1.  It is NOT'',
     +   '' being changed!'')')
      end if
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
        if (idouble .eq. 1) then
          mg=mg+1
          write(6,'('' WARNING:  changing mglevg to '',i5)') mg
        else
          if (mg .gt. 1) then
            mg=mg-1
            write(6,'('' WARNING:  changing mglevg to '',i5)') mg
          end if
        end if
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
        read(2,*) num,ng,i1,i2,i3,i4,i5,i6,is1,is2
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) i2=i2*2-1
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i4=i4*2-1
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) i5=i5*2-1
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) i6=i6*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' jsta not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' iend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' jend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i6/2) .eq. (float(i6)/2.)) then
                write(6,'('' kend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i6=(i6+1)/2
              end if
            end if
          end if
        write(3,'(2i7,i14,7i7)') num,ng,i1,i2,i3,i4,i5,i6,is1,is2
      enddo
c  Limits2:
      read(2,'(a80)') name
      write(3,'('' NUMBER   GRID     :    ISTA   JSTA   KSTA  '',
     + '' IEND   JEND   KEND  ISVA1  ISVA2'')')
      do n=1,nbli
        read(2,*) num,ng,i1,i2,i3,i4,i5,i6,is1,is2
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) i2=i2*2-1
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i4=i4*2-1
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) i5=i5*2-1
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) i6=i6*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' jsta not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
            end if
            if (i3 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' ksta not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' iend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' jend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i6/2) .eq. (float(i6)/2.)) then
                write(6,'('' kend not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i6=(i6+1)/2
              end if
            end if
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
        read(2,*) ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) i5=i5*2-1
            end if
            if (i7 .ne. 0) then
              if (kc .eq. 1) i7=i7*2-1
            end if
            if (i8 .ne. 0) then
              if (kc .eq. 1) i8=i8*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jsta/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' jend/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i7 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i7/2) .eq. (float(i7)/2.)) then
                write(6,'('' ksta/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i7=(i7+1)/2
              end if
            end if
            if (i8 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i8/2) .eq. (float(i8)/2.)) then
                write(6,'('' kend/PLOT3D not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i8=(i8+1)/2
              end if
            end if
          end if
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
        read(2,*) ig,ip,i1,i2,i3,i4,i5,i6,i7,i8,i9
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) i5=i5*2-1
            end if
            if (i7 .ne. 0) then
              if (kc .eq. 1) i7=i7*2-1
            end if
            if (i8 .ne. 0) then
              if (kc .eq. 1) i8=i8*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jsta/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
            if (i5 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' jend/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i7 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i7/2) .eq. (float(i7)/2.)) then
                write(6,'('' ksta/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i7=(i7+1)/2
              end if
            end if
            if (i8 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i8/2) .eq. (float(i8)/2.)) then
                write(6,'('' kend/PRINT not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i8=(i8+1)/2
              end if
            end if
          end if
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
        read(2,*) ig,i1,i2,i3,i4,i5,i6,iw,inorm
          if (idouble .eq. 1) then
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i1=i1*2-1
              end if
            end if
            if (i2 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) i2=i2*2-1
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) i3=i3*2-1
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) i4=i4*2-1
            end if
            if (i5 .ne. 0) then
              if (kc .eq. 1) i5=i5*2-1
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) i6=i6*2-1
            end if
          else
            if (i1 .ne. 0) then
              if (i2d .ne. 1) then
              if (ic .eq. 1) then
              if (float(i1/2) .eq. (float(i1)/2.)) then
                write(6,'('' ista/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i1=(i1+1)/2
              end if
              end if
            end if
            if (i2 .ne. 0) then
              if (ic .eq. 1) then
              if (i2d .ne. 1) then
              if (float(i2/2) .eq. (float(i2)/2.)) then
                write(6,'('' iend/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i2=(i2+1)/2
              end if
              end if
            end if
            if (i3 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i3/2) .eq. (float(i3)/2.)) then
                write(6,'('' jsta/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i3=(i3+1)/2
              end if
            end if
            if (i4 .ne. 0) then
              if (jc .eq. 1) then
              if (float(i4/2) .eq. (float(i4)/2.)) then
                write(6,'('' jend/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i4=(i4+1)/2
              end if
            end if
            if (i5 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i5/2) .eq. (float(i5)/2.)) then
                write(6,'('' ksta/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i5=(i5+1)/2
              end if
            end if
            if (i6 .ne. 0) then
              if (kc .eq. 1) then
              if (float(i6/2) .eq. (float(i6)/2.)) then
                write(6,'('' kend/NCS not MGable, zone='',i5)') n
                write(6,'('' new input file only partially written'')')
                stop
              end if
              i6=(i6+1)/2
              end if
            end if
          end if
        write(3,'(3i7,i9,i7,i9,3i7)') ig,i1,i2,i3,i4,i5,i6,iw,inorm
      enddo
      end if
c  Mov grid stuff - trans:
      read(2,'(a80)',end=999) name
      if (name .eq. '') then
        write(6,'(/,'' WARNING:  input file has a blank line after'',
     +   '' control surface info...'')')
        write(6,'(''           any subsequent lines will be ignored'')')
        goto 999
      end if
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ntrans
      write(3,'(i8)') ntrans
      read(2,'(a80)') name
      write(3,'(a80)') name
      if (ntrans .ne. 0) then
        read(2,*) xref
        write(3,'(f8.4)') xref
      end if
      read(2,'(a80)') name
      write(3,'('' GRID ITRANS     RFREQ      XMAG      YMAG'',
     +  ''      ZMAG'')')
      if (ntrans .ne. 0) then
        do n=1,ntrans
          read(2,*) ig,it,fr,x,y,z
          if (abs(x).gt.99. .or. abs(y).gt.99. .or. abs(z).gt.99. .or.
     +        abs(fr).gt.99.) then
            write(6,'('' Check NTRANS data section to make sure the'',
     +       '' numbers fit with given format!'')')
          end if
          write(3,'(i5,i7,4f10.5)') ig,it,fr,x,y,z
        enddo
      end if
      read(2,'(a80)') name
      write(3,'('' GRID     DXMAX     DYMAX     DZMAX'')')
      if (ntrans .ne. 0) then
        do n=1,ntrans
          read(2,*) ig,x,y,z
          write(3,'(i5,3f10.4)') ig,x,y,z
        enddo
      end if
c  Mov grid stuff - rotation:
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ntrans
      write(3,'(i8)') ntrans
      read(2,'(a80)') name
      write(3,'(a80)') name
      if (ntrans .ne. 0) then
        read(2,*) xref
        write(3,'(f8.4)') xref
      end if
      read(2,'(a80)') name
      write(3,'('' GRID IROTAT     RFREQ    THXMAG    THYMAG   '',
     +  '' THZMAG    XORIG    YORIG    ZORIG'')')
      if (ntrans .ne. 0) then
        do n=1,ntrans
          read(2,*) ig,ir,fr,thx,thy,thz,x,y,z
          if (abs(x).gt.99. .or. abs(y).gt.99. .or. abs(z).gt.99. .or.
     +        abs(thx).gt.99. .or. abs(thy).gt.99. .or. 
     +        abs(thz).gt.99. .or. abs(fr).gt.99) then
            write(6,'('' Check NROTAT data section to make sure the'',
     +       '' numbers fit with given format!'')')
          end if
          write(3,'(i5,i7,4f10.5,3f9.4)') ig,ir,fr,thx,thy,thz,x,y,z
        enddo
      end if
      read(2,'(a80)') name
      write(3,'('' GRID    THXMAX    THYMAX    THZMAX'')')
      if (ntrans .ne. 0) then
        do n=1,ntrans
          read(2,*) ig,x,y,z
          write(3,'(i5,3f10.4)') ig,x,y,z
        enddo
      end if
c  Dynamic patch data
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) ninter
      write(3,'(i10)') ninter
      read(2,'(a80)') name
      write(3,'(a80)') name
      do n=1,ninter
        read(2,*) i1,i2,i3,i4,i5,i6,i7,i8,i9
        write(3,'(2i6,7i9)') i1,i2,i3,i4,i5,i6,i7,i8,i9
      enddo
c
      write(6,'('' WARNING:  currently, no changes are made to dynamic'',
     + '' patch input data'')')
c
      do n=1,ninter
      read(2,'(a80)') name
      write(3,'(a80)') name
      read(2,*) i1,i2,i3,i4,i5,i6,nfb
      write(3,'(2i6,5i9)') i1,i2,i3,i4,i5,i6,nfb
      do m=1,nfb
        read(2,'(a80)') name
        write(3,'(a80)') name
        read(2,*) i1,i2,i3,i4,i5,f1,f2
        write(3,'(i12,4i9,2f9.3)') i1,i2,i3,i4,i5,f1,f2
        read(2,'(a80)') name
        write(3,'(a80)') name
        read(2,*) x,y,z,dx,dy,dz
        write(3,'(f21.3,5f9.3)') x,y,z,dx,dy,dz
      enddo
      enddo
c
 999  continue
      write(6,'(/,'' new CFL3D input file successfully output'',
     + '' to '',a80)') file2
      write(6,'(/,'' MOST FILE NAMES IN NEW INPUT DECK ARE UNCHANGED'',
     + '' - BE SURE TO ALTER THEM IF DESIRED'')')
      stop
      end
