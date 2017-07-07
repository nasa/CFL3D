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
      program initialize_field
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     Creates a CFL3D V6 restart file with user-specified initial conditions.
c     This can be used when "freestream" initial conditions (which CFL3D does
c     by default) are not desired.
c     After creating the restart file with this tool, run CFL3D with IREST
c     equal to 1 or -1.
c     This tool does not allow setting of initial turbulence quantities.
c     It also sets ivisc=0, so that, upon restart, CFL3D will compute min distance
c     function for turbulent cases.
c
c     If using this alone (not in conjunction with cfl3d makefile):
c     f90 -64 -r8 initialize_field.f umalloc_r.o -o initialize_field
c*****************************************************************************
c
      parameter (nblock=10000)
c
      dimension rms(1),clw(1),cdw(1),cdpw(1),
     .          cdvw(1),cxw(1),cyw(1),czw(1),
     .          cmxw(1),cmyw(1),cmzw(1),
     .          fmdotw(1),cftmomw(1),cftpw(1),
     .          cftvw(1),cfttotw(1),
     .          rmstr1(1),rmstr2(1),nneg1(1),
     .          nneg2(1)
      dimension titlw(20)
      dimension id(nblock),jd(nblock),kd(nblock)
      character file1*80,file2*80
      character*1 dum1
c
      gamma=1.4
c
      write(6,'('' what is corresponding input file name?'')')
      read(5,'(a80)') file1
      open(2,file=file1,form='formatted',status='old')
      write(6,'('' input new V6 unformatted restart file name'',
     . '' desired:'')')
      read(5,'(a80)') file2
      open(3,file=file2,form='unformatted',status='unknown')
      write(6,'('' write ghost values (1=yes, 0=no)'',
     +  '' (1=CFL3Ds default):'')')
      read(5,*) iwghost
c
      ncycmax=1
      do n=1,ncycmax
        rms(n)=1.
        clw(n)=1.
        cdw(n)=1.
        cdpw(n)=1.
        cdvw(n)=1.
        cxw(n)=1.
        cyw(n)=1.
        czw(n)=1.
        cmxw(n)=1. 
        cmyw(n)=1.
        cmzw(n)=1.
        fmdotw(n)=1.
        cftmomw(n)=1.
        cftpw(n)=1.
        cftvw(n)=1.
        cfttotw(n)=1.
        rmstr1(n)=1.
        rmstr2(n)=1.
        nneg1(n)=0
        nneg2(n)=0
      enddo
c 
c   Read input file
      do n=1,14
        read(2,*)
      enddo
      read(2,'(a1)') dum1
      if (dum1 .eq. '>') then
        do n=1,500
          read(2,'(a1)') dum1
          if (dum1 .eq. '<') goto 1002
        enddo
        write(6,'('' Error, too many lines (>500) of keyword input'')')
        stop
 1002   continue
      else
        backspace(2)
      end if
      read(2,*)
      read(2,*)
      read(2,*) xmachw,alphw,betaw,reuew
      reuew=reuew*1.e6
      do n=1,5
        read(2,*)
      enddo
      read(2,*) ngrid
      nblk=abs(ngrid)
      if (nblk .gt. nblock) then
        write(6,'('' need to increase nblock to '',i5)') nblk
        stop
      end if
      read(2,*)
      do n=1,nblk
        read(2,*)
      enddo
      read(2,*)
      do n=1,nblk
        read(2,*) id(n),jd(n),kd(n)
      enddo
c
      write(6,'('' Are you using mesh sequencing (1=yes)?'')')
      read(5,*) imesh
      if(imesh .eq. 1) then
        write(6,'('' How many levels down are you starting'',
     .   '' (mseq-1)?'')')
        read(5,*) idown
        do n=1,nblk
          do jj=1,idown
            if(id(n) .ne. 2) then
              id(n)=id(n)/2+1        
            end if
            jd(n)=jd(n)/2+1        
            kd(n)=kd(n)/2+1        
          enddo
        enddo
      end if
c
      write(6,'('' input 0=input once for all zones'')')
      write(6,'('' input 1=input individually for each zone'')')
      read(5,*) izon
c
      iskip = 1
      ntr=1
      time=0.
      do n=1,20
        titlw(n)=0.
      enddo
c
c     Do over all the blocks
      do 9897 nrty=1,nblk
      it=id(nrty)
      jt=jd(nrty)
      kt=kd(nrty)
      if (izon .eq. 1) then
        write(6,'('' for block '',i5,'', idim,jdim,kdim='',3i5)') nrty,
     +    it,jt,kt
      end if
      write(3) titlw,xmachw,jt,kt,it,alphw,reuew,ntr,time
c
c     Convergence data (residual,force coefficients, mass flow, etc.)
c
      if (iskip.gt.0) then
        write(3) (rms(n),     n=1,ntr),(clw(n),     n=1,ntr),
     .           (cdw(n),     n=1,ntr),(cdpw(n),    n=1,ntr),
     .           (cdvw(n),    n=1,ntr),(cxw(n),     n=1,ntr),
     .           (cyw(n),     n=1,ntr),(czw(n),     n=1,ntr),
     .           (cmxw(n),    n=1,ntr),(cmyw(n),    n=1,ntr),
     .           (cmzw(n),    n=1,ntr),(fmdotw(n),  n=1,ntr),
     .           (cftmomw(n), n=1,ntr),(cftpw(n),   n=1,ntr),
     .           (cftvw(n),   n=1,ntr),(cfttotw(n), n=1,ntr)
      end if
c
c     Primative variables (rho,u,v,w,p)
c
      jdim1=jt-1
      kdim1=kt-1
      idim1=it-1
c
      call writeq(jdim1,kdim1,idim1,xmachw,alphw,betaw,
     +  gamma,nrty,iwghost,izon,rho0,u0,v0,w0,p0)
c
c     Turbulence quantities
c
      iv1=0
      iv2=0
      iv3=0
      write(3) iv1,iv2,iv3
c
c     Turbulence model convergence data
c
      if (iskip.gt.0) then
         write(3) (rmstr1(n),n=1,ntr),(rmstr2(n),n=1,ntr),
     .            (nneg1(n), n=1,ntr),(nneg2(n), n=1,ntr)
      end if
c
c
      iskip=0
      if (izon .eq. 0) then
        write(6,'('' for block '',i5,'', idim,jdim,kdim='',3i5)') nrty,
     +    it,jt,kt
      end if
 9897 continue
c
      write(6,'('' New V6 restart file successfully written'',
     . '' for '',i5,'' blocks'',/)') nblk
      if (iwghost .eq. 1) then
        write(6,'('' ghost values defaulted to local interior'',
     .   '' values'')')
      end if
c
      stop
      end
c
c   *************************************************************
      subroutine writeq(jdim1,kdim1,idim1,xmachw,alphw,betaw,
     +  gamma,nrty,iwghost,izon,rho0,u0,v0,w0,p0)
c
      integer stats
c
      allocatable :: q(:,:,:,:)
      allocatable :: qi0(:,:,:,:)
      allocatable :: qj0(:,:,:,:)
      allocatable :: qk0(:,:,:,:)
c
c     allocate memory
c
      memuse = 0
      allocate( q(jdim1,kdim1,idim1,5), stat=stats )
      call umalloc_r(jdim1*kdim1*idim1*5,0,'q',memuse,stats)
      allocate( qi0(jdim1,kdim1,5,4), stat=stats )
      call umalloc_r(jdim1*kdim1*5*4,0,'qi0',memuse,stats)
      allocate( qj0(kdim1,idim1,5,4), stat=stats )
      call umalloc_r(kdim1*idim1*5*4,0,'qj0',memuse,stats)
      allocate( qk0(jdim1,idim1,5,4), stat=stats )
      call umalloc_r(jdim1*idim1*5*4,0,'qk0',memuse,stats)
c
      if (izon .eq. 0) then
c   do once only for 1st zone, all other zones the same
c   range is automatically everywhere in this case
        if (nrty .eq. 1) then
          write(6,'('' input 0 for freestream, '',
     +     ''1 for other:'')')
          read(5,*) iset
          if (iset .eq. 0) then
            rho0=1.
            u0     = xmachw*cos(alphw)*cos(betaw)
            w0     = xmachw*sin(alphw)*cos(betaw)
            v0     = -xmachw*sin(betaw)
            p0     = 1./gamma
          else
            write(6,'('' input desired rho,u,v,w,p:'')')
            read(5,*) rho0,u0,v0,w0,p0
          end if
        else
          continue
        end if
        do i=1,idim1
        do j=1,jdim1
        do k=1,kdim1
          q(j,k,i,1)=rho0
          q(j,k,i,2)=u0
          q(j,k,i,3)=v0
          q(j,k,i,4)=w0
          q(j,k,i,5)=p0
        enddo
        enddo
        enddo
      else
c
c   do individually for each zone
      write(6,'('' For block '',i5,'' input 0 for freestream, '',
     + ''1 for other:'')') nrty
      read(5,*) iset
c   Default = freestream
        rho0=1.
        u0     = xmachw*cos(alphw)*cos(betaw)
        w0     = xmachw*sin(alphw)*cos(betaw)
        v0     = -xmachw*sin(betaw)
        p0     = 1./gamma
        do i=1,idim1
        do j=1,jdim1
        do k=1,kdim1
          q(j,k,i,1)=rho0
          q(j,k,i,2)=u0
          q(j,k,i,3)=v0
          q(j,k,i,4)=w0
          q(j,k,i,5)=p0
        enddo
        enddo
        enddo
      if(iset .eq. 1) then
 1055   continue
        write(6,'('' input desired rho,u,v,w,p:'')')
        read(5,*) rho0,u0,v0,w0,p0
        write(6,'('' input index range ilo,ihi,jlo,jhi,klo,khi'',
     +    ''(0 = min or max):'')')
        write(6,'(''    (input cell center values... e.g., max ihi,'',
     +    ''jhi, khi should be'')')
        write(6,'(''    one less than idim, jdim, kdim)'')')
        read(5,*) ilo,ihi,jlo,jhi,klo,khi
        if(ilo .le. 0) ilo=1
        if(jlo .le. 0) jlo=1
        if(klo .le. 0) klo=1
        if(ihi .le. 0) ihi=idim1
        if(ihi .gt. idim1) then
          write(6,'('' ihi cannnot excede idim1, set to '',i5)') idim1
          ihi=idim1
        end if
        if(jhi .le. 0) jhi=jdim1
        if(jhi .gt. jdim1) then
          write(6,'('' jhi cannnot excede jdim1, set to '',i5)') jdim1
          jhi=jdim1
        end if
        if(khi .le. 0) khi=kdim1
        if(khi .gt. kdim1) then
          write(6,'('' khi cannnot excede kdim1, set to '',i5)') kdim1
          khi=kdim1
        end if
        do i=ilo,ihi
        do j=jlo,jhi
        do k=klo,khi
          q(j,k,i,1)=rho0
          q(j,k,i,2)=u0
          q(j,k,i,3)=v0
          q(j,k,i,4)=w0
          q(j,k,i,5)=p0
        enddo
        enddo
        enddo
        write(6,'('' More index ranges? (1=yes) (if you dont specify'',
     +   '' it will be set to freestream)'')')
        read(5,*) iset2
        if(iset2 .eq. 1) goto 1055
      end if
      end if
c
      write(3) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
c
      if (iwghost .eq. 1) then
c   Default for ghost values = interior value next to it
        do m=1,2
        do l=1,5
        do j=1,jdim1
        do k=1,kdim1
          qi0(j,k,l,m)=q(j,k,1,l)
        enddo
        enddo
        do i=1,idim1
        do k=1,kdim1
          qj0(k,i,l,m)=q(1,k,i,l)
        enddo
        enddo
        do i=1,idim1
        do j=1,jdim1
          qk0(j,i,l,m)=q(j,1,i,l)
        enddo
        enddo
        enddo
        enddo
c
        do m=3,4
        do l=1,5
        do j=1,jdim1
        do k=1,kdim1
          qi0(j,k,l,m)=q(j,k,idim1,l)
        enddo
        enddo
        do i=1,idim1
        do k=1,kdim1
          qj0(k,i,l,m)=q(jdim1,k,i,l)
        enddo
        enddo
        do i=1,idim1
        do j=1,jdim1
          qk0(j,i,l,m)=q(j,kdim1,i,l)
        enddo
        enddo
        enddo
        enddo
        write(3) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),l=1,5),m=1,4),
     .           ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),m=1,4),
     .           ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),m=1,4)
      end if
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(q)
         deallocate(qi0)
         deallocate(qj0)
         deallocate(qk0)
      end if
c
      return
      end
