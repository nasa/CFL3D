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
      program p3d_to_cfl3drst
c
c     $Id$
c
c***********************************************************************
c     Purpose:
c     Creates a CFL3D V6 restart file from a plot3d Q-file.
c     Must be from a grid that is the same size as the desired restart file,
c     and the plot3d q-file must be at GRIDPOINTS.
c     This tool does not allow setting of initial turbulence quantities.
c     It also sets ivisc=0, so that upon restart, CFL3D computes min distance
c     function for turbulent cases.
c     If it is converting a 2-D plot3d q-file, it assumes that the results
c     are in the x-z plane!!!
c
c     If using this alone (not in conjunction with cfl3d makefile):
c     f90 -64 -r8 p3d_to_cfl3drst.f umalloc_r.o -o p3d_to_cfl3drst
c***********************************************************************
c
      parameter (nblock=10000)
c
      dimension rms(1),clw(1),cdw(1),cdpw(1),
     .          cdvw(1),cxw(1),cyw(1),czw(1),
     .          cmxw(1),cmyw(1),cmzw(1),
     .          fmdotw(1),cftmomw(1),cftpw(1),
     .          cftvw(1),cfttotw(1),
     .          rmstr1(1),rmstr2(1)
      dimension titlw2(20)
c
      dimension nneg1(1),nneg2(1)
      dimension id(nblock),jd(nblock),kd(nblock)
      character file1*80,file2*80
c
      gamma=1.4
c
      write(6,'('' what is plot3d q-file name?'')')
      read(5,'(a80)') file1
      write(6,'(''   note: q values must be STANDARD in nondim:'')')
      write(6,'(''   i.e.: rho/rhoinf, rho*u/(rhoinf*ainf),'',
     + '' rho*E/(rhoinf*ainf**2)'',/)')
      write(6,'('' input 1 to read single-precision q-file'')')
      write(6,'('' input 0 to read double-precision q-file:'')')
      read(5,*) isingle
#if defined ASN_P3D
      write(6,'('' input 1 if need to call asn (on cray)'',
     + '', 0 otherwise:'')')
      read(5,*) iasn
      if (iasn .eq. 1) then
        call asnfile(file1, '-F f77 -N ieee',IER)
      end if
#endif
      write(6,'('' input 0=3-D q-file, 1=2-D q-file'')')
      read(5,*) i2d
      open(13,file=file1,form='unformatted',status='old')
      write(6,'('' input new V6 unformatted restart file name:'')')
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
c   Read q-file
c
      read(13) nblk
      write(6,'('' number of zones='',i5)') nblk
      if(nblk .gt. nblock) then
        write(6,'('' need to increase nblock to '',i5)') nblk
        stop
      end if
      if (i2d .eq. 0) then
        read(13) (id(n),jd(n),kd(n),n=1,nblk)
      else
        read(13) (jd(n),kd(n),n=1,nblk)
        do n=1,nblk
          id(n)=2
        enddo
      end if
c
      imax  = 1
      jmax  = 1
      kmax  = 1
      do n=1,nblk
         write(6,'('' read in zone '',i5,'', i,j,k='',3i6)')
     +     n,id(n),jd(n),kd(n)
         imax  = max(imax,id(n))
         jmax  = max(jmax,jd(n))
         kmax  = max(kmax,kd(n))
      end do
c
      call writeq(imax,jmax,kmax,nblk,isingle,i2d,
     + id,jd,kd,nblock,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     + cmxw,cmyw,cmzw,fmdotw,cftmomw,cftpw,cftvw,
     + cfttotw,rmstr1,rmstr2,nneg1,nneg2,iwghost)
c
      write(6,'('' successful completion'')')
      stop
      end
c
c *****************************************************************
      subroutine writeq(imax,jmax,kmax,nblk,isingle,i2d,
     + id,jd,kd,nblock,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     + cmxw,cmyw,cmzw,fmdotw,cftmomw,cftpw,cftvw,
     + cfttotw,rmstr1,rmstr2,nneg1,nneg2,iwghost)
c
      integer stats
c
      real*4 xmachw,alphw,reuew,time
c
      dimension cdpw(1)
      dimension cdvw(1)
      dimension cdw(1)
      dimension cftmomw(1)
      dimension cftpw(1)
      dimension cfttotw(1)
      dimension cftvw(1)
      dimension clw(1)
      dimension cmxw(1)
      dimension cmyw(1)
      dimension cmzw(1)
      dimension cxw(1)
      dimension cyw(1)
      dimension czw(1)
      dimension fmdotw(1)
      dimension id(nblock)
      dimension idsav(nblock)
      dimension jd(nblock)
      dimension jdsav(nblock)
      dimension kd(nblock)
      dimension kdsav(nblock)
      dimension nneg1(1)
      dimension nneg2(1)
      dimension rms(1)
      dimension rmstr1(1)
      dimension rmstr2(1)
      dimension titlw2(20)
c
      real(4), allocatable :: q(:,:,:,:,:)
      allocatable :: q2(:,:,:,:,:)
      allocatable :: qi0(:,:,:,:)
      allocatable :: qj0(:,:,:,:)
      allocatable :: qk0(:,:,:,:)
c
c     allocate memory
c
      memuse = 0
      allocate( q(imax,jmax,kmax,5,nblk), stat=stats )
      call umalloc_r(imax*jmax*kmax*5*nblk,0,'q',memuse,stats)
      allocate( q2(imax,jmax,kmax,5,nblk), stat=stats )
      call umalloc_r(imax*jmax*kmax*5*nblk,0,'q2',memuse,stats)
      allocate( qi0(jmax,kmax,5,4), stat=stats )
      call umalloc_r(jmax*kmax*5*4,0,'qi0',memuse,stats)
      allocate( qj0(kmax,imax,5,4), stat=stats )
      call umalloc_r(kmax*imax*5*4,0,'qj0',memuse,stats)
      allocate( qk0(jmax,imax,5,4), stat=stats )
      call umalloc_r(jmax*imax*5*4,0,'qk0',memuse,stats)
c
      do 60 n=1,nblk
        if (isingle .eq. 1) then
          read(13) xmachw,alphw,reuew,time
          xmachw2=xmachw
          alphw2=alphw
          reuew2=reuew
          time2=time
        else
          read(13) xmachw2,alphw2,reuew2,time2
        end if
        write(6,'('' n='',i5,'': M,id,jd,kd='',f12.5,3i5)')
     +   n,xmachw2,id(n),jd(n),kd(n)
        if(i2d .eq. 0) then
          if (isingle .eq. 1) then
            read(13) ((((q(i,j,k,l,n),i=1,id(n)),j=1,jd(n)),
     .      k=1,kd(n)),l=1,5)
            do l=1,5
            do i=1,id(n)
            do j=1,jd(n)
            do k=1,kd(n)
              q2(i,j,k,l,n)=q(i,j,k,l,n)
            enddo
            enddo
            enddo
            enddo
          else
            read(13) ((((q2(i,j,k,l,n),i=1,id(n)),j=1,jd(n)),
     .      k=1,kd(n)),l=1,5)
          end if
        else
          if (isingle .eq. 1) then
            read(13) (((q(1,j,k,l,n),j=1,jd(n)),k=1,kd(n)),l=1,4)
            do l=1,5
            do i=1,id(n)
            do j=1,jd(n)
            do k=1,kd(n)
              q2(i,j,k,l,n)=q(i,j,k,l,n)
            enddo
            enddo
            enddo
            enddo
          else
            read(13) (((q2(1,j,k,l,n),j=1,jd(n)),k=1,kd(n)),l=1,4)
          end if
          do j=1,jd(n)
          do k=1,kd(n)
            q2(1,j,k,5,n)=q2(1,j,k,4,n)
            q2(1,j,k,4,n)=q2(1,j,k,3,n)
            q2(1,j,k,3,n)=0.
            q2(2,j,k,5,n)=q2(1,j,k,5,n)
            q2(2,j,k,4,n)=q2(1,j,k,4,n)
            q2(2,j,k,3,n)=0.
            q2(2,j,k,1,n)=q2(1,j,k,1,n)
            q2(2,j,k,2,n)=q2(1,j,k,2,n)
          enddo
          enddo
        end if
c
c      convert plot3d quantities to primative
c
        gamma=1.4
        do 1269 j=1,jd(n)
        do 1269 k=1,kd(n)
        do 1269 i=1,id(n)
          rho = q2(i,j,k,1,n)
          u   = q2(i,j,k,2,n)/rho
          v   = q2(i,j,k,3,n)/rho
          w   = q2(i,j,k,4,n)/rho
          vsq = u*u + v*v + w*w
          pp   = (gamma-1.)*(q2(i,j,k,5,n)-.5*rho*vsq)
          csq = gamma*pp/rho
          dmsq= vsq/csq
          q2(i,j,k,1,n) = rho
          q2(i,j,k,2,n) = u
          q2(i,j,k,3,n) = v
          q2(i,j,k,4,n) = w
          q2(i,j,k,5,n) = pp
1269    continue
        idsav(n)=id(n)
        jdsav(n)=jd(n)
        kdsav(n)=kd(n)
60    continue
      write(6,'('' Are you using mesh sequencing (1=yes)?'')')
      read(5,*) imesh
      if(imesh .eq. 1) then
        write(6,'('' How many levels down are you starting?'')')
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
      iskip = 1
      ntr=1
      do n=1,20
        titlw2(n)=0.
      enddo
c
      do 9897 nrty=1,nblk
      it=id(nrty)
      jt=jd(nrty)
      kt=kd(nrty)
      write(6,'('' for block '',i5,'', idim,jdim,kdim='',3i5)') nrty,
     +  it,jt,kt
      write(3) titlw2,xmachw2,jt,kt,it,alphw2,reuew2,ntr,time2
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
      jdim1=jdsav(nrty)-1
      kdim1=kdsav(nrty)-1
      idim1=idsav(nrty)-1
c
c
      if (iwghost .eq. 1) then
c   Default for ghost values = interior value next to it
        do m=1,2
        do l=1,5
        i=1
        do j=1,jdim1
        do k=1,kdim1
      qi0(j,k,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        j=1
        do i=1,idim1
        do k=1,kdim1
      qj0(k,i,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        k=1
        do i=1,idim1
        do j=1,jdim1
      qk0(j,i,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        enddo
        enddo
c
        do m=3,4
        do l=1,5
        i=idim1
        do j=1,jdim1
        do k=1,kdim1
      qi0(j,k,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        j=jdim1
        do i=1,idim1
        do k=1,kdim1
      qj0(k,i,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        k=kdim1
        do i=1,idim1
        do j=1,jdim1
      qk0(j,i,l,m)=0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
        enddo
        enddo
        enddo
        enddo
      end if
c
      if (imesh .eq. 1) then
      is=idown+1
      write(3) ((((0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
     + ,j=1,jdim1,is),k=1,kdim1,is),i=1,idim1,is),l=1,5)
      if (iwghost .eq. 1) then
        write(3) ((((qi0(j,k,l,m),j=1,jdim1,is),k=1,kdim1,is),l=1,5),
     +   m=1,4),
     +           ((((qj0(k,i,l,m),k=1,kdim1,is),i=1,idim1,is),l=1,5),
     +   m=1,4),
     +           ((((qk0(j,i,l,m),j=1,jdim1,is),i=1,idim1,is),l=1,5),
     +   m=1,4)
      end if
      else
      write(3) ((((0.125*(q2(i  ,j  ,k  ,l,nrty)+q2(i+1,j  ,k  ,l,nrty)+
     +                    q2(i  ,j+1,k  ,l,nrty)+q2(i  ,j  ,k+1,l,nrty)+
     +                    q2(i+1,j+1,k  ,l,nrty)+q2(i+1,j  ,k+1,l,nrty)+
     +                    q2(i  ,j+1,k+1,l,nrty)+q2(i+1,j+1,k+1,l,nrty))
     + ,j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
      if (iwghost .eq. 1) then
        write(3) ((((qi0(j,k,l,m),j=1,jdim1),k=1,kdim1),l=1,5),
     +   m=1,4),
     +           ((((qj0(k,i,l,m),k=1,kdim1),i=1,idim1),l=1,5),
     +   m=1,4),
     +           ((((qk0(j,i,l,m),j=1,jdim1),i=1,idim1),l=1,5),
     +   m=1,4)
      end if
      end if
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
 9897 continue
c
      write(6,'('' New V6 restart file successfully written'',
     . '' for '',i5,'' blocks'',/)') nblk
      if (iwghost .eq. 1) then
        write(6,'('' ghost values defaulted to local interior'',
     .   '' values'')')
      end if
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(q)
         deallocate(q2)
         deallocate(qi0)
         deallocate(qj0)
         deallocate(qk0)
      end if
c
      stop
      end
