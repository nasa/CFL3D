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
c     f90 -64 -r8 v6_restart_mod.f umalloc_r.o ufree_r.o -o v6_restart_mod
c     On CRAY:
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 -c umalloc_c.F
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 -c ufree_c.F
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 -c ufree_r.F
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 -c umalloc_r.F
c ar ruc libs.a umalloc_c.o ufree_c.o umalloc_r.o ufree_r.o
c f90 -DASN_P3D -DCRAY_TIME -DWKSTN_OFF -Oscalar3,vector3 v6_restart_mod.f -o v6_restart_mod libs.a
c*****************************************************************************
c
      dimension titlw(20)
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
          end if
          if (iunfi .eq. 1) then
          read(2,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          else
          read(2,*,end=9898) (((dum,j=1,jdim1),k=1,kdim1),i=1,idim1)
          end if
          if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4) then
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
     . ihalf,nblk,imaxc,jmaxc,kmaxc,ic,jc,kc,titlw)
      write(6,'(/,'' successful program completion'')')
      stop
 9898 write(6,'('' Error, end of file reached... are you sure you'')')
      write(6,'('' input the correct number of zones?'')')
      stop
      end
c
c *********************************************************************
      subroutine readwrite(idim,jdim,kdim,ncycmax,iunfi,iunfo,irghost,
     . ihalf,nblk,idimc,jdimc,kdimc,ic,jc,kc,titlw)
c
      dimension q(jdim,kdim,idim,5)
      dimension qc(jdimc,kdimc,idimc,5)
      dimension work(jdimc,kdimc,idimc,5)
      dimension qi0(jdim+1,kdim+1,5,4),qj0(kdim+1,idim,5,4),
     .          qk0(jdim+1,idim,5,4)
      dimension qi0c(jdimc+1,kdimc+1,5,4),qj0c(kdimc+1,idimc,5,4),
     .          qk0c(jdimc+1,idimc,5,4)
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),cdpw(ncycmax),
     .          cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          fmdotw(ncycmax),cftmomw(ncycmax),cftpw(ncycmax),
     .          cftvw(ncycmax),cfttotw(ncycmax),
     .          rmstr1(ncycmax),rmstr2(ncycmax)
      dimension nneg1(ncycmax),nneg2(ncycmax)
      dimension titlw(20)
c
      pointer (ip_q,q),(ip_qc,qc),(ip_qi0,qi0),(ip_qj0,qj0),(ip_qk0,qk0)
      pointer (ip_qi0c,qi0c),(ip_qj0c,qj0c),(ip_qk0c,qk0c)
      pointer (ip_rms,rms),(ip_clw,clw),(ip_cdw,cdw),(ip_cdpw,cdpw),
     .        (ip_cdvw,cdvw),(ip_cxw,cxw),(ip_cyw,cyw),(ip_czw,czw),
     .        (ip_cmxw,cmxw),(ip_cmyw,cmyw),(ip_cmzw,cmzw),
     .        (ip_fmdotw,fmdotw),(ip_cftmomw,cftmomw),(ip_cftpw,cftpw),
     .        (ip_cftvw,cftvw),(ip_cfttotw,cfttotw),
     .        (ip_rmstr1,rmstr1),(ip_rmstr2,rmstr2)
      pointer (ip_work,work)
      pointer (ip_nneg1,nneg1),(ip_nneg2,nneg2)
c
c     allocate memory
c
      icall  = 0
      memuse = 0
      call umalloc_r(ip_q,jdim*kdim*idim*5,0,icall,memuse)
      call umalloc_r(ip_qc,jdimc*kdimc*idimc*5,0,icall,memuse)
      call umalloc_r(ip_qi0,(jdim+1)*(kdim+1)*20,0,icall,memuse)
      call umalloc_r(ip_qj0,(kdim+1)*idim*20,0,icall,memuse)
      call umalloc_r(ip_qk0,(jdim+1)*idim*20,0,icall,memuse)
      call umalloc_r(ip_qi0c,(jdimc+1)*(kdimc+1)*20,0,icall,memuse)
      call umalloc_r(ip_qj0c,(kdimc+1)*idimc*20,0,icall,memuse)
      call umalloc_r(ip_qk0c,(jdimc+1)*idimc*20,0,icall,memuse)
      call umalloc_r(ip_rms,ncycmax,0,icall,memuse)
      call umalloc_r(ip_clw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cdw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cdpw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cdvw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cxw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cyw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_czw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cmxw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cmyw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cmzw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_fmdotw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cftmomw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cftpw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cftvw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_cfttotw,ncycmax,0,icall,memuse)
      call umalloc_r(ip_rmstr1,ncycmax,0,icall,memuse)
      call umalloc_r(ip_rmstr2,ncycmax,0,icall,memuse)
      call umalloc_r(ip_work,jdimc*kdimc*idimc*5,0,icall,memuse)
      call umalloc_r(ip_nneg1,ncycmax,1,icall,memuse)
      call umalloc_r(ip_nneg2,ncycmax,1,icall,memuse)
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
      if (iunfi .eq. 1) then
      read(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
      else
      read(2,*) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,5)
      end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,5,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,5)
      else
c       Don't coarsen
        do l=1,5
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
     +  l=1,5)
      else
      write(3,*) ((((qc(j,k,i,l),j=1,jdim1c),k=1,kdim1c),i=1,idim1c),
     +  l=1,5)
      end if
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
      end if
c
c   vist3d data:
      if (iv1.ge.2 .or. iv2.ge.2 .or. iv3.ge.2) then
         if (iunfi .eq. 1) then
         read(2) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         else
         read(2,*) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,1,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,1)
      else
c       Don't coarsen
        do j=1,jdim1
        do k=1,kdim1
        do i=1,idim1
          qc(j,k,i,1)=q(j,k,i,1)
        enddo
        enddo
        enddo
      end if
         if (iunfo .eq. 1) then
         write(3) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
         else
         write(3,*) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),
     .     i=1,idim1c)
         end if
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
         if (iunfi .eq. 1) then
         read(2) ((((q(j,k,i,m),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .           m=1,2)
         else
         read(2,*) ((((q(j,k,i,m),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .           m=1,2)
         end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,2,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,2)
      else
c       Don't coarsen
        do l=1,2
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
         write(3) ((((qc(j,k,i,m),j=1,jdim1c),k=1,kdim1c),
     +           i=1,idim1c),m=1,2)
         else
         write(3,*) ((((qc(j,k,i,m),j=1,jdim1c),k=1,kdim1c),
     +           i=1,idim1c),m=1,2)
         end if
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
      end if
c
c   smin data:
         if (iunfi .eq. 1) then
         read(2) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         else
         read(2,*) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,1,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,1)
      else
c       Don't coarsen
        do j=1,jdim1
        do k=1,kdim1
        do i=1,idim1
          qc(j,k,i,1)=q(j,k,i,1)
        enddo
        enddo
        enddo
      end if
         if (iunfo .eq. 1) then
         write(3) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
         else
         write(3,*) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
         end if
         if (iv1.eq.4 .or. iv2.eq.4 .or. iv3.eq.4) then
c   xjb,xkb,blnum data
         if (iunfi .eq. 1) then
           read(2) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
           read(2) (((q(j,k,i,2),j=1,jdim1),k=1,kdim1),i=1,idim1)
           read(2) (((q(j,k,i,3),j=1,jdim1),k=1,kdim1),i=1,idim1)
         else
           read(2,*) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
           read(2,*) (((q(j,k,i,2),j=1,jdim1),k=1,kdim1),i=1,idim1)
           read(2,*) (((q(j,k,i,3),j=1,jdim1),k=1,kdim1),i=1,idim1)
         end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,3,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,3)
      else
c       Don't coarsen
        do l=1,3
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
           write(3) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           write(3) (((qc(j,k,i,2),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           write(3) (((qc(j,k,i,3),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           else
           write(3,*) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           write(3,*) (((qc(j,k,i,2),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           write(3,*) (((qc(j,k,i,3),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           end if
         end if
c   cmuv data
         if (iv1.eq.8 .or. iv2.eq.8 .or. iv3.eq.8  .or.
     .       iv1.eq.9 .or. iv2.eq.9 .or. iv3.eq.9  .or.
     .       iv1.eq.13.or. iv2.eq.13.or. iv3.eq.13 .or.
     .       iv1.eq.14.or. iv2.eq.14.or. iv3.eq.14) then
         if (iunfi .eq. 1) then
           read(2) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         else
           read(2,*) (((q(j,k,i,1),j=1,jdim1),k=1,kdim1),i=1,idim1)
         end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,1,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,1)
      else
c       Don't coarsen
        do j=1,jdim1
        do k=1,kdim1
        do i=1,idim1
          qc(j,k,i,1)=q(j,k,i,1)
        enddo
        enddo
        enddo
      end if
           if (iunfo .eq. 1) then
           write(3) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),i=1,idim1c)
           else
           write(3,*) (((qc(j,k,i,1),j=1,jdim1c),k=1,kdim1c),
     .       i=1,idim1c)
           end if
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
           if (iunfi .eq. 1) then
           read(2) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .             l=1,5)
           else
           read(2,*) ((((q(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .             l=1,5)
           end if
      if (ihalf .eq. 1) then
c       Make finer:
        call makefin(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,idim1c,jdim1c,kdim1c,ic,jc,kc,5,
     +    work)
      else if (ihalf .eq. 0) then
c       Coarsen
        call makecoar(q,qc,idim,jdim,kdim,idimc,jdimc,kdimc,
     +    idim1,jdim1,kdim1,ic,jc,kc,it,5)
      else
c       Don't coarsen
        do l=1,5
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
           write(3) ((((qc(j,k,i,l),j=1,jdim1c),k=1,kdim1c),
     +      i=1,idim1c),l=1,5)
           else
           write(3,*) ((((qc(j,k,i,l),j=1,jdim1c),k=1,kdim1c),
     +      i=1,idim1c),l=1,5)
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
         call ufree_r(ip_q)
         call ufree_r(ip_qc)
         call ufree_r(ip_qi0)
         call ufree_r(ip_qj0)
         call ufree_r(ip_qk0)
         call ufree_r(ip_qi0c)
         call ufree_r(ip_qj0c)
         call ufree_r(ip_qk0c)
         call ufree_r(ip_rms)
         call ufree_r(ip_clw)
         call ufree_r(ip_cdw)
         call ufree_r(ip_cdpw)
         call ufree_r(ip_cdvw)
         call ufree_r(ip_cxw)
         call ufree_r(ip_cyw)
         call ufree_r(ip_czw)
         call ufree_r(ip_cmxw)
         call ufree_r(ip_cmyw)
         call ufree_r(ip_cmzw)
         call ufree_r(ip_fmdotw)
         call ufree_r(ip_cftmomw)
         call ufree_r(ip_cftpw)
         call ufree_r(ip_cftvw)
         call ufree_r(ip_cfttotw)
         call ufree_r(ip_rmstr1)
         call ufree_r(ip_rmstr2)
         call ufree_r(ip_work)
         call ufree_r(ip_nneg1)
         call ufree_r(ip_nneg2)
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
