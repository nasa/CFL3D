      program Get_FD
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
      dimension rms1(ncycarg),clw1(ncycarg),
     .          cdw1(ncycarg),cdpw1(ncycarg),
     .          cdvw1(ncycarg),cxw1(ncycarg),
     .          cyw1(ncycarg),czw1(ncycarg),
     .          cmxw1(ncycarg),cmyw1(ncycarg),
     .          cmzw1(ncycarg),fmdotw1(ncycarg),
     .          cftmomw1(ncycarg),cftpw1(ncycarg),
     .          cftvw1(ncycarg),cfttotw1(ncycarg)
      dimension rms2(ncycarg),clw2(ncycarg),
     .          cdw2(ncycarg),cdpw2(ncycarg),
     .          cdvw2(ncycarg),cxw2(ncycarg),
     .          cyw2(ncycarg),czw2(ncycarg),
     .          cmxw2(ncycarg),cmyw2(ncycarg),
     .          cmzw2(ncycarg),fmdotw2(ncycarg),
     .          cftmomw2(ncycarg),cftpw2(ncycarg),
     .          cftvw2(ncycarg),cfttotw2(ncycarg)
      dimension titlw(20)
c
      character*80 restfile1,restfile2,fdfile
c
      common /files/ restfile1,restfile2,fdfile
c
      pointer (ip_rms1,rms1),(ip_clw1,clw1),(ip_cdw1,cdw1),
     .        (ip_cdpw1,cdpw1),(ip_cdvw1,cdvw1),(ip_cxw1,cxw1),
     .        (ip_cyw1,cyw1),(ip_czw1,czw1),(ip_cmxw1,cmxw1),
     .        (ip_cmyw1,cmyw1),(ip_cmzw1,cmzw1),(ip_fmdotw1,fmdotw1),
     .        (ip_cftmomw1,cftmomw1),(ip_cftpw1,cftpw1),
     .        (ip_cftvw1,cftvw1),(ip_cfttotw1,cfttotw1)
      pointer (ip_rms2,rms2),(ip_clw2,clw2),(ip_cdw2,cdw2),
     .        (ip_cdpw2,cdpw2),(ip_cdvw2,cdvw2),(ip_cxw2,cxw2),
     .        (ip_cyw2,cyw2),(ip_czw2,czw2),(ip_cmxw2,cmxw2),
     .        (ip_cmyw2,cmyw2),(ip_cmzw2,cmzw2),(ip_fmdotw2,fmdotw2),
     .        (ip_cftmomw2,cftmomw2),(ip_cftpw2,cftpw2),
     .        (ip_cftvw2,cftvw2),(ip_cfttotw2,cfttotw2)
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
      icall  = 0
      memuse = 0
      call umalloc_r(ip_rms1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_clw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdpw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdvw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cxw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cyw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_czw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmxw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmyw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmzw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_fmdotw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftmomw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftpw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftvw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cfttotw1,ncycarg,0,icall,memuse)
      call umalloc_r(ip_rms2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_clw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdpw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cdvw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cxw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cyw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_czw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmxw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmyw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cmzw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_fmdotw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftmomw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftpw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cftvw2,ncycarg,0,icall,memuse)
      call umalloc_r(ip_cfttotw2,ncycarg,0,icall,memuse)
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
      write(8,'(''DERIVATIVES VIA CENTRAL DIFFERENCES, '',
     .          ''WITH STEP SIZE = '',e12.4)') stepsize
      write(8,'('' '')')
      write(8,'(''dcl/ddv = '',e13.5)') 
     .            (clw1(ntr1)-clw2(ntr2))/stepsize/2.
      write(8,'('' '')')
      write(8,'(''dcd/ddv = '',e13.5)') 
     .            (cdw1(ntr1)-cdw2(ntr2))/stepsize/2.
      write(8,'('' '')')
      write(8,'(''dcy/ddv = '',e13.5)') 
     .            (cyw1(ntr1)-cyw2(ntr2))/stepsize/2.
      write(8,'('' '')')
      write(8,'(''dcmy/ddv = '',e13.5)') 
     .            (cmyw1(ntr1)-cmyw2(ntr2))/stepsize/2.
      write(8,'('' '')')
      write(8,'(''dcz/ddv = '',e13.5)')
     .            (czw1(ntr1)-czw2(ntr2))/stepsize/2.
      write(8,'('' '')')
      write(8,'(''dcmz/ddv = '',e13.5)')
     .            (cmzw1(ntr1)-cmzw2(ntr2))/stepsize/2.
      write(8,'('' '')')
c
c     output convergence history
c
      write(8,'('' '')')
      write(8,'(''derivative convergence history'')')
      write(8,'('' '')')
      if (ialph.eq.0) then
         write(8,'(''  iter     residual      dcl/ddv      dcd/ddv'',
     .             ''      dcy/ddv     dcmy/ddv'')')
      else
         write(8,'(''  iter     residual      dcl/ddv      dcd/ddv'',
     .             ''      dcz/ddv     dcmz/ddv'')')
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
      else
         do n=1,ntr1
            write(8,'(i6,5e13.5)') n,1.,
     .                             (clw1(n)-clw2(n))/stepsize/2.,
     .                             (cdw1(n)-cdw2(n))/stepsize/2.,
     .                             (czw1(n)-czw2(n))/stepsize/2.,
     .                             (cmzw1(n)-cmzw2(n))/stepsize/2.
        end do
      end if
c
c     free memory
c
      call ufree_r(ip_rms1)
      call ufree_r(ip_clw1)
      call ufree_r(ip_cdw1)
      call ufree_r(ip_cdpw1)
      call ufree_r(ip_cdvw1)
      call ufree_r(ip_cxw1)
      call ufree_r(ip_cyw1)
      call ufree_r(ip_czw1)
      call ufree_r(ip_cmxw1)
      call ufree_r(ip_cmyw1)
      call ufree_r(ip_cmzw1)
      call ufree_r(ip_fmdotw1)
      call ufree_r(ip_cftmomw1)
      call ufree_r(ip_cftpw1)
      call ufree_r(ip_cftvw1)
      call ufree_r(ip_cfttotw1)
      call ufree_r(ip_rms2)
      call ufree_r(ip_clw2)
      call ufree_r(ip_cdw2)
      call ufree_r(ip_cdpw2)
      call ufree_r(ip_cdvw2)
      call ufree_r(ip_cxw2)
      call ufree_r(ip_cyw2)
      call ufree_r(ip_czw2)
      call ufree_r(ip_cmxw2)
      call ufree_r(ip_cmyw2)
      call ufree_r(ip_cmzw2)
      call ufree_r(ip_fmdotw2)
      call ufree_r(ip_cftmomw2)
      call ufree_r(ip_cftpw2)
      call ufree_r(ip_cftvw2)
      call ufree_r(ip_cfttotw2)
c
      return
      end
