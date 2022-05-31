      parameter (jmax=297,kmax=297,nblmx=3,nexpcp=140)
c
      dimension x(jmax,kmax),y(jmax,kmax),q(jmax,kmax,4)
      dimension jdim(nblmx),kdim(nblmx)
      dimension cp(jmax),cpexp(nexpcp),xexp(nexpcp)
c
      open (unit=3,file='plot3dg.bin',
     .form='unformatted',status='old')
      open (unit=4,file='plot3dq.bin',
     .form='unformatted',status='old')
      open (unit=7,file='2822.cpexp',form='formatted',
     .status='old')
      open (unit=10,file='cp_vs_x.dat',
     .form='formatted',status='unknown')
c
      read(3) ngrid
      read(3) (jdim(n),kdim(n),n=1,ngrid) 
      jd = jdim(1)
      kd = kdim(1)
c
c     check dimensions...note only block 1 actually stored
c
      if (jd.gt.jmax .or. kd.gt.kmax) then
         write(6,*)'stopping...jmax,kmax = ',jmax,kmax
         write(6,*)'block 1....jd,kd     = ',jd,kd
         stop
      end if
c
      read(3)  ((x(j,k),j=1,jd),k=1,kd),
     .         ((y(j,k),j=1,jd),k=1,kd),
     .         ((iblk,j=1,jd),k=1,kd)
      read(4) ngrid
      read(4) (jdim(n),kdim(n),n=1,ngrid)
      read(4) fsmach,alpha,reue,time
      read(4) (((q(j,k,l),j=1,jd),k=1,kd),l=1,4)
c
      rhoe  = 1.
      tinf  = 460.
      c2b   = 198.6/tinf
      c2bp  = c2b+1.
      gamma = 1.4
      gamm1 = gamma-1
c
c     convert to primitive variables
c
      do j=1,jd
         do k=1,kd
           rho = q(j,k,1)
           uu  = q(j,k,2)/rho
           vv  = q(j,k,3)/rho
           pp  = gamma*(gamma-1.)*(q(j,k,4)- 0.5*rho*(uu**2+vv**2))
           q(j,k,2) = uu
           q(j,k,3) = vv
           q(j,k,4) = pp
         end do
      end do
c
c     generate Cp data
c
      do j=1,jd
         pp    = q(j,1,4)
         cp(j) = (pp-1.)/gamma/(.5*fsmach**2)
      end do
c
      jte1 = 41
      jte2 = 217
c
c     output -Cp
c
      write(10,'(''zone'')')
      write(10,'(''# x/c     Cp'')')
      do j=jte1,jte2
         write(10,*) x(j,1),-cp(j)
      end do
c
      read(7,*)
      read(7,*)
      read(7,*)
      read(7,*) ncpdat
      if (ncpdat.gt.nexpcp) then
         write(6,*)'stopping...must increase parameter nexpcp to ',
     .   ncpdat
         stop
      end if
      write(10,'(''zone'')')
      write(10,'(''# EXP data'')')
      write(10,'(''# x/c     Cp'')')
      do j=1,ncpdat
         read(7,*) xexp(j),cpexp(j)
         write(10,*) xexp(j),-cpexp(j)
      end do
c
      stop
      end 
