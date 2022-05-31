      parameter (jmax=297,kmax=297,nblmx=3,nexpcp=40,nexpvel=10)
c
      dimension x(jmax,kmax),y(jmax,kmax),q(jmax,kmax,4)
      dimension jdim(nblmx),kdim(nblmx)
      dimension cp(jmax),cpexp(nexpcp),xexp(nexpcp),etaexp(nexpvel)
      dimension u1(kmax),y1(kmax),v1(kmax),t1(kmax)
      dimension ymaxexp(nexpvel)
c
      write(6,*)'enter 0 to process standard-grid solution'
      write(6,*)'enter 1 to process chimera-grid solution'
      read(5,*) ixmera
      if (ixmera.gt.0) then
         open (unit=3,file='plot3dg.bin_xmera',
     .   form='unformatted',status='old')
         open (unit=4,file='plot3dq.bin_xmera',
     .   form='unformatted',status='old')
      else
         open (unit=3,file='plot3dg.bin_standard',
     .   form='unformatted',status='old')
         open (unit=4,file='plot3dq.bin_standard',
     .   form='unformatted',status='old')
      end if
      open (unit=7,file='4412.cpexp',form='formatted',
     .status='old')
      open (unit=8,file='4412.velexp',form='formatted',
     .status='old')
c
c     stations for velocity profiles
c
      nvel = 7
      if (nvel.gt.nexpvel) then
         write(6,*)'stopping...must increase parameter nexpvel to ',nvel
         stop
      end if
      etaexp(1) = 0.620
      etaexp(2) = 0.675
      etaexp(3) = 0.731
      etaexp(4) = 0.786
      etaexp(5) = 0.842
      etaexp(6) = 0.897
      etaexp(7) = 0.953
c
      if (ixmera.gt.0) then
         open (unit=10,file='cp_vs_x_xmera.dat',
     .   form='formatted',status='unknown')
         open (unit=11,file='u_vs_y_xmera.dat',
     .   form='formatted',status='unknown')
      else
         open (unit=10,file='cp_vs_x_standard.dat',
     .   form='formatted',status='unknown')
         open (unit=11,file='u_vs_y_standard.dat',
     .   form='formatted',status='unknown')
      end if
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
      if (ixmera.gt.0) then
         jte1 = 25
         jte2 = 201
      else
         jte1 = 41
         jte2 = 217
      end if
c
c     output -Cp
c
      write(10,'(''zone'')')
      if (ixmera.gt.0) then
         write(10,'(''# CFD data - chimera grid'')')
      else
         write(10,'(''# CFD data - standard grid'')')
      end if
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
c     generate velocity profile data
c
      rhoe=1. 
      tinf=460. 
      c2b=198.6/tinf
      c2bp=c2b+1. 
      gamma=1.4 
c 
c     note: grid, q data read for block 1 only
c     jhalf dictates which half of block 1 data is used
c     to search for the x/c values corresponding to the
c     experimental stations
c          = 0 use data at all j stations
c            1 use "lower half" data, i.e. j=1,jdim/2
c            2 use "upper half" data, i.e. j=jdim/2,jdim
c 
      jhalf = 2
c
      jd = jdim(1)
      kd = kdim(1)
      if (jhalf.eq.1) then
         js = 1
         je = jd/2
      else if (jhalf.eq.2) then
         js = jd/2
         je = jd
      else
         js = 1
         je = jd
      end if
c
c     read/write experimental profiles. The exp velocity data is 
c     multiplied by 0.925 to convert it to u/uinf
c
      nskip = 14
      do n=1,nskip
         read(8,*)
      end do
      do n=1,nvel
      write(11,'(''zone'')')
      write(11,'(''#Experimental Profiles, x/c = '',f6.4)') etaexp(n)
      write(11,'(''#U/Uref  y/n'')')
         ymaxexp(n) = 0.
         read(8,*) kdpts
         do nn=1,kdpts
            read(8,*) velexp,yexp
            write(11,'(f10.6,2x,f10.6)') velexp*0.925,yexp
            if (yexp.gt.ymaxexp(n)) ymaxexp(n) = yexp
         end do
      end do
c 
c     Computational profiles
c
      do 300 n=1,nvel
c
      write(11,'(''zone'')')
      write(11,'(''#Computational Profiles, x/c = '',f6.4)') etaexp(n)
      write(11,'(''#U/Uref  y/n'')')
c
c     determine nearest grid location from which to extract
c     velocity profiles, and output them
c
      jnd=1
      xmr=abs(x(1,1)-etaexp(n))
      do j=js,je
         if(abs(x(j,1)-etaexp(n)) .lt. xmr) then
           jnd=j
           xmr=abs(x(j,1)-etaexp(n))
         end if
      end do
      write(6,'(''for profile station '',i2,'' j(grid)'',i4,
     .'' x(grid) ='',f8.4,'' xexp ='',f8.4)') n,jnd,x(jnd,1),etaexp(n)
c
      j=jnd 
      ynn=y(j+1,1)-y(j-1,1) 
      xnn=x(j+1,1)-x(j-1,1) 
      ynorm=sqrt(xnn*xnn+ynn*ynn) 
      xn=xnn/ynorm
      yn=ynn/ynorm
c
c     re will equal ue*x/vnu
c
      ccord  = 1.0 
      ue=fsmach
      vnu=ue/reue
      re=sqrt(ue*x(j,1)/vnu) 
      d99=ccord 
      qsqr=q(j,1,4)/q(j,1,1) 
      fmu1=c2bp*(qsqr**1.5)/(c2b+qsqr) 
      do k=1,kd 
         u1(k)=(q(j,k,2)*xn+q(j,k,3)*yn)/ue 
         y1(k)=(sqrt((y(j,k)-y(j,1))**2+(x(j,k)-x(j,1))**2))/d99
         v1(k)=(-q(j,k,2)*yn+q(j,k,3)*xn)*re/ue 
         t1(k)=q(j,k,4)/q(j,k,1)
      end do
      cf1=fmu1*(u1(2)-u1(1))/(0.5*y1(2)*rhoe*reue)
c
c     output computational profiles
c
      do k=1,kd
         write(11,'(f10.6,2x,f10.6)') u1(k),y1(k)
c        if (y1(k).gt.ymaxexp(n)) go to 299
      end do
 299  continue
c
 300  continue
c 
      if (ixmera.eq.0) then
         write(6,*) 'tecplot file for velocity: u_vs_y_standard.dat'
         write(6,*) 'tecplot file for       Cp: cp_vs_x_standard.dat'
      else
         write(6,*) 'tecplot file for velocity: u_vs_y_xmera.dat'
         write(6,*) 'tecplot file for       Cp: cp_vs_x_xmera.dat'
      end if
c 
      stop
      end 
