c
c      gets velocity and pressure data for transonic duct solution 
c      for comparison with experimental data
c
      parameter (id=161,jd=97,kd=2,nbl=10,nvmax=10,nst=20)
c
      dimension x(id,jd,kd,3),q(id,jd,kd,nvmax),
     .idim(nbl),jdim(nbl),kdim(nbl),xprofile(nst),vel(nst,jd),
     .y(nst,jd),ptop(id),pbot(id),vel_exp(nst,jd),y_exp(nst,jd),
     .xbot_exp(id),xtop_exp(id),ptop_exp(id),pbot_exp(id),
     .x_exp(nst),nyval(nst)
c
      character*60 gfile,qfile
      character*80 string
c
      gamma = 1.4
c
c     dimensionalization values
c
      xref   = 0.1444 !ft
      astar  = 1029.8 !ft/sec
      pstar  = 1488.7 !lb/ft**2
      ptotal = 2819.5 !lb/ft**2
c
c     locations for experimental velocity profiles
c
      nsta = 4
      xprofile(1) = 2.882
      xprofile(2) = 4.611
      xprofile(3) = 6.340
      xprofile(4) = 7.493
c
c     write(6,*) 'name of plot3d grid file to read from? '
c     read(5,'(a60)') gfile
c     write(6,*) 'name of plot3d solution file to read from? '
c     read(5,'(a60)') qfile
c
      gfile='plot3dg.bin'
      qfile='plot3dq.bin'
c
      open(unit=3,file=gfile,form='unformatted',status='old')
      open(unit=4,file=qfile,form='unformatted',status='old')
c
c     read grid file to get number of blocks
c
      read(3) ngd
      if (ngd.gt.nbl) then
         write(6,*)'stopping...must set parameter nmax to at least ',
     .   nbl
         stop
      end if
c
c     loop over all blocks
c
      do 100 n=1,ngd
c
c     read grid file
c
      if (n.eq.1) read(3) (idim(nn),jdim(nn),nn=1,ngd)
      kdim(n)=1
      if (n.eq.1) then
         idmax = id
         jdmax = jd
         kdmax = kd
         do 5 nn=1,ngd
         if (idim(nn).gt.idmax) idmax = idim(nn)
         if (jdim(nn).gt.jdmax) jdmax = jdim(nn)
         if (kdim(nn).gt.kdmax) kdmax = kdim(nn)
 5       continue
         if (idmax.gt.id .or. jdmax.gt.jd .or. kdmax.gt.kd) then
            if (idmax.gt.id) then
               write(6,*)'stopping: must set parameter id to ',
     .                   'at least ',idmax
            end if
            if (jdmax.gt.jd) then
               write(6,*)'stopping: must set parameter jd to ',
     .                   'at least ',jdmax
            end if
            if (kdmax.gt.kd) then
               write(6,*)'stopping: must set parameter kd to ',
     .                   'at least ',kdmax
            end if
            stop
         end if
      end if
      ni=idim(n)
      nj=jdim(n)
      nk=kdim(n)
      read(3) (((x(i,j,k,1),i=1,ni),j=1,nj),k=1,nk),
     .        (((x(i,j,k,2),i=1,ni),j=1,nj),k=1,nk),
     .        (((idum,i=1,ni),j=1,nj),k=1,nk)
c
c     read solution file
c
      if (n.eq.1) read(4) ngd
      if (n.eq.1) read(4) (idim(nn),jdim(nn),nn=1,ngd)
      read(4) fsmach,alpha,re,time
      ni=idim(n)
      nj=jdim(n)
      nk=kdim(n)
      read(4) ((((q(i,j,k,l),i=1,ni),j=1,nj),k=1,nk),l=1,4)
c
c     convert plot3d quantities to primitive
c
      do 110 j=1,jdim(n)
      do 110 k=1,kdim(n)
      do 110 i=1,idim(n)
      rho  = q(i,j,k,1)
      u    = q(i,j,k,2)/rho
      v    = q(i,j,k,3)/rho
      vsq  = u*u + v*v
      p    = (gamma-1.)*(q(i,j,k,4)-.5*rho*vsq)
      csq  = gamma*p/rho
      dmsq = vsq/csq
      mach = sqrt(dmsq)
      p0  = p*(1.+.5*(gamma-1.)*dmsq)**(gamma/(gamma-1.))
      q(i,j,k,2) = u
      q(i,j,k,3) = v
      q(i,j,k,4) = p
      q(i,j,k,5) = mach
      q(i,j,k,6) = p0
 110  continue
c
c     interpolate velocity profiles to experimental x-stations
c
      do 200 ns=1,nsta
         iloc = -1
         do 210 i=1,idim(n)-1
         if(x(i,1,1,1) .lt. xprofile(ns) .and. 
     .      x(i+1,1,1,1) .ge. xprofile(ns)) iloc = i
 210     continue
         if (iloc.eq.-1) then
            write(6,*)'stopping...cannot find x location in ',
     .      'grid for data station ',n
            stop
         end if
         x2 = x(iloc+1,1,1,1)
         x1 = x(iloc,1,1,1)
         do 220 j=1,jdim(n)
         vel(ns,j) = q(iloc,j,1,2) + (xprofile(ns)-x1)/(x2-x1)*
     .               (q(iloc+1,j,1,2)-q(iloc,j,1,2))
         y(ns,j)   = x(iloc,j,1,2) + (xprofile(ns)-x1)/(x2-x1)*
     .               (x(iloc+1,j,1,2)-x(iloc,j,1,2))
c        make velocity dimensional
         vel(ns,j) = vel(ns,j)*astar*0.3048 !0.3048 converts ft/s to m/s
 220     continue 
c        make y dimensionless wrt local duct height
         ymax = y(ns,jdim(n))
         do j=1,jdim(n)
            y(ns,j) = y(ns,j)/ymax
         end do
 200  continue
c
c     top and bottom wall pressures (normalized with total pressure)
c     q(4) = p  normalized so that q(4) = 1/gamma at cfl3d reference
c     conditions
c
      do 230 i=1,idim(n)
      ptop(i) = gamma*q(i,jdim(n),1,4)*pstar/ptotal
      pbot(i) = gamma*q(i,1,1,4)*pstar/ptotal
 230  continue
c
 100  continue
c
c     read in exp. data files (tecplot), and write out to new files so that
c     the computational data can be appended to it.
c
      write(6,*)'input 0 for weak-shock case, 1 for strong-shock case'
      read(5,*) icase
      if (icase.eq.0) then
         open(unit=10,file='sajwvel.data',form='formatted',
     .   status='old')
         open(unit=11,file='sajwpres.data',form='formatted',
     .   status='old')
      else
         open(unit=10,file='sajsvel.data',form='formatted',
     .   status='old')
         open(unit=11,file='sajspres.data',form='formatted',
     .   status='old')
      end if
c
      if (icase.eq.0) then
         open(unit=8,file='u_vs_y_ws.dat',form='formatted',
     .   status='unknown')
         open(unit=9,file='p_vs_x_ws.dat',form='formatted',
     .   status='unknown')
      else
         open(unit=8,file='u_vs_y_ss.dat',form='formatted',
     .   status='unknown')
         open(unit=9,file='p_vs_x_ss.dat',form='formatted',
     .   status='unknown')
      end if
c
      do nn=1,9999
         read(10,'(a80)',end=998) string
         write(8,'(a80)') string
      end do
 998  continue
      do nn=1,9999
         read(11,'(a80)',end=999) string
         write(9,'(a80)') string
      end do
 999  continue
c
c     write velocity profile data file (tecplot format)
c
      write(8,'(''Zone T="Computation", F=Point'')')
      do j=1,jdim(1)
         write(8,'(i3,f8.4,f8.2,f8.4,f8.2,f8.4,f8.2,f8.4,f8.2)') 
     .   j,(y(ns,j),vel(ns,j),ns=1,nsta)
      end do
c
c     write wall pressure data file (tecplot format)
c
      write(9,'(''Zone T="Computation", F=Point'')')
      do i=1,idim(1)
         write(9,*) x(i,jdim(1),1,1),pbot(i),x(i,jdim(1),1,1),ptop(i)
      end do
c
      write(6,*)'done'
      if (icase.eq.0) then
         write(6,*) 'tecplot file for velocity: u_vs_y_ws.dat'
         write(6,*) 'tecplot file for pressure: p_vs_x_ws.dat'
      else
         write(6,*) 'tecplot file for velocity: u_vs_y_ss.dat'
         write(6,*) 'tecplot file for pressure: p_vs_x_ss.dat'
      end if

      stop
      end 
