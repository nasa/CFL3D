      program geterror
c
c   get the error to investigate order of temporal accuracy
c   (hardwired specifically for this case)
c
      dimension cl(5),cd(5),dt(5)
c
c   first, write dt vs. error to show slopes
c   1st order:
      open(11,file='slope1.dat',form='formatted',status='unknown')
      dt1=1.e-3
      error1=1.e-6
      dt2=1.
      error2=1.e-3
      write(11,'(2e17.7)') dt1,error1
      write(11,'(2e17.7)') dt2,error2
c   2nd order:
      open(12,file='slope2.dat',form='formatted',status='unknown')
      dt1=1.e-3
      error1=1.e-6
      dt2=1.
      error2=1.
      write(12,'(2e17.7)') dt1,error1
      write(12,'(2e17.7)') dt2,error2
c
c   next, read the final results from each of the runs with different dt
      open(13,file='cfl3d.res_0.2',form='formatted',status='old')
      open(14,file='cfl3d.res_0.1',form='formatted',status='old')
      open(15,file='cfl3d.res_0.05',form='formatted',status='old')
      open(16,file='cfl3d.res_0.025',form='formatted',status='old')
      open(17,file='cfl3d.res_0.0125',form='formatted',status='old')
      dt(1)=0.2
      dt(2)=0.1
      dt(3)=0.05
      dt(4)=0.025
      dt(5)=0.0125
      do n=13,17
        read(n,*)
        read(n,*)
        read(n,*)
        read(n,*)
        read(n,*) num
        do m=1,num
          read(n,*) idum,dum,cl(n-12),cd(n-12),dum,dum
        enddo
      enddo
c
c   use finest 2 results with Richardson extrapolation (assuming 2nd order)
c   to get an estimate for results with an infinitely small dt
      clexact = 4.*cl(5)/3. - cl(4)/3.
      cdexact = 4.*cd(5)/3. - cd(4)/3.
c
c   finally, write out results
      open(18,file='clerror.dat',form='formatted',status='unknown')
      open(19,file='cderror.dat',form='formatted',status='unknown')
      do n=1,5
        write(18,'(2e17.7)') dt(n),abs(clexact-cl(n))
        write(19,'(2e17.7)') dt(n),abs(cdexact-cd(n))
      enddo
      write(6,'('' output to clerror.dat and cderror.dat'')')
      write(6,'('' (slope curves for reference in slope1.dat'',
     +  '' and slope2.dat)'')')
      stop
      end
