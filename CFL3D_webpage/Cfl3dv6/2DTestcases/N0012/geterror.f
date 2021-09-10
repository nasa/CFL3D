      program geterror
c
c   get the error to investigate order of spatial accuracy
c   (hardwired specifically for this case)
c
      dimension cd(5),dh(5)
c
c   first, write dh vs. error to show slopes
c   1st order:
      open(11,file='slope1.dat',form='formatted',status='unknown')
      dh1=1.e-3
      error1=1.e-6
      dh2=1.
      error2=1.e-3
      write(11,'(''variables="dh","error"'')')
      write(11,'(2e17.7)') dh1,error1
      write(11,'(2e17.7)') dh2,error2
c   2nd order:
      open(12,file='slope2.dat',form='formatted',status='unknown')
      dh1=1.e-3
      error1=1.e-6
      dh2=1.
      error2=1.
      write(12,'(''variables="dh","error"'')')
      write(12,'(2e17.7)') dh1,error1
      write(12,'(2e17.7)') dh2,error2
c
c   next, read the final results from each of the runs with different dh
      open(13,file='n0012_65t0.res',form='formatted',status='old')
      open(14,file='n0012_129t0.res',form='formatted',status='old')
      open(15,file='n0012_257t0.res',form='formatted',status='old')
      open(16,file='n0012_513t0.res',form='formatted',status='old')
      open(17,file='n0012_1025t0.res',form='formatted',status='old')
      dh(1)=sqrt(1./(65.*33.))
      dh(2)=sqrt(1./(129.*65.))
      dh(3)=sqrt(1./(257.*129.))
      dh(4)=sqrt(1./(513.*257.))
      dh(5)=sqrt(1./(1025.*513.))
      do n=13,17
        read(n,*)
        read(n,*)
        read(n,*)
        read(n,*)
        read(n,*) num
        do m=1,num
          read(n,*) idum,dum,dum,cd(n-12),dum,dum
        enddo
      enddo
c
c   use finest 2 results with Richardson extrapolation (assuming 2nd order)
c   to get an estimate for results with an infinitely small dh
      cdexact = 4.*cd(5)/3. - cd(4)/3.
      write(6,'('' extrapolated cdexact is '',e18.8)') cdexact
c
c   finally, write out results
      open(19,file='cderror.dat',form='formatted',status='unknown')
      write(19,'(''variables="dh","error","percent"'')')
      do n=1,5
        write(19,'(3e17.7)') dh(n),abs(cdexact-cd(n)),
     +    abs(cdexact-cd(n))/cdexact
      enddo
      write(6,'('' output to cderror.dat'')')
      write(6,'('' (slope curves for reference in slope1.dat'',
     +  '' and slope2.dat)'')')
      stop
      end
