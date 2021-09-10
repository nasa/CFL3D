      program getalpha
c
c   read cfl3d.res file and create alpha vs cl and cm plot
c   this code is currently HARDWIRED for this particular
c   pitching airfoil case
c
      pi=acos(-1.)
      open(2,file='cfl3d.res',form='formatted',status='old')
      do n=1,5
        read(2,*)
      enddo
c
c  iss = number of steady-state iterations:
      iss=1500
c
      do n=1,iss
        read(2,*)
      enddo
c
      open(3,file='alpha.dat',form='formatted',status='unknown')
      write(3,'(''variables="alpha, deg","C_L","C_M"'')')
c  input parameters for this particular case, to extract alpha
c  from iteration number:
      alpha0=4.86
      alpha1=2.44
      zkr=.01547
      dt=0.4
c
      do n=1,100000
        read(2,*,end=999) it,res,cl,cd,cy,cm
        time=float(it-iss)*dt
        alpha=alpha0+alpha1*sin(2.*pi*zkr*time)
        write(3,'(3e15.5)') alpha,cl,cm
      enddo
 999  continue
      write(6,'('' output to alpha.dat'')')
      stop
      end
