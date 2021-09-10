      program getcp_bump
c 
c  Reads cp from cfl3d.prout (from bump run using bumpperiodic.inp)
c  This program is HARDWIRED for the specific bump
c  testcase bumpperiodic.inp
c
      dimension x(181),cp(181)
c
      open(2,file='cfl3d.prout',form='formatted',status='old')
      open(3,file='cp.dat',form='formatted',status='unknown')
c  read front matter that we do not need
      do n=1,10
        read(2,*)
      enddo
c  read cp
      write(3,'(''variables="x","C_p"'')')
      write(3,'(''zone,t="CFD"'')')
      do n=1,181
        read(2,*) idum,idum,idum,x(n),dum,dum,dum,dum,dum,dum,dum,
     +   dum,cp(n)
c  write results
        write(3,'(2e15.5)') x(n),cp(n)
      enddo
c
      write(6,'('' tecplot file output to cp.dat'')')
      stop
      end
