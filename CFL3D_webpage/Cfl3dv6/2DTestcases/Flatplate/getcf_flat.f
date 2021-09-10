      program getcf_flat
c 
c  Reads cf from cfl3d.prout (from flat plate run using grdflat5.inp)
c  and compares to theory (eqn 6-118 from White, 1974).
c  This program is HARDWIRED for the specific flatplate
c  testcase grdflat5.inp
c
      dimension x(63),cf(63),cftheory(63)
c
      open(2,file='cfl3d.prout',form='formatted',status='old')
      open(3,file='cf.dat',form='formatted',status='unknown')
c  read front matter that we do not need
      do n=1,77
        read(2,*)
      enddo
c  Reynolds number for this case is 6 million, based on plate length
      re=6.e6
c  read cf
      write(3,'(''variables="x","c_f"'')')
      write(3,'(''zone,t="CFD"'')')
      do n=1,63
        read(2,*) idum,idum,idum,x(n),dum,dum,dum,dum,dum,cf(n)
c  write results only for x > 0
        if (x(n) .gt. 0.) then
          write(3,'(2e15.5)') x(n),cf(n)
        end if
      enddo
      write(3,'(''zone,t="theory"'')')
      do n=1,63
        if (x(n) .gt. 0.) then
c  get rex
          rex=re*x(n)
c  apply White formula
          cftheory(n)=0.025*(rex**(-1./7.))
c  write results
          write(3,'(2e15.5)') x(n),cftheory(n)
        end if
      enddo
c
      write(6,'('' tecplot file output to cf.dat'')')
      stop
      end
